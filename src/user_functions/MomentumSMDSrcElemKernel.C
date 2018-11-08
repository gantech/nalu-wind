/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "user_functions/MomentumSMDSrcElemKernel.h"
#include "AlgTraits.h"
#include "master_element/MasterElement.h"
#include "SolutionOptions.h"
#include "TimeIntegrator.h"

// template and scratch space
#include "BuildTemplates.h"
#include "ScratchViews.h"

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra {
namespace nalu {

template<typename AlgTraits>
MomentumSMDSrcElemKernel<AlgTraits>::MomentumSMDSrcElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ElemDataRequests& dataPreReqs,
  bool lumped)
  : Kernel(),
    ipNodeMap_(sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_)->ipNodeMap()),
    cur_time_(0.0),
    u_infty_(2.0),
    A_(1.0),
    sigma_(0.5),
    alpha_(1.0),
    omega_(1.0),
    M_(1.0),
    C_(1.0),
    K_(10.0),
    mu_(1e-5)
{
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  coordinates_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

  MasterElement* meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_);

  if ( lumped ) {
    get_scv_shape_fn_data<AlgTraits>([&](double* ptr){meSCV->shifted_shape_fcn(ptr);}, v_shape_function_);
  }
  else {
    get_scv_shape_fn_data<AlgTraits>([&](double* ptr){meSCV->shape_fcn(ptr);}, v_shape_function_);
  }

  // add master elements
  dataPreReqs.add_cvfem_volume_me(meSCV);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_master_element_call(SCV_VOLUME, CURRENT_COORDINATES);
}

template<typename AlgTraits>
MomentumSMDSrcElemKernel<AlgTraits>::~MomentumSMDSrcElemKernel()
{}

template<typename AlgTraits>
void
MomentumSMDSrcElemKernel<AlgTraits>::setup(const TimeIntegrator& timeIntegrator)
{
    cur_time_ = timeIntegrator.get_current_time();
}

template<typename AlgTraits>
void
MomentumSMDSrcElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType**>& lhs,
  SharedMemView<DoubleType*>& rhs,
  ScratchViews<DoubleType>& scratchViews)
{

  // Forcing nDim = 3 instead of using AlgTraits::nDim_ here to avoid compiler
  // warnings when this template is instantiated for 2-D topologies.
  NALU_ALIGNED DoubleType w_scvCoords[3];

  SharedMemView<DoubleType**>& v_coordinates = scratchViews.get_scratch_view_2D(*coordinates_);
  SharedMemView<DoubleType*>& v_scv_volume = scratchViews.get_me_views(CURRENT_COORDINATES).scv_volume;

  DoubleType sigma2 = sigma_ * sigma_;
  DoubleType oneOverSigma2 = 1.0/sigma2 ;
  DoubleType oneOverSigma5 = oneOverSigma2 * oneOverSigma2 / sigma_ ;
  DoubleType sinomegat = stk::math::sin(omega_ * cur_time_);
  DoubleType sin2omegat = stk::math::sin(2.0 * omega_ * cur_time_);
  DoubleType cosomegat = stk::math::cos(omega_ * cur_time_);

  DoubleType us = u_infty_ * (1.0 + alpha_ * sinomegat/sigma_ * stk::math::exp( - alpha_ * sinomegat * alpha_ * sinomegat / (sigma_*sigma_)) );
  DoubleType aoa = stk::math::atan( -alpha_ * omega_*cosomegat/us );
  DoubleType aoa_deg = aoa * 180.0/M_PI;
  DoubleType cl_polyfit[16];
  cl_polyfit[0] = -4.76967067e-31;
  cl_polyfit[1] = -2.00718968e-29;
  cl_polyfit[2] = 6.07390811e-26;
  cl_polyfit[3] = 2.36792535e-24;
  cl_polyfit[4] = -3.12374084e-21;
  cl_polyfit[5] = -1.11110307e-19; 
  cl_polyfit[6] = 8.31416678e-17;   
  cl_polyfit[7] = 2.63597939e-15;
  cl_polyfit[8] = -1.22021153e-12; 
  cl_polyfit[9] = -3.31592104e-11;  
  cl_polyfit[10] = 9.78048623e-09;
  cl_polyfit[11] = 2.10958481e-07;
  cl_polyfit[12] = -4.05844281e-05;
  cl_polyfit[13] = -5.79150921e-04;
  cl_polyfit[14] = 7.25422559e-02;
  cl_polyfit[15] = 4.47575641e-01;
  DoubleType cl = cl_polyfit[0];
  for (int it=0; it < 15; it++)
      cl = cl + (cl * aoa_deg) + cl_polyfit[it+1];

  DoubleType fsiForce = M_PI * aoa * (us*us + alpha_*alpha_*omega_*omega_*cosomegat*cosomegat) * stk::math::cos(aoa);
  
  for (int ip=0; ip < AlgTraits::numScvIp_; ++ip) {

    // nearest node to ip
    const int nearestNode = ipNodeMap_[ip];

    // zero out
    for ( int j =0; j < AlgTraits::nDim_; ++j )
        w_scvCoords[j] = 0.0;

    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
        const DoubleType r = v_shape_function_(ip,ic);
        for ( int j = 0; j < AlgTraits::nDim_; ++j )
            w_scvCoords[j] += r*v_coordinates(ic,j);
    }

    DoubleType x = w_scvCoords[0];
    DoubleType y = w_scvCoords[1];

    DoubleType yrels = y - alpha_ * sinomegat;
    DoubleType expfn = stk::math::exp(-(x*x + y*y)*oneOverSigma2);
    DoubleType expfn_fsi = stk::math::exp(-(x*x + yrels*yrels)*oneOverSigma2);



    // Compute RHS contributions
    const DoubleType scV = v_scv_volume(ip);
    const int nnNdim = nearestNode * AlgTraits::nDim_;

    //Static case - No FSI    
//    rhs(nnNdim + 0) += (expfn*u_infty_*(sigma2*u_infty_*x*(-(expfn*sigma_) + 4.*y) + 4.*mu_*y*(-2.*sigma2 + x*x + y*y))) * oneOverSigma5 * scV ;

//    rhs(nnNdim + 1) += (expfn*u_infty_*(-1.*expfn*sigma2*sigma_*u_infty_*y + 2.*sigma2*u_infty_*(-x*x + y*y) - 4.*mu_*x*(-2.*sigma2 + x*x + y*y))) * oneOverSigma5 * scV;

    //Oscillating case - Still no FSI - What worked earlier
//    rhs(nnNdim + 0) += (u_infty_*((-1.*sigma2*sigma_*u_infty_*x) * expfn - 8.*mu_*sigma2*y + 4.*sigma2*u_infty_*x*y + 4.*mu_*x*x*y + 4.*mu_*y*y*y + alpha_*omega_*sigma2*(1.*sigma2 - 2.*y*y)*cosomegat + 12.*alpha_*alpha_*mu_*y*sinomegat * sinomegat - 4.*alpha_*alpha_*alpha_*mu_*sinomegat * sinomegat * sinomegat + 2.*alpha_*alpha_*omega_*sigma2*y*sin2omegat + alpha_*sinomegat*(8.*mu_*sigma2 - 4.*sigma2*u_infty_*x - 4.*mu_*x*x - 12.*mu_*y*y - 1.*alpha_*alpha_*omega_*sigma2*sin2omegat))) * expfn * oneOverSigma5 * scV;

//    rhs(nnNdim + 1) += (u_infty_*(8.*mu_*sigma2*x - 2.*sigma2*u_infty_*x*x - 4.*mu_*x*x*x - (1.*sigma2*sigma_*u_infty_*y) * expfn + 2.*sigma2*u_infty_*y*y - 4.*mu_*x*y*y + 2.*alpha_*omega_*sigma2*x*y*cosomegat + alpha_*((1.*sigma2*sigma_*u_infty_) * expfn - 4.*sigma2*u_infty_*y + 8.*mu_*x*y)*sinomegat + alpha_*alpha_*(2.*sigma2*u_infty_ - 4.*mu_*x)*sinomegat * sinomegat - 1.*alpha_*alpha_*omega_*sigma2*x*sin2omegat)) * expfn * oneOverSigma5 * scV;

    //Oscillating case - Still no FSI - Flipped vorticity
    // rhs(nnNdim + 0) += (u_infty_*((-1.*sigma2*sigma_*u_infty_*x) * expfn + 8.*mu_*sigma2*y - 4.*mu_*x*x*y - 4.*mu_*y*y*y + alpha_*omega_*sigma2*(-1.*sigma2 + 2.*y*y)*cosomegat - 12.*alpha_*alpha_*mu_*y*sinomegat * sinomegat + 4.*alpha_*alpha_*alpha_*mu_*sinomegat * sinomegat * sinomegat - 2.*alpha_*alpha_*omega_*sigma2*y*sin2omegat + alpha_*sinomegat*(-8.*mu_*sigma2 + 4.*mu_*x*x + 12.*mu_*y*y + 1.*alpha_*alpha_*omega_*sigma2*sin2omegat))) * expfn * oneOverSigma5 * scV;

    // rhs(nnNdim + 1) += (u_infty_*(-2.*sigma2*sigma2*u_infty_ - 8.*mu_*sigma2*x + 2.*sigma2*u_infty_*x*x + 4.*mu_*x*x*x - (1.*sigma2*sigma_*u_infty_*y) * expfn + 2.*sigma2*u_infty_*y*y + 4.*mu_*x*y*y - 2.*alpha_*omega_*sigma2*x*y*cosomegat + alpha_*((1.*sigma2*sigma_*u_infty_) * expfn - 4.*sigma2*u_infty_*y - 8.*mu_*x*y)*sinomegat + alpha_*alpha_*(2.*sigma2*u_infty_ + 4.*mu_*x)*sinomegat * sinomegat + 1.*alpha_*alpha_*omega_*sigma2*x*sin2omegat)) * expfn * oneOverSigma5 * scV;

    // Velocity field is not oscillating - Flipped vorticity
    rhs(nnNdim + 0) += u_infty_*(-((sigma2*sigma_*u_infty_*x)*expfn) + (8.*mu_*sigma2*y) + (mu_*y*(-4.*x*x - 4.*y*y)) ) * expfn * oneOverSigma5 * scV;

    rhs(nnNdim + 1) += u_infty_*(-((sigma2*sigma_*u_infty_*y)*expfn) + (-2.*sigma2*sigma2*u_infty_ + 4.*mu_*x*x*x + 4.*mu_*x*y*y + sigma2*(-8.*mu_*x + 2.*u_infty_*x*x + 2.*u_infty_*y*y) ) ) * expfn * oneOverSigma5 * scV;
        
    //Oscillating case - With FSI
    rhs(nnNdim + 1) += fsiForce * expfn_fsi * scV / (sigma_ * sigma_ * M_PI);

    //Didn't work
//    rhs(nnNdim + 0) += (A_*expfn*oneOverSigma2*(-2*A_*x + 3.*(-(A_*expfn) + u_infty_)*x - 2.*alpha_*omega_*stk::math::cos(omega_*cur_time_)*yrels - (2.*A_*expfn*x*yrels*yrels)*oneOverSigma2 + (2.*(-(A_*expfn) + u_infty_)*x*yrels*yrels)*oneOverSigma2 - mu_*(4.0 - (4.0*x*x)*oneOverSigma2 - (4*yrels*yrels)*oneOverSigma2 - (2*(x + (2*x*yrels*yrels)*oneOverSigma2))/3.)))*scV ;

//    rhs(nnNdim + 1) += (A_*expfn*oneOverSigma2*(1.*alpha_*omega_*x*stk::math::cos(omega_*cur_time_) - 2*A_*yrels - 1.*(-(A_*expfn) + u_infty_)*yrels + (2.*(-(A_*expfn) + u_infty_)*x*x*yrels)*oneOverSigma2 - (2.*alpha_*omega_*x*stk::math::cos(omega_*cur_time_)*yrels*yrels)*oneOverSigma2 - (4.*A_*expfn*x*x*yrels*yrels*yrels)*oneOverSigma4 - mu_*((12.0*x*yrels)*oneOverSigma2 - (4.0*x*x*x*yrels)*oneOverSigma4 - (4.0*x*yrels*yrels*yrels)*oneOverSigma4 - (2.0*(x + (2.0*x*yrels*yrels)*oneOverSigma2))/3.)))*scV;
    //+ (alpha_*sigma2*(C_*omega_*stk::math::cos(omega_*cur_time_) + (K_ - M_*omega_*omega_)*stk::math::sin(omega_*cur_time_)))/A_ 
  }
}

INSTANTIATE_KERNEL(MomentumSMDSrcElemKernel);

}  // nalu
}  // sierra
