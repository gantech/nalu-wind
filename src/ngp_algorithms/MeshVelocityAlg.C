// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#include "ngp_algorithms/MeshVelocityAlg.h"
#include "BuildTemplates.h"
#include "master_element/MasterElement.h"
#include "master_element/MasterElementFactory.h"
#include "ngp_algorithms/ViewHelper.h"
#include "ngp_utils/NgpLoopUtils.h"
#include "ngp_utils/NgpFieldOps.h"
#include "Realm.h"
#include "ScratchViews.h"
#include "SolutionOptions.h"
#include "utils/StkHelpers.h"

namespace sierra {
namespace nalu {

template<typename AlgTraits>
MeshVelocityAlg<AlgTraits>::MeshVelocityAlg(
  Realm& realm,
  stk::mesh::Part* part
) : Algorithm(realm, part),
    elemData_(realm.meta_data()),
    modelCoords_(get_field_ordinal(realm.meta_data(), "coordinates")),
    currentCoords_(get_field_ordinal(realm.meta_data(), "current_coordinates")),
    meshDispNp1_(
      get_field_ordinal(realm.meta_data(), "mesh_displacement", stk::mesh::StateNP1)),
    meshDispN_(
      get_field_ordinal(realm.meta_data(), "mesh_displacement", stk::mesh::StateN)),
    faceVelMag_(
      get_field_ordinal(realm.meta_data(), "face_velocity_mag", stk::topology::ELEM_RANK)),
    sweptVolumeNp1_(
      get_field_ordinal(
        realm.meta_data(), "swept_face_volume", stk::mesh::StateNP1, stk::topology::ELEM_RANK)),
    sweptVolumeN_(
      get_field_ordinal(
        realm.meta_data(), "swept_face_volume", stk::mesh::StateN, stk::topology::ELEM_RANK)),
    meSCS_(MasterElementRepo::get_surface_master_element<AlgTraits>())
{
  elemData_.add_cvfem_surface_me(meSCS_);

  elemData_.add_coordinates_field(modelCoords_, AlgTraits::nDim_, MODEL_COORDINATES);
  elemData_.add_coordinates_field(currentCoords_, AlgTraits::nDim_, CURRENT_COORDINATES);
  elemData_.add_element_field(faceVelMag_, AlgTraits::numScsIp_);
  elemData_.add_element_field(sweptVolumeN_, AlgTraits::numScsIp_);
  elemData_.add_element_field(sweptVolumeNp1_, AlgTraits::numScsIp_);
  elemData_.add_gathered_nodal_field(meshDispNp1_, AlgTraits::nDim_);
  elemData_.add_gathered_nodal_field(meshDispN_, AlgTraits::nDim_);
  
  elemData_.add_master_element_call(SCS_AREAV, CURRENT_COORDINATES);

  meSCS_->general_shape_fcn(AlgTraits::numScsIp_, isoParCoords_, isoCoordsShapeFcn_);
}

template<typename AlgTraits>
void MeshVelocityAlg<AlgTraits>::execute()
{
  using ElemSimdDataType = sierra::nalu::nalu_ngp::ElemSimdData<ngp::Mesh>;
  const auto& meshInfo = realm_.mesh_info();
  const auto& meta = meshInfo.meta();
  const DoubleType dt = realm_.get_time_step();
  const DoubleType gamma1 = realm_.get_gamma1();
  const DoubleType gamma2 = realm_.get_gamma2();
  const DoubleType gamma3 = realm_.get_gamma3();  
  const auto& ngpMesh = meshInfo.ngp_mesh();
  const auto& fieldMgr = meshInfo.ngp_field_manager();
  auto faceVel = fieldMgr.template get_field<double>(faceVelMag_);
  auto ngpSweptVol = fieldMgr.template get_field<double>(sweptVolumeNp1_);
  const auto faceVelOps = nalu_ngp::simd_elem_field_updater(ngpMesh, faceVel);
  const auto sweptVolOps = nalu_ngp::simd_elem_field_updater(ngpMesh, ngpSweptVol);

  const auto modelCoordsID = modelCoords_;
  const auto currentCoordsID = currentCoords_;
  const auto meshDispNp1ID = meshDispNp1_;
  const auto meshDispNID = meshDispN_;
  const auto faceVelID = faceVelMag_;
  const auto sweptVolNID = sweptVolumeN_;
  const auto sweptVolNp1ID = sweptVolumeNp1_;
  const auto* meSCS = meSCS_;

  const stk::mesh::Selector sel = meta.locally_owned_part()
    & stk::mesh::selectUnion(partVec_)
    & !(realm_.get_inactive_selector());

  const std::string algName = "compute_mesh_vel_" + std::to_string(AlgTraits::topo_);
  nalu_ngp::run_elem_algorithm(
    algName, meshInfo, stk::topology::ELEM_RANK, elemData_, sel,
    KOKKOS_LAMBDA(ElemSimdDataType& edata) {
      auto& scrView = edata.simdScrView;
      const auto& meViews = scrView.get_me_views(CURRENT_COORDINATES);
      const auto& v_areav = meViews.scs_areav;

      const auto& mCoords = scrView.get_scratch_view_2D(modelCoordsID);
      const auto& cCoords = scrView.get_scratch_view_2D(currentCoordsID);
      const auto& dispNp1 = scrView.get_scratch_view_2D(meshDispNp1ID);
      const auto& dispN = scrView.get_scratch_view_2D(meshDispNID);
      const auto& sweptVolN = scrView.get_scratch_view_1D(sweptVolNID);

      DoubleType dx[19][AlgTraits::nDim_];

      for (int i=0; i < 19; i++) {
          for (int j=0; j < AlgTraits::nDim_; j++)
              dx[i][j]= 0.0;
          for (int k=0; k < AlgTraits::nodesPerElement_; k++ ) {
              const DoubleType r = isoCoordsShapeFcn_[i*AlgTraits::nodesPerElement_ + k];
              for (int j=0; j < AlgTraits::nDim_; j++)
                  dx[i][j] += r * (dispNp1(k,j) - dispN(k,j)) ;
          }
      }


      DoubleType ws_coords_n[AlgTraits::nDim_ * AlgTraits::nodesPerElement_];
      DoubleType ws_scs_area_n[AlgTraits::numScsIp_ * AlgTraits::nDim_];
      double t_coords_n[AlgTraits::nDim_ * AlgTraits::nodesPerElement_];
      double t_scs_area_n[AlgTraits::numScsIp_ * AlgTraits::nDim_];
      for (int k=0; k < AlgTraits::nodesPerElement_; k++) {
          for (int j=0; j < AlgTraits::nDim_; j++)
              ws_coords_n[k*AlgTraits::nDim_ + j] = mCoords(k,j) + dispN(k,j);
      }
      {
        double scs_error = 0.0;
        for (int is=0; is < edata.numSimdElems; ++is) {
          for (int i=0; i < AlgTraits::nDim_*AlgTraits::nodesPerElement_; ++i)
            t_coords_n[i] = stk::simd::get_data(ws_coords_n[i], is);
          meSCS_->determinant(1, &t_coords_n[0], &t_scs_area_n[0], &scs_error);
          for (int i=0; i < AlgTraits::nDim_*AlgTraits::nodesPerElement_; ++i)
            stk::simd::set_data(ws_scs_area_n[i], is, t_scs_area_n[i]);
        }
      }

      for (int ip =0; ip < AlgTraits::numScsIp_; ++ip) {

          DoubleType dx_cg[AlgTraits::nDim_];
          DoubleType vdmvb[AlgTraits::nDim_];
          DoubleType vcmva[AlgTraits::nDim_];
          DoubleType vrhs[AlgTraits::nDim_];
          
          const int na = scsFaceNodeMap_[ip][0];
          const int nb = scsFaceNodeMap_[ip][1];
          const int nc = scsFaceNodeMap_[ip][2];
          const int nd = scsFaceNodeMap_[ip][3];

          for (int j=0; j < AlgTraits::nDim_; j++) {
              dx_cg[j] = 0.25 * ( dx[na][j] + dx[nb][j] + dx[nc][j] + dx[nd][j] );
              vdmvb[j] = dx[nd][j] - dx[nb][j];
              vcmva[j] = dx[nc][j] - dx[na][j];
          }

          vrhs[0] = 0.5 * ( ws_scs_area_n[ip*AlgTraits::nDim_+0] + v_areav(ip,0) )
              - dt * ( vdmvb[1] * vcmva[2] - vdmvb[2] * vcmva[1] ) / 6.0;
          vrhs[1] = 0.5 * ( ws_scs_area_n[ip*AlgTraits::nDim_+1] + v_areav(ip,1) )
              - dt * ( vdmvb[2] * vcmva[0] - vdmvb[0] * vcmva[2] ) / 6.0;
          vrhs[2] = 0.5 * ( ws_scs_area_n[ip*AlgTraits::nDim_+2] + v_areav(ip,2) )
              - dt * ( vdmvb[0] * vcmva[1] - vdmvb[1] * vcmva[0] ) / 6.0;
          
          const DoubleType tmp = dx_cg[0] * vrhs[0] + dx_cg[1] * vrhs[1] + dx_cg[2] * vrhs[2];
          sweptVolOps(edata, ip) = tmp;
          faceVelOps(edata, ip) =  ( gamma1 * tmp + (gamma1 + gamma2)  * sweptVolN(ip) );
          
      }

      
    });
}

INSTANTIATE_KERNEL(MeshVelocityAlg)

}  // nalu
}  // sierra
