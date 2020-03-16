// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

/*
    Scale-adaptive-simulation - k-omega-SST-SAS model.
    References:
        Egorov, Y., & Menter F.R. (2008).
        Development and Application of SST-SAS Model in the DESIDER Project.
        Advances in Hybrid RANS-LES Modelling,
        Notes on Num. Fluid Mech. And Multidisciplinary Design,
        Volume 97, 261-270.

    The default model coefficients are
        Cmu = betaStar 0.09
        Cs             0.11;
        kappa          0.41;
        zeta2          3.51;
        sigmaPhi       2.0/3.0;
        C              2;
        delta          cubeRootVol;
    
    Other references for the SAS model may have different coefficients.

        Menter, F.R., Egorov, Y. 
        The Scale-Adaptive Simulation Method for Unsteady Turbulent 
        Flow Predictions. Part 1: Theory and Model Description. 
        Flow Turbulence Combust 85, 113–138 (2010). 
        https://doi.org/10.1007/s10494-010-9264-5

*/ 
#ifndef SDRSSTSASNODEKERNEL_H
#define SDRSSTSASNODEKERNEL_H

#include "node_kernels/NodeKernel.h"
#include "FieldTypeDef.h"

#include "stk_mesh/base/BulkData.hpp"
#include "stk_ngp/Ngp.hpp"

namespace sierra {
namespace nalu {

class Realm;

class SDRSSTSASNodeKernel : public NGPNodeKernel<SDRSSTSASNodeKernel>
{
public:
  SDRSSTSASNodeKernel(const stk::mesh::MetaData&);

  KOKKOS_FORCEINLINE_FUNCTION
  SDRSSTSASNodeKernel() = default;

  KOKKOS_FUNCTION
  virtual ~SDRSSTSASNodeKernel() = default;

  virtual void setup(Realm&) override;

  KOKKOS_FUNCTION
  virtual void execute(
    NodeKernelTraits::LhsType&,
    NodeKernelTraits::RhsType&,
    const stk::mesh::FastMeshIndex&) override;

private:
  ngp::Field<double> tke_;
  ngp::Field<double> sdr_;
  ngp::Field<double> density_;
  ngp::Field<double> tvisc_;
  ngp::Field<double> dudx_;
  ngp::Field<double> laplacianU_;
  ngp::Field<double> dkdx_;
  ngp::Field<double> dwdx_;
  ngp::Field<double> dualNodalVolume_;
  ngp::Field<double> fOneBlend_;

  unsigned tkeID_             {stk::mesh::InvalidOrdinal};
  unsigned sdrID_             {stk::mesh::InvalidOrdinal};
  unsigned densityID_         {stk::mesh::InvalidOrdinal};
  unsigned tviscID_           {stk::mesh::InvalidOrdinal};
  unsigned dudxID_            {stk::mesh::InvalidOrdinal};
  unsigned laplacianUID_      {stk::mesh::InvalidOrdinal};
  unsigned dkdxID_            {stk::mesh::InvalidOrdinal};
  unsigned dwdxID_            {stk::mesh::InvalidOrdinal};
  unsigned dualNodalVolumeID_ {stk::mesh::InvalidOrdinal};
  unsigned fOneBlendID_       {stk::mesh::InvalidOrdinal};

  NodeKernelTraits::DblType betaStar_;
  NodeKernelTraits::DblType tkeProdLimitRatio_;
  NodeKernelTraits::DblType sigmaWTwo_;
  NodeKernelTraits::DblType betaOne_;
  NodeKernelTraits::DblType betaTwo_;
  NodeKernelTraits::DblType gammaOne_;
  NodeKernelTraits::DblType gammaTwo_;
  NodeKernelTraits::DblType kappa_;
  NodeKernelTraits::DblType sasZetaTwo_;
  NodeKernelTraits::DblType sasSigmaPhi_;
  NodeKernelTraits::DblType sasC_;
  NodeKernelTraits::DblType sasCs_;
  NodeKernelTraits::DblType relaxFac_;
  NodeKernelTraits::DblType oneOvercMuPower0p25_;

  const int nDim_;
};

}  // nalu
}  // sierra


#endif /* SDRSSTSASNODEKERNEL_H */
