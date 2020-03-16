// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//


#include "node_kernels/SDRSSTSASNodeKernel.h"
#include "Realm.h"
#include "SolutionOptions.h"
#include "SimdInterface.h"
#include "utils/StkHelpers.h"
#include "utils/ComputeVectorLaplacian.h"

#include "stk_mesh/base/MetaData.hpp"

namespace sierra {
namespace nalu {

SDRSSTSASNodeKernel::SDRSSTSASNodeKernel(
  const stk::mesh::MetaData& meta
) : NGPNodeKernel<SDRSSTSASNodeKernel>(),
    tkeID_(get_field_ordinal(meta, "turbulent_ke")),
    sdrID_(get_field_ordinal(meta, "specific_dissipation_rate")),
    densityID_(get_field_ordinal(meta, "density")),
    tviscID_(get_field_ordinal(meta, "turbulent_viscosity")),
    dudxID_(get_field_ordinal(meta, "dudx")),
    laplacianUID_(get_field_ordinal(meta,"laplacianU")),
    dkdxID_(get_field_ordinal(meta, "dkdx")),
    dwdxID_(get_field_ordinal(meta, "dwdx")),
    dualNodalVolumeID_(get_field_ordinal(meta, "dual_nodal_volume")),
    fOneBlendID_(get_field_ordinal(meta, "sst_f_one_blending")),
    nDim_(meta.spatial_dimension())
{}

void
SDRSSTSASNodeKernel::setup(Realm& realm)
{
  const auto& fieldMgr = realm.ngp_field_manager();

  tke_             = fieldMgr.get_field<double>(tkeID_);
  sdr_             = fieldMgr.get_field<double>(sdrID_);
  density_         = fieldMgr.get_field<double>(densityID_);
  tvisc_           = fieldMgr.get_field<double>(tviscID_);
  dudx_            = fieldMgr.get_field<double>(dudxID_);
  laplacianU_      = fieldMgr.get_field<double>(laplacianUID_);
  dkdx_            = fieldMgr.get_field<double>(dkdxID_);
  dwdx_            = fieldMgr.get_field<double>(dwdxID_);
  dualNodalVolume_ = fieldMgr.get_field<double>(dualNodalVolumeID_);
  fOneBlend_       = fieldMgr.get_field<double>(fOneBlendID_);

  const std::string dofName = "specific_dissipation_rate";
  relaxFac_ = realm.solutionOptions_->get_relaxation_factor(dofName);

  // Update turbulence model constants
  betaStar_ = realm.get_turb_model_constant(TM_betaStar);
  tkeProdLimitRatio_ = realm.get_turb_model_constant(TM_tkeProdLimitRatio);
  sigmaWTwo_ = realm.get_turb_model_constant(TM_sigmaWTwo);
  betaOne_ = realm.get_turb_model_constant(TM_betaOne);
  betaTwo_ = realm.get_turb_model_constant(TM_betaTwo);
  gammaOne_ = realm.get_turb_model_constant(TM_gammaOne);
  gammaTwo_ = realm.get_turb_model_constant(TM_gammaTwo);
  kappa_ = realm.get_turb_model_constant(TM_kappa);
  sasZetaTwo_ = realm.get_turb_model_constant(TM_sas_zetaTwo);
  sasSigmaPhi_ = realm.get_turb_model_constant(TM_sas_sigmaPhi);
  sasC_ = realm.get_turb_model_constant(TM_sas_C);
  sasCs_ = realm.get_turb_model_constant(TM_sas_Cs);
  oneOvercMuPower0p25_ = 1.0/sqrt(sqrt(betaStar_));
}

void
SDRSSTSASNodeKernel::execute(
  NodeKernelTraits::LhsType& lhs,
  NodeKernelTraits::RhsType& rhs,
  const stk::mesh::FastMeshIndex& node)
{
  using DblType = NodeKernelTraits::DblType;

  const DblType tke       = tke_.get(node, 0);
  const DblType sdr       = sdr_.get(node, 0);
  const DblType density   = density_.get(node, 0);
  const DblType tvisc     = tvisc_.get(node, 0);
  const DblType dVol      = dualNodalVolume_.get(node, 0);
  const DblType fOneBlend = fOneBlend_.get(node, 0);

  DblType magSqrGradK = 0.0;
  DblType magSqrGradOmega = 0.0;
  DblType magLaplU = 0.0;
  DblType S2 = 0.0;
  DblType Pk = 0.0;
  DblType crossDiff = 0.0;
  for (int i=0; i < nDim_; ++i) {
    magLaplU += laplacianU_.get(node,i) * laplacianU_.get(node,i);
    magSqrGradK += dkdx_.get(node,i) * dkdx_.get(node,i);
    magSqrGradOmega += dwdx_.get(node,i) * dwdx_.get(node,i);
    crossDiff += dkdx_.get(node, i) * dwdx_.get(node, i);
    const int offset = nDim_ * i;
    for (int j=0; j < nDim_; ++j) {
      const auto dudxij = dudx_.get(node, offset+j);
      S2 += dudxij * (dudxij + dudx_.get(node, j*nDim_ + i));
    }
  }
  magLaplU = sqrt(magLaplU);
  Pk = S2 * tvisc;

  const DblType Dk = betaStar_ * density * sdr * tke;

  // Clip production term
  Pk = stk::math::min(tkeProdLimitRatio_ * Dk, Pk);

  // Blend constants for SDR
  const DblType omf1 = (1.0 - fOneBlend);
  const DblType beta = fOneBlend * betaOne_ + omf1 * betaTwo_;
  const DblType gamma = fOneBlend * gammaOne_ + omf1 * gammaTwo_;
  const DblType sigmaD = 2.0 * omf1 * sigmaWTwo_;

  const DblType small = 1.0e-18;
  const DblType L = oneOvercMuPower0p25_ * stk::math::sqrt(tke) / sdr;
  const DblType Delta = stk::math::cbrt(dVol);
  const DblType Lnuk = stk::math::max(kappa_ * sqrt(S2) / (magLaplU + small),  sasCs_ * stk::math::sqrt(kappa_ * sasZetaTwo_ / (beta/betaStar_ - gamma)) * Delta );
  const DblType Qsas = density * stk::math::max(sasZetaTwo_ * kappa_ * S2 * (L * L / (Lnuk * Lnuk)) - 2.0 * sasC_ * tke / sasSigmaPhi_ * stk::math::max(magSqrGradOmega/(sdr * sdr), magSqrGradK * magSqrGradK / (tke * tke)) , 0.0);
      
  // Production term with appropriate clipping of tvisc
  const DblType Pw = gamma * density * Pk / stk::math::max(tvisc, 1.0e-16);
  const DblType Dw = beta * density * sdr * sdr;
  const DblType Sw = sigmaD * density * crossDiff / sdr;

  rhs(0) += (Pw - Dw + Sw + Qsas) * dVol;
  lhs(0, 0) += (2.0 * beta * density * sdr + stk::math::max(Sw / sdr, 0.0)) *
               dVol;
}

} // namespace nalu
}  // sierra
