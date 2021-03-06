/*------------------------------------------------------------------------*/
/*  Copyright 2019 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "edge_kernels/ContinuityEdgeSolverAlg.h"
#include "utils/StkHelpers.h"

namespace sierra {
namespace nalu {

ContinuityEdgeSolverAlg::ContinuityEdgeSolverAlg(
  Realm& realm,
  stk::mesh::Part* part,
  EquationSystem* eqSystem
) : AssembleEdgeSolverAlgorithm(realm, part, eqSystem)
{
  const auto& meta = realm.meta_data();

  coordinates_ = get_field_ordinal(meta, realm.get_coordinates_name());
  const std::string velField = realm.does_mesh_move()? "velocity_rtm" : "velocity";
  velocityRTM_ = get_field_ordinal(meta, velField);
  densityNp1_ = get_field_ordinal(meta, "density", stk::mesh::StateNP1);
  pressure_ = get_field_ordinal(meta, "pressure");
  Gpdx_ = get_field_ordinal(meta, "dpdx");
  edgeAreaVec_ = get_field_ordinal(meta, "edge_area_vector", stk::topology::EDGE_RANK);
  Udiag_ = get_field_ordinal(meta, "momentum_diag");
}

void
ContinuityEdgeSolverAlg::execute()
{
  const int ndim = realm_.meta_data().spatial_dimension();

  // Non-orthogonal correction factor for continuity equation system
  const std::string dofName = "pressure";
  const DblType nocFac
    = (realm_.get_noc_usage(dofName) == true) ? 1.0 : 0.0;

  // Classic Nalu projection timescale
  const DblType dt = realm_.get_time_step();
  const DblType gamma1 = realm_.get_gamma1();
  const DblType tauScale = dt / gamma1;

  // Interpolation option for rho*U
  const DblType interpTogether = realm_.get_mdot_interp();
  const DblType om_interpTogether = (1.0 - interpTogether);

  // STK ngp::Field instances for capture by lambda
  const auto& fieldMgr = realm_.ngp_field_manager();
  const auto coordinates = fieldMgr.get_field<double>(coordinates_);
  const auto velocity = fieldMgr.get_field<double>(velocityRTM_);
  const auto Gpdx = fieldMgr.get_field<double>(Gpdx_);
  const auto density = fieldMgr.get_field<double>(densityNp1_);
  const auto pressure = fieldMgr.get_field<double>(pressure_);
  const auto udiag = fieldMgr.get_field<double>(Udiag_);
  const auto edgeAreaVec = fieldMgr.get_field<double>(edgeAreaVec_);

  run_algorithm(
    realm_.bulk_data(),
    KOKKOS_LAMBDA(
      ShmemDataType& smdata,
      const stk::mesh::FastMeshIndex& edge,
      const stk::mesh::FastMeshIndex& nodeL,
      const stk::mesh::FastMeshIndex& nodeR)
    {
      // Scratch work array for edgeAreaVector
      NALU_ALIGNED DblType av[nDimMax_];
      // Populate area vector work array
      for (int d=0; d < ndim; ++d)
        av[d] = edgeAreaVec.get(edge, d);

      const DblType pressureL = pressure.get(nodeL, 0);
      const DblType pressureR = pressure.get(nodeR, 0);

      const DblType densityL = density.get(nodeL, 0);
      const DblType densityR = density.get(nodeR, 0);

      const DblType udiagL = udiag.get(nodeL, 0);
      const DblType udiagR = udiag.get(nodeR, 0);

      const DblType projTimeScale = 0.5 * (1.0/udiagL + 1.0/udiagR);
      const DblType rhoIp = 0.5 * (densityL + densityR);

      DblType axdx = 0.0;
      DblType asq = 0.0;
      for (int d=0; d < ndim; ++d) {
        const DblType dxj = coordinates.get(nodeR, d) - coordinates.get(nodeL, d);
        asq += av[d] * av[d];
        axdx += av[d] * dxj;
      }
      const DblType inv_axdx = 1.0 / axdx;

      DblType tmdot = -projTimeScale * (pressureR - pressureL) * asq * inv_axdx;
      for (int d = 0; d < ndim; ++d) {
        const DblType dxj =
          coordinates.get(nodeR, d) - coordinates.get(nodeL, d);
        // non-orthogonal correction
        const DblType kxj = av[d] - asq * inv_axdx * dxj;
        const DblType rhoUjIp = 0.5 * (densityR * velocity.get(nodeR, d) +
                                       densityL * velocity.get(nodeL, d));
        const DblType ujIp = 0.5 * (velocity.get(nodeR, d) + velocity.get(nodeL, d));
        const DblType GjIp = 0.5 * (Gpdx.get(nodeR, d) / udiagR +
                                    Gpdx.get(nodeL, d) / udiagL);
        tmdot += (interpTogether * rhoUjIp +
                  om_interpTogether * rhoIp * ujIp + GjIp) * av[d]
          - kxj * GjIp * nocFac;
      }
      tmdot /= tauScale;
      const DblType lhsfac = -asq * inv_axdx * projTimeScale / tauScale;

      // Left node entries
      smdata.lhs(0, 0) = -lhsfac;
      smdata.lhs(0, 1) = +lhsfac;
      smdata.rhs(0) = -tmdot;

      // Right node entries
      smdata.lhs(1, 0) = +lhsfac;
      smdata.lhs(1, 1) = -lhsfac;
      smdata.rhs(1) = tmdot;
    });
}

}  // nalu
}  // sierra
