// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//


#include "ngp_algorithms/TurbViscKsgsAlg.h"
#include "ngp_utils/NgpLoopUtils.h"
#include "ngp_utils/NgpTypes.h"
#include "ngp_utils/NgpFieldManager.h"
#include "Realm.h"
#include "SolutionOptions.h"
#include "utils/StkHelpers.h"

#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/NgpMesh.hpp"

namespace sierra{
namespace nalu{

TurbViscKsgsAlg::TurbViscKsgsAlg(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType* tvisc
) : Algorithm(realm, part),
    tviscField_(tvisc),
    tke_(get_field_ordinal(realm.meta_data(), "turbulent_ke")),
    density_(get_field_ordinal(realm.meta_data(), "density")),
    tvisc_(tvisc->mesh_meta_data_ordinal()),
    dualNodalVolume_(get_field_ordinal(realm.meta_data(), "dual_nodal_volume")),
    stabLscale_(get_field_ordinal(realm.meta_data(), "stab_lscale")),
    dhdx_(get_field_ordinal(realm.meta_data(), "dhdx")),
    specificHeat_(get_field_ordinal(realm.meta_data(), "specific_heat")),
    cmuEps_(realm.get_turb_model_constant(TM_cmuEps)),
    nDim_(realm.meta_data().spatial_dimension())
{
    const SolutionOptions& solnOpts = *realm.solutionOptions_;
    beta_ = solnOpts.thermalExpansionCoeff_;
    const std::vector<double>& solnOptsGravity = solnOpts.get_gravity_vector(nDim_);
    for (int i = 0; i < nDim_; i++)
        gravity_[i] = solnOptsGravity[i];
}

void
TurbViscKsgsAlg::execute()
{
  using Traits = nalu_ngp::NGPMeshTraits<stk::mesh::NgpMesh>;

  const auto& meta = realm_.meta_data();

  stk::mesh::Selector sel = (
    meta.locally_owned_part() | meta.globally_shared_part())
    & stk::mesh::selectField(*tviscField_);

  const auto& meshInfo = realm_.mesh_info();
  const auto ngpMesh = meshInfo.ngp_mesh();
  const auto& fieldMgr = meshInfo.ngp_field_manager();
  auto tke = fieldMgr.get_field<double>(tke_);
  auto density = fieldMgr.get_field<double>(density_);
  auto dualNodalVolume = fieldMgr.get_field<double>(dualNodalVolume_);
  auto tvisc = fieldMgr.get_field<double>(tvisc_);
  auto stabLscale = fieldMgr.get_field<double>(stabLscale_);
  auto dhdx = fieldMgr.get_field<double>(dhdx_);
  auto specificHeat = fieldMgr.get_field<double>(specificHeat_);

  tke.sync_to_device();
  density.sync_to_device();
  dualNodalVolume.sync_to_device();
  tvisc.sync_to_device();
  stabLscale.sync_to_device();
  dhdx.sync_to_device();

  const DblType invDim = 1.0 / static_cast<double>(meta.spatial_dimension());
  const DblType cmuEps = cmuEps_;

  nalu_ngp::run_entity_algorithm(
    "TurbViscKsgsAlg",
    ngpMesh, stk::topology::NODE_RANK, sel,
    KOKKOS_LAMBDA(const Traits::MeshIndex& meshIdx) {

      const DblType filter = stk::math::pow(dualNodalVolume.get(meshIdx, 0),
                                            invDim);
      const DblType strat =
          beta_ *
          (dhdx.get(meshIdx, 0)*gravity_[0]
           + dhdx.get(meshIdx, 1)*gravity_[1]
           + dhdx.get(meshIdx, 2)*gravity_[2])/specificHeat.get(meshIdx, 0);
      if (strat > 1e-10)
          stabLscale.get(meshIdx, 0) = stk::math::min(filter, 0.75 *stk::math::sqrt(tke.get(meshIdx, 0)/ strat));
      else
          stabLscale.get(meshIdx, 0) = filter;

      tvisc.get(meshIdx, 0) = cmuEps*density.get(meshIdx, 0)*std::sqrt(tke.get(meshIdx, 0))*stabLscale.get(meshIdx,0);
    });
  tvisc.modify_on_device();
}

} // namespace nalu
} // namespace Sierra
