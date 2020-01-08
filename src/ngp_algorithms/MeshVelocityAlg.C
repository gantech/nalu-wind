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
    meshDisp_(get_field_ordinal(realm.meta_data(), "mesh_displacement")),
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

  elemData_.add_master_element_call(SCS_AREAV, CURRENT_COORDINATES);

  meSCS_->general_shape_fcn(AlgTraits::numScsIp_, isoParCoords_, isoCoordsShapeFcn_);
}

template<typename AlgTraits>
void MeshVelocityAlg<AlgTraits>::execute()
{
  using ElemSimdDataType = sierra::nalu::nalu_ngp::ElemSimdData<ngp::Mesh>;
  const auto& meshInfo = realm_.mesh_info();
  const auto& meta = meshInfo.meta();
  const auto& ngpMesh = meshInfo.ngp_mesh();
  const auto& fieldMgr = meshInfo.ngp_field_manager();
  auto faceVel = fieldMgr.template get_field<double>(faceVelMag_);
  auto ngpSweptVol = fieldMgr.template get_field<double>(sweptVolumeNp1_);
  const auto faceVelOps = nalu_ngp::simd_elem_field_updater(ngpMesh, faceVel);
  const auto sweptVolOps = nalu_ngp::simd_elem_field_updater(ngpMesh, ngpSweptVol);

  const auto modelCoordsID = modelCoords_;
  const auto currentCoordsID = currentCoords_;
  const auto meshDispID = meshDisp_;
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

      const auto& cCoords = scrView.get_scratch_view_2D(currentCoordsID);
      const auto& disp = scrView.get_scratch_view_2D(meshDispID);
    });
}

INSTANTIATE_KERNEL(MeshVelocityAlg)

}  // nalu
}  // sierra
