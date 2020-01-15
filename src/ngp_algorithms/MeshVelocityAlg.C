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
#include "master_element/Hex8GeometryFunctions.h"
#include "ngp_algorithms/ViewHelper.h"
#include "ngp_utils/NgpLoopUtils.h"
#include "ngp_utils/NgpFieldOps.h"
#include "Realm.h"
#include "ScratchViews.h"
#include "SolutionOptions.h"
#include "utils/StkHelpers.h"

#include <cmath>

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
  meSCS_->general_shape_fcn(19, isoParCoords_, isoCoordsShapeFcn_);
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
  const auto& ngpMesh = meshInfo.ngp_mesh();
  const auto& fieldMgr = meshInfo.ngp_field_manager();
  auto faceVel = fieldMgr.template get_field<double>(faceVelMag_);
  auto ngpSweptVol = fieldMgr.template get_field<double>(sweptVolumeNp1_);
  const auto faceVelOps = nalu_ngp::simd_elem_field_updater(ngpMesh, faceVel);
  const auto sweptVolOps = nalu_ngp::simd_elem_field_updater(ngpMesh, ngpSweptVol);

  const auto modelCoordsID = modelCoords_;
  const auto meshDispNp1ID = meshDispNp1_;
  const auto meshDispNID = meshDispN_;
  const auto sweptVolNID = sweptVolumeN_;

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
      const auto& dispNp1 = scrView.get_scratch_view_2D(meshDispNp1ID);
      const auto& dispN = scrView.get_scratch_view_2D(meshDispNID);
      const auto& sweptVolN = scrView.get_scratch_view_1D(sweptVolNID);

      /* DoubleType dx[19][AlgTraits::nDim_]; */
      DoubleType scs_coords_n[19][AlgTraits::nDim_];
      DoubleType scs_coords_np1[19][AlgTraits::nDim_];
      
      for (int i=0; i < 19; i++) {
          for (int j=0; j < AlgTraits::nDim_; j++) {
              /* dx[i][j]= 0.0; */
              scs_coords_n[i][j] = 0.0;
              scs_coords_np1[i][j] = 0.0;
          }
          for (int k=0; k < AlgTraits::nodesPerElement_; k++ ) {
              const DoubleType r = isoCoordsShapeFcn_[i*AlgTraits::nodesPerElement_ + k];
              for (int j=0; j < AlgTraits::nDim_; j++) {
                  /* dx[i][j] += r * (dispNp1(k,j) - dispN(k,j)) ; */
                  scs_coords_n[i][j] += r * (mCoords(k,j) + dispN(k,j));
                  scs_coords_np1[i][j] += r * (mCoords(k,j) + dispNp1(k,j));
              }
          }
          /* std::cerr << "scs_coords_n " << i << " = " << scs_coords_n[i][0] << "," */
          /*           << scs_coords_n[i][1] << "," */
          /*           << scs_coords_n[i][2] << std::endl; */
      }

      /* DoubleType ws_coords_n[AlgTraits::nDim_ * AlgTraits::nodesPerElement_]; */
      /* DoubleType ws_scs_area_n[AlgTraits::numScsIp_ * AlgTraits::nDim_]; */
      /* double t_coords_n[AlgTraits::nDim_ * AlgTraits::nodesPerElement_]; */
      /* double t_scs_area_n[AlgTraits::numScsIp_ * AlgTraits::nDim_]; */
      /* for (int k=0; k < AlgTraits::nodesPerElement_; k++) { */
      /*     for (int j=0; j < AlgTraits::nDim_; j++) */
      /*         ws_coords_n[k*AlgTraits::nDim_ + j] = mCoords(k,j) + dispN(k,j); */
      /* } */

      /* { */
      /*   double scs_error = 0.0; */
      /*   for (int is=0; is < edata.numSimdElems; ++is) { */
      /*     for (int i=0; i < AlgTraits::nDim_*AlgTraits::nodesPerElement_; ++i) */
      /*       t_coords_n[i] = stk::simd::get_data(ws_coords_n[i], is); */
      /*     meSCS_->determinant(1, &t_coords_n[0], &t_scs_area_n[0], &scs_error); */
      /*     for (int i=0; i < AlgTraits::nDim_*AlgTraits::numScsIp_; ++i) */
      /*       stk::simd::set_data(ws_scs_area_n[i], is, t_scs_area_n[i]); */
      /*   } */
      /* } */

      /* const int elem_faces[6][4] = { */
      /*                                      {2,6,5,1}, */
      /*                                      {3,7,6,2}, */
      /*                                      {0,4,7,3}, */
      /*                                      {1,5,4,0}, */
      /*                                      {7,4,5,6}, */
      /*                                      {1,0,3,2} */
      /* }; */

      /* DoubleType elem_face_mesh_vel[6]; */
      /* DoubleType div_mesh_vel = 0.0; */
      
      /* for (int iface=0; iface < 6; ++iface) { */

      /*     DoubleType dx_cg[AlgTraits::nDim_]; */
      /*     DoubleType vdmvb[AlgTraits::nDim_]; */
      /*     DoubleType vamvc[AlgTraits::nDim_]; */
      /*     DoubleType vrhs[AlgTraits::nDim_]; */

      /*     DoubleType xamc_n[AlgTraits::nDim_]; */
      /*     DoubleType xdmb_n[AlgTraits::nDim_]; */
      /*     DoubleType xamc_np1[AlgTraits::nDim_]; */
      /*     DoubleType xdmb_np1[AlgTraits::nDim_]; */
          
      /*     DoubleType area_n[AlgTraits::nDim_]; */
      /*     DoubleType area_np1[AlgTraits::nDim_]; */

      /*     const int na = elem_faces[iface][0]; */
      /*     const int nb = elem_faces[iface][1]; */
      /*     const int nc = elem_faces[iface][2]; */
      /*     const int nd = elem_faces[iface][3]; */

      /*     for (int j=0; j < AlgTraits::nDim_; j++) { */
      /*         dx_cg[j] = 0.25 * ( (dispNp1(na,j)-dispN(na,j)) + (dispNp1(nb,j)-dispN(nb,j)) + (dispNp1(nc,j)-dispN(nc,j)) + (dispNp1(nd,j)-dispN(nd,j)) ); */
      /*         vdmvb[j] = (dispNp1(nd,j)-dispN(nd,j)) - (dispNp1(nb,j)-dispN(nb,j)); */
      /*         vamvc[j] = (dispNp1(na,j)-dispN(na,j)) - (dispNp1(nc,j)-dispN(nc,j)); */
      /*         xamc_n[j] = (mCoords(na,j) + dispN(na,j) - mCoords(nc,j) - dispN(nc,j)); */
      /*         xamc_np1[j] = (mCoords(na,j) + dispNp1(na,j) - mCoords(nc,j) - dispNp1(nc,j)); */
      /*         xdmb_n[j] = (mCoords(nd,j) + dispN(nd,j) - mCoords(nb,j) - dispN(nb,j)); */
      /*         xdmb_np1[j] = (mCoords(nd,j) + dispNp1(nd,j) - mCoords(nb,j) - dispNp1(nb,j)); */
      /*     } */
      /*     area_n[0] = 0.5 * (xdmb_n[1]*xamc_n[2] - xdmb_n[2]*xamc_n[1]); */
      /*     area_n[1] = 0.5 * (xdmb_n[2]*xamc_n[0] - xdmb_n[0]*xamc_n[2]); */
      /*     area_n[2] = 0.5 * (xdmb_n[0]*xamc_n[1] - xdmb_n[1]*xamc_n[0]); */
      /*     area_np1[0] = 0.5 * (xdmb_np1[1]*xamc_np1[2] - xdmb_np1[2]*xamc_np1[1]); */
      /*     area_np1[1] = 0.5 * (xdmb_np1[2]*xamc_np1[0] - xdmb_np1[0]*xamc_np1[2]); */
      /*     area_np1[2] = 0.5 * (xdmb_np1[0]*xamc_np1[1] - xdmb_np1[1]*xamc_np1[0]); */

      /*     elem_face_mesh_vel[iface] = dx_cg[0] * (0.5 * (area_n[0] + area_np1[0]) - (vdmvb[1]*vamvc[2] - vdmvb[2]*vamvc[1]) / 12.0 )  + */
      /*                                 dx_cg[1] * (0.5 * (area_n[1] + area_np1[1]) - (vdmvb[2]*vamvc[1] - vdmvb[1]*vamvc[2]) / 12.0 ) + */
      /*                                 dx_cg[2] * (0.5 * (area_n[2] + area_np1[2]) - (vdmvb[0]*vamvc[0] - vdmvb[0]*vamvc[0]) / 12.0); */
      /*     /\* std::cerr << "area_1 = " << area_n[0] << "," << area_n[1] << "," << area_n[2] << std::endl ; *\/ */
      /*     /\* std::cerr << "area_2 = " << area_np1[0] << "," << area_np1[1] << "," << area_np1[2] << std::endl ; *\/ */
      /*     /\* std::cerr << "dx_cg = " << dx_cg[0] << "," << dx_cg[1] << "," << dx_cg[2] << std::endl ; *\/ */
          
      /*     /\* std::cerr << "elem_face_mesh_vel[" << iface << "] = " << elem_face_mesh_vel[iface] << std::endl ; *\/ */

      /*     DoubleType face_svc[8][3]; */
      /*     for (int j=0; j < 3; j++) { */
              
      /*         face_svc[0][j] = mCoords(na,j) + dispN(na,j); */
      /*         face_svc[1][j] = mCoords(nb,j) + dispN(nb,j); */
      /*         face_svc[2][j] = mCoords(nc,j) + dispN(nc,j); */
      /*         face_svc[3][j] = mCoords(nd,j) + dispN(nd,j); */
      /*         face_svc[4][j] = mCoords(na,j) + dispNp1(na,j); */
      /*         face_svc[5][j] = mCoords(nb,j) + dispNp1(nb,j); */
      /*         face_svc[6][j] = mCoords(nc,j) + dispNp1(nc,j); */
      /*         face_svc[7][j] = mCoords(nd,j) + dispNp1(nd,j); */
              
      /*     } */
      /*     DoubleType tmp = hex_volume_grandy(face_svc); */
      /*     // std::cerr << "Grandy vol calc = " << tmp << std::endl; */
      /*     div_mesh_vel += tmp; //elem_face_mesh_vel[iface]; */
      /* } */
      /* double n_div_mesh_vel = stk::simd::get_data(div_mesh_vel,0); */
      /* if( std::abs(n_div_mesh_vel) > 1e-15) */
      /*     std::cerr << "Div Mesh Vel = " << div_mesh_vel << std::endl; */

      for (int ip =0; ip < AlgTraits::numScsIp_; ++ip) {

          /* DoubleType dx_cg[AlgTraits::nDim_]; */
          /* DoubleType vdmvb[AlgTraits::nDim_]; */
          /* DoubleType vamvc[AlgTraits::nDim_]; */
          /* DoubleType vrhs[AlgTraits::nDim_]; */

          /* DoubleType xcma_n[AlgTraits::nDim_]; */
          /* DoubleType xdmb_n[AlgTraits::nDim_]; */
          /* DoubleType xcma_np1[AlgTraits::nDim_];           */
          /* DoubleType xdmb_np1[AlgTraits::nDim_]; */
          
          /* DoubleType area_n[AlgTraits::nDim_]; */
          /* DoubleType area_np1[AlgTraits::nDim_]; */

          const int na = scsFaceNodeMap_[ip][0];
          const int nb = scsFaceNodeMap_[ip][1];
          const int nc = scsFaceNodeMap_[ip][2];
          const int nd = scsFaceNodeMap_[ip][3];

          /* for (int j=0; j < AlgTraits::nDim_; j++) { */
          /*     dx_cg[j] = 0.25 * ( dx[na][j] + dx[nb][j] + dx[nc][j] + dx[nd][j] ); */
          /*     vdmvb[j] = dx[nd][j] - dx[nb][j]; */
          /*     vamvc[j] = dx[na][j] - dx[nc][j]; */
          /*     xcma_n[j] = scs_coords_n[na][j] - scs_coords_n[nc][j]; */
          /*     xdmb_n[j] = scs_coords_n[nd][j] - scs_coords_n[nb][j]; */
          /*     xcma_np1[j] = scs_coords_np1[nc][j] - scs_coords_np1[na][j]; */
          /*     xdmb_np1[j] = scs_coords_np1[nd][j] - scs_coords_np1[nb][j]; */
          /* } */

          DoubleType scs_vol_coords[8][3];

          for (int j=0; j < AlgTraits::nDim_; j++) {
              scs_vol_coords[0][j] = scs_coords_n[na][j];
              scs_vol_coords[1][j] = scs_coords_n[nb][j];
              scs_vol_coords[2][j] = scs_coords_n[nc][j];
              scs_vol_coords[3][j] = scs_coords_n[nd][j];
              scs_vol_coords[4][j] = scs_coords_np1[na][j];
              scs_vol_coords[5][j] = scs_coords_np1[nb][j];
              scs_vol_coords[6][j] = scs_coords_np1[nc][j];
              scs_vol_coords[7][j] = scs_coords_np1[nd][j];
          }

          DoubleType tmp = hex_volume_grandy(scs_vol_coords);
          sweptVolOps(edata,ip) = tmp;

          /* area_n[0] = 0.5*(xcma_n[1]*xdmb_n[2] - xcma_n[2]*xdmb_n[1]); */
          /* area_n[1] = 0.5*(xcma_n[2]*xdmb_n[0] - xcma_n[0]*xdmb_n[2]); */
          /* area_n[2] = 0.5*(xcma_n[0]*xdmb_n[1] - xcma_n[1]*xdmb_n[0]); */
          /* area_np1[0] = 0.5*(xcma_np1[1]*xdmb_np1[2] - xcma_np1[2]*xdmb_np1[1]); */
          /* area_np1[1] = 0.5*(xcma_np1[2]*xdmb_np1[0] - xcma_np1[0]*xdmb_np1[2]); */
          /* area_np1[2] = 0.5*(xcma_np1[0]*xdmb_np1[1] - xcma_np1[1]*xdmb_np1[0]); */

          /* std::cerr << "ip = " << ip << ", "; */
          /* std::cerr << "dx = " << dx[na][0] << "," << dx[na][1] << "," << dx[na][2] << "," */
          /*           << dx[nb][0] << "," << dx[nb][1] << "," << dx[nb][2] << "," */
          /*           << dx[nc][0] << "," << dx[nc][1] << "," << dx[nc][2] << "," */
          /*           << dx[nd][0] << "," << dx[nd][1] << "," << dx[nd][2] << std::endl; */
          /* std::cerr << "old area n = " << ws_scs_area_n[ip*AlgTraits::nDim_+0] << "," << ws_scs_area_n[ip*AlgTraits::nDim_+1] << "," << ws_scs_area_n[ip*AlgTraits::nDim_+2] << std::endl; */
          /* std::cerr << "old area np1 = " << v_areav(ip,0) << "," << v_areav(ip,1) << "," << v_areav(ip,2) << std::endl; */

          /* vrhs[0] = 0.5 * ( ws_scs_area_n[ip*AlgTraits::nDim_+0] + v_areav(ip,0) ) */
          /*     - ( vdmvb[1] * vamvc[2] - vdmvb[2] * vamvc[1] ) / 12.0; */
          /* vrhs[1] = 0.5 * ( ws_scs_area_n[ip*AlgTraits::nDim_+1] + v_areav(ip,1) ) */
          /*     - ( vdmvb[2] * vamvc[0] - vdmvb[0] * vamvc[2] ) / 12.0; */
          /* vrhs[2] = 0.5 * ( ws_scs_area_n[ip*AlgTraits::nDim_+2] + v_areav(ip,2) ) */
          /*     - ( vdmvb[0] * vamvc[1] - vdmvb[1] * vamvc[0] ) / 12.0; */

          /* std::cerr << "Cross product terms = " << ( vdmvb[1] * vamvc[2] - vdmvb[2] * vamvc[1] ) / 12.0 */
          /*           << "," << ( vdmvb[2] * vamvc[0] - vdmvb[0] * vamvc[2] ) / 12.0 */
          /*           << "," << ( vdmvb[0] * vamvc[1] - vdmvb[1] * vamvc[0] ) / 12.0 */
          /*           << std::endl ; */
          
          /* std::cerr << "new area n = " << area_n[0] << "," << area_n[1] << "," << area_n[2] << std::endl; */
          /* std::cerr << "new area np1 = " << area_np1[0] << "," << area_np1[1] << "," << area_np1[2] << std::endl; */
          /* /\* vrhs[0] = 0.5 * ( area_n[0] + area_np1[0] ) *\/ */
          /* /\*     - dt * ( vdmvb[1] * vamvc[2] - vdmvb[2] * vamvc[1] ) / 12.0; *\/ */
          /* /\* vrhs[1] = 0.5 * ( area_n[1] + area_np1[1] ) *\/ */
          /* /\*     - dt * ( vdmvb[2] * vamvc[0] - vdmvb[0] * vamvc[2] ) / 12.0; *\/ */
          /* /\* vrhs[2] = 0.5 * ( area_n[2] + area_np1[2] ) *\/ */
          /* /\*     - dt * ( vdmvb[0] * vamvc[1] - vdmvb[1] * vamvc[0] ) / 12.0; *\/ */

          /* const DoubleType tmp = dx_cg[0] * vrhs[0] + dx_cg[1] * vrhs[1] + dx_cg[2] * vrhs[2]; */
          /* sweptVolOps(edata, ip) = tmp; */
          faceVelOps(edata, ip) =  ( gamma1 * tmp + (gamma1 + gamma2)  * sweptVolN(ip) )/dt;
          
      }

    });
}

INSTANTIATE_KERNEL(MeshVelocityAlg)

}  // nalu
}  // sierra
