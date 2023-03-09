// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#include <aero/actuator/ActuatorFunctorsVG.h>
#include <aero/actuator/UtilitiesActuator.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <NaluEnv.h>
#include <FieldTypeDef.h>
#include "utils/LinearInterpolation.h"
#include <cmath>
#include <string>
#include <ostream>

namespace sierra {
namespace nalu {

InterpActuatorDensityVG::InterpActuatorDensityVG(
  ActuatorBulkVG& actBulk, stk::mesh::BulkData& stkBulk)
  : actBulk_(actBulk),
    stkBulk_(stkBulk),
    coordinates_(stkBulk_.mesh_meta_data().get_field<VectorFieldType>(
      stk::topology::NODE_RANK, "coordinates")),
    density_(stkBulk_.mesh_meta_data().get_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "density"))
{
  actBulk_.density_.sync_host();
  actBulk_.density_.modify_host();
}

void
InterpActuatorDensityVG::operator()(int index) const
{
  auto rho = actBulk_.density_.view_host();
  auto localCoord = actBulk_.localCoords_;

  if (actBulk_.pointIsLocal_(index)) {

    stk::mesh::Entity elem = stkBulk_.get_entity(
      stk::topology::ELEMENT_RANK, actBulk_.elemContainingPoint_(index));

    const int nodesPerElem = stkBulk_.num_nodes(elem);

    // just allocate for largest expected size (hex27)
    double ws_coordinates[81], ws_density[81];

    // Check to make sure the size is sufficient
    ThrowAssert(81 >= 3 * nodesPerElem);

    actuator_utils::gather_field(
      3, &ws_coordinates[0], *coordinates_, stkBulk_.begin_nodes(elem),
      nodesPerElem);

    actuator_utils::gather_field_for_interp(
      1, &ws_density[0], *density_, stkBulk_.begin_nodes(elem), nodesPerElem);

    actuator_utils::interpolate_field(
      1, elem, stkBulk_, &(localCoord(index, 0)), &ws_density[0],
      &(rho(index)));
    rho(index) /= actBulk_.localParallelRedundancy_(index);
  }
}

void
ActVGWriteToFile(
  ActuatorBulkVG& actBulk, const ActuatorMetaVG& actMeta)
{
  if (!actMeta.has_output_file_)
    return;
  std::string filename = actMeta.output_filenames_[actBulk.localTurbineId_];
  ActDualViewHelper<ActuatorFixedMemSpace> helper;
  auto vel = helper.get_local_view(actBulk.velocity_);
  auto force = helper.get_local_view(actBulk.actuatorForce_);
  auto density = helper.get_local_view(actBulk.density_);
  const int offset = actBulk.turbIdOffset_.h_view(actBulk.localTurbineId_);

  if (actBulk.localTurbineId_ == NaluEnv::self().parallel_rank()) {
    std::ofstream outFile;
    // ThrowErrorIf(NaluEnv::self().parallel_rank()!=0);

    outFile.open(filename, std::ios_base::app);
    const int stop =
      offset + actMeta.numPointsTurbine_.h_view(actBulk.localTurbineId_);

    for (int index = offset; index < stop; ++index) {
      const int i = index - offset;
      // write cached stuff from earlier computations
      outFile << actBulk.output_cache_[i];
      outFile << vel(index, 0) << ", " << vel(index, 1) << ", " << vel(index, 2)
              << ", ";
      outFile << force(index, 0) << ", " << force(index, 1) << ", "
              << force(index, 2) << ", ";
      outFile << density(index) << std::endl;
      actBulk.output_cache_[i].clear();
    }
    outFile.close();
  }
}

ActVGAssignVel::ActVGAssignVel(ActuatorBulkVG& actBulk)
  : velocity_(helper_.get_local_view(actBulk.velocity_)),
    density_(helper_.get_local_view(actBulk.density_)),
    points_(helper_.get_local_view(actBulk.pointCentroid_)),
    offset_(helper_.get_local_view(actBulk.turbIdOffset_)),
    debug_output_(actBulk.debug_output_),
    turbId_(actBulk.localTurbineId_)
{
}

void
ActVGAssignVel::operator()(int index) const
{

  const int pointId = index - offset_(turbId_);
  auto vel = Kokkos::subview(velocity_, index, Kokkos::ALL);
  auto rho = Kokkos::subview(density_, index);

  // Use this to double check the velocities and point positions
  auto point = Kokkos::subview(points_, index, Kokkos::ALL);
  if (debug_output_)
    NaluEnv::self().naluOutput()
      << "Blade " << turbId_ // LCCOUT
      << " pointId: " << pointId << std::scientific << std::setprecision(5)
      << " point: " << point(0) << " " << point(1) << " " << point(2) << " "
      << " vel: " << vel(0) << " " << vel(1) << " " << vel(2) << " "
      << " rho: " << *rho.data() << std::endl;
  // Do nothing otherwise
}

void
ActVGComputeForce(
  ActuatorBulkVG& actBulk, const ActuatorMetaVG& actMeta)
{

  ActDualViewHelper<ActuatorFixedMemSpace> helper;
  helper.touch_dual_view(actBulk.actuatorForce_);

  auto density = helper.get_local_view(actBulk.density_);
  auto velocity = helper.get_local_view(actBulk.velocity_);
  auto force = helper.get_local_view(actBulk.actuatorForce_);
  auto offset = helper.get_local_view(actBulk.turbIdOffset_);

  auto areas = helper.get_local_view(actMeta.areas_);
  auto bvec = helper.get_local_view(actMeta.bvec_);
  auto tvec = helper.get_local_view(actMeta.tvec_);
  auto nvec = helper.get_local_view(actMeta.nvec_);

  const double cvg = actMeta.Cvg_;

  const int turbId = actBulk.localTurbineId_;

  const int debug_output = actBulk.debug_output_;

  //Loop 'index' over actuator points of all the vgs
  Kokkos::parallel_for(
    "ActVGComputeForce", actBulk.local_range_policy(),
    ACTUATOR_LAMBDA(int index) {
      const int localId = index - offset(turbId);

      auto pforce = Kokkos::subview(force, index, Kokkos::ALL);

      auto bvecl = Kokkos::subview(bvec, turbId, Kokkos::ALL);
      auto tvecl = Kokkos::subview(tvec, turbId, Kokkos::ALL);
      auto nvecl = Kokkos::subview(nvec, turbId, Kokkos::ALL);

      auto vel = Kokkos::subview(velocity, index, Kokkos::ALL);

      double rho = density(index);  //Scalar arrays don't need Kokkos::subview
      double area = areas(index);

      // Magnitude of velocity
      double velmag_sq =
          vel(0) * vel(0) + vel(1) * vel(1) + vel(2) * vel(2);
      double velmag =  sqrt(velmag_sq);

      // Following eq. 1 in Troldborg et al.
      double alpha =
          (vel(0) * tvecl(0) + vel(1) * tvecl(1) + vel(2) * tvecl(2)) *
          (vel(0) * nvecl(0) + vel(1) * nvecl(1) + vel(2) * nvecl(2)) /
          velmag_sq;

      // Magnitude of force
      double fmag = cvg * rho * area * velmag_sq * alpha;

      // Force direction following convention in Troldborg et al.
      double el[3];
      el[0] = (vel(1) * bvecl(2) - vel(2) * bvecl(1)) / velmag;
      el[1] = (vel(2) * bvecl(0) - vel(0) * bvecl(2)) / velmag;
      el[2] = (vel(0) * bvecl(1) - vel(1) * bvecl(0)) / velmag;

      // Set the pointForce
      pforce(0) = -fmag * el[0];
      pforce(1) = -fmag * el[1];
      pforce(2) = -fmag * el[2];

      if (debug_output)
        NaluEnv::self().naluOutput()
          << "Blade " << turbId // LCCOUT
          << " pointId: " << localId << std::setprecision(5)
          << " cvg = " << cvg
          << " alpha: " << alpha
          << " rho = " << rho
          << " vel: " << vel(0) << ", " << vel(1) << ", " << vel(2)
          << " area = " << areas(index)
          << " nvecl = " << nvecl(0) << ", " << nvecl(1) << ", " << nvecl(2)
          << " tvecl = " << tvecl(0) << ", " << tvecl(1) << ", " << tvecl(2)
          << " bvecl = " << bvecl(0) << ", " << bvecl(1) << ", " << bvecl(2)
          << " fmag = " << fmag
          << " el = " << el[0] << ", " << el[1] << ", " << el[2]
          << std::endl;

    });

  actuator_utils::reduce_view_on_host(force);
}

void
ActVGComputeThrustInnerLoop::operator()(
  const uint64_t,
  const double*,
  double* sourceTerm,
  const double,
  const double scvIp) const
{

  auto offsets = actBulk_.turbIdOffset_.view_host();

  if (NaluEnv::self().parallel_rank() < actBulk_.num_vgs_) {
    int turbId = NaluEnv::self().parallel_rank();
    auto vgforce = Kokkos::subview(actBulk_.vg_force_, turbId, Kokkos::ALL);

    double force_term[3];

    for (size_t i = 0; i < 3; i++) {
      force_term[i] = sourceTerm[i] * scvIp;
      vgforce(i) += force_term[i];
    }
  }
}

void
ActVGSpreadForceWhProjInnerLoop::preloop()
{
  actBulk_.actuatorForce_.sync_host();
}

void
ActVGSpreadForceWhProjInnerLoop::operator()(
  const uint64_t pointId,
  const double* nodeCoords,
  double* sourceTerm,
  const double dual_vol,
  const double scvIp) const
{

  auto pointCoords =
    Kokkos::subview(actBulk_.pointCentroid_.view_host(), pointId, Kokkos::ALL);

  auto pointForce =
    Kokkos::subview(actBulk_.actuatorForce_.view_host(), pointId, Kokkos::ALL);

  auto epsilon =
    Kokkos::subview(actBulk_.epsilon_.view_host(), pointId, Kokkos::ALL);

  /* auto orientation = Kokkos::subview( */
  /*   actBulk_.orientationTensor_.view_host(), pointId, Kokkos::ALL); */

  double distance[3] = {0, 0, 0};
  // double projectedDistance[3] = {0, 0, 0};
  // double projectedForce[3] = {0, 0, 0};

  actuator_utils::compute_distance(
    3, nodeCoords, pointCoords.data(), &distance[0]);

  /* // transform distance from Cartesian to blade coordinate system */
  /* for (size_t i = 0; i < 3; i++) { */
  /*   for (size_t j = 0; j < 3; j++) { */
  /*     projectedDistance[i] += distance[j] * orientation(i + j * 3); */
  /*   } */
  /* } */

  const double gauss = actuator_utils::Gaussian_projection(
    3, &distance[0], epsilon.data());

  /* for (size_t j = 0; j < 3; j++) { */
  /*   projectedForce[j] = gauss * pointForce(j); */
  /* } */

  for (size_t j = 0; j < 3; j++) {
    sourceTerm[j] += gauss * pointForce[j] * scvIp / dual_vol;
  }
}

} /* namespace nalu */
} /* namespace sierra */
