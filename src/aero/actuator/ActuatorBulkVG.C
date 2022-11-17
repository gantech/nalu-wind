// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#include <aero/actuator/ActuatorBulkVG.h>
#include <aero/actuator/UtilitiesActuator.h>
#include <NaluEnv.h>

namespace sierra {
namespace nalu {

ActuatorMetaVG::ActuatorMetaVG(const int num_force_pts, const ActuatorMeta& actMeta)
  : ActuatorMeta(actMeta),
    isotropicGaussian_(false),
    num_force_pts_(num_force_pts),
    Cvg_(4.0),
    areas_("areas", numberOfActuators_ * num_force_pts),
    centers_("areas", numberOfActuators_ * num_force_pts),
    bvec_("bvec", numberOfActuators_),
    tvec_("tvec", numberOfActuators_),
    nvec_("nvec", numberOfActuators_),
    output_filenames_(numberOfActuators_),
    has_output_file_(false)
{
    //Areas and Centers not initialized yet
}

ActuatorBulkVG::ActuatorBulkVG(const ActuatorMetaVG& actMeta)
  : ActuatorBulk(actMeta),
    density_("actDensity", actMeta.numPointsTotal_),
    alpha_("actAngleOfAttack", actMeta.numPointsTotal_),
    vg_force_("vgForce", actMeta.numberOfActuators_),
    num_force_pts_(actMeta.num_force_pts_),
    assignedProc_("assignedProcBulk", actMeta.numberOfActuators_),
    num_vgs_(actMeta.numberOfActuators_),
    debug_output_(actMeta.debug_output_),
    output_cache_(
      actMeta.has_output_file_
        ? actMeta.numPointsTurbine_.h_view(localTurbineId_)
        : 0)

{
  // Allocate blades to turbines
  const int nProcs = NaluEnv::self().parallel_size();
  const int nVG = actMeta.numberOfActuators_;
  const int intDivision = nVG / nProcs;
  const int remainder = actMeta.numberOfActuators_ % nProcs;

  if (actMeta.debug_output_)
    NaluEnv::self().naluOutputP0()
      << " nProcs: " << nProcs << " nVG:  " << nVG
      << " intDiv: " << intDivision << " remain: " << remainder
      << std::endl; // LCCOUT

  if (remainder && intDivision) // this doesn't work for nProcs=1
    throw std::runtime_error(" ERRORXX: more blades than ranks");
  if (nVG > nProcs)
    throw std::runtime_error(" ERROR: more blades than ranks");

  for (int i = 0; i < nVG; i++) {
    assignedProc_.h_view(i) = i;
    NaluEnv::self().naluOutputP0()
      << " VG#: " << i << " Proc#: " << assignedProc_.h_view(i)
      << std::endl;
  }

  // Double check offsets
  if (actMeta.debug_output_)
    for (int i = 0; i < actMeta.numberOfActuators_; ++i) {
      NaluEnv::self().naluOutputP0()
        << "Offset blade: " << i << " " << turbIdOffset_.h_view(i)
        << " num_force_pts: " << num_force_pts_
        << std::endl; // LCCOUT
    }
  init_epsilon(actMeta);
  init_points(actMeta);
  add_output_headers(actMeta);
  NaluEnv::self().naluOutputP0()
    << "Done ActuatorBulkVG Init " << std::endl; // LCCOUT
}

void
ActuatorBulkVG::add_output_headers(const ActuatorMetaVG& actMeta)
{
  if (!actMeta.has_output_file_)
    return;
  std::string filename = actMeta.output_filenames_[localTurbineId_];

  if (localTurbineId_ == NaluEnv::self().parallel_rank()) {
    std::ofstream outFile;

    outFile.open(filename, std::ios_base::out);
    outFile << "pointId,alpha,forceX,forceY,forceZ,density\n";
    outFile.close();
  }
}

void
ActuatorBulkVG::init_epsilon(const ActuatorMetaVG& actMeta)
{
  // set epsilon and radius

  epsilon_.modify_host();
  searchRadius_.modify_host();

  const int nVGs = actMeta.n_vgs_;
  for (int iVG = 0; iVG < nVGs; iVG++) {
    // LCC test this for non-isotropic
    if (NaluEnv::self().parallel_rank() == assignedProc_.h_view(iVG)) {
      const size_t numForcePts = actMeta.num_force_pts_;
      const int offset = turbIdOffset_.h_view(iVG);
      auto areas = actMeta.areas_;
      for (int np = 0; np < numForcePts; np++) {
        auto epsilonLocal =
          Kokkos::subview(epsilon_.view_host(), np + offset, Kokkos::ALL);

        for (int i = 0; i < 3; i++) {
          // Define the optimal epsilon
            epsilonLocal(i) = std::sqrt(areas.h_view(np+offset));
        }
        // The radius of the searching. This is given in terms of
        //   the maximum of epsilon.x/y/z/.
        //
        // This is the length where the value of the Gaussian becomes
        // 0.1 % (1.0 / .001 = 1000) of the value at the center of the Gaussian
        searchRadius_.h_view(np + offset) =  epsilonLocal(0) * sqrt(log(1.e3));

      } // loop over np
    }
  } // loop over iVG

  actuator_utils::reduce_view_on_host(epsilon_.view_host());
  actuator_utils::reduce_view_on_host(searchRadius_.view_host());
  epsilon_.sync_host();
  searchRadius_.sync_host();
}

// Initializes the point coordinates
void
ActuatorBulkVG::init_points(const ActuatorMetaVG& actMeta)
{
  pointCentroid_.modify_host();

  const int nVGs = actMeta.n_vgs_;
  for (int iVG = 0; iVG < nVGs; iVG++) {
    if (NaluEnv::self().parallel_rank() == assignedProc_.h_view(iVG)) {
      const size_t numForcePts = actMeta.num_force_pts_;
      const int offset = turbIdOffset_.h_view(iVG);
      const double denom = (double)numForcePts;

      // set every pointCentroid
      for (int np = 0; np < numForcePts; np++) {
        auto pointLocal =
          Kokkos::subview(pointCentroid_.view_host(), np + offset, Kokkos::ALL);
        auto centersLocal =
            Kokkos::subview(actMeta.centers_.view_host(), np + offset, Kokkos::ALL);

        for (int i = 0; i < 3; i++)
            pointLocal(i) = centersLocal(i);

        if (actMeta.debug_output_)
          NaluEnv::self().naluOutput()
            << "VG " << iVG // LCCOUT
            << " pointId: " << np << std::scientific << std::setprecision(5)
            << " point: " << pointLocal(0) << " " << pointLocal(1) << " "
            << pointLocal(2) << std::endl;

      } // loop over np
    }
  } // loop over iVG
  actuator_utils::reduce_view_on_host(pointCentroid_.view_host());
  pointCentroid_.sync_host();
}

Kokkos::RangePolicy<ActuatorFixedExecutionSpace>
ActuatorBulkVG::local_range_policy()
{
  auto rank = NaluEnv::self().parallel_rank();
  if (rank < num_vgs_) {
    const int offset = turbIdOffset_.h_view(rank);
    const int size = num_force_pts_;
    return Kokkos::RangePolicy<ActuatorFixedExecutionSpace>(
      offset, offset + size);
  } else {
    return Kokkos::RangePolicy<ActuatorFixedExecutionSpace>(0, 0);
  }
}

void
ActuatorBulkVG::zero_actuator_views()
{
  dvHelper_.touch_dual_view(actuatorForce_);
  dvHelper_.touch_dual_view(velocity_);
  dvHelper_.touch_dual_view(density_);
  Kokkos::deep_copy(dvHelper_.get_local_view(actuatorForce_), 0.0);
  Kokkos::deep_copy(dvHelper_.get_local_view(velocity_), 0.0);
  Kokkos::deep_copy(dvHelper_.get_local_view(density_), 0.0);

#ifdef ENABLE_ACTSIMPLE_PTMOTION
  dvHelper_.touch_dual_view(pointCentroid_);
  Kokkos::deep_copy(dvHelper_.get_local_view(pointCentroid_), 0.0);
#endif
}

} // namespace nalu
} // namespace sierra
