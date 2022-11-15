// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#ifndef ACTUATORBULKVG_H_
#define ACTUATORBULKVG_H_

#include <aero/actuator/ActuatorBulk.h>

namespace sierra {
namespace nalu {

struct ActuatorMetaVG : public ActuatorMeta
{
  ActuatorMetaVG(const int num_force_pts, const ActuatorMeta& actMeta);
  virtual ~ActuatorMetaVG() {}

  // HOST ONLY
  bool isotropicGaussian_;

  // Stuff for the simple blade
  bool debug_output_;
  bool useSpreadActuatorForce;
  std::size_t num_force_pts_;
  std::size_t n_vgs_;
  double Cvg_;            // Constant

  ActScalarDblDv areas_;         // Areas of actuator points
  ActVectorDblDv centers_;       // Coordinates of actuator points

  ActVectorDblDv bvec_;         // Vector pointing upwards for the VG
  ActVectorDblDv tvec_;         // Vector pointing along the bottom
                                  // for the VG
  ActVectorDblDv nvec_;         // Vector pointing perpendicular to the plane
                                  // for the VG

  // Kokkos quantities
  // for the blade definitions
  std::size_t max_num_force_pts_;

  std::vector<std::string> output_filenames_;
  bool has_output_file_;
};

struct ActuatorBulkVG : public ActuatorBulk
{
  ActuatorBulkVG(const ActuatorMetaVG& actMeta);
  virtual ~ActuatorBulkVG() {}

  Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace> local_range_policy();

  void init_epsilon(const ActuatorMetaVG& actMeta);
  void init_points(const ActuatorMetaVG& actMeta);
  void init_orientation(const ActuatorMetaVG& actMeta);
  void add_output_headers(const ActuatorMetaVG& actMeta);
  virtual void zero_actuator_views();

  ActScalarDblDv density_;
  ActScalarDblDv alpha_;

  ActFixVectorDbl vg_force_;

  // Stuff for VGs
  const size_t num_force_pts_;
  ActScalarIntDv assignedProc_;
  const int num_vgs_;
  const bool debug_output_;

  ActDualViewHelper<ActuatorMemSpace> dvHelper_;
  std::vector<std::string> output_cache_;
};

} // namespace nalu
} // namespace sierra

#endif /* ACTUATORBULKVG_H_ */
