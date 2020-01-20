// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#include "utils/UnitTestGCL.h"

TEST_F(GCLTest, rigid_rotation_elem)
{
  if (bulk_.parallel_size() > 1) return;

  const std::string meshDims = "3x3x3|offset:0,65,0";
  const bool secondOrder = false;
  const double deltaT = 0.003; // approx 0.25 deg motion for given omega
  const std::string mesh_motion =
    "mesh_motion:                                                          \n"
    "  - name: interior                                                    \n"
    "    frame: non_inertial                                               \n"
    "    mesh_parts: [ block_1 ]                                           \n"
    "    motion:                                                           \n"
    "      - type: rotation                                                \n"
    "        omega: 1.5707963267948966                                     \n"
    "        axis: [1.0, 0.0, 0.0]                                         \n"
    "        centroid: [0.0, 0.0, 0.0]                                     \n";

  fill_mesh_and_init_fields(meshDims);
  init_time_integrator(secondOrder, deltaT);
  register_algorithms(mesh_motion);
  init_states();
  compute_div_mesh_vel();
  compute_dvoldt();
  compute_absolute_error();
}

TEST_F(GCLTest, rigid_rotation_edge)
{
  if (bulk_.parallel_size() > 1) return;

  realm_.realmUsesEdges_ = true;
  const std::string meshDims = "3x3x3|offset:0,65,0";
  const bool secondOrder = false;
  const double deltaT = 0.003; // approx 0.25 deg motion for given omega
  const std::string mesh_motion =
    "mesh_motion:                                                          \n"
    "  - name: interior                                                    \n"
    "    frame: non_inertial                                               \n"
    "    mesh_parts: [ block_1 ]                                           \n"
    "    motion:                                                           \n"
    "      - type: rotation                                                \n"
    "        omega: 1.5707963267948966                                     \n"
    "        axis: [1.0, 0.0, 0.0]                                         \n"
    "        centroid: [0.0, 0.0, 0.0]                                     \n";

  fill_mesh_and_init_fields(meshDims);
  init_time_integrator(secondOrder, deltaT);
  register_algorithms(mesh_motion);
  init_states();
  compute_div_mesh_vel();
  compute_dvoldt();
  compute_absolute_error();
}

TEST_F(GCLTest, rigid_scaling_elem)
{
  if (bulk_.parallel_size() > 1) return;

  const std::string meshDims = "3x3x3|offset:0,65,0";
  const bool secondOrder = false;
  const double deltaT = 0.003; // approx 0.25 deg motion for given omega
  const std::string mesh_motion =
    "mesh_motion:                                                          \n"
    "  - name: interior                                                    \n"
    "    frame: non_inertial                                               \n"
    "    mesh_parts: [ block_1 ]                                           \n"
    "    motion:                                                           \n"
    "      - type: scaling                                                 \n"
    "        rate: [1.0, 1.0, 1.0]                                         \n"
    "        centroid: [0.0, 0.0, 0.0]                                     \n";

  fill_mesh_and_init_fields(meshDims);
  init_time_integrator(secondOrder, deltaT);
  register_algorithms(mesh_motion);
  init_states();
  compute_div_mesh_vel();
  compute_dvoldt();
  compute_relative_error();
}

TEST_F(GCLTest, rigid_scaling_edge)
{
  if (bulk_.parallel_size() > 1) return;

  realm_.realmUsesEdges_ = true;
  const std::string meshDims = "3x3x3|offset:0,65,0";
  const bool secondOrder = false;
  const double deltaT = 0.003; // approx 0.25 deg motion for given omega
  const std::string mesh_motion =
    "mesh_motion:                                                          \n"
    "  - name: interior                                                    \n"
    "    frame: non_inertial                                               \n"
    "    mesh_parts: [ block_1 ]                                           \n"
    "    motion:                                                           \n"
    "      - type: scaling                                                 \n"
    "        rate: [1.0, 1.0, 1.0]                                         \n"
    "        centroid: [0.0, 0.0, 0.0]                                     \n";

  fill_mesh_and_init_fields(meshDims);
  init_time_integrator(secondOrder, deltaT);
  register_algorithms(mesh_motion);
  init_states();
  compute_div_mesh_vel();
  compute_dvoldt();
  compute_relative_error();
}


TEST_F(GCLTest, rigid_translation)
{
  if (bulk_.parallel_size() > 1) return;

  const std::string meshDims = "3x3x3|offset:0,65,0";
  const bool secondOrder = false;
  const double deltaT = 0.003; // approx 0.25 deg motion for given omega
  const std::string mesh_motion =
    "mesh_motion:                                                          \n"
    "  - name: interior                                                    \n"
    "    frame: non_inertial                                               \n"
    "    mesh_parts: [ block_1 ]                                           \n"
    "    motion:                                                           \n"
    "      - type: translation                                             \n"
    "        velocity: [1.0, 0.0, 0.0]                                     \n";

  fill_mesh_and_init_fields(meshDims);
  init_time_integrator(secondOrder, deltaT);
  register_algorithms(mesh_motion);
  init_states();
  compute_div_mesh_vel();
  compute_dvoldt();
  compute_absolute_error();
}
