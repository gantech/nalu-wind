#include <gtest/gtest.h>
#include <limits>

#include "mesh_motion/AirfoilSMD.h"
#include "UnitTestRealm.h"

namespace {

const double testTol = 1e-14;

} // namespace

TEST(meshMotion, airfoil_smd_update)
{
  // create a yaml node describing Airfoil SMD
  const std::string smd_info = "axis: [0.0, 0.0, -1.0] \n"
                              "centroid: [0.25, 0.0, 0.0] \n"
                              "mass_matrix: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0] \n"
                              "stiffness_matrix: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0] \n"
                              "damping_matrix: [0.1, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.1] \n"
                              "x_init: [0.0, 0.0, 2.0] \n"
                              "loads_scale: 1.0 \n"
                              "centroid: [0.3,0.5,0.0] \n";

  YAML::Node smd_node = YAML::Load(smd_info);

  // initialize the mesh rotation class
  sierra::nalu::AirfoilSMD smd_class(smd_node);
  smd_class.setup(0.001); // Set timestep here
  
  // TODO: Make appropriate forces and moments to test here
  vs::Vector f_np1;
  vs::Vector m_np1;
  smd_class.update_timestep(f_np1, m_np1);

  vs::Vector trans_disp = smd_class.get_trans_disp();
  double rot_disp = smd_class.get_rot_disp();
  vs::Vector trans_vel = smd_class.get_trans_vel();
  vs::Vector rot_vel = smd_class.get_rot_vel();

  // TODO: Make appropriate expected values of the states at t=n+1
  const double gold_x0 = 0.0;
  const double gold_x1 = 0.0;
  const double gold_x2 = 0.0;

  const double gold_xdot0 = 0.0;
  const double gold_xdot1 = 0.0;
  const double gold_xdot2 = 0.0;
  
  EXPECT_NEAR(trans_disp[0], gold_x0, testTol);
  EXPECT_NEAR(trans_disp[1], gold_x1, testTol);
  EXPECT_NEAR(trans_disp[2], 0.0, testTol);
  EXPECT_NEAR(rot_disp, gold_x2, testTol);
  EXPECT_NEAR(trans_vel[0], gold_xdot0, testTol);
  EXPECT_NEAR(trans_vel[1], gold_xdot1, testTol);
  EXPECT_NEAR(trans_vel[2], 0.0, testTol);
  EXPECT_NEAR(rot_vel[0], 0.0, testTol);
  EXPECT_NEAR(rot_vel[1], 0.0, testTol);
  EXPECT_NEAR(rot_vel[2], rot_vel[2], testTol);

}

