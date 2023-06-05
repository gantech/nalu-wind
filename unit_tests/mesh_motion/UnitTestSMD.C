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
                              "mass_matrix: [7000.0, -90.0, 60.0, -90.0, 6400.0, 1800.0, 60.0, 1800.0, 7800] \n"
                              "stiffness_matrix: [131000.0, 7700.0, -97000.0, 7700.0, 66000.0, -31000.0, -97000.0, -31000.0, 4800000.0] \n"
                              "damping_matrix: [300.0, 8.0, -70.0, 8.0, 200.0, 20.0, -70.0, 20.0, 3700.0] \n"
                              "x_init: [0.1, 0.2, 0.01] \n"
                              "loads_scale: 0.75 \n"
                              "centroid: [0.3,0.5,0.0] \n"
                              "alpha : -0.05 \n";

  YAML::Node smd_node = YAML::Load(smd_info);

  // initialize the mesh rotation class
  sierra::nalu::AirfoilSMD smd_class(smd_node);
  smd_class.setup(0.01); // Set timestep here
  
  // Make appropriate forces and moments to test here
  vs::Vector f_np1;
  vs::Vector m_np1;

  // Meaningless values are passed for force z and moments x,y
  f_np1.x() = 25000.0;   f_np1.y() = 35000.0;   f_np1.z() = 7000000.0;
  m_np1.x() = 7000000.0; m_np1.y() = 7000000.0; m_np1.z() = 50000.0;

  smd_class.update_timestep(f_np1, m_np1);

  vs::Vector trans_disp = smd_class.get_trans_disp();
  double rot_disp = smd_class.get_rot_disp();
  vs::Vector trans_vel = smd_class.get_trans_vel();
  vs::Vector rot_vel = smd_class.get_rot_vel();

  // appropriate expected values of the states at t=n+1
  const double gold_x0 = 0.10004032391557534;
  const double gold_x1 = 0.20007471912761296;
  const double gold_x2 = 0.010036213720668577;

  const double gold_xdot0 = 0.00804649562500858;
  const double gold_xdot1 = 0.014909939296912205;
  const double gold_xdot2 = 0.007226320677630027;
  
  EXPECT_NEAR(trans_disp[0], gold_x0, testTol);
  EXPECT_NEAR(trans_disp[1], gold_x1, testTol);
  EXPECT_NEAR(trans_disp[2], 0.0, testTol);
  EXPECT_NEAR(rot_disp, gold_x2, testTol);
  EXPECT_NEAR(trans_vel[0], gold_xdot0, testTol);
  EXPECT_NEAR(trans_vel[1], gold_xdot1, testTol);
  EXPECT_NEAR(trans_vel[2], 0.0, testTol);
  EXPECT_NEAR(rot_vel[0], 0.0, testTol);
  EXPECT_NEAR(rot_vel[1], 0.0, testTol);
  EXPECT_NEAR(rot_vel[2], gold_xdot2, testTol);

}


TEST(meshMotion, airfoil_smd_update2)
{
  // create a yaml node describing Airfoil SMD
  const std::string smd_info = "axis: [0.0, 0.0, -1.0] \n"
                              "centroid: [0.25, 0.0, 0.0] \n"
                              "mass_matrix: [7000.0, -90.0, 60.0, -90.0, 6400.0, 1800.0, 60.0, 1800.0, 7800] \n"
                              "stiffness_matrix: [1310.0, 77.0, -970.0, 77.0, 660.0, -310.0, -970.0, -310.0, 48000.0] \n"
                              "damping_matrix: [300.0, 8.0, -70.0, 8.0, 200.0, 20.0, -70.0, 20.0, 3700.0] \n"
                              "force_transform_matrix: [62.0, 0.0, 10.0, 0.0, 61.0, 0.0,  1.0, -6.0, 60.0]\n"
                              "x_init: [0.1, 0.2, 0.01] \n"
                              "xdot_init: [0.2, 0.3, -0.05] \n"
                              "a_init: [0.3, 0.2, 0.0] \n"
                              "f_init: [200.0, 400.0, 50.0] \n"
                              "loads_scale: 0.75 \n"
                              "centroid: [0.3,0.5,0.0] \n"
                              "alpha : -0.05 \n";


  YAML::Node smd_node = YAML::Load(smd_info);

  // initialize the mesh rotation class
  sierra::nalu::AirfoilSMD smd_class(smd_node);
  smd_class.setup(0.01); // Set timestep here
  
  // Make appropriate forces and moments to test here
  vs::Vector f_np1;
  vs::Vector m_np1;

  // Meaningless values are passed for force z and moments x,y
  f_np1.x() = 250.0;   f_np1.y() = 350.0;   f_np1.z() = 70000.0;
  m_np1.x() = 70000.0; m_np1.y() = 70000.0; m_np1.z() = 500.0;

  smd_class.update_timestep(f_np1, m_np1);

  vs::Vector trans_disp = smd_class.get_trans_disp();
  double rot_disp = smd_class.get_rot_disp();
  vs::Vector trans_vel = smd_class.get_trans_vel();
  vs::Vector rot_vel = smd_class.get_rot_vel();

  // appropriate expected values of the states at t=n+1
  const double gold_x0 = 0.10208538238772817;
  const double gold_x1 = 0.20307565383807263;
  const double gold_x2 = 0.009577049930923145;

  const double gold_xdot0 = 0.2170445580954166;
  const double gold_xdot1 = 0.31510099263127095;
  const double gold_xdot2 = -0.03462495709470137;
  
  EXPECT_NEAR(trans_disp[0], gold_x0, testTol);
  EXPECT_NEAR(trans_disp[1], gold_x1, testTol);
  EXPECT_NEAR(trans_disp[2], 0.0, testTol);
  EXPECT_NEAR(rot_disp, gold_x2, testTol);
  EXPECT_NEAR(trans_vel[0], gold_xdot0, testTol);
  EXPECT_NEAR(trans_vel[1], gold_xdot1, testTol);
  EXPECT_NEAR(trans_vel[2], 0.0, testTol);
  EXPECT_NEAR(rot_vel[0], 0.0, testTol);
  EXPECT_NEAR(rot_vel[1], 0.0, testTol);
  EXPECT_NEAR(rot_vel[2], gold_xdot2, testTol);

}

TEST(meshMotion, airfoil_smd_predict)
{
  // create a yaml node describing Airfoil SMD
  const std::string smd_info = "axis: [0.0, 0.0, -1.0] \n"
                              "centroid: [0.25, 0.0, 0.0] \n"
                              "mass_matrix: [7000.0, -90.0, 60.0, -90.0, 6400.0, 1800.0, 60.0, 1800.0, 7800] \n"
                              "stiffness_matrix: [1310.0, 77.0, -970.0, 77.0, 660.0, -310.0, -970.0, -310.0, 48000.0] \n"
                              "damping_matrix: [300.0, 8.0, -70.0, 8.0, 200.0, 20.0, -70.0, 20.0, 3700.0] \n"
                              "force_transform_matrix: [62.0, 0.0, 10.0, 0.0, 61.0, 0.0,  1.0, -6.0, 60.0]\n"
                              "x_init: [0.1, 0.2, 0.01] \n"
                              "xdot_init: [1.0, -2.0, 0.05] \n"
                              "a_init: [-0.5, -7.0, 1.0] \n"
                              "x_nm1: [0.0, 0.0, 0.0] \n"
                              "xdot_nm1: [1.1, -1.4, -0.05] \n"
                              "a_nm1: [0.5, -8.0, 2.0] \n"
                              "f_init: [200.0, 400.0, 50.0] \n"
                              "loads_scale: 0.75 \n"
                              "centroid: [0.3,0.5,0.0] \n"
                              "alpha : -0.05 \n";


  YAML::Node smd_node = YAML::Load(smd_info);

  // initialize the mesh rotation class
  sierra::nalu::AirfoilSMD smd_class(smd_node);
  smd_class.setup(0.01); // Set timestep here
  
  // Predict next state
  smd_class.predict_states();

  vs::Vector trans_disp = smd_class.get_trans_disp();
  double rot_disp = smd_class.get_rot_disp();
  vs::Vector trans_vel = smd_class.get_trans_vel();
  vs::Vector rot_vel = smd_class.get_rot_vel();

  // appropriate expected values of the states at t=n+1
  const double gold_x0 = 0.1095;
  const double gold_x1 = 0.177;
  const double gold_x2 = 0.011;

  const double gold_xdot0 = 0.99;
  const double gold_xdot1 = -2.065;
  const double gold_xdot2 = 0.055;
  
  EXPECT_NEAR(trans_disp[0], gold_x0, testTol);
  EXPECT_NEAR(trans_disp[1], gold_x1, testTol);
  EXPECT_NEAR(trans_disp[2], 0.0, testTol);
  EXPECT_NEAR(rot_disp, gold_x2, testTol);
  EXPECT_NEAR(trans_vel[0], gold_xdot0, testTol);
  EXPECT_NEAR(trans_vel[1], gold_xdot1, testTol);
  EXPECT_NEAR(trans_vel[2], 0.0, testTol);
  EXPECT_NEAR(rot_vel[0], 0.0, testTol);
  EXPECT_NEAR(rot_vel[1], 0.0, testTol);
  EXPECT_NEAR(rot_vel[2], gold_xdot2, testTol);

}

TEST(meshMotion, airfoil_smd_advance)
{
  // create a yaml node describing Airfoil SMD
  const std::string smd_info = "axis: [0.0, 0.0, -1.0] \n"
                              "centroid: [0.25, 0.0, 0.0] \n"
                              "mass_matrix: [7000.0, -90.0, 60.0, -90.0, 6400.0, 1800.0, 60.0, 1800.0, 7800] \n"
                              "stiffness_matrix: [1310.0, 77.0, -970.0, 77.0, 660.0, -310.0, -970.0, -310.0, 48000.0] \n"
                              "damping_matrix: [300.0, 8.0, -70.0, 8.0, 200.0, 20.0, -70.0, 20.0, 3700.0] \n"
                              "force_transform_matrix: [62.0, 0.0, 10.0, 0.0, 61.0, 0.0,  1.0, -6.0, 60.0]\n"
                              "x_init: [0.1, 0.2, 0.01] \n"
                              "xdot_init: [0.2, 0.3, -0.05] \n"
                              "a_init: [0.3, 0.2, 0.0] \n"
                              "f_init: [200.0, 400.0, 50.0] \n"
                              "loads_scale: 0.75 \n"
                              "centroid: [0.3,0.5,0.0] \n"
                              "alpha : -0.05 \n";

  YAML::Node smd_node = YAML::Load(smd_info);

  // initialize the mesh rotation class
  sierra::nalu::AirfoilSMD smd_class(smd_node);
  smd_class.setup(0.01); // Set timestep here
  
  // Check two updates before/after advance to make sure advance 
  // does not copy pointers in a way that breaks test

  // Forces for first update
  vs::Vector f_np1;
  vs::Vector m_np1;

  // Forces for second update
  vs::Vector f_np2;
  vs::Vector m_np2;


  // Update
  f_np1.x() = 250.0;     f_np1.y() = 350.0;     f_np1.z() = 7000000.0;
  m_np1.x() = 7000000.0; m_np1.y() = 7000000.0; m_np1.z() = 500.0;

  smd_class.update_timestep(f_np1, m_np1);

  // Advance
  smd_class.advance_timestep();

  // Update
  f_np2.x() = 350.0;     f_np2.y() = 450.0;     f_np2.z() = 7000000.0;
  m_np2.x() = 7000000.0; m_np2.y() = 7000000.0; m_np2.z() = 700.0;

  smd_class.update_timestep(f_np2, m_np2);

  // Check final
  vs::Vector trans_disp = smd_class.get_trans_disp();
  double rot_disp = smd_class.get_rot_disp();
  vs::Vector trans_vel = smd_class.get_trans_vel();
  vs::Vector rot_vel = smd_class.get_rot_vel();

  // appropriate expected values of the states at t=n+1
  const double gold_x0 = 0.10443046438429568;
  const double gold_x1 = 0.20636932163370228;
  const double gold_x2 = 0.009410229273121231;

  const double gold_xdot0 = 0.25195734751525295;
  const double gold_xdot1 = 0.34362641777071373;
  const double gold_xdot2 = 0.0012428410464164694;
  
  EXPECT_NEAR(trans_disp[0], gold_x0, testTol);
  EXPECT_NEAR(trans_disp[1], gold_x1, testTol);
  EXPECT_NEAR(trans_disp[2], 0.0, testTol);
  EXPECT_NEAR(rot_disp, gold_x2, testTol);
  EXPECT_NEAR(trans_vel[0], gold_xdot0, testTol);
  EXPECT_NEAR(trans_vel[1], gold_xdot1, testTol);
  EXPECT_NEAR(trans_vel[2], 0.0, testTol);
  EXPECT_NEAR(rot_vel[0], 0.0, testTol);
  EXPECT_NEAR(rot_vel[1], 0.0, testTol);
  EXPECT_NEAR(rot_vel[2], gold_xdot2, testTol);

}

