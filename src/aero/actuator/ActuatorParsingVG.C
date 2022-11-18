// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//
#include <aero/actuator/ActuatorBulk.h>
#include <aero/actuator/ActuatorBulkVG.h>
#include <NaluParsing.h>
#include <aero/actuator/ActuatorParsingVG.h>
#include <aero/actuator/ActuatorParsing.h>
#include <NaluEnv.h>

namespace sierra {
namespace nalu {

namespace {

/** Resizes a vector to size N
 *
 * - If vec is of size 1, resizes a vector to size N (assuming all
 *   constant values).
 * - If vec is of size N, do nothing.
 *
 * Similar to std::vector::resize(), except this returns helpful
 * error message if vec is not size 1 or N.  Useful when validating
 * user input.
 *
 */
std::vector<double>
extend_double_vector(std::vector<double> vec, const unsigned N)
{
  if ((vec.size() != 1) && (vec.size() != N))
    throw std::runtime_error("Vector is not of size 1 or " + std::to_string(N));
  if (vec.size() == 1) { // Extend the vector to size N
    std::vector<double> newvec(N, vec[0]);
    return newvec;
  }
  if (vec.size() == N)
    return vec;
  return vec; // Should not get here
}
} // namespace

ActuatorMetaVG
actuator_VG_parse(const YAML::Node& y_node, const ActuatorMeta& actMeta)
{

  NaluEnv::self().naluOutputP0()
    << "In actuator_VG_parse() " << std::endl; // LCCOUT

  const YAML::Node y_actuator = y_node["actuator"];
  ThrowErrorMsgIf(
    !y_actuator, "actuator argument is "
                 "missing from yaml node passed to actuator_VG_parse");

  size_t num_force_pts_vg;
  get_required(y_actuator, "num_force_pts_vg", num_force_pts_vg);
  ActuatorMetaVG actMetaVG(num_force_pts_vg, actMeta);

  get_if_present(y_actuator, "c_vg", actMetaVG.Cvg_);

  // Load the debug option
  const YAML::Node debug_output = y_actuator["debug_output"];
  if (debug_output)
    actMetaVG.debug_output_ = debug_output.as<bool>();
  else
    actMetaVG.debug_output_ = false;

  size_t n_vgs;
  get_required(y_actuator, "n_vgs", n_vgs);
  actMetaVG.n_vgs_ = n_vgs;
  actMetaVG.numPointsTotal_ = n_vgs * num_force_pts_vg;

  if (actMetaVG.n_vgs_ > 0) {

    for (unsigned iVG = 0; iVG < n_vgs; iVG++) {

      actMetaVG.numPointsTurbine_.h_view(iVG) = num_force_pts_vg;

      const YAML::Node cur_vg = y_actuator["VG" + std::to_string(iVG)];
      get_if_present_no_default(
        cur_vg, "output_file_name", actMetaVG.output_filenames_[iVG]);
      if (
        !actMetaVG.output_filenames_[iVG].empty() &&
        NaluEnv::self().parallel_rank() == (int)iVG) {
        actMetaVG.has_output_file_ = true;
      }

      if (actMetaVG.debug_output_)
        NaluEnv::self().naluOutputP0()
          << "Reading VG: " << iVG << " num_force_pts_vg: "
          << actMetaVG.numPointsTurbine_.h_view(iVG)
          << std::endl; // LCCOUT

      std::string vg_filename;
      get_required(cur_vg, "data_file", vg_filename);
      YAML::Node vg = YAML::LoadFile(vg_filename.c_str());

      //Items to read from VG file

      epsilon_parsing(iVG, cur_vg, actMetaVG);

      // Handle VG properties
      std::vector<double> vec_tmp(3);
      get_required(vg, "bvec", vec_tmp);
      for (size_t j = 0; j < 3; j++)
        actMetaVG.bvec_.h_view(iVG, j) = vec_tmp[j];

      get_required(vg, "tvec", vec_tmp);
      for (size_t j = 0; j < 3; j++)
          actMetaVG.tvec_.h_view(iVG, j) = vec_tmp[j];

      get_required(vg, "nvec", vec_tmp);
      for (size_t j = 0; j < 3; j++)
          actMetaVG.nvec_.h_view(iVG, j) = vec_tmp[j];

      YAML::Node centers = vg["centers"];
      YAML::Node areas = vg["areas"];
      for (size_t j = 0; j < num_force_pts_vg; j++) {
          vec_tmp = centers[j].as<std::vector<double>>();
          for (size_t k = 0; k < 3; k++)
              actMetaVG.centers_.h_view(iVG*num_force_pts_vg + j, k) = vec_tmp[k];
          actMetaVG.areas_.h_view(iVG*num_force_pts_vg + j) =
              areas[j].as<double>();
      }

      // output directions
      if (actMetaVG.debug_output_)
        NaluEnv::self().naluOutputP0() << "VG: " << iVG << std::endl;

    } // End loop over VGs
  } else {
    throw std::runtime_error("Number of VGs <= 0 ");
  }

  if (actMetaVG.debug_output_) {
    NaluEnv::self().naluOutputP0()
      << " actMetaVG.numPointsTotal_ = " << actMetaVG.numPointsTotal_
      << std::endl; // LCCOUT
    NaluEnv::self().naluOutputP0()
      << "Done actuator_VG_parse()" << std::endl; // LCCOUT
  }
  return actMetaVG;
}

} // namespace nalu
} // namespace sierra
