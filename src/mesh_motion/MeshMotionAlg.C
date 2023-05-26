#include "mesh_motion/MeshMotionAlg.h"

#include "mesh_motion/FrameMoving.h"

#include "NaluParsing.h"

#include <cassert>
#include <iostream>

namespace sierra {
namespace nalu {

MeshMotionAlg::MeshMotionAlg(stk::mesh::BulkData& bulk, const YAML::Node& node)
{
  load(bulk, node);

  set_deformation_flag();
}

void
MeshMotionAlg::load(stk::mesh::BulkData& bulk, const YAML::Node& node)
{
  // get motion information for entire mesh
  const int num_groups = node.size();

  int num_groups_smd = 0;
  int num_groups_mesh_motion = 0;
  for (int i = 0; i < num_groups; i++) {
    // extract current motion group info
    const auto& ginfo = node[i];
    bool enable_smd=false;
    get_if_present(ginfo, "enable_smd", enable_smd, enable_smd);
    std::cerr << "Enable SMD = " << enable_smd << std::endl;
    if ( enable_smd )
      num_groups_smd += 1;
    else
      num_groups_mesh_motion += 1;
  }

  movingFrameVec_.resize(num_groups_mesh_motion);
  smdFrameVec_.resize(num_groups_smd);

  int i_mm = 0;
  int i_smd = 0;
  for (int i = 0; i < num_groups; i++) {
    // extract current motion group info
    const auto& ginfo = node[i];
    bool enable_smd=false;
    get_if_present(ginfo, "enable_smd", enable_smd, enable_smd);
    if (enable_smd) {
      smdFrameVec_[i_smd].reset(new FrameSMD(bulk, ginfo));
      i_smd++;
    }
    else {
      movingFrameVec_[i_mm].reset(new FrameMoving(bulk, ginfo));
      i_mm++;
    }
  }
}

void
MeshMotionAlg::set_deformation_flag()
{
  for (size_t i = 0; i < movingFrameVec_.size(); i++)
    if (movingFrameVec_[i]->is_deforming())
      isDeforming_ = true;

  for (size_t i = 0; i < smdFrameVec_.size(); i++)
      if (smdFrameVec_[i]->is_deforming())
          isDeforming_ = true;
}

void
MeshMotionAlg::initialize(const double time)
{
  if (isInit_)
    throw std::runtime_error("MeshMotionAlg::initialize(): Re-initialization "
                             "of MeshMotionAlg not valid");

  for (size_t i = 0; i < movingFrameVec_.size(); i++) {
    movingFrameVec_[i]->setup();
    smdFrameVec_[i]->setup();

    // update coordinates and velocity
    movingFrameVec_[i]->update_coordinates_velocity(time);
    smdFrameVec_[i]->update_coordinates_velocity(time);
  }

  isInit_ = true;
}

void
MeshMotionAlg::restart_reinit(const double time)
{
  if (isInit_) {
    isInit_ = false;
    initialize(time);
  } else {
    throw std::runtime_error(
      "MeshMotionAlg::restart_reinit(): Re-initialization of MeshMotionAlg for "
      "restart should be called after initialize");
  }
}

void
MeshMotionAlg::execute(const double time)
{
  for (size_t i = 0; i < movingFrameVec_.size(); i++)
    movingFrameVec_[i]->update_coordinates_velocity(time);
  for (size_t i = 0; i < smdFrameVec_.size(); i++)
    smdFrameVec_[i]->update_coordinates_velocity(time);

}

void
MeshMotionAlg::post_compute_geometry()
{
  for (size_t i = 0; i < movingFrameVec_.size(); i++)
    movingFrameVec_[i]->post_compute_geometry();
  for (size_t i = 0; i < smdFrameVec_.size(); i++)
    smdFrameVec_[i]->post_compute_geometry();
}

stk::mesh::PartVector
MeshMotionAlg::get_partvec()
{
  stk::mesh::PartVector fpartVec;
  for (size_t i = 0; i < movingFrameVec_.size(); i++) {
    stk::mesh::PartVector fPartVec = movingFrameVec_[i]->get_partvec();
    for (auto p : fPartVec)
      fpartVec.push_back(p);
  }
  return fpartVec;
}

void
MeshMotionAlg::predict_states_smd()
{
  for (size_t i = 0; i < smdFrameVec_.size(); i++)
    smdFrameVec_[i]->predict_states();
}

void
MeshMotionAlg::update_timestep_smd()
{
  for (size_t i = 0; i < smdFrameVec_.size(); i++)
    smdFrameVec_[i]->update_timestep();
}

void
MeshMotionAlg::advance_timestep_smd()
{
  for (size_t i = 0; i < smdFrameVec_.size(); i++)
    smdFrameVec_[i]->advance_timestep();
}
    
} // namespace nalu
} // namespace sierra
