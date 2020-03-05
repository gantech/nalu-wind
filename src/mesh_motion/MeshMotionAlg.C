#include "mesh_motion/MeshMotionAlg.h"

#include "mesh_motion/FrameInertial.h"
#include "mesh_motion/FrameNonInertial.h"
#include "mesh_motion/FrameOpenFAST.h"

#include "NaluParsing.h"

#include <cassert>
#include <iostream>

namespace sierra{
namespace nalu{

MeshMotionAlg::MeshMotionAlg(
  stk::mesh::BulkData& bulk,
  const YAML::Node& node,
  OpenfastFSI* openfast)
{
  load(bulk, node, openfast);
}

void MeshMotionAlg::load(
  stk::mesh::BulkData& bulk,
  const YAML::Node& node,
  OpenfastFSI* openfast)
{
  // get motion information for entire mesh
  const int num_groups = node.size();
  frameVec_.resize(num_groups);

  std::cout << "MeshMotionAlg: Num groups = " << num_groups << std::endl ;

  // temporary vector to store frame names
  std::vector<std::string> frameNames(num_groups);

  for (int i=0; i < num_groups; i++) {

    // extract current motion group info
    const auto& ginfo = node[i];

    // get name of motion group
    frameNames[i] = ginfo["name"].as<std::string>();

    // get frame definition of motion group
    std::string frame;
    get_required(ginfo, "frame", frame);

    if( frame == "inertial" )
      frameVec_[i].reset(new FrameInertial(bulk, ginfo));
    else if( frame == "non_inertial" )
      frameVec_[i].reset(new FrameNonInertial(bulk, ginfo));
    else
      throw std::runtime_error("MeshMotion: Invalid frame type: " + frame);

    // get the reference frame index if it exists
    if(ginfo["reference"])
    {
      std::string refFrameName = ginfo["reference"].as<std::string>();

      auto it = std::find(frameNames.begin(), frameNames.end(), refFrameName);

      if( it ==  frameNames.end() )
        throw std::runtime_error("MeshMotion: Invalid reference frame: " + refFrameName);

      refFrameMap_[i] = frameVec_[std::distance(frameNames.begin(), it)];
    }
  }

  if (openfast != NULL) {
      int nTurbinesGlob = openfast->get_nTurbinesGlob();
      frameVec_.resize(num_groups + nTurbinesGlob);

      for (auto iTurb=0; iTurb < nTurbinesGlob; iTurb++) {
          fsiTurbine *fsiTurbineData = openfast->get_fsiTurbineData(iTurb);
          YAML::Node node; //Empty node
          if (fsiTurbineData != NULL) { //Could be a turbine handled through actuator line or something
              frameVec_[num_groups+iTurb].reset(new FrameOpenFAST(bulk, node, fsiTurbineData));
          } else {
              frameVec_[num_groups+iTurb].reset(new FrameInertial(bulk, node));
          }
      }
  }

}

void MeshMotionAlg::initialize( const double time )
{
  if(isInit_)
    throw std::runtime_error("MeshMotionAlg::initialize(): Re-initialization of MeshMotionAlg not valid");

  std::cout << "FrameVec size = " << frameVec_.size() << std::endl ;
  for (size_t i=0; i < frameVec_.size(); i++)
  {
    frameVec_[i]->setup();

    // set reference frame if they exist
    if( refFrameMap_.find(i) != refFrameMap_.end() )
    {
      MotionBase::TransMatType ref_frame = refFrameMap_[i]->get_inertial_frame();
      frameVec_[i]->set_ref_frame(ref_frame);
    }

    // update coordinates and velocity
    frameVec_[i]->update_coordinates_velocity(time);
  }

  isInit_ = true;
}

void MeshMotionAlg::execute(const double time)
{
  for (size_t i=0; i < frameVec_.size(); i++) {

    if( !frameVec_[i]->is_inertial() )
      frameVec_[i]->update_coordinates_velocity(time);
  }
}

void MeshMotionAlg::post_compute_geometry()
{
  for (size_t i=0; i < frameVec_.size(); i++)
    frameVec_[i]->post_compute_geometry();
}

} // nalu
} // sierra
