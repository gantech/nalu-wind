#ifndef MESHMOTIONALG_H
#define MESHMOTIONALG_H

#include "FrameMoving.h"
#include "FrameSMD.h"

namespace sierra {
namespace nalu {

class MeshMotionAlg
{
public:
  // TODO fsi data needs to be supplied to mesh motion/or triggered
  MeshMotionAlg(std::shared_ptr<stk::mesh::BulkData> bulk, const YAML::Node&);

  ~MeshMotionAlg() {}

  void initialize(const double, std::shared_ptr<stk::mesh::BulkData> bulk);

  void restart_reinit(const double, std::shared_ptr<stk::mesh::BulkData> bulk);

  void execute(const double);

  void post_compute_geometry();

  void predict_states_smd();

  void update_timestep_smd();

  void advance_timestep_smd();

  stk::mesh::PartVector get_partvec();

  bool is_deforming() { return isDeforming_; }

  bool is_smd() { return isSMD_; }

private:
  MeshMotionAlg() = delete;
  MeshMotionAlg(const MeshMotionAlg&) = delete;

  void load(std::shared_ptr<stk::mesh::BulkData>, const YAML::Node&);

  void set_deformation_flag();

  /** Moving frame vector
   *
   *  Vector of moving frames
   *  Size is the number of groups under mesh_motion in input file
   */
  std::vector<std::shared_ptr<FrameMoving>> movingFrameVec_;

  /** Spring-Mass-Damper frame vector
   *
   *  Vector of frames connected to a spring-mass-damper
   *  Size is the number of groups under mesh_motion in input file
   */
  std::vector<std::shared_ptr<FrameSMD>> smdFrameVec_;

  //! flag to guard against multiple invocations of initialize()
  bool isInit_ = false;

  //! flag to denote if mesh deformation exists
  bool isDeforming_ = false;

  //! flag to denote if any SMD frames are active
  bool isSMD_ = false;
};

} // namespace nalu
} // namespace sierra

#endif /* MESHMOTIONALG_H */
