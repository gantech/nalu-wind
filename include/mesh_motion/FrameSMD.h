#ifndef FRAMESMD_H
#define FRAMESMD_H

#include "NgpMotion.h"
#include "FrameBase.h"

// stk base header files
#include "stk_mesh/base/CoordinateSystems.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/MetaData.hpp"

#include "mesh_motion/SMD.h"

namespace YAML {
class Node;
}

namespace sierra {
namespace nalu {

class FrameSMD : public FrameBase
{
public:
  FrameSMD(std::shared_ptr<stk::mesh::BulkData>, const YAML::Node&);

  virtual ~FrameSMD();

  virtual void setup(std::shared_ptr<stk::mesh::BulkData> bulk);

  stk::mesh::PartVector get_partvec() { return partVec_; };

  bool is_deforming() { return isDeforming_; }

  virtual void update_coordinates_velocity(const double time);

  void post_compute_geometry();

  void predict_states();

  void update_timestep();

  void advance_timestep();

private:
  FrameSMD() = delete;
  FrameSMD(const FrameSMD&) = delete;

  void load(const YAML::Node&);

  /** Spring-Mass-Damper vector
   *
   *  A vector of size number of motion groups
   */
  std::vector<std::unique_ptr<SMD>> smd_;

};

} // namespace nalu
} // namespace sierra

#endif /* FRAMESMD_H */
