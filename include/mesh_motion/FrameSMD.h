#ifndef FRAMESMD_H
#define FRAMESMD_H

#include "NgpMotion.h"
#include "FrameBase.h"
#include "mesh_motion/SMD.h"
#include "aero/aero_utils/CalcLoads.h"

// stk base header files
#include "stk_mesh/base/CoordinateSystems.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/MetaData.hpp"


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

  virtual void setup(std::shared_ptr<stk::mesh::BulkData> bulk) {}

  void setup(const double dt, std::shared_ptr<stk::mesh::BulkData> bulk);
    
  virtual void initialize();

  stk::mesh::PartVector get_partvec() { return partVec_; };

  bool is_deforming() { return isDeforming_; }

  virtual void update_coordinates_velocity(const double time);

  void post_compute_geometry();

  void predict_states();

  void update_timestep(double cur_time);

  void advance_timestep(const double cur_time);

private:
  FrameSMD() = delete;
  FrameSMD(const FrameSMD&) = delete;

  void load(const YAML::Node&);

  double ramp_function(double position, const double start, const double end);

  /** Spring-Mass-Damper vector
   *
   *  A vector of size number of motion groups
   */
  std::vector<std::unique_ptr<SMD>> smd_;

  // Pointer to Algorithm that calculates loads on the surfaces of the SMD
  std::unique_ptr<CalcLoads> calc_loads_;

  // Scale loads by this factor when transfering to SMD
  double loads_scale_{1.0};

  // Have a transition function to go from rigid body motion to no motion of the mesh
  // Transition starts at ramp_lower_ and ends at ramp_upper_
  // Default values set based on initial SMD Airfoil simulations with chord=1.0
  double mesh_ramp_lower_{30.0};
  double mesh_ramp_upper_{100.0};

  // Default is that there is no ramp in for the load.
  double load_ramp_lower_{-1.0};
  double load_ramp_upper_{0.0};
    
};

} // namespace nalu
} // namespace sierra

#endif /* FRAMESMD_H */
