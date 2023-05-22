#ifndef FRAMESMD_H
#define FRAMESMD_H

#include "NgpMotion.h"
#include "FrameBase.h"

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
  FrameSMD(stk::mesh::BulkData&, const YAML::Node&);

  virtual ~FrameSMD();

  virtual void setup();

  stk::mesh::PartVector get_partvec() { return partVec_; };

  bool is_deforming() { return isDeforming_; }

  virtual void update_coordinates_velocity(const double time);
    
  void post_compute_geometry();
   
private:
  FrameSMD() = delete;
  FrameSMD(const FrameSMD&) = delete;

  void load(const YAML::Node&);
    
};

} // namespace nalu
} // namespace sierra

#endif /* FRAMESMD_H */
