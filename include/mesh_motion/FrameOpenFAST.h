#ifndef FRAMEOPENFAST_H
#define FRAMEOPENFAST_H

#include "FrameBase.h"
#include "FSIturbine.h"

#include "yaml-cpp/yaml.h"

#include <cassert>
#include <float.h>

namespace sierra{
namespace nalu{

class FrameOpenFAST : public FrameBase
{
public:
  FrameOpenFAST(
    stk::mesh::BulkData& bulk,
    YAML::Node node,
    fsiTurbine* fsiturbinedata
  ) : FrameBase(bulk,node,false),
      fsiTurbineData_(fsiturbinedata)
  {}

  virtual ~FrameOpenFAST() {}

  //Redefine setup to empty function
  void setup() {};

  void update_coordinates_velocity(const double time);

private:
    FrameOpenFAST() = delete;
    FrameOpenFAST(const FrameOpenFAST&) = delete;

    fsiTurbine* fsiTurbineData_;
    /** Compute transformation matrix
     *
     * @return 4x4 matrix representing composite addition of motions
     */
    MotionBase::TransMatType compute_transformation(
      const double,
      const double*);
};

} // nalu
} // sierra

#endif /* FRAMEOPENFAST_H */
