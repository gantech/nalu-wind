#ifndef MOTIONOSCILLATIONCYLINDERRBF_H
#define MOTIONOSCILLATIONCYLINDERRBF_H

#include "NgpMotion.h"
#include <Kokkos_Core.hpp>

namespace sierra {
namespace kynema_ugf {

class MotionOscillationCylinderRBF : public NgpMotionKernel<MotionOscillationCylinderRBF>
{
public:
  MotionOscillationCylinderRBF(const YAML::Node&);

  MotionOscillationCylinderRBF() = default;

  virtual ~MotionOscillationCylinderRBF() = default;

  /** Function to compute motion-specific transformation matrix
   *
   * @param[in] time Current time
   * @param[in] xyz  Coordinates
   * @return Transformation matrix
   */
  KOKKOS_FUNCTION
  virtual mm::TransMatType
  build_transformation(const double& time, const mm::ThreeDVecType& xyz);

  /** Function to compute motion-specific velocity
   *
   * @param[in]  time      Current time
   * @param[in]  compTrans Transformation matrix
   *                       including all motions
   * @param[in]  mxyz      Model coordinates
   * @param[in]  cxyz      Transformed coordinates
   * @return Velocity vector associated with coordinates
   */
  KOKKOS_FUNCTION
  virtual mm::ThreeDVecType compute_velocity(
    const double& time,
    const mm::TransMatType& compTrans,
    const mm::ThreeDVecType& mxyz,
    const mm::ThreeDVecType& cxyz);

private:
  void load(const YAML::Node&);

  // Cylinder parameters
  double cylinderRadius_{1.0};
  double cylinderHeight_{4.0};
  double amplitudeFactor_{0.3};
  double frequency_{1.0};
  double rigidRadius_{1.0};
  double outerBlendRadius_{8.0};

  // Domain extents (min and max for x, y, z)
  double domainExtents_[6]{-8.0, -8.0, 0.0, 8.0, 8.0, 4.0};
};

} // namespace kynema_ugf
} // namespace sierra

#endif /* MOTIONOSCILLATIONCYLINDERRBF_H */
