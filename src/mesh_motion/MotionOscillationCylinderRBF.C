#include "mesh_motion/MotionOscillationCylinderRBF.h"

#include <KynemaUGFEnv.h>
#include <KynemaUGFParsing.h>
#include <cmath>
#include <stdexcept>

namespace sierra {
namespace kynema_ugf {

namespace {

KOKKOS_INLINE_FUNCTION
double
radial_blend_factor(
  const double radius, const double rigidRadius, const double outerBlendRadius)
{
  if (radius <= rigidRadius)
    return 1.0;

  if (radius >= outerBlendRadius)
    return 0.0;

  const double xi =
    (radius - rigidRadius) / (outerBlendRadius - rigidRadius);

  // Smooth monotone cubic blend: 1 at xi=0, 0 at xi=1.
  return 1.0 - 3.0 * xi * xi + 2.0 * xi * xi * xi;
}

}

MotionOscillationCylinderRBF::MotionOscillationCylinderRBF(const YAML::Node& node)
  : NgpMotionKernel<MotionOscillationCylinderRBF>()
{
  load(node);
  isDeforming_ = true;
}

void
MotionOscillationCylinderRBF::load(const YAML::Node& node)
{
  // Parse time bounds
  get_if_present(node, "start_time", startTime_, startTime_);
  startTime_ = startTime_ - DBL_EPSILON;

  get_if_present(node, "end_time", endTime_, endTime_);
  endTime_ = endTime_ + DBL_EPSILON;

  // Parse cylinder parameters
  get_if_present(node, "cylinder_radius", cylinderRadius_, cylinderRadius_);
  get_if_present(node, "cylinder_height", cylinderHeight_, cylinderHeight_);
  get_required(node, "frequency", frequency_);
  get_if_present(node, "amplitude_factor", amplitudeFactor_, amplitudeFactor_);

  // Parse domain extents if provided
  if (node["domain_extents"]) {
    for (int i = 0; i < 6; ++i) {
      domainExtents_[i] = node["domain_extents"][i].as<double>();
    }
  }

  get_if_present(node, "rigid_radius", rigidRadius_, cylinderRadius_);

  const double defaultOuterBlendRadius = stk::math::min(
    stk::math::min(stk::math::abs(domainExtents_[0]), stk::math::abs(domainExtents_[3])),
    stk::math::min(stk::math::abs(domainExtents_[1]), stk::math::abs(domainExtents_[4])));
  outerBlendRadius_ = defaultOuterBlendRadius;
  get_if_present(
    node, "outer_blend_radius", outerBlendRadius_, outerBlendRadius_);

  if (rigidRadius_ < cylinderRadius_) {
    throw std::runtime_error(
      "MotionOscillationCylinderRBF: rigid_radius must be >= cylinder_radius.");
  }

  if (outerBlendRadius_ <= rigidRadius_) {
    throw std::runtime_error(
      "MotionOscillationCylinderRBF: outer_blend_radius must be > rigid_radius.");
  }

  KynemaUGFEnv::self().kynema_ugfOutputP0()
    << "MotionOscillationCylinderRBF initialized with analytic x-y radial blend"
    << ", frequency=" << frequency_ << ", amplitude_factor="
    << amplitudeFactor_ << ", rigid_radius=" << rigidRadius_
    << ", outer_blend_radius=" << outerBlendRadius_
    << std::endl;
}

KOKKOS_FUNCTION
mm::TransMatType
MotionOscillationCylinderRBF::build_transformation(
  const double& time, const mm::ThreeDVecType& xyz)
{
  mm::TransMatType transMat;

  if (time < startTime_) {
    return transMat;
  }

  // Compute deformation amplitude at current time
  double motionTime = (time < endTime_) ? time : endTime_;
  double angle = 2.0 * M_PI * frequency_ * (motionTime - startTime_);
  double amplitudeAtTime =
    amplitudeFactor_ * cylinderRadius_ * 2.0 * stk::math::sin(angle);

  // Check if point is within deformation domain
  if (xyz[0] < domainExtents_[0] || xyz[0] > domainExtents_[3] ||
      xyz[1] < domainExtents_[1] || xyz[1] > domainExtents_[4] ||
      xyz[2] < domainExtents_[2] || xyz[2] > domainExtents_[5]) {
    return transMat; // Outside domain, no deformation
  }

  const double radius = stk::math::sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]);
  const double blend =
    radial_blend_factor(radius, rigidRadius_, outerBlendRadius_);
  const double deformation_y = amplitudeAtTime * blend;

  // Set y-displacement in transformation matrix (only y-direction)
  transMat[1 * mm::matSize + 3] = deformation_y;

  return transMat;
}

KOKKOS_FUNCTION
mm::ThreeDVecType
MotionOscillationCylinderRBF::compute_velocity(
  const double& time,
  const mm::TransMatType& /* compTrans */,
  const mm::ThreeDVecType& /* mxyz */,
  const mm::ThreeDVecType& cxyz)
{
  if (time < startTime_ || time > endTime_) {
    return mm::ThreeDVecType{0, 0, 0};
  }

  // Check if point is within deformation domain
  if (cxyz[0] < domainExtents_[0] || cxyz[0] > domainExtents_[3] ||
      cxyz[1] < domainExtents_[1] || cxyz[1] > domainExtents_[4] ||
      cxyz[2] < domainExtents_[2] || cxyz[2] > domainExtents_[5]) {
    return mm::ThreeDVecType{0, 0, 0}; // Outside domain, no velocity
  }

  // Compute velocity amplitude using analytical derivative
  // d/dt[A * sin(2*pi*f*t)] = A * 2*pi*f * cos(2*pi*f*t)
  double motionTime = (time < endTime_) ? time : endTime_;
  double angle = 2.0 * M_PI * frequency_ * (motionTime - startTime_);
  double velocityAmplitude =
    amplitudeFactor_ * cylinderRadius_ * 2.0 * 2.0 * M_PI * frequency_ *
    stk::math::cos(angle);

  const double radius = stk::math::sqrt(cxyz[0] * cxyz[0] + cxyz[1] * cxyz[1]);
  const double blend =
    radial_blend_factor(radius, rigidRadius_, outerBlendRadius_);
  const double velocity_y = velocityAmplitude * blend;

  // Velocity only in y-direction
  return mm::ThreeDVecType{0, velocity_y, 0};
}

} // namespace kynema_ugf
} // namespace sierra
