#ifndef MOTIONAIRFOILSMDKERNEL_H
#define MOTIONAIRFOILSMDKERNEL_H

#include <aero/fsi/CalcLoads.h>

#include "NgpMotion.h"

namespace sierra{
namespace nalu{

class MotionAirfoilSMDKernel : public NgpMotionKernel<MotionAirfoilSMDKernel>
{
public:
  MotionAirfoilSMDKernel(const YAML::Node&);

  MotionAirfoilSMDKernel() = default;

  virtual ~MotionAirfoilSMDKernel() = default;

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

protected:
  void load(const YAML::Node&);

  double get_cur_angle(const double time);

  double get_cur_ang_vel(const double time);
  
  mm::ThreeDVecType axis_{0.0, 0.0, 1.0};

  double omega_{0.0};
  double amplt_{0.0};
  double ang_vel_{0.0};
  double phase_{0.0};
  double angle_{0.0};

// private:
//   // Pointer to Algorithm that calculates loads on the surfaces of the Turbine
//   std::unique_ptr<CalcLoads> calc_loads_;
    
};

} // nalu
} //sierra

#endif /* MOTIONAIRFOILSMDKERNEL_H */
