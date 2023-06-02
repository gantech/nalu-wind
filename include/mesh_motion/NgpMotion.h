#ifndef NGPMOTION_H
#define NGPMOTION_H

#include <algorithm>
#include <cassert>
#include <cfloat>

#include "NGPInstance.h"
#include "ngp_utils/NgpTypes.h"
#include "vs/vector.h"

namespace YAML {
class Node;
}

namespace sierra {
namespace nalu {

namespace mm {

static constexpr int matSize = nalu_ngp::NDimMax + 1;

struct TransMatType
{
  double Mat_[matSize * matSize];

  // initialize matrix to be an identity matrix
  KOKKOS_FORCEINLINE_FUNCTION
  constexpr TransMatType()
    : Mat_{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}
  {
  }

  // initialize matrix using all components
  KOKKOS_FORCEINLINE_FUNCTION
  constexpr TransMatType(
    const double& x1,
    const double& x2,
    const double& x3,
    const double& x4,
    const double& y1,
    const double& y2,
    const double& y3,
    const double& y4,
    const double& z1,
    const double& z2,
    const double& z3,
    const double& z4,
    const double& c1,
    const double& c2,
    const double& c3,
    const double& c4)
    : Mat_{x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, c1, c2, c3, c4}
  {
  }

  // bracket access operators
  KOKKOS_FORCEINLINE_FUNCTION
  double& operator[](int pos) { return Mat_[pos]; }
  KOKKOS_FORCEINLINE_FUNCTION
  const double& operator[](int pos) const { return Mat_[pos]; }

  // zeroed out matrix
  KOKKOS_FORCEINLINE_FUNCTION
  static constexpr TransMatType zero()
  {
    return TransMatType{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  }

  // Identity matrix
  KOKKOS_FORCEINLINE_FUNCTION
  static constexpr TransMatType I()
  {
    return TransMatType{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  }
};

struct ThreeDVecType
{
  double Vec_[nalu_ngp::NDimMax];

  // initialize empty vector
  KOKKOS_FORCEINLINE_FUNCTION
  constexpr ThreeDVecType() : Vec_{0, 0, 0} {}

  // initialize vector using all components
  KOKKOS_FORCEINLINE_FUNCTION
  ThreeDVecType(const double& x, const double& y, const double& z)
    : Vec_{x, y, z}
  {
  }

  // bracket access operators
  KOKKOS_FORCEINLINE_FUNCTION
  double& operator[](int pos) { return Vec_[pos]; }
  KOKKOS_FORCEINLINE_FUNCTION
  const double& operator[](int pos) const { return Vec_[pos]; }
};

} // namespace mm

class NgpMotion
{
public:
  NgpMotion() = default;

  virtual ~NgpMotion() = default;

  virtual NgpMotion* create_on_device() = 0;

  virtual void free_on_device() = 0;

  /** Function to compute motion-specific transformation matrix
   *
   * @param[in] time Current time
   * @param[in] xyz  Coordinates
   * @return Transformation matrix
   */
  KOKKOS_FUNCTION
  virtual mm::TransMatType
  build_transformation(const double& time, const mm::ThreeDVecType& xyz) = 0;


  /** Function to compute motion-specific transformation matrix
   *
   * @param[in] time Current time
   * @param[in] trans_disp Translation displacements
   * @param[in] origin Origin of rotation
   * @param[in] axis Axis of rotation
   * @param[in] angle Angle of rotation in radians
   * @return Transformation matrix
   */
  KOKKOS_FUNCTION
  mm::TransMatType
  build_transformation(
    const double& time,
    const vs::Vector trans_disp,
    const vs::Vector origin,
    const vs::Vector axis,
    const double rot_angle)
  {
    mm::TransMatType transMat;
    
    // Build matrix for translating object to cartesian origin
    transMat[0 * mm::matSize + 3] = -origin[0];
    transMat[1 * mm::matSize + 3] = -origin[1];
    transMat[2 * mm::matSize + 3] = -origin[2];
  
    // Build matrix for rotating object
    // compute magnitude of axis around which to rotate
    double mag = 0.0;
    for (int d = 0; d < nalu_ngp::NDimMax; d++)
      mag += axis[d] * axis[d];
    mag = stk::math::sqrt(mag);
  
    // build quaternion based on angle and axis of rotation
    const double cosang = stk::math::cos(0.5 * rot_angle);
    const double sinang = stk::math::sin(0.5 * rot_angle);
    const double q0 = cosang;
    const double q1 = sinang * axis[0] / mag;
    const double q2 = sinang * axis[1] / mag;
    const double q3 = sinang * axis[2] / mag;
  
    // rotation matrix based on quaternion
    mm::TransMatType tempMat;
    // 1st row
    tempMat[0 * mm::matSize + 0] = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
    tempMat[0 * mm::matSize + 1] = 2.0 * (q1 * q2 - q0 * q3);
    tempMat[0 * mm::matSize + 2] = 2.0 * (q0 * q2 + q1 * q3);
    // 2nd row
    tempMat[1 * mm::matSize + 0] = 2.0 * (q1 * q2 + q0 * q3);
    tempMat[1 * mm::matSize + 1] = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
    tempMat[1 * mm::matSize + 2] = 2.0 * (q2 * q3 - q0 * q1);
    // 3rd row
    tempMat[2 * mm::matSize + 0] = 2.0 * (q1 * q3 - q0 * q2);
    tempMat[2 * mm::matSize + 1] = 2.0 * (q0 * q1 + q2 * q3);
    tempMat[2 * mm::matSize + 2] = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
  
    // composite addition of motions in current group
    transMat = add_motion(tempMat, transMat);
  
    // Build matrix for translating object back to its origin
    tempMat = mm::TransMatType::I();
    tempMat[0 * mm::matSize + 3] = origin[0];
    tempMat[1 * mm::matSize + 3] = origin[1];
    tempMat[2 * mm::matSize + 3] = origin[2];
  
    // composite addition of motions
    return add_motion(tempMat, transMat);
  }
    

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
    const mm::ThreeDVecType& cxyz) = 0;

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
  mm::ThreeDVecType
  compute_velocity(
    const double& time,
    const mm::ThreeDVecType& mxyz,
    const vs::Vector origin,
    const vs::Vector trans_vel,
    const vs::Vector rot_vel)
  {
    mm::ThreeDVecType vel;

    // compute relative coords and vector omega (dimension 3) for general cross
    // product
    mm::ThreeDVecType rel_coord;
    for (int d = 0; d < nalu_ngp::NDimMax; d++)
      rel_coord[d] = mxyz[d] - origin[d];

    // v = vtrans + \omega \cross \x
    vel[0] = trans_vel[0] + rot_vel[1] * rel_coord[2] - rot_vel[2] * rel_coord[1];
    vel[1] = trans_vel[0] + rot_vel[2] * rel_coord[0] - rot_vel[0] * rel_coord[2];
    vel[2] = trans_vel[0] + rot_vel[0] * rel_coord[1] - rot_vel[1] * rel_coord[0];
    return vel;
  }
      
    
  /** Composite addition of motions
   *
   * @param[in] motionL Left matrix in composite transformation of matrices
   * @param[in] motionR Right matrix in composite transformation of matrices
   * @return    4x4 matrix representing composite addition of motions
   */
  KOKKOS_FORCEINLINE_FUNCTION
  mm::TransMatType
  add_motion(const mm::TransMatType& motionL, const mm::TransMatType& motionR)
  {
    mm::TransMatType tempMat = mm::TransMatType::zero();
    for (int r = 0; r < mm::matSize; r++) {
      for (int c = 0; c < mm::matSize; c++) {
        for (int k = 0; k < mm::matSize; k++) {
          tempMat[r * mm::matSize + c] +=
            motionL[r * mm::matSize + k] * motionR[k * mm::matSize + c];
        }
      }
    }
    return tempMat;
  }

  void set_computed_centroid(const mm::ThreeDVecType& centroid)
  {
    for (int d = 0; d < nalu_ngp::NDimMax; ++d)
      origin_[d] = centroid[d];
  }

  bool is_deforming() { return isDeforming_; }

protected:
  /** Centroid
   *
   * A 3x1 vector storing the centroid as computed
   * to a collection of parts or as defined in the input file
   */
  mm::ThreeDVecType origin_;

  double startTime_{0.0};
  double endTime_{DBL_MAX};

  bool isDeforming_ = false;
};

template <typename T>
class NgpMotionKernel : public NgpMotion
{
public:
  NgpMotionKernel() = default;

  virtual ~NgpMotionKernel() = default;

  virtual NgpMotion* create_on_device() final
  {
    free_on_device();
    deviceCopy_ = nalu_ngp::create<T>(*dynamic_cast<T*>(this));
    return deviceCopy_;
  }

  virtual void free_on_device() final
  {
    if (deviceCopy_ != nullptr) {
      nalu_ngp::destroy<T>(dynamic_cast<T*>(deviceCopy_));
      deviceCopy_ = nullptr;
    }
  }

  T* device_copy() const { return deviceCopy_; }

private:
  T* deviceCopy_{nullptr};
};

} // namespace nalu
} // namespace sierra

#endif /* NGPMOTION_H */
