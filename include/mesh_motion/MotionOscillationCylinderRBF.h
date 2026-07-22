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

  void generate_control_points_on_device();
  void assemble_rbf_matrix_on_device();
  void solve_rbf_weights_on_device();

  // Cylinder parameters
  double cylinderRadius_{1.0};
  double cylinderHeight_{4.0};
  double amplitudeFactor_{0.3};
  double frequency_{1.0};
  double rbfBasisParameter_{1.0};

  // Domain extents (min and max for x, y, z)
  double domainExtents_[6]{-8.0, -8.0, 0.0, 8.0, 8.0, 4.0};

  // Number of control points per z-plane on cylinder surface
  int numCylinderControlPointsPerPlane_{36};

  // Number of z-planes for control points (fixed to 5 planes)
  int numControlPlanes_{5};

  // Total number of control points (cylinder + hexahedron)
  int totalControlPoints_{0};

  // Split totals for indexing and RHS assembly
  int totalCylinderControlPoints_{0};
  int totalHexControlPoints_{0};

  // Control points (numCP x 3): each row is [x, y, z]
  Kokkos::View<double**> controlPoints_;

  // RBF matrix (numCP x numCP)
  Kokkos::View<double**> rbfMatrix_;

  // RBF weights for y-direction (numCP x 1)
  Kokkos::View<double*> rbfWeights_;
};

} // namespace kynema_ugf
} // namespace sierra

#endif /* MOTIONOSCILLATIONCYLINDERRBF_H */
