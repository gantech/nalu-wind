#include "mesh_motion/MotionOscillationCylinderRBF.h"

#include <KynemaUGFEnv.h>
#include <KynemaUGFParsing.h>
#include <cmath>
#include <stdexcept>

namespace sierra {
namespace kynema_ugf {

MotionOscillationCylinderRBF::MotionOscillationCylinderRBF(const YAML::Node& node)
  : NgpMotionKernel<MotionOscillationCylinderRBF>()
{
  load(node);
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
  get_if_present(node, "rbf_basis_parameter", rbfBasisParameter_, rbfBasisParameter_);
  get_if_present(
    node, "num_cylinder_control_points", numCylinderControlPoints_,
    numCylinderControlPoints_);

  // Parse domain extents if provided
  if (node["domain_extents"]) {
    for (int i = 0; i < 6; ++i) {
      domainExtents_[i] = node["domain_extents"][i].as<double>();
    }
  }

  // Calculate total control points
  // Cylinder surface: numCylinderControlPoints
  // Hexahedron: 8 corners + 12 edge midpoints + 6 face centers = 26 points
  totalControlPoints_ = numCylinderControlPoints_ + 26;

  // Allocate Kokkos views on device
  controlPoints_ = Kokkos::View<double**>("controlPoints", totalControlPoints_, 3);
  rbfMatrix_ = Kokkos::View<double**>("rbfMatrix", totalControlPoints_, totalControlPoints_);
  rbfWeights_ = Kokkos::View<double*>("rbfWeights", totalControlPoints_);

  // Generate control points on device
  generate_control_points_on_device();

  // Assemble RBF matrix on device
  assemble_rbf_matrix_on_device();

  // Solve for RBF weights on device (time-independent)
  solve_rbf_weights_on_device();

  KynemaUGFEnv::self().kynema_ugfOutputP0()
    << "MotionOscillationCylinderRBF initialized with " << totalControlPoints_
    << " control points, frequency=" << frequency_ << ", amplitude_factor="
    << amplitudeFactor_ << std::endl;
}

// Kernel to generate control points on the cylinder surface and hexahedron boundary
void
MotionOscillationCylinderRBF::generate_control_points_on_device()
{
  auto controlPoints = controlPoints_;
  auto cylinderRadius = cylinderRadius_;
  auto cylinderHeight = cylinderHeight_;
  auto numCylCP = numCylinderControlPoints_;
  auto xMin = domainExtents_[0];
  auto yMin = domainExtents_[1];
  auto zMin = domainExtents_[2];
  auto xMax = domainExtents_[3];
  auto yMax = domainExtents_[4];
  auto zMax = domainExtents_[5];

  // Lambda to generate points on device
  auto generate_lambda = [=]() {
    // Generate cylinder surface control points (uniformly spaced in azimuthal direction)
    int cpPerEnd = numCylCP / 2;
    for (int i = 0; i < cpPerEnd; ++i) {
      double theta = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(cpPerEnd);
      double x = cylinderRadius * std::cos(theta);
      double y = cylinderRadius * std::sin(theta);

      // Point at z = 0
      controlPoints(i, 0) = x;
      controlPoints(i, 1) = y;
      controlPoints(i, 2) = 0.0;

      // Point at z = height
      controlPoints(cpPerEnd + i, 0) = x;
      controlPoints(cpPerEnd + i, 1) = y;
      controlPoints(cpPerEnd + i, 2) = cylinderHeight;
    }

    // Generate hexahedron control points
    int hexCPStart = numCylCP;

    // 8 corners
    int cornerIdx = hexCPStart;
    double corners[8][3] = {
      {xMin, yMin, zMin}, {xMax, yMin, zMin}, {xMax, yMax, zMin}, {xMin, yMax, zMin},
      {xMin, yMin, zMax}, {xMax, yMin, zMax}, {xMax, yMax, zMax}, {xMin, yMax, zMax}};
    for (int i = 0; i < 8; ++i) {
      controlPoints(cornerIdx + i, 0) = corners[i][0];
      controlPoints(cornerIdx + i, 1) = corners[i][1];
      controlPoints(cornerIdx + i, 2) = corners[i][2];
    }

    // 12 edge midpoints
    int edgeIdx = hexCPStart + 8;
    // Bottom face edges (z = zMin)
    controlPoints(edgeIdx + 0, 0) = (xMin + xMax) / 2.0;
    controlPoints(edgeIdx + 0, 1) = yMin;
    controlPoints(edgeIdx + 0, 2) = zMin;
    controlPoints(edgeIdx + 1, 0) = xMax;
    controlPoints(edgeIdx + 1, 1) = (yMin + yMax) / 2.0;
    controlPoints(edgeIdx + 1, 2) = zMin;
    controlPoints(edgeIdx + 2, 0) = (xMin + xMax) / 2.0;
    controlPoints(edgeIdx + 2, 1) = yMax;
    controlPoints(edgeIdx + 2, 2) = zMin;
    controlPoints(edgeIdx + 3, 0) = xMin;
    controlPoints(edgeIdx + 3, 1) = (yMin + yMax) / 2.0;
    controlPoints(edgeIdx + 3, 2) = zMin;

    // Top face edges (z = zMax)
    controlPoints(edgeIdx + 4, 0) = (xMin + xMax) / 2.0;
    controlPoints(edgeIdx + 4, 1) = yMin;
    controlPoints(edgeIdx + 4, 2) = zMax;
    controlPoints(edgeIdx + 5, 0) = xMax;
    controlPoints(edgeIdx + 5, 1) = (yMin + yMax) / 2.0;
    controlPoints(edgeIdx + 5, 2) = zMax;
    controlPoints(edgeIdx + 6, 0) = (xMin + xMax) / 2.0;
    controlPoints(edgeIdx + 6, 1) = yMax;
    controlPoints(edgeIdx + 6, 2) = zMax;
    controlPoints(edgeIdx + 7, 0) = xMin;
    controlPoints(edgeIdx + 7, 1) = (yMin + yMax) / 2.0;
    controlPoints(edgeIdx + 7, 2) = zMax;

    // Vertical edges (connecting bottom to top)
    controlPoints(edgeIdx + 8, 0) = xMin;
    controlPoints(edgeIdx + 8, 1) = yMin;
    controlPoints(edgeIdx + 8, 2) = (zMin + zMax) / 2.0;
    controlPoints(edgeIdx + 9, 0) = xMax;
    controlPoints(edgeIdx + 9, 1) = yMin;
    controlPoints(edgeIdx + 9, 2) = (zMin + zMax) / 2.0;
    controlPoints(edgeIdx + 10, 0) = xMax;
    controlPoints(edgeIdx + 10, 1) = yMax;
    controlPoints(edgeIdx + 10, 2) = (zMin + zMax) / 2.0;
    controlPoints(edgeIdx + 11, 0) = xMin;
    controlPoints(edgeIdx + 11, 1) = yMax;
    controlPoints(edgeIdx + 11, 2) = (zMin + zMax) / 2.0;

    // 6 face centers
    int faceIdx = hexCPStart + 20;
    // Bottom face
    controlPoints(faceIdx + 0, 0) = (xMin + xMax) / 2.0;
    controlPoints(faceIdx + 0, 1) = (yMin + yMax) / 2.0;
    controlPoints(faceIdx + 0, 2) = zMin;
    // Top face
    controlPoints(faceIdx + 1, 0) = (xMin + xMax) / 2.0;
    controlPoints(faceIdx + 1, 1) = (yMin + yMax) / 2.0;
    controlPoints(faceIdx + 1, 2) = zMax;
    // Front face
    controlPoints(faceIdx + 2, 0) = (xMin + xMax) / 2.0;
    controlPoints(faceIdx + 2, 1) = yMin;
    controlPoints(faceIdx + 2, 2) = (zMin + zMax) / 2.0;
    // Back face
    controlPoints(faceIdx + 3, 0) = (xMin + xMax) / 2.0;
    controlPoints(faceIdx + 3, 1) = yMax;
    controlPoints(faceIdx + 3, 2) = (zMin + zMax) / 2.0;
    // Left face
    controlPoints(faceIdx + 4, 0) = xMin;
    controlPoints(faceIdx + 4, 1) = (yMin + yMax) / 2.0;
    controlPoints(faceIdx + 4, 2) = (zMin + zMax) / 2.0;
    // Right face
    controlPoints(faceIdx + 5, 0) = xMax;
    controlPoints(faceIdx + 5, 1) = (yMin + yMax) / 2.0;
    controlPoints(faceIdx + 5, 2) = (zMin + zMax) / 2.0;
  };

  // Execute on host (simple point generation doesn't require device)
  generate_lambda();
}

// Kernel to assemble RBF matrix with exponential basis function
void
MotionOscillationCylinderRBF::assemble_rbf_matrix_on_device()
{
  auto controlPoints = controlPoints_;
  auto rbfMatrix = rbfMatrix_;
  auto rbfBasisParam = rbfBasisParameter_;
  int numCP = totalControlPoints_;

  // Parallel for to fill RBF matrix
  Kokkos::parallel_for(
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {numCP, numCP}),
    KOKKOS_LAMBDA(const int i, const int j) {
      // Compute distance between control points i and j
      double dx = controlPoints(i, 0) - controlPoints(j, 0);
      double dy = controlPoints(i, 1) - controlPoints(j, 1);
      double dz = controlPoints(i, 2) - controlPoints(j, 2);
      double dist_sq = dx * dx + dy * dy + dz * dz;

      // Exponential RBF basis function
      rbfMatrix(i, j) = std::exp(-rbfBasisParam * dist_sq);
    });

  Kokkos::fence();
}

// Simple LU decomposition solver on device
void
MotionOscillationCylinderRBF::solve_rbf_weights_on_device()
{
  int numCP = totalControlPoints_;

  // Create RHS vector for y-direction deformation
  // Cylinder surface control points: deformation = 0.3 * diameter * sin(2*pi*f*t)
  // Since we're solving time-independent, we solve for unit amplitude = 1.0
  // Then scale by actual amplitude when applying in build_transformation
  Kokkos::View<double*> rhs("rhs", numCP);

  auto numCylCP = numCylinderControlPoints_;
  auto controlPoints = controlPoints_;

  // Set RHS: 1.0 for cylinder surface points, 0.0 for hexahedron boundary points
  Kokkos::parallel_for(
    numCP, KOKKOS_LAMBDA(const int i) {
      if (i < numCylCP) {
        rhs(i) = 1.0; // Cylinder surface points get unit deformation
      } else {
        rhs(i) = 0.0; // Hexahedron boundary points get zero deformation
      }
    });

  Kokkos::fence();

  // Copy RBF matrix to host for LU decomposition (simple serial solver)
  auto rbfMatrix_host =
    Kokkos::create_mirror_view(Kokkos::HostSpace(), rbfMatrix_);
  Kokkos::deep_copy(rbfMatrix_host, rbfMatrix_);

  auto rhs_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), rhs);
  Kokkos::deep_copy(rhs_host, rhs);

  auto rbfWeights_host =
    Kokkos::create_mirror_view(Kokkos::HostSpace(), rbfWeights_);

  // LU decomposition and forward/backward substitution
  // Copy RHS to weights for in-place solution
  for (int i = 0; i < numCP; ++i) {
    rbfWeights_host(i) = rhs_host(i);
  }

  // Perform LU decomposition with partial pivoting
  std::vector<int> pivots(numCP);
  for (int i = 0; i < numCP; ++i) {
    pivots[i] = i;
  }

  // Forward elimination with partial pivoting
  for (int k = 0; k < numCP; ++k) {
    // Find pivot
    int pivot_row = k;
    double max_val = std::abs(rbfMatrix_host(k, k));
    for (int i = k + 1; i < numCP; ++i) {
      if (std::abs(rbfMatrix_host(i, k)) > max_val) {
        max_val = std::abs(rbfMatrix_host(i, k));
        pivot_row = i;
      }
    }

    if (std::abs(max_val) < 1e-14) {
      throw std::runtime_error(
        "MotionOscillationCylinderRBF: RBF matrix is singular!");
    }

    // Swap rows in matrix, RHS, and pivots
    if (pivot_row != k) {
      for (int j = k; j < numCP; ++j) {
        std::swap(rbfMatrix_host(k, j), rbfMatrix_host(pivot_row, j));
      }
      std::swap(rbfWeights_host(k), rbfWeights_host(pivot_row));
      std::swap(pivots[k], pivots[pivot_row]);
    }

    // Eliminate below
    for (int i = k + 1; i < numCP; ++i) {
      double factor = rbfMatrix_host(i, k) / rbfMatrix_host(k, k);
      rbfMatrix_host(i, k) = 0.0;
      for (int j = k + 1; j < numCP; ++j) {
        rbfMatrix_host(i, j) -= factor * rbfMatrix_host(k, j);
      }
      rbfWeights_host(i) -= factor * rbfWeights_host(k);
    }
  }

  // Back substitution
  for (int i = numCP - 1; i >= 0; --i) {
    double sum = 0.0;
    for (int j = i + 1; j < numCP; ++j) {
      sum += rbfMatrix_host(i, j) * rbfWeights_host(j);
    }
    rbfWeights_host(i) = (rbfWeights_host(i) - sum) / rbfMatrix_host(i, i);
  }

  // Copy solution back to device
  Kokkos::deep_copy(rbfWeights_, rbfWeights_host);

  KynemaUGFEnv::self().kynema_ugfOutputP0()
    << "MotionOscillationCylinderRBF: RBF weights solved successfully"
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

  // Compute RBF interpolation for y-direction deformation
  double deformation_y = 0.0;
  for (int i = 0; i < totalControlPoints_; ++i) {
    double dx = xyz[0] - controlPoints_(i, 0);
    double dy = xyz[1] - controlPoints_(i, 1);
    double dz = xyz[2] - controlPoints_(i, 2);
    double dist_sq = dx * dx + dy * dy + dz * dz;

    // Exponential RBF basis function
    double phi_i = stk::math::exp(-rbfBasisParameter_ * dist_sq);

    // Accumulate weighted contribution
    deformation_y += rbfWeights_(i) * phi_i * amplitudeAtTime;
  }

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

  // Compute RBF interpolation for y-direction velocity
  double velocity_y = 0.0;
  for (int i = 0; i < totalControlPoints_; ++i) {
    double dx = cxyz[0] - controlPoints_(i, 0);
    double dy = cxyz[1] - controlPoints_(i, 1);
    double dz = cxyz[2] - controlPoints_(i, 2);
    double dist_sq = dx * dx + dy * dy + dz * dz;

    // Exponential RBF basis function
    double phi_i = stk::math::exp(-rbfBasisParameter_ * dist_sq);

    // Accumulate weighted contribution
    velocity_y += rbfWeights_(i) * phi_i * velocityAmplitude;
  }

  // Velocity only in y-direction
  return mm::ThreeDVecType{0, velocity_y, 0};
}

} // namespace kynema_ugf
} // namespace sierra
