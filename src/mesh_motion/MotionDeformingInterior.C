#include "mesh_motion/MotionDeformingInterior.h"

#include <NaluParsing.h>
#include "utils/ComputeVectorDivergence.h"

// stk_mesh/base/fem
#include <stk_mesh/base/FieldBLAS.hpp>

#include <cmath>

namespace sierra{
namespace nalu{

MotionDeformingInterior::MotionDeformingInterior(
  stk::mesh::MetaData& meta,
  const YAML::Node& node)
  : MotionBase()
{
  load(node);

  // declare divergence of mesh velocity for this motion
  ScalarFieldType *divV = &(meta.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "div_mesh_velocity"));
  stk::mesh::put_field_on_mesh(*divV, meta.universal_part(), nullptr);
  stk::mesh::field_fill(0.0, *divV);
}

void MotionDeformingInterior::load(const YAML::Node& node)
{
  // perturb start and end times with a small value for
  // accurate comparison with floats
  double eps = std::numeric_limits<double>::epsilon();

  get_if_present(node, "start_time", startTime_, startTime_);
  startTime_ = startTime_-eps;

  get_if_present(node, "end_time", endTime_, endTime_);
  endTime_ = endTime_+eps;

  // get lower bounds of deforming part of mesh
  if( !node["xyz_min"] )
    NaluEnv::self().naluOutputP0() << "MotionDeformingInterior: Need to define lower bounds of mesh that deform" << std::endl;
  xyzMin_ = node["xyz_min"].as<ThreeDVecType>();

  // get upper bounds of deforming part of mesh
  if( !node["xyz_max"] )
    NaluEnv::self().naluOutputP0() << "MotionDeformingInterior: Need to define upper bounds of mesh that deform" << std::endl;
  xyzMax_ = node["xyz_max"].as<ThreeDVecType>();

  // get amplitude it was defined
  if( node["amplitude"] )
    amplitude_ = node["amplitude"].as<ThreeDVecType>();

  // get amplitude it was defined
  if( node["frequency"] )
    frequency_ = node["frequency"].as<ThreeDVecType>();

  // get origin based on if it was defined
  if( node["centroid"] )
    origin_ = node["centroid"].as<ThreeDVecType>();
}

void MotionDeformingInterior::build_transformation(
  const double time,
  const double* xyz)
{
  if(time < (startTime_)) return;

  double motionTime = (time < endTime_)? time : endTime_;

  scaling_mat(motionTime,xyz);
}

void MotionDeformingInterior::scaling_mat(
  const double time,
  const double* xyz)
{
  reset_mat(transMat_);

  // return identity matrix if point is outside bounds
  if( (xyz[0] <= xyzMin_[0]) || (xyz[0] >= xyzMax_[0]) ||
      (xyz[1] <= xyzMin_[1]) || (xyz[1] >= xyzMax_[1]) ||
      (xyz[2] <= xyzMin_[2]) || (xyz[2] >= xyzMax_[2])  )
    return;

  double eps = std::numeric_limits<double>::epsilon();

  double radius_x = std::abs(xyz[0]-origin_[0]);
  double radius_y = std::abs(xyz[1]-origin_[1]);
  double radius_z = std::abs(xyz[2]-origin_[2]);

  double curr_radius_x = radius_x + amplitude_[0]*(1 - std::cos(2*M_PI*frequency_[0]*time));
  double curr_radius_y = radius_y + amplitude_[1]*(1 - std::cos(2*M_PI*frequency_[1]*time));
  double curr_radius_z = radius_z + amplitude_[2]*(1 - std::cos(2*M_PI*frequency_[2]*time));

  double scaling_x = curr_radius_x/radius_x; if(radius_x <= eps) scaling_x = 1.0;
  double scaling_y = curr_radius_y/radius_y; if(radius_y <= eps) scaling_y = 1.0;
  double scaling_z = curr_radius_z/radius_z; if(radius_z <= eps) scaling_z = 1.0;

  // Build matrix for translating object to cartesian origin
  transMat_[0][3] = -origin_[0];
  transMat_[1][3] = -origin_[1];
  transMat_[2][3] = -origin_[2];

  // Build matrix for scaling object
  TransMatType currTransMat = {};

  currTransMat[0][0] = scaling_x;
  currTransMat[1][1] = scaling_y;
  currTransMat[2][2] = scaling_z;
  currTransMat[3][3] = 1.0;

  // composite addition of motions in current group
  transMat_ = add_motion(currTransMat,transMat_);

  // Build matrix for translating object back to its origin
  reset_mat(currTransMat);
  currTransMat[0][3] = origin_[0];
  currTransMat[1][3] = origin_[1];
  currTransMat[2][3] = origin_[2];

  // composite addition of motions
  transMat_ = add_motion(currTransMat,transMat_);
}

MotionBase::ThreeDVecType MotionDeformingInterior::compute_velocity(
  const double /* time */,
  const TransMatType&  /* compTrans */,
  const double* /* mxyz */,
  const double* /* cxyz */ )
{
  ThreeDVecType vel = {};
  return vel;
}

void MotionDeformingInterior::post_compute_geometry(
  stk::mesh::BulkData& bulk,
  stk::mesh::PartVector& partVec,
  stk::mesh::PartVector& partVecBc,
  bool& computedMeshVelDiv)
{
  return;
}

} // nalu
} // sierra
