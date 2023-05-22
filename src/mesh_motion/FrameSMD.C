#include "mesh_motion/FrameSMD.h"

#include "FieldTypeDef.h"
#include "mesh_motion/MotionAirfoilSMDKernel.h"
#include "NaluParsing.h"
#include "ngp_utils/NgpLoopUtils.h"
#include "ngp_utils/NgpReducers.h"
#include "ngp_utils/NgpTypes.h"
#include "utils/ComputeVectorDivergence.h"

// stk_mesh/base/fem
#include <stk_mesh/base/FieldBLAS.hpp>
#include "stk_mesh/base/GetNgpMesh.hpp"

namespace sierra {
namespace nalu {

FrameSMD::FrameSMD(stk::mesh::BulkData& bulk, const YAML::Node& node)
  : FrameBase(bulk)
{
  load(node);

  // set deformation flag based on motions in the frame
  for (auto& mm : motionKernels_)
    if (mm->is_deforming())
      isDeforming_ = true;
}

FrameSMD::~FrameSMD()
{
  // Release the device pointers if any
  for (auto& kern : motionKernels_) {
    kern->free_on_device();
  }
}

void
FrameSMD::load(const YAML::Node& node)
{
  // get any part names associated with current motion group
  populate_part_vec(node);

  // check if centroid needs to be computed
  get_if_present(node, "compute_centroid", computeCentroid_, computeCentroid_);

  if (node["motion"]) {
    // extract the motions in the current group
    const auto& motions = node["motion"];

    const int num_motions = motions.size();
    motionKernels_.resize(num_motions);

    // create the classes associated with every motion in current group
    for (int i = 0; i < num_motions; i++) {

      // get the motion definition for i-th transformation
      const auto& motion_def = motions[i];

      // motion type should always be defined by the user
      std::string type;
      get_required(motion_def, "smd_type", type);

      // determine type of mesh motion based on user definition in input file
      if (type == "airfoil_smd")
        motionKernels_[i].reset(
          new MotionAirfoilSMDKernel(motion_def));
      else
        throw std::runtime_error(
          "FrameSMD: Invalid mesh motion type: " + type);

    } // end for loop - i index
  }
}

void
FrameSMD::setup()
{
  // compute and set centroid if requested
  if (computeCentroid_) {
    mm::ThreeDVecType computedCentroid;
    compute_centroid_on_parts(computedCentroid);
    set_computed_centroid(computedCentroid);
  }
}

void
FrameSMD::update_coordinates_velocity(const double time)
{
  assert(partVec_.size() > 0);

  // create NGP view of motion kernels
  const size_t numKernels = motionKernels_.size();
  auto ngpKernels = nalu_ngp::create_ngp_view<NgpMotion>(motionKernels_);

  // define mesh entities
  const int nDim = meta_.spatial_dimension();
  const auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk_);
  const stk::mesh::EntityRank entityRank = stk::topology::NODE_RANK;

  // get the parts in the current motion frame
  stk::mesh::Selector sel =
    stk::mesh::selectUnion(partVec_) &
    (meta_.locally_owned_part() | meta_.globally_shared_part());

  // get the field from the NGP mesh
  stk::mesh::NgpField<double> modelCoords =
    stk::mesh::get_updated_ngp_field<double>(
      *meta_.get_field<VectorFieldType>(entityRank, "coordinates"));
  stk::mesh::NgpField<double> currCoords =
    stk::mesh::get_updated_ngp_field<double>(
      *meta_.get_field<VectorFieldType>(entityRank, "current_coordinates"));
  stk::mesh::NgpField<double> displacement =
    stk::mesh::get_updated_ngp_field<double>(
      *meta_.get_field<VectorFieldType>(entityRank, "mesh_displacement"));
  stk::mesh::NgpField<double> meshVelocity =
    stk::mesh::get_updated_ngp_field<double>(
      *meta_.get_field<VectorFieldType>(entityRank, "mesh_velocity"));

  // sync fields to device
  modelCoords.sync_to_device();
  currCoords.sync_to_device();
  displacement.sync_to_device();
  meshVelocity.sync_to_device();

  // always reset velocity field
  nalu_ngp::run_entity_algorithm(
    "FrameSMD_reset_velocity", ngpMesh, entityRank, sel,
    KOKKOS_LAMBDA(
      const nalu_ngp::NGPMeshTraits<stk::mesh::NgpMesh>::MeshIndex& mi) {
      for (int d = 0; d < nDim; ++d)
        meshVelocity.get(mi, d) = 0.0;
    });

  // NGP for loop to update coordinates and velocity
  nalu_ngp::run_entity_algorithm(
    "FrameSMD_update_coordinates_velocity", ngpMesh, entityRank, sel,
    KOKKOS_LAMBDA(
      const nalu_ngp::NGPMeshTraits<stk::mesh::NgpMesh>::MeshIndex& mi) {
      // temporary current and model coords for a generic 2D and 3D
      // implementation
      mm::ThreeDVecType mX;
      mm::ThreeDVecType cX;

      // copy over model coordinates and reset velocity
      for (int d = 0; d < nDim; ++d)
        mX[d] = modelCoords.get(mi, d);

      // initialize composite transformation matrix
      mm::TransMatType compTransMat;

      // create composite transformation matrix based off of all motions
      for (size_t i = 0; i < numKernels; ++i) {
        NgpMotion* kernel = ngpKernels(i);

        // build and get transformation matrix
        mm::TransMatType currTransMat = kernel->build_transformation(time, mX);

        // composite addition of motions in current group
        compTransMat = kernel->add_motion(currTransMat, compTransMat);
      }

      // perform matrix multiplication between transformation matrix
      // and old coordinates to obtain current coordinates
      for (int d = 0; d < nDim; ++d) {
        currCoords.get(mi, d) = compTransMat[d * mm::matSize + 0] * mX[0] +
                                compTransMat[d * mm::matSize + 1] * mX[1] +
                                compTransMat[d * mm::matSize + 2] * mX[2] +
                                compTransMat[d * mm::matSize + 3];

        displacement.get(mi, d) =
          currCoords.get(mi, d) - modelCoords.get(mi, d);
      } // end for loop - d index

      // copy over current coordinates
      for (int d = 0; d < nDim; ++d)
        cX[d] = currCoords.get(mi, d);

      // compute velocity vector on current node resulting from all
      // motions in current motion frame
      for (size_t i = 0; i < numKernels; ++i) {
        NgpMotion* kernel = ngpKernels(i);

        // evaluate velocity associated with motion
        mm::ThreeDVecType mm_vel =
          kernel->compute_velocity(time, compTransMat, mX, cX);

        for (int d = 0; d < nDim; ++d)
          meshVelocity.get(mi, d) += mm_vel[d];
      } // end for loop - mm
    }); // end NGP for loop

  // Mark fields as modified on device
  currCoords.modify_on_device();
  displacement.modify_on_device();
  meshVelocity.modify_on_device();
}

void
FrameSMD::post_compute_geometry()
{
  for (auto& mm : motionKernels_) {
    if (!mm->is_deforming())
      continue;

    // compute divergence of mesh velocity
    ScalarFieldType* meshDivVelocity = meta_.get_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "div_mesh_velocity");
    GenericFieldType* faceVelMag = meta_.get_field<GenericFieldType>(
      stk::topology::ELEMENT_RANK, "face_velocity_mag");

    if (faceVelMag == NULL) {
      faceVelMag = meta_.get_field<GenericFieldType>(
        stk::topology::EDGE_RANK, "edge_face_velocity_mag");
      compute_edge_scalar_divergence(
        bulk_, partVec_, partVecBc_, faceVelMag, meshDivVelocity);
    } else {
      compute_scalar_divergence(
        bulk_, partVec_, partVecBc_, faceVelMag, meshDivVelocity);
    }

    // Mesh velocity divergence is not motion-specific and
    // is computed for the aggregated mesh velocity
    break;
  }
}

} // namespace nalu
} // namespace sierra
