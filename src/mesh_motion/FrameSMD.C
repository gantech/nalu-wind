#include "mesh_motion/FrameSMD.h"

#include "FieldTypeDef.h"
#include "mesh_motion/MotionAirfoilSMDKernel.h"
#include <NaluEnv.h>
#include "NaluParsing.h"
#include "ngp_utils/NgpLoopUtils.h"
#include "ngp_utils/NgpReducers.h"
#include "ngp_utils/NgpTypes.h"
#include "utils/ComputeVectorDivergence.h"
#include "mesh_motion/SMD.h"
#include "mesh_motion/AirfoilSMD.h"
#include <ngp_utils/NgpFieldManager.h>
#include "ngp_utils/NgpMeshInfo.h"

// stk_mesh/base/fem
#include <stk_mesh/base/FieldBLAS.hpp>
#include "stk_mesh/base/GetNgpMesh.hpp"

namespace sierra {
namespace nalu {

FrameSMD::FrameSMD(std::shared_ptr<stk::mesh::BulkData> bulk, const YAML::Node& node)
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

    if (num_motions > 1) {
        throw std::runtime_error("FrameSMD: num_motions is more than 1. Only one motion is supported");
    }

    motionKernels_.resize(num_motions);
    smd_.resize(num_motions);

    // create the classes associated with every motion in current group
    for (int i = 0; i < num_motions; i++) {

      // get the motion definition for i-th transformation
      const auto& motion_def = motions[i];

      get_if_present(motion_def, "loads_scale", loads_scale_);

      get_if_present(motion_def, "mesh_transition_start", mesh_ramp_lower_);
      get_if_present(motion_def, "mesh_transition_end", mesh_ramp_upper_);

      if (mesh_ramp_lower_ >= mesh_ramp_upper_) {
          throw std::runtime_error("FrameSMD: Mesh transition should start at a lower value that it ends. "
                                   "Requires mesh_transition_start < mesh_transition_end.");
      }
      
      get_if_present(motion_def, "load_transition_start", load_ramp_lower_);
      get_if_present(motion_def, "load_transition_end", load_ramp_upper_);

      if (load_ramp_lower_ >= load_ramp_upper_) {
          throw std::runtime_error("FrameSMD: Load fade in required to start at earlier time than it ends. "
                                   "Requires load_transition_start < load_transition_end.");
      }

      // motion type should always be defined by the user
      std::string type;
      get_required(motion_def, "type", type);

      // determine type of mesh motion based on user definition in input file
      if (type == "airfoil_smd") {
        motionKernels_[i].reset(
          new MotionAirfoilSMDKernel(motion_def));
        smd_[i].reset(
          new AirfoilSMD(motion_def));
      } else
        throw std::runtime_error(
          "FrameSMD: Invalid mesh motion type: " + type);

      if ( !bulk_->parallel_rank())
        smd_[i]->prepare_nc_file();

    } // end for loop - i index
  }
}

void
FrameSMD::setup(const double dt, std::shared_ptr<stk::mesh::BulkData> bulk)
{
  // compute and set centroid if requested
  if (computeCentroid_) {
    mm::ThreeDVecType computedCentroid;
    compute_centroid_on_parts(computedCentroid);
    set_computed_centroid(computedCentroid);
  }

  calc_loads_ = std::make_unique<CalcLoads>(partVecBc_);
  calc_loads_->setup(bulk_);

  for (auto& i_smd : smd_) 
    i_smd->setup(dt);
}

void FrameSMD::initialize()
{
  calc_loads_->initialize();
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
  const auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulk_);
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
  stk::mesh::NgpField<double> ndtw =
      stk::mesh::get_updated_ngp_field<double>(
          *meta_.get_field<ScalarFieldType>(entityRank, "minimum_distance_to_wall"));

  // sync fields to device
  modelCoords.sync_to_device();
  currCoords.sync_to_device();
  displacement.sync_to_device();
  meshVelocity.sync_to_device();
  ndtw.sync_to_device();

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

        vs::Vector trans_disp = smd_[i]->get_trans_disp();
        const double rot_angle = smd_[i]->get_rot_disp();

        vs::Vector axis = smd_[i]->get_rot_axis();
        vs::Vector origin = smd_[i]->get_origin();
        // build and get transformation matrix
        mm::TransMatType currTransMat = kernel->build_transformation(time, trans_disp, origin, axis, rot_angle);

        // composite addition of motions in current group
        compTransMat = kernel->add_motion(currTransMat, compTransMat);
      }

      double wdist = ndtw.get(mi,0);
      double odist = stk::math::sqrt(modelCoords.get(mi,0) * modelCoords.get(mi,0) +
                                     modelCoords.get(mi,1) * modelCoords.get(mi,1) +
                                     modelCoords.get(mi,2) * modelCoords.get(mi,2));
      double mesh_ramp_func = ramp_function(wdist, mesh_ramp_lower_, mesh_ramp_upper_);
      
      // perform matrix multiplication between transformation matrix
      // and old coordinates to obtain current coordinates
      for (int d = 0; d < nDim; ++d) {
        double cur_coord = compTransMat[d * mm::matSize + 0] * mX[0] +
                           compTransMat[d * mm::matSize + 1] * mX[1] +
                           compTransMat[d * mm::matSize + 2] * mX[2] +
                           compTransMat[d * mm::matSize + 3];

        displacement.get(mi, d) =
            (cur_coord - modelCoords.get(mi, d)) * mesh_ramp_func;
        
        currCoords.get(mi, d) = modelCoords.get(mi, d) + displacement.get(mi, d);

      } // end for loop - d index

      // compute velocity vector on current node resulting from all
      // motions in current motion frame
      for (size_t i = 0; i < numKernels; ++i) {
        NgpMotion* kernel = ngpKernels(i);


        vs::Vector trans_vel = smd_[i]->get_trans_vel();
        vs::Vector rot_vel = smd_[i]->get_rot_vel();
        vs::Vector origin = smd_[i]->get_origin();
        // evaluate velocity associated with motion
        mm::ThreeDVecType mm_vel =
            kernel->compute_velocity(time, mX, origin, trans_vel, rot_vel);

        for (int d = 0; d < nDim; ++d)
            meshVelocity.get(mi, d) += mm_vel[d] * mesh_ramp_func;
        //* stk::math::tanh( -0.5 * (1.0 - ndtw.get(mi,0)-20.0)/40.0 );
        //* exp(- (ndtw.get(mi, 0)-10.0) * (ndtw.get(mi, 0)-10.0)/100.0);
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
        *bulk_, partVec_, partVecBc_, faceVelMag, meshDivVelocity);
    } else {
      compute_scalar_divergence(
        *bulk_, partVec_, partVecBc_, faceVelMag, meshDivVelocity);
    }

    const auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulk_);
    const stk::mesh::EntityRank entityRank = stk::topology::NODE_RANK;

    stk::mesh::NgpField<double> divmeshvel =  
        stk::mesh::get_updated_ngp_field<double>(
            *meta_.get_field<ScalarFieldType>(entityRank, "div_mesh_velocity"));
    ScalarFieldType* dnv = meta_.get_field<ScalarFieldType>(entityRank, "dual_nodal_volume");
    stk::mesh::NgpField<double> dnvNp1 =  
        stk::mesh::get_updated_ngp_field<double>(dnv->field_of_state(stk::mesh::StateNP1));
    stk::mesh::NgpField<double> dnvN =  
        stk::mesh::get_updated_ngp_field<double>(dnv->field_of_state(stk::mesh::StateN));
    stk::mesh::NgpField<double> dnvNm1 =  
        stk::mesh::get_updated_ngp_field<double>(dnv->field_of_state(stk::mesh::StateNM1));
    dnvNp1.sync_to_device();
    dnvN.sync_to_device();
    dnvNm1.sync_to_device();    
    divmeshvel.sync_to_device();
    double gamma1 = 3.0/2.0;
    double gamma2 = -2.0;
    double gamma3 = 1.0/2.0;
    double dt = 0.001;
    
    // get the parts in the current motion frame
    stk::mesh::Selector sel = meta_.locally_owned_part() | meta_.globally_shared_part();
    // NGP for loop to update coordinates and velocity
    nalu_ngp::run_entity_algorithm(
        "FrameSMD_compute_div_mesh_vel", ngpMesh, entityRank, sel,
        KOKKOS_LAMBDA(
            const nalu_ngp::NGPMeshTraits<stk::mesh::NgpMesh>::MeshIndex& mi) {

            double volume_contrib = (gamma1 * dnvNp1.get(mi,0) +
                                     gamma2 * dnvN.get(mi,0) +
                                     gamma3 * dnvNm1.get(mi,0) ) / dt;
            divmeshvel.get(mi,0) -= volume_contrib;
    }); // end NGP for loop
    //divmeshvel.modify_on_device();
        
    // Mesh velocity divergence is not motion-specific and
    // is computed for the aggregated mesh velocity
    break;
  }
}

void
FrameSMD::predict_states()
{
  for (auto& i_smd : smd_)
    i_smd->predict_states();
}

void
FrameSMD::update_timestep(double cur_time)
{
  calc_loads_->execute();
  for (auto& i_smd : smd_) {
    // Calc 6DOF forces here and pass to
    vs::Vector fnp1;
    vs::Vector mnp1;
    calc_loads_->calc_force_moment(i_smd->get_origin(), fnp1, mnp1);

    double load_ramp_func = ramp_function(cur_time, load_ramp_lower_, load_ramp_upper_);

    i_smd->update_timestep(loads_scale_ * (1 - load_ramp_func) * fnp1, 
                           loads_scale_ * (1 - load_ramp_func) * mnp1);
  }
}

void
FrameSMD::advance_timestep(const double cur_time)
{
  for (auto& i_smd : smd_) {
    i_smd->advance_timestep();
    
    if ( !bulk_->parallel_rank())
      i_smd->write_nc_def_loads(cur_time);
  }
  
}    

double 
FrameSMD::ramp_function(double position, const double start, const double end)
{
  double ramp_func = 1.0;

  if (position < start) {
      ramp_func = 1.0;
  } else if (position < end) {
      ramp_func = 1.0 - 3.0 * stk::math::pow((position-start)/(end-start), 2) 
                      + 2.0 * stk::math::pow((position-start)/(end-start), 3);
  } else {
      ramp_func = 0.0;
  }

  return ramp_func;
}

} // namespace nalu
} // namespace sierra
