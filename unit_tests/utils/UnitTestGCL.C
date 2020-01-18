// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#include <limits>

#include "kernels/UnitTestKernelUtils.h"
#include "UnitTestHelperObjects.h"

#include "AlgTraits.h"
#include "mesh_motion/MeshMotionAlg.h"
#include "ngp_algorithms/GeometryInteriorAlg.h"
#include "ngp_algorithms/GeometryBoundaryAlg.h"
#include "ngp_algorithms/WallFuncGeometryAlg.h"
#include "ngp_algorithms/GeometryAlgDriver.h"
#include "ngp_algorithms/NodalGradAlgDriver.h"
#include "ngp_algorithms/NodalGradElemAlg.h"
#include "ngp_algorithms/NodalGradBndryElemAlg.h"
#include "utils/StkHelpers.h"

namespace {

class GCLTest : public ::testing::Test
{
public:
  GCLTest()
    : naluObj_(),
      realm_(naluObj_.create_realm()),
      meta_(realm_.meta_data()),
      bulk_(realm_.bulk_data()),
      geomAlgDriver_(realm_),
      nodalGradAlgDriver_(realm_, "dvdx"),
      currCoords_(
        &meta_.declare_field<VectorFieldType>(
          stk::topology::NODE_RANK, "current_coordinates", numStates_)),
      dualVol_(
        &meta_.declare_field<ScalarFieldType>(
          stk::topology::NODE_RANK, "dual_nodal_volume", numStates_)),
      elemVol_(
        &meta_.declare_field<ScalarFieldType>(
          stk::topology::ELEM_RANK, "element_volume")),
      edgeAreaVec_(
        &meta_.declare_field<VectorFieldType>(
          stk::topology::EDGE_RANK, "edge_area_vector")),
      exposedAreaVec_(
        &meta_.declare_field<GenericFieldType>(
          meta_.side_rank(), "exposed_area_vector")),
      meshDisp_(
        &meta_.declare_field<VectorFieldType>(
          stk::topology::NODE_RANK, "mesh_displacement", numStates_)),
      meshVel_(
        &meta_.declare_field<VectorFieldType>(
          stk::topology::NODE_RANK, "mesh_velocity", numStates_)),
      dMeshVeldx_(
        &meta_.declare_field<GenericFieldType>(stk::topology::NODE_RANK, "dvdx")),
      divMeshVel_(
        &meta_.declare_field<ScalarFieldType>
        (stk::topology::NODE_RANK, "div_mesh_velocity")),
      dVoldt_(
        &meta_.declare_field<ScalarFieldType>
        (stk::topology::NODE_RANK, "dvol_dt"))
  {
    realm_.timeIntegrator_ = naluObj_.sim_.timeIntegrator_;
    stk::mesh::put_field_on_mesh(
      *currCoords_, meta_.universal_part(), spatialDim_, nullptr);
    stk::mesh::put_field_on_mesh(
      *dualVol_, meta_.universal_part(), 1, nullptr);
    stk::mesh::put_field_on_mesh(
      *elemVol_, meta_.universal_part(), 1, nullptr);
    stk::mesh::put_field_on_mesh(
      *edgeAreaVec_, meta_.universal_part(), spatialDim_, nullptr);
    stk::mesh::put_field_on_mesh(
      *exposedAreaVec_, meta_.universal_part(),
      spatialDim_ * sierra::nalu::AlgTraitsQuad4::numScsIp_, nullptr);
    stk::mesh::put_field_on_mesh(
      *meshDisp_, meta_.universal_part(), spatialDim_, nullptr);
    stk::mesh::put_field_on_mesh(
      *meshVel_, meta_.universal_part(), spatialDim_, nullptr);
    stk::mesh::put_field_on_mesh(
      *dMeshVeldx_, meta_.universal_part(), spatialDim_ * spatialDim_, nullptr);
    stk::mesh::put_field_on_mesh(
      *divMeshVel_, meta_.universal_part(), 1, nullptr);
    stk::mesh::put_field_on_mesh(
      *dVoldt_, meta_.universal_part(), 1, nullptr);
  }

  virtual ~GCLTest() = default;

  void fill_mesh_and_init_fields(
    const std::string meshSize,
    bool doPerturb=true,
    bool generateSidesets=true)
  {
    std::string meshSpec = "generated:" + meshSize;
    if (generateSidesets)
      meshSpec += "|sideset:xXyYzZ";
    unit_test_utils::fill_hex8_mesh(meshSpec, bulk_);
    if (doPerturb)
      unit_test_utils::perturb_coord_hex_8(bulk_);

    partVec_ = {meta_.get_part("block_1")};
    coordinates_ = static_cast<const VectorFieldType*>(
      meta_.coordinate_field());
    EXPECT_TRUE(coordinates_ != nullptr);

    stk::mesh::create_edges(bulk_, meta_.universal_part());
  }

  void init_time_integrator(bool secondOrder=true, double timeStep=0.1)
  {
    auto& timeInt = *realm_.timeIntegrator_;
    timeInt.secondOrderTimeAccurate_ = secondOrder;
    timeInt.timeStepN_ = timeStep;
    timeInt.timeStepNm1_ = timeStep;
    timeInt.timeStepCount_ = 2;
    if (timeInt.secondOrderTimeAccurate_)
      timeInt.compute_gamma();
  }

  void register_algorithms(const std::string& motion_options)
  {
    // Force creation of edge area vector
    realm_.realmUsesEdges_ = true;
    // Force mesh motion logic everywhere
    realm_.solutionOptions_->meshMotion_ = true;
    const YAML::Node motionNode = YAML::Load(motion_options);
    realm_.meshMotionAlg_.reset(
      new sierra::nalu::MeshMotionAlg(bulk_, motionNode["mesh_motion"]));

    const bool useShifted = false;
    auto* part = meta_.get_part("surface_1");
    geomAlgDriver_.register_elem_algorithm<sierra::nalu::GeometryInteriorAlg>(
      sierra::nalu::INTERIOR, partVec_[0], "geometry");
    nodalGradAlgDriver_.register_elem_algorithm<sierra::nalu::VectorNodalGradElemAlg>(
      sierra::nalu::INTERIOR, partVec_[0], "nodal_grad_dvdx",
      meshVel_, dMeshVeldx_, useShifted);

    for (auto* surfPart: part->subsets()) {
      geomAlgDriver_.register_face_algorithm<sierra::nalu::GeometryBoundaryAlg>(
        sierra::nalu::BOUNDARY, surfPart, "geometry");
      nodalGradAlgDriver_.register_face_algorithm<sierra::nalu::VectorNodalGradBndryElemAlg>(
        sierra::nalu::WALL, surfPart, "nodal_grad_dvdx",
        meshVel_, dMeshVeldx_, useShifted);
    }
  }

  void compute_mesh_velocity()
  {
    const stk::mesh::Selector sel = meta_.universal_part();
    const double dt = realm_.get_time_step();
    const double gamma1 = realm_.get_gamma1();
    const double gamma2 = realm_.get_gamma2();
    const double gamma3 = realm_.get_gamma3();

    const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);

    const auto& dispNm1 = meshDisp_->field_of_state(stk::mesh::StateNM1);
    const auto& dispN = meshDisp_->field_of_state(stk::mesh::StateN);
    const auto& dispNp1 = meshDisp_->field_of_state(stk::mesh::StateNP1);

    std::cerr << "Timestep info: dt = " << dt
              << " gamma1 = " << gamma1
              << " gamma2 = " << gamma2
              << " gamma3 = " << gamma3 << std::endl;
    for (auto* b: bkts) {
      const double* dxNm1 = stk::mesh::field_data(dispNm1, *b);
      const double* dxN = stk::mesh::field_data(dispN, *b);
      const double* dxNp1 = stk::mesh::field_data(dispNp1, *b);
      double* mVel = stk::mesh::field_data(*meshVel_, *b);

      for (size_t in=0; in < b->size(); ++in) {
        size_t offset = in * spatialDim_;
        for (unsigned d=0; d < spatialDim_; ++d)
          mVel[offset + d] = (gamma1 * dxNp1[offset + d] +
                              gamma2 * dxN[offset + d] +
                              gamma3 * dxNm1[offset + d]) / dt;
      }
    }
  }

  void compute_div_mesh_vel()
  {
    const stk::mesh::Selector sel = get_selector();
    const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);

    for (auto* b: bkts) {
      const double* dVol = stk::mesh::field_data(*dualVol_, *b);
      const double* dvdx = stk::mesh::field_data(*dMeshVeldx_, *b);
      double* divV = stk::mesh::field_data(*divMeshVel_, *b);

      for (size_t in=0; in < b->size(); ++in) {
        size_t offset = in * spatialDim_ * spatialDim_;
        double sum = 0.0;

        for (unsigned d=0; d < spatialDim_; ++d)
          sum += dvdx[offset + spatialDim_ * d + d];
        divV[in] = sum * dVol[in];

        auto nelems = b->num_elements(in);
        if (nelems == 8)
          std::cerr << "div(meshVel) for interior node = " << divV[in] << std::endl;
      }
    }
  }

  void compute_dvoldt()
  {
    const stk::mesh::Selector sel = get_selector();
    const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);

    const auto& dVolNm1 = dualVol_->field_of_state(stk::mesh::StateNM1);
    const auto& dVolN = dualVol_->field_of_state(stk::mesh::StateN);
    const auto& dVolNp1 = dualVol_->field_of_state(stk::mesh::StateNP1);
    const double dt = realm_.get_time_step();
    const double gamma1 = realm_.get_gamma1();
    const double gamma2 = realm_.get_gamma2();
    const double gamma3 = realm_.get_gamma3();

    double minVal = std::numeric_limits<double>::max();
    double maxVal = std::numeric_limits<double>::lowest();
    for (auto* b: bkts) {
      const double* dvNm1 = stk::mesh::field_data(dVolNm1, *b);
      const double* dvN = stk::mesh::field_data(dVolN, *b);
      const double* dvNp1 = stk::mesh::field_data(dVolNp1, *b);
      double* dvdt = stk::mesh::field_data(*dVoldt_, *b);

      for (size_t in=0; in < b->size(); ++in) {
        dvdt[in] = (gamma1 * dvNp1[in] + gamma2 * dvN[in] + gamma3 * dvNm1[in]) / dt;
        minVal = std::min(dvdt[in], minVal);
        maxVal = std::max(dvdt[in], maxVal);
      }
    }
    std::cerr << "dVol/dt: min = " << minVal << " max = " << maxVal << std::endl;
  }

  void compute_error()
  {
    const stk::mesh::Selector sel = get_selector();
    const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);

    double minErr = std::numeric_limits<double>::max();
    double maxErr = std::numeric_limits<double>::lowest();
    for (auto *b: bkts) {
      const double* divV = stk::mesh::field_data(*divMeshVel_, *b);
      const double* dvdt = stk::mesh::field_data(*dVoldt_, *b);

      for (size_t in=0; in < b->size(); in++) {
        const double err = dvdt[in] - divV[in];
        minErr = std::min(minErr, err);
        maxErr = std::max(maxErr, err);
      }
    }
    std::cerr << "Error: min = " << minErr << " max = " << maxErr << std::endl;
  }

  /** Selector to loop over interior nodes only
   */
  stk::mesh::Selector get_selector()
  {
    stk::mesh::PartVector bdyParts = { meta_.get_part("surface_1")};
    return (
      meta_.universal_part() & !stk::mesh::selectUnion(bdyParts));
  }

  void init_states()
  {
    const double deltaT = realm_.get_time_step();
    auto& motionAlg = *realm_.meshMotionAlg_;
    motionAlg.initialize(0.0);
    for (int it=0; it < numStates_; ++it) {
      realm_.swap_states();
      motionAlg.execute(it * deltaT);
      geomAlgDriver_.execute();
    }
  }

  YAML::Node doc_;
  YAML::Node realmNode_;
  unit_test_utils::NaluTest naluObj_;
  sierra::nalu::Realm& realm_;

  static constexpr int numStates_{3};
  const unsigned spatialDim_{3};
  stk::mesh::MetaData& meta_;
  stk::mesh::BulkData& bulk_;
  sierra::nalu::GeometryAlgDriver geomAlgDriver_;
  sierra::nalu::VectorNodalGradAlgDriver nodalGradAlgDriver_;
  stk::mesh::PartVector partVec_;

  const VectorFieldType* coordinates_{nullptr};
  VectorFieldType* currCoords_{nullptr};
  ScalarFieldType* dualVol_{nullptr};
  ScalarFieldType* elemVol_{nullptr};
  VectorFieldType* edgeAreaVec_{nullptr};
  GenericFieldType* exposedAreaVec_{nullptr};
  VectorFieldType* meshDisp_{nullptr};
  VectorFieldType* meshVel_{nullptr};
  GenericFieldType* dMeshVeldx_{nullptr};
  ScalarFieldType* divMeshVel_{nullptr};
  ScalarFieldType* dVoldt_{nullptr};
};

} // namespace

TEST_F(GCLTest, rigid_rotation)
{
  if (bulk_.parallel_size() > 1) return;

  const std::string meshDims = "2x2x2|offset:0,65,0";
  const bool secondOrder = true;
  const double deltaT = 0.003; // approx 0.25 deg motion for given omega
  const std::string mesh_motion =
    "mesh_motion:                                                          \n"
    "  - name: interior                                                    \n"
    "    frame: non_inertial                                               \n"
    "    mesh_parts: [ block_1 ]                                           \n"
    "    motion:                                                           \n"
    "      - type: rotation                                                \n"
    "        omega: 1.5707963267948966                                     \n"
    "        axis: [1.0, 0.0, 0.0]                                         \n"
    "        centroid: [0.0, 0.0, 0.0]                                     \n";

  fill_mesh_and_init_fields(meshDims);
  init_time_integrator(secondOrder, deltaT);
  register_algorithms(mesh_motion);
  init_states();
  compute_mesh_velocity();
  nodalGradAlgDriver_.execute();
  compute_div_mesh_vel();
  compute_dvoldt();
  compute_error();
}

TEST_F(GCLTest, rigid_translation)
{
  if (bulk_.parallel_size() > 1) return;

  const std::string meshDims = "2x2x2|offset:0,65,0";
  const bool secondOrder = true;
  const double deltaT = 0.003; // approx 0.25 deg motion for given omega
  const std::string mesh_motion =
    "mesh_motion:                                                          \n"
    "  - name: interior                                                    \n"
    "    frame: non_inertial                                               \n"
    "    mesh_parts: [ block_1 ]                                           \n"
    "    motion:                                                           \n"
    "      - type: translation                                             \n"
    "        velocity: [1.0, 0.0, 0.0]                                     \n";

  fill_mesh_and_init_fields(meshDims);
  init_time_integrator(secondOrder, deltaT);
  register_algorithms(mesh_motion);
  init_states();
  compute_mesh_velocity();
  nodalGradAlgDriver_.execute();
  compute_div_mesh_vel();
  compute_dvoldt();
  compute_error();
}
