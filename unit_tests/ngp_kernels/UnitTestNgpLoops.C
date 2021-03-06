/*------------------------------------------------------------------------*/
/*  Copyright 2019 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "gtest/gtest.h"
#include "UnitTestUtils.h"

#include "ngp_utils/NgpLoopUtils.h"
#include "ngp_utils/NgpFieldOps.h"

#include <cmath>

class NgpLoopTest : public ::testing::Test
{
public:
  NgpLoopTest()
    : meta(3),
      bulk(meta, MPI_COMM_WORLD),
      partVec(),
      density(&meta.declare_field<ScalarFieldType>(
                 stk::topology::NODE_RANK, "density")),
      pressure(&meta.declare_field<ScalarFieldType>(
                 stk::topology::NODE_RANK, "pressure")),
      velocity(&meta.declare_field<VectorFieldType>(
                 stk::topology::NODE_RANK, "velocity")),
      mdotEdge(&meta.declare_field<ScalarFieldType>(
                 stk::topology::EDGE_RANK, "mass_flow_rate")),
      massFlowRate(&meta.declare_field<GenericFieldType>(
                     stk::topology::ELEM_RANK, "mass_flow_rate_scs"))
  {
    const double ten = 10.0;
    const double zero = 0.0;
    const double oneVec[3] = {1.0, 1.0, 1.0};
    sierra::nalu::HexSCS hex8SCS;
    stk::mesh::put_field_on_mesh(*density, meta.universal_part(), 1, &ten);
    stk::mesh::put_field_on_mesh(*pressure, meta.universal_part(), 1, &zero);
    stk::mesh::put_field_on_mesh(*velocity, meta.universal_part(), 3, oneVec);
    stk::mesh::put_field_on_mesh(
      *massFlowRate, meta.universal_part(), hex8SCS.num_integration_points(),
      &zero);
    stk::mesh::put_field_on_mesh(*mdotEdge, meta.universal_part(), 1, &zero);
  }

  ~NgpLoopTest() = default;

  void fill_mesh_and_init_fields(const std::string& meshSpec = "generated:2x2x2")
  {
    unit_test_utils::fill_hex8_mesh(meshSpec, bulk);
    partVec = { meta.get_part("block_1")};

    coordField = static_cast<const VectorFieldType*>(meta.coordinate_field());
    EXPECT_TRUE(coordField != nullptr);
  }

  stk::mesh::MetaData meta;
  stk::mesh::BulkData bulk;
  stk::mesh::PartVector partVec;
  const VectorFieldType* coordField{nullptr};
  ScalarFieldType* density{nullptr};
  ScalarFieldType* pressure{nullptr};
  VectorFieldType* velocity{nullptr};
  ScalarFieldType* mdotEdge{nullptr};
  GenericFieldType* massFlowRate{nullptr};
};

void basic_node_loop(
  const stk::mesh::BulkData& bulk,
  ScalarFieldType& pressure)
{
  using Traits = sierra::nalu::nalu_ngp::NGPMeshTraits<ngp::Mesh>;
  const double presSet = 4.0;

  const auto& meta = bulk.mesh_meta_data();
  stk::mesh::Selector sel = meta.universal_part();
  ngp::Mesh ngpMesh(bulk);
  ngp::Field<double> ngpPressure(bulk, pressure);

  sierra::nalu::nalu_ngp::run_entity_algorithm(
    ngpMesh, stk::topology::NODE_RANK, sel,
    KOKKOS_LAMBDA(const typename Traits::MeshIndex& meshIdx) {
      ngpPressure.get(meshIdx, 0) = presSet;
    });

  ngpPressure.modify_on_device();
  ngpPressure.sync_to_host();

  // Do checks
  {
    const auto& bkts = bulk.get_buckets(stk::topology::NODE_RANK, sel);
    const double tol = 1.0e-16;
    for (const auto* b: bkts) {
      for (const auto node: *b) {
        const double* pres = stk::mesh::field_data(pressure, node);
        EXPECT_NEAR(presSet, pres[0], tol);
      }
    }
  }
}

void basic_elem_loop(
  const stk::mesh::BulkData& bulk,
  ScalarFieldType& pressure,
  GenericFieldType& massFlowRate)
{
  const double flowRate = 4.0;
  const double presSet = 10.0;

  ngp::Mesh ngpMesh(bulk);
  ngp::Field<double> ngpMassFlowRate(bulk, massFlowRate);
  ngp::Field<double> ngpPressure(bulk, pressure);

  const auto& meta = bulk.mesh_meta_data();
  stk::mesh::Selector sel = meta.universal_part();

  sierra::nalu::nalu_ngp::run_elem_algorithm(
    ngpMesh, stk::topology::ELEMENT_RANK, sel,
    KOKKOS_LAMBDA(const sierra::nalu::nalu_ngp::EntityInfo<ngp::Mesh>& einfo) {
      ngpMassFlowRate.get(einfo.meshIdx, 0) = flowRate;

      const auto& nodes = einfo.entityNodes;
      const int numNodes = nodes.size();
      for (int i=0; i < numNodes; ++i)
        ngpPressure.get(ngpMesh, nodes[i], 0) = presSet;
    });

  ngpMassFlowRate.modify_on_device();
  ngpMassFlowRate.sync_to_host();
  ngpPressure.modify_on_device();
  ngpPressure.sync_to_host();


  {
    const auto& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK,  sel);
    const double tol = 1.0e-16;
    for (const stk::mesh::Bucket* b : elemBuckets)
    {
      for (stk::mesh::Entity elem : *b)
      {
        const double* flowRateData = stk::mesh::field_data(massFlowRate, elem);
        EXPECT_NEAR(flowRate, *flowRateData, tol);

        const stk::mesh::Entity* nodes = bulk.begin_nodes(elem);
        const unsigned numNodes = bulk.num_nodes(elem);
        for (unsigned n = 0; n < numNodes; ++n)
        {
          const double* pres = stk::mesh::field_data(pressure, nodes[n]);
          EXPECT_NEAR(presSet, pres[0], tol);
        }
      }
    }
  }
}

void basic_edge_loop(
  const stk::mesh::BulkData& bulk,
  ScalarFieldType& pressure,
  ScalarFieldType& mdotEdge)
{
  const double flowRate = 4.0;
  const double presSet = 10.0;

  ngp::Mesh ngpMesh(bulk);
  ngp::Field<double> ngpMassFlowRate(bulk, mdotEdge);
  ngp::Field<double> ngpPressure(bulk, pressure);

  const auto& meta = bulk.mesh_meta_data();
  stk::mesh::Selector sel = meta.universal_part();

  sierra::nalu::nalu_ngp::run_edge_algorithm(
    ngpMesh, sel,
    KOKKOS_LAMBDA(const sierra::nalu::nalu_ngp::EntityInfo<ngp::Mesh>& einfo) {
      ngpMassFlowRate.get(einfo.meshIdx, 0) = flowRate;

      const auto& nodes = einfo.entityNodes;
      const int numNodes = nodes.size();
      for (int i=0; i < numNodes; ++i)
        ngpPressure.get(ngpMesh, nodes[i], 0) = presSet;
    });

  {
    const auto& edgeBuckets = bulk.get_buckets(stk::topology::EDGE_RANK,  sel);
    const double tol = 1.0e-16;
    for (const stk::mesh::Bucket* b : edgeBuckets)
    {
      for (stk::mesh::Entity edge : *b)
      {
        const double* flowRateData = stk::mesh::field_data(mdotEdge, edge);
        EXPECT_NEAR(flowRate, *flowRateData, tol);

        const stk::mesh::Entity* nodes = bulk.begin_nodes(edge);
        const unsigned numNodes = bulk.num_nodes(edge);
        for (unsigned n = 0; n < numNodes; ++n)
        {
          const double* pres = stk::mesh::field_data(pressure, nodes[n]);
          EXPECT_NEAR(presSet, pres[0], tol);
        }
      }
    }
  }
}

void elem_loop_scratch_views(
  const stk::mesh::BulkData& bulk,
  ScalarFieldType& pressure,
  VectorFieldType& velocity)
{
  using Traits = sierra::nalu::nalu_ngp::NGPMeshTraits<ngp::Mesh>;
  using Hex8Traits = sierra::nalu::AlgTraitsHex8;
  using ElemSimdData = sierra::nalu::nalu_ngp::ElemSimdData<ngp::Mesh>;
  typedef Kokkos::DualView<double*, Kokkos::LayoutRight, sierra::nalu::DeviceSpace> DoubleTypeView;

  const auto& meta = bulk.mesh_meta_data();
  sierra::nalu::ElemDataRequests dataReq(meta);
  auto meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element<
    sierra::nalu::AlgTraitsHex8>();
  dataReq.add_cvfem_volume_me(meSCV);

  auto* coordsField = bulk.mesh_meta_data().coordinate_field();
  dataReq.add_coordinates_field(*coordsField, 3, sierra::nalu::CURRENT_COORDINATES);
  dataReq.add_gathered_nodal_field(velocity, 3);
  dataReq.add_gathered_nodal_field(pressure, 1);
  dataReq.add_master_element_call(
    sierra::nalu::SCV_VOLUME, sierra::nalu::CURRENT_COORDINATES);
  dataReq.add_master_element_call(
    sierra::nalu::SCV_SHIFTED_SHAPE_FCN, sierra::nalu::CURRENT_COORDINATES);

  sierra::nalu::nalu_ngp::MeshInfo<> meshInfo(bulk);
  stk::mesh::Selector sel = meta.universal_part();

  const unsigned velID = velocity.mesh_meta_data_ordinal();
  const unsigned presID = pressure.mesh_meta_data_ordinal();

  // Field updates
  Traits::DblType xVel = 1.0;
  Traits::DblType yVel = 2.0;
  Traits::DblType zVel = 3.0;

  const auto ngpMesh = meshInfo.ngp_mesh();
  const auto& fieldMgr = meshInfo.ngp_field_manager();
  auto ngpVel = fieldMgr.get_field<double>(velID);

  const int numNodes = 8;
  DoubleTypeView volCheck("scv_volume", numNodes);
  Kokkos::deep_copy(volCheck.h_view, 0.0);
  volCheck.template modify<typename DoubleTypeView::host_mirror_space>();
  volCheck.template sync<typename DoubleTypeView::execution_space>();

  sierra::nalu::nalu_ngp::run_elem_algorithm(
    meshInfo, stk::topology::ELEM_RANK, dataReq, sel,
    KOKKOS_LAMBDA(ElemSimdData & edata) {
      const auto ngpVelOp = sierra::nalu::nalu_ngp::simd_nodal_field_updater(
        ngpMesh, ngpVel, edata);

      Traits::DblType test = 0.0;
      auto& scrViews = edata.simdScrView;
      auto& v_pres = scrViews.get_scratch_view_1D(presID);
      auto& v_vel = scrViews.get_scratch_view_2D(velID);
      auto& scv_vol =
        scrViews.get_me_views(sierra::nalu::CURRENT_COORDINATES).scv_volume;

      test += v_vel(0, 0) + v_pres(0) * scv_vol(0);

      for (int i=0; i < numNodes; ++i) {
        volCheck.d_view(i) = stk::simd::get_data(scv_vol(i), 0);

        // Scatter SIMD value to nodes
        ngpVelOp(i, 0) = xVel;
        ngpVelOp(i, 1) = yVel;
        ngpVelOp(i, 2) = zVel;
      }
    });

  ngpVel.modify_on_device();
  ngpVel.sync_to_host();

  volCheck.modify<DoubleTypeView::execution_space>();
  volCheck.sync<DoubleTypeView::host_mirror_space>();

  for (int i=0; i < numNodes; ++i)
    EXPECT_NEAR(volCheck.h_view(i), 0.125, 1.0e-12);

  {
    const double xVel = 1.0;
    const double yVel = 2.0;
    const double zVel = 3.0;
    const auto& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK,  sel);
    const double tol = 1.0e-16;
    for (const stk::mesh::Bucket* b : elemBuckets)
    {
      for (stk::mesh::Entity elem : *b) {
        const stk::mesh::Entity* nodes = bulk.begin_nodes(elem);
        for (int i=0; i < Hex8Traits::nodesPerElement_; ++i) {
          const double* velptr = stk::mesh::field_data(velocity, nodes[i]);
          EXPECT_NEAR(velptr[0], xVel, tol);
          EXPECT_NEAR(velptr[1], yVel, tol);
          EXPECT_NEAR(velptr[2], zVel, tol);
        }
      }
    }
  }
}
void calc_mdot_elem_loop(
  const stk::mesh::BulkData& bulk,
  ScalarFieldType& density,
  VectorFieldType& velocity,
  GenericFieldType& massFlowRate)
{
  using Traits = sierra::nalu::nalu_ngp::NGPMeshTraits<ngp::Mesh>;
  using Hex8Traits = sierra::nalu::AlgTraitsHex8;
  using ElemSimdData = sierra::nalu::nalu_ngp::ElemSimdData<ngp::Mesh>;

  const auto& meta = bulk.mesh_meta_data();
  sierra::nalu::ElemDataRequests dataReq(meta);
  auto meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element<Hex8Traits>();
  dataReq.add_cvfem_surface_me(meSCS);

  auto* coordsField = bulk.mesh_meta_data().coordinate_field();
  dataReq.add_coordinates_field(*coordsField, 3, sierra::nalu::CURRENT_COORDINATES);
  dataReq.add_gathered_nodal_field(velocity, 3);
  dataReq.add_gathered_nodal_field(density, 1);
  dataReq.add_master_element_call(
    sierra::nalu::SCS_AREAV, sierra::nalu::CURRENT_COORDINATES);
  dataReq.add_master_element_call(
    sierra::nalu::SCS_SHIFTED_SHAPE_FCN, sierra::nalu::CURRENT_COORDINATES);

  sierra::nalu::nalu_ngp::MeshInfo<> meshInfo(bulk);
  stk::mesh::Selector sel = meta.universal_part();

  const unsigned velID = velocity.mesh_meta_data_ordinal();
  const unsigned rhoID = density.mesh_meta_data_ordinal();
  const auto mdotID = massFlowRate.mesh_meta_data_ordinal();
  const auto ngpMesh = meshInfo.ngp_mesh();
  const auto& fieldMgr = meshInfo.ngp_field_manager();
  ngp::Field<double> ngpMdot = fieldMgr.get_field<double>(mdotID);

  sierra::nalu::nalu_ngp::run_elem_algorithm(
    meshInfo, stk::topology::ELEM_RANK, dataReq, sel,
    KOKKOS_LAMBDA(ElemSimdData& edata) {
      // SIMD Element field operation handler
      const auto mdotOps = sierra::nalu::nalu_ngp::simd_elem_field_updater(
        ngpMesh, ngpMdot, edata);

      NALU_ALIGNED Traits::DblType rhoU[Hex8Traits::nDim_];
      auto& scrViews = edata.simdScrView;
      auto& v_rho = scrViews.get_scratch_view_1D(rhoID);
      auto& v_vel = scrViews.get_scratch_view_2D(velID);
      auto& meViews = scrViews.get_me_views(sierra::nalu::CURRENT_COORDINATES);
      auto& v_area = meViews.scs_areav;
      auto& v_shape_fcn = meViews.scs_shifted_shape_fcn;

      for (int ip = 0; ip < Hex8Traits::numScsIp_; ++ip) {
        for (int d=0; d < Hex8Traits::nDim_; ++d)
          rhoU[d] = 0.0;

        for (int ic=0; ic < Hex8Traits::nodesPerElement_; ++ic) {
          const auto r = v_shape_fcn(ip, ic);
          for (int d=0; d < Hex8Traits::nDim_; ++d)
            rhoU[d] += r * v_rho(ic) * v_vel(ic, d);
        }

        Traits::DblType tmdot = 0.0;
        for (int d=0; d < Hex8Traits::nDim_; ++d)
          tmdot += rhoU[d] * v_area(ip, d);

        mdotOps(ip) = tmdot;
      }
    });

  ngpMdot.modify_on_device();
  ngpMdot.sync_to_host();

  {
    const double flowRate = 2.5;
    const auto& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK,  sel);
    const double tol = 1.0e-16;
    for (const stk::mesh::Bucket* b : elemBuckets)
    {
      for (stk::mesh::Entity elem : *b)
      {
        const double* flowRateData = stk::mesh::field_data(massFlowRate, elem);
        for (int i=0; i < Hex8Traits::numScsIp_; ++i)
          EXPECT_NEAR(flowRate, std::abs(flowRateData[i]), tol);
      }
    }
  }
}

TEST_F(NgpLoopTest, NGP_basic_node_loop)
{
  fill_mesh_and_init_fields("generated:2x2x2");

  basic_node_loop(bulk, *pressure);
}

TEST_F(NgpLoopTest, NGP_basic_elem_loop)
{
  fill_mesh_and_init_fields("generated:2x2x2");

  basic_elem_loop(bulk, *pressure, *massFlowRate);
}

TEST_F(NgpLoopTest, NGP_basic_edge_loop)
{
  fill_mesh_and_init_fields("generated:2x2x2");

  basic_edge_loop(bulk, *pressure, *mdotEdge);
}

TEST_F(NgpLoopTest, NGP_elem_loop_scratch_views)
{
  fill_mesh_and_init_fields("generated:2x2x2");

  elem_loop_scratch_views(bulk, *pressure, *velocity);
}

TEST_F(NgpLoopTest, NGP_calc_mdot_elem_loop)
{
  fill_mesh_and_init_fields("generated:2x2x2");

  calc_mdot_elem_loop(bulk, *density, *velocity, *massFlowRate);
}
