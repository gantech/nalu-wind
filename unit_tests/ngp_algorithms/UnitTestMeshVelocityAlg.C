// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//


#include "kernels/UnitTestKernelUtils.h"
#include "UnitTestHelperObjects.h"

#include "AlgTraits.h"
#include "ngp_algorithms/GeometryInteriorAlg.h"
#include "ngp_algorithms/GeometryBoundaryAlg.h"
#include "ngp_algorithms/WallFuncGeometryAlg.h"
#include "ngp_algorithms/GeometryAlgDriver.h"
#include "ngp_algorithms/MeshVelocityAlg.h"
#include "utils/StkHelpers.h"

#include "mesh_motion/MotionRotation.h"

#include "UnitTestRealm.h"
#include "UnitTestUtils.h"

namespace {
    
    std::vector<double> transform(
        const sierra::nalu::MotionBase::TransMatType& transMat,
        const double* xyz )
    {
        std::vector<double> transCoord(3,0.0);

        // perform matrix multiplication between transformation matrix
        // and original coordinates to obtain transformed coordinates
        for (int d = 0; d < sierra::nalu::MotionBase::threeDVecSize; d++) {
            transCoord[d] = transMat[d][0]*xyz[0]
                +transMat[d][1]*xyz[1]
                +transMat[d][2]*xyz[2]
                +transMat[d][3];
        }

        return transCoord;
    }

    namespace hex8_golds {
        namespace ngp_mesh_velocity {
            static constexpr double swept_vol[12] = {
                -1.385467742974961e-17, -0.006239588540426781,
                1.385467742974961e-17, -0.00623958854042678,
                -1.385467742974961e-17, -0.018718765621280307,
                1.3854677429749644e-17, -0.0187187656212803,
                0.006239588540426745, 0.006239588540426742,
                0.018718765621280272, 0.018718765621280272,
            };

            static constexpr double face_vel_mag[12] = {
                -2.0782016144624415e-16, -0.06239588540426796,
                2.0782016144624415e-16, -0.06239588540426792,
                -2.0782016144624415e-16, -0.18718765621280323,
                2.0782016144624464e-16, -0.1871876562128031,
                0.062395885404267354, 0.06239588540426732,
                0.18718765621280262, 0.18718765621280262,
            };
        }
    }


}

TEST_F(TestKernelHex8Mesh, NGP_mesh_velocity)
{
  // Only execute for 1 processor runs
  if (bulk_.parallel_size() > 1) return;

  // declare relevant fields
  dnvField_ = &(meta_.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume", 3));
  stk::mesh::put_field_on_mesh(*dnvField_, meta_.universal_part(), nullptr);

  VectorFieldType* meshDisp_ = &(meta_.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "mesh_displacement", 3));
  stk::mesh::put_field_on_mesh(*meshDisp_, meta_.universal_part(), nullptr);

  VectorFieldType* cCoords_ = &(meta_.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "current_coordinates"));
  stk::mesh::put_field_on_mesh(*cCoords_, meta_.universal_part(), nullptr);

  const auto& meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(stk::topology::HEX_8);
  GenericFieldType* sweptVolume_ = &(meta_.declare_field<GenericFieldType>(stk::topology::ELEM_RANK, "swept_face_volume", 3));
  stk::mesh::put_field_on_mesh(*sweptVolume_, meta_.universal_part(), meSCS->num_integration_points(), nullptr);

  GenericFieldType* faceVelMag_ = &(meta_.declare_field<GenericFieldType>(stk::topology::ELEM_RANK, "face_velocity_mag", 2));
  stk::mesh::put_field_on_mesh(*faceVelMag_, meta_.universal_part(), meSCS->num_integration_points(), nullptr);

  fill_mesh_and_init_fields();

  unit_test_utils::HelperObjects helperObjs(
    bulk_, stk::topology::HEX_8, 1, partVec_[0]);

  sierra::nalu::TimeIntegrator timeIntegrator;
  timeIntegrator.timeStepN_ = 0.1;
  timeIntegrator.timeStepNm1_ = 0.1;
  timeIntegrator.currentTime_ = 0.2;
  timeIntegrator.gamma1_ = 1.5;
  timeIntegrator.gamma2_ = -2.0;
  timeIntegrator.gamma3_ = 0.5;
  helperObjs.realm.timeIntegrator_ = &timeIntegrator;
  // Force computation of edge area vector
  helperObjs.realm.realmUsesEdges_ = true;
  helperObjs.realm.solutionOptions_->meshMotion_ = true;
  
  sierra::nalu::GeometryAlgDriver geomAlgDriver(helperObjs.realm);
  geomAlgDriver.register_elem_algorithm<sierra::nalu::GeometryInteriorAlg>(
    sierra::nalu::INTERIOR, partVec_[0], "geometry");
  geomAlgDriver.execute();

  sierra::nalu::MeshVelocityAlg<sierra::nalu::AlgTraitsHex8> mvAlg(helperObjs.realm, partVec_[0]);

  //First set the mesh displacement corresponding to rotation about x-axis
  // create a yaml node describing rotation
  const std::string rotInfo =
      "omega: 1.0              \n"
      "axis: [1.0,0.0,0.0]     \n"
      "centroid: [0.0,0.0,0.0] \n"
      ;
  YAML::Node rotNode = YAML::Load(rotInfo);
  sierra::nalu::MotionRotation rotClass(rotNode);
  VectorFieldType *meshDispNp1 = &(meshDisp_->field_of_state(stk::mesh::StateNP1));
  VectorFieldType *meshDispN = &(meshDisp_->field_of_state(stk::mesh::StateN));
  VectorFieldType *meshDispNm1 = &(meshDisp_->field_of_state(stk::mesh::StateNM1));
  {
      stk::mesh::Selector sel = meta_.universal_part();
      const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);
      for (const auto* b: bkts) {
          for (const auto node: *b) {
              double* dispNp1 = stk::mesh::field_data(*meshDispNp1, node);
              double* dispN = stk::mesh::field_data(*meshDispN, node);
              double* dispNm1 = stk::mesh::field_data(*meshDispNm1, node);
              double* ccoord = stk::mesh::field_data(*cCoords_, node);
              double* mcoord = stk::mesh::field_data(*coordinates_, node);
              dispNm1[0] = 0.0;
              dispNm1[1] = 0.0;
              dispNm1[2] = 0.0;

              rotClass.build_transformation(0.1,mcoord);
              std::vector<double> rot_xyz = transform(rotClass.get_trans_mat(), mcoord);
              dispN[0] = rot_xyz[0] - mcoord[0];
              dispN[1] = rot_xyz[1] - mcoord[1];
              dispN[2] = rot_xyz[2] - mcoord[2];
              
              rotClass.build_transformation(0.2,mcoord);
              rot_xyz = transform(rotClass.get_trans_mat(), mcoord);
              dispNp1[0] = rot_xyz[0] - mcoord[0];
              dispNp1[1] = rot_xyz[1] - mcoord[1];
              dispNp1[2] = rot_xyz[2] - mcoord[2];
              
              ccoord[0] = rot_xyz[0];
              ccoord[1] = rot_xyz[1];
              ccoord[2] = rot_xyz[2];
          }
      }
  }

  geomAlgDriver.execute();

  GenericFieldType *sweptVolN = &(sweptVolume_->field_of_state(stk::mesh::StateN));
  
  stk::mesh::Selector sel = meta_.universal_part();  
  {
      const auto& bkts = bulk_.get_buckets(stk::topology::ELEM_RANK, sel);
      for (const auto* b: bkts) {
          double *sv = stk::mesh::field_data(*sweptVolN, *b, 0);
          sv[0] = 0.0 ;
          sv[1] = -0.006239588540426752 ;
          sv[2] = 0.0 ;
          sv[3] = -0.006239588540426754 ;
          sv[4] = 0.0 ;
          sv[5] = -0.018718765621280272 ;
          sv[6] = 0.0 ;
          sv[7] = -0.018718765621280272 ;
          sv[8] = 0.006239588540426764 ;
          sv[9] = 0.006239588540426764 ;
          sv[10] = 0.018718765621280286 ;
          sv[11] = 0.018718765621280286 ;
      }
  }
  
  mvAlg.execute();
  
  const double tol = 1.0e-16;
  namespace gold_values = ::hex8_golds::ngp_mesh_velocity;
  {
    const auto& bkts = bulk_.get_buckets(stk::topology::ELEM_RANK, sel);
    int counter=0;
    for (const auto* b: bkts) {
        const double *sv = stk::mesh::field_data(*sweptVolume_, *b, 0);
        const double *fvm = stk::mesh::field_data(*faceVelMag_, *b, 0);        
        for (int i=0; i < 12; i++) {
            EXPECT_NEAR(gold_values::swept_vol[i], sv[i], tol);
            counter++;
            EXPECT_NEAR(gold_values::face_vel_mag[i], fvm[i], tol);
            counter++;
        }
    }
    EXPECT_EQ(counter,24);
  }
  
}
