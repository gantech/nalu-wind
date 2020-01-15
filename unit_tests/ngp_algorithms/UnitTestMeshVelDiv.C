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
#include "utils/ComputeVectorDivergence.h"

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

}

TEST_F(MeshVelocityKernelHex8Mesh, NGP_mesh_vel_div)
{
  // Only execute for 1 processor runs
  if (bulk_.parallel_size() > 1) return;

  MeshVelocityKernelHex8Mesh::fill_mesh_and_init_fields("generated:2x2x2",true,false);

  unit_test_utils::HelperObjects helperObjs(
    bulk_, stk::topology::HEX_8, 1, partVec_[0]);

  sierra::nalu::TimeIntegrator timeIntegrator;
  timeIntegrator.timeStepN_ = 0.1;
  timeIntegrator.timeStepNm1_ = 0.1;
  timeIntegrator.currentTime_ = 0.2;
  timeIntegrator.gamma1_ = 1.0;
  timeIntegrator.gamma2_ = -1.0;
  timeIntegrator.gamma3_ = 0.0;
  helperObjs.realm.timeIntegrator_ = &timeIntegrator;
  // Force computation of edge area vector
  helperObjs.realm.realmUsesEdges_ = true;
  helperObjs.realm.solutionOptions_->meshMotion_ = true;
  
  sierra::nalu::GeometryAlgDriver geomAlgDriver(helperObjs.realm);
  geomAlgDriver.register_elem_algorithm<sierra::nalu::GeometryInteriorAlg>(
    sierra::nalu::INTERIOR, partVec_[0], "geometry");

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
  mvAlg.execute();
  stk::mesh::PartVector bndyPartVec;
  sierra::nalu::compute_scalar_divergence(bulk_, partVec_, bndyPartVec, faceVelMag_, divMeshVelField_);

  const double tol = 1.0e-14;
  double full_dnv_mdv = 1e15;
  bool foundNode = false;
  {
    stk::mesh::Selector sel = meta_.universal_part();
    const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);
    for (const auto* b: bkts) {
      for (const auto node: *b) {
        const double* meshVelDiv = stk::mesh::field_data(*divMeshVelField_,node);
        if (bulk_.num_elements(node) == 8) {
          full_dnv_mdv = *meshVelDiv;
          foundNode = true;
          break;
        }
        if (foundNode)
          break;
      }
      if (foundNode)
        break;
    }
  }
  EXPECT_NEAR(0.0, full_dnv_mdv, tol);
  
}
