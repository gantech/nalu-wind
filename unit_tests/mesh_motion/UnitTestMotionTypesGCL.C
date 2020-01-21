#include <gtest/gtest.h>
#include <limits>

#include "mesh_motion/MeshMotionAlg.h"
#include "mesh_motion/MotionRotation.h"
#include "mesh_motion/MotionScaling.h"
#include "mesh_motion/MotionTranslation.h"
#include "mesh_motion/MotionWaves.h"

#include "ComputeGeometryInteriorAlgorithm.h"
#include "ComputeGeometryBoundaryAlgorithm.h"
#include "Realm.h"
#include "SolutionOptions.h"
#include "TimeIntegrator.h"
#include "utils/ComputeVectorDivergence.h"

#include <stk_mesh/base/FieldParallel.hpp>

#include "UnitTestRealm.h"
#include "UnitTestUtils.h"

#include <string>

namespace {
  // create a yaml node describing translation
  const std::string mInfo =
    "mesh_motion:																	 		 \n"
    "  																						 		 \n"
    "  - name: waterwaves 										     		 \n"
    "    mesh_parts: [ block_1 ]											 \n"
    "    frame: non_inertial									     		 \n"
    "    motion:																	 		 \n"
 		"      - type: water_waves										 		 \n"
    "        wave_motion_model: Sinusoidal_full_domain \n"
    "        amplitude: 0.0159 												 \n"
    "        waveperiod: 8. 													 \n"
    "        wavelength: 1.0 													 \n"
    ;
  
  const YAML::Node meshMotionNode = YAML::Load(mInfo);
  const YAML::Node frames				  = meshMotionNode["mesh_motion"];
   
  // Tolerance
  const double testTol = 1e-12;
}

TEST(utils, test_mesh_motion_gcl)
{
  
  // create realm
  unit_test_utils::NaluTest naluObj;
  sierra::nalu::Realm& realm = naluObj.create_realm();  
  realm.solutionOptions_->meshMotion_ =true;

  sierra::nalu::TimeIntegrator timeIntegrator;
  timeIntegrator.secondOrderTimeAccurate_ = false;
  realm.timeIntegrator_ = &timeIntegrator;
  // register mesh motion fields and initialize coordinate fields
  realm.register_nodal_fields( &(realm.meta_data().universal_part()) );
  realm.init_current_coordinates();


  // create mesh motion algorithm class
  const YAML::Node meshMotionNode = YAML::Load(mInfo);
  std::unique_ptr<sierra::nalu::MeshMotionAlg> meshMotionAlg;
  meshMotionAlg.reset(
    new sierra::nalu::MeshMotionAlg( realm.bulk_data(), meshMotionNode["mesh_motion"] ));

  // get relevant fields
  VectorFieldType* modelCoords = realm.meta_data().get_field<VectorFieldType>(
    stk::topology::NODE_RANK, "coordinates");
  VectorFieldType* currCoords = realm.meta_data().get_field<VectorFieldType>(
    stk::topology::NODE_RANK, "current_coordinates");
  VectorFieldType* meshVelocity = realm.meta_data().get_field<VectorFieldType>(
    stk::topology::NODE_RANK, "mesh_velocity");
  // define auxillary variable to calculate the divergence of the vector field
  ScalarFieldType *div_mesh_velocity = &(realm.meta_data().declare_field<ScalarFieldType>(
		stk::topology::NODE_RANK, "div_mesh_vector"));
  stk::mesh::put_field_on_mesh(*div_mesh_velocity, realm.meta_data().universal_part(), nullptr);
  ScalarFieldType *duaNdlVol = &(realm.meta_data().declare_field<ScalarFieldType>(
		stk::topology::NODE_RANK, "dual_nodal_volume"));
  stk::mesh::put_field_on_mesh(*duaNdlVol, realm.meta_data().universal_part(), nullptr);
  ScalarFieldType& elemVol = realm.meta_data().declare_field<ScalarFieldType>(
    stk::topology::ELEMENT_RANK, "element_volume");  
  stk::mesh::put_field_on_mesh(elemVol, realm.meta_data().universal_part(), nullptr);
  
  // create mesh
  const std::string meshSpec("generated:144x4x1|bbox:0,0,0,6,5,1|sideset:xXyYzZ|show");
  unit_test_utils::fill_hex8_mesh(meshSpec, realm.bulk_data());

  const sierra::nalu::MasterElement* meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(stk::topology::QUAD_4);
  const int numScsIp = meFC->num_integration_points();
  GenericFieldType *exposedAreaVec = &(realm.meta_data().declare_field<GenericFieldType>(realm.meta_data().side_rank(), "exposed_area_vector"));

  // get the parts in the current motion frame
  stk::mesh::Selector sel = stk::mesh::Selector(realm.meta_data().universal_part())
    & (realm.meta_data().locally_owned_part() | realm.meta_data().globally_shared_part());
  const auto& bkts = realm.bulk_data().get_buckets(stk::topology::NODE_RANK, sel);

  /////////////////////////////////////////////////////////////
  // initialize and execute mesh motion algorithm
  const double currTime = 0.0;
  meshMotionAlg->initialize(currTime);
  // Compute divergence of the mesh velocity vector here
  // create dual volumes
  sierra::nalu::ComputeGeometryInteriorAlgorithm geomAlg(realm, &(realm.meta_data().universal_part()));
  geomAlg.execute();
  stk::mesh::parallel_sum(realm.bulk_data(), {duaNdlVol});

  sierra::nalu::ComputeGeometryBoundaryAlgorithm bndyGeomAlg(realm, realm.meta_data().get_part("surface_1"));
  bndyGeomAlg.execute();
   
  stk::mesh::PartVector partVec;
  partVec.push_back( &(realm.meta_data().universal_part()) );

  stk::mesh::PartVector bndyPartVec;
  bndyPartVec.push_back( realm.meta_data().get_part("surface_1") );

  sierra::nalu::compute_vector_divergence( realm.bulk_data(),
                                           partVec, bndyPartVec,
                                           meshVelocity, div_mesh_velocity );
  // check values
  for (auto b: bkts) {
    for (size_t in=0; in < b->size(); in++) {

      auto node = (*b)[in]; // mesh node and NOT YAML node

      EXPECT_NEAR(0., 0., testTol);

    } // end for loop - in index
  } // end for loop - bkts
   
  
}

