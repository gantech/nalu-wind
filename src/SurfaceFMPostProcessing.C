/*------------------------------------------------------------------------*/
/*  Copyright 2019 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <SurfaceFMPostProcessing.h>

#include <nalu_make_unique.h>
#include <Realm.h>
#include <PostProcessingData.h>
#include <SurfaceForceAndMomentAlgorithm.h>
#include <SurfaceForceAndMomentAlgorithmDriver.h>
#include <SurfaceForceAndMomentWallFunctionAlgorithm.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

// stk_topo
#include <stk_topology/topology.hpp>

namespace sierra {
namespace nalu {


SurfaceFMPostProcessing::SurfaceFMPostProcessing(
    Realm& realm)
    : realm_(realm)
{}

void SurfaceFMPostProcessing::register_surface_pp(
    const PostProcessingData &theData)
{
  stk::mesh::MetaData &meta = realm_.meta_data();

  stk::mesh::PartVector partVector;
  std::vector<std::string> targetNames = theData.targetNames_;
  for ( size_t in = 0; in < targetNames.size(); ++in) {
    stk::mesh::Part *targetPart = meta.get_part(targetNames[in]);
    if ( NULL == targetPart ) {
      NaluEnv::self().naluOutputP0() << "SurfacePP: can not find part with name: " << targetNames[in];
    }
    else {
      // found the part
      const std::vector<stk::mesh::Part*> & mesh_parts = targetPart->subsets();
      for( std::vector<stk::mesh::Part*>::const_iterator i = mesh_parts.begin();
          i != mesh_parts.end(); ++i )
      {
        stk::mesh::Part * const part = *i ;
        if ( !(meta.side_rank() == part->primary_entity_rank()) ) {
          NaluEnv::self().naluOutputP0() << "SurfacePP: part is not a face: " << targetNames[in];
        }
        partVector.push_back(part);
      }
    }
  }

  const std::string thePhysics = theData.physics_;

  if ( NULL == driverAlg_ )
      driverAlg_ = make_unique<SurfaceForceAndMomentAlgorithmDriver>(realm_);

  // register nodal fields in common
  VectorFieldType *pressureForce =  &(meta.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "pressure_force"));
  stk::mesh::put_field_on_mesh(*pressureForce, stk::mesh::selectUnion(partVector), meta.spatial_dimension(), nullptr);
  VectorFieldType *tauWall =  &(meta.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "tau_wall"));
  stk::mesh::put_field_on_mesh(*tauWall, stk::mesh::selectUnion(partVector), nullptr);
  ScalarFieldType *yplus =  &(meta.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "yplus"));
  stk::mesh::put_field_on_mesh(*yplus, stk::mesh::selectUnion(partVector), nullptr);
 
  // force output for these variables
  realm_.augment_output_variable_list(pressureForce->name());
  realm_.augment_output_variable_list(tauWall->name());
  realm_.augment_output_variable_list(yplus->name());

  if ( thePhysics == "surface_force_and_moment" ) {
    ScalarFieldType *assembledArea =  &(meta.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "assembled_area_force_moment"));
    stk::mesh::put_field_on_mesh(*assembledArea, stk::mesh::selectUnion(partVector), nullptr);
    SurfaceForceAndMomentAlgorithm *ppAlg
      = new SurfaceForceAndMomentAlgorithm(
          realm_, partVector, theData.outputFileName_, theData.frequency_,
          theData.parameters_, realm_.realmUsesEdges_);
    driverAlg_->algVec_.push_back(ppAlg);
  }
  else if ( thePhysics == "surface_force_and_moment_wall_function" ) {
    ScalarFieldType *assembledArea =  &(meta.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "assembled_area_force_moment_wf"));
    stk::mesh::put_field_on_mesh(*assembledArea, stk::mesh::selectUnion(partVector), nullptr);
    SurfaceForceAndMomentWallFunctionAlgorithm *ppAlg
      = new SurfaceForceAndMomentWallFunctionAlgorithm(
          realm_, partVector, theData.outputFileName_, theData.frequency_,
          theData.parameters_, realm_.realmUsesEdges_);
    driverAlg_->algVec_.push_back(ppAlg);
  }
  
}

void SurfaceFMPostProcessing::execute()
{
    driverAlg_->execute();
}


}
}
