/*------------------------------------------------------------------------*/
/*  Copyright 2019 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <SurfaceFMPostProcessing.h>

#include <nalu_make_unique.h>
#include <FieldFunctions.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <master_element/MasterElement.h>
#include <NaluEnv.h>
#include <PostProcessingData.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// stk_topo
#include <stk_topology/topology.hpp>

// stk util
#include <stk_util/parallel/ParallelReduce.hpp>

// basic c++
#include <fstream>

namespace sierra {
namespace nalu {


SurfaceFMPostProcessing::SurfaceFMPostProcessing(
    Realm& realm)
    : realm_(realm),
      yplusCrit_(11.63),
      elog_(9.8),
      kappa_(realm.get_turb_model_constant(TM_kappa)),
      pressureForce_(NULL),
      tauWall_(NULL),
      yplus_(NULL),
      assembledArea_(NULL),
      velocity_(NULL),
      bcVelocity_(NULL),
      wallFrictionVelocityBip_(NULL),
      wallNormalDistanceBip_(NULL)
{
    stk::mesh::MetaData & meta = realm_.meta_data();
    pressure_ = meta.get_field<ScalarFieldType>(
        stk::topology::NODE_RANK, "pressure");
    coordinates_ = meta.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, realm_.get_coordinates_name());
    density_ = meta.get_field<ScalarFieldType>(
        stk::topology::NODE_RANK, "density");
    // extract viscosity name
    const std::string viscName = realm_.is_turbulent()
        ? "effective_viscosity_u" : "viscosity";
    viscosity_ = meta.get_field<ScalarFieldType>(
        stk::topology::NODE_RANK, viscName);
    dudx_ = meta.get_field<GenericFieldType>(
        stk::topology::NODE_RANK, "dudx");
    exposedAreaVec_ = meta.get_field<GenericFieldType>(
        meta.side_rank(), "exposed_area_vector");
    velocity_ = meta.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "velocity");
    
}

void SurfaceFMPostProcessing::register_surface_pp(
    const PostProcessingData &theData)
{
  stk::mesh::MetaData &meta = realm_.meta_data();

  stk::mesh::PartVector partVector;
  std::vector<std::string> targetNames = theData.targetNames_;
  for ( size_t in = 0; in < targetNames.size(); ++in) {
    stk::mesh::Part *targetPart = meta.get_part(targetNames[in]);
    if ( NULL == targetPart ) {
      NaluEnv::self().naluOutputP0() <<
          "SurfacePP: can not find part with name: " << targetNames[in];
    }
    else {
      // found the part
      auto& mesh_parts = targetPart->subsets();
      for( auto part: mesh_parts )
      {
        if ( !(meta.side_rank() == part->primary_entity_rank()) ) {
          NaluEnv::self().naluOutputP0() << "SurfacePP: part is not a face: "
                                         << targetNames[in];
        }
        partVector.push_back(part);
        allPartVector_.push_back(part);
      }
    }
  }

  const std::string thePhysics = theData.physics_;

  // register nodal fields in common
  pressureForce_ = &(meta.declare_field<VectorFieldType>(
                         stk::topology::NODE_RANK, "pressure_force"));
  stk::mesh::put_field_on_mesh(*pressureForce_,
                               stk::mesh::selectUnion(partVector),
                               meta.spatial_dimension(), nullptr);
  tauWall_ =  &(meta.declare_field<VectorFieldType>(
                    stk::topology::NODE_RANK, "tau_wall"));
  stk::mesh::put_field_on_mesh(*tauWall_,
                               stk::mesh::selectUnion(partVector), nullptr);
  yplus_ =  &(meta.declare_field<ScalarFieldType>(
                  stk::topology::NODE_RANK, "yplus"));
  stk::mesh::put_field_on_mesh(*yplus_,
                               stk::mesh::selectUnion(partVector), nullptr);
  assembledArea_ =  &(meta.declare_field<ScalarFieldType>(
                          stk::topology::NODE_RANK,
                          "assembled_area_force_moment"));
  stk::mesh::put_field_on_mesh(*assembledArea_,
                               stk::mesh::selectUnion(partVector), nullptr);
 
  // force output for these variables
  realm_.augment_output_variable_list(pressureForce_->name());
  realm_.augment_output_variable_list(tauWall_->name());
  realm_.augment_output_variable_list(yplus_->name());

  std::array<double, 3> centroidCoords{ {
              theData.parameters_[0],
              theData.parameters_[1],
              theData.parameters_[2]} };
  
  if ( thePhysics == "surface_force_and_moment" ) {
    SurfaceFMData new_sfm_data = {partVector,
                                  theData.outputFileName_,
                                  centroidCoords,
                                  theData.frequency_,
                                  false};
    surfaceFMData_.push_back(new_sfm_data);
  }
  else if ( thePhysics == "surface_force_and_moment_wall_function" ) {
    SurfaceFMData new_sfm_data = {partVector,
                                  theData.outputFileName_,
                                  centroidCoords,
                                  theData.frequency_,
                                  true};
    surfaceFMData_.push_back(new_sfm_data);
    bcVelocity_ = meta.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "wall_velocity_bc");
    wallFrictionVelocityBip_ = meta.get_field<GenericFieldType>(
        meta.side_rank(), "wall_friction_velocity_bip");
    wallNormalDistanceBip_ = meta.get_field<GenericFieldType>(
        meta.side_rank(), "wall_normal_distance_bip");
  }

  create_file(theData.outputFileName_);
  
}


void SurfaceFMPostProcessing::create_file(
    std::string fileName)
{

    // deal with file name and banner
    if ( NaluEnv::self().parallel_rank() == 0 ) {
        std::ofstream myfile;
        myfile.open(fileName.c_str());
        myfile << std::setw(16) 
               << "Time" << std::setw(16) 
               << "Fpx"  << std::setw(16)
               << "Fpy" << std::setw(16)
               << "Fpz" << std::setw(16) 
               << "Fvx"  << std::setw(16)
               << "Fvy" << std::setw(16)
               << "Fvz" << std::setw(16) 
               << "Mtx"  << std::setw(16)
               << "Mty" << std::setw(16)
               << "Mtz" << std::setw(16) 
               << "Y+min" << std::setw(16)
               << "Y+max"<< std::endl;
        myfile.close();
    }

}
void SurfaceFMPostProcessing::execute()
{
    // zero fields
    zero_fields();

    if (realm_.has_mesh_motion())
        calc_assembled_area(allPartVector_);

    for (auto& sfm_data: surfaceFMData_) {
        if(sfm_data.wallFunction)
            calc_surface_force_wallfn(sfm_data);
        else
            calc_surface_force(sfm_data);

    }
    
    // parallel assembly
    parallel_assemble_fields();
    
}

//======================
// assemble area
//======================
void SurfaceFMPostProcessing::calc_assembled_area(
    stk::mesh::PartVector partvec)
{
  // common
  stk::mesh::BulkData & bulk = realm_.bulk_data();
  stk::mesh::MetaData & meta = realm_.meta_data();
  const int nDim = meta.spatial_dimension();

  // define some common selectors
  stk::mesh::Selector sel = meta.locally_owned_part()
      &stk::mesh::selectUnion(partvec);

  auto& face_buckets = realm_.get_buckets( meta.side_rank(), sel );
  for ( auto b: face_buckets ) {
      // face master element
    MasterElement *meFC =
        MasterElementRepo::get_surface_master_element(b->topology());
    const int numScsBip = meFC->numIntPoints_;

    // mapping from ip to nodes for this ordinal;
    // face perspective (use with face_node_relations)
    const int *faceIpNodeMap = meFC->ipNodeMap();

    const size_t length   = b->size();
    for ( size_t k = 0 ; k < length ; ++k ) {

      // get face
        stk::mesh::Entity face = (*b)[k];
      // face node relations
      stk::mesh::Entity const * face_node_rels = bulk.begin_nodes(face);
      // pointer to face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec_, face);

      for ( int ip = 0; ip < numScsBip; ++ip ) {
        // offsets
        const int offSetAveraVec = ip*nDim;
        // nearest node mapping to this ip
        const int localFaceNode = faceIpNodeMap[ip];
        // extract nodal fields
        stk::mesh::Entity node = face_node_rels[localFaceNode];
        double *assembledArea = stk::mesh::field_data(*assembledArea_, node );
        // aMag
        double aMag = 0.0;
        for ( int j = 0; j < nDim; ++j)
          aMag += areaVec[offSetAveraVec+j]*areaVec[offSetAveraVec+j];
        aMag = std::sqrt(aMag);

        // assemble nodal quantities
        *assembledArea += aMag;
      }
    }
  }

  // parallel assemble
  std::vector<const stk::mesh::FieldBase*> fields;
  fields.push_back(assembledArea_);
  const std::vector<const stk::mesh::FieldBase*>& const_fields = fields;
  stk::mesh::parallel_sum(bulk, const_fields);

  // periodic assemble
  if ( realm_.hasPeriodic_) {
      // fields are not defined at all slave/master node pairs
      const bool bypassFieldCheck = false; 
      realm_.periodic_field_update(assembledArea_, 1, bypassFieldCheck);
  }
    
}

void SurfaceFMPostProcessing::calc_surface_force(
    SurfaceFMData & sfm_data)
{

  // check to see if this is a valid step to process output file
  const int timeStepCount = realm_.get_time_step_count();
  const bool processMe =
      (timeStepCount % sfm_data.frequency_) == 0 ? true : false;

  // do not waste time here
  if ( !processMe )
    return;

  // common
  stk::mesh::BulkData & bulk = realm_.bulk_data();
  stk::mesh::MetaData & meta = realm_.meta_data();
  const int nDim = meta.spatial_dimension();

  // set min and max values
  double yplusMin = 1.0e8;
  double yplusMax = -1.0e8;

  // nodal fields to gather
  std::vector<double> ws_pressure;
  std::vector<double> ws_density;
  std::vector<double> ws_viscosity;

  // master element
  std::vector<double> ws_face_shape_function;

  // deal with state
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

  // define vector of parent topos; should always be UNITY in size
  std::vector<stk::topology> parentTopo;

  const double currentTime = realm_.get_current_time();

  // local force and moment; i.e., to be assembled
  std::array<double,9> l_force_moment;

  // work force, moment and radius; i.e., to be pushed to cross_product()
  std::array<double,3> ws_p_force;
  std::array<double,3> ws_v_force;
  std::array<double,3> ws_t_force;
  std::array<double,3> ws_tau;
  std::array<double,3> ws_moment;
  std::array<double,3> ws_radius;

  // will need surface normal
  std::array<double,3> ws_normal;

  // define some common selectors
  stk::mesh::Selector sel = meta.locally_owned_part()
    &stk::mesh::selectUnion(sfm_data.partVector_);
  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta.side_rank(), sel );
  for ( auto b: face_buckets ) {
    // face master element
    MasterElement *meFC =
        MasterElementRepo::get_surface_master_element(b->topology());
    const int nodesPerFace = meFC->nodesPerElement_;
    const int numScsBip = meFC->numIntPoints_;

    // mapping from ip to nodes for this ordinal;
    // face perspective (use with face_node_relations)
    const int *faceIpNodeMap = meFC->ipNodeMap();

    // extract connected element topology
    b->parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
    ThrowAssert ( parentTopo.size() == 1 );
    stk::topology theElemTopo = parentTopo[0];

    // extract master element for this element topo
    MasterElement *meSCS =
        MasterElementRepo::get_surface_master_element(theElemTopo);

    // algorithm related; element
    ws_pressure.resize(nodesPerFace);
    ws_density.resize(nodesPerFace);
    ws_viscosity.resize(nodesPerFace);
    ws_face_shape_function.resize(numScsBip*nodesPerFace);
    
    // shape functions
    if ( realm_.realmUsesEdges_ )
        meFC->shifted_shape_fcn(ws_face_shape_function.data());
    else
        meFC->shape_fcn(ws_face_shape_function.data());

    const size_t length = b->size();
    for ( size_t k = 0 ; k < length ; ++k ) {

      // get face
      stk::mesh::Entity face = (*b)[k];
      // face node relations
      stk::mesh::Entity const * face_node_rels = bulk.begin_nodes(face);
      //======================================
      // gather nodal data off of face
      //======================================
      for ( int ni = 0; ni < nodesPerFace; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];
        // gather scalars
        ws_pressure[ni] = *stk::mesh::field_data(*pressure_, node);
        ws_density[ni] = *stk::mesh::field_data(densityNp1, node);
        ws_viscosity[ni] = *stk::mesh::field_data(*viscosity_, node);
      }

      // pointer to face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec_, face);
      // extract the connected element to this exposed face;
      // should be single in size!
      const stk::mesh::Entity* face_elem_rels = bulk.begin_elements(face);
      ThrowAssert( bulk.num_elements(face) == 1 );
      // get element; its face ordinal number
      stk::mesh::Entity element = face_elem_rels[0];
      const int face_ordinal = bulk.begin_element_ordinals(face)[0];
      // get the relations off of element
      stk::mesh::Entity const * elem_node_rels = bulk.begin_nodes(element);

      for ( int ip = 0; ip < numScsBip; ++ip ) {
        // offsets
        const int offSetAveraVec = ip*nDim;
        const int offSetSF_face = ip*nodesPerFace;
        const int localFaceNode = faceIpNodeMap[ip];
        const int opposingNode = meSCS->opposingNodes(face_ordinal,ip);

        // interpolate to bip
        double pBip = 0.0;
        double rhoBip = 0.0;
        double muBip = 0.0;
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = ws_face_shape_function[offSetSF_face+ic];
          pBip += r*ws_pressure[ic];
          rhoBip += r*ws_density[ic];
          muBip += r*ws_viscosity[ic];
        }

        // extract nodal fields
        stk::mesh::Entity node = face_node_rels[localFaceNode];
        const double * coord = stk::mesh::field_data(*coordinates_, node );
        const double *duidxj = stk::mesh::field_data(*dudx_, node );
        double *pressureForce = stk::mesh::field_data(*pressureForce_, node );
        double *tauWall = stk::mesh::field_data(*tauWall_, node );
        double *yplus = stk::mesh::field_data(*yplus_, node );
        const double assembledArea =
            *stk::mesh::field_data(*assembledArea_, node );

        // divU and aMag
        double divU = 0.0;
        double aMag = 0.0;
        for ( int j = 0; j < nDim; ++j) {
          divU += duidxj[j*nDim+j];
          aMag += areaVec[offSetAveraVec+j]*areaVec[offSetAveraVec+j];
        }
        aMag = std::sqrt(aMag);

        // normal
        for ( int i = 0; i < nDim; ++i ) {
          ws_normal[i] = areaVec[offSetAveraVec+i]/aMag;
        }

        // load radius; assemble force -sigma_ij*njdS and compute tau_ij njDs
        for ( int i = 0; i < nDim; ++i ) {
          const double ai = areaVec[offSetAveraVec+i];
          ws_radius[i] = coord[i] - sfm_data.centroidCoords_[i];
          // set forces
          ws_v_force[i] = 2.0/3.0*muBip*divU*ai;
          tauWall[i] += 2.0/3.0*muBip*divU*ai;
          ws_p_force[i] = pBip*ai;
          pressureForce[i] += pBip*ai;
          double dflux = 0.0;
          double tauijNj = 0.0;
          const int offSetI = nDim*i;
          for ( int j = 0; j < nDim; ++j ) {
            const int offSetTrans = nDim*j+i;
            dflux += -muBip*(duidxj[offSetI+j] +
                             duidxj[offSetTrans])*areaVec[offSetAveraVec+j];
            tauijNj += -muBip*(duidxj[offSetI+j] +
                               duidxj[offSetTrans])*ws_normal[j];
          }
          // accumulate viscous force and set tau for component i
          ws_v_force[i] += dflux;
          tauWall[i] += dflux;
          ws_tau[i] = tauijNj;
        }
        
        // compute total force and tangential tau
        double tauTangential = 0.0;
        for ( int i = 0; i < nDim; ++i ) {
          ws_t_force[i] = ws_p_force[i] + ws_v_force[i];
          double tauiTangential = (1.0-ws_normal[i]*ws_normal[i])*ws_tau[i];
          for ( int j = 0; j < nDim; ++j ) {
            if ( i != j )
              tauiTangential -= ws_normal[i]*ws_normal[j]*ws_tau[j];
          }
          tauTangential += tauiTangential*tauiTangential;
        }

        // assemble nodal quantities;
        // scaled by area for L2 lumped nodal projection
        const double areaFac = aMag/assembledArea;

        cross_product(&ws_t_force[0], &ws_moment[0], &ws_radius[0]);

        // assemble force and moment
        for ( int j = 0; j < 3; ++j ) {
          l_force_moment[j] += ws_p_force[j];
          l_force_moment[j+3] += ws_v_force[j];
          l_force_moment[j+6] += ws_moment[j];
        }

        //==================
        // deal with yplus
        //==================

        // left and right nodes; right is on the face; left is the opposing node
        stk::mesh::Entity nodeL = elem_node_rels[opposingNode];
        stk::mesh::Entity nodeR = face_node_rels[localFaceNode];

        // extract nodal fields
        const double * coordL = stk::mesh::field_data(*coordinates_, nodeL );
        const double * coordR = stk::mesh::field_data(*coordinates_, nodeR );

        // determine yp; ~nearest opposing edge normal distance to wall
        double ypBip = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double nj = ws_normal[j];
          const double ej = coordR[j] - coordL[j];
          ypBip += nj*ej*nj*ej;
        }
        ypBip = std::sqrt(ypBip);

        const double tauW = std::sqrt(tauTangential);
        const double uTau = std::sqrt(tauW/rhoBip);
        const double yplusBip = rhoBip*ypBip/muBip*uTau;

        // nodal field
        *yplus += yplusBip*areaFac;

        // min and max
        yplusMin = std::min(yplusMin, yplusBip);
        yplusMax = std::max(yplusMax, yplusBip);

      }
    }
  }

  // parallel assemble and output
  double g_force_moment[9] = {};
  stk::ParallelMachine comm = NaluEnv::self().parallel_comm();
  
  // Parallel assembly of L2
  stk::all_reduce_sum(comm, &l_force_moment[0], &g_force_moment[0], 9);
  
  // min/max
  double g_yplusMin = 0.0, g_yplusMax = 0.0;
  stk::all_reduce_min(comm, &yplusMin, &g_yplusMin, 1);
  stk::all_reduce_max(comm, &yplusMax, &g_yplusMax, 1);
  
  // deal with file name and banner
  if ( NaluEnv::self().parallel_rank() == 0 ) {
      std::ofstream myfile;
      myfile.open(sfm_data.outputFileName_.c_str(), std::ios_base::app);
      myfile << std::setprecision(6) 
             << std::setw(16) 
             << currentTime << std::setw(16) 
             << g_force_moment[0] << std::setw(16)
             << g_force_moment[1] << std::setw(16)
             << g_force_moment[2] << std::setw(16)
             << g_force_moment[3] << std::setw(16)
             << g_force_moment[4] << std::setw(16)
             << g_force_moment[5] <<  std::setw(16)
             << g_force_moment[6] << std::setw(16)
             << g_force_moment[7] << std::setw(16)
             << g_force_moment[8] <<  std::setw(16)
             << g_yplusMin << std::setw(16) << g_yplusMax << std::endl;
      myfile.close();
  }
  
}


void SurfaceFMPostProcessing::calc_surface_force_wallfn(
    SurfaceFMData & sfm_data)
{

  // check to see if this is a valid step to process output file
  const int timeStepCount = realm_.get_time_step_count();
  const bool processMe =
      (timeStepCount % sfm_data.frequency_) == 0 ? true : false;

  // do not waste time here
  if ( !processMe )
    return;

  stk::mesh::BulkData & bulk = realm_.bulk_data();
  stk::mesh::MetaData & meta = realm_.meta_data();

  const int nDim = meta.spatial_dimension();

  // set min and max values
  double yplusMin = 1.0e8;
  double yplusMax = -1.0e8;

  // bip values
  std::array<double,3> uBip;
  std::array<double,3> uBcBip;
  std::array<double,3> unitNormal;

  // tangential work array
  std::array<double,3> uiTangential;
  std::array<double,3> uiBcTangential;

  // nodal fields to gather
  std::vector<double> ws_velocityNp1;
  std::vector<double> ws_bcVelocity;
  std::vector<double> ws_pressure;
  std::vector<double> ws_density;
  std::vector<double> ws_viscosity;

  // master element
  std::vector<double> ws_face_shape_function;

  // deal with state
  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);
  const double currentTime = realm_.get_current_time();

  // local force and MomentWallFunction; i.e., to be assembled
  std::array<double,9> l_force_moment;

  // work force, MomentWallFunction
  // and radius; i.e., to be pused to cross_product()
  std::array<double,3> ws_p_force;
  std::array<double,3> ws_v_force;
  std::array<double,3> ws_t_force;
  std::array<double,3> ws_moment;
  std::array<double,3> ws_radius;

  // define some common selectors
  stk::mesh::Selector sel = meta.locally_owned_part()
    &stk::mesh::selectUnion(sfm_data.partVector_);
  auto& face_buckets =
    realm_.get_buckets( meta.side_rank(), sel );
  for ( auto b: face_buckets) {
    // face master element
    MasterElement *meFC =
        MasterElementRepo::get_surface_master_element(b->topology());
    const int nodesPerFace = meFC->nodesPerElement_;
    const int numScsBip = meFC->numIntPoints_;

    // mapping from ip to nodes for this ordinal
    const int *faceIpNodeMap = meFC->ipNodeMap();

    // algorithm related; element
    ws_velocityNp1.resize(nodesPerFace*nDim);
    ws_bcVelocity.resize(nodesPerFace*nDim);
    ws_pressure.resize(nodesPerFace);
    ws_density.resize(nodesPerFace);
    ws_viscosity.resize(nodesPerFace);
    ws_face_shape_function.resize(numScsBip*nodesPerFace);

    // shape functions
    if ( realm_.realmUsesEdges_ )
        meFC->shifted_shape_fcn(ws_face_shape_function.data());
    else
        meFC->shape_fcn(ws_face_shape_function.data());

    const size_t length   = b->size();
    for ( size_t k = 0 ; k < length ; ++k ) {

      // get face
      stk::mesh::Entity face = (*b)[k];
      // face node relations
      stk::mesh::Entity const * face_node_rels = bulk.begin_nodes(face);
      //======================================
      // gather nodal data off of face
      //======================================
      for ( int ni = 0; ni < nodesPerFace; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];
        // gather scalars
        ws_pressure[ni]    = *stk::mesh::field_data(*pressure_, node);
        ws_density[ni]    = *stk::mesh::field_data(densityNp1, node);
        ws_viscosity[ni] = *stk::mesh::field_data(*viscosity_, node);
        // gather vectors
        double * uNp1 = stk::mesh::field_data(velocityNp1, node);
        double * uBc = stk::mesh::field_data(*bcVelocity_, node);
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
          ws_velocityNp1[offSet+j] = uNp1[j];
          ws_bcVelocity[offSet+j] = uBc[j];
        }
      }

      // pointer to face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec_, face);
      const double *wallNormalDistanceBip =
          stk::mesh::field_data(*wallNormalDistanceBip_, face);
      const double *wallFrictionVelocityBip =
          stk::mesh::field_data(*wallFrictionVelocityBip_, face);

      for ( int ip = 0; ip < numScsBip; ++ip ) {

        // offsets
        const int offSetAveraVec = ip*nDim;
        const int offSetSF_face = ip*nodesPerFace;
        const int localFaceNode = faceIpNodeMap[ip];
        
        // zero out vector quantities; squeeze in aMag
        double aMag = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          uBip[j] = 0.0;
          uBcBip[j] = 0.0;
          const double axj = areaVec[offSetAveraVec+j];
          aMag += axj*axj;
        }
        aMag = std::sqrt(aMag);

        // interpolate to bip
        double pBip = 0.0;
        double rhoBip = 0.0;
        double muBip = 0.0;
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = ws_face_shape_function[offSetSF_face+ic];
          pBip += r*ws_pressure[ic];
          rhoBip += r*ws_density[ic];
          muBip += r*ws_viscosity[ic];
          const int offSetFN = ic*nDim;
          for ( int j = 0; j < nDim; ++j ) {
            uBip[j] += r*ws_velocityNp1[offSetFN+j];
            uBcBip[j] += r*ws_bcVelocity[offSetFN+j];
          }
        }

        // form unit normal
        for ( int j = 0; j < nDim; ++j ) {
          unitNormal[j] = areaVec[offSetAveraVec+j]/aMag;
        }

        // determine tangential velocity
        double uTangential = 0.0;
        double uknk = 0.0;
        double ukBcnk = 0.0;
        for (int i=0; i < nDim; ++i) {
            uknk += uBip[i] * unitNormal[i];
            ukBcnk += uBcBip[i] * unitNormal[i];
        }
        for ( int i = 0; i < nDim; ++i ) {
            // save off tangential components and augment magnitude
            uiTangential[i] = uBip[i] - uknk * unitNormal[i];
            uiBcTangential[i] = uBcBip[i] - ukBcnk * unitNormal[i];
            uTangential += (uiTangential[i] - uiBcTangential[i]) *
                (uiTangential[i] - uiBcTangential[i]);
        }
        uTangential = std::sqrt(uTangential);

        // extract bip data
        const double yp = wallNormalDistanceBip[ip];
        const double utau= wallFrictionVelocityBip[ip];

        // determine yplus
        const double yplusBip = rhoBip*yp*utau/muBip;

        // min and max
        yplusMin = std::min(yplusMin, yplusBip);
        yplusMax = std::max(yplusMax, yplusBip);

        double lambda = muBip/yp*aMag;
        if ( yplusBip > yplusCrit_)
          lambda = rhoBip*kappa_*utau/std::log(elog_*yplusBip)*aMag;

        // extract nodal fields
        stk::mesh::Entity node = face_node_rels[localFaceNode];
        const double * coord = stk::mesh::field_data(*coordinates_, node );
        double *pressureForce = stk::mesh::field_data(*pressureForce_, node );
        double *tauWall = stk::mesh::field_data(*tauWall_, node );
        double *yplus = stk::mesh::field_data(*yplus_, node );
        const double assembledArea =
            *stk::mesh::field_data(*assembledArea_, node );

        // load radius; assemble force -sigma_ij*njdS
        for ( int i = 0; i < nDim; ++i ) {
          const double ai = areaVec[offSetAveraVec+i];
          ws_radius[i] = coord[i] - sfm_data.centroidCoords_[i];
          const double uDiff = uiTangential[i] - uiBcTangential[i];
          ws_p_force[i] = pBip*ai;
          ws_v_force[i] = lambda*uDiff;
          ws_t_force[i] = ws_p_force[i] + ws_v_force[i];
          pressureForce[i] += ws_p_force[i];
          tauWall[i] += ws_v_force[i];
        }

        cross_product(&ws_t_force[0], &ws_moment[0], &ws_radius[0]);

        // assemble for and moment
        for ( size_t j = 0; j < 3; ++j ) {
          l_force_moment[j] += ws_p_force[j];
          l_force_moment[j+3] += ws_v_force[j];
          l_force_moment[j+6] += ws_moment[j];
        }
        
        // deal with yplus
        *yplus += yplusBip*aMag/assembledArea;

      }
    }
  }

  // parallel assemble and output
  double g_force_moment[9] = {};
  stk::ParallelMachine comm = NaluEnv::self().parallel_comm();
  
  // Parallel assembly of L2
  stk::all_reduce_sum(comm, &l_force_moment[0], &g_force_moment[0], 9);
  
  // min/max
  double g_yplusMin = 0.0, g_yplusMax = 0.0;
  stk::all_reduce_min(comm, &yplusMin, &g_yplusMin, 1);
  stk::all_reduce_max(comm, &yplusMax, &g_yplusMax, 1);
  
  // deal with file name and banner
  if ( NaluEnv::self().parallel_rank() == 0 ) {
      std::ofstream myfile;
      myfile.open(sfm_data.outputFileName_.c_str(), std::ios_base::app);
      myfile << std::setprecision(6) 
             << std::setw(16) 
             << currentTime << std::setw(16) 
             << g_force_moment[0] << std::setw(16)
             << g_force_moment[1] << std::setw(16)
             << g_force_moment[2] << std::setw(16)
             << g_force_moment[3] << std::setw(16)
             << g_force_moment[4] << std::setw(16)
             << g_force_moment[5] << std::setw(16)
             << g_force_moment[6] << std::setw(16)
             << g_force_moment[7] << std::setw(16)
             << g_force_moment[8] << std::setw(16)
             << g_yplusMin << std::setw(16) << g_yplusMax << std::endl;
      myfile.close();
  }

}

//--------------------------------------------------------------------------
//-------- zero_fields -------------------------------------------------------
//--------------------------------------------------------------------------
void
SurfaceFMPostProcessing::zero_fields()
{

  // common
  stk::mesh::BulkData & bulk = realm_.bulk_data();
  stk::mesh::MetaData & meta = realm_.meta_data();

  // zero fields
  field_fill( meta, bulk, 0.0, *pressureForce_, realm_.get_activate_aura());
  field_fill( meta, bulk, 0.0, *tauWall_, realm_.get_activate_aura());
  field_fill( meta, bulk, 0.0, *yplus_, realm_.get_activate_aura());
  field_fill( meta, bulk, 0.0, *assembledArea_, realm_.get_activate_aura());

}

//--------------------------------------------------------------------------
//-------- parralel_assemble_fields ----------------------------------------
//--------------------------------------------------------------------------
void
SurfaceFMPostProcessing::parallel_assemble_fields()
{

  stk::mesh::BulkData & bulk = realm_.bulk_data();
  stk::mesh::MetaData & meta = realm_.meta_data();
  const size_t nDim = meta.spatial_dimension();

  stk::mesh::parallel_sum(bulk, {pressureForce_, tauWall_, yplus_});

  // periodic assemble
  if ( realm_.hasPeriodic_) {
    // fields are not defined at all slave/master node pairs
    const bool bypassFieldCheck = false; 
    realm_.periodic_field_update(pressureForce_, nDim, bypassFieldCheck);
    realm_.periodic_field_update(tauWall_, 1, bypassFieldCheck);
    realm_.periodic_field_update(yplus_, 1, bypassFieldCheck);
  }

}

//--------------------------------------------------------------------------
//-------- cross_product ----------------------------------------------------
//--------------------------------------------------------------------------
void
SurfaceFMPostProcessing::cross_product(
    double *force,
    double *cross,
    double *rad)
{
    cross[0] =   rad[1]*force[2] - rad[2]*force[1];
    cross[1] = -(rad[0]*force[2] - rad[2]*force[0]);
    cross[2] =   rad[0]*force[1] - rad[1]*force[0];
}


}
}
