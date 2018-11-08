/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifdef NALU_USES_SMD

#include <ActuatorLineSMD.h>
#include <FieldTypeDef.h>
#include <NaluParsing.h>
#include <NaluEnv.h>
#include <Realm.h>
#include <Simulation.h>
#include <nalu_make_unique.h>

// master elements
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

// stk_search
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/IdentProc.hpp>

// basic c++
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <cmath>

namespace sierra{
namespace nalu{


// constructor
ActuatorLineSMDInfo::ActuatorLineSMDInfo():
    ActuatorInfo()
{
  // nothing to do
}


// destructor
ActuatorLineSMDInfo::~ActuatorLineSMDInfo()
{
  // nothing to do
}


// constructor
ActuatorLineSMDPointInfo::ActuatorLineSMDPointInfo(
  size_t globTurbId,
  Point centroidCoords,
  double searchRadius,
  Coordinates epsilon
  )
  : ActuatorPointInfo(
       centroidCoords,
       searchRadius,
       1.0e16,
       stk::mesh::Entity()
    ),
    globSMDId_(globTurbId),
    epsilon_(epsilon)
{
  // nothing to do
}


// destructor
ActuatorLineSMDPointInfo::~ActuatorLineSMDPointInfo()
{
  // nothing to do
}


// constructor
ActuatorLineSMD::ActuatorLineSMD(
  Realm &realm,
  const YAML::Node &node)
  : Actuator(realm, node)
{
  Actuator::load(node);
  // load the data
  load(node);

}


//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ActuatorLineSMD::~ActuatorLineSMD()
{
  p_smd.end(); // Call destructors in p_smd
}


// Multiply the point force by the weight at this element location.
void
ActuatorLineSMD::compute_node_force_given_weight(
  const int &nDim,
  const double &g,
  const double *pointForce,
  double *nodeForce)
{

  for ( int j = 0; j < nDim; ++j )
    nodeForce[j] = pointForce[j]*g;
}

/**
 * This method calculates the isotropic Gaussian projection of width epsilon of
 * a unit body force at the actuator point to another point at a distance *dis*
 * \f[
 * g(dis) = \frac{1}{\pi^{3/2}} \epsilon^3} e^{-\left( dis/ \epsilon \right)^2}
 * \f]
*/
double
ActuatorLineSMD::isotropic_Gaussian_projection(
  const int &nDim,
  const double &dis,
  const Coordinates &epsilon)
{
  // Compute the force projection weight at this location using an
  // isotropic Gaussian.
  double g;
  const double pi = acos(-1.0);
  if ( nDim == 2 )
    g = (1.0 / (pow(epsilon.x_,2.0) * pi)) * exp(-pow((dis/epsilon.x_),2.0));
  else
    g = (1.0 / (pow(epsilon.x_,3.0) * pow(pi,1.5))) * exp(-pow((dis/epsilon.x_),2.0));

  return g;
}


//--------------------------------------------------------------------------
//-------- load ------------------------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineSMD::load(
  const YAML::Node & y_node)
{
  // check for any data probes
  const YAML::Node y_actuatorLine = y_node["actuator"];
  if (y_actuatorLine) {
    NaluEnv::self().naluOutputP0() << "ActuatorLineSMD::load" << std::endl;

    // search specifications
    std::string searchMethodName = "na";
    get_if_present(y_actuatorLine, "search_method", searchMethodName, searchMethodName);

    // determine search method for this pair
    if ( searchMethodName == "boost_rtree" )
      searchMethod_ = stk::search::BOOST_RTREE;
    else if ( searchMethodName == "stk_kdtree" )
      searchMethod_ = stk::search::KDTREE;
    else
      NaluEnv::self().naluOutputP0() << "ActuatorLineSMD::search method not declared; will use stk_kdtree" << std::endl;

    // extract the set of from target names; each spec is homogeneous in this respect
    const YAML::Node searchTargets = y_actuatorLine["search_target_part"];
    if (searchTargets.Type() == YAML::NodeType::Scalar) {
      searchTargetNames_.resize(1);
      searchTargetNames_[0] = searchTargets.as<std::string>() ;
    }
    else {
      searchTargetNames_.resize(searchTargets.size());
      for (size_t i=0; i < searchTargets.size(); ++i) {
        searchTargetNames_[i] = searchTargets[i].as<std::string>() ;
      }
    }

    // Populate object of inputs class to SMD

    get_if_present(y_actuatorLine, "dry_run", i_smd.dry_run, false);
    get_if_present(y_actuatorLine, "debug", i_smd.debug, false);
    get_required(y_actuatorLine, "t_start", i_smd.t_start);
    get_required(y_actuatorLine, "n_every_checkpoint", i_smd.n_every_checkpoint);
    get_required(y_actuatorLine, "dt_smd", i_smd.dt);
    get_required(y_actuatorLine, "t_max", i_smd.t_max); // t_max is the total duration to which you want to run SMD.
    get_required(y_actuatorLine, "predictor", i_smd.predictor);

    // Only one SMD for now
    actuatorInfo_.emplace_back(new ActuatorLineSMDInfo());
    auto actuatorLineInfo = dynamic_cast<ActuatorLineSMDInfo*>(actuatorInfo_.back().get());

    const YAML::Node epsilon = y_actuatorLine["epsilon"];
    if ( epsilon )
      actuatorLineInfo->epsilon_ = epsilon.as<Coordinates>() ;
    else
      throw std::runtime_error("ActuatorLineSMD: lacking epsilon vector");

  }
}

/** Called after load, but before initialize. The mesh isn't loaded yet. For now, this function only
    checks that the Nalu time step is an integral multiple of the SMD time step
*/
void
ActuatorLineSMD::setup()
{
  // objective: declare the part, register coordinates; must be before populate_mesh()

  double dtNalu = realm_.get_time_step_from_file();
  tStepRatio_ = dtNalu/i_smd.dt ;
  if (std::abs(dtNalu - tStepRatio_ * i_smd.dt) < 0.001) {// TODO: Fix arbitrary number 0.001
    NaluEnv::self().naluOutputP0() << "Time step ratio  dtNalu/dtSMD: " << tStepRatio_ << std::endl ;
  } else {
    throw std::runtime_error("ActuatorLineSMD: Ratio of Nalu's time step is not an integral multiple of SMD time step");
  }
  i_smd.n_substeps = tStepRatio_;
  p_smd.setInputs(i_smd);

}

/** This function searches for the processor containing the base point of the SMD system and allocates the SMD
 *  to that processor. It does this through a stk::coarse_search of bounding boxes around the processor domains.
*/
void
ActuatorLineSMD::allocateSMDToProc()
{

    p_smd.setProcNo(0, 0);

}

/** This method allocates the turbines to processors, initializes the SMD, populates the  map of
 *  ActuatorLinePointInfo with the actuator point and determines the elements associated with each
 *  actuator point.
*/
void
ActuatorLineSMD::initialize()
{

  allocateSMDToProc();
  p_smd.init() ;
  
  if (NaluEnv::self().parallel_rank() == p_smd.get_procNo(0)) {
      std::vector<double> ws_pointGasVelocity(2);
      ws_pointGasVelocity[0] = 2.0; ws_pointGasVelocity[1] = 0.0; 
      p_smd.setVelocity_n(ws_pointGasVelocity, 0, 0);
      p_smd.setVelocity_nm1(ws_pointGasVelocity, 0, 0);
      p_smd.setVelocity_np1(ws_pointGasVelocity, 0, 0);
  }
  p_smd.solution0();
  
  update(); // Update location of actuator points, ghosting etc.
}


/**
 * This method should be called whenever the actuator points have moved and does the following:
 *
 * + creates a new map of actuator points in ActuatorLinePointInfoMap,
 * + searches the element bounding boxes for the elements within the search radius of each
 *   actuator point,
 * + identifies the elements to be ghosted to the processor controlling the smd,
 * + identifies the bestElem_ that contains each actuator point.
*/
void
ActuatorLineSMD::update()
{
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  // initialize need to ghost and elems to ghost
  needToGhostCount_ = 0;
  elemsToGhost_.clear();

  // clear actuatorLinePointInfoMap_
  actuatorPointInfoMap_.clear();

  bulkData.modification_begin();

  if ( actuatorGhosting_ == NULL) {
    // create new ghosting
    std::string theGhostName = "nalu_actuator_line_ghosting";
    actuatorGhosting_ = &bulkData.create_ghosting( theGhostName );
  }
  else {
    bulkData.destroy_ghosting(*actuatorGhosting_);
  }

  bulkData.modification_end();

  // clear some of the search info
  boundingSphereVec_.clear();
  boundingElementBoxVec_.clear();
  searchKeyPair_.clear();

  // set all of the candidate elements in the search target names
  populate_candidate_elements();

  // create the ActuatorLineSMDPointInfo
  create_actuator_line_point_info_map();

  // coarse search
  determine_elems_to_ghost();

  // manage ghosting
  manage_ghosting();

  // complete filling in the set of elements connected to the centroid
  complete_search();

}

void ActuatorLineSMD:: predict_struct_time_step() {

    if ( ! p_smd.isDryRun() )
        p_smd.update_states_driver_time_step();
        
}

// predict the state of OpenFAST at time zero
void ActuatorLineSMD::init_predict_struct_states() {

    if ( ! p_smd.isDryRun() ) {

        if ( p_smd.isTimeZero() ) {
            p_smd.solution0();
        }
    }
}

void ActuatorLineSMD:: advance_struct_time_step() {

    if ( ! p_smd.isDryRun() )
        p_smd.advance_to_next_driver_time_step();
}

void ActuatorLineSMD::sample_vel() {

  // meta/bulk data and nDim
  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  const int nDim = metaData.spatial_dimension();

  // extract fields
  VectorFieldType *coordinates
    = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  VectorFieldType *velocity = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  ScalarFieldType *density
    = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  ScalarFieldType *dualNodalVolume
    = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  
  // fixed size scratch
  std::vector<double> ws_pointGasVelocity(nDim);
  std::vector<double> ws_elemCentroid(nDim);
  std::vector<double> ws_pointForce(nDim);
  std::vector<double> ws_elemForce(nDim);
  double ws_pointGasDensity;

  // parallel communicate data to the ghosted elements; again can communicate points to element ranks
  if ( NULL != actuatorGhosting_ ) {
    std::vector< const stk::mesh::FieldBase *> ghostFieldVec;
    // fields that are needed
    ghostFieldVec.push_back(coordinates);
    ghostFieldVec.push_back(velocity);
    ghostFieldVec.push_back(dualNodalVolume);
    //    ghostFieldVec.push_back(viscosity);
    stk::mesh::communicate_field_data(*actuatorGhosting_, ghostFieldVec);
  }

  // loop over map and get velocity at points
  std::map<size_t, ActuatorLineSMDPointInfo *>::iterator iterPoint;
  int np=0;
  for (auto&& iterPoint : actuatorPointInfoMap_){

    // actuator line info object of interest
    auto infoObject = dynamic_cast<ActuatorLineSMDPointInfo*>(iterPoint.second.get());
    if( infoObject==NULL){
        throw std::runtime_error("Object in ActuatorPointInfo is not the correct type.  Should be ActuatorLineSMDPointInfo.");
    }

    //==========================================================================
    // extract the best element; compute drag given this velocity, property, etc
    // this point drag value will be used by all other elements below
    //==========================================================================
    stk::mesh::Entity bestElem = infoObject->bestElem_;
    int nodesPerElement = bulkData.num_nodes(bestElem);

    // resize some work vectors
    resize_std_vector(nDim, ws_coordinates_, bestElem, bulkData);
    resize_std_vector(nDim, ws_velocity_, bestElem, bulkData);
    //    resize_std_vector(1, ws_viscosity_, bestElem, bulkData);
    resize_std_vector(1, ws_density_, bestElem, bulkData);
    
    // gather nodal data to element nodes; both vector and scalar; coords are used in determinant calc
    gather_field(nDim, &ws_coordinates_[0], *coordinates, bulkData.begin_nodes(bestElem),
                 nodesPerElement);

    NaluEnv::self().naluOutput() << "Best element co-ordinates " << std::endl;
    NaluEnv::self().naluOutput() << ws_coordinates_[0] << " " << ws_coordinates_[1] << std::endl ;
    NaluEnv::self().naluOutput() << ws_coordinates_[2] << " " << ws_coordinates_[3] << std::endl ;
    NaluEnv::self().naluOutput() << ws_coordinates_[4] << " " << ws_coordinates_[5] << std::endl ;
    NaluEnv::self().naluOutput() << ws_coordinates_[6] << " " << ws_coordinates_[7] << std::endl ;

    NaluEnv::self().naluOutput() << "BestX" << infoObject->bestX_ << std::endl ;
    
    gather_field_for_interp(nDim, &ws_velocity_[0], *velocity, bulkData.begin_nodes(bestElem),
                            nodesPerElement);
    
    NaluEnv::self().naluOutput() << "Best element velocities " << std::endl;
    NaluEnv::self().naluOutput() << ws_velocity_[0] << " " << ws_velocity_[1] << std::endl ;
    NaluEnv::self().naluOutput() << ws_velocity_[2] << " " << ws_velocity_[3] << std::endl ;
    NaluEnv::self().naluOutput() << ws_velocity_[4] << " " << ws_velocity_[5] << std::endl ;
    NaluEnv::self().naluOutput() << ws_velocity_[6] << " " << ws_velocity_[7] << std::endl ;

    NaluEnv::self().naluOutput() << "IsoParCoords" << std::endl ;
    NaluEnv::self().naluOutput() << infoObject->isoParCoords_[0] << " " << infoObject->isoParCoords_[1] << " " << infoObject->isoParCoords_[2] << " " << infoObject->isoParCoords_[3] << std::endl ;
    
    gather_field_for_interp(1, &ws_density_[0], *density, bulkData.begin_nodes(bestElem),
                            nodesPerElement);

    // interpolate velocity
    interpolate_field(nDim, bestElem, bulkData, infoObject->isoParCoords_.data(),
                      &ws_velocity_[0], ws_pointGasVelocity.data());

    // interpolate density
    interpolate_field(1, bestElem, bulkData, infoObject->isoParCoords_.data(),
                      &ws_density_[0], &ws_pointGasDensity);

    NaluEnv::self().naluOutput() << "Velocity at point " << ws_pointGasVelocity[0] << ", " << ws_pointGasVelocity[1] << std::endl;
        
    p_smd.setVelocity_np1(ws_pointGasVelocity, np, infoObject->globSMDId_);
    np = np + 1;

  }
    
}


/** This function is called at each time step. This assembles the source terms
 *  in the momentum equation for Nalu.
*/
void
ActuatorLineSMD::execute()
{

  p_smd.extrapolateStatesVelDriver(); //Predict the velocity and states at time step 'n+1'
  
  update(); //Move the actuator point to the predicted location and do search, ghosting etc. 

  //Now get the force and spread to the nodes and assemble source term
  
  // meta/bulk data and nDim
  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  const int nDim = metaData.spatial_dimension();

  // extract fields
  VectorFieldType *coordinates
    = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  VectorFieldType *velocity = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  VectorFieldType *actuator_source
    = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "actuator_source");
  ScalarFieldType *actuator_source_lhs
    = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "actuator_source_lhs");
  ScalarFieldType *g
    = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "g");
  ScalarFieldType *dualNodalVolume
    = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

  // fixed size scratch
  std::vector<double> ws_pointGasVelocity(nDim);
  std::vector<double> ws_elemCentroid(nDim);
  std::vector<double> ws_pointForce(nDim);
  std::vector<double> ws_elemForce(nDim);

  // zero out source term; do this manually since there are custom ghosted entities
  stk::mesh::Selector s_nodes = stk::mesh::selectField(*actuator_source);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    double * actSrc = stk::mesh::field_data(*actuator_source, b);
    double * actSrcLhs = stk::mesh::field_data(*actuator_source_lhs, b);
    double * gF = stk::mesh::field_data(*g, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      actSrcLhs[k] = 0.0;
      gF[k] = 0.0;
      const int offSet = k*3;
      for ( int j = 0; j < nDim; ++j ) {
        actSrc[offSet+j] = 0.0;
      }
    }
  }

  // parallel communicate data to the ghosted elements; again can communicate points to element ranks
  if ( NULL != actuatorGhosting_ ) {
    std::vector< const stk::mesh::FieldBase *> ghostFieldVec;
    // fields that are needed
    ghostFieldVec.push_back(coordinates);
    ghostFieldVec.push_back(velocity);
    ghostFieldVec.push_back(dualNodalVolume);
    //    ghostFieldVec.push_back(viscosity);
    stk::mesh::communicate_field_data(*actuatorGhosting_, ghostFieldVec);
  }

  // loop over map and assemble source terms
  std::map<size_t, ActuatorLineSMDPointInfo *>::iterator iterPoint;
  int np=0;
  for (auto&& iterPoint : actuatorPointInfoMap_){
    
    // actuator line info object of interest
    auto infoObject = dynamic_cast<ActuatorLineSMDPointInfo*>(iterPoint.second.get());
    if( infoObject==NULL){
        throw std::runtime_error("Object in ActuatorPointInfo is not the correct type.  Should be ActuatorLineSMDPointInfo.");
    }

    //==========================================================================
    // extract the best element; compute drag given this velocity, property, etc
    // this point drag value will be used by all other elements below
    //==========================================================================
    p_smd.getForce_np1(ws_pointForce, np, infoObject->globSMDId_);
    NaluEnv::self().naluOutput() << "Predicte force at n+1 is  " << ws_pointForce[1] << std::endl ;
    
    // get the vector of elements
    std::set<stk::mesh::Entity> nodeVec = infoObject->nodeVec_;

    spread_actuator_force_to_node_vec(nDim, nodeVec, ws_pointForce, &(infoObject->centroidCoords_[0]), *coordinates, *actuator_source, *dualNodalVolume, infoObject->epsilon_);

    np=np+1;
  }

  // parallel assemble (contributions from ghosted and locally owned)
  const std::vector<const stk::mesh::FieldBase*> sumFieldVec(1, actuator_source);
  stk::mesh::parallel_sum_including_ghosts(bulkData, sumFieldVec);

}

// Creates bounding boxes around the subdomain of each processor
void
ActuatorLineSMD::populate_candidate_procs()
{
  stk::mesh::MetaData & metaData = realm_.meta_data();

  const int nDim = metaData.spatial_dimension();

  // fields
  VectorFieldType *coordinates = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // point data structures
  //  Point minCorner, maxCorner;
  std::vector<Point> minCorner(1), maxCorner(1);
  std::vector<Point> gMinCorner, gMaxCorner;

  // initialize max and min
  for (int j = 0; j < nDim; ++j ) {
    minCorner[0][j] = +1.0e16;
    maxCorner[0][j] = -1.0e16;
  }

  // extract part
  stk::mesh::PartVector searchParts;
  for ( size_t k = 0; k < searchTargetNames_.size(); ++k ) {
    stk::mesh::Part *thePart = metaData.get_part(searchTargetNames_[k]);
    if ( NULL != thePart )
      searchParts.push_back(thePart);
    else
      throw std::runtime_error("ActuatorLineSMD: Part is null" + searchTargetNames_[k]);
  }

  // selector and bucket loop
  stk::mesh::Selector s_locally_owned = metaData.locally_owned_part()
    &stk::mesh::selectUnion(searchParts);

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_locally_owned );

  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
        ib != elem_buckets.end() ; ++ib ) {

    stk::mesh::Bucket & b = **ib;

    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get element
      stk::mesh::Entity node = b[k];

      // pointers to real data
      const double * coords = stk::mesh::field_data(*coordinates, node );

      // check max/min
      for ( int j = 0; j < nDim; ++j ) {
	minCorner[0][j] = std::min(minCorner[0][j], coords[j]);
	maxCorner[0][j] = std::max(maxCorner[0][j], coords[j]);
      }

    }
  }

  stk::parallel_vector_concat(NaluEnv::self().parallel_comm(), minCorner, gMinCorner);
  stk::parallel_vector_concat(NaluEnv::self().parallel_comm(), maxCorner, gMaxCorner);

  for(int j = 0; j < NaluEnv::self().parallel_size(); j++) {
    // setup ident
    stk::search::IdentProc<uint64_t,int> theIdent(j, j);

    // create the bounding point box and push back
    boundingElementBox theBox(Box(gMinCorner[j],gMaxCorner[j]), theIdent);
    boundingProcBoxVec_.push_back(theBox);

  }

}

//--------------------------------------------------------------------------
//-------- create_actuator_line_point_info_map -----------------------------
//--------------------------------------------------------------------------
void
ActuatorLineSMD::create_actuator_line_point_info_map() {

  stk::mesh::MetaData & metaData = realm_.meta_data();
  const int nDim = metaData.spatial_dimension();

  size_t np = 0;

  for ( size_t iTurb = 0; iTurb < actuatorInfo_.size(); ++iTurb ) {

    const auto actuatorLineInfo = dynamic_cast<ActuatorLineSMDInfo*>(actuatorInfo_[iTurb].get());
    if(actuatorLineInfo==NULL){
        throw std::runtime_error("Object in ActuatorInfo is not correct type.  Should be ActuatorLineSMDInfo.");
    }

    int processorId = p_smd.get_procNo(iTurb);
    if ( processorId == NaluEnv::self().parallel_rank() ) {

      // define a point that will hold the centroid
      Point centroidCoords;

      // scratch array for coordinates and dummy array for velocity
      std::vector<double> currentCoords (3,0.0);

      // loop over all points for this turbine
      const int numForcePts = 1;

      if (! p_smd.isDryRun() ) {
	for(int iNode = 0; iNode < numForcePts; iNode++) {
	  stk::search::IdentProc<uint64_t,int> theIdent(np, NaluEnv::self().parallel_rank());

	  // set model coordinates from SMD
	  // move the coordinates
	  p_smd.getCoordinates_np1(currentCoords, 0, 0);
          NaluEnv::self().naluOutput() << "Predicted location at n+1 is  " << currentCoords[1] << std::endl ;
          double actCoord = std::sin(realm_.get_current_time());
          NaluEnv::self().naluOutput() << "Actual location at n+1 is supposed to be " << actCoord << std::endl ;
          
          double searchRadius = actuatorLineInfo->epsilon_.x_ * 2000.0;

	  for ( int j = 0; j < nDim; ++j )
	    centroidCoords[j] = currentCoords[j];

	  // create the bounding point sphere and push back
	  boundingSphere theSphere( Sphere(centroidCoords, searchRadius), theIdent);
	  boundingSphereVec_.push_back(theSphere);

	  // create the point info and push back to map
	  actuatorPointInfoMap_.insert(std::make_pair(np, make_unique<ActuatorLineSMDPointInfo>
                                                      (
                                                          iTurb, centroidCoords,
                                                          searchRadius, actuatorLineInfo->epsilon_
                                                      ))
              );

	  np=np+1;
	}

      }
      else {
	NaluEnv::self().naluOutput() << "Proc " << NaluEnv::self().parallel_rank() << " glob iTurb " << iTurb << std::endl ;
      }

    }

  }

}


//--------------------------------------------------------------------------
//-------- compute_elem_centroid -------------------------------------------
//--------------------------------------------------------------------------
void
ActuatorLineSMD::compute_elem_centroid(
  const int &nDim,
  double *elemCentroid,
  const int & nodesPerElement)
{
  // zero
  for ( int j = 0; j < nDim; ++j )
    elemCentroid[j] = 0.0;

  // assemble
  for ( int ni = 0; ni < nodesPerElement; ++ni ) {
    for ( int j=0; j < nDim; ++j ) {
      elemCentroid[j] += ws_coordinates_[ni*nDim+j]/nodesPerElement;
    }
  }
}


// Spread actuator force to nodes
void
ActuatorLineSMD::spread_actuator_force_to_node_vec(
  const int &nDim,
  std::set<stk::mesh::Entity>& nodeVec,
  const std::vector<double>& actuator_force,
  const double * actuator_node_coordinates,
  const stk::mesh::FieldBase & coordinates,
  stk::mesh::FieldBase & actuator_source,
  const stk::mesh::FieldBase & dual_nodal_volume,
  const Coordinates & epsilon
)
{
    std::vector<double> ws_nodeForce(nDim);

    std::set<stk::mesh::Entity>::iterator iNode;
    for (iNode = nodeVec.begin(); iNode != nodeVec.end(); ++iNode ) {

      stk::mesh::Entity node = *iNode;
      const double * node_coords = (double*)stk::mesh::field_data(coordinates, node );
      // compute distance
      const double distance = compute_distance(nDim, node_coords, actuator_node_coordinates);
      // project the force to this node with projection function
      double gA = isotropic_Gaussian_projection(nDim, distance, epsilon);
      compute_node_force_given_weight(nDim, gA, &actuator_force[0], &ws_nodeForce[0]);
      double * sourceTerm = (double*)stk::mesh::field_data(actuator_source, node );
      for ( int j=0; j < nDim; ++j ) sourceTerm[j] += ws_nodeForce[j];

    }

}

std::string ActuatorLineSMD::get_class_name(){
    return "ActuatorLineSMD";
}

} // namespace nalu
} // namespace Sierra

#endif
