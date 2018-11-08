/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifdef NALU_USES_SMD

/** @file ActuatorLineSMD.h
 *  @brief A class to couple Nalu with SMD for actuator line simulation of a spring-mass-damper system
 *
 */

#ifndef ActuatorLineSMD_h
#define ActuatorLineSMD_h

#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include "Actuator.h"

// SMD C++ API
#include "smd.H"

namespace sierra{
namespace nalu{

class Realm;

/** Class that holds all of the information relevant to each spring-mass-damper system
 *
 *
 */
//
class ActuatorLineSMDInfo : public ActuatorInfo{
public:
  ActuatorLineSMDInfo();
  ~ActuatorLineSMDInfo();
  Coordinates epsilon_; ///< The Gaussian spreading width in (chordwise, spanwise, thickness) directions
};

/** Class that holds all of the search action for each actuator point
 *
 *
 */
//
class ActuatorLineSMDPointInfo : public ActuatorPointInfo{
public:
  ActuatorLineSMDPointInfo(
			    size_t globSMDId, Point centroidCoords, double searchRadius, Coordinates epsilon);
  ~ActuatorLineSMDPointInfo();
  size_t globSMDId_; ///< Global spring-mass-damper system number.
  Coordinates epsilon_; ///< The Gaussian spreading width in (chordwise, spanwise, thickness) directions for this actuator point.
};

/** The ActuatorLineSMD class couples Nalu with the third party library SMD for actuator line simulations of a spring-mass-damper system
 *
 * SMD is a library that models a point spring-mass-damper system. The effect of the turbine
 * on the flow field is modeled using the actuator line approach. The force exerted by the wind
 * turbine on the flow field is lumpled into a set of body forces at a discrete set of actuator
 * points. This class spreads the the body force at each actuator point using a Gaussian function.

 * 1) During the load phase - the turbine data from the yaml file is read and stored in an
 *    object of the smd::smdInputs class

 * 2) During the initialize phase - The processor containing the hub of each turbine is found
 *    through a search and assigned to be the one controlling SMD for that turbine. All
 *    processors controlling > 0 turbines initialize FAST, populate the map of ActuatorLinePointInfo
 *    and initialize element searches for all the actuator points associated with the turbines.
 *
 * 3) Elements are ghosted to the owning point rank. We tried the opposite approach of
 *    ghosting the actuator points to the processor owning the elements. The second approach
 *    was found to peform poorly compared to the first method.
 *
 * 4) A time lagged simple FSI model is used to interface Nalu with the turbine model:
 *    + The velocity at time step at time step 'n' is sampled at the actuator points and sent
 *       to SMD
 *    + SMD advances the spring-mass-damper system upto the next Nalu time step 'n+1'
 *    + The body forces at the actuator points are converted to the source terms of the momentum
 *      equation to advance Nalu to the next time step 'n+1'.
 *
 * 5) During the execute phase called every time step, we sample the velocity at the spring-mass-damper
 *    point and pass it to SMD. The spring-mass-damper system is advanced upto Nalu's
 *    next time step to get the body forces at the actuator point. We then iterate over the
 *    ActuatorLinePointInfoMap (now only single point) to assemble source terms.
 *
 *    actuator:
 *     type: ActLineSMD
 *     search_method: boost_rtree
 *     search_target_part: Unspecified-2-HEX
 *     dry_run:  False
 *     debug:    False
 *     t_max:    5.0
 *     n_every_checkpoint: 100
 *     epsilon: [ 5.0, 5.0, 5.0 ]
 */

class ActuatorLineSMD: public Actuator {
 public:

  ActuatorLineSMD(
    Realm &realm,
    const YAML::Node &node);
  ~ActuatorLineSMD();

  // load all of the options
  void load(
    const YAML::Node & node) override;

  // load the options for each turbine
  void readSMDData(int iTurb, smd::smdInputs & fi, YAML::Node turbNode);

  // setup part creation and nodal field registration (before populate_mesh())
  void setup() override;

  // allocate smd to processor containing base location
  void allocateSMDToProc() ;

  // Allocate SMD to to processor0, initialize SMD and get location of actuator points
  void initialize() override;

  // setup part creation and nodal field registration (after populate_mesh())
  void update();

  // determine processor bounding box in the mesh
  void populate_candidate_procs();

  // fill in the map that will hold point and ghosted elements
  void create_actuator_line_point_info_map();

  // predict initial structural model states
  void init_predict_struct_states();
  
  // predict the state of the structural model at the next time step
  void predict_struct_time_step();

  // firmly advance the state of the structural model to the next time step
  void advance_struct_time_step();
  
  // sample velocity at the actuator points and send to the structural model
  void sample_vel();
  
  // populate nodal field and output norms (if appropriate)
  void execute() override;

  // centroid of the element
  void compute_elem_centroid(
    const int &nDim,
    double *elemCentroid,
    const int &nodesPerElement);

  // compute the body force at an element given a
  // projection weighting.
  void compute_node_force_given_weight(
    const int &nDim,
    const double &g,
    const double *pointForce,
    double *nodeForce);

  // isotropic Gaussian projection function.
  double isotropic_Gaussian_projection(
    const int &nDim,
    const double &dis,
    const Coordinates &epsilon);

  // finally, perform the assembly

  void spread_actuator_force_to_node_vec(
    const int &nDim,
    std::set<stk::mesh::Entity>& nodeVec,
    const std::vector<double>& actuator_force,
    const double * actuator_node_coordinates,
    const stk::mesh::FieldBase & coordinates,
    stk::mesh::FieldBase & actuator_source,
    const stk::mesh::FieldBase & dual_nodal_volume,
    const Coordinates & epsilon);
      
  int tStepRatio_;  ///< Ratio of Nalu time step to SMD time step (dtNalu/dtSMD) - Should be an integral number

  // bounding box data types for stk_search
  std::vector<boundingSphere> boundingHubSphereVec_; ///< bounding box around the hub point of each turbine
  std::vector<boundingElementBox> boundingProcBoxVec_; ///< bounding box around all the nodes residing locally on each processor

  std::string get_class_name() override;

  smd::smdInputs i_smd; ///< Object to hold input information for SMD
  smd::smd p_smd; ///< SMD API handle

  std::vector<std::vector<double>> thrust;
  std::vector<std::vector<double>> torque;

};


} // namespace nalu
} // namespace Sierra

#endif

#endif
