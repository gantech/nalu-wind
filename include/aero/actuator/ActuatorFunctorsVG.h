// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#ifndef ACTUATORFUNCTORSVG_H_
#define ACTUATORFUNCTORSVG_H_

#include <aero/actuator/ActuatorTypes.h>
#include <aero/actuator/ActuatorBulkVG.h>
#include <aero/actuator/ActuatorFunctors.h>
#include <NaluParsedTypes.h>
#include <NaluEnv.h>

namespace sierra {
namespace nalu {

struct InterpActuatorDensity
{
  using execution_space = ActuatorFixedExecutionSpace;

  InterpActuatorDensity(
    ActuatorBulkVG& actBulk, stk::mesh::BulkData& stkBulk);

  void operator()(int index) const;

  ActuatorBulkVG& actBulk_;
  stk::mesh::BulkData& stkBulk_;
  VectorFieldType* coordinates_;
  ScalarFieldType* density_;
};

void ActVGWriteToFile(
  ActuatorBulkVG& actBulk, const ActuatorMetaVG& actMeta);

struct ActVGAssignVel
{
  using execution_space = ActuatorFixedExecutionSpace;

  ActVGAssignVel(ActuatorBulkVG& actBulk);
  void operator()(int index) const;

  ActDualViewHelper<ActuatorFixedMemSpace> helper_;
  ActFixVectorDbl velocity_;
  ActFixScalarDbl density_;
  ActFixVectorDbl points_;
  ActFixScalarInt offset_;
  const int debug_output_;
  const int turbId_;
};

void ActVGComputeForce(
  ActuatorBulkVG& actBulk, const ActuatorMetaVG& actMeta);

struct ActVGComputeThrustInnerLoop
{

  ActVGComputeThrustInnerLoop(ActuatorBulkVG& actBulk)
    : actBulk_(actBulk)
  {
  }
  void operator()(
    const uint64_t pointId,
    const double* nodeCoords,
    double* sourceTerm,
    const double dualNodalVolume,
    const double scvIp) const;
  void preloop() {}

  ActuatorBulkVG& actBulk_;
};

struct ActVGSpreadForceWhProjInnerLoop
{
  ActVGSpreadForceWhProjInnerLoop(ActuatorBulkVG& actBulk)
    : actBulk_(actBulk)
  {
  }

  void operator()(
    const uint64_t pointId,
    const double* nodeCoords,
    double* sourceTerm,
    const double dualNodalVolume,
    const double scvIp) const;
  void preloop();

  ActuatorBulkVG& actBulk_;
};

using ActVGComputeThrust = GenericLoopOverCoarseSearchResults<
  ActuatorBulkVG,
  ActVGComputeThrustInnerLoop>;

using ActVGSpreadForceWhProjection = GenericLoopOverCoarseSearchResults<
  ActuatorBulkVG,
  ActVGSpreadForceWhProjInnerLoop>;

} /* namespace nalu */
} /* namespace sierra */

#endif /* ACTUATORFUNCTORSVG_H_ */
