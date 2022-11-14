// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#ifndef ACTUATOREXECUTORSVG_H_
#define ACTUATOREXECUTORSVG_H_

#include <aero/actuator/ActuatorBulkVG.h>
#include <aero/actuator/ActuatorFunctorsVG.h>
#include <aero/actuator/UtilitiesActuator.h>
#include <aero/actuator/ActuatorExecutor.h>
#include <memory>

namespace sierra {
namespace nalu {

class ActuatorLineVG : public ActuatorExecutor
{
public:
  ActuatorLineVG(
    const ActuatorMetaVG& actMeta,
    ActuatorBulkVG& actBulk,
    stk::mesh::BulkData& stkBulk);

  virtual ~ActuatorLineVG(){};
  void operator()() final;

private:
  const ActuatorMetaVG& actMeta_;
  ActuatorBulkVG& actBulk_;
  stk::mesh::BulkData& stkBulk_;
  const int numActPoints_;
  const bool useSpreadActuatorForce_;
};

} /* namespace nalu */
} /* namespace sierra */

#endif /* ACTUATORLINEVG_H_ */
