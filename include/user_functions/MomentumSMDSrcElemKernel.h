/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MOMENTUMSMDSRCELEMKERNEL_H
#define MOMENTUMSMDSRCELEMKERNEL_H

#include "kernel/Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class TimeIntegrator;
class SolutionOptions;
class MasterElement;
class ElemDataRequests;

/** CMM source term for MMS test for Spring-Mass-Damper system using actuator lines for momentum equation (velocity DOF)
 */
template<typename AlgTraits>
class MomentumSMDSrcElemKernel: public Kernel
{
public:
  MomentumSMDSrcElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    ElemDataRequests&,
    bool lumped);

  virtual ~MomentumSMDSrcElemKernel();

  /** Perform pre-timestep work for the computational kernel
   */
  virtual void setup(const TimeIntegrator&);

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  MomentumSMDSrcElemKernel() = delete;

  VectorFieldType *coordinates_{nullptr};

  DoubleType cur_time_;
  const DoubleType u_infty_;
  const DoubleType A_;
  const DoubleType sigma_;
  const DoubleType alpha_;
  const DoubleType omega_;
  const DoubleType M_;
  const DoubleType C_;
  const DoubleType K_;
  const DoubleType mu_;

  const int* ipNodeMap_;

  // scratch space
  Kokkos::View<DoubleType[AlgTraits::numScvIp_][AlgTraits::nodesPerElement_]> v_shape_function_ { "v_shape_func" };
};

}  // nalu
}  // sierra

#endif /* MOMENTUMSMDSRCELEMKERNEL_H */
