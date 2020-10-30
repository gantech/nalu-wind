// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//


#ifndef TurbViscKsgsAlg_h
#define TurbViscKsgsAlg_h

#include "Algorithm.h"
#include "FieldTypeDef.h"

#include "node_kernels/NodeKernel.h"
#include "stk_mesh/base/Ngp.hpp"
#include "stk_mesh/base/NgpField.hpp"
#include "stk_mesh/base/Types.hpp"

namespace sierra{
namespace nalu{

class Realm;

class TurbViscKsgsAlg : public Algorithm
{
public:
  using DblType = double;

  TurbViscKsgsAlg(
    Realm &realm,
    stk::mesh::Part* part,
    ScalarFieldType* tvisc);

  virtual ~TurbViscKsgsAlg() = default;

  virtual void execute() override;

private:
  ScalarFieldType* tviscField_ {nullptr};
  unsigned tke_  {stk::mesh::InvalidOrdinal};
  unsigned density_  {stk::mesh::InvalidOrdinal};
  unsigned tvisc_  {stk::mesh::InvalidOrdinal};
  unsigned dualNodalVolume_  {stk::mesh::InvalidOrdinal};
  unsigned stabLscale_ {stk::mesh::InvalidOrdinal};
  unsigned dhdx_ {stk::mesh::InvalidOrdinal};
  unsigned specificHeat_   {stk::mesh::InvalidOrdinal};

  NALU_ALIGNED NodeKernelTraits::DblType gravity_[NodeKernelTraits::NDimMax];
  NodeKernelTraits::DblType beta_{0.0};

  const DblType cmuEps_;
  const int    nDim_;

};

} // namespace nalu
} // namespace Sierra

#endif
