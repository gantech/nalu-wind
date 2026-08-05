// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#ifndef KYNEMAFMBINTERFACE_H
#define KYNEMAFMBINTERFACE_H

#include <memory>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra {
namespace kynema_ugf {

class KynemaFMBInterface
{
public:
  virtual ~KynemaFMBInterface() = default;

  virtual void
  setup(double dtKynemaUGF, std::shared_ptr<stk::mesh::BulkData> bulk) = 0;

  virtual void initialize(int restartFreqKynemaUGF, double curTime) = 0;

  virtual void map_displacements(double currentTime, bool updateCurCoor) = 0;

  virtual void
  advance_struct_timestep(const double currentTime, const double dT) = 0;

  virtual void map_loads(const double currentTime) = 0;

  virtual stk::mesh::PartVector get_mesh_blocks() const = 0;
};

} // namespace kynema_ugf
} // namespace sierra

#endif