/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <SurfaceForceAndMomentAlgorithmDriver.h>
#include <Algorithm.h>
#include <AlgorithmDriver.h>
#include <FieldFunctions.h>
#include <FieldTypeDef.h>
#include <Realm.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

//==========================================================================
// Class Definition
//==========================================================================
// SurfaceForceAndMomentAlgorithmDriver - Drives fluids post processing algs
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SurfaceForceAndMomentAlgorithmDriver::SurfaceForceAndMomentAlgorithmDriver(
  Realm &realm)
  : AlgorithmDriver(realm)
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
SurfaceForceAndMomentAlgorithmDriver::~SurfaceForceAndMomentAlgorithmDriver()
{
  std::vector<Algorithm *>::iterator iter, iter_end;
  iter_end = algVec_.end();
  for(iter = algVec_.begin(); iter != iter_end; ++iter)
    delete *iter;
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
SurfaceForceAndMomentAlgorithmDriver::execute()
{

  // pre-work; basically, allow algs to load up nodal area(s)
  for ( size_t k = 0; k < algVec_.size(); ++k )
    algVec_[k]->pre_work();

  // assemble area
  parallel_assemble_area();

  // execute
  for ( size_t k = 0; k < algVec_.size(); ++k )
    algVec_[k]->execute();
  
  // parallel assembly
  parallel_assemble_fields();

}


} // namespace nalu
} // namespace Sierra
