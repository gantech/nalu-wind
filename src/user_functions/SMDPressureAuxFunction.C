/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/SMDPressureAuxFunction.h>
#include <algorithm>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

SMDPressureAuxFunction::SMDPressureAuxFunction() :
  AuxFunction(0,1),
  u_infty_(2.0),
  sigma_(0.5),
  alpha_(1.0),
  omega_(1.0)
{
  // does nothing
}

void
SMDPressureAuxFunction::do_evaluate(
  const double *coords,
  const double time,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  const double ys = 0.0; //alpha_ * std::sin(omega_ * time);
  const double oneOverSigma2 =1.0/(sigma_*sigma_);

  for(unsigned p=0; p < numPoints; ++p) {

    const double x = coords[0];
    const double y = coords[1];

    fieldPtr[0] = u_infty_*u_infty_*(1.0 - ((y-ys)/sigma_) * std::exp( -(x*x + (y-ys)*(y-ys)) * oneOverSigma2) );

    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

SMDPressureGradientAuxFunction::SMDPressureGradientAuxFunction(
   const unsigned beginPos,
   const unsigned endPos) :
   AuxFunction(beginPos, endPos),
    u_infty_(2.0),
    sigma_(0.5),
    alpha_(1.0),
    omega_(1.0)
{
    // does nothing
}
    
void
SMDPressureGradientAuxFunction::do_evaluate(
  const double *coords,
  const double time,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
    const double ys = 0.0; //alpha_ * std::sin(omega_ * time);
    const double oneOverSigma = 1.0/sigma_;
    const double oneOverSigma2 =1.0/(sigma_*sigma_);
    
    for(unsigned p=0; p < numPoints; ++p) {
        
        const double x = coords[0];
        const double y = coords[1];

        const double expfn = std::exp( -(x*x + (y-ys)*(y-ys)) * oneOverSigma2);
        fieldPtr[0] = 2.0 * u_infty_*u_infty_ * x*(y-ys) * expfn * oneOverSigma2*oneOverSigma;
        fieldPtr[1] = u_infty_*u_infty_ * expfn * oneOverSigma * ( -1.0 + 2.0*(y-ys)*(y-ys)*oneOverSigma2 );
        
        fieldPtr += fieldSize;
        coords += spatialDimension;
    }
}
    
} // namespace nalu
} // namespace Sierra
