#include "AirfoilSMD.h"

namespace sierra {
namespace nalu {

AirfoilSMD::AirfoilSMD()
{
  // Initialize all matrices here - We'll figure out how to get them from user
  // input later
  M.xx = 1.0;
  M.yy = 1.0;
  M.zz = 1.0;

}

void
AirfoilSMD::predict_states() {

  // Create a simple extrapolator from x_nm1 and x_n to x_np1 here

}


void
AirfoilSMD::predict_time_step(vs::Vector F_np1) {

  // Implement generalized alpha, or RK time integrator scheme here
  // Go from x_nm1, x_n to x_np1


}


}
}
