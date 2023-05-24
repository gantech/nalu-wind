#include "mesh_motion/AirfoilSMD.h"

namespace sierra {
namespace nalu {

AirfoilSMD::AirfoilSMD()
{
  // Initialize all matrices here - We'll figure out how to get them from user
  // input later
  M.xx() = 1.0;
  M.yy() = 1.0;
  M.zz() = 1.0;

}

void
AirfoilSMD::predict_states() {

  // Create a simple extrapolator from x_nm1 and x_n to x_np1 here

  double dt = 0.01;

  // Prediction based solely on time n
  x_np1 = x_n + dt * v_n + (0.5*dt*dt)*a_n;

  // TODO: Consider second order schemes for prediction

}


void
AirfoilSMD::predict_time_step(vs::Vector F_np1) {

  //TODO: Verify that the matrix inverse is correct for a few cases.


  // Implement generalized alpha, or RK time integrator scheme here
  // Go from x_nm1, x_n to x_np1
  double dt = 0.01;

  // Generalized Alpha Method parameters
  double beta = (1.0 - alpha)*(1.0 - alpha)/4.0;
  double gamma = (1.0 - 2.0*alpha)/2.0;

  // Formulate update as left * a_np1 = right
  vs::Tensor Left;
  vs::Vector right;

  Left = M + ((1 + alpha)*dt*gamma)*C + ((1+alpha)*dt*dt*beta)*K;

  right = C & (-1.0*(v_n + (1 + alpha)*dt*(1-gamma)*a_n))
          + (K & (-1.0*(x_n + (1 + alpha)*dt*v_n + (1 + alpha)*0.5*dt*dt*(1 - 2*beta)*a_n)))
          + (1 + alpha) * F_np1 - alpha * f_n;
   
   f_np1 = F_np1;

  // Solve the matrix problem to get a_np1
  a_np1 = Left.inv() & right;

  x_np1 = x_n + dt*v_n + 0.5*dt*dt*((1 - 2*beta)*a_n + 2*beta*a_np1);
  v_np1 = v_n + dt*((1 - gamma)*a_n + gamma*a_np1);
}


}
}
