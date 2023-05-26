#include "mesh_motion/AirfoilSMD.h"
#include "NaluParsing.h"

namespace sierra {
namespace nalu {

AirfoilSMD::AirfoilSMD(const YAML::Node& node)
  : SMD()
{
  std::vector<double> mass;
  get_required(node, "mass_matrix", mass);
  M_.xx() = mass[0]; M_.xy() = mass[1]; M_.xz() = mass[2];
  M_.yx() = mass[3]; M_.yy() = mass[4]; M_.yz() = mass[5];
  M_.zx() = mass[6]; M_.zy() = mass[7]; M_.zz() = mass[8];

  std::vector<double> stiff;
  get_required(node, "stiffness_matrix", stiff);
  K_.xx() = stiff[0]; K_.xy() = stiff[1]; K_.xz() = stiff[2];
  K_.yx() = stiff[3]; K_.yy() = stiff[4]; K_.yz() = stiff[5];
  K_.zx() = stiff[6]; K_.zy() = stiff[7]; K_.zz() = stiff[8];

  std::vector<double> damp;
  get_required(node, "damping_matrix", damp);
  C_.xx() = stiff[0]; C_.xy() = stiff[1]; C_.xz() = stiff[2];
  C_.yx() = stiff[3]; C_.yy() = stiff[4]; C_.yz() = stiff[5];
  C_.zx() = stiff[6]; C_.zy() = stiff[7]; C_.zz() = stiff[8];

  std::vector<double> x_init;
  get_required(node, "x_init", x_init);
  x_n_.x() = x_init[0]; x_n_.y() = x_init[1]; x_n_.z() = x_init[2]; 

  std::vector<double> xdot_init;
  get_if_present(node, "xdot_init", xdot_init);
  xdot_n_.x() = xdot_init[0]; xdot_n_.y() = xdot_init[1]; xdot_n_.z() = xdot_init[2]; 

  std::vector<double> x_nm1;
  get_if_present(node, "x_nm1", x_nm1);
  x_nm1_.x() = x_nm1[0]; x_nm1_.y() = x_nm1[1]; x_nm1_.z() = x_nm1[2]; 

  std::vector<double> xdot_nm1;
  get_if_present(node, "xdot_nm1", xdot_nm1);
  xdot_nm1_.x() = xdot_nm1[0]; xdot_nm1_.y() = xdot_nm1[1]; xdot_nm1_.z() = xdot_nm1[2];

  std::vector<double> f_init;
  get_if_present(node, "f_init", f_init);
  f_n_.x() = f_init[0]; f_n_.y() = f_init[1]; f_n_.z() = f_init[2]; 

  std::vector<double> f_nm1;
  get_if_present(node, "f_nm1", f_nm1);
  f_nm1_.x() = f_nm1[0]; f_nm1_.y() = f_nm1[1]; f_nm1_.z() = f_nm1[2]; 
  
}

void
AirfoilSMD::predict_states() {

  // Create a simple extrapolator from x_nm1 and x_n to x_np1 here

  double dt = 0.01;

  // Prediction based solely on time n
  x_np1_ = x_n_ + dt * v_n_ + (0.5*dt*dt)*a_n_;

  // TODO: Consider second order schemes for prediction

}

void
AirfoilSMD::advance_timestep() {
  x_nm1_ = x_n_;
  x_n_ = x_np1_;
  f_nm1_ = f_n_;
  f_n_ = f_np1_;
  xdot_nm1_ = xdot_n_;
  xdot_n_ = xdot_np1_;
}


void
AirfoilSMD::update_timestep(vs::Vector F_np1, vs::Vector M_np1) {

  //TODO: Verify that the matrix inverse is correct for a few cases.

  // Implement generalized alpha, or RK time integrator scheme here
  // Go from x_nm1, x_n to x_np1
  double dt = 0.01;

  // Generalized Alpha Method parameters
  double beta = (1.0 - alpha_)*(1.0 - alpha_)/4.0;
  double gamma = (1.0 - 2.0*alpha_)/2.0;

  // Formulate update as left * a_np1 = right
  vs::Tensor Left;
  vs::Vector right;

  Left = M_ + ((1 + alpha_)*dt*gamma)*C_ + ((1+alpha_)*dt*dt*beta)*K_;

  right = C_ & (-1.0*(v_n_ + (1 + alpha_)*dt*(1-gamma)*a_n_))
          + (K_ & (-1.0*(x_n_ + (1 + alpha_)*dt*v_n_ + (1 + alpha_)*0.5*dt*dt*(1 - 2*beta)*a_n_)))
          + (1 + alpha_) * F_np1 - alpha_ * f_n_;
   
   f_np1_ = F_np1;

  // Solve the matrix problem to get a_np1
  a_np1_ = Left.inv() & right;

  x_np1_ = x_n_ + dt*v_n_ + 0.5*dt*dt*((1 - 2*beta)*a_n_ + 2*beta*a_np1_);
  v_np1_ = v_n_ + dt*((1 - gamma)*a_n_ + gamma*a_np1_);
}


}
}
