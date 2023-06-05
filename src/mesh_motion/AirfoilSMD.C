#include "mesh_motion/AirfoilSMD.h"
#include "NaluParsing.h"

#include "netcdf.h"

#include <sstream>

namespace sierra {
namespace nalu {

inline void
check_nc_error(int code, std::string msg)
{
  if (code != 0)
    throw std::runtime_error("AirfoilSMD:: NetCDF error: " + msg);
}
    
AirfoilSMD::AirfoilSMD(const YAML::Node& node)
  : SMD(node)
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
  C_.xx() = damp[0]; C_.xy() = damp[1]; C_.xz() = damp[2];
  C_.yx() = damp[3]; C_.yy() = damp[4]; C_.yz() = damp[5];
  C_.zx() = damp[6]; C_.zy() = damp[7]; C_.zz() = damp[8];

  std::vector<double> transform(9,0.0);
  transform[0] = 1.0; transform[4] = 1.0; transform[8] = 1.0;
  get_if_present(node, "force_transform_matrix", transform);
  T_.xx() = transform[0]; T_.xy() = transform[1]; T_.xz() = transform[2];
  T_.yx() = transform[3]; T_.yy() = transform[4]; T_.yz() = transform[5];
  T_.zx() = transform[6]; T_.zy() = transform[7]; T_.zz() = transform[8];

  std::vector<double> x_init;
  get_required(node, "x_init", x_init);
  x_n_.x() = x_init[0]; x_n_.y() = x_init[1]; x_n_.z() = x_init[2]; 

  std::vector<double> xdot_init(3,0.0);
  get_if_present(node, "xdot_init", xdot_init);
  xdot_n_.x() = xdot_init[0]; xdot_n_.y() = xdot_init[1]; xdot_n_.z() = xdot_init[2]; 

  std::vector<double> a_init(3,0.0); //acceleration
  get_if_present(node, "a_init", a_init);
  a_n_.x() = a_init[0]; a_n_.y() = a_init[1]; a_n_.z() = a_init[2]; 

  std::vector<double> x_nm1(3,0.0);
  get_if_present(node, "x_nm1", x_nm1);
  x_nm1_.x() = x_nm1[0]; x_nm1_.y() = x_nm1[1]; x_nm1_.z() = x_nm1[2]; 

  std::vector<double> xdot_nm1(3,0.0);
  get_if_present(node, "xdot_nm1", xdot_nm1);
  xdot_nm1_.x() = xdot_nm1[0]; xdot_nm1_.y() = xdot_nm1[1]; xdot_nm1_.z() = xdot_nm1[2];

  std::vector<double> a_nm1(3,0.0); //acceleration
  get_if_present(node, "a_nm1", a_nm1);
  a_nm1_.x() = a_nm1[0]; a_nm1_.y() = a_nm1[1]; a_nm1_.z() = a_nm1[2]; 

  std::vector<double> f_init(3,0.0);
  get_if_present(node, "f_init", f_init);
  f_n_.x() = f_init[0]; f_n_.y() = f_init[1]; f_n_.z() = f_init[2]; 

  std::vector<double> f_nm1(3,0.0);
  get_if_present(node, "f_nm1", f_nm1);
  f_nm1_.x() = f_nm1[0]; f_nm1_.y() = f_nm1[1]; f_nm1_.z() = f_nm1[2]; 

  std::vector<double> origin;
  get_required(node, "centroid", origin);
  origin_.x() = origin[0]; origin_.y() = origin[1]; origin_.z() = origin[2];

  get_if_present(node, "alpha", alpha_);
  if ( (alpha_ > 0.0) || (alpha_ < -0.33333) )
    throw std::runtime_error("AirfoilSMD:: alpha should be '-0.333 =< alpha =< 0.0'. Instead alpha is  " + std::to_string(alpha_) );
  
}

vs::Vector
AirfoilSMD::get_trans_disp()
{
  vs::Vector disp;
  disp[0] = x_np1_[0];
  disp[1] = x_np1_[1];
  disp[2] = 0.0;
  return disp;
}

double
AirfoilSMD::get_rot_disp()
{
  return x_np1_[2];
}

vs::Vector
AirfoilSMD::get_trans_vel()
{
  vs::Vector vel;
  vel[0] = xdot_np1_[0];
  vel[1] = xdot_np1_[1];
  vel[2] = 0.0;
  return vel;
}
    
vs::Vector
AirfoilSMD::get_rot_vel()
{
  vs::Vector vel;
  vel[0] = 0.0;
  vel[1] = 0.0;
  vel[2] = xdot_np1_[2];
  return vel;
}

vs::Vector
AirfoilSMD::get_rot_axis()
{
  vs::Vector axis;
  axis[0] = 0.0;
  axis[1] = 0.0;
  axis[2] = 1.0;
  return axis;
}
    
void
AirfoilSMD::setup(const double dt) {
  dt_ = dt;
}

void
AirfoilSMD::predict_states() {

  // Create a simple extrapolator from x_nm1 and x_n to x_np1 here

  double dt = 0.01;

  // Second order predictor
  x_np1_ = x_n_ + dt*(1.5*xdot_n_ - 0.5*xdot_nm1_);

  xdot_np1_ = xdot_n_ + dt*(1.5*a_n_ - 0.5*a_nm1_);


}

void
AirfoilSMD::advance_timestep() {

  x_nm1_ = x_n_;
  x_n_ = x_np1_;

  f_nm1_ = f_n_;
  f_n_ = f_np1_;

  xdot_nm1_ = xdot_n_;
  xdot_n_ = xdot_np1_;

  a_nm1_ = a_n_;
  a_n_ = a_np1_;

  tstep_ += 1;
}


void
AirfoilSMD::update_timestep(vs::Vector F_np1, vs::Vector M_np1) {

  // Implement generalized alpha, or RK time integrator scheme here
  // Go from x_nm1, x_n to x_np1
  double dt = 0.01;

  // Generalized Alpha Method parameters
  double beta = (1.0 - alpha_)*(1.0 - alpha_)/4.0;
  double gamma = (1.0 - 2.0*alpha_)/2.0;

  // temporary vector of forces before transform
  vs::Vector temp_fnp1;

  // Formulate update as left * a_np1_ = right
  vs::Tensor Left;
  vs::Vector right;

  Left = M_ + (((1 + alpha_)*dt*gamma)*C_) + (((1+alpha_)*dt*dt*beta)*K_);

  // Create force vector from appropriate force and moment entries
  temp_fnp1[0] = F_np1[0];
  temp_fnp1[1] = F_np1[1];
  temp_fnp1[2] = M_np1[2];

  f_np1_ = temp_fnp1;

  right = (C_ & (-1.0*(xdot_n_ + (1 + alpha_)*dt*(1-gamma)*a_n_)))
          + (K_ & (-1.0*(x_n_ + (1 + alpha_)*dt*xdot_n_ + (1 + alpha_)*0.5*dt*dt*(1 - 2*beta)*a_n_)))
          + (T_ & ( ((1 + alpha_) * f_np1_) -( alpha_ * f_n_) ));
   

  // Solve the matrix problem to get a_np1_
  a_np1_ = Left.inv() & right;

  x_np1_ = x_n_ + dt*xdot_n_ + 0.5*dt*dt*((1 - 2*beta)*a_n_ + 2*beta*a_np1_);
  xdot_np1_ = xdot_n_ + dt*((1 - gamma)*a_n_ + gamma*a_np1_);
}

void
AirfoilSMD::prepare_nc_file() {
  
  int ncid, n_dim, n_tsteps, varid, ierr;
  
  // Create the file
  std::string filename = "af_smd_deflloads.nc";
  ierr = nc_create(filename.c_str(), NC_CLOBBER, &ncid);
  check_nc_error(ierr, "nc_create");
  
  // Define dimensions
  ierr = nc_def_dim(ncid, "n_dim", 3, &n_dim);
  ierr = nc_def_dim(ncid, "n_tsteps", NC_UNLIMITED, &n_tsteps);

  const std::vector<int> matrix_dims{n_dim, n_dim};
  const std::vector<int> state_dims{n_tsteps, n_dim};
  const std::vector<int> vector_dims{n_dim};
    
  ierr = nc_def_var(ncid, "time", NC_DOUBLE, 1, &n_tsteps, &varid);
  nc_var_ids_["time"] = varid;

  ierr = nc_def_var(ncid, "mass_matrix", NC_DOUBLE, 2, matrix_dims.data(), &varid);
  nc_var_ids_["mass_matrix"] = varid;

  ierr = nc_def_var(ncid, "stiffness_matrix", NC_DOUBLE, 2, matrix_dims.data(), &varid);
  nc_var_ids_["stiffness_matrix"] = varid;

  ierr = nc_def_var(ncid, "damping_matrix", NC_DOUBLE, 2, matrix_dims.data(), &varid);
  nc_var_ids_["damping_matrix"] = varid;

  ierr = nc_def_var(ncid, "origin", NC_DOUBLE, 1, vector_dims.data(), &varid);
  nc_var_ids_["origin"] = varid;
  
  ierr = nc_def_var(ncid, "x", NC_DOUBLE, 2, state_dims.data(), &varid);
  nc_var_ids_["x"] = varid;

  ierr = nc_def_var(ncid, "xdot", NC_DOUBLE, 2, state_dims.data(), &varid);
  nc_var_ids_["xdot"] = varid;

  ierr = nc_def_var(ncid, "f", NC_DOUBLE, 2, state_dims.data(), &varid);
  nc_var_ids_["f"] = varid;

  //! Indicate that we are done defining variables, ready to write data
  ierr = nc_enddef(ncid);
  check_nc_error(ierr, "nc_enddef");

  
  {
    std::vector<size_t> count_dim{3};
    std::vector<size_t> start_dim{0};
    ierr = nc_put_vara_double(
        ncid, nc_var_ids_["origin"], start_dim.data(), count_dim.data(),
        origin_.data());
  }

  {
    std::vector<size_t> count_dim{1,3};
    {
    std::vector<size_t> start_dim{0,0};
    ierr = nc_put_vara_double(
      ncid, nc_var_ids_["mass_matrix"], start_dim.data(), count_dim.data(),
      M_.x().data());
    ierr = nc_put_vara_double(
      ncid, nc_var_ids_["stiffness_matrix"], start_dim.data(), count_dim.data(),
      K_.x().data());
    ierr = nc_put_vara_double(
      ncid, nc_var_ids_["damping_matrix"], start_dim.data(), count_dim.data(),
      C_.x().data());
    }
    {
    std::vector<size_t> start_dim{1,0};
    ierr = nc_put_vara_double(
      ncid, nc_var_ids_["mass_matrix"], start_dim.data(), count_dim.data(),
      M_.y().data());
    ierr = nc_put_vara_double(
      ncid, nc_var_ids_["stiffness_matrix"], start_dim.data(), count_dim.data(),
      K_.y().data());
    ierr = nc_put_vara_double(
      ncid, nc_var_ids_["damping_matrix"], start_dim.data(), count_dim.data(),
      C_.y().data());
    }
    {
    std::vector<size_t> start_dim{2,0};
    ierr = nc_put_vara_double(
      ncid, nc_var_ids_["mass_matrix"], start_dim.data(), count_dim.data(),
      M_.z().data());
    ierr = nc_put_vara_double(
      ncid, nc_var_ids_["stiffness_matrix"], start_dim.data(), count_dim.data(),
      K_.z().data());
    ierr = nc_put_vara_double(
      ncid, nc_var_ids_["damping_matrix"], start_dim.data(), count_dim.data(),
      C_.z().data());
    }
  }
  {
    size_t count0 = 1;
    double cur_time = 0.0;
    ierr = nc_put_vara_double(ncid, nc_var_ids_["time"], &tstep_, &count0, &cur_time);
    std::vector<size_t> count_dim{1,3};
    {
      std::vector<size_t> start_dim{0,0};
      ierr = nc_put_vara_double(
          ncid, nc_var_ids_["x"], start_dim.data(), count_dim.data(),
          x_n_.data());
      ierr = nc_put_vara_double(
          ncid, nc_var_ids_["xdot"], start_dim.data(), count_dim.data(),
          xdot_n_.data());
      ierr = nc_put_vara_double(
          ncid, nc_var_ids_["f"], start_dim.data(), count_dim.data(),
          f_n_.data());
    }
  }
  
  ierr = nc_close(ncid);
  check_nc_error(ierr, "nc_close");
}

void
AirfoilSMD::write_nc_def_loads(const double cur_time)
{

  int ncid, ierr;
  std::string filename = "af_smd_deflloads.nc";
  ierr = nc_open(filename.c_str(), NC_WRITE, &ncid);
  check_nc_error(ierr, "nc_open");
  ierr = nc_enddef(ncid);

  std::cerr << "tstep = " << tstep_ << std::endl;
  {
    size_t count0 = 1;
    ierr = nc_put_vara_double(ncid, nc_var_ids_["time"], &tstep_, &count0, &cur_time);
    std::vector<size_t> count_dim{1,3};
    {
      std::vector<size_t> start_dim{tstep_,0};
      ierr = nc_put_vara_double(
          ncid, nc_var_ids_["x"], start_dim.data(), count_dim.data(),
          x_n_.data());
      ierr = nc_put_vara_double(
          ncid, nc_var_ids_["xdot"], start_dim.data(), count_dim.data(),
          xdot_n_.data());
      ierr = nc_put_vara_double(
          ncid, nc_var_ids_["f"], start_dim.data(), count_dim.data(),
          f_n_.data());
    }
  }

  ierr = nc_close(ncid);
  check_nc_error(ierr, "nc_close");
    
}

}
}
