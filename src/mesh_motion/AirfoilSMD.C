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
  C_.xx() = stiff[0]; C_.xy() = stiff[1]; C_.xz() = stiff[2];
  C_.yx() = stiff[3]; C_.yy() = stiff[4]; C_.yz() = stiff[5];
  C_.zx() = stiff[6]; C_.zy() = stiff[7]; C_.zz() = stiff[8];

  std::vector<double> x_init;
  get_required(node, "x_init", x_init);
  x_n_.x() = x_init[0]; x_n_.y() = x_init[1]; x_n_.z() = x_init[2]; 

  std::vector<double> xdot_init(3,0.0);
  get_if_present(node, "xdot_init", xdot_init);
  xdot_n_.x() = xdot_init[0]; xdot_n_.y() = xdot_init[1]; xdot_n_.z() = xdot_init[2]; 

  std::vector<double> x_nm1(3,0.0);
  get_if_present(node, "x_nm1", x_nm1);
  x_nm1_.x() = x_nm1[0]; x_nm1_.y() = x_nm1[1]; x_nm1_.z() = x_nm1[2]; 

  std::vector<double> xdot_nm1(3,0.0);
  get_if_present(node, "xdot_nm1", xdot_nm1);
  xdot_nm1_.x() = xdot_nm1[0]; xdot_nm1_.y() = xdot_nm1[1]; xdot_nm1_.z() = xdot_nm1[2];

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

void
AirfoilSMD::setup(const double dt) {
  dt_ = dt;
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

  tstep_ += 1;
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
   
   f_np1_[0] = F_np1[0];
   f_np1_[1] = F_np1[1];
   f_np1_[2] = M_np1[2];

  // Solve the matrix problem to get a_np1
  a_np1_ = Left.inv() & right;

  x_np1_ = x_n_ + dt*v_n_ + 0.5*dt*dt*((1 - 2*beta)*a_n_ + 2*beta*a_np1_);
  v_np1_ = v_n_ + dt*((1 - gamma)*a_n_ + gamma*a_np1_);
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
