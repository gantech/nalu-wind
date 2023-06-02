#ifndef AirfoilSMD_H
#define AirfoilSMD_H

#include "yaml-cpp/yaml.h"

#include "mesh_motion/SMD.h"
#include "vs/vector.h"
#include "vs/tensor.h"

namespace sierra {
namespace nalu {


/** Spring-Mass-Damper system to predict structural response of an airfoil

   3-DOF system with flap, edge and twist


*/
class AirfoilSMD: public SMD
{
public:

    AirfoilSMD(const YAML::Node& node);

    virtual ~AirfoilSMD() {}

    void setup(double dt);

    void predict_states();

    void update_timestep(vs::Vector F_np1, vs::Vector M_np1);

    void advance_timestep();

    vs::Vector get_origin() { return origin_; }

    //! Prepare netcdf file to write deflections and loads
    void prepare_nc_file();
    
    //! Write deflections and loads to netcdf file
    void
    write_nc_def_loads(const double cur_time);

    vs::Vector get_trans_disp();

    double get_rot_disp();
    
    vs::Vector get_trans_vel();

    vs::Vector get_rot_vel();

    vs::Vector get_rot_axis();

private:

    AirfoilSMD(const AirfoilSMD&) = delete;

    double alpha_;

    vs::Tensor M_;
    vs::Tensor C_;
    vs::Tensor K_;
    
    vs::Vector x_np1_;
    vs::Vector x_n_;
    vs::Vector x_nm1_;

    vs::Vector xdot_np1_;    
    vs::Vector xdot_n_;    
    vs::Vector xdot_nm1_;
    
    vs::Vector v_np1_;
    vs::Vector v_n_;
    vs::Vector v_nm1_;

    vs::Vector a_np1_;
    vs::Vector a_n_;
    vs::Vector a_nm1_;

    vs::Vector f_np1_;
    vs::Vector f_n_;
    vs::Vector f_nm1_;

    vs::Vector origin_;

    double dt_{-1.0};
    
    size_t tstep_ {0};

    //! Map of `{variableName : netCDF_ID}` obtained from the NetCDF C interface
    std::unordered_map<std::string, int> nc_var_ids_;

};

} // namespace nalu
} // namespace sierra

#endif /* AirfoilSMD_H */
