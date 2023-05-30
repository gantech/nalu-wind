#ifndef SMD_H
#define SMD_H

#include "yaml-cpp/yaml.h"

#include "vs/vector.h"
#include "vs/tensor.h"

namespace sierra {
namespace nalu {


/** Spring-Mass-Damper system to predict structural response of an object

   Dummy class for a generic 6-DOF system
**/
class SMD
{
public:

    SMD(const YAML::Node& /* node */) {}

    virtual ~SMD() {};

    virtual void predict_states() = 0;

    virtual void update_timestep(vs::Vector F_np1, vs::Vector M_np1) = 0;

    virtual void advance_timestep() = 0;

    virtual void prepare_nc_file() = 0;

    virtual void write_nc_def_loads(const double cur_time) = 0;

    virtual vs::Vector get_origin() = 0;

private:

    SMD(const SMD&) = delete;
};

} // namespace nalu
} // namespace sierra

#endif /* SMD_H */
