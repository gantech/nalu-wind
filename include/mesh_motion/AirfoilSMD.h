#ifndef AirfoilSMD_H
#define AirfoilSMD_H

#include "vs/vector.h"
#include "vs/tensor.h"

namespace sierra {
namespace nalu {


/** Spring-Mass-Damper system to predict structural response of an airfoil

   3-DOF system with flap, edge and twist


*/
class AirfoilSMD
{
public:

    AirfoilSMD();

    virtual ~AirfoilSMD();

    void predict_states();

    void predict_time_step(vs::Vector F_np1);

    namespace sierra {
        namespace nalu {
private:

    AirfoilSMD(const AirfoilSMD&) = delete;


    vs::Tensor M;
    vs::Tensor C;
    vs::Tensor K;
    vs::Vector x_np1;
    vs::Vector x_n;
    vs::Vector x_nm1;

    vs::Vector f_nm1;
    vs::Vector f_n;
    vs::Vector f_np1;

};

} // namespace nalu
} // namespace sierra

#endif /* AirfoilSMD_H */
