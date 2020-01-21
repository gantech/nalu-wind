#include "mesh_motion/MotionWaves.h"

#include <NaluParsing.h>
#include "utils/ComputeVectorDivergence.h"

// stk_mesh/base/fem
#include <stk_mesh/base/FieldBLAS.hpp>

#include <cmath>

namespace sierra{
namespace nalu{

MotionWaves::MotionWaves(
	stk::mesh::MetaData& meta,
	const YAML::Node& node)
	: MotionBase()
{
  load(node);

  // declare divergence of mesh velocity for this motion
  ScalarFieldType *divV = &(meta.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "div_mesh_velocity"));
  stk::mesh::put_field_on_mesh(*divV, meta.universal_part(), nullptr);
  stk::mesh::field_fill(0.0, *divV);
}

void MotionWaves::load(const YAML::Node& node)
{
  // Get type of input prescribed wave
	
	get_if_present(node,"wave_motion_model", waveModel_, waveModel_);

  if (waveModel_=="Sinusoidal_full_domain"){
  get_if_present(node, "amplitude", amplitude_, amplitude_);
  get_if_present(node, "waveperiod", waveperiod_, waveperiod_);	
	get_if_present(node, "wavelength",wavelength_,wavelength_); 
	}
  else if (waveModel_=="Linear_prescribed"){
  get_if_present(node, "amplitude", amplitude_, amplitude_);
  get_if_present(node, "waveperiod", waveperiod_, waveperiod_);	
	get_if_present(node, "wavelength",wavelength_,wavelength_); 
	}
  else if (waveModel_=="StokesSecondOrder_prescribed"){
  get_if_present(node, "amplitude", amplitude_, amplitude_);
  get_if_present(node, "waveperiod", waveperiod_, waveperiod_);	
	get_if_present(node, "wavelength",wavelength_,wavelength_);  
  }
  else if (waveModel_=="StokesThirdOrder_prescribed"){
  
  }
	else if (waveModel_ == "HOS"){

	}
	else {
    throw std::runtime_error("invalid wave_motion model specified ");
	} 	     
  double eps = std::numeric_limits<double>::epsilon();
	
  // Time parameters
  get_if_present(node, "start_time", startTime_, startTime_);
  startTime_ = startTime_-eps;
  get_if_present(node, "end_time", endTime_, endTime_);
  endTime_ = endTime_+eps;
	get_if_present(node, "sea_level_z",sealevelz_,sealevelz_); 
  // get sea-level based on if it was defined

}
void MotionWaves::build_transformation(
  const double time,
  const double* xyz)
{
  if(time < (startTime_)) return;

  double motionTime = (time < endTime_)? time : endTime_;

	double k=2.*M_PI/wavelength_;
  double omega=2.*M_PI/waveperiod_;
  double phase=k*xyz[0]-omega*motionTime;

	ThreeDVecType curr_disp={};
	curr_disp[0]=0.;
	curr_disp[1]=0.;
	
	
	if(waveModel_== "Sinusoidal_full_domain"){
	curr_disp[2]=sealevelz_+amplitude_*std::cos(phase);
	}
	else if(waveModel_== "Linear_prescribed"){
	curr_disp[2]=sealevelz_+amplitude_*std::cos(phase)*std::exp(-0.1*xyz[2]/amplitude_);
	}
  else if (waveModel_ == "StokesSecondOrder_prescribed"){
  curr_disp[2]=sealevelz_+amplitude_*(std::cos(phase)+k*amplitude_*(3-omega*omega)/(4.*omega*omega*omega)*std::cos(2*phase));
  }
	else if (waveModel_ == "StokesThirdOrder_prescribed"){
	curr_disp[2]=sealevelz_*((1.0-1.0/16.0*(k*amplitude_)*(k*amplitude_))*std::cos(phase)+1.0/2.0*(k*amplitude_)*std::cos(2.*phase)+3./8.*(k*amplitude_)*(k*amplitude_)*std::cos(3*phase));
	}
	else if (waveModel_ =="HOS"){
	
	}
	else {
    throw std::runtime_error("invalid wave_motion model specified ");
	}
	translation_mat(curr_disp);
}

void MotionWaves::translation_mat(const ThreeDVecType& curr_disp)
{
  reset_mat(transMat_);

  // Build matrix for translating object
  transMat_[0][3] = curr_disp[0];
  transMat_[1][3] = curr_disp[1];
  transMat_[2][3] = curr_disp[2];
}

MotionBase::ThreeDVecType MotionWaves::compute_velocity(
  const double time,
  const TransMatType&  /* compTrans */,
  const double* mxyz,
  const double* /* cxyz */ )
{
  ThreeDVecType vel = {};

  if( (time < startTime_) || (time > endTime_) ) return vel;
  
  double motionTime = (time < endTime_)? time : endTime_;

	double k=2.*M_PI/wavelength_;
  double omega=2.*M_PI/waveperiod_;
  double phase=k*mxyz[0]-omega*motionTime;
  double VerticalWaveVelocity;
  double HorizontalWaveVelocity;
	
  if(waveModel_== "Sinusoidal_full_domain"){
  VerticalWaveVelocity = amplitude_*omega*std::sin(phase);  
  HorizontalWaveVelocity = 0.;
  }
  else if(waveModel_== "Linear_prescribed"){
  VerticalWaveVelocity = amplitude_*omega*std::sin(phase)*std::exp(-0.1*mxyz[2]/amplitude_);  
  HorizontalWaveVelocity = amplitude_*omega*std::cos(phase);
  }
  else if (waveModel_ == "StokesSecondOrder_prescribed"){
  VerticalWaveVelocity=0.;
  HorizontalWaveVelocity=0.;	
	}
	else if (waveModel_ == "StokesThirdOrder_prescribed"){
  VerticalWaveVelocity=0.;
  HorizontalWaveVelocity=0.;	
  }
	else if (waveModel_ =="HOS"){
  VerticalWaveVelocity=0.;
  HorizontalWaveVelocity=0.;		
	}
	else {
    throw std::runtime_error("invalid wave_motion model specified ");
	}

  double eps = std::numeric_limits<double>::epsilon();
  
 
  if(mxyz[2] < sealevelz_ + eps){
  vel[0] = HorizontalWaveVelocity;
  }else{
  vel[0] = 0.;
  }
	vel[1] = 0.;
	vel[2] = VerticalWaveVelocity ;
  
  return vel;
}

void MotionWaves::post_compute_geometry(
  stk::mesh::BulkData& bulk,
  stk::mesh::PartVector& partVec,
  stk::mesh::PartVector& partVecBc,
  bool& computedMeshVelDiv)
{
  if(computedMeshVelDiv) return;

  // compute divergence of mesh velocity
  VectorFieldType* meshVelocity = bulk.mesh_meta_data().get_field<VectorFieldType>(
    stk::topology::NODE_RANK, "mesh_velocity");

  ScalarFieldType* meshDivVelocity = bulk.mesh_meta_data().get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "div_mesh_velocity");

  compute_vector_divergence(bulk, partVec, partVecBc, meshVelocity, meshDivVelocity, true);
  computedMeshVelDiv = true;
}


} // nalu
} // sierra

