#ifndef MODALSHAPEANALYSIS_H
#define MODALSHAPEANALYSIS_H

#include "aero/fsi/FSIturbine.h"
#include "yaml-cpp/yaml.h"

#include <array>
#include <memory>

namespace sierra {

namespace nalu {


class ModeShapeAnalysis
{
public:
  ModeShapeAnalysis(const YAML::Node&);
  virtual ~ModeShapeAnalysis() = default;

  void setup(double dtNalu, std::shared_ptr<stk::mesh::BulkData> bulk);

  void initialize(int restartFreqNalu, double curTime);

  void map_displacements(double, bool);

  void predict_struct_states();

  void predict_struct_timestep(const double curTime);

  void advance_struct_timestep(const double curTime);

  void compute_div_mesh_velocity();

  fsiTurbine* get_fsiTurbineData()
  {
    return fsiTurbineData_.get();
  }

  bool get_meshmotion() { return mesh_motion_; }

  void map_loads(const int tStep, const double curTime);

  void set_rotational_displacement(
    std::array<double, 3> axis, double omega, double curTime);
  void end_openfast();

  double total_openfastfsi_execution_time() { return 0.0; }
  double total_nalu_fsi_execution_time() { return naluTimer_.second; }

private:
  ModeShapeAnalysis() = delete;
  ModeShapeAnalysis(const ModeShapeAnalysis&) = delete;

  void load(const YAML::Node&);

  void get_displacements(double);

  void compute_mapping();

  void send_loads(const double curTime);
  void timer_start(std::pair<double, double>& timer);
  void timer_stop(std::pair<double, double>& timer);

  std::shared_ptr<stk::mesh::BulkData> bulk_;

  std::vector<std::string> partNames_;

  std::unique_ptr<fsiTurbine> fsiTurbineData_;

  std::string ncFileName_ ;
  // Mode shape in local (blade root) frame at finite element nodes
  std::vector<std::array<double, 6>> modeShape_;
  // Mode shape phase in local (blade root) frame at finite element nodes
  std::vector<std::array<double, 6>> modeShapePhase_;
  // Number of finite element nodes
  size_t nFEnds_;
  // Interpolation matrix from finite element nodes to quadrature points
  std::vector<std::vector<double>> interpMatrix_;
  // Reference position of blade nodes in local (blade root) frame
  std::vector<vs::Vector> refPosLoc_;
  std::vector<vs::Vector> refOrientLoc_;


  double t_start_; //When to start the mode oscillations

  double modeFreq_;

  bool mesh_motion_;

  int tStep_{0}; // Time step count

  double dt_{-1.0}; // Store nalu-wind step

  std::pair<double, double> naluTimer_{
    0.0, 0.0}; // store time taken in openfast calls

  int writeFreq_{
    30}; // Frequency to write line loads and deflections to netcdf file


};

} // namespace nalu

} // namespace sierra

#endif /* MODALSHAPEANALYSIS_H */
