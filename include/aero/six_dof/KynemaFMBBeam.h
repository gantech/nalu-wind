// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#ifndef KYNEMAFMBBEAM_H
#define KYNEMAFMBBEAM_H

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "aero/six_dof/KynemaFMBInterface.h"
#include "KokkosInterface.h"
#include "FieldTypeDef.h"
#include "aero/fsi/CalcLoadsAssembled.h"
#include "yaml-cpp/yaml.h"

#include <interfaces/blade/blade_interface.hpp>
#include <interfaces/components/beam_builder.hpp>

namespace sierra {
namespace kynema_ugf {

struct BeamSectionData
{
  double s{0.0};
  std::array<std::array<double, 6>, 6> mass_matrix{};
  std::array<std::array<double, 6>, 6> stiffness_matrix{};
};

template <typename Space>
struct BeamDataT
{
  using scalar_view_type = Kokkos::View<double*, Space>;
  using pos_view_type = Kokkos::View<double* [7], Space>;
  using load_view_type = Kokkos::View<double* [6], Space>;

  scalar_view_type bary_weights; // Barycentric weights for interpolation nodes
  scalar_view_type node_xi; // Non-dimensional location of each node along the beam
  pos_view_type pos;       // First 3 are position, last 4 are quaternion orientation
  pos_view_type disp;      // First 3 are translational displacement, last 4 are quaternion rotations
  load_view_type loads;     // First 3 are forces, last 3 are moments

  BeamDataT() = default;

  explicit BeamDataT(const std::size_t nNodes, const std::string& labelPrefix = "beam")
  {
    resize(nNodes, labelPrefix);
  }

  void resize(const std::size_t nNodes, const std::string& labelPrefix = "beam")
  {
    bary_weights = scalar_view_type(labelPrefix + "_bary_weights", nNodes);
    node_xi = scalar_view_type(labelPrefix + "_node_xi", nNodes);
    pos = pos_view_type(labelPrefix + "_pos", nNodes);
    disp = pos_view_type(labelPrefix + "_disp", nNodes);
    loads = load_view_type(labelPrefix + "_loads", nNodes);
  }

  std::size_t size() const { return node_xi.extent(0); }

};

using BeamDataDevice = BeamDataT<MemSpace>;
using BeamDataHost = BeamDataT<HostSpace>;
using BeamData = BeamDataDevice;

struct BeamBody
{
  bool use_restart_data = false;
  std::string restart_file_name = "beam.restart";
  std::string output_file_name = "";   // per-step text log (mirrors SixDof output_file_name)
  std::string output_file_path = "";   // kynema NetCDF output path
  bool write_output = false;

  bool enable_static_solve = false;
  double time_step = 0.01;
  double rho_inf{0.0};
  int number_of_nonlinear_iterations = 5;
  double absolute_error_tolerance = 1.0e-5;
  double relative_error_tolerance = 1.0e-3;

  int element_order = 10;
  int section_refinement = 0;
  bool prescribed_root_motion = false;
  size_t n_nodes = 0;

  BeamData beam_data;
  BeamDataHost beam_data_host;

  kynema_fmb::interfaces::components::BeamInput::QuadratureRule quadrature_rule =
    kynema_fmb::interfaces::components::BeamInput::QuadratureRule::GaussLobatto;
  kynema_fmb::interfaces::components::ReferenceAxisOrientation
    reference_axis_orientation =
      kynema_fmb::interfaces::components::ReferenceAxisOrientation::X;
  kynema_fmb::interfaces::components::ReferenceAxisOrientation
    section_orientation =
      kynema_fmb::interfaces::components::ReferenceAxisOrientation::X;

  std::vector<std::array<double, 3>> key_points_location;
  std::vector<double> key_points_twist;
  std::vector<double> key_points_s;
  std::vector<BeamSectionData> sections;

  std::array<double, 3> root_location = {0.0, 0.0, 0.0};
  std::array<double, 3> tip_location = {0.0, 0.0, 0.0};

  std::vector<std::string> forcing_surface_names;
  std::vector<std::string> moving_mesh_block_names;
  stk::mesh::PartVector forcing_surfaces;
  stk::mesh::PartVector moving_mesh_blocks;

  std::shared_ptr<kynema_fmb::interfaces::BladeInterface> kynema_interface =
    nullptr;
  std::shared_ptr<stk::mesh::BulkData> bulk = nullptr;
  GenericFieldType* total_force = nullptr;
  std::shared_ptr<CalcLoadsAssembled> calc_loads = nullptr;

};

class KynemaFMBBeam : public KynemaFMBInterface
{
public:
  explicit KynemaFMBBeam(const YAML::Node& node);
  ~KynemaFMBBeam() override = default;

  void
  setup(double dtKynemaUGF, std::shared_ptr<stk::mesh::BulkData> bulk) override;

  void initialize(int restartFreqKynemaUGF, double curTime) override;

  void map_displacements(double currentTime, bool updateCurCoor) override;

  void
  advance_struct_timestep(const double currentTime, const double dT) override;

  void map_loads(const double currentTime) override;

  stk::mesh::PartVector get_mesh_blocks() const override;

private:
  KynemaFMBBeam() = delete;
  KynemaFMBBeam(const KynemaFMBBeam&) = delete;

  void load(const YAML::Node& node);
  void load_beam(const YAML::Node& body_node, const YAML::Node& common_node);

  void setup_beam(
    BeamBody& beam,
    const double dtKynemaUGF,
    std::shared_ptr<stk::mesh::BulkData> bulk);

  void map_loads_beam(BeamBody& beam);
  void map_displacements_beam(BeamBody& beam, bool updateCur);
  void compute_mapping_beam(BeamBody& beam);

  kynema_fmb::interfaces::components::BeamInput::QuadratureRule
  parse_quadrature_rule(const std::string& name) const;

  kynema_fmb::interfaces::components::ReferenceAxisOrientation
  parse_orientation(const std::string& name) const;

  std::array<std::array<double, 6>, 6> parse_matrix6(const YAML::Node& node) const;

  std::shared_ptr<stk::mesh::BulkData> bulk_;
  std::array<double, 3> gravity_ = {0.0, 0.0, 0.0};
  int number_of_bodies_{0};
  int restart_frequency_{0};
  double dt_{-1.0};

  std::vector<BeamBody> beam_bodies_;
};

} // namespace kynema_ugf
} // namespace sierra

#endif