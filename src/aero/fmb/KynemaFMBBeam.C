// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#include "KynemaUGFEnv.h"
#include "aero/fmb/KynemaFMBBeam.h"
#include "aero/fsi/MapLoad.h"
#include <KynemaUGFParsing.h>

#include <interfaces/blade/blade_interface_builder.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>

#include "stk_mesh/base/GetNgpMesh.hpp"
#include "ngp_utils/NgpLoopUtils.h"
#include "ngp_utils/NgpTypes.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>

namespace sierra {
namespace kynema_ugf {

KynemaFMBBeam::KynemaFMBBeam(const YAML::Node& node)
{
  load(node);
}

kynema_fmb::interfaces::components::BeamInput::QuadratureRule
KynemaFMBBeam::parse_quadrature_rule(const std::string& name) const
{
  using Rule = kynema_fmb::interfaces::components::BeamInput::QuadratureRule;
  if (name == "GaussLegendre")
    return Rule::GaussLegendre;
  if (name == "GaussLobatto" || name == "Trapezoidal")
    return Rule::GaussLobatto;

  throw std::runtime_error(
    "Unsupported quadrature_rule for beam: " + name +
    ". Supported: GaussLegendre, GaussLobatto");
}

kynema_fmb::interfaces::components::ReferenceAxisOrientation
KynemaFMBBeam::parse_orientation(const std::string& name) const
{
  using Orientation =
    kynema_fmb::interfaces::components::ReferenceAxisOrientation;
  if (name == "X")
    return Orientation::X;
  if (name == "Z")
    return Orientation::Z;

  throw std::runtime_error(
    "Unsupported beam orientation: " + name + ". Supported: X, Z");
}

std::array<std::array<double, 6>, 6>
KynemaFMBBeam::parse_matrix6(const YAML::Node& node) const
{
  if (!node.IsSequence() || node.size() != 6)
    throw std::runtime_error("Matrix must contain 6 rows");

  std::array<std::array<double, 6>, 6> matrix{};
  for (std::size_t r = 0; r < 6; ++r) {
    if (!node[r].IsSequence() || node[r].size() != 6)
      throw std::runtime_error("Each matrix row must contain 6 values");
    for (std::size_t c = 0; c < 6; ++c)
      matrix[r][c] = node[r][c].as<double>();
  }
  return matrix;
}

void
KynemaFMBBeam::load(const YAML::Node& node)
{
  const int ndim = 3;
  get_required(node, "number_of_bodies", number_of_bodies_);

  if (node["gravity"]) {
    for (int idim = 0; idim < ndim; ++idim)
      gravity_[idim] = node["gravity"][idim].as<double>();
  }

  for (int ibody = 0; ibody < number_of_bodies_; ++ibody) {
    const auto body_key = "Body" + std::to_string(ibody);
    if (!node[body_key]) {
      throw std::runtime_error(
        "Node for " + body_key + " not present or correct in input file");
    }

    auto body_node = node[body_key];
    if (!body_node["type"]) {
      throw std::runtime_error("Beam body is missing required key: type");
    }

    const std::string body_type = body_node["type"].as<std::string>();
    if (body_type != "beam") {
      throw std::runtime_error(
        "unrecognized body type for beam interface. Expected type=beam");
    }

    load_beam(body_node, node);
  }
}

void
KynemaFMBBeam::load_beam(
  const YAML::Node& body_node,
  const YAML::Node& common_node)
{
  BeamBody new_body;
  BeamInterface new_iface;

  const YAML::Node empty_map(YAML::NodeType::Map);

  const YAML::Node common_config = common_node["config"];
  const YAML::Node common_solution =
    common_config ? common_config["solution"] : empty_map;
  const YAML::Node common_outputs =
    common_config ? common_config["outputs"] : empty_map;
  const YAML::Node common_blade =
    common_config ? common_config["blade"] : empty_map;

  const YAML::Node body_config = body_node["config"];
  const YAML::Node body_outputs =
    body_config ? body_config["outputs"] : empty_map;
  const YAML::Node body_blade = body_config ? body_config["blade"] : empty_map;

  // Keep high-level solution config common across beams.
  get_if_present(common_solution, "enable_static_solve", new_iface.enable_static_solve, false);
  get_required(common_solution, "time_step", new_iface.time_step);
  get_if_present(common_solution, "damping_factor", new_iface.rho_inf, 0.0);
  get_if_present(
    common_solution,
    "max_nonlinear_iterations",
    new_iface.number_of_nonlinear_iterations,
    5);
  get_if_present(
    common_solution,
    "absolute_error_tolerance",
    new_iface.absolute_error_tolerance,
    1.0e-5);
  get_if_present(
    common_solution,
    "relative_error_tolerance",
    new_iface.relative_error_tolerance,
    1.0e-3);

  // Allow per-beam output path override.
  get_if_present(common_outputs, "write_output", new_iface.write_output, false);
  get_if_present(common_outputs, "output_file_path", new_iface.output_file_path, std::string(""));
  get_if_present_no_default(body_outputs, "write_output", new_iface.write_output);
  get_if_present_no_default(body_outputs, "output_file_path", new_iface.output_file_path);

  get_if_present(common_blade, "element_order", new_iface.element_order, 10);
  get_if_present(common_blade, "section_refinement", new_iface.section_refinement, 0);
  get_if_present(
    common_blade,
    "prescribed_root_motion",
    new_iface.prescribed_root_motion,
    false);

  std::string quadrature_rule_name("GaussLobatto");
  get_if_present(common_blade, "quadrature_rule", quadrature_rule_name, std::string("GaussLobatto"));
  new_iface.quadrature_rule = parse_quadrature_rule(quadrature_rule_name);

  // Allow per-beam orientation overrides.
  std::string reference_axis_orientation_name("X");
  get_if_present(
    common_blade,
    "reference_axis_orientation",
    reference_axis_orientation_name,
    std::string("X"));
  get_if_present_no_default(
    body_blade,
    "reference_axis_orientation",
    reference_axis_orientation_name);
  new_iface.reference_axis_orientation =
    parse_orientation(reference_axis_orientation_name);

  std::string section_orientation_name("X");
  get_if_present(
    common_blade,
    "section_orientation",
    section_orientation_name,
    std::string("X"));
  get_if_present_no_default(
    body_blade,
    "section_orientation",
    section_orientation_name);
  new_iface.section_orientation = parse_orientation(section_orientation_name);

  get_if_present_no_default(body_node, "output_file_name", new_iface.output_file_name);
  get_if_present_no_default(body_node, "use_restart_data", new_iface.use_restart_data);
  get_if_present_no_default(body_node, "restart_file_name", new_iface.restart_file_name);

  const YAML::Node key_points_location_node = body_node["key_points"]["location"];
  const YAML::Node key_points_twist_node = body_node["key_points"]["twist"];
  const YAML::Node section_s_node = body_node["beam_properties"]["stations"];
  const YAML::Node mass_matrices_node =
    body_node["beam_properties"]["mass_matrices"];
  const YAML::Node stiffness_matrices_node =
    body_node["beam_properties"]["stiffness_matrices"];

  if (
    !key_points_location_node || !key_points_twist_node || !section_s_node ||
    !mass_matrices_node || !stiffness_matrices_node) {
    throw std::runtime_error(
      "Beam body is missing one or more required keys: key_points/beam_properties");
  }

  if (key_points_location_node.size() < 2) {
    throw std::runtime_error(
      "Need at least two key-point locations for each beam body");
  }
  if (key_points_location_node.size() != key_points_twist_node.size()) {
    throw std::runtime_error(
      "key_points.location and key_points.twist size mismatch");
  }

  new_iface.key_points_location.reserve(key_points_location_node.size());
  new_iface.key_points_twist.reserve(key_points_twist_node.size());

  for (const auto& location_node : key_points_location_node) {
    if (!location_node.IsSequence() || location_node.size() != 3) {
      throw std::runtime_error(
        "Each key point location must be a 3-value sequence");
    }
    new_iface.key_points_location.push_back(
      {location_node[0].as<double>(), location_node[1].as<double>(),
       location_node[2].as<double>()});
  }

  for (const auto& twist_node : key_points_twist_node)
    new_iface.key_points_twist.push_back(twist_node.as<double>());

  new_iface.key_points_s.assign(new_iface.key_points_location.size(), 0.0);
  double cumulative_length = 0.0;
  for (std::size_t i = 1; i < new_iface.key_points_location.size(); ++i) {
    const auto dx =
      new_iface.key_points_location[i][0] - new_iface.key_points_location[i - 1][0];
    const auto dy =
      new_iface.key_points_location[i][1] - new_iface.key_points_location[i - 1][1];
    const auto dz =
      new_iface.key_points_location[i][2] - new_iface.key_points_location[i - 1][2];
    cumulative_length += std::sqrt(dx * dx + dy * dy + dz * dz);
    new_iface.key_points_s[i] = cumulative_length;
  }
  if (cumulative_length > 0.0) {
    for (auto& s : new_iface.key_points_s)
      s /= cumulative_length;
  }

  if (
    section_s_node.size() != mass_matrices_node.size() ||
    section_s_node.size() != stiffness_matrices_node.size()) {
    throw std::runtime_error(
      "beam_properties stations/mass_matrices/stiffness_matrices size mismatch");
  }

  new_iface.sections.reserve(section_s_node.size());
  for (std::size_t i = 0; i < section_s_node.size(); ++i) {
    BeamSectionData section;
    section.s = section_s_node[i].as<double>();
    section.mass_matrix = parse_matrix6(mass_matrices_node[i]);
    section.stiffness_matrix = parse_matrix6(stiffness_matrices_node[i]);
    new_iface.sections.push_back(section);
  }

  new_iface.root_location = new_iface.key_points_location.front();
  new_iface.tip_location = new_iface.key_points_location.back();

  if (body_node["forcing_surfaces"]) {
    for (std::size_t isurf = 0; isurf < body_node["forcing_surfaces"].size();
         ++isurf) {
      new_body.forcing_surface_names.emplace_back(
        body_node["forcing_surfaces"][isurf].as<std::string>());
    }
  }

  if (body_node["moving_mesh_blocks"]) {
    for (std::size_t iblock = 0; iblock < body_node["moving_mesh_blocks"].size();
         ++iblock) {
      new_body.moving_mesh_block_names.emplace_back(
        body_node["moving_mesh_blocks"][iblock].as<std::string>());
    }
  }

  beam_bodies_.emplace_back(std::move(new_body));
  beam_interfaces_.emplace_back(std::move(new_iface));
}

void
KynemaFMBBeam::setup_beam(
  BeamBody& beam,
  BeamInterface& iface,
  const double dtKynemaUGF,
  std::shared_ptr<stk::mesh::BulkData> bulk)
{
  (void)dtKynemaUGF;

  kynema_fmb::interfaces::BladeInterfaceBuilder builder;
  auto& solution = builder.Solution();
  if (iface.enable_static_solve)
    solution.EnableStaticSolve();
  else
    solution.EnableDynamicSolve();

  // solution.SetGravity(gravity_)
  solution
    .SetTimeStep(iface.time_step)
    .SetDampingFactor(iface.rho_inf)
    .SetMaximumNonlinearIterations(iface.number_of_nonlinear_iterations)
    .SetAbsoluteErrorTolerance(iface.absolute_error_tolerance)
    .SetRelativeErrorTolerance(iface.relative_error_tolerance);

  if (iface.write_output && !iface.output_file_path.empty()) {
    builder.Outputs().SetOutputFilePath(iface.output_file_path);
  }

  builder.Blade()
    .SetElementOrder(static_cast<std::size_t>(iface.element_order))
    .SetQuadratureRule(iface.quadrature_rule)
    .SetSectionRefinement(static_cast<std::size_t>(iface.section_refinement))
    .PrescribedRootMotion(iface.prescribed_root_motion);

  for (std::size_t i = 0; i < iface.key_points_s.size(); ++i)
    builder.Blade().AddRefAxisTwist(iface.key_points_s[i], iface.key_points_twist[i]);

  for (std::size_t i = 0; i < iface.key_points_s.size(); ++i)
    builder.Blade().AddRefAxisPoint(
      iface.key_points_s[i], iface.key_points_location[i], iface.reference_axis_orientation);

  for (const auto& section : iface.sections)
    builder.Blade().AddSection(
      section.s, section.mass_matrix, section.stiffness_matrix, iface.section_orientation);

  std::array<double, 7> root_position = {0.0, 0.0, 0.0, 0.70710678, 0, -0.70710678, 0};
  builder.Blade().SetRootPosition(root_position);

  iface.kynema_interface =
    std::make_shared<kynema_fmb::interfaces::BladeInterface>(builder.Build());
  beam.n_nodes = iface.kynema_interface->Blade().nodes.size();
  beam.bulk = bulk;
  beam.beam_data.resize(beam.n_nodes);
  beam.beam_data_host.resize(beam.n_nodes);

  auto& meta = bulk->mesh_meta_data();
  beam.total_force = meta.get_field<double>(stk::topology::NODE_RANK, "tforce");
  if (beam.total_force == nullptr)
    beam.total_force = &(meta.declare_field<double>(stk::topology::NODE_RANK, "tforce"));

  for (const auto& surface_name : beam.forcing_surface_names) {
    stk::mesh::Part* part = meta.get_part(surface_name);
    if (part == nullptr)
      throw std::runtime_error("Unable to find forcing surface part: " + surface_name);
    beam.forcing_surfaces.push_back(part);
    stk::mesh::put_field_on_mesh(
      *beam.total_force, *part, meta.spatial_dimension(), nullptr);
  }

  auto beam_xi = &(meta.declare_field<double>(stk::topology::NODE_RANK, "beam_xi"));

  for (const auto& block_name : beam.moving_mesh_block_names) {
    stk::mesh::Part* part = meta.get_part(block_name);
    if (part == nullptr)
      throw std::runtime_error("Unable to find moving mesh block part: " + block_name);
    beam.moving_mesh_blocks.push_back(part);
    stk::mesh::put_field_on_mesh(*beam_xi, *part, 1, nullptr);
  }

  beam.calc_loads = std::make_shared<CalcLoadsAssembled>(beam.forcing_surfaces);
  beam.calc_loads->setup(bulk);
}

void
KynemaFMBBeam::setup(
  double dtKynemaUGF,
  std::shared_ptr<stk::mesh::BulkData> bulk)
{
  bulk_ = bulk;
  dt_ = dtKynemaUGF;
  for (std::size_t ibeam = 0; ibeam < beam_bodies_.size(); ++ibeam)
    setup_beam(beam_bodies_[ibeam], beam_interfaces_[ibeam], dtKynemaUGF, bulk);
}

void
KynemaFMBBeam::initialize(int restartFreqKynemaUGF, double curTime)
{
  restart_frequency_ = restartFreqKynemaUGF;

  for (std::size_t ibeam = 0; ibeam < beam_bodies_.size(); ++ibeam) {
    auto& beam = beam_bodies_[ibeam];
    auto& iface = beam_interfaces_[ibeam];
    if (!iface.kynema_interface)
      continue;

    if (iface.use_restart_data) {
      KynemaUGFEnv::self().kynema_ugfOutputP0()
        << "Beam restart read is not yet implemented for beam body " << ibeam
        << std::endl;
    }

    get_displacements(beam, iface);
    compute_mapping_beam(beam);

    auto& tip_node = iface.kynema_interface->Blade().nodes.back();
    tip_node.loads[0] = 0.0005;
    tip_node.loads[1] = 0.0005;
    tip_node.loads[2] = 0.0;
    tip_node.loads[3] = 0.0;
    tip_node.loads[4] = 0.0;
    tip_node.loads[5] = 0.0;

    if (iface.enable_static_solve) {
      const bool converged = iface.kynema_interface->Step();
      if (!converged) {
        KynemaUGFEnv::self().kynema_ugfOutputP0()
          << "Initial static solve did not converge for beam body " << ibeam
          << std::endl;
      }
    }
  }

  map_displacements(curTime, false);

  if (curTime < 1e-10 && bulk_ != nullptr) {
    auto& meta = bulk_->mesh_meta_data();
    const VectorFieldType* meshDisp =
      meta.get_field<double>(stk::topology::NODE_RANK, "mesh_displacement");
    const VectorFieldType* meshVel =
      meta.get_field<double>(stk::topology::NODE_RANK, "mesh_velocity");

    if (meshDisp == nullptr || meshVel == nullptr)
      return;

    const VectorFieldType* meshDispNp1 =
      &(meshDisp->field_of_state(stk::mesh::StateNP1));
    VectorFieldType* meshDispN = &(meshDisp->field_of_state(stk::mesh::StateN));
    VectorFieldType* meshDispNm1 =
      &(meshDisp->field_of_state(stk::mesh::StateNM1));
    const VectorFieldType* meshVelNp1 =
      &(meshVel->field_of_state(stk::mesh::StateNP1));

    meshDisp->sync_to_host();
    meshVel->sync_to_host();
    meshDispNp1->sync_to_host();
    meshDispN->sync_to_host();
    meshDispNm1->sync_to_host();
    meshVelNp1->sync_to_host();

    stk::mesh::Selector sel = meta.universal_part();
    const auto& bkts = bulk_->get_buckets(stk::topology::NODE_RANK, sel);
    for (const auto* b : bkts) {
      for (const auto node : *b) {
        const double* velNp1 = stk::mesh::field_data(*meshVelNp1, node);
        const double* dispNp1 = stk::mesh::field_data(*meshDispNp1, node);
        double* dispN = stk::mesh::field_data(*meshDispN, node);
        double* dispNm1 = stk::mesh::field_data(*meshDispNm1, node);
        for (std::size_t i = 0; i < 3; ++i) {
          dispN[i] = dispNp1[i] - dt_ * velNp1[i];
          dispNm1[i] = dispNp1[i] - 2.0 * dt_ * velNp1[i];
        }
      }
    }
    meshDispN->modify_on_host();
    meshDispNm1->modify_on_host();
  }
}


void
KynemaFMBBeam::get_positions(BeamBody& beam, BeamInterface& iface) 
{

  if (!iface.kynema_interface)
    return;

  auto node_xi_vec = iface.kynema_interface->Blade().node_xi;

  for (std::size_t i = 0; i < node_xi_vec.size(); ++i) {
    beam.beam_data_host.node_xi(i) = node_xi_vec[i];
    for (std::size_t j = 0; j < 7; ++j) {
      beam.beam_data_host.pos(i, j) = iface.kynema_interface->Blade().nodes[i].position[j];
      KynemaUGFEnv::self().kynema_ugfOutputP0() << "nodes[" << i << "].position[" << j << "] = " << iface.kynema_interface->Blade().nodes[i].position[j] << std::endl;
    }
  }
  
}


void
KynemaFMBBeam::get_displacements(BeamBody& beam, BeamInterface& iface) 
{

  if (!iface.kynema_interface || beam.moving_mesh_blocks.empty())
    return;

  constexpr int kMaxInterpolationNodes = 32;
  auto n_nodes = beam.n_nodes;
  auto& nodes = iface.kynema_interface->Blade().nodes;
  for (std::size_t i = 0; i < n_nodes; ++i) {
    for (std::size_t j = 0; j < 7; ++j) {
      beam.beam_data_host.disp(i, j) = nodes[i].displacement[j];
    }
    for (std::size_t j = 0; j < 6; ++j)
      beam.beam_data_host.vel(i, j) = nodes[i].velocity[j];
  }

}



void
KynemaFMBBeam::map_displacements(double currentTime, bool updateCurCoor)
{
  (void)currentTime;
  for (std::size_t ibeam = 0; ibeam < beam_bodies_.size(); ++ibeam) {
    get_displacements(beam_bodies_[ibeam], beam_interfaces_[ibeam]);
    map_displacements_beam(beam_bodies_[ibeam], updateCurCoor);
  }
}

void
KynemaFMBBeam::map_loads()
{
  for (std::size_t ibeam = 0; ibeam < beam_bodies_.size(); ++ibeam) {
    auto& beam = beam_bodies_[ibeam];
    auto& iface = beam_interfaces_[ibeam];

    map_loads_beam(beam);

    iface.loads[BeamInterface::lts::NP1].resize(beam.n_nodes);
    for (std::size_t inode = 0; inode < beam.n_nodes; ++inode) {
      for (int idim = 0; idim < 6; ++idim)
        iface.loads[BeamInterface::lts::NP1][inode][idim] =
          beam.beam_data_host.loads(inode, idim);
    }
  }
}

void
KynemaFMBBeam::advance_struct_timestep(const double currentTime, const double dT)
{
  for (std::size_t ibeam = 0; ibeam < beam_bodies_.size(); ++ibeam) {
    auto& beam = beam_bodies_[ibeam];
    auto& iface = beam_interfaces_[ibeam];
    if (!iface.kynema_interface)
      continue;

    (void)dT;
    const bool converged = iface.kynema_interface->Step();


    if (!converged) {
      KynemaUGFEnv::self().kynema_ugfOutputP0()
        << "Kynema did not converge for beam body " << ibeam << std::endl;
    }

    if (
      !iface.output_file_name.empty() &&
      KynemaUGFEnv::self().parallel_rank() == 0) {
      auto& nodes = iface.kynema_interface->Blade().nodes;
      if (!nodes.empty()) {
        std::ofstream outfile(iface.output_file_name, std::ios::app);
        const auto& root = nodes.front();
        const auto& tip = nodes.back();
        outfile << currentTime;
        for (int i = 0; i < 7; ++i)
          outfile << " " << root.position[i];
        for (int i = 0; i < 7; ++i)
          outfile << " " << tip.position[i];
        outfile << std::endl;
      }
    }
  }
}

void
KynemaFMBBeam::predict_loads() 
{
  for (std::size_t ibeam = 0; ibeam < beam_bodies_.size(); ++ibeam) {
    auto& beam = beam_bodies_[ibeam];
    auto& iface = beam_interfaces_[ibeam];
    const auto n_nodes = beam.n_nodes;

    iface.loads[BeamInterface::lts::NP1].resize(n_nodes);
    iface.loads[BeamInterface::lts::N].resize(n_nodes);
    iface.loads[BeamInterface::lts::NM1].resize(n_nodes);

    for (std::size_t inode = 0; inode < n_nodes; ++inode) {
      for (int idim = 0; idim < 6; ++idim) {
        iface.loads[BeamInterface::lts::NP1][inode][idim] =
          2.0 * iface.loads[BeamInterface::lts::N][inode][idim] -
          iface.loads[BeamInterface::lts::NM1][inode][idim];
      }
    }
  }
}

void 
KynemaFMBBeam::send_loads(const double /* currentTime*/)
{
  for (std::size_t ibeam = 0; ibeam < beam_bodies_.size(); ++ibeam) {
    auto& beam = beam_bodies_[ibeam];
    auto& iface = beam_interfaces_[ibeam];
    if (!iface.kynema_interface)
      continue;

    auto& nodes = iface.kynema_interface->Blade().nodes;
    const auto n_nodes = beam.n_nodes;

    iface.loads[BeamInterface::lts::NP1].resize(n_nodes);

    for (std::size_t inode = 0; inode < n_nodes; ++inode) {
      for (int idim = 0; idim < 6; ++idim)
        nodes[inode].loads[idim] = iface.loads[BeamInterface::lts::NP1][inode][idim];
    }
  }
}

void 
KynemaFMBBeam::finalize_struct_timestep()
{
  for (std::size_t ibeam = 0; ibeam < beam_bodies_.size(); ++ibeam) {
    auto& beam = beam_bodies_[ibeam];
    auto& iface = beam_interfaces_[ibeam];
    const auto n_nodes = beam.n_nodes;

    iface.loads[BeamInterface::lts::NP1].resize(n_nodes);
    iface.loads[BeamInterface::lts::N].resize(n_nodes);
    iface.loads[BeamInterface::lts::NM1].resize(n_nodes);

    for (std::size_t inode = 0; inode < n_nodes; ++inode) {
      for (int idim = 0; idim < 6; ++idim) {
        iface.loads[BeamInterface::lts::NM1][inode][idim] =
          iface.loads[BeamInterface::lts::N][inode][idim];
        iface.loads[BeamInterface::lts::N][inode][idim] =
          iface.loads[BeamInterface::lts::NP1][inode][idim];
      }
    }
  }
}

stk::mesh::PartVector
KynemaFMBBeam::get_mesh_blocks() const
{
  stk::mesh::PartVector all_mesh_blocks;
  for (const auto& beam : beam_bodies_) {
    for (const auto& block : beam.moving_mesh_blocks)
      all_mesh_blocks.push_back(block);
  }
  return all_mesh_blocks;
}

} // namespace kynema_ugf
} // namespace sierra
