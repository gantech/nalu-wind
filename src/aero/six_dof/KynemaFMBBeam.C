// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#include "KynemaUGFEnv.h"
#include "aero/six_dof/KynemaFMBBeam.h"
#include "aero/six_dof/KynemaFMBBeamUtils.h"
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
  BeamBody beam;

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
  get_if_present(common_solution, "enable_static_solve", beam.enable_static_solve, false);
  get_required(common_solution, "time_step", beam.time_step);
  get_if_present(common_solution, "damping_factor", beam.rho_inf, 0.0);
  get_if_present(
    common_solution,
    "max_nonlinear_iterations",
    beam.number_of_nonlinear_iterations,
    5);
  get_if_present(
    common_solution,
    "absolute_error_tolerance",
    beam.absolute_error_tolerance,
    1.0e-5);
  get_if_present(
    common_solution,
    "relative_error_tolerance",
    beam.relative_error_tolerance,
    1.0e-3);

  // Allow per-beam output path override.
  get_if_present(common_outputs, "write_output", beam.write_output, false);
  get_if_present(common_outputs, "output_file_path", beam.output_file_path, std::string(""));
  get_if_present_no_default(body_outputs, "write_output", beam.write_output);
  get_if_present_no_default(body_outputs, "output_file_path", beam.output_file_path);

  get_if_present(common_blade, "element_order", beam.element_order, 10);
  get_if_present(common_blade, "section_refinement", beam.section_refinement, 0);
  get_if_present(
    common_blade,
    "prescribed_root_motion",
    beam.prescribed_root_motion,
    false);

  std::string quadrature_rule_name("GaussLobatto");
  get_if_present(common_blade, "quadrature_rule", quadrature_rule_name, std::string("GaussLobatto"));
  beam.quadrature_rule = parse_quadrature_rule(quadrature_rule_name);

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
  beam.reference_axis_orientation =
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
  beam.section_orientation = parse_orientation(section_orientation_name);

  get_if_present_no_default(body_node, "output_file_name", beam.output_file_name);
  get_if_present_no_default(body_node, "use_restart_data", beam.use_restart_data);
  get_if_present_no_default(body_node, "restart_file_name", beam.restart_file_name);

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

  beam.key_points_location.reserve(key_points_location_node.size());
  beam.key_points_twist.reserve(key_points_twist_node.size());

  for (const auto& location_node : key_points_location_node) {
    if (!location_node.IsSequence() || location_node.size() != 3) {
      throw std::runtime_error(
        "Each key point location must be a 3-value sequence");
    }
    beam.key_points_location.push_back(
      {location_node[0].as<double>(), location_node[1].as<double>(),
       location_node[2].as<double>()});
  }

  for (const auto& twist_node : key_points_twist_node)
    beam.key_points_twist.push_back(twist_node.as<double>());

  beam.key_points_s.assign(beam.key_points_location.size(), 0.0);
  double cumulative_length = 0.0;
  for (std::size_t i = 1; i < beam.key_points_location.size(); ++i) {
    const auto dx =
      beam.key_points_location[i][0] - beam.key_points_location[i - 1][0];
    const auto dy =
      beam.key_points_location[i][1] - beam.key_points_location[i - 1][1];
    const auto dz =
      beam.key_points_location[i][2] - beam.key_points_location[i - 1][2];
    cumulative_length += std::sqrt(dx * dx + dy * dy + dz * dz);
    beam.key_points_s[i] = cumulative_length;
  }
  if (cumulative_length > 0.0) {
    for (auto& s : beam.key_points_s)
      s /= cumulative_length;
  }

  if (
    section_s_node.size() != mass_matrices_node.size() ||
    section_s_node.size() != stiffness_matrices_node.size()) {
    throw std::runtime_error(
      "beam_properties stations/mass_matrices/stiffness_matrices size mismatch");
  }

  beam.sections.reserve(section_s_node.size());
  for (std::size_t i = 0; i < section_s_node.size(); ++i) {
    BeamSectionData section;
    section.s = section_s_node[i].as<double>();
    section.mass_matrix = parse_matrix6(mass_matrices_node[i]);
    section.stiffness_matrix = parse_matrix6(stiffness_matrices_node[i]);
    beam.sections.push_back(section);
  }

  beam.root_location = beam.key_points_location.front();
  beam.tip_location = beam.key_points_location.back();

  if (body_node["forcing_surfaces"]) {
    for (std::size_t isurf = 0; isurf < body_node["forcing_surfaces"].size();
         ++isurf) {
      beam.forcing_surface_names.emplace_back(
        body_node["forcing_surfaces"][isurf].as<std::string>());
    }
  }

  if (body_node["moving_mesh_blocks"]) {
    for (std::size_t iblock = 0; iblock < body_node["moving_mesh_blocks"].size();
         ++iblock) {
      beam.moving_mesh_block_names.emplace_back(
        body_node["moving_mesh_blocks"][iblock].as<std::string>());
    }
  }

  beam_bodies_.emplace_back(std::move(beam));
}

void
KynemaFMBBeam::setup_beam(
  BeamBody& beam,
  const double dtKynemaUGF,
  std::shared_ptr<stk::mesh::BulkData> bulk)
{
  (void)dtKynemaUGF;

  kynema_fmb::interfaces::BladeInterfaceBuilder builder;
  auto& solution = builder.Solution();
  if (beam.enable_static_solve)
    solution.EnableStaticSolve();
  else
    solution.EnableDynamicSolve();

  // solution.SetGravity(gravity_)
  solution
    .SetTimeStep(beam.time_step)
    .SetDampingFactor(beam.rho_inf)
    .SetMaximumNonlinearIterations(beam.number_of_nonlinear_iterations)
    .SetAbsoluteErrorTolerance(beam.absolute_error_tolerance)
    .SetRelativeErrorTolerance(beam.relative_error_tolerance);

  if (beam.write_output && !beam.output_file_path.empty()) {
    builder.Outputs().SetOutputFilePath(beam.output_file_path);
  }

  builder.Blade()
    .SetElementOrder(static_cast<std::size_t>(beam.element_order))
    .SetQuadratureRule(beam.quadrature_rule)
    .SetSectionRefinement(static_cast<std::size_t>(beam.section_refinement))
    .PrescribedRootMotion(beam.prescribed_root_motion);

  for (std::size_t i = 0; i < beam.key_points_s.size(); ++i)
    builder.Blade().AddRefAxisTwist(beam.key_points_s[i], beam.key_points_twist[i]);

  for (std::size_t i = 0; i < beam.key_points_s.size(); ++i)
    builder.Blade().AddRefAxisPoint(
      beam.key_points_s[i], beam.key_points_location[i], beam.reference_axis_orientation);

  for (const auto& section : beam.sections)
    builder.Blade().AddSection(
      section.s, section.mass_matrix, section.stiffness_matrix, beam.section_orientation);

  std::array<double, 7> root_position = {0.0, 0.0, 0.0, 0.70710678, 0, -0.70710678, 0};
  builder.Blade().SetRootPosition(root_position);

  beam.kynema_interface =
    std::make_shared<kynema_fmb::interfaces::BladeInterface>(builder.Build());
  beam.n_nodes = beam.kynema_interface->Blade().nodes.size();
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
  for (auto& beam : beam_bodies_)
    setup_beam(beam, dtKynemaUGF, bulk);
}

void
KynemaFMBBeam::initialize(int restartFreqKynemaUGF, double curTime)
{
  restart_frequency_ = restartFreqKynemaUGF;

  for (std::size_t ibeam = 0; ibeam < beam_bodies_.size(); ++ibeam) {
    auto& beam = beam_bodies_[ibeam];
    if (!beam.kynema_interface)
      continue;

    if (beam.use_restart_data) {
      KynemaUGFEnv::self().kynema_ugfOutputP0()
        << "Beam restart read is not yet implemented for beam body " << ibeam
        << std::endl;
    }

    compute_mapping_beam(beam);

    auto& tip_node = beam.kynema_interface->Blade().nodes.back();
    tip_node.loads[0] = 0.0005;
    tip_node.loads[1] = 0.0005;
    tip_node.loads[2] = 0.0;
    tip_node.loads[3] = 0.0;
    tip_node.loads[4] = 0.0;
    tip_node.loads[5] = 0.0;

    if (beam.enable_static_solve) {
      const bool converged = beam.kynema_interface->Step();
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
KynemaFMBBeam::compute_mapping_beam(BeamBody& beam)
{
  if (!beam.kynema_interface)
    return;

  auto node_xi_vec = beam.kynema_interface->Blade().node_xi;

  for (std::size_t i = 0; i < node_xi_vec.size(); ++i) {
    beam.beam_data_host.node_xi(i) = node_xi_vec[i];
    for (std::size_t j = 0; j < 7; ++j) {
      beam.beam_data_host.pos(i, j) = beam.kynema_interface->Blade().nodes[i].position[j];
      KynemaUGFEnv::self().kynema_ugfOutputP0() << "nodes[" << i << "].position[" << j << "] = " << beam.kynema_interface->Blade().nodes[i].position[j] << std::endl;
    }

  }
  Kokkos::deep_copy(beam.beam_data.node_xi, beam.beam_data_host.node_xi);
  Kokkos::deep_copy(beam.beam_data.pos, beam.beam_data_host.pos);

  constexpr int kMaxInterpolationNodes = 32;
  auto node_xi_host = beam.beam_data_host.node_xi;
  auto n_nodes = beam.n_nodes;
  // Precompute barycentric weights for interpolation on xi-space nodes
  ComputeBarycentricWeights(
      node_xi_host.data(), static_cast<int>(n_nodes), beam.beam_data_host.bary_weights.data()
  );
  Kokkos::deep_copy(beam.beam_data.bary_weights, beam.beam_data_host.bary_weights);

  auto& meta = beam.bulk->mesh_meta_data();
  const auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*(beam.bulk));
  const stk::mesh::EntityRank entityRank = stk::topology::NODE_RANK;

  // get the parts in the current motion frame
  stk::mesh::Selector sel =
    stk::mesh::selectUnion(beam.moving_mesh_blocks) &
    (meta.locally_owned_part() | meta.globally_shared_part());
  // get the field from the NGP mesh
  stk::mesh::NgpField<double> modelCoords =
    stk::mesh::get_updated_ngp_field<double>(
      *meta.get_field<double>(entityRank, "coordinates"));
  stk::mesh::NgpField<double> beam_xi =
    stk::mesh::get_updated_ngp_field<double>(
      *meta.get_field<double>(entityRank, "beam_xi"));
  modelCoords.sync_to_device();
  beam_xi.sync_to_device();

  auto positions = beam.beam_data.pos;
  auto node_xi = beam.beam_data.node_xi;
  auto bary_weights = beam.beam_data.bary_weights;

  kynema_ugf_ngp::run_entity_algorithm(
    "KynemaFMBBeam_compute_mapping", ngpMesh, entityRank, sel,
    KOKKOS_LAMBDA(
      const kynema_ugf_ngp::NGPMeshTraits<stk::mesh::NgpMesh>::MeshIndex& mi) {

      const double query_point[3] = {modelCoords(mi, 0), modelCoords(mi, 1), modelCoords(mi, 2)};
      double scratch_weights[kMaxInterpolationNodes];
      double scratch_dweights[kMaxInterpolationNodes];
      double closest_position[3] = {0.0, 0.0, 0.0};
      double xi = 0.0;
      double dist2 = 0.0;

      FindClosestPointOnBlade(
          query_point,
          node_xi.data(),
          positions.data(),
          n_nodes,
          bary_weights.data(),
          scratch_weights,
          scratch_dweights,
          xi,
          closest_position,
          dist2
      );

      beam_xi.get(mi,0) = xi;
    });

  beam_xi.modify_on_device();

}

void
KynemaFMBBeam::map_displacements_beam(BeamBody& beam, bool updateCur)
{
  if (!beam.kynema_interface || beam.moving_mesh_blocks.empty())
    return;


  constexpr int kMaxInterpolationNodes = 32;
  auto n_nodes = beam.n_nodes;
  auto& nodes = beam.kynema_interface->Blade().nodes;
  for (std::size_t i = 0; i < n_nodes; ++i) {
    for (std::size_t j = 0; j < 7; ++j) {
      // TODO: Also get velocity
      beam.beam_data_host.disp(i, j) = nodes[i].displacement[j];
    }
  }
  Kokkos::deep_copy(beam.beam_data.disp, beam.beam_data_host.disp);


  auto& meta = beam.bulk->mesh_meta_data();
  const auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*(beam.bulk));
  const stk::mesh::EntityRank entityRank = stk::topology::NODE_RANK;

  // get the parts in the current motion frame
  stk::mesh::Selector sel =
    stk::mesh::selectUnion(beam.moving_mesh_blocks) &
    (meta.locally_owned_part() | meta.globally_shared_part());
  // get the field from the NGP mesh
  stk::mesh::NgpField<double> modelCoords =
    stk::mesh::get_updated_ngp_field<double>(
      *meta.get_field<double>(entityRank, "coordinates"));
  stk::mesh::NgpField<double> curCoords =
    stk::mesh::get_updated_ngp_field<double>(
      *meta.get_field<double>(entityRank, "current_coordinates"));
  stk::mesh::NgpField<double> meshDisp =
    stk::mesh::get_updated_ngp_field<double>(
      *meta.get_field<double>(entityRank, "mesh_displacement"));
  stk::mesh::NgpField<double> meshVel =
    stk::mesh::get_updated_ngp_field<double>(
      *meta.get_field<double>(entityRank, "mesh_velocity"));
  stk::mesh::NgpField<double> beam_xi =
    stk::mesh::get_updated_ngp_field<double>(
      *meta.get_field<double>(entityRank, "beam_xi"));

  modelCoords.sync_to_device();
  curCoords.sync_to_device();
  meshDisp.sync_to_device();
  meshVel.sync_to_device();
  beam_xi.sync_to_device();

  auto displacements = beam.beam_data.disp;
  auto node_xi = beam.beam_data.node_xi;
  auto positions = beam.beam_data.pos;
  auto bary_weights = beam.beam_data.bary_weights;

  kynema_ugf_ngp::run_entity_algorithm(
    "KynemaFMBBeam_map_displacements", ngpMesh, entityRank, sel,
    KOKKOS_LAMBDA(
      const kynema_ugf_ngp::NGPMeshTraits<stk::mesh::NgpMesh>::MeshIndex& mi) {

      double scratch_weights[kMaxInterpolationNodes];
      double pos[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      double disp[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      const double xi = beam_xi.get(mi, 0);
      double qp[3] = {modelCoords(mi, 0), modelCoords(mi, 1), modelCoords(mi, 2)};

      InterpolateFieldAtPoint(
          xi,
          node_xi.data(),
          positions.data(),
          n_nodes,
          7,
          bary_weights.data(),
          scratch_weights,
          pos
      );
      InterpolateFieldAtPoint(
          xi,
          node_xi.data(),
          displacements.data(),
          n_nodes,
          7,
          bary_weights.data(),
          scratch_weights,
          disp
      );

      double rel_pos_g[3] = {qp[0] - pos[0], qp[1] - pos[1], qp[2] - pos[2]};
      double rel_pos_l[3] = {0.0, 0.0, 0.0};
      RotateVectorByQuaternionInv(pos + 3, rel_pos_g, rel_pos_l);

      double tmp_disp[3] = {0.0, 0.0, 0.0};
      double rot_disp[3] = {0.0, 0.0, 0.0};
      RotateVectorByQuaternion(disp + 3, rel_pos_l, tmp_disp);
      RotateVectorByQuaternion(pos + 3, tmp_disp, rot_disp);

      meshDisp.get(mi, 0) = disp[0] + rot_disp[0] - rel_pos_g[0];
      meshDisp.get(mi, 1) = disp[1] + rot_disp[1] - rel_pos_g[1];
      meshDisp.get(mi, 2) = disp[2] + rot_disp[2] - rel_pos_g[2];

      if (updateCur) {
        curCoords.get(mi, 0) = modelCoords(mi, 0) + meshDisp.get(mi, 0);
        curCoords.get(mi, 1) = modelCoords(mi, 1) + meshDisp.get(mi, 1);
        curCoords.get(mi, 2) = modelCoords(mi, 2) + meshDisp.get(mi, 2);
      }
    });

  meshDisp.modify_on_device();
  curCoords.modify_on_device();
}

void
KynemaFMBBeam::map_displacements(double currentTime, bool updateCurCoor)
{
  (void)currentTime;
  for (auto& beam : beam_bodies_)
    map_displacements_beam(beam, updateCurCoor);
}

void
KynemaFMBBeam::map_loads_beam(BeamBody& beam)
{
  if (!beam.kynema_interface || !beam.calc_loads || beam.forcing_surfaces.empty())
    return;

  // First calculate the assembled forces on the beam surfaces
  beam.calc_loads->initialize();
  beam.calc_loads->execute();


  auto& meta = beam.bulk->mesh_meta_data();
  const auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*(beam.bulk));
  const stk::mesh::EntityRank entityRank = stk::topology::NODE_RANK;
    
  // get the parts in the current motion frame
  stk::mesh::Selector sel =
    stk::mesh::selectUnion(beam.forcing_surfaces) &
    (meta.locally_owned_part() | meta.globally_shared_part());
  stk::mesh::NgpField<double> modelCoords =
    stk::mesh::get_updated_ngp_field<double>(
      *meta.get_field<double>(entityRank, "coordinates"));    
  stk::mesh::NgpField<double> tforce =
    stk::mesh::get_updated_ngp_field<double>(
      *meta.get_field<double>(entityRank, "tforce"));
  stk::mesh::NgpField<double> beam_xi =
    stk::mesh::get_updated_ngp_field<double>(
      *meta.get_field<double>(entityRank, "beam_xi"));      

  tforce.sync_to_device();
  beam_xi.sync_to_device();

  constexpr int kMaxInterpolationNodes = 32;
  auto n_nodes = beam.n_nodes;
  auto bary_weights = beam.beam_data.bary_weights;
  auto node_xi = beam.beam_data.node_xi;
  auto positions = beam.beam_data.pos;

  for (std::size_t i = 0; i < n_nodes; ++i) {
    for (std::size_t j = 0; j < 6; ++j) {
      beam.beam_data_host.loads(i, j) = 0.0;
    }
  }
  Kokkos::deep_copy(beam.beam_data.loads, beam.beam_data_host.loads);
  auto loads = beam.beam_data.loads;

  kynema_ugf_ngp::run_entity_algorithm(
    "KynemaFMBBeam_map_loads", ngpMesh, entityRank, sel,
    KOKKOS_LAMBDA(
      const kynema_ugf_ngp::NGPMeshTraits<stk::mesh::NgpMesh>::MeshIndex& mi) {

      const double query_point[3] = {modelCoords(mi, 0), modelCoords(mi, 1), modelCoords(mi, 2)};
      double scratch_weights[kMaxInterpolationNodes];
      double closest_position[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      double xi = beam_xi.get(mi, 0);
      const double force[3] = {tforce(mi, 0), tforce(mi, 1), tforce(mi, 2)};
      double moment[3] = {0.0, 0.0, 0.0};

      // Get the closest position on the beam to the query point and compute the relative position vector
      // Relative position vector is expected to not change with beam deformation
      InterpolateFieldAtPoint(
          xi,
          node_xi.data(),
          positions.data(),
          n_nodes,
          7,
          bary_weights.data(),
          scratch_weights,
          closest_position
      );

      double rel_pos[3] = {
        query_point[0] - closest_position[0], 
        query_point[1] - closest_position[1], 
        query_point[2] - closest_position[2]};
      
      // Transferring force from beam surface to beam node will add moment
      CrossProduct3(rel_pos, force, moment);

      // Compute all shape functions for the current xi value to distribute the force and moment to all beam nodes
      LagrangePolynomialInterpWeights(
          xi,
          node_xi.data(),
          bary_weights.data(),
          n_nodes,
          scratch_weights
      );

      // Multiply the force and moment by the shape function weights and add to the loads on each beam node
      for (std::size_t i = 0; i < n_nodes; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
          Kokkos::atomic_add(&loads(i,j), scratch_weights[i] * force[j]);
          Kokkos::atomic_add(&loads(i,j+3), scratch_weights[i] * moment[j]);
        }
      }

    });

  Kokkos::fence();
  Kokkos::deep_copy(beam.beam_data_host.loads, beam.beam_data.loads);

  // Sum the loads across all ranks to get the total load on each beam node
  MPI_Allreduce(
    MPI_IN_PLACE, beam.beam_data_host.loads.data(), 6 * beam.n_nodes, MPI_DOUBLE, MPI_SUM,
    bulk_->parallel());

}

void
KynemaFMBBeam::map_loads(const double currentTime)
{
  for (auto& beam : beam_bodies_) {
    map_loads_beam(beam);

  // // Now apply the loads to the Kynema interface
  // beam.kynema_interface->Blade().ClearLoads();
  // for (std::size_t i = 0; i < n_nodes; ++i) {
  //   for (std::size_t j = 0; j < 6; ++j) {
  //     beam.kynema_interface->Blade().nodes[i].loads[j] =
  //       beam.beam_data_host.loads(i, j);
  //   }
  // }

    auto& tip_node = beam.kynema_interface->Blade().nodes.back();
    tip_node.loads[0] = 0.0005;
    tip_node.loads[1] = 0.0005;
    tip_node.loads[2] = 0.0;
    tip_node.loads[3] = 0.0;
    tip_node.loads[4] = 0.0;
    tip_node.loads[5] = 0.0;
  }
}

void
KynemaFMBBeam::advance_struct_timestep(const double currentTime, const double dT)
{
  for (std::size_t ibeam = 0; ibeam < beam_bodies_.size(); ++ibeam) {
    auto& beam = beam_bodies_[ibeam];
    if (!beam.kynema_interface)
      continue;

    (void)dT;
    const bool converged = beam.kynema_interface->Step();


    if (!converged) {
      KynemaUGFEnv::self().kynema_ugfOutputP0()
        << "Kynema did not converge for beam body " << ibeam << std::endl;
    }

    if (
      !beam.output_file_name.empty() &&
      KynemaUGFEnv::self().parallel_rank() == 0) {
      auto& nodes = beam.kynema_interface->Blade().nodes;
      if (!nodes.empty()) {
        std::ofstream outfile(beam.output_file_name, std::ios::app);
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
