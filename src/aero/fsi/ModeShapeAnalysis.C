// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#include "aero/fsi/ModeShapeAnalysis.h"
#include "aero/fsi/FSIturbine.h"
#include "aero/aero_utils/WienerMilenkovic.h"
#include <NaluParsing.h>

#include "netcdf.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>

namespace sierra {

namespace nalu {

inline void
check_nc_error(int code, std::string msg)
{
  if (code != 0)
    throw std::runtime_error("FSIturbine:: NetCDF error: " + msg);
}

ModeShapeAnalysis::ModeShapeAnalysis(const YAML::Node& node)
  : mesh_motion_(true)
{
  load(node);
}

void
ModeShapeAnalysis::end_openfast()
{
}

void
ModeShapeAnalysis::load(const YAML::Node& node)
{

  fsiTurbineData_ = std::make_unique<fsiTurbine>(0, node);
  get_required(node, "turbine_base_pos", turbineBasePos_);

  get_required(node, "nc_file_name", ncFileName_);

  int ncid, ierr;
  ierr = nc_open(ncFileName_.c_str(), NC_NOCLOBBER, &ncid);
  //check_nc_error(ierr, "nc_open");

  int dimid;
  size_t n_bld_nds;

  ierr = nc_inq_dimid(ncid, "n_bld_nds", &dimid);
  ierr = nc_inq_dimlen(ncid, dimid, &n_bld_nds);
  fsiTurbineData_->params_.numBlades = 1;
  fsiTurbineData_->params_.nBRfsiPtsBlade.resize(1);
  fsiTurbineData_->params_.nBRfsiPtsBlade[0] = n_bld_nds;
  fsiTurbineData_->params_.nBRfsiPtsTwr = 0;

  if (node["mode"]) {
    get_required(node["mode"], "freq", modeFreq_);
    get_required(node["mode"], "shape", modeShape_);
    NaluEnv::self().naluOutputP0() << "Mode shape at freq " << modeFreq_  << " is :" << std::endl;
    for (int i = 0; i < n_bld_nds; i++) {
        NaluEnv::self().naluOutputP0() << modeShape_[i][0] << ", " << modeShape_[i][1] << ", " << modeShape_[i][2] << ", "
                                       << modeShape_[i][3] << ", "  << modeShape_[i][4] << ", "  << modeShape_[i][5] << std::endl;
    }
  } else {
      throw std::runtime_error(
          "mode is required in mode_shape_analysis with both freq and shape entries");
  }

  ierr = nc_close(ncid);
  // check_nc_error(ierr, "nc_close");
  std::cerr << "Only dealing with 1 blade for now corresponding to blade 0 " << std::endl;
}

void
ModeShapeAnalysis::setup(double dtNalu, std::shared_ptr<stk::mesh::BulkData> bulk)
{
  bulk_ = bulk;
  dt_ = dtNalu;

  fsiTurbineData_->setup(bulk_);
}

void
ModeShapeAnalysis::initialize(int restartFreqNalu, double curTime)
{

  // TODO: Check here on the processor containing the turbine that the number of
  // blades on the turbine is the same as the number of blade parts specified in
  // the input file.

  // TODO: In the documentation, mention that the CFD mesh must always be
  // created for the turbine in the reference position defined in OpenFAST, i.e.
  // with blade1 pointing up and the other blades following it in order as the
  // turbine rotates clockwise facing downwind. If the mesh is not created this
  // way, the mapping won't work. Any non-zero initial azimuth and/or initial
  // yaw must be only specified in the OpenFAST input file and the mesh will
  // automatically be deformed after calling solution0. Requiring the initial
  // CFD mesh to be in the reference configuration may not always work if the
  // mesh domain and initial yaw setting does not align with the reference
  // configuration. May be this is isn't an issue because unlike AeroDyn,
  // ExtLoads does not create the initial mesh independent of the
  // ElastoDyn/BeamDyn. May be ExtLoads already has the correct yaw and azimuth
  // setting from OpenFAST after the init call. In this case, the CFD mesh must
  // start in the correct azimuth and yaw configuration. In which case, the
  // initial yaw and azimuth must be obtained from OpenFAST and the mesh around
  // the turbine must be deformed through rigid body motion first before
  // starting any mapping.

  // TODO: Get parameters here
  // FAST.get_turbineParams(i, fsiTurbineData_->params_);


  fsiTurbineData_->initialize();

  compute_mapping();

  map_displacements(curTime, false);

  if (curTime < 1e-10) {

    auto& meta = bulk_->mesh_meta_data();

    const VectorFieldType* meshDisp =
      meta.get_field<double>(stk::topology::NODE_RANK, "mesh_displacement");
    const VectorFieldType* meshVel =
      meta.get_field<double>(stk::topology::NODE_RANK, "mesh_velocity");

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
        for (size_t i = 0; i < 3; i++) {
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
ModeShapeAnalysis::compute_mapping()
{


        size_t nBldPts = fsiTurbineData_->params_.nBRfsiPtsBlade[0];

        int ncid, ierr;
        ierr = nc_open(ncFileName_.c_str(), NC_NOCLOBBER, &ncid);
        //check_nc_error(ierr, "nc_open");

        std::vector<double> tmpArray;
        int nc_var_id;

        {
        nc_inq_varid(ncid, "bld_ref_pos", &nc_var_id);
        tmpArray.resize(nBldPts);
        std::vector<size_t> count_dim{1, 1, nBldPts};
        for (size_t iDim = 0; iDim < 3; iDim++) {
            std::vector<size_t> start_dim{0, iDim, 0};
            ierr = nc_get_vara_double(
              ncid, nc_var_id, start_dim.data(), count_dim.data(),
              tmpArray.data());
            for (size_t i = 0; i < nBldPts; i++)
                fsiTurbineData_->brFSIdata_.bld_ref_pos[i * 6 + iDim] = tmpArray[i];
        }
        nc_inq_varid(ncid, "bld_ref_orient", &nc_var_id);
        for (size_t iDim = 0; iDim < 3; iDim++) {
            std::vector<size_t> start_dim{0, iDim, 0};
            ierr = nc_get_vara_double(
                                      ncid, nc_var_id, start_dim.data(), count_dim.data(),
                                      tmpArray.data());
            for (size_t i = 0; i < nBldPts; i++)
                fsiTurbineData_->brFSIdata_.bld_ref_pos[i * 6 + 3 + iDim] = tmpArray[i];
        }

        for (size_t i = 0; i < nBldPts; i++) {
            NaluEnv::self().naluOutputP0() << fsiTurbineData_->brFSIdata_.bld_ref_pos[i * 6 + 0] << ", "
                                           << fsiTurbineData_->brFSIdata_.bld_ref_pos[i * 6 + 1] << ", "
                                           << fsiTurbineData_->brFSIdata_.bld_ref_pos[i * 6 + 2] << ", "
                                           << fsiTurbineData_->brFSIdata_.bld_ref_pos[i * 6 + 3] << ", "
                                           << fsiTurbineData_->brFSIdata_.bld_ref_pos[i * 6 + 4] << ", "
                                           << fsiTurbineData_->brFSIdata_.bld_ref_pos[i * 6 + 5] << std::endl;

        }

        }

        {
        nc_inq_varid(ncid, "bld_root_ref_pos", &nc_var_id);
        tmpArray.resize(3);
        std::vector<size_t> count_dim{1,3};
        std::vector<size_t> start_dim{0,0};
        ierr = nc_get_vara_double(
                                  ncid, nc_var_id, start_dim.data(), count_dim.data(),
                                  tmpArray.data());
        for (size_t i = 0; i < 3; i++)
            fsiTurbineData_->brFSIdata_.bld_root_ref_pos[i] = tmpArray[i];

        nc_inq_varid(ncid, "hub_ref_orient", &nc_var_id);
        ierr = nc_get_vara_double(
                                  ncid, nc_var_id, start_dim.data(), count_dim.data(),
                                  tmpArray.data());
        for (size_t i = 0; i < 3; i++)
            fsiTurbineData_->brFSIdata_.bld_root_ref_pos[i+3] = tmpArray[i];

        nc_inq_varid(ncid, "hub_ref_pos", &nc_var_id);
        ierr = nc_get_var_double(
                                 ncid, nc_var_id,
                                 tmpArray.data());
        for (size_t i = 0; i < 3; i++)
            fsiTurbineData_->brFSIdata_.hub_ref_pos[i] = tmpArray[i];

        nc_inq_varid(ncid, "hub_ref_orient", &nc_var_id);
        ierr = nc_get_var_double(
                                 ncid, nc_var_id,
                                 tmpArray.data());
        for (size_t i = 0; i < 3; i++)
            fsiTurbineData_->brFSIdata_.hub_ref_pos[i+3] = tmpArray[i];
        }

        ierr = nc_close(ncid);

        /* FAST.getTowerRefPositions( */
        /*   fsiTurbineData_->brFSIdata_.twr_ref_pos.data(), i); */
        /* FAST.getBladeRefPositions( */
        /*   fsiTurbineData_->brFSIdata_.bld_ref_pos.data(), i); */
        /* FAST.getBladeRootRefPositions( */
        /*   fsiTurbineData_->brFSIdata_.bld_root_ref_pos.data(), i); */
        /* FAST.getHubRefPosition( */
        /*   fsiTurbineData_->brFSIdata_.hub_ref_pos.data(), i); */
        /* FAST.getNacelleRefPosition( */
        /*   fsiTurbineData_->brFSIdata_.nac_ref_pos.data(), i); */
        /* FAST.getBladeRloc(fsiTurbineData_->brFSIdata_.bld_rloc.data(), i); */
        /* FAST.getBladeChord(fsiTurbineData_->brFSIdata_.bld_chord.data(), i); */


      // int iError = MPI_Bcast(
      //   fsiTurbineData_->brFSIdata_.twr_ref_pos.data(),
      //   (fsiTurbineData_->params_.nBRfsiPtsTwr) * 6, MPI_DOUBLE, 0,
      //   bulk_->parallel());
      // int nTotBldNodes = fsiTurbineData_->params_.nTotBRfsiPtsBlade;
      // int nBlades = fsiTurbineData_->params_.numBlades;
      // iError = MPI_Bcast(
      //   fsiTurbineData_->brFSIdata_.bld_ref_pos.data(), nTotBldNodes * 6,
      //   MPI_DOUBLE, 0, bulk_->parallel());
      // iError = MPI_Bcast(
      //   fsiTurbineData_->brFSIdata_.bld_root_ref_pos.data(), nBlades * 6,
      //   MPI_DOUBLE, 0, bulk_->parallel());
      // iError = MPI_Bcast(
      //   fsiTurbineData_->brFSIdata_.hub_ref_pos.data(), 6, MPI_DOUBLE,
      //   0, bulk_->parallel());
      // iError = MPI_Bcast(
      //   fsiTurbineData_->brFSIdata_.nac_ref_pos.data(), 6, MPI_DOUBLE,
      //   0, bulk_->parallel());
      // iError = MPI_Bcast(
      //   fsiTurbineData_->brFSIdata_.bld_rloc.data(), nTotBldNodes,
      //   MPI_DOUBLE, 0, bulk_->parallel());
      // No need to bcast chord
      fsiTurbineData_->computeMapping();
      fsiTurbineData_->computeLoadMapping();

}

void
ModeShapeAnalysis::predict_struct_states()
{
}

void
ModeShapeAnalysis::predict_struct_timestep(const double curTime)
{
  send_loads(curTime);
}

void
ModeShapeAnalysis::advance_struct_timestep(const double curTime)
{

  tStep_ += 1;

  // int nTurbinesGlob = FAST.get_nTurbinesGlob();
  // for (int i=0; i < nTurbinesGlob; i++) {
  //     if(fsiTurbineData_[i] != nullptr)
  //         fsiTurbineData_[i]->write_nc_def_loads(tStep_, curTime);
  // }
}

void
ModeShapeAnalysis::send_loads(const double curTime)
{

  fsiTurbineData_->mapLoads();

  int nTotBldNodes = fsiTurbineData_->params_.nTotBRfsiPtsBlade;
  if (bulk_->parallel_rank() == 0) {
    int iError = MPI_Reduce(
      MPI_IN_PLACE, fsiTurbineData_->brFSIdata_.twr_ld.data(),
      (fsiTurbineData_->params_.nBRfsiPtsTwr) * 6, MPI_DOUBLE, MPI_SUM,
      0, bulk_->parallel());
    iError = MPI_Reduce(
      MPI_IN_PLACE, fsiTurbineData_->brFSIdata_.bld_ld.data(),
      nTotBldNodes * 6, MPI_DOUBLE, MPI_SUM, 0, bulk_->parallel());
  } else {
    int iError = MPI_Reduce(
      fsiTurbineData_->brFSIdata_.twr_ld.data(), NULL,
      (fsiTurbineData_->params_.nBRfsiPtsTwr) * 6, MPI_DOUBLE, MPI_SUM,
      0, bulk_->parallel());
    iError = MPI_Reduce(
      fsiTurbineData_->brFSIdata_.bld_ld.data(), NULL,
      (nTotBldNodes) * 6, MPI_DOUBLE, MPI_SUM, 0, bulk_->parallel());
  }

}

void
ModeShapeAnalysis::get_displacements(double current_time)
{
    // TODO:: Set displacements here from mode shapes

    double sinomegat = stk::math::sin(2.0 * 3.14159265358979323846 * modeFreq_ * current_time);
    size_t n_bld_nds = fsiTurbineData_->params_.nBRfsiPtsBlade[0];
    vs::Vector wm_ref;
    vs::Vector wm_def;
    vs::Vector wm_tmp;
    vs::Vector wm_bld_root;
    vs::Vector wm_final;
    for (int j = 0; j < 3; j++)
        wm_bld_root[j] = -fsiTurbineData_->brFSIdata_.bld_root_ref_pos[3+j];
    for (size_t i = 0; i < n_bld_nds; i++) {
        for (int j = 0; j < 3; j++)
            fsiTurbineData_->brFSIdata_.bld_def[i*6+j] = modeShape_[i][j] * sinomegat;
        for (int j = 0; j < 3; j++) {
            wm_ref[j] = fsiTurbineData_->brFSIdata_.bld_ref_pos[i*6+3+j] ;
            wm_def[j] = modeShape_[i][j+3] * sinomegat;
        }
        wm_final = wmp::compose(wm_def, wm_ref);
        //wm_final = wmp::compose(wm_bld_root, wm_tmp);

        for (int j = 0; j < 3; j++)
            fsiTurbineData_->brFSIdata_.bld_def[i*6+3+j] = wm_final[j];
    }
}

void
ModeShapeAnalysis::compute_div_mesh_velocity()
{
  timer_start(naluTimer_);
  fsiTurbineData_->compute_div_mesh_velocity();
  timer_stop(naluTimer_);
}

void
ModeShapeAnalysis::set_rotational_displacement(
  std::array<double, 3> axis, double omega, double curTime)
{

  fsiTurbineData_->setRotationDisplacement(axis, omega, curTime);
}

void
ModeShapeAnalysis::map_displacements(double current_time, bool updateCurCoor)
{

  timer_start(naluTimer_);
  get_displacements(current_time);

  stk::mesh::Selector sel;
  fsiTurbineData_->mapDisplacements(current_time);
  sel &= stk::mesh::selectUnion(fsiTurbineData_->getPartVec());

  if (updateCurCoor) {
    auto& meta = bulk_->mesh_meta_data();
    const VectorFieldType* modelCoords =
      meta.get_field<double>(stk::topology::NODE_RANK, "coordinates");
    VectorFieldType* curCoords =
      meta.get_field<double>(stk::topology::NODE_RANK, "current_coordinates");
    VectorFieldType* displacement =
      meta.get_field<double>(stk::topology::NODE_RANK, "mesh_displacement");

    modelCoords->sync_to_host();
    curCoords->sync_to_host();
    displacement->sync_to_host();

    const auto& bkts = bulk_->get_buckets(stk::topology::NODE_RANK, sel);
    for (const auto* b : bkts) {
      for (const auto node : *b) {
        for (size_t in = 0; in < b->size(); in++) {

          double* cc = stk::mesh::field_data(*curCoords, node);
          double* mc = stk::mesh::field_data(*modelCoords, node);
          double* cd = stk::mesh::field_data(*displacement, node);

          for (int j = 0; j < 3; ++j) {
            cc[j] = mc[j] + cd[j];
          }
        }
      }
    }

    curCoords->modify_on_host();
    curCoords->sync_to_device();
  }
  timer_stop(naluTimer_);
}

void
ModeShapeAnalysis::map_loads(const int tStep, const double curTime)
{
  timer_start(naluTimer_);
  fsiTurbineData_->mapLoads();
  int nTotBldNodes = fsiTurbineData_->params_.nTotBRfsiPtsBlade;
  if (bulk_->parallel_rank() == 0) {
    int iError = MPI_Reduce(
      MPI_IN_PLACE, fsiTurbineData_->brFSIdata_.twr_ld.data(),
      (fsiTurbineData_->params_.nBRfsiPtsTwr) * 6, MPI_DOUBLE, MPI_SUM,
      0, bulk_->parallel());
    iError = MPI_Reduce(
      MPI_IN_PLACE, fsiTurbineData_->brFSIdata_.bld_ld.data(),
      nTotBldNodes * 6, MPI_DOUBLE, MPI_SUM, 0, bulk_->parallel());
  } else {
    int iError = MPI_Reduce(
      fsiTurbineData_->brFSIdata_.twr_ld.data(), NULL,
      (fsiTurbineData_->params_.nBRfsiPtsTwr) * 6, MPI_DOUBLE, MPI_SUM,
      0, bulk_->parallel());
    iError = MPI_Reduce(
      fsiTurbineData_->brFSIdata_.bld_ld.data(), NULL,
      (nTotBldNodes) * 6, MPI_DOUBLE, MPI_SUM, 0, bulk_->parallel());
  }
  fsiTurbineData_->write_nc_def_loads(tStep, curTime);

  timer_stop(naluTimer_);
}

void
ModeShapeAnalysis::timer_start(std::pair<double, double>& timer)
{
  timer.first = NaluEnv::self().nalu_time();
}

void
ModeShapeAnalysis::timer_stop(std::pair<double, double>& timer)
{
  timer.first = NaluEnv::self().nalu_time() - timer.first;
  timer.second += timer.first;
}

} // namespace nalu

} // namespace sierra
