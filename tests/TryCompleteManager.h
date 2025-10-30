#pragma once
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TMultiGraph.h>

#include "SANDKalmanFilter.h" 
#include "SANDProcessTracklets.h" 
#include "SANDTrackerUtils.h"  
#include "utils.h" 

//--- TTree STEP---
// Unità DA CONTROLLARE: z,x,y in mm; p in MeV; angoli in rad; var_x/var_y in m^2
struct StepInfos {
  Int_t   run{}, event{}, track_id{}, step_idx{}; /*orientation{}*/; // 0=Vertical (meas X), 1=Horizontal (meas Y)
  Double_t z{};                 // z [mm]
  Double_t x_true{}, y_true{}; // truth pos (mm)
  Double_t x_pred{}, y_pred{}; // KF predicted (mm)
  Double_t x_filt{}, y_filt{}; // KF filtered  (mm)
  Double_t x_smooth{}, y_smooth{};    // KF smoothed  (mm)
  Double_t x_meas{}, y_meas{}; // misura
  Double_t p_true{},p_smooth{}; // MeV
  Double_t innov_pos{}, innov_ang{}, chi2{};
  Double_t var_x{}, var_y{};  
};

// ------------ TTree TRACKS ------------
// Prende solo upstream o downstream point
struct TrackInfos {
  Int_t   run{}, event{}, track_id{}, n_steps{};
  Double_t p_init_true{};      // MeV (da particle.initial_mom)
  Double_t p_true_last{};      // MeV (truth all’ultimo step)
  Double_t p_smooth_last{};    // MeV (smoothed all’ultimo step)
};

struct StepsTree { TFile* file{nullptr}; TTree* tree{nullptr}; StepInfos info{}; };
struct TracksTree{ TTree* tree{nullptr}; TrackInfos info{}; };

void create_kf_trees(const std::string& filename, StepsTree& steps, TracksTree& tracks);
void close_kf_trees(StepsTree& steps, TracksTree& tracks);

void tryCompleteManager(
      sand_reco::kf::utils::TrackletMap& z_to_tracklets,
      SParticleInfo& particle,
      TMultiGraph* mg, TMultiGraph* mgx,
      StepsTree* stepsTree,
      TracksTree* tracksTree,
      int run, int event);