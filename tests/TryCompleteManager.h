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
// Unità DA CONTROLLARE  E RISCRIVERE: z,x,y in mm; p in MeV; angoli in rad; var_x/var_y in m^2
struct StepInfos {
  Int_t   run{}, event{}, track_id{}, step_idx{}; /*orientation{}*/; // 0=Vertical (meas X), 1=Horizontal (meas Y)
  Double_t z{}; //mm
  //---true---
  Double_t x_true{};
  Double_t y_true{};
  Double_t invR_true{};
  Double_t tanL_true{};
  Double_t phi_true{};
  //---smooth---
  Double_t x_smooth{};
  Double_t y_smooth{};
  Double_t invR_smooth{};
  Double_t tanL_smooth{};
  Double_t phi_smooth{};
  
  Double_t sigma_x_smooth{};
  Double_t sigma_y_smooth{};
  Double_t sigma_invR_smooth{};
  Double_t sigma_tanL_smooth{};
  Double_t sigma_phi_smooth{};

  Double_t x_meas{}, y_meas{};
  Double_t x_pred{}, y_pred{}; 
  Double_t x_filt{}, y_filt{};

  Double_t p_true{}; //MeV
  Double_t p_smooth{}; //MeV
  Double_t innov_pos{}, innov_ang{}, chi2{};

};

// ------------ TTree TRACKS ------------
struct TrackInfos {
  Int_t   run{}, event{}, track_id{}, n_steps{};

  Double_t x_true_seed{};
  Double_t y_true_seed{};
  Double_t invR_true_seed{};
  Double_t tanL_true_seed{};
  Double_t phi_true_seed{};
  
  Double_t x_smooth_seed{};
  Double_t y_smooth_seed{};
  Double_t invR_smooth_seed{};
  Double_t tanL_smooth_seed{};
  Double_t phi_smooth_seed{};
  
  Double_t p_true_first{};      // MeV (da particle.initial_mom)
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