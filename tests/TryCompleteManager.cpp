#include "TryCompleteManager.h"
#include <cmath>
#include <stdexcept>

using namespace sand_reco::kf;

void create_kf_trees(const std::string& filename, StepsTree& steps, TracksTree& tracks) {
  steps.file = TFile::Open(filename.c_str(), "RECREATE");
  if (!steps.file || steps.file->IsZombie()) throw std::runtime_error("Cannot create ROOT file: " + filename);

  // STEPS
  steps.tree = new TTree("steps", "KF per-step plots");
  auto& s = steps.info;
  steps.tree->Branch("run",     &s.run);         
  steps.tree->Branch("event",   &s.event);
  steps.tree->Branch("track_id",&s.track_id);
  steps.tree->Branch("step_idx",&s.step_idx);
  steps.tree->Branch("z",       &s.z);

  steps.tree->Branch("x_true",    &s.x_true);
  steps.tree->Branch("y_true",    &s.y_true);
  steps.tree->Branch("invR_true", &s.invR_true);
  steps.tree->Branch("tanL_true", &s.tanL_true);
  steps.tree->Branch("phi_true",  &s.phi_true);
  
  steps.tree->Branch("x_smooth",    &s.x_smooth);
  steps.tree->Branch("y_smooth",    &s.y_smooth);
  steps.tree->Branch("invR_smooth", &s.invR_smooth);
  steps.tree->Branch("tanL_smooth", &s.tanL_smooth);
  steps.tree->Branch("phi_smooth",  &s.phi_smooth);
  
  steps.tree->Branch("sigma_x_smooth",    &s.sigma_x_smooth);
  steps.tree->Branch("sigma_y_smooth",    &s.sigma_y_smooth);
  steps.tree->Branch("sigma_invR_smooth", &s.sigma_invR_smooth);
  steps.tree->Branch("sigma_tanL_smooth", &s.sigma_tanL_smooth);
  steps.tree->Branch("sigma_phi_smooth",  &s.sigma_phi_smooth);

  steps.tree->Branch("x_meas",  &s.x_meas);
  steps.tree->Branch("y_meas",  &s.y_meas);
  steps.tree->Branch("x_pred",  &s.x_pred);
  steps.tree->Branch("y_pred",  &s.y_pred);
  steps.tree->Branch("x_filt",  &s.x_filt);
  steps.tree->Branch("y_filt",  &s.y_filt);

  steps.tree->Branch("p_true",  &s.p_true);
  steps.tree->Branch("p_smooth",&s.p_smooth);

  steps.tree->Branch("innov_pos", &s.innov_pos);
  steps.tree->Branch("innov_ang", &s.innov_ang);
  steps.tree->Branch("chi2",      &s.chi2);

  // TRACKS
  tracks.tree = new TTree("tracks", "KF per-track summary");
  auto& t = tracks.info;
  tracks.tree->Branch("run",          &t.run); 
  tracks.tree->Branch("event",        &t.event);
  tracks.tree->Branch("track_id",     &t.track_id);
  tracks.tree->Branch("n_steps",      &t.n_steps);
  tracks.tree->Branch("p_true_first",  &t.p_true_first);
  tracks.tree->Branch("p_true_last",  &t.p_true_last);
  tracks.tree->Branch("p_smooth_last",&t.p_smooth_last);

}

void close_kf_trees(StepsTree& steps, TracksTree& tracks) {
  if (!steps.file) return;
  steps.file->cd();
  if (steps.tree)  steps.tree->Write();
  if (tracks.tree) tracks.tree->Write();
  steps.file->Write();
  steps.file->Close();
  steps.tree = nullptr; tracks.tree = nullptr; steps.file = nullptr;
}

void tryCompleteManager(
  sand_reco::kf::utils::TrackletMap& z_to_tracklets,
 SParticleInfo& particle,
  TMultiGraph* /*mg*/, TMultiGraph* /*mgx*/,
  StepsTree* stepsTree,
  TracksTree* tracksTree,
  int run, int event)
{
  Manager manager;
  manager.initFromMC(&z_to_tracklets, particle);
  manager.run();

  const auto& track = manager.getTrack();
  const auto& vsteps = track.getSteps();
  const int nSteps = static_cast<int>(vsteps.size());
  if (nSteps == 0) return;

  // ---- TRACK  ----
  const auto initial_state = sand_reco::kf::utils::getStateVector(
      particle.initial_mom * 1E-3, particle.initial_pos * 1E-3, particle.charge);
  const double p_true_first = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(
      initial_state.radius(), initial_state.tanLambda());

  const auto& last = vsteps.back();
  const auto true_state_last = sand_reco::kf::utils::getStateVector(
      last.getTrueMomentum() * 1E-3, last.getTruePosition() * 1E-3, particle.charge);
  const double p_true_last = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(
      true_state_last.radius(), true_state_last.tanLambda());

  const auto smoothed_last_sv = last.getStage(TrackStep::TrackStateStage::kSmoothing).getStateVector();
  const double p_smooth_last = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(
      smoothed_last_sv.radius(), smoothed_last_sv.tanLambda());

  if (tracksTree && tracksTree->tree) {
    auto& tr = tracksTree->info;
    tr.run = run;
    tr.event = event;
    tr.track_id = particle.id;
    tr.n_steps = nSteps;
    tr.p_true_first = p_true_first;
    tr.p_true_last = p_true_last;
    tr.p_smooth_last = p_smooth_last;
    tracksTree->tree->Fill();
  }

  // ---- STEPS ----
  if (!(stepsTree && stepsTree->tree)) return;
  auto& infos = stepsTree->info;

  for (int i=0;i<nSteps;++i) {
    const auto& st = vsteps[i];
    infos.run = run;
    infos.event = event;
    infos.track_id = particle.id;
    infos.step_idx = i;
    //infos.orientation = Orienttin(st.getOrientation());
    infos.z = st.getZ();

    // STATES
    const auto true_state = sand_reco::kf::utils::getStateVector(
        st.getTrueMomentum()*1E-3, st.getTruePosition()*1E-3, particle.charge); //truth from tracklet 
    const auto& prediction = st.getStage(TrackStep::TrackStateStage::kPrediction).getStateVector();
    const auto& filtering = st.getStage(TrackStep::TrackStateStage::kFiltering).getStateVector();
    const auto& smoothing = st.getStage(TrackStep::TrackStateStage::kSmoothing).getStateVector();
    
    //---truth (from trk)--- 
    infos.x_true = st.getTruePosition().X();
    infos.y_true = st.getTruePosition().Y();
    infos.invR_true = true_state.signedInverseRadius();
    infos.tanL_true = true_state.tanLambda();
    infos.phi_true = true_state.phi();
    //---smooth---
    infos.x_smooth = smoothing.x()*1000.0 ;
    infos.y_smooth = smoothing.y()*1000.0 ;
    infos.invR_smooth = smoothing.signedInverseRadius();
    infos.tanL_smooth = smoothing.tanLambda();
    infos.phi_smooth = smoothing.phi();

    const auto& smooth_cov_matrix = st.getStage(sand_reco::kf::TrackStep::TrackStateStage::kSmoothing).getStateCovMatrix();
    infos.sigma_x_smooth     = (smooth_cov_matrix.GetNrows() > 0 && smooth_cov_matrix.GetNcols() > 0) ? std::sqrt(smooth_cov_matrix(0,0)) * 1000 : NAN;
    infos.sigma_y_smooth     = (smooth_cov_matrix.GetNrows() > 1 && smooth_cov_matrix.GetNcols() > 1) ? std::sqrt(smooth_cov_matrix(1,1)) * 1000 : NAN;
    infos.sigma_invR_smooth  = (smooth_cov_matrix.GetNrows() > 2 && smooth_cov_matrix.GetNcols() > 2) ? std::sqrt(smooth_cov_matrix(2,2)) : NAN;
    infos.sigma_tanL_smooth   = (smooth_cov_matrix.GetNrows()> 3 && smooth_cov_matrix.GetNcols() > 3) ? std::sqrt(smooth_cov_matrix(3,3)) : NAN;
    infos.sigma_phi_smooth   = (smooth_cov_matrix.GetNrows() > 4 && smooth_cov_matrix.GetNcols() > 4) ? std::sqrt(smooth_cov_matrix(4,4)) : NAN;

    //---measurement--- (truth from trk + smearing)    
    infos.x_meas = st.getX();
    infos.y_meas = st.getY();

    infos.x_pred = prediction.x()*1000.0;
    infos.y_pred = prediction.y()*1000.0;
    infos.x_filt  = filtering.x()*1000.0;
    infos.y_filt  = filtering.y()*1000.0;

    // momentum
    infos.p_true   = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(true_state.radius(), true_state.tanLambda());
    infos.p_smooth = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(smoothing.radius(), smoothing.tanLambda());

    // Innovation e chi2
    const auto& innovation = st.getInnovation();
    infos.innov_pos = (innovation.size()>0) ? innovation[0] : NAN;
    infos.innov_ang = (innovation.size()>1) ? innovation[1] : NAN;
    infos.chi2           = st.getChi2();

    stepsTree->tree->Fill();
  }
}
