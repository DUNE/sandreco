#include <TVector3.h>
#include <TTree.h>
#include <TFile.h>
#include <TMarker.h>
#include <TArrow.h>
#include <TDatabasePDG.h>
#include <TGraph.h>
#include <TMultiGraph.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <random>
#include <limits>

#include "SANDGeoManager.h"
#include "SANDTrackletFinder.h"
#include "SANDProcessTracklets.h"
#include "SANDTrackerClusterCollection.h"
#include "SANDTrackerDigitCollection.h"
#include "SANDKalmanFilter.h"
#include "SANDTrackerUtils.h"
#include "utils.h"

struct PlotContainer {
  TH1D meas_x_res;
  TH1D meas_y_res;
  TH1D smooth_x_res;
  TH1D smooth_y_res;
  TH1D smooth_p_res;
  TH1D smooth_tan_res;
  TH1D smooth_phi_res;
  TH1D x_pull_smooth;
  TH1D y_pull_smooth;
  TH1D r_inv_pull_smooth;
  TH1D tan_pull_smooth;
  TH1D phi_pull_smooth;
  TH1D p_res_step;
  TH2D p_init_vs_reco_smooth;
  TH2D p_first_vs_reco_smooth;
  TH2D p_diff_vs_points;
  TH1D h_gpos_distribution;
  TH1D h_gang_distribution;
  TH1D deltaP_over_P;
  TH1D chi2;

  PlotContainer()
    : meas_x_res("meas_x_res", "meas_x_res; meas_x_res; (x_{meas} - x_{true}) [mm]; Entries", 200, -2.0, 2.0),
      meas_y_res("meas_y_res", "meas_y_res; meas_y_res; (y_{meas} - y_{true}) [mm]; Entries", 200, -2.0, 2.0),
      smooth_x_res("smooth_x_res", "smooth_x_res; (x_{smoothed} - x_{true}) [mm]; Entries", 200, -2.0, 2.0),
      smooth_y_res("smooth_y_res", "smooth_y_res; (y_{smoothed} - y_{true}) [mm]; Entries", 200, -2.0, 2.0),
      smooth_p_res("smooth_p_res", "smooth_p_res; (p_{smoothed} - p_{true}) [GeV]; Entries", 200, -200, 200),
      smooth_tan_res("smooth_tan_res", "smooth_tan_res; (tan_{smoothed} - tan_{true}); Entries", 400, -0.2, 0.2),
      smooth_phi_res("smooth_phi_res", "smooth_phi_res; (phi_{smoothed} - phi_{true}) [rad]; Entries", 400, -0.2, 0.2),
      x_pull_smooth("x_pull_smooth", "(smoothed_x - true_x)/#sigma_{x}; (smoothed_x - true_x)/#sigma_{x}; Entries", 200, -5.0, 5.0),
      y_pull_smooth("y_pull_smooth", "(smoothed_y - true_y)/#sigma_{y}; (smoothed_y - true_y)/#sigma_{y}; Entries", 200, -5.0, 5.0),
      r_inv_pull_smooth("r_inv_pull_smooth", "(smoothed_r_inv - true_r_inv)/#sigma_{r_inv}; (smoothed_r_inv - true_r_inv)/#sigma_{r_inv}; Entries", 200, -5.0, 5.0),
      tan_pull_smooth("tan_pull_smooth", "(smoothed_tan - true_tan)/#sigma_{tan}; (smoothed_tan - true_tan)/#sigma_{tan}; Entries", 200, -5.0, 5.0),
      phi_pull_smooth("phi_pull_smooth", "(smoothed_phi - true_phi)/#sigma_{phi}; (smoothed_phi - true_phi)/#sigma_{phi}; Entries", 200, -5.0, 5.0),
      p_res_step("p_res_step", "p_{smooth} - p_{true} per step; MeV; Entries", 200, -500.0, 500.0),
      p_init_vs_reco_smooth("p_init_vs_reco_smooth", "Initial true p vs last smoothed reco p; p_{true}^{init} [MeV]; p_{reco}^{smooth,last} [MeV]", 200, 0, 5000, 200, 0, 5000),
      p_first_vs_reco_smooth("p_first_vs_reco_smooth", "True p  at first hit vs last smoothed reco p; p_{true}^{init} [MeV]; p_{reco}^{smooth,last} [MeV]", 200, 0, 5000, 200, 0, 5000),
      p_diff_vs_points("p_diff_vs_points", "DeltaP vs nPoints p; nPoints; (p_{true} - p{smoothed}) [MeV];", 200, 0, 200, 200, -200, 200),
      h_gpos_distribution("h_gpos_distribution", "Innovation on position", 100, -3, 3),
      h_gang_distribution("h_gang_distribution", "Innovation on direction", 100, -3, 3),
      deltaP_over_P("deltaP_over_P", "DeltaP vs True P", 100, -1, 1),
      chi2("chi2", "chi2", 1000, 0, 100000)
  {}
};


static void printChecks(
    int track_ID,
    const sand_reco::kf::utils::TrackletMap& z_to_tracklets,
    const sand_reco::kf::Track& track)
{

  //# Planes with at least one measurement
  //# Total number of measurements per track 
  size_t planes_with_meas = 0;
  size_t total_meas = 0;
  for (const auto& kv : z_to_tracklets) {
    const auto& v = kv.second;
    if (!v.empty()) {
      planes_with_meas++;
      total_meas += v.size();
    }
  }

  //# Total steps
  const size_t n_steps = track.getSteps().size();
  //# True tracklet
  size_t n_true_tracklet= 0;
  //# Tracklet used as measurements
  size_t n_measured_tracklet = 0;
  
  for (const auto& s : track.getSteps()) {
    const TMatrixD& m = s.getMeasurement();
    if (m.GetNrows() >= 2 && m.GetNcols() >= 1) {
      ++n_measured_tracklet;
    }

    const TVector3& ptrue = s.getTrueMomentum();
    const TVector3& xtrue = s.getTruePosition();
    if (std::isfinite(ptrue.X()) && std::isfinite(ptrue.Y()) && std::isfinite(ptrue.Z()) &&
        std::isfinite(xtrue.X()) && std::isfinite(xtrue.Y()) && std::isfinite(xtrue.Z())) {
      ++n_true_tracklet;
    }
  }


  std::cout
    << "Traccia " << track_ID << ":\n"
    << "  # piani con almeno una misura   = " << planes_with_meas << "\n"
    << "  # misure totali (event)         = " << total_meas << "\n"
    << "  # step (traccia)                = " << n_steps << "\n"
    << "  # true tracklet (traccia)       = " << n_true_tracklet<< "\n"
    << "  # measurement (traccia)         = " << n_measured_tracklet << "\n"
    << std::flush;
}



  void tryCompleteManager(sand_reco::kf::utils::TrackletMap z_to_tracklets, SParticleInfo particle, 
                          TMultiGraph* mg, TMultiGraph* mgx, PlotContainer& plot_container)
  
  {

  sand_reco::kf::Manager manager;
  manager.initFromMC(&z_to_tracklets, particle);
  manager.run();
  
  const auto& track = manager.getTrack();
  printChecks(particle.id, z_to_tracklets, track);
  
  if (track.getSteps().size() > 0) { // was commented, with > 3
    auto initial_state = sand_reco::kf::utils::getStateVector(particle.initial_mom * 1E-3, particle.initial_pos * 1E-3, particle.charge);
    auto initial_mom = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(initial_state.radius(), initial_state.tanLambda());

    auto last_step = track.getSteps().back(); //crash if empty due to the .back().
    auto smoothed_state = last_step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kSmoothing).getStateVector();
    auto smoothed_mom = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(smoothed_state.radius(), smoothed_state.tanLambda());
    auto true_state = sand_reco::kf::utils::getStateVector(last_step.getTrueMomentum(), last_step.getTruePosition(), particle.charge);
    auto true_mom = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(true_state.radius(), true_state.tanLambda());

    if (std::isfinite(initial_mom) && std::isfinite(smoothed_mom)) {
      plot_container.p_init_vs_reco_smooth.Fill(initial_mom, smoothed_mom);
    }
    if (std::isfinite(true_mom) && std::isfinite(smoothed_mom)) {
      plot_container.p_first_vs_reco_smooth.Fill(true_mom, smoothed_mom);
    }

    plot_container.p_diff_vs_points.Fill( track.getSteps().size(), smoothed_mom - true_mom);

    std::cout << "Initial Momentum " << initial_mom << std::endl;
    std::cout << "Smoothed Momentum " << smoothed_mom << std::endl;
    std::cout << "True Momentum " << true_mom << std::endl;
    std::cout << "track.getSteps() " << track.getSteps().size() << std::endl;
    
    TGraph* yz_predicted = new TGraph(track.getSteps().size());
    TGraph* yz_filtered = new TGraph(track.getSteps().size());
    TGraph* yz_smoothed = new TGraph(track.getSteps().size());
    TGraph* yz_measured = new TGraph(track.getSteps().size());
    TGraph* xz_predicted = new TGraph(track.getSteps().size());
    TGraph* xz_filtered = new TGraph(track.getSteps().size());
    TGraph* xz_smoothed = new TGraph(track.getSteps().size());
    TGraph* xz_measured = new TGraph(track.getSteps().size());

    for (uint i = 0; i < track.getSteps().size(); i++) {
      auto step = track.getSteps()[i];

      // -------------------------------------------------------------------------------------
      // Plots trajectories reconstructed and KF stages
      // -------------------------------------------------------------------------------------
      auto prediction = step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kPrediction).getStateVector();
      auto filtering  = step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kFiltering).getStateVector();
      auto smoothing  = step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kSmoothing).getStateVector();
      
      yz_predicted->SetPoint(i, step.getZ(), prediction.y()*1000 );
      yz_filtered->SetPoint(i, step.getZ() , filtering.y()*1000);
      yz_smoothed->SetPoint(i, step.getZ() , smoothing.y()*1000);
      yz_measured->SetPoint(i, step.getZ() , step.getY());
      xz_predicted->SetPoint(i, step.getZ(), prediction.x()*1000 );
      xz_filtered->SetPoint(i, step.getZ() , filtering.x()*1000);
      xz_smoothed->SetPoint(i, step.getZ() , smoothing.x()*1000);
      xz_measured->SetPoint(i, step.getZ() , step.getX());
      

      // -------------------------------------------------------------------------------------
      // Plots for pull-tests
      // -------------------------------------------------------------------------------------
      // Parameters from true tracklet
      const TVector3& true_pos_from_trk = step.getTruePosition();   // mm
      const TVector3& true_mom_from_trk = step.getTrueMomentum();
      auto true_step_state = sand_reco::kf::utils::getStateVector(true_mom_from_trk, true_pos_from_trk, particle.charge);
      auto true_step_mom = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(
                                                  true_step_state.radius(), true_step_state.tanLambda());
      
      // Parameters from measurement vector (which is the true with a smearing)
      const auto& meas = step.getMeasurement();   // TMatrixD 2x1
      const double meas_x = step.getX();
      const double meas_y = step.getY();

      // Parameters from reconstructed trajectory of KF
      const double smooth_x = smoothing.x();
      const double smooth_y = smoothing.y();
      const double smooth_p = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(
                                                  smoothing.radius(), smoothing.tanLambda());

      
      // plot
      // Position residuals
      plot_container.meas_x_res.Fill(meas_x - true_pos_from_trk.X());  // mm
      plot_container.meas_y_res.Fill(meas_y - true_pos_from_trk.Y());  // mm

      const double res_smooth_x = smooth_x - true_pos_from_trk.X();  // mm
      const double res_smooth_y = smooth_y - true_pos_from_trk.Y();  // mm
      const double res_smooth_p = (smooth_p - true_step_mom);  // GeV
      const double res_smooth_r_inv = (smoothing.signedInverseRadius() - true_step_state.signedInverseRadius());  // 1 / m
      const double res_smooth_tan = smoothing.tanLambda() - true_step_state.tanLambda();  // 
      const double res_smooth_phi = smoothing.phi() - true_step_state.phi();  // rad
      plot_container.smooth_x_res.Fill(res_smooth_x);
      plot_container.smooth_y_res.Fill(res_smooth_y);
      plot_container.smooth_p_res.Fill(res_smooth_p);
      plot_container.smooth_tan_res.Fill(res_smooth_tan);
      plot_container.smooth_phi_res.Fill(res_smooth_phi);
      
      plot_container.deltaP_over_P.Fill((smooth_p - true_step_mom) / true_step_mom );

      // Momentum residuals
      const auto& smooth_cov_matrix = step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kSmoothing).getStateCovMatrix();
      const double sigma_x = (smooth_cov_matrix.GetNrows() > 0 && smooth_cov_matrix.GetNcols() > 0) ? std::sqrt(smooth_cov_matrix(0,0)) * 1000.0 : 0.0;
      const double sigma_y = (smooth_cov_matrix.GetNrows() > 1 && smooth_cov_matrix.GetNcols() > 1) ? std::sqrt(smooth_cov_matrix(1,1)) * 1000.0 : 0.0;

      if (x_pull_smooth && sigma_x > 0.0 && std::isfinite(sigma_x))
        x_pull_smooth->Fill(res_smooth_x / sigma_x);
      const double sigma_phi   = (smooth_cov_matrix.GetNrows() > 4 && smooth_cov_matrix.GetNcols() > 4) ? std::sqrt(smooth_cov_matrix(4,4)) : 0.0;

      if (sigma_x > 0.0 && std::isfinite(sigma_x))
        plot_container.x_pull_smooth.Fill(res_smooth_x / sigma_x);

      if (sigma_y > 0.0 && std::isfinite(sigma_y))
       plot_container.y_pull_smooth.Fill(res_smooth_y / sigma_y);

      if (sigma_r_inv > 0.0 && std::isfinite(sigma_r_inv))
       plot_container.r_inv_pull_smooth.Fill(res_smooth_r_inv / sigma_r_inv);

      if (sigma_tan > 0.0 && std::isfinite(sigma_tan))
       plot_container.tan_pull_smooth.Fill(res_smooth_tan / sigma_tan);

      if (sigma_phi > 0.0 && std::isfinite(sigma_phi))
       plot_container.phi_pull_smooth.Fill(res_smooth_phi / sigma_phi);
      
      // Innovation test
      auto& innovation = step.getInnovation();
      if (innovation.empty()) {
        continue;
      }
        
        plot_container.h_gpos_distribution.Fill(innovation[0]);
        plot_container.h_gang_distribution.Fill(innovation[1]);
        
        plot_container.chi2.Fill(step.getChi2());


          // if (meas.GetNrows() < 2 || meas.GetNcols() < 1) continue;
          // const double meas_pos  = meas(0,0); 
          // const double meas_ang  = meas(1,0); 
          
          // const auto& orientation = step.getOrientation();
          // if (orientation == sand_reco::kf::Orientation::kVertical) {
          //     // Piano x
          //     const double res_x = true_pos.X() - meas_pos;
          //     smeared_meas_x_res->Fill(res_x);
          //   } else {
          //     // Piano y
          //     const double res_y = true_pos.Y() - meas_pos;
          //     smeared_meas_y_res->Fill(res_y);
          //   }
          
        }
        
    yz_predicted->SetLineColor(3);
    yz_predicted->SetMarkerStyle(3);
    mg->Add(yz_predicted);
    yz_filtered->SetLineColor(4);
    yz_filtered->SetMarkerStyle(4);
    mg->Add(yz_filtered);
    yz_smoothed->SetLineColor(6);
    yz_smoothed->SetMarkerStyle(5);
    mg->Add(yz_smoothed);
    yz_measured->SetLineColor(2);
    yz_measured->SetMarkerStyle(2);
    mg->Add(yz_measured);

    xz_predicted->SetLineColor(3);
    xz_predicted->SetMarkerStyle(3);
    mgx->Add(xz_predicted);
    xz_filtered->SetLineColor(4);
    xz_filtered->SetMarkerStyle(4);
    mgx->Add(xz_filtered);
    xz_smoothed->SetLineColor(6);
    xz_smoothed->SetMarkerStyle(5);
    mgx->Add(xz_smoothed);
    xz_measured->SetLineColor(2);
    xz_measured->SetMarkerStyle(2);
    mgx->Add(xz_measured);
  }

  return;
}

void processEventWithKF(SANDGeoManager* sand_geo, TG4Event* mc_event, std::vector<dg_wire>* digits, PlotContainer& plot_container)

{
  // unknown parameters
  int p[9] = {100, -2000, 2000, 100, -4000, -0, 100, 22500, 26000};

  // fill ID vs index map of DigitCollection
  sand_reco::tracker::DigitCollection::fillMap(digits);

  // get the vector of digits
  auto digit_vec =  sand_reco::tracker::DigitCollection::getDigits();

  // check if the vector of digits is empty
  if (digit_vec.empty()) {
    return;
  }

  // get the name of the tracker
  std::string tracker_name = sand_reco::tracker::DigitCollection::getDigits().begin()->det;

  // build the clusters starting from digits and return them as a collection
  sand_reco::tracker::ClusterCollection clusters(sand_geo, digit_vec, sand_reco::tracker::ClusterCollection::ClusteringMethod::kCellAdjacency);
  
  // map
  // key  : z of the module, 
  // value: vector of tracklet reconstructed in the module
  std::map<double, std::vector<Tracklet>> z_to_tracklets;

  // init the SAND Tracker Utils
  SANDTrackerUtils::init(sand_geo->getTGeoManager());

  // random number generator
  TRandom3 rand(0);

  // loop over tracker modules
  for (const auto& container:clusters.getContainers()) {
    // loop over clusters
    for (const auto& cluster_in_container:container->getClusters()) {

      TVector3 first_point;
      TVector3 last_point;
      TVector3 first_point_momentum;
      TVector3 last_point_momentum;
      double min_z = 10e8;
      double max_z = -10e8;
      // loop over digits in each clusters
      for (uint d = 0; d < cluster_in_container.getDigits().size(); d++) {
        auto digit = sand_reco::tracker::DigitCollection::getDigit(cluster_in_container.getDigits()[d]);
        
        // find digit with min and max z
        if (digit.z > max_z) {
          max_z = digit.z;
          last_point = TVector3(digit.x, digit.y, digit.z);
          last_point_momentum = TVector3(digit.px, digit.py, digit.pz);
        }
        if (digit.z < min_z) {
          min_z = digit.z;
          first_point = TVector3(digit.x, digit.y, digit.z);
          first_point_momentum = TVector3(digit.px, digit.py, digit.pz);
        }
      }

      // eval position, direction and momentum of the particle at the z of the module
      auto true_tracklet = getTrueTrackletOfCluster(first_point, last_point, first_point_momentum, last_point_momentum, cluster_in_container.getZ());
      TVector3 true_pos = true_tracklet.pos_;
      TVector3 true_dir = true_tracklet.dir_;
      TVector3 true_mom = true_tracklet.mom_;

      // eval yz nd xz angles
      double true_theta_yz = atan(true_dir.Y() / true_dir.Z());
      double true_theta_xz = atan(true_dir.X() / true_dir.Z());
      if (true_theta_xz > M_PI_2) true_theta_xz -= M_PI;

      // smear tracket to simulate the measurement
      Tracklet measurement_from_true_tracklet;
      measurement_from_true_tracklet.x = true_pos.X()  + rand.Gaus() * SANDTrackerUtils::getSigmaPositionMeasurement() * 1E3;
      measurement_from_true_tracklet.y = true_pos.Y()  + rand.Gaus() * SANDTrackerUtils::getSigmaPositionMeasurement() * 1E3;
      measurement_from_true_tracklet.theta_xz = true_theta_xz + rand.Gaus(0, SANDTrackerUtils::getSigmaAngleMeasurement());
      measurement_from_true_tracklet.theta_yz = true_theta_yz + rand.Gaus(0, SANDTrackerUtils::getSigmaAngleMeasurement());
      for (uint d = 0; d < cluster_in_container.getDigits().size(); d++) {
        auto digit = sand_reco::tracker::DigitCollection::getDigit(cluster_in_container.getDigits()[d]);
        measurement_from_true_tracklet.digits.push_back(digit);
      }
      // std::cout << true_pos.X() << std::endl;
      measurement_from_true_tracklet.true_pos_ = true_pos;
      measurement_from_true_tracklet.true_dir_ = true_dir;
      measurement_from_true_tracklet.true_mom_ = true_mom;
      
      z_to_tracklets[cluster_in_container.getZ()].push_back(measurement_from_true_tracklet);
    }
  }
  
  if (z_to_tracklets.empty()) {
    return;
  }

  EDEPTree tree;
  tree.InizializeFromEdep(*mc_event, sand_geo->getTGeoManager());
  
  std::vector<EDEPTrajectory> primaryTrj;
  tree.Filter(std::back_insert_iterator<std::vector<EDEPTrajectory>>(primaryTrj), 
    [](const EDEPTrajectory& trj) { return trj.GetParentId() == -1;} );

  TDatabasePDG pdg_db;
  std::vector<SParticleInfo> particleInfos;
  std::map<double, std::vector<TVectorD>> z_to_best_tracklet;

  std::vector<int> indeces;
  int ii = -1;
  for (auto trj:primaryTrj) {
    ii++;

    if (trj.GetHitMap().find(string_to_component[tracker_name]) == trj.GetHitMap().end()) {
      continue;
    }
    
    if (trj.GetTrajectoryPoints().find(string_to_component[tracker_name]) == trj.GetTrajectoryPoints().end()) {
      continue;
    }

    auto particle = pdg_db.GetParticle(trj.GetPDGCode());

    if (!particle) {
      continue;
    }

    if (particle->Mass() == 0 || particle->Charge() == 0) {
      continue;
    }

    SParticleInfo pi;
    pi.pdg_code = trj.GetPDGCode();
    pi.id       = trj.GetId();
    pi.mass     = particle->Mass();
    pi.charge   = particle->Charge() / 3;

    double max_z = 0;
    bool to_be_reconstructed = false;

    for (auto& point : trj.GetTrajectoryPoints().at(string_to_component[tracker_name])) {
      if (point.GetPosition().Z() > max_z && point.GetMomentum().Z() > 100) {
        max_z = point.GetPosition().Z();
        pi.pos = point.GetPosition().Vect();
        pi.mom = point.GetMomentum();
        to_be_reconstructed = true;
      }
    }

    if (!to_be_reconstructed) continue;

    double sigma_pos = rand.Gaus(0, SANDTrackerUtils::getSigmaPositionMeasurement() * 1E3);
    double sigma_mom = 0.05;

    double x_smeared = rand.Gaus(pi.pos.X(), sigma_pos);
    double y_smeared = rand.Gaus(pi.pos.Y(), sigma_pos);
    double px_smeared = pi.mom.X() * rand.Gaus(1, sigma_mom);
    double py_smeared = pi.mom.Y() * rand.Gaus(1, sigma_mom);
    double pz_smeared = pi.mom.Z() * rand.Gaus(1, sigma_mom);

    pi.pos = TVector3(x_smeared, y_smeared, pi.pos.Z());
    pi.mom = TVector3(px_smeared, py_smeared, pz_smeared);
    pi.initial_pos = trj.GetTrajectoryPoints().at(string_to_component[tracker_name])[0].GetPosition().Vect();
    pi.initial_mom = trj.GetTrajectoryPoints().at(string_to_component[tracker_name])[0].GetMomentum();
    particleInfos.push_back(pi);
    indeces.push_back(ii);
  }

   
  
  int nParticles = particleInfos.size();
  
  if (nParticles == 0) {
    std::cerr << "no particles to be reconstructed...process aborted"
              << std::endl;
    return;
  }

  for (int ip = 0; ip < (int)indeces.size(); ip++) {
    std::string name_mg = "YZ_" + std::to_string(indeces[ip]);
    TMultiGraph* mg = new TMultiGraph(name_mg.c_str(), name_mg.c_str());
    TGraph* yz_true = new TGraph(primaryTrj[indeces[ip]].GetTrajectoryPoints().at(string_to_component[tracker_name]).size());
    std::string name_mgx = "XZ_" + std::to_string(indeces[ip]);
    TMultiGraph* mgx = new TMultiGraph(name_mgx.c_str(), name_mgx.c_str());
    TGraph* xz_true = new TGraph(primaryTrj[indeces[ip]].GetTrajectoryPoints().at(string_to_component[tracker_name]).size());

    for (uint i = 0; i <  primaryTrj[indeces[ip]].GetTrajectoryPoints().at(string_to_component[tracker_name]).size(); i++){
      auto point = primaryTrj[indeces[ip]].GetTrajectoryPoints().at(string_to_component[tracker_name])[i];
       yz_true->SetPoint(i, point.GetPosition().Z() , point.GetPosition().Y());
       xz_true->SetPoint(i, point.GetPosition().Z() , point.GetPosition().X());
    }


    tryCompleteManager(z_to_tracklets, particleInfos[ip], mg, mgx, plot_container);

    std::string title = name_mg + "; z [mm]; y [mm]";
    mg->SetTitle(title.c_str());
    yz_true->SetMarkerStyle(4);
    mg->Add(yz_true);
    mg->Write();
    title = name_mgx + "; z [mm]; x [mm]";
    mgx->SetTitle(title.c_str());
    xz_true->SetMarkerStyle(4);
    mgx->Add(xz_true);
    mgx->Write();
  }

}

int main(int argc, char* argv[])
{
  gStyle->SetOptStat(0);  

  TFile f(argv[1], "READ");
  TGeoManager* geo = 0;
  geo = (TGeoManager*)f.Get("EDepSimGeometry");


  // MC info tree
  TTree* t_h = (TTree*)f.Get("EDepSimEvents");
  TG4Event* ev = new TG4Event;
  t_h->SetBranchAddress("Event", &ev);
  
  TFile f_d(argv[2], "READ");
  TTree* t = (TTree*)f_d.Get("tDigit");

  std::vector<dg_wire>* digits = 0;
  t->SetBranchAddress("dg_wire", &digits);

  bool plots = false;
  
  TFile* pull_test = new TFile("pull_test.root", "RECREATE");
  PlotContainer plot_container;
  
  SANDGeoManager sand_geo;
  sand_geo.init(geo);
  
  std::string geometry;
  if (geo->FindVolumeFast("STTtracker_PV")) {
    geometry = "STT";
  } else if (geo->FindVolumeFast("SANDtracker_PV")) {
    geometry = "DRIFT";
  } 
  sand_geo.fillAdjacentCells(geometry);

  TFile* h_out = new TFile("h_out.root", "RECREATE");

  int nev = t->GetEntries();
  for (int i = 0; i < nev; i++) {
    t_h->GetEntry(i);
    t->GetEntry(i);

    if (!plots) {
      pull_test->cd();
      processEventWithKF(&sand_geo, ev, digits, plot_container);
    }

    // if (plots) {
    //   h_out->cd();
    //   int p[9] = {100, -2000, 2000, 100, -4000, -0, 100, 22500, 26000};

    //   sand_reco::tracker::DigitCollection::fillMap(digits);
    //   sand_reco::tracker::ClusterCollection clusters(&sand_geo, sand_reco::tracker::DigitCollection::getDigits(), sand_reco::tracker::ClusterCollection::ClusteringMethod::kCellAdjacency);
    //   auto digit_vec =  sand_reco::tracker::DigitCollection::getDigits();
    //   std::string tracker_name = digit_vec.begin()->det;

    //   TCanvas* canvas_cluster = new TCanvas("canvas_cluster","canvas_cluster",2000,1000);
    //   canvas_cluster->Divide(2,1);
      
    //   TH2D* h_cluster_yz = new TH2D("h","h", p[6],p[7], p[8], p[3],p[4], p[5]);
    //   TH2D* h_cluster_xz = new TH2D("h","h", p[6],p[7], p[8], p[0],p[1], p[2]);
    //   canvas_cluster->cd(1);
    //   h_cluster_yz->Draw();
    //   canvas_cluster->cd(2);
    //   h_cluster_xz->Draw();
    //   canvas_cluster->Print("clu.pdf(","pdf");

    //   std::map<double, std::vector<TVectorD>> z_to_tracklets;

    //   int color = 2;
    //   for (const auto& container:clusters.getContainers()) {
    //     int gg = 0;
    //     for (const auto& cluster_in_container:container->getClusters()) {
    //       gg++;
    //       if (color > 9) color = 2;
          
    //       EDEPTree tree;
    //       tree.InizializeFromEdep(*ev, sand_geo.getTGeoManager());
          
    //       std::vector<EDEPTrajectory> primaryTrj;
    //       tree.Filter(std::back_insert_iterator<std::vector<EDEPTrajectory>>(primaryTrj), 
    //         [](const EDEPTrajectory& trj) { return trj.GetParentId() == -1;} );

    //       for (auto trj:primaryTrj) {
            
    //         if (trj.GetTrajectoryPoints().find(string_to_component[tracker_name]) == trj.GetTrajectoryPoints().end()) {
    //           continue;
    //         }
            
    //         for (auto& point : trj.GetTrajectoryPoints().at(string_to_component[tracker_name])) {
    //           TEllipse* pt_yz = new TEllipse(point.GetPosition().Z(), point.GetPosition().Y(), 1);
    //           TEllipse* pt_xz = new TEllipse(point.GetPosition().Z(), point.GetPosition().X(), 1);
    //           pt_yz->SetFillStyle(0);
    //           pt_yz->SetLineWidth(1);
    //           pt_yz->SetLineColor(1);
    //           pt_xz->SetFillStyle(0);
    //           pt_xz->SetLineWidth(1);
    //           pt_xz->SetLineColor(1);
    //           canvas_cluster->cd(1);
    //           pt_yz->Draw();
    //           canvas_cluster->cd(2);
    //           pt_xz->Draw();
    //         }
    //       }

    //       TVector3 first_point;
    //       TVector3 last_point;
    //       double min_z = 10e8;
    //       double max_z = -10e8;
    //       for (uint d = 0; d < cluster_in_container.getDigits().size(); d++) {
    //         auto digit = sand_reco::tracker::DigitCollection::getDigit(cluster_in_container.getDigits()[d]);
            
    //         if (digit.z > max_z) {
    //           max_z = digit.z;
    //           last_point = TVector3(digit.x, digit.y, digit.z);
    //         }
    //         if (digit.z < min_z) {
    //           min_z = digit.z;
    //           first_point = TVector3(digit.x, digit.y, digit.z);
    //         }
    //       }

    //       auto true_tracklet = getTrueTrackletOfCluster(first_point, last_point, TVector3(0,0,0),  TVector3 (0,0,0), cluster_in_container.getZ());
    //       double z_start = cluster_in_container.getZ();
    //       // Draw tracklets
    //       TVector2 start_true_tracklet_yz(z_start, true_tracklet[0].Y());
    //       TVector2 start_true_tracklet_xz(z_start, true_tracklet[0].X());
    //       double zy_end = z_start + 5 * cos(atan(true_tracklet[1].Y() / true_tracklet[1].Z()));
    //       double zx_end = z_start + 5 * cos(atan(true_tracklet[1].X() / true_tracklet[1].Z()));
    //       double y_end = true_tracklet[0].Y() + 5 * sin(atan(true_tracklet[1].Y() / true_tracklet[1].Z()));
    //       double x_end = true_tracklet[0].X() + 5 * sin(atan(true_tracklet[1].X() / true_tracklet[1].Z()));
    //       TVector2 end_true_tracklet_yz(zy_end, y_end);
    //       TVector2 end_true_tracklet_xz(zx_end, x_end);
          
    //       TLine* line_yz_true_tracklet = new TLine(start_true_tracklet_yz.X(), start_true_tracklet_yz.Y(), end_true_tracklet_yz.X(), end_true_tracklet_yz.Y());
    //       TLine* line_xz_true_tracklet = new TLine(start_true_tracklet_xz.X(), start_true_tracklet_xz.Y(), end_true_tracklet_xz.X(), end_true_tracklet_xz.Y());
    //       line_yz_true_tracklet->SetLineColor(color + 1);
    //       line_yz_true_tracklet->SetLineWidth(3);
    //       line_xz_true_tracklet->SetLineColor(color + 1);
    //       line_xz_true_tracklet->SetLineWidth(3);
          
    //       canvas_cluster->cd(1);
    //       line_yz_true_tracklet->Draw();
    //       canvas_cluster->cd(2);
    //       line_xz_true_tracklet->Draw();
        


    //       for (auto digit:digit_vec) {
    //         auto cell = sand_geo.getCellInfo(sand_geometry::tracker::CellID(digit.did));
    //         auto h = cell->second.getSize().h;
    //         auto w = cell->second.getSize().w;
    //         TBox* box_yz = new TBox(cell->second.getWire().getCenter().Z() - h/2., cell->second.getWire().getCenter().Y() - w/2., cell->second.getWire().getCenter().Z() + h/2., cell->second.getWire().getCenter().Y() + w/2.);
    //         TBox* box_xz = new TBox(cell->second.getWire().getCenter().Z() - h/2., cell->second.getWire().getCenter().X() - w/2., cell->second.getWire().getCenter().Z() + h/2., cell->second.getWire().getCenter().X() + w/2.);
    //         box_yz->SetFillStyle(0);
    //         box_yz->SetLineColor(1);
    //         box_yz->SetLineWidth(1);
    //         box_xz->SetFillStyle(0);
    //         box_xz->SetLineColor(1);
    //         box_xz->SetLineWidth(1);
    //         canvas_cluster->cd(1);
    //         box_yz->Draw();
    //         canvas_cluster->cd(2);
    //         box_xz->Draw();

    //         TMarker* mark_yz = new TMarker(cell->second.getWire().getCenter().Z(), cell->second.getWire().getCenter().Y(), 5);
    //         mark_yz->SetMarkerColor(1);
    //         mark_yz->SetMarkerSize(0.5);
    //         canvas_cluster->cd(1);
    //         mark_yz->Draw();

    //         TMarker* mark_xz = new TMarker(cell->second.getWire().getCenter().Z(), cell->second.getWire().getCenter().X(), 5);
    //         mark_xz->SetMarkerColor(1);
    //         mark_xz->SetMarkerSize(0.5);
    //         canvas_cluster->cd(2);
    //         mark_xz->Draw();


    //         for (auto& kk:digit.hindex) {
    //           const TG4HitSegment& hseg = ev->SegmentDetectors[digit.det].at(kk);
    //           TLine* l_yz = new TLine(hseg.Start.Z(), hseg.Start.Y(), hseg.Stop.Z(), hseg.Stop.Y());
    //           TLine* l_xz = new TLine(hseg.Start.Z(), hseg.Start.X(), hseg.Stop.Z(), hseg.Stop.X());
    //           l_yz->SetLineColor(1);
    //           l_xz->SetLineColor(1);
    //           canvas_cluster->cd(1);
    //           l_yz->Draw();
    //           canvas_cluster->cd(2);
    //           l_xz->Draw();
    //         }
    //       }

    //       std::vector<sand_reco::tracker::DigitID> digits_cluster = cluster_in_container.getDigits();

    //       for (uint d = 0; d < digits_cluster.size(); d++) {
    //         canvas_cluster->cd();

    //         auto digit = sand_reco::tracker::DigitCollection::getDigit(digits_cluster[d]);
    //         auto cell = sand_geo.getCellInfo(sand_geometry::tracker::CellID(digit.did));

            
    //         // Draw cells of cluster
    //         auto h = cell->second.getSize().h;
    //         auto w = cell->second.getSize().w;

    //         TBox* box_yz = new TBox(cell->second.getWire().getCenter().Z() - h/2., cell->second.getWire().getCenter().Y() - w/2., cell->second.getWire().getCenter().Z() + h/2., cell->second.getWire().getCenter().Y() + w/2.);
    //         TBox* box_xz = new TBox(cell->second.getWire().getCenter().Z() - h/2., cell->second.getWire().getCenter().X() - w/2., cell->second.getWire().getCenter().Z() + h/2., cell->second.getWire().getCenter().X() + w/2.);
    //         box_yz->SetFillStyle(0);
    //         box_yz->SetLineColor(color);
    //         box_yz->SetLineWidth(1);
    //         box_xz->SetFillStyle(0);
    //         box_xz->SetLineColor(color);
    //         box_xz->SetLineWidth(1);
    //         canvas_cluster->cd(1);
    //         box_yz->Draw();
    //         canvas_cluster->cd(2);
    //         box_xz->Draw();
            
    //         // Draw true drift time of digits in cluster
    //         TEllipse* el_yz = new TEllipse(cell->second.getWire().getCenter().Z(), cell->second.getWire().getCenter().Y(), 
    //                               sand_reco::stt::wire_radius + cell->second.getDriftVelocity() * digit.drift_time);
    //         TEllipse* el_xz = new TEllipse(cell->second.getWire().getCenter().Z(), cell->second.getWire().getCenter().X(), 
    //                               sand_reco::stt::wire_radius + cell->second.getDriftVelocity() * digit.drift_time);
    //         el_yz->SetFillStyle(0);
    //         el_yz->SetLineWidth(1);
    //         el_yz->SetLineColor(1);
    //         el_xz->SetFillStyle(0);
    //         el_xz->SetLineWidth(1);
    //         el_xz->SetLineColor(1);
    //         canvas_cluster->cd(1);
    //         el_yz->Draw();
    //         canvas_cluster->cd(2);
    //         el_xz->Draw();

    //         // Draw hit segments for the cluster
    //         for (auto& kk:digit.hindex) {
    //           const TG4HitSegment& hseg = ev->SegmentDetectors[digit.det].at(kk);
    //           TLine* l_yz = new TLine(hseg.Start.Z(), hseg.Start.Y(), hseg.Stop.Z(), hseg.Stop.Y());
    //           TLine* l_xz = new TLine(hseg.Start.Z(), hseg.Start.X(), hseg.Stop.Z(), hseg.Stop.X());
    //           l_yz->SetLineColor(color);
    //           l_xz->SetLineColor(color);
    //           canvas_cluster->cd(1);
    //           l_yz->Draw();
    //           canvas_cluster->cd(2);
    //           l_xz->Draw();
    //         }
    //       }

    //       color++; 
    //       canvas_cluster->Write();
    //       canvas_cluster->Print("clu.pdf","pdf");
    //       canvas_cluster->Clear();

    //       canvas_cluster->Divide(2,1);
    //       canvas_cluster->cd(1);
    //       h_cluster_yz->Draw();
    //       canvas_cluster->cd(2);
    //       h_cluster_xz->Draw();

    //     }
    //   }
    //   canvas_cluster->Print("clu.pdf)","pdf");

    // }
  }

  if (!plots) {
    pull_test->cd();
    plot_container.meas_x_res.Write();
    plot_container.meas_y_res.Write();
    plot_container.smooth_x_res.Write();
    plot_container.smooth_y_res.Write();
    plot_container.smooth_p_res.Write();
    plot_container.smooth_tan_res.Write();
    plot_container.smooth_phi_res.Write();
    plot_container.x_pull_smooth.Write();
    plot_container.y_pull_smooth.Write();
    plot_container.r_inv_pull_smooth.Write();
    plot_container.tan_pull_smooth.Write();
    plot_container.phi_pull_smooth.Write();
    plot_container.p_res_step.Write();
    plot_container.p_init_vs_reco_smooth.Write();
    plot_container.p_first_vs_reco_smooth.Write();
    plot_container.p_diff_vs_points.Write();
    plot_container.h_gpos_distribution.Write();
    plot_container.h_gang_distribution.Write();
    plot_container.chi2.Write();
    plot_container.deltaP_over_P.Write();
    pull_test->Close();
  }
}
