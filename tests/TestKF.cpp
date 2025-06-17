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
#include "utils.h"

//#include "EDEPTree.h"

void tryCompleteManager(sand_reco::kf::TrackletMap z_to_tracklets, SParticleInfo particle, TH1D* h_gpos_distribution, TH1D* h_gang_distribution, TMultiGraph* mg) {
  sand_reco::kf::Manager manager;
  manager.initFromMC(&z_to_tracklets, particle);
  manager.run();

  auto track = manager.getTrack();
  if (track.getSteps().size() > 3) {
    std::cout << track.getSteps().size() << std::endl;
    auto last_step = track.getSteps().back();
    auto reco_state =
          last_step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kSmoothing).getStateVector();
    auto reco_mom = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(
                                  reco_state.radius(), reco_state.tanLambda());

    std::cout << "Initial Smoothed Reco Momentum " << reco_mom << std::endl;
    
    TGraph* yz_predicted = new TGraph(track.getSteps().size());
    TGraph* yz_filtered = new TGraph(track.getSteps().size());
    TGraph* yz_smoothed = new TGraph(track.getSteps().size());
    TGraph* yz_measured = new TGraph(track.getSteps().size());
    
    int i = 0;
    for (auto& step : track.getSteps()) {
      auto prediction = step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kPrediction).getStateVector();
      auto filtering = step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kFiltering).getStateVector();
      auto smoothing =  step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kSmoothing).getStateVector();
      
      yz_predicted->SetPoint(i, step.getZ(), prediction.y()*1000 );
      yz_filtered->SetPoint(i, step.getZ() , filtering.y()*1000);
      yz_smoothed->SetPoint(i, step.getZ() , smoothing.y()*1000);
      yz_measured->SetPoint(i, step.getZ() , step.getY());
      i++;
  
      auto tanLambda = step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kFiltering).getStateVector().tanLambda();
      
      auto& innovation = step.getInnovation();
      if (innovation.empty()) {
        continue;
      }
  
      h_gpos_distribution->Fill(innovation[0]);
      h_gang_distribution->Fill(innovation[1]);
    
      
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
    }
  }

  return;
}

void processEventWithKF(SANDGeoManager* sand_geo, TG4Event* mc_event, std::vector<dg_wire>* digits, TH1D* h_gpos_distribution,TH1D* h_gang_distribution)
{
  
  int p[9] = {100, -2000, 2000, 100, -4000, -1000, 100, 23800, 26000};

  sand_reco::tracker::DigitCollection::fillMap(digits);
  auto digit_map =  sand_reco::tracker::DigitCollection::getDigits();
  if (sand_reco::tracker::DigitCollection::getDigits().empty()) {
    return;
  }
  std::string tracker_name = sand_reco::tracker::DigitCollection::getDigits().begin()->det;
  sand_reco::tracker::ClusterCollection clusters(sand_geo, sand_reco::tracker::DigitCollection::getDigits(), sand_reco::tracker::ClusterCollection::ClusteringMethod::kCellAdjacency);
  
  TrackletFinder traklet_finder;
  traklet_finder.setVolumeParameters(p);
  traklet_finder.setSigmaPosition(0.2);
  traklet_finder.setSigmaAngle(0.02);

  std::map<double, std::vector<TVectorD>> z_to_tracklets;

  SANDTrackerUtils::init(sand_geo->getTGeoManager());

  int gg = 0;
  for (const auto& container:clusters.getContainers()) {
    for (const auto& cluster_in_container:container->getClusters()) {
      // std::cout << 100 * gg / container->getClusters().size() << std::endl;
      // gg++;
      // if (gg == 100) break;
      // if(cluster_in_container.getZ() < 25650) continue;

      traklet_finder.setCells(cluster_in_container);
      auto minima = traklet_finder.findTracklets();

      double z_start = cluster_in_container.getZ();
      for (uint trk = 0; trk < minima.size(); trk++) {
        if (minima[trk][4] < 1E-2) {
          z_to_tracklets[cluster_in_container.getZ()].push_back(minima[trk]);
        }
      }
      traklet_finder.clear();
    }
  }
  
  for (auto z : z_to_tracklets) std::cout << z.first << std::endl;
  
  int sum = 0;
  for (auto el:z_to_tracklets) {
    sum += el.second.size();
  }
  if (sum == 0) {
    return;
  }

  EDEPTree tree;
  tree.InizializeFromEdep(*mc_event, sand_geo->getTGeoManager());
  
  std::vector<EDEPTrajectory> primaryTrj;
  tree.Filter(std::back_insert_iterator<std::vector<EDEPTrajectory>>(primaryTrj), 
    [](const EDEPTrajectory& trj) { return trj.GetParentId() == -1;} );

  TDatabasePDG pdg_db;
  std::vector<SParticleInfo> particleInfos;
  TRandom3 rand(0);
  std::map<double, std::vector<TVectorD>> z_to_best_tracklet;

  double sigma_pos = 0;
  double sigma_mom = 0;
  for (auto trj:primaryTrj) {

    if (trj.GetHitMap().find(string_to_component[tracker_name]) == trj.GetHitMap().end()) {
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
    pi.mass = particle->Mass();
    pi.charge = particle->Charge() / 3;

    // std::cout << pi.mass << " " << pi.charge << " " << pi.pdg_code << std::endl;
    // std::cout << trj.GetHitMap().at(string_to_component[tracker_name]).size() << std::endl;

    double max_z = 0;
    bool to_be_reconstructed = false;
    for (auto& point : trj.GetTrajectoryPoints().at(string_to_component[tracker_name])) {
      // std::cout << point.GetPosition().Z() << " " << point.GetMomentum().Mag() << " " << point.GetMomentum().Z() << std::endl;
      if (point.GetPosition().Z() > max_z && point.GetMomentum().Z() > 100) {
        max_z = point.GetPosition().Z();
        pi.pos = point.GetPosition().Vect();
        pi.mom = point.GetMomentum();
        to_be_reconstructed = true;
      }
    }

    if (!to_be_reconstructed) continue;
    double x_smeared = rand.Gaus(pi.pos.X(), sigma_pos);
    double y_smeared = rand.Gaus(pi.pos.Y(), sigma_pos);
    double px_smeared = pi.mom.X() * rand.Gaus(1, sigma_mom);
    double py_smeared = pi.mom.Y() * rand.Gaus(1, sigma_mom);
    double pz_smeared = pi.mom.Z() * rand.Gaus(1, sigma_mom);

    pi.pos = TVector3(x_smeared, y_smeared, pi.pos.Z());
    pi.mom = TVector3(px_smeared, py_smeared, pz_smeared);
    particleInfos.push_back(pi);


    std::cout << "Initial Momentum " << trj.GetInitialMomentum().Vect().Mag() << std::endl;
    std::cout << "Selected Momentum " << pi.mom.Mag() << " " << pi.mom.Z() << std::endl;

    std::map<double, std::vector<TVector3>> z_to_interpolated_tracklets = getInterpolatedZ(trj.GetTrajectoryPoints().at(string_to_component[tracker_name]), z_to_tracklets);
    z_to_best_tracklet = findBestTracklet(z_to_tracklets, z_to_interpolated_tracklets);
  }
  
  int nParticles = particleInfos.size();
  
  if (nParticles == 0) {
    std::cerr << "no particles to be reconstructed...process aborted"
              << std::endl;
    return;
  }

  for (int ip = 0; ip < nParticles; ip++) {
    std::string name_mg = "YZ_" + std::to_string(ip);
    TMultiGraph* mg = new TMultiGraph(name_mg.c_str(), name_mg.c_str());
    TGraph* yz_true = new TGraph(primaryTrj[ip].GetTrajectoryPoints().at(string_to_component[tracker_name]).size());

    for (uint i = 0; i <  primaryTrj[ip].GetTrajectoryPoints().at(string_to_component[tracker_name]).size(); i++){
      auto point = primaryTrj[ip].GetTrajectoryPoints().at(string_to_component[tracker_name])[i];
       yz_true->SetPoint(i, point.GetPosition().Z() , point.GetPosition().Y());
    }

    bool use_interpolated = true;
    if (use_interpolated) {
      tryCompleteManager(z_to_best_tracklet, particleInfos[ip], h_gpos_distribution, h_gang_distribution, mg);
    } else {
      tryCompleteManager(z_to_tracklets, particleInfos[ip], h_gpos_distribution, h_gang_distribution, mg);
    }

    mg->SetTitle("YZ view; z [mm]; y [mm]");
    yz_true->SetMarkerStyle(4);
    mg->Add(yz_true);
    mg->Write();
  }
}

int main(int argc, char* argv[])
{
  gStyle->SetOptStat(0);  

  TFile f(argv[2], "READ");
  TGeoManager* geo = 0;
  geo = (TGeoManager*)f.Get("EDepSimGeometry");


  // MC info tree
  TTree* t_h = (TTree*)f.Get("EDepSimEvents");
  TG4Event* ev = new TG4Event;
  t_h->SetBranchAddress("Event", &ev);
  
  TFile f_d(argv[3], "READ");
  TTree* t = (TTree*)f_d.Get("tDigit");

  std::vector<dg_wire>* digits = 0;
  t->SetBranchAddress("dg_wire", &digits);

    
  TFile* innovation_test = new TFile("innovation_test.root", "RECREATE");
  TH1D* h_gpos_distribution = new TH1D("h_gpos_distribution", "Innovation", 100, -3, 3);
  TH1D* h_gang_distribution = new TH1D("h_gang_distribution", "Innovation", 100, -3, 3);
  SANDGeoManager sand_geo;
  sand_geo.init(geo);
  
  std::string geometry;
  if (geo->FindVolumeFast("STTtracker_PV")) {
    geometry = "STT";
  } else if (geo->FindVolumeFast("SANDtracker_PV")) {
    geometry = "DRIFT";
  } 
  sand_geo.fillAdjacentCells(geometry);

  for (int i = 0; i < 20; i++) {
    t_h->GetEntry(i);
    t->GetEntry(i);

    processEventWithKF(&sand_geo, ev, digits, h_gpos_distribution, h_gang_distribution);
  }
  h_gpos_distribution->Write();
  h_gang_distribution->Write();
  innovation_test->Close();
}
