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

void tryCompleteManager(sand_reco::kf::TrackletMap z_to_tracklets, SParticleInfo particle, TH1D* h_gpos_distribution, TH1D* h_gang_distribution, TMultiGraph* mg, TMultiGraph* mgx) {
  sand_reco::kf::Manager manager;
  manager.initFromMC(&z_to_tracklets, particle);
  manager.run();

  auto track = manager.getTrack();
   if (track.getSteps().size() > 0) { // was commented, with > 3
    std::cout << track.getSteps().size() << std::endl;
    auto last_step = track.getSteps().back(); //crash if empty due to the .back().
    auto reco_state =
          last_step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kSmoothing).getStateVector();
    auto reco_mom = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(
                                  reco_state.radius(), reco_state.tanLambda());

    std::cout << "Initial Smoothed Reco Momentum " << reco_mom << std::endl;
    
    TGraph* yz_predicted = new TGraph(track.getSteps().size());
    TGraph* yz_filtered = new TGraph(track.getSteps().size());
    TGraph* yz_smoothed = new TGraph(track.getSteps().size());
    TGraph* yz_measured = new TGraph(track.getSteps().size());
    TGraph* xz_predicted = new TGraph(track.getSteps().size());
    TGraph* xz_filtered = new TGraph(track.getSteps().size());
    TGraph* xz_smoothed = new TGraph(track.getSteps().size());
    TGraph* xz_measured = new TGraph(track.getSteps().size());
    
    int i = 0;
    for (auto& step : track.getSteps()) {
      auto prediction = step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kPrediction).getStateVector();
      auto filtering = step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kFiltering).getStateVector();
      auto smoothing =  step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kSmoothing).getStateVector();
      
      yz_predicted->SetPoint(i, step.getZ(), prediction.y()*1000 );
      yz_filtered->SetPoint(i, step.getZ() , filtering.y()*1000);
      yz_smoothed->SetPoint(i, step.getZ() , smoothing.y()*1000);
      yz_measured->SetPoint(i, step.getZ() , step.getY());
      xz_predicted->SetPoint(i, step.getZ(), prediction.x()*1000 );
      xz_filtered->SetPoint(i, step.getZ() , filtering.x()*1000);
      xz_smoothed->SetPoint(i, step.getZ() , smoothing.x()*1000);
      xz_measured->SetPoint(i, step.getZ() , step.getX());
      i++;
      
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

void processEventWithKF(SANDGeoManager* sand_geo, TG4Event* mc_event, std::vector<dg_wire>* digits, TH1D* h_gpos_distribution,TH1D* h_gang_distribution,
                        TH1D* h_x_diff, TH1D* h_y_diff, TH1D* h_theta_x_diff, TH1D* h_theta_y_diff)
{
  
  int p[9] = {100, -2000, 2000, 100, -4000, -0, 100, 22500, 26000};

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

      TVector3 first_point;
      TVector3 last_point;
      double min_z = 10e8;
      double max_z = -10e8;
      for (uint d = 0; d < cluster_in_container.getDigits().size(); d++) {
        auto digit = sand_reco::tracker::DigitCollection::getDigit(cluster_in_container.getDigits()[d]);
        
        if (digit.z > max_z) {
          max_z = digit.z;
          last_point = TVector3(digit.x, digit.y, digit.z);
        }
        if (digit.z < min_z) {
          min_z = digit.z;
          first_point = TVector3(digit.x, digit.y, digit.z);
        }
      }
      auto true_tracklet = getTrueTrackletOfCluster(first_point, last_point, cluster_in_container.getZ());

      traklet_finder.setCells(cluster_in_container);
      auto minima = traklet_finder.findTracklets();

      double best_score = 10e8;
      TVectorD best_tracklet(9);

                
      double z_start = cluster_in_container.getZ();
      for (uint trk = 0; trk < minima.size(); trk++) {
        // if (minima[trk][4] < 1E-2) {
          double score = getScore(minima[trk], true_tracklet);
          if (score < best_score) {
            best_tracklet = minima[trk];
            best_score = score;
          }
          // z_to_tracklets[cluster_in_container.getZ()].push_back(minima[trk]);
        // }
      }
      z_to_tracklets[cluster_in_container.getZ()].push_back(best_tracklet);

      h_x_diff->Fill(best_tracklet[0] - true_tracklet[0].X());
      h_y_diff->Fill(best_tracklet[1] - true_tracklet[0].Y());
      h_theta_x_diff->Fill(best_tracklet[2] - atan(true_tracklet[1].Z() / true_tracklet[1].X()));
      h_theta_y_diff->Fill(best_tracklet[3] - atan(true_tracklet[1].Y() / true_tracklet[1].Z()));

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

    std::cout << "HIT " << trj.GetHitMap().at(string_to_component[tracker_name]).size() << std::endl;

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

    for (auto t:z_to_interpolated_tracklets) {
      std::cout << "Interpolated z: " << t.first 
                << ", does the same z from tracklets exist?: " <<  (z_to_tracklets.find(t.first) != z_to_tracklets.end())  << std::endl;
    }

    std::vector<double> z_difference = computeZDistance(trj.GetTrajectoryPoints().at(string_to_component[tracker_name]), z_to_tracklets);
    
    for(auto z:z_difference){
      std::cout << "distnza fra due z è: " << z << std::endl;
    }
  
    
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
    std::string name_mgx = "XZ_" + std::to_string(ip);
    TMultiGraph* mgx = new TMultiGraph(name_mgx.c_str(), name_mgx.c_str());
    TGraph* xz_true = new TGraph(primaryTrj[ip].GetTrajectoryPoints().at(string_to_component[tracker_name]).size());

    for (uint i = 0; i <  primaryTrj[ip].GetTrajectoryPoints().at(string_to_component[tracker_name]).size(); i++){
      auto point = primaryTrj[ip].GetTrajectoryPoints().at(string_to_component[tracker_name])[i];
       yz_true->SetPoint(i, point.GetPosition().Z() , point.GetPosition().Y());
       xz_true->SetPoint(i, point.GetPosition().Z() , point.GetPosition().X());
    }

    bool use_interpolated = false;
    if (use_interpolated) {
      tryCompleteManager(z_to_best_tracklet, particleInfos[ip], h_gpos_distribution, h_gang_distribution, mg, mgx);
    } else {
      tryCompleteManager(z_to_tracklets, particleInfos[ip], h_gpos_distribution, h_gang_distribution, mg, mgx);
    }

    mg->SetTitle("YZ view; z [mm]; y [mm]");
    yz_true->SetMarkerStyle(4);
    mg->Add(yz_true);
    mg->Write();
    mgx->SetTitle("XZ view; z [mm]; x [mm]");
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

  bool plots = true;
  TFile* innovation_test = new TFile("innovation_test.root", "RECREATE");
  TH1D* h_gpos_distribution = new TH1D("h_gpos_distribution", "Innovation", 100, -3, 3);
  TH1D* h_gang_distribution = new TH1D("h_gang_distribution", "Innovation", 100, -3, 3);
  TH1D* h_x_diff = new TH1D("h_x_diff", "h_x_diff", 1000, -3, 3);
  TH1D* h_y_diff = new TH1D("h_y_diff", "h_y_diff", 1000, -3, 3);
  TH1D* h_theta_x_diff = new TH1D("h_theta_x_diff", "h_theta_x_diff", 1000, -3, 3);
  TH1D* h_theta_y_diff = new TH1D("h_theta_y_diff", "h_theta_y_diff", 1000, -3, 3);
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

  for (int i = 2; i < 3; i++) {
    t_h->GetEntry(i);
    t->GetEntry(i);

    if (!plots) {
      innovation_test->cd();
      processEventWithKF(&sand_geo, ev, digits, h_gpos_distribution, h_gang_distribution, h_x_diff, h_y_diff, h_theta_x_diff, h_theta_y_diff);
    }

    if (plots) {
      h_out->cd();
      int p[9] = {100, -2000, 2000, 100, -4000, -0, 100, 22500, 26000};

      sand_reco::tracker::DigitCollection::fillMap(digits);
      sand_reco::tracker::ClusterCollection clusters(&sand_geo, sand_reco::tracker::DigitCollection::getDigits(), sand_reco::tracker::ClusterCollection::ClusteringMethod::kCellAdjacency);
      auto digit_map =  sand_reco::tracker::DigitCollection::getDigits();
      std::string tracker_name = digit_map.begin()->det;
      
      TrackletFinder traklet_finder;
      traklet_finder.setVolumeParameters(p);
      traklet_finder.setSigmaPosition(0.2);
      traklet_finder.setSigmaAngle(0.02);
      

      TCanvas* canvas_cluster = new TCanvas("canvas_cluster","canvas_cluster",2000,1000);
      canvas_cluster->Divide(2,1);
      
      TH2D* h_cluster_yz = new TH2D("h","h", p[6],p[7], p[8], p[3],p[4], p[5]);
      TH2D* h_cluster_xz = new TH2D("h","h", p[6],p[7], p[8], p[0],p[1], p[2]);
      canvas_cluster->cd(1);
      h_cluster_yz->Draw();
      canvas_cluster->cd(2);
      h_cluster_xz->Draw();
      canvas_cluster->Print("clu.pdf(","pdf");

      std::map<double, std::vector<TVectorD>> z_to_tracklets;

      int color = 2;
      for (const auto& container:clusters.getContainers()) {
        int gg = 0;
        for (const auto& cluster_in_container:container->getClusters()) {
          std::cout << (double)gg / container->getClusters().size() * 100 << std::endl;
          gg++;
          if (color > 9) color = 2;
          
          TVector3 first_point;
          TVector3 last_point;
          double min_z = 10e8;
          double max_z = -10e8;
          for (uint d = 0; d < cluster_in_container.getDigits().size(); d++) {
            auto digit = sand_reco::tracker::DigitCollection::getDigit(cluster_in_container.getDigits()[d]);
            
            if (digit.z > max_z) {
              max_z = digit.z;
              last_point = TVector3(digit.x, digit.y, digit.z);
            }
            if (digit.z < min_z) {
              min_z = digit.z;
              first_point = TVector3(digit.x, digit.y, digit.z);
            }
          }

          auto true_tracklet = getTrueTrackletOfCluster(first_point, last_point, cluster_in_container.getZ());
          
          
          traklet_finder.setCells(cluster_in_container);
          auto minima = traklet_finder.findTracklets();

          // Draw tracklets
          if (minima.size() != 0) {
            canvas_cluster->cd();
            std::sort(minima.begin(), minima.end(),
            [](TVectorD v1, TVectorD v2){ return v1[4] < v2[4];});
            double z_start = cluster_in_container.getZ();

            double best_score = 10e8;
            TVectorD best_tracklet(9);

            for (uint trk = 0; trk < minima.size(); trk++) {
              if (minima[trk][4] < 1E-2) {
                
                double score = getScore(minima[trk], true_tracklet);
                std::cout << "SCORE: " << score << std::endl;
                if (score < best_score) {
                  best_tracklet = minima[trk];
                  best_tracklet.Print();
                  best_score = score;
                }
                
                z_to_tracklets[cluster_in_container.getZ()].push_back(minima[trk]);

                TVector2 start_tracklet_yz(z_start, minima[trk][1]);
                TVector2 start_tracklet_xz(z_start, minima[trk][0]);
                double zy_end = z_start + 5 * cos(minima[trk][3]);
                double zx_end = z_start + 5 * sin(minima[trk][2]);
                double y_end = minima[trk][1] + 5 * sin(minima[trk][3]);
                double x_end = minima[trk][0] + 5 * cos(minima[trk][2]);
                TVector2 end_tracklet_yz(zy_end, y_end);
                TVector2 end_tracklet_xz(zx_end, x_end);
                
                TLine* line_yz_tracklet = new TLine(start_tracklet_yz.X(), start_tracklet_yz.Y(), end_tracklet_yz.X(), end_tracklet_yz.Y());
                TLine* line_xz_tracklet = new TLine(start_tracklet_xz.X(), start_tracklet_xz.Y(), end_tracklet_xz.X(), end_tracklet_xz.Y());
                line_yz_tracklet->SetLineColor(color);
                line_yz_tracklet->SetLineWidth(1);
                line_xz_tracklet->SetLineColor(color);
                line_xz_tracklet->SetLineWidth(1);
                
                canvas_cluster->cd(1);
                line_yz_tracklet->Draw();
                canvas_cluster->cd(2);
                line_xz_tracklet->Draw();
              }
            }

            TVector2 start_tracklet_yz(z_start, best_tracklet[1]);
            TVector2 start_tracklet_xz(z_start, best_tracklet[0]);
            double zy_end = z_start + 5 * cos(best_tracklet[3]);
            double zx_end = z_start + 5 * sin(best_tracklet[2]);
            double y_end = best_tracklet[1] + 5 * sin(best_tracklet[3]);
            double x_end = best_tracklet[0] + 5 * cos(best_tracklet[2]);
            TVector2 end_tracklet_yz(zy_end, y_end);
            TVector2 end_tracklet_xz(zx_end, x_end);
            
            TLine* line_yz_tracklet = new TLine(start_tracklet_yz.X(), start_tracklet_yz.Y(), end_tracklet_yz.X(), end_tracklet_yz.Y());
            TLine* line_xz_tracklet = new TLine(start_tracklet_xz.X(), start_tracklet_xz.Y(), end_tracklet_xz.X(), end_tracklet_xz.Y());
            line_yz_tracklet->SetLineColor(color);
            line_yz_tracklet->SetLineWidth(3);
            line_xz_tracklet->SetLineColor(color);
            line_xz_tracklet->SetLineWidth(3);
            
            canvas_cluster->cd(1);
            line_yz_tracklet->Draw();
            canvas_cluster->cd(2);
            line_xz_tracklet->Draw();

            TVector2 start_true_tracklet_yz(z_start, true_tracklet[0].Y());
            TVector2 start_true_tracklet_xz(z_start, true_tracklet[0].X());
            zy_end = z_start + 5 * cos(atan(true_tracklet[1].Y() / true_tracklet[1].Z()));
            zx_end = z_start + 5 * cos(atan(true_tracklet[1].X() / true_tracklet[1].Z()));
            y_end = true_tracklet[0].Y() + 5 * sin(atan(true_tracklet[1].Y() / true_tracklet[1].Z()));
            x_end = true_tracklet[0].X() + 5 * sin(atan(true_tracklet[1].X() / true_tracklet[1].Z()));
            TVector2 end_true_tracklet_yz(zy_end, y_end);
            TVector2 end_true_tracklet_xz(zx_end, x_end);
            
            TLine* line_yz_true_tracklet = new TLine(start_true_tracklet_yz.X(), start_true_tracklet_yz.Y(), end_true_tracklet_yz.X(), end_true_tracklet_yz.Y());
            TLine* line_xz_true_tracklet = new TLine(start_true_tracklet_xz.X(), start_true_tracklet_xz.Y(), end_true_tracklet_xz.X(), end_true_tracklet_xz.Y());
            line_yz_true_tracklet->SetLineColor(color + 1);
            line_yz_true_tracklet->SetLineWidth(3);
            line_xz_true_tracklet->SetLineColor(color + 1);
            line_xz_true_tracklet->SetLineWidth(3);
            
            canvas_cluster->cd(1);
            line_yz_true_tracklet->Draw();
            canvas_cluster->cd(2);
            line_xz_true_tracklet->Draw();
          }

          bool ok = false;
          for (uint trk = 0; trk < minima.size(); trk++) {
            if (minima[trk][4] < 1E-2) {
              ok = true;
              break;
            }
          }  
          if(!ok) {
            traklet_finder.clear();
            continue;
          }

          auto digitId_to_drift_time = traklet_finder.getDigitToDriftTimeMap();
          
          // Draw cells of all digits
          // for (auto digit:digit_map) {
          //   auto cell = sand_geo.getCellInfo(sand_geometry::tracker::CellID(digit.did));
          //   auto h = cell->second.getSize().h;
          //   auto w = cell->second.getSize().w;


          //   TVector3 r = cell->second.getWire().getDirection();
          //   TVector3 leftend = cell->second.getWire().getReadoutPoint();

          //   TVector3 AP = TVector3(digit.x, digit.y, digit.z) - leftend; 
          //   double t_prime = AP.Dot(r) / r.Mag2();
          //   t_prime = std::max(0.0, std::min(1.0, t_prime));
          //   TVector3 position_along_wire = leftend + t_prime * r;

          //   TBox* box_yz = new TBox(position_along_wire.Z() - h/2., position_along_wire.Y() - w/2., position_along_wire.Z() + h/2., position_along_wire.Y() + w/2.);
          //   TBox* box_xz = new TBox(position_along_wire.Z() - h/2., position_along_wire.X() - w/2., position_along_wire.Z() + h/2., position_along_wire.X() + w/2.);
          //   box_yz->SetFillStyle(0);
          //   box_yz->SetLineColor(1);
          //   box_yz->SetLineWidth(1);
          //   box_xz->SetFillStyle(0);
          //   box_xz->SetLineColor(1);
          //   box_xz->SetLineWidth(1);
          //   canvas_cluster->cd(1);
          //   box_yz->Draw();
          //   canvas_cluster->cd(2);
          //   box_xz->Draw();

          //   TMarker* mark_yz = new TMarker(position_along_wire.Z(), position_along_wire.Y(), 5);
          //   mark_yz->SetMarkerColor(1);
          //   mark_yz->SetMarkerSize(0.5);
          //   canvas_cluster->cd(1);
          //   mark_yz->Draw();

          //   TMarker* mark_xz = new TMarker(position_along_wire.Z(), position_along_wire.X(), 5);
          //   mark_xz->SetMarkerColor(1);
          //   mark_xz->SetMarkerSize(0.5);
          //   canvas_cluster->cd(2);
          //   mark_xz->Draw();
          // }


          for (auto digit:digit_map) {
            auto cell = sand_geo.getCellInfo(sand_geometry::tracker::CellID(digit.did));
            auto h = cell->second.getSize().h;
            auto w = cell->second.getSize().w;
            TBox* box_yz = new TBox(cell->second.getWire().getCenter().Z() - h/2., cell->second.getWire().getCenter().Y() - w/2., cell->second.getWire().getCenter().Z() + h/2., cell->second.getWire().getCenter().Y() + w/2.);
            TBox* box_xz = new TBox(cell->second.getWire().getCenter().Z() - h/2., cell->second.getWire().getCenter().X() - w/2., cell->second.getWire().getCenter().Z() + h/2., cell->second.getWire().getCenter().X() + w/2.);
            box_yz->SetFillStyle(0);
            box_yz->SetLineColor(1);
            box_yz->SetLineWidth(1);
            box_xz->SetFillStyle(0);
            box_xz->SetLineColor(1);
            box_xz->SetLineWidth(1);
            canvas_cluster->cd(1);
            box_yz->Draw();
            canvas_cluster->cd(2);
            box_xz->Draw();

            TMarker* mark_yz = new TMarker(cell->second.getWire().getCenter().Z(), cell->second.getWire().getCenter().Y(), 5);
            mark_yz->SetMarkerColor(1);
            mark_yz->SetMarkerSize(0.5);
            canvas_cluster->cd(1);
            mark_yz->Draw();

            TMarker* mark_xz = new TMarker(cell->second.getWire().getCenter().Z(), cell->second.getWire().getCenter().X(), 5);
            mark_xz->SetMarkerColor(1);
            mark_xz->SetMarkerSize(0.5);
            canvas_cluster->cd(2);
            mark_xz->Draw();
          }

          std::vector<sand_reco::tracker::DigitID> digits_cluster = cluster_in_container.getDigits();
          // for (uint d = 0; d < digits_cluster.size(); d++) {
          //   canvas_cluster->cd();

          //   auto digit = sand_reco::tracker::DigitCollection::getDigit(digits_cluster[d]);
          //   auto cell = sand_geo.getCellInfo(sand_geometry::tracker::CellID(digit.did));

            
          //   // Draw cells of cluster
          //   auto h = cell->second.getSize().h;
          //   auto w = cell->second.getSize().w;
          //   TVector3 r = cell->second.getWire().getDirection();
          //   TVector3 leftend = cell->second.getWire().getReadoutPoint();

          //   TVector3 AP = TVector3(digit.x, digit.y, digit.z) - leftend; 
          //   double t_prime = AP.Dot(r) / r.Mag2();
          //   t_prime = std::max(0.0, std::min(1.0, t_prime));
          //   TVector3 position_along_wire = leftend + t_prime * r;

          //   TBox* box_yz = new TBox(position_along_wire.Z() - h/2., position_along_wire.Y() - w/2., position_along_wire.Z() + h/2., position_along_wire.Y() + w/2.);
          //   TBox* box_xz = new TBox(position_along_wire.Z() - h/2., position_along_wire.X() - w/2., position_along_wire.Z() + h/2., position_along_wire.X() + w/2.);
          //   box_yz->SetFillStyle(0);
          //   box_yz->SetLineColor(color);
          //   box_yz->SetLineWidth(1);
          //   box_xz->SetFillStyle(0);
          //   box_xz->SetLineColor(color);
          //   box_xz->SetLineWidth(1);
          //   canvas_cluster->cd(1);
          //   box_yz->Draw();
          //   canvas_cluster->cd(2);
          //   box_xz->Draw();
            

          //   // Draw reco drift time of digits in cluster
          //   TEllipse* el_yz_comp = new TEllipse(position_along_wire.Z(), position_along_wire.Y(), 
          //                         sand_reco::stt::wire_radius + cell->second.getDriftVelocity() * digitId_to_drift_time[digits_cluster[d]]);
          //   TEllipse* el_xz_comp = new TEllipse(position_along_wire.Z(), position_along_wire.X(), 
          //                         sand_reco::stt::wire_radius + cell->second.getDriftVelocity() * digitId_to_drift_time[digits_cluster[d]]);
          //   el_yz_comp->SetFillStyle(0);
          //   el_yz_comp->SetLineColor(color);
          //   el_yz_comp->SetLineWidth(1);
          //   el_xz_comp->SetFillStyle(0);
          //   el_xz_comp->SetLineColor(color);
          //   el_xz_comp->SetLineWidth(1);
          //   canvas_cluster->cd(1);
          //   el_yz_comp->Draw();
          //   canvas_cluster->cd(2);
          //   el_xz_comp->Draw();
            
          //   // Draw true drift time of digits in cluster
          //   TEllipse* el_yz = new TEllipse(position_along_wire.Z(), position_along_wire.Y(), 
          //                         sand_reco::stt::wire_radius + cell->second.getDriftVelocity() * digit.drift_time);
          //   TEllipse* el_xz = new TEllipse(position_along_wire.Z(), position_along_wire.X(), 
          //                         sand_reco::stt::wire_radius + cell->second.getDriftVelocity() * digit.drift_time);
          //   el_yz->SetFillStyle(0);
          //   el_yz->SetLineWidth(1);
          //   el_yz->SetLineColor(1);
          //   el_xz->SetFillStyle(0);
          //   el_xz->SetLineWidth(1);
          //   el_xz->SetLineColor(1);
          //   canvas_cluster->cd(1);
          //   el_yz->Draw();
          //   canvas_cluster->cd(2);
          //   el_xz->Draw();

          //   // Draw hit segments for the cluster
          //   for (auto& kk:digit.hindex) {
          //     const TG4HitSegment& hseg = ev->SegmentDetectors[digit.det].at(kk);
          //     TLine* l_yz = new TLine(hseg.Start.Z(), hseg.Start.Y(), hseg.Stop.Z(), hseg.Stop.Y());
          //     TLine* l_xz = new TLine(hseg.Start.Z(), hseg.Start.X(), hseg.Stop.Z(), hseg.Stop.X());
          //     l_yz->SetLineColor(1);
          //     l_xz->SetLineColor(1);
          //     canvas_cluster->cd(1);
          //     l_yz->Draw();
          //     canvas_cluster->cd(2);
          //     l_xz->Draw();
          //   }
          // }

          for (uint d = 0; d < digits_cluster.size(); d++) {
            canvas_cluster->cd();

            auto digit = sand_reco::tracker::DigitCollection::getDigit(digits_cluster[d]);
            auto cell = sand_geo.getCellInfo(sand_geometry::tracker::CellID(digit.did));

            
            // Draw cells of cluster
            auto h = cell->second.getSize().h;
            auto w = cell->second.getSize().w;

            TBox* box_yz = new TBox(cell->second.getWire().getCenter().Z() - h/2., cell->second.getWire().getCenter().Y() - w/2., cell->second.getWire().getCenter().Z() + h/2., cell->second.getWire().getCenter().Y() + w/2.);
            TBox* box_xz = new TBox(cell->second.getWire().getCenter().Z() - h/2., cell->second.getWire().getCenter().X() - w/2., cell->second.getWire().getCenter().Z() + h/2., cell->second.getWire().getCenter().X() + w/2.);
            box_yz->SetFillStyle(0);
            box_yz->SetLineColor(color);
            box_yz->SetLineWidth(1);
            box_xz->SetFillStyle(0);
            box_xz->SetLineColor(color);
            box_xz->SetLineWidth(1);
            canvas_cluster->cd(1);
            box_yz->Draw();
            canvas_cluster->cd(2);
            box_xz->Draw();
            

            // Draw reco drift time of digits in cluster
            TEllipse* el_yz_comp = new TEllipse(cell->second.getWire().getCenter().Z(), cell->second.getWire().getCenter().Y(), 
                                  sand_reco::stt::wire_radius + cell->second.getDriftVelocity() * digitId_to_drift_time[digits_cluster[d]]);
            TEllipse* el_xz_comp = new TEllipse(cell->second.getWire().getCenter().Z(), cell->second.getWire().getCenter().X(), 
                                  sand_reco::stt::wire_radius + cell->second.getDriftVelocity() * digitId_to_drift_time[digits_cluster[d]]);
            el_yz_comp->SetFillStyle(0);
            el_yz_comp->SetLineColor(color);
            el_yz_comp->SetLineWidth(1);
            el_xz_comp->SetFillStyle(0);
            el_xz_comp->SetLineColor(color);
            el_xz_comp->SetLineWidth(1);
            canvas_cluster->cd(1);
            el_yz_comp->Draw();
            canvas_cluster->cd(2);
            el_xz_comp->Draw();
            
            // Draw true drift time of digits in cluster
            TEllipse* el_yz = new TEllipse(cell->second.getWire().getCenter().Z(), cell->second.getWire().getCenter().Y(), 
                                  sand_reco::stt::wire_radius + cell->second.getDriftVelocity() * digit.drift_time);
            TEllipse* el_xz = new TEllipse(cell->second.getWire().getCenter().Z(), cell->second.getWire().getCenter().X(), 
                                  sand_reco::stt::wire_radius + cell->second.getDriftVelocity() * digit.drift_time);
            el_yz->SetFillStyle(0);
            el_yz->SetLineWidth(1);
            el_yz->SetLineColor(1);
            el_xz->SetFillStyle(0);
            el_xz->SetLineWidth(1);
            el_xz->SetLineColor(1);
            canvas_cluster->cd(1);
            el_yz->Draw();
            canvas_cluster->cd(2);
            el_xz->Draw();

            // Draw hit segments for the cluster
            for (auto& kk:digit.hindex) {
              const TG4HitSegment& hseg = ev->SegmentDetectors[digit.det].at(kk);
              TLine* l_yz = new TLine(hseg.Start.Z(), hseg.Start.Y(), hseg.Stop.Z(), hseg.Stop.Y());
              TLine* l_xz = new TLine(hseg.Start.Z(), hseg.Start.X(), hseg.Stop.Z(), hseg.Stop.X());
              l_yz->SetLineColor(1);
              l_xz->SetLineColor(1);
              canvas_cluster->cd(1);
              l_yz->Draw();
              canvas_cluster->cd(2);
              l_xz->Draw();
            }
          }

          color++; 
          canvas_cluster->Write();
          canvas_cluster->Print("clu.pdf","pdf");
          canvas_cluster->Clear();

          canvas_cluster->Divide(2,1);
          canvas_cluster->cd(1);
          h_cluster_yz->Draw();
          canvas_cluster->cd(2);
          h_cluster_xz->Draw();
          traklet_finder.clear();

        }
      }
      canvas_cluster->Print("clu.pdf)","pdf");








    }
  }


  if (!plots) {
    innovation_test->cd();
    h_gpos_distribution->Write();
    h_gang_distribution->Write();
    h_x_diff->Write();
    h_y_diff->Write();
    h_theta_x_diff->Write();
    h_theta_y_diff->Write();
    innovation_test->Close();
  }
}
