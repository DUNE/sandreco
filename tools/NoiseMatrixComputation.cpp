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
#include "SANDTrackerClusterCollection.h"
#include "SANDTrackerDigitCollection.h"
#include "SANDKalmanFilter.h"
#include "utils.h"

#include "EDEPTree.h"

struct TrackletInfo {
  double x_trj, theta_x_trj, tan_x_trj, y_trj, theta_y_trj, tan_y_trj; //True measurements
  double x_trk, theta_x_trk, tan_x_trk, y_trk, theta_y_trk, tan_y_trk; //Predicted measurements
  double delta_x, delta_theta_x, delta_tan_x, delta_y, delta_theta_y, delta_tan_y; //Residuals
};


double ComputeStd(const std::vector<double>& values, double mean) {
  double squared_difference = 0.0;
  for (double v : values) {
      squared_difference += (v - mean) * (v - mean);
  }
  return std::sqrt(squared_difference / values.size());
}




void ProcessTracklets(SANDGeoManager* sand_geo, TG4Event* mc_event, std::vector<dg_wire>* digits, TH1D* h_D, TH1D* h_x, TH1D* h_y, TH1D* h_theta_x, TH1D* h_theta_y, TrackletInfo& info, TTree* tracklet_info)
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
  traklet_finder.setSigmaAngle(0.2);

  std::map<double, std::vector<TVectorD>> z_to_tracklets;

  SANDTrackerUtils::init(sand_geo->getTGeoManager());

  int gg = 0;
  for (const auto& container:clusters.getContainers()) {
    for (const auto& cluster_in_container:container->getClusters()) {
      // std::cout << 100 * gg / container->GetClusters().size() << std::endl;
      // gg++;
      // if (gg == 100) break;
      // if(cluster_in_container.GetZ() < 25650) continue;

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

  double sigma_pos = 0;
  double sigma_mom = 0;
  std::vector<double> p_true;

  for (auto trj:primaryTrj) {

    SParticleInfo pi;
    pi.pdg_code = trj.GetPDGCode();
    pi.id       = trj.GetId();
    auto particle = pdg_db.GetParticle(pi.pdg_code);
    if (!particle) continue;
    pi.mass = particle->Mass();
    pi.charge = particle->Charge() / 3;
  
    double x_smeared = rand.Gaus(trj.GetTrajectoryPoints().at(string_to_component[tracker_name]).back().GetPosition().Vect().X(), sigma_pos);
    double y_smeared = rand.Gaus(trj.GetTrajectoryPoints().at(string_to_component[tracker_name]).back().GetPosition().Vect().Y(), sigma_pos);
    double px_smeared = trj.GetTrajectoryPoints().at(string_to_component[tracker_name]).back().GetMomentum().X() * rand.Gaus(1, sigma_mom);
    double py_smeared = trj.GetTrajectoryPoints().at(string_to_component[tracker_name]).back().GetMomentum().Y() * rand.Gaus(1, sigma_mom);
    double pz_smeared = trj.GetTrajectoryPoints().at(string_to_component[tracker_name]).back().GetMomentum().Z() * rand.Gaus(1, sigma_mom);

    pi.pos = TVector3(x_smeared, y_smeared, trj.GetTrajectoryPoints().at(string_to_component[tracker_name]).back().GetPosition().Vect().Z());
    pi.mom = TVector3(px_smeared, py_smeared, pz_smeared);
    // std::cout << trj.GetTrajectoryPoints().size() << std::endl;
    // pi.pos = trj.GetTrajectoryPoints().at(string_to_component[tracker_name]).back().GetPosition().Vect();
    // pi.mom = trj.GetTrajectoryPoints().at(string_to_component[tracker_name]).back().GetMomentum();

    particleInfos.push_back(pi);
    double initial_momentum = trj.GetInitialMomentum().Vect().Mag();
    p_true.push_back(initial_momentum);
    
    std::cout << "Initial Momentum " << trj.GetInitialMomentum().Vect().Mag() << std::endl;
    std::cout << "x" << trj.GetTrajectoryPoints().at(string_to_component[tracker_name]).back().GetPosition().Vect().X() << std::endl;

    //Find the closest z coordinates of the trajectory to the traklet
    for (const auto& z : z_to_tracklets) {
      double z_trk = z.first;
      uint i_min = 0;
      double z_min = 10E8;

      for(uint i = 0; i < trj.GetTrajectoryPoints().at(string_to_component[tracker_name]).size(); i++ ){ 
        auto point = trj.GetTrajectoryPoints().at(string_to_component[tracker_name]).at(i);
        double z_trj = point.GetPosition().Vect().Z();
        double z_distance = fabs(z_trk - z_trj);

        if(z_distance < z_min){
          z_min = z_distance;
          i_min = i;
        }
      }
    
      int i_second_min = 10E8;
      if(i_min == trj.GetTrajectoryPoints().at(string_to_component[tracker_name]).size() - 1 ){
        i_second_min = i_min - 1;
      } else if (i_min == 0) {
        i_second_min = i_min + 1;
      } else {
        auto prev_point = trj.GetTrajectoryPoints().at(string_to_component[tracker_name]).at(i_min - 1);
        auto next_point = trj.GetTrajectoryPoints().at(string_to_component[tracker_name]).at(i_min + 1);
        double prev_z = prev_point.GetPosition().Vect().Z();
        double next_z = next_point.GetPosition().Vect().Z();

        double prev_z_distance = fabs(z_trk - prev_z);
        double next_z_distance = fabs(z_trk - next_z);

        if (prev_z_distance < next_z_distance) {
          i_second_min = i_min - 1;
        } else {
          i_second_min = i_min +   1;
        }
      }

      //Interpolation of z coordinate
      auto point_min        = trj.GetTrajectoryPoints().at(string_to_component[tracker_name])[i_min];
      auto point_second_min = trj.GetTrajectoryPoints().at(string_to_component[tracker_name])[i_second_min];
      
      auto first_point = point_min;
      auto second_point = point_second_min;
      if (point_min.GetPosition().Z() < point_second_min.GetPosition().Z()) {
        first_point  = point_min;
        second_point = point_second_min;
      } else {
        first_point  = point_second_min;
        second_point = point_min;
      }

      double z_diff = z_trk - first_point.GetPosition().Z();
      double alpha = z_diff / fabs(second_point.GetPosition().Z() - first_point.GetPosition().Z());

      TVector3 interpolated_pos = first_point.GetPosition().Vect() * (1 - alpha) + second_point.GetPosition().Vect() * alpha;
      TVector3 interpolated_mom = first_point.GetMomentum() * (1 - alpha) + second_point.GetMomentum() * alpha;

      //Select the best traklet based on position and direction
      auto point = trj.GetTrajectoryPoints().at(string_to_component[tracker_name]).at(i_min);
      auto p_trj =  point.GetMomentum();
      

      TVector3 p_trj_dir = interpolated_mom.Unit();
      TVector3 best_trj_point = interpolated_pos;
      TVectorD best_tracklet(z.second[0].GetNrows());

      std::vector<double> position_errors, direction_errors;
      
      double best_score = 1e8;
      for (const auto& tracklet : z.second) {
        double x_trk = tracklet[0];
        double y_trk = tracklet[1];
        double theta_xz = tracklet[2];
        double theta_yz = tracklet[3];

        //Find the closest (x,y)
        // double position_distance =  sqrt(pow(x_trk - x_trj, 2) + pow(y_trk - y_trj, 2));
        double position_distance = sqrt(pow(x_trk - best_trj_point.X(), 2) + pow(y_trk - best_trj_point.Y(), 2));
        
        //Find the best direction
        double px_trk = cos(theta_xz);
        double py_trk = sin(theta_yz);
        double pz_trk = sqrt(1 - px_trk * px_trk - py_trk * py_trk);

        TVector3 p_trk_dir(px_trk, py_trk, pz_trk);    
        
        //to be precise this is the cos of the angle between the trajectory and the traklet              
        double direction = p_trk_dir.Dot(p_trj_dir); 

        position_errors.push_back(position_distance);
        direction_errors.push_back(direction);

        double score = position_distance / 200E-3 + acos(direction) / 0.2;
        if (score < best_score) {
          best_score = score;
          std::cout << "z: " << z.first<< " position distance  "  
                    << position_distance 
                    << " angular distance  "  
                    << direction
                    << " SCORE " 
                    << best_score 
                    << std::endl;
        }
      }

      double mean_errors_position = std::accumulate(position_errors.begin(), position_errors.end(), 0.0) / position_errors.size();
      double mean_errors_direction = std::accumulate(direction_errors.begin(), direction_errors.end(), 0.0) / direction_errors.size();

      double sigma_position = ComputeStd(position_errors, mean_errors_position);
      double sigma_direction = ComputeStd(direction_errors, mean_errors_direction);

      best_score = 1e8;
      double best_direction =  -1;
      double best_position_distance = 1e8;
      for (size_t i = 0; i < position_errors.size(); i++) {
        double pos_sigmas = position_errors[i] / 200E-3;
        double ang_sigmas = acos(direction_errors[i]) / 0.2;
        // if (pos_sigmas < 3 && ang_sigmas < 3) {
          double score = pos_sigmas + ang_sigmas;
          if (score < best_score) {
              best_score = score;
              best_tracklet = z.second[i];
          }
        // }
      }

      if (best_score < 15) {
        double theta_best_xz = best_tracklet[2];
        double theta_best_yz = best_tracklet[3];   
        double px_best_trk = cos(theta_best_xz);
        double py_best_trk = sin(theta_best_yz);
        std::cout << theta_best_xz << " " << theta_best_yz << std::endl;
        double pz_best_trk = sqrt(1 - px_best_trk * px_best_trk - py_best_trk * py_best_trk);
        TVector3 p_best_trk_dir(px_best_trk, py_best_trk, pz_best_trk);    
            
        best_direction = p_best_trk_dir.Dot(p_trj_dir); 

        std::cout << " min position distance  "  
                  << sqrt(pow(best_tracklet[0] - best_trj_point.X(), 2) + pow(best_tracklet[1] - best_trj_point.Y(), 2)) 
                  << " min angular distance  "  
                  << best_direction
                  << " SCORE " 
                  << best_score 
                  << std::endl;

        //Tracklet infos
        info.x_trk = best_tracklet[0];
        info.theta_x_trk = best_tracklet[2];
        info.tan_x_trk = tan(best_tracklet[2]);
        info.y_trk = best_tracklet[1];
        info.theta_y_trk = best_tracklet[3];
        info.tan_y_trk = tan(best_tracklet[3]);
        //Trajectory infos
        info.x_trj = best_trj_point.X();
        info.theta_x_trj = atan2(p_trj_dir.Z(), p_trj_dir.X());
        info.tan_x_trj = tan(atan2(p_trj_dir.Z(), p_trj_dir.X()));
        info.y_trj = best_trj_point.Y();
        info.theta_y_trj = std::fmod(atan(p_trj_dir.Y() / p_trj_dir.Z()), M_PI);
        info.tan_y_trj = tan(std::fmod(atan(p_trj_dir.Y() / p_trj_dir.Z()), M_PI));  
        //Residuals
        info.delta_x = best_tracklet[0] - best_trj_point.X();        
        info.delta_theta_x = best_tracklet[2] - atan2(p_trj_dir.Z(), p_trj_dir.X());
        info.delta_tan_x = tan(best_tracklet[2] - atan2(p_trj_dir.Z(), p_trj_dir.X()));
        info.delta_y = best_tracklet[1] - best_trj_point.Y();    
        info.delta_theta_y = best_tracklet[3] - std::fmod(atan(p_trj_dir.Y() / p_trj_dir.Z()), M_PI);      
        info.delta_tan_y = tan(best_tracklet[3] - std::fmod(atan(p_trj_dir.Y() / p_trj_dir.Z()), M_PI));
        
        tracklet_info->Fill();
   
        h_D->Fill(best_score);
        h_x->Fill(best_tracklet[0] - best_trj_point.X());
        h_y->Fill(best_tracklet[1] - best_trj_point.Y());
        h_theta_x->Fill(best_tracklet[2] - atan2(p_trj_dir.Z(), p_trj_dir.X()));
        h_theta_y->Fill(best_tracklet[3] - std::fmod(atan(p_trj_dir.Y() / p_trj_dir.Z()), M_PI));
       
      }
    
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
    for (uint i = 0; i <  primaryTrj[ip].GetTrajectoryPoints().at(string_to_component[tracker_name]).size(); i++){
      auto point = primaryTrj[ip].GetTrajectoryPoints().at(string_to_component[tracker_name])[i];
       yz_true->SetPoint(i, point.GetPosition().Z() , point.GetPosition().Y());
    }
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

  TFile* tools_test = new TFile("tools_test.root", "RECREATE");
  // TH1D* h_p_reco = new TH1D("h_p_reco", "Reconstrcuted momentum", 500, -2, 2);
  TH1D* h_D = new TH1D("h_D", "Distribution of D ;D score;Entries", 200, 0, 200);
  TH1D* h_x = new TH1D("h_x", ";x_{trk} - x_{trj} [mm];Entries", 100, -5, 5);
  TH1D* h_y = new TH1D("h_y", ";y_{trk} - y_{trj} [mm];Entries", 100, -5, 5);
  TH1D* h_theta_x = new TH1D("h_theta_x", ";#theta^{trk}_{xz} - #theta^{trj}_{xz} [rad];Entries", 200, -3, 3);
  TH1D* h_theta_y = new TH1D("h_theta_y", ";#theta^{trk}_{yz} - #theta^{trj}_{yz} [rad];Entries", 200, -2, 2);

  TrackletInfo info;

  TTree* tracklet_info = new TTree("tracklet_info", "Tracklets info");
  tracklet_info->Branch("x_trk", &info.x_trk);
  tracklet_info->Branch("theta_x_trk", &info.theta_x_trk);
  tracklet_info->Branch("tan_x_trk", &info.tan_x_trk);
  tracklet_info->Branch("y_trk", &info.y_trk);
  tracklet_info->Branch("theta_y_trk", &info.theta_y_trk);
  tracklet_info->Branch("tan_y_trk", &info.tan_y_trk);

  tracklet_info->Branch("x_trj", &info.x_trj);
  tracklet_info->Branch("theta_x_trj", &info.theta_x_trj);
  tracklet_info->Branch("tan_x_trj", &info.tan_x_trj);
  tracklet_info->Branch("y_trj", &info.y_trj);
  tracklet_info->Branch("theta_y_trj", &info.theta_y_trj);
  tracklet_info->Branch("tan_y_trj", &info.tan_y_trj);

  tracklet_info->Branch("delta_x", &info.delta_x);
  tracklet_info->Branch("delta_theta_x", &info.delta_theta_x);
  tracklet_info->Branch("delta_tan_x", &info.delta_tan_x);
  tracklet_info->Branch("delta_y", &info.delta_y);
  tracklet_info->Branch("delta_theta_y", &info.delta_theta_y);
  tracklet_info->Branch("delta_tan_y", &info.delta_tan_y);

  SANDGeoManager sand_geo;
  sand_geo.init(geo);

  std::string geometry;
  if (geo->FindVolumeFast("STTtracker_PV")) {
    geometry = "STT";
  } else if (geo->FindVolumeFast("SANDtracker_PV")) {
    geometry = "DRIFT";
  } 
  sand_geo.fillAdjacentCells(geometry);

  int nev = t_h->GetEntries();
  for (int i = 0; i < 200; i++) {
    t_h->GetEntry(i);
    t->GetEntry(i); 

    tools_test->cd();

    ProcessTracklets(&sand_geo, ev, digits, h_D, h_x, h_y, h_theta_x, h_theta_y, info, tracklet_info);
  }
   
  h_D->Write("h_D", TObject::kOverwrite);
  h_x->Write("h_x", TObject::kOverwrite);
  h_y->Write("h_y", TObject::kOverwrite);
  h_theta_x->Write("h_theta_x", TObject::kOverwrite);
  h_theta_y->Write("h_theta_y", TObject::kOverwrite);
  tracklet_info->Write("tracklet_info", TObject::kOverwrite);
  tools_test->Close();
}


