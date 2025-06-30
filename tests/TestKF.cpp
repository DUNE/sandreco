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

//#include "EDEPTree.h"

void tryCompleteManager(sand_reco::kf::TrackletMap z_to_tracklets, SParticleInfo particle, TH1D* h_gpos_distribution, TH1D* h_gang_distribution, TH1D* x_res, TH1D* y_res, TH1D* theta_y_res, TH1D* theta_x_res, TH1D* mom_res, TMultiGraph* mg, TMultiGraph* mgx) {
  sand_reco::kf::Manager manager;
  manager.initFromMC(&z_to_tracklets, particle);
  manager.run();

  auto track = manager.getTrack();
  if (track.getSteps().size() > 0) { // was commented, with > 3
    auto last_step = track.getSteps().back(); //crash if empty due to the .back().
    auto reco_state =
          last_step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kSmoothing).getStateVector();
    auto reco_mom = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(
                                  reco_state.radius(), reco_state.tanLambda());

    auto initial_state = sand_reco::kf::utils::getStateVector(particle.initial_mom * 1E-3, particle.initial_pos * 1E-3, particle.charge);
    auto initial_mom = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(initial_state.radius(), initial_state.tanLambda());
 

    std::cout << "Initial Momentum " << initial_mom << std::endl;
    std::cout << "Initial Smoothed Reco Momentum " << reco_mom << std::endl;
    
    TGraph* yz_predicted = new TGraph(track.getSteps().size());
    TGraph* yz_filtered = new TGraph(track.getSteps().size());
    TGraph* yz_smoothed = new TGraph(track.getSteps().size());
    TGraph* yz_measured = new TGraph(track.getSteps().size());
    TGraph* xz_predicted = new TGraph(track.getSteps().size());
    TGraph* xz_filtered = new TGraph(track.getSteps().size());
    TGraph* xz_smoothed = new TGraph(track.getSteps().size());
    TGraph* xz_measured = new TGraph(track.getSteps().size());

    x_res->Fill(initial_state.x() - reco_state.x()); 
    y_res->Fill(initial_state.y() - reco_state.y()); 
    theta_y_res->Fill(initial_state.phi() - reco_state.phi()); 
    theta_x_res->Fill(initial_state.tanLambda() - reco_state.tanLambda()); 
    mom_res->Fill(initial_mom - reco_mom); 


    
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
  }

  return;
}

void processEventWithKF(SANDGeoManager* sand_geo, TG4Event* mc_event, std::vector<dg_wire>* digits, TH1D* h_gpos_distribution,TH1D* h_gang_distribution,
                        TH1D* h_x_diff, TH1D* h_y_diff, TH1D* h_theta_x_diff, TH1D* h_theta_y_diff, TH1D* x_res, TH1D* y_res, TH1D* theta_y_res, TH1D* theta_x_res, TH1D* mom_res)
{
  
  int p[9] = {100, -2000, 2000, 100, -4000, -0, 100, 22500, 26000};

  sand_reco::tracker::DigitCollection::fillMap(digits);
  auto digit_map =  sand_reco::tracker::DigitCollection::getDigits();
  if (sand_reco::tracker::DigitCollection::getDigits().empty()) {
    return;
  }
  std::string tracker_name = sand_reco::tracker::DigitCollection::getDigits().begin()->det;
  sand_reco::tracker::ClusterCollection clusters(sand_geo, sand_reco::tracker::DigitCollection::getDigits(), sand_reco::tracker::ClusterCollection::ClusteringMethod::kCellAdjacency);
  
  std::map<double, std::vector<TVectorD>> z_to_tracklets;

  SANDTrackerUtils::init(sand_geo->getTGeoManager());

  for (const auto& container:clusters.getContainers()) {
    for (const auto& cluster_in_container:container->getClusters()) {

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
      
      TVector3 true_pos = true_tracklet[0];
      TVector3 true_dir = true_tracklet[1];     
      double true_theta_yz = atan(true_dir.Y() / true_dir.Z());
      double true_theta_xz = atan(true_dir.X() / true_dir.Z());
      if (true_theta_xz > M_PI_2) true_theta_xz -= M_PI;

      TVectorD measurement_from_true_tracklet(4);
      measurement_from_true_tracklet(0) = true_pos.X();
      measurement_from_true_tracklet(1) = true_pos.Y();
      measurement_from_true_tracklet(2) = true_theta_xz;
      measurement_from_true_tracklet(3) = true_theta_yz;

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
  TRandom3 rand(0);
  std::map<double, std::vector<TVectorD>> z_to_best_tracklet;

  double sigma_pos = 0;
  double sigma_mom = 0;
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


    tryCompleteManager(z_to_tracklets, particleInfos[ip], h_gpos_distribution, h_gang_distribution, x_res, y_res, theta_y_res, theta_x_res, mom_res, mg, mgx);

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

  bool plots = false;
  TFile* innovation_test = new TFile("innovation_test.root", "RECREATE");
  TH1D* h_gpos_distribution = new TH1D("h_gpos_distribution", "Innovation", 100, -3, 3);
  TH1D* h_gang_distribution = new TH1D("h_gang_distribution", "Innovation", 100, -3, 3);
  TH1D* h_x_diff = new TH1D("h_x_diff", "h_x_diff", 1000, -3, 3);
  TH1D* h_y_diff = new TH1D("h_y_diff", "h_y_diff", 1000, -3, 3);
  TH1D* h_theta_x_diff = new TH1D("h_theta_x_diff", "h_theta_x_diff", 1000, -3, 3);
  TH1D* h_theta_y_diff = new TH1D("h_theta_y_diff", "h_theta_y_diff", 1000, -3, 3);
  TH1D* x_res = new TH1D("x_res", "x_res", 1000, -0.1, 0.1);
  TH1D* y_res = new TH1D("y_res", "y_res", 1000, -0.1, 0.1);
  TH1D* theta_y_res = new TH1D("theta_y_res", "theta_y_res", 1000, -1, 1);
  TH1D* theta_x_res = new TH1D("theta_x_res", "theta_x_res", 1000, -1, 1);
  TH1D* mom_res = new TH1D("mom_res", "mom_res", 1000, -1000, 1000);
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

  for (int i = 0; i < 500; i++) {
    t_h->GetEntry(i);
    t->GetEntry(i);

    if (!plots) {
      innovation_test->cd();
      processEventWithKF(&sand_geo, ev, digits, h_gpos_distribution, h_gang_distribution, h_x_diff, h_y_diff, h_theta_x_diff, h_theta_y_diff,  x_res,  y_res,  theta_y_res,  theta_x_res,  mom_res);
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
          double z_start = cluster_in_container.getZ();
          // Draw tracklets
          TVector2 start_true_tracklet_yz(z_start, true_tracklet[0].Y());
          TVector2 start_true_tracklet_xz(z_start, true_tracklet[0].X());
          double zy_end = z_start + 5 * cos(atan(true_tracklet[1].Y() / true_tracklet[1].Z()));
          double zx_end = z_start + 5 * cos(atan(true_tracklet[1].X() / true_tracklet[1].Z()));
          double y_end = true_tracklet[0].Y() + 5 * sin(atan(true_tracklet[1].Y() / true_tracklet[1].Z()));
          double x_end = true_tracklet[0].X() + 5 * sin(atan(true_tracklet[1].X() / true_tracklet[1].Z()));
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

          std::vector<sand_reco::tracker::DigitID> digits_cluster = cluster_in_container.getDigits();

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
              l_yz->SetLineColor(color);
              l_xz->SetLineColor(color);
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
    x_res->Write();
    y_res->Write();
    theta_y_res->Write();
    theta_x_res->Write();
    mom_res->Write();
    innovation_test->Close();
  }
}
