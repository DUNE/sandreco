#include <TVector3.h>
#include <TTree.h>
#include <TFile.h>
#include <TMarker.h>
#include <TArrow.h>
#include <TDatabasePDG.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <random>

#include "SANDGeoManager.h"
#include "SANDTrackletFinder.h"
#include "SANDTrackerClusterCollection.h"
#include "SANDTrackerDigitCollection.h"
#include "SANDKalmanFilter.h"
#include "utils.h"

#include "EDEPTree.h"

void tryCompleteManager(sand_reco::kf::TrackletMap z_to_tracklets, SParticleInfo particleInfo) {
  sand_reco::kf::Manager manager;
  manager.initFromMC(&z_to_tracklets, particleInfo);
  manager.run();

  auto track = manager.getTrack();
  if (track.getSteps().size() > 3) {
    std::cout << track.getSteps().size() << std::endl;
    auto step = track.getSteps().back();
    auto reco_state =
          step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kSmoothing).getStateVector();
    auto reco_mom = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(
                                  reco_state.radius(), reco_state.tanLambda());

    std::cout << "Initial Smoothed Reco Momentum " << reco_mom << std::endl;
    

    
    int i = 0;
    for (auto& step : track.getSteps()) {
      auto prediction = step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kPrediction).getStateVector();
      auto filtering = step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kFiltering).getStateVector();
      auto smoothing =  step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kSmoothing).getStateVector();
      
      
    }

  return;
}

void processEventWithKF(SANDGeoManager* sand_geo, TG4Event* mc_event, std::vector<dg_wire>* digits)
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

    particleInfos.push_back(pi);

    std::cout << "Initial Momentum " << trj.GetInitialMomentum().Vect().Mag() << std::endl;
    std::cout << "Selected Momentum " << pi.mom.Mag() << " " << pi.mom.Z() << std::endl;
  }
  int nParticles = particleInfos.size();
  
  if (nParticles == 0) {
    std::cerr << "no particles to be reconstructed...process aborted"
              << std::endl;
    return;
  }

  for (int ip = 0; ip < nParticles; ip++) {
    tryCompleteManager(z_to_tracklets, particleInfos[ip]);
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

    
  TFile* h_out = new TFile("h_out.root", "RECREATE");
  TH1D*  h_res = new TH1D("h_res", "h_res", 1000,-100,100);
  TH1D*  h_minima1000 = new TH1D("minima1000", "minima1000", 1000,0,100000);
  TH1D*  h_minima_100 = new TH1D("minima100", "minima100", 1000,0,100);
  TH1D*  h_minima_0_1 = new TH1D("minima0.1", "minima0.1", 1000,0,0.1);
  TH1D*  h_minima_0_0001 = new TH1D("minima0.0001", "minima0.0001", 1000,0,0.0001);

  SANDGeoManager sand_geo;
  sand_geo.init(geo);

  for (int i = 0; i < 20; i++) {
    t_h->GetEntry(i);
    t->GetEntry(i);

    processEventWithKF(&sand_geo, ev, digits);

    continue;
    int p[9] = {100, -2000, 2000, 100, -4000, -1000, 100, 23800, 26000};

    sand_reco::tracker::DigitCollection::fillMap(digits);
    sand_reco::tracker::ClusterCollection clusters(&sand_geo, sand_reco::tracker::DigitCollection::getDigits(), sand_reco::tracker::ClusterCollection::ClusteringMethod::kCellAdjacency);
    auto digit_map =  sand_reco::tracker::DigitCollection::getDigits();
    
    TrackletFinder traklet_finder;
    traklet_finder.setVolumeParameters(p);
    traklet_finder.setSigmaPosition(0.2);
    traklet_finder.setSigmaAngle(0.2);
    

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
        // if (gg == 500) break;
        if (color > 9) color = 2;
        

        traklet_finder.setCells(cluster_in_container);
        auto minima = traklet_finder.findTracklets();
        // // Draw tracklets
        // if (minima.size() != 0) {
        //   canvas_cluster->cd();
        //   std::sort(minima.begin(), minima.end(),
        //             [](TVectorD v1, TVectorD v2){ return v1[4] < v2[4];});
        //   double z_start = cluster_in_container.getZ();
        //   for (uint trk = 0; trk < minima.size(); trk++) {
        //     h_minima1000->Fill(minima[trk][4]);
        //     h_minima_100->Fill(minima[trk][4]);
        //     h_minima_0_1->Fill(minima[trk][4]);
        //     h_minima_0_0001->Fill(minima[trk][4]);
            
        //     if (minima[trk][4] < 1E-2) {
        //       // std::cout << minima[trk][0] << " " << minima[trk][2] << std::endl;
              
        //       z_to_tracklets[cluster_in_container.getZ()].push_back(minima[trk]);

        //       TVector2 start_tracklet_yz(z_start, minima[trk][1]);
        //       TVector2 start_tracklet_xz(z_start, minima[trk][0]);
        //       double z_end = z_start + 5 * cos(minima[trk][3]);
        //       double y_end = minima[trk][1] + 5 * sin(minima[trk][3]);
        //       double x_end = minima[trk][0] + 5 * cos(minima[trk][2]);
        //       TVector2 end_tracklet_yz(z_end, y_end);
        //       TVector2 end_tracklet_xz(z_end, x_end);
              
        //       TLine* line_yz_tracklet = new TLine(start_tracklet_yz.X(), start_tracklet_yz.Y(), end_tracklet_yz.X(), end_tracklet_yz.Y());
        //       TLine* line_xz_tracklet = new TLine(start_tracklet_xz.X(), start_tracklet_xz.Y(), end_tracklet_xz.X(), end_tracklet_xz.Y());
        //       line_yz_tracklet->SetLineColor(color);
        //       line_yz_tracklet->SetLineWidth(1);
        //       line_xz_tracklet->SetLineColor(color);
        //       line_xz_tracklet->SetLineWidth(1);
              
        //       canvas_cluster->cd(1);
        //       line_yz_tracklet->Draw();
        //       canvas_cluster->cd(2);
        //       line_xz_tracklet->Draw();
        //     }
        //   }
        // }

        bool ok = false;
        for (uint trk = 0; trk < minima.size(); trk++) {
            if (minima[trk][4] < 1E-4) {
              ok = true;
              break;
            }
        }  
        if(!ok) {
          traklet_finder.clear();
          continue;
        }

        auto digitId_to_drift_time = traklet_finder.getDigitToDriftTimeMap();
        // // Draw cells of all digits
        // for (auto digit:digit_map) {
        //   auto cell = sand_geo.getCellInfo(sand_geometry::tracker::CellID(digit.did));
        //   double h,w;
        //   cell->second.size(w,h);


        //   TVector3 r = cell->second.wire().getDirection();
        //   TVector3 leftend = cell->second.wire().getReadoutPoint();

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


        std::vector<sand_reco::tracker::DigitID> digits_cluster = cluster_in_container.getDigits();
        for (uint d = 0; d < digits_cluster.size(); d++) {
          // canvas_cluster->cd();

          auto digit = sand_reco::tracker::DigitCollection::getDigit(digits_cluster[d]);
          // auto cell = sand_geo.getCellInfo(sand_geometry::tracker::CellID(digit.did));

          
          // // Draw cells of cluster
          // double h,w;
          // cell->second.size(w,h);
          // TVector3 r = cell->second.wire().getDirection();
          // TVector3 leftend = cell->second.wire().getReadoutPoint();

          // TVector3 AP = TVector3(digit.x, digit.y, digit.z) - leftend; 
          // double t_prime = AP.Dot(r) / r.Mag2();
          // t_prime = std::max(0.0, std::min(1.0, t_prime));
          // TVector3 position_along_wire = leftend + t_prime * r;

          // TBox* box_yz = new TBox(position_along_wire.Z() - h/2., position_along_wire.Y() - w/2., position_along_wire.Z() + h/2., position_along_wire.Y() + w/2.);
          // TBox* box_xz = new TBox(position_along_wire.Z() - h/2., position_along_wire.X() - w/2., position_along_wire.Z() + h/2., position_along_wire.X() + w/2.);
          // box_yz->SetFillStyle(0);
          // box_yz->SetLineColor(color);
          // box_yz->SetLineWidth(1);
          // box_xz->SetFillStyle(0);
          // box_xz->SetLineColor(color);
          // box_xz->SetLineWidth(1);
          // canvas_cluster->cd(1);
          // box_yz->Draw();
          // canvas_cluster->cd(2);
          // box_xz->Draw();
          

          // // Draw reco drift time of digits in cluster
          // TEllipse* el_yz_comp = new TEllipse(position_along_wire.Z(), position_along_wire.Y(), 
          //                       sand_reco::stt::wire_radius + cell->second.driftVelocity() * digitId_to_drift_time[digits_cluster[d]]);
          // TEllipse* el_xz_comp = new TEllipse(position_along_wire.Z(), position_along_wire.X(), 
          //                       sand_reco::stt::wire_radius + cell->second.driftVelocity() * digitId_to_drift_time[digits_cluster[d]]);
          // el_yz_comp->SetFillStyle(0);
          // el_yz_comp->SetLineColor(color);
          // el_yz_comp->SetLineWidth(1);
          // el_xz_comp->SetFillStyle(0);
          // el_xz_comp->SetLineColor(color);
          // el_xz_comp->SetLineWidth(1);
          // canvas_cluster->cd(1);
          // el_yz_comp->Draw();
          // canvas_cluster->cd(2);
          // el_xz_comp->Draw();
          
          // // Draw true drift time of digits in cluster
          // TEllipse* el_yz = new TEllipse(position_along_wire.Z(), position_along_wire.Y(), 
          //                       sand_reco::stt::wire_radius + cell->second.driftVelocity() * digit.drift_time);
          // TEllipse* el_xz = new TEllipse(position_along_wire.Z(), position_along_wire.X(), 
          //                       sand_reco::stt::wire_radius + cell->second.driftVelocity() * digit.drift_time);
          // el_yz->SetFillStyle(0);
          // el_yz->SetLineWidth(1);
          // el_yz->SetLineColor(1);
          // el_xz->SetFillStyle(0);
          // el_xz->SetLineWidth(1);
          // el_xz->SetLineColor(1);
          // canvas_cluster->cd(1);
          // el_yz->Draw();
          // canvas_cluster->cd(2);
          // el_xz->Draw();

          // // Draw hit segments for the cluster
          // for (auto& kk:digit.hindex) {
          //   const TG4HitSegment& hseg = ev->SegmentDetectors[digit.det].at(kk);
          //   TLine* l_yz = new TLine(hseg.Start.Z(), hseg.Start.Y(), hseg.Stop.Z(), hseg.Stop.Y());
          //   TLine* l_xz = new TLine(hseg.Start.Z(), hseg.Start.X(), hseg.Stop.Z(), hseg.Stop.X());
          //   l_yz->SetLineColor(1);
          //   l_xz->SetLineColor(1);
          //   canvas_cluster->cd(1);
          //   l_yz->Draw();
          //   canvas_cluster->cd(2);
          //   l_xz->Draw();
          // }
          
          // std::cout << digitId_to_drift_time[digits_cluster[d]] << " " << digit.drift_time << std::endl;
          h_res->Fill(digitId_to_drift_time[digits_cluster[d]] - digit.drift_time);
        }
        color++; 
        // canvas_cluster->Write();
        // canvas_cluster->Print("clu.pdf","pdf");
        // canvas_cluster->Clear();

        // canvas_cluster->Divide(2,1);
        // canvas_cluster->cd(1);
        // h_cluster_yz->Draw();
        // canvas_cluster->cd(2);
        // h_cluster_xz->Draw();
        traklet_finder.clear();

      }
    }
    // canvas_cluster->Print("clu.pdf)","pdf");
    // h_minima1000->Write();
    // h_minima_100->Write();
    // h_minima_0_1->Write();
    // h_minima_0_0001->Write();

    int sum = 0;
    for (auto el:z_to_tracklets) {
      std::cout << "At z = " << el.first << " there are " << el.second.size() << " tracklets" << std::endl;
      sum += el.second.size();
    }
    std::cout << "Total tracklets: " << sum << std::endl;




    TCanvas* canvas_digitization = new TCanvas("canvas_digitization","canvas_digitization",2000,1000);
    canvas_digitization->Divide(2,1);
    TH2D* h_digitization_yz = new TH2D("h","h", p[6],p[7], p[8], p[3],p[4], p[5]);
    TH2D* h_digitization_xz = new TH2D("h","h", p[6],p[7], p[8], p[0],p[1], p[2]);
    canvas_digitization->cd(1);
    h_digitization_yz->Draw();
    canvas_digitization->cd(2);
    h_digitization_xz->Draw();
    for (const auto& digit:digit_map) {
      auto cell = sand_geo.getCellInfo(sand_geometry::tracker::CellID(digit.did));
      

      // Draw cells of cluster
      auto cell_size = cell->second.getSize();
      double h = cell_size.h;
      double w = cell_size.w;
      TVector3 r = cell->second.getWire().getDirection();
      TVector3 leftend = cell->second.getWire().getReadoutPoint();

      TVector3 AP = TVector3(digit.x, digit.y, digit.z) - leftend; 
      double t_prime = AP.Dot(r) / r.Mag2();
      t_prime = std::max(0.0, std::min(1.0, t_prime));
      TVector3 position_along_wire = leftend + t_prime * r;


      TEllipse* el_yz = new TEllipse(position_along_wire.Z(), position_along_wire.Y(), sand_reco::stt::wire_radius + cell->second.getDriftVelocity() * digit.drift_time);
      TEllipse* el_xz = new TEllipse(position_along_wire.Z(), position_along_wire.X(), sand_reco::stt::wire_radius + cell->second.getDriftVelocity() * digit.drift_time);
      TBox* box_yz = new TBox(position_along_wire.Z() - h/2., position_along_wire.Y() - w/2., position_along_wire.Z() + h/2., position_along_wire.Y() + w/2.);
      TBox* box_xz = new TBox(position_along_wire.Z() - h/2., position_along_wire.X() - w/2., position_along_wire.Z() + h/2., position_along_wire.X() + w/2.);

      canvas_digitization->cd(1);
      el_yz->SetFillStyle(0);
      el_yz->Draw();
      
      box_yz->SetFillStyle(0);
      box_yz->SetLineWidth(1);
      canvas_digitization->cd(1);
      box_yz->Draw();

      canvas_digitization->cd(2);
      el_xz->SetFillStyle(0);
      el_xz->Draw();
      
      box_xz->SetFillStyle(0);
      box_xz->SetLineWidth(1);
      box_xz->Draw();

      for (auto& hi:digit.hindex) {
        const TG4HitSegment& hseg = ev->SegmentDetectors[digit.det].at(hi);
        TLine* l_yz = new TLine(hseg.Start.Z(), hseg.Start.Y(), hseg.Stop.Z(), hseg.Stop.Y());
        canvas_digitization->cd(1);
        l_yz->Draw();
        TLine* l_xz = new TLine(hseg.Start.Z(), hseg.Start.X(), hseg.Stop.Z(), hseg.Stop.X());
        canvas_digitization->cd(2);
        l_xz->Draw();
      }

      TMarker* mark_yz = new TMarker(position_along_wire.Z(), position_along_wire.Y(), 5);
      mark_yz->SetMarkerColor(1);
      mark_yz->SetMarkerSize(0.5);
      canvas_digitization->cd(1);
      mark_yz->Draw();

      TMarker* mark_xz = new TMarker(position_along_wire.Z(), position_along_wire.X(), 5);
      mark_xz->SetMarkerColor(1);
      mark_xz->SetMarkerSize(0.5);
      canvas_digitization->cd(2);
      mark_xz->Draw();
    }
    canvas_digitization->SaveAs("./c2D.png");
    canvas_digitization->SaveAs("./c2D.C");



  }
  h_res->Write();

}