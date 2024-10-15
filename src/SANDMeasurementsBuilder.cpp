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

void TryCompleteManager(TrackletMap z_to_tracklets, SParticleInfo particleInfo) {
  SANDKalmanFilterManager manager;
  manager.InitFromMC(&z_to_tracklets, particleInfo);
  manager.Run();

  auto track = manager.GetTrack();

  auto step = track.GetSteps().back();
  auto reco_state =
        step.GetStage(SANDKFTrackStep::SANDKFTrackStateStage::kSmoothing).GetStateVector();
  auto reco_mom = SANDTrackerUtils::GetMomentumInMeVFromRadiusInMM(
                                reco_state.Radius(), reco_state.TanLambda());

  std::cout << "Initial Smoothed Reco Momentum " << reco_mom << std::endl;

  return;
}

void ProcessEventWithKF(SANDGeoManager* sand_geo, TG4Event* mc_event, std::vector<dg_wire>* digits)
{
  
  int p[9] = {100, -2000, 2000, 100, -4000, -1000, 100, 23800, 26000};

  SANDTrackerDigitCollection::FillMap(digits);
  auto digit_map =  SANDTrackerDigitCollection::GetDigits();
  if (SANDTrackerDigitCollection::GetDigits().empty()) {
    return;
  }
  std::string tracker_name = SANDTrackerDigitCollection::GetDigits().begin()->det;
  SANDTrackerClusterCollection clusters(sand_geo, SANDTrackerDigitCollection::GetDigits(), SANDTrackerClusterCollection::ClusteringMethod::kCellAdjacency);
  
  TrackletFinder traklet_finder;
  traklet_finder.SetVolumeParameters(p);
  traklet_finder.SetSigmaPosition(0.2);
  traklet_finder.SetSigmaAngle(0.2);

  std::map<double, std::vector<TVectorD>> z_to_tracklets;

  SANDTrackerUtils::Init(sand_geo->GetTGeoManager());

  int gg = 0;
  for (const auto& container:clusters.GetContainers()) {
    for (const auto& cluster_in_container:container->GetClusters()) {
      // std::cout << 100 * gg / container->GetClusters().size() << std::endl;
      // gg++;
      // if (gg == 100) break;
      // if(cluster_in_container.GetZ() < 25650) continue;

      traklet_finder.SetCells(cluster_in_container);
      auto minima = traklet_finder.FindTracklets();
      double z_start = cluster_in_container.GetZ();
      for (uint trk = 0; trk < minima.size(); trk++) {
        if (minima[trk][4] < 1E-2) {
          z_to_tracklets[cluster_in_container.GetZ()].push_back(minima[trk]);
        }
      }
      traklet_finder.Clear();
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
  tree.InizializeFromEdep(*mc_event, sand_geo->GetTGeoManager());
  
  std::vector<EDEPTrajectory> primaryTrj;
  tree.Filter(std::back_insert_iterator<std::vector<EDEPTrajectory>>(primaryTrj), 
    [](const EDEPTrajectory& trj) { return trj.GetParentId() == -1;} );

  TDatabasePDG pdg_db;
  std::vector<SParticleInfo> particleInfos;
  for (auto trj:primaryTrj) {
    SParticleInfo pi;
    pi.pdg_code = trj.GetPDGCode();
    pi.id       = trj.GetId();
    auto particle = pdg_db.GetParticle(pi.pdg_code);
    if (!particle) continue;
    pi.mass = particle->Mass();
    pi.charge = particle->Charge() / 3;

    pi.pos = trj.GetTrajectoryPoints().at(string_to_component[tracker_name]).back().GetPosition().Vect();
    pi.mom = trj.GetTrajectoryPoints().at(string_to_component[tracker_name]).back().GetMomentum();
    particleInfos.push_back(pi);

    std::cout << "Initial Momentum " << trj.GetInitialMomentum().Vect().Mag() << std::endl;
  }

  int nParticles = particleInfos.size();
  
  if (nParticles == 0) {
    std::cerr << "no particles to be reconstructed...process aborted"
              << std::endl;
    return;
  }

  for (int ip = 0; ip < nParticles; ip++) {
    TryCompleteManager(z_to_tracklets, particleInfos[ip]);
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

  for (int i = 1; i < 2; i++) {
    t_h->GetEntry(i);
    t->GetEntry(i);

    ProcessEventWithKF(&sand_geo, ev, digits);


    int p[9] = {100, -2000, 2000, 100, -4000, -1000, 100, 23800, 26000};

    SANDTrackerDigitCollection::FillMap(digits);
    SANDTrackerClusterCollection clusters(&sand_geo, SANDTrackerDigitCollection::GetDigits(), SANDTrackerClusterCollection::ClusteringMethod::kCellAdjacency);
    auto digit_map =  SANDTrackerDigitCollection::GetDigits();
    
    TrackletFinder traklet_finder;
    traklet_finder.SetVolumeParameters(p);
    traklet_finder.SetSigmaPosition(0.2);
    traklet_finder.SetSigmaAngle(0.2);
    

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
    for (const auto& container:clusters.GetContainers()) {
      int gg = 0;
      for (const auto& cluster_in_container:container->GetClusters()) {
        std::cout << (double)gg / container->GetClusters().size() * 100 << std::endl;
        gg++;
        if (gg == 500) break;
        if (color > 9) color = 2;
        

        traklet_finder.SetCells(cluster_in_container);
        auto minima = traklet_finder.FindTracklets();
        // Draw tracklets
        if (minima.size() != 0) {
          canvas_cluster->cd();
          std::sort(minima.begin(), minima.end(),
                    [](TVectorD v1, TVectorD v2){ return v1[4] < v2[4];});
          double z_start = cluster_in_container.GetZ();
          for (uint trk = 0; trk < minima.size(); trk++) {
            h_minima1000->Fill(minima[trk][4]);
            h_minima_100->Fill(minima[trk][4]);
            h_minima_0_1->Fill(minima[trk][4]);
            h_minima_0_0001->Fill(minima[trk][4]);
            
            if (minima[trk][4] < 1E-2) {
              // std::cout << minima[trk][0] << " " << minima[trk][2] << std::endl;
              
              z_to_tracklets[cluster_in_container.GetZ()].push_back(minima[trk]);

              TVector2 start_tracklet_yz(z_start, minima[trk][1]);
              TVector2 start_tracklet_xz(z_start, minima[trk][0]);
              double z_end = z_start + 5 * cos(minima[trk][3]);
              double y_end = minima[trk][1] + 5 * sin(minima[trk][3]);
              double x_end = minima[trk][0] + 5 * cos(minima[trk][2]);
              TVector2 end_tracklet_yz(z_end, y_end);
              TVector2 end_tracklet_xz(z_end, x_end);
              
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
        }

        auto digitId_to_drift_time = traklet_finder.GetDigitToDriftTimeMap();
        
        // Draw cells of all digits
        for (auto digit:digit_map) {
          auto cell = sand_geo.get_cell_info(SANDTrackerCellID(digit.did));
          double h,w;
          cell->second.size(w,h);


          TVector3 r = cell->second.wire().getDirection();
          TVector3 leftend = cell->second.wire().getReadoutPoint();

          TVector3 AP = TVector3(digit.x, digit.y, digit.z) - leftend; 
          double t_prime = AP.Dot(r) / r.Mag2();
          t_prime = std::max(0.0, std::min(1.0, t_prime));
          TVector3 position_along_wire = leftend + t_prime * r;

          TBox* box_yz = new TBox(position_along_wire.Z() - h/2., position_along_wire.Y() - w/2., position_along_wire.Z() + h/2., position_along_wire.Y() + w/2.);
          TBox* box_xz = new TBox(position_along_wire.Z() - h/2., position_along_wire.X() - w/2., position_along_wire.Z() + h/2., position_along_wire.X() + w/2.);
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
        }

        std::vector<SANDTrackerDigitID> digits_cluster = cluster_in_container.GetDigits();
        for (uint d = 0; d < digits_cluster.size(); d++) {
          canvas_cluster->cd();

          auto digit = SANDTrackerDigitCollection::GetDigit(digits_cluster[d]);
          auto cell = sand_geo.get_cell_info(SANDTrackerCellID(digit.did));

          
          // Draw cells of cluster
          double h,w;
          cell->second.size(w,h);
          TVector3 r = cell->second.wire().getDirection();
          TVector3 leftend = cell->second.wire().getReadoutPoint();

          TVector3 AP = TVector3(digit.x, digit.y, digit.z) - leftend; 
          double t_prime = AP.Dot(r) / r.Mag2();
          t_prime = std::max(0.0, std::min(1.0, t_prime));
          TVector3 position_along_wire = leftend + t_prime * r;

          TBox* box_yz = new TBox(position_along_wire.Z() - h/2., position_along_wire.Y() - w/2., position_along_wire.Z() + h/2., position_along_wire.Y() + w/2.);
          TBox* box_xz = new TBox(position_along_wire.Z() - h/2., position_along_wire.X() - w/2., position_along_wire.Z() + h/2., position_along_wire.X() + w/2.);
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
          TEllipse* el_yz_comp = new TEllipse(position_along_wire.Z(), position_along_wire.Y(), 
                                sand_reco::stt::wire_radius + cell->second.driftVelocity() * digitId_to_drift_time[digits_cluster[d]]);
          TEllipse* el_xz_comp = new TEllipse(position_along_wire.Z(), position_along_wire.X(), 
                                sand_reco::stt::wire_radius + cell->second.driftVelocity() * digitId_to_drift_time[digits_cluster[d]]);
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
          TEllipse* el_yz = new TEllipse(position_along_wire.Z(), position_along_wire.Y(), 
                                sand_reco::stt::wire_radius + cell->second.driftVelocity() * digit.drift_time);
          TEllipse* el_xz = new TEllipse(position_along_wire.Z(), position_along_wire.X(), 
                                sand_reco::stt::wire_radius + cell->second.driftVelocity() * digit.drift_time);
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
          
          // h_res->Fill(digitId_to_drift_time[d] - digit.drift_time);
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
        traklet_finder.Clear();

      }
    }
    canvas_cluster->Print("clu.pdf)","pdf");
    h_minima1000->Write();
    h_minima_100->Write();
    h_minima_0_1->Write();
    h_minima_0_0001->Write();

    // int sum = 0;
    // for (auto el:z_to_tracklets) {
    //   std::cout << "At z = " << el.first << " there are " << el.second.size() << " tracklets" << std::endl;
    //   sum += el.second.size();
    // }
    // std::cout << "Total tracklets: " << sum << std::endl;




    // TCanvas* canvas_digitization = new TCanvas("canvas_digitization","canvas_digitization",2000,1000);
    // canvas_digitization->Divide(2,1);
    // TH2D* h_digitization_yz = new TH2D("h","h", p[6],p[7], p[8], p[3],p[4], p[5]);
    // TH2D* h_digitization_xz = new TH2D("h","h", p[6],p[7], p[8], p[0],p[1], p[2]);
    // canvas_digitization->cd(1);
    // h_digitization_yz->Draw();
    // canvas_digitization->cd(2);
    // h_digitization_xz->Draw();
    // for (const auto& digit:digit_map) {
    //   auto cell = sand_geo.get_cell_info(SANDTrackerCellID(digit.did));
          
    //   TEllipse* el_yz = new TEllipse(position_along_wire.Z(), position_along_wire.Y(), sand_reco::stt::wire_radius + cell->second.driftVelocity() * digit.drift_time);
    //   TEllipse* el_xz = new TEllipse(position_along_wire.Z(), position_along_wire.X(), sand_reco::stt::wire_radius + cell->second.driftVelocity() * digit.drift_time);
    //   double h,w;
    //   cell->second.size(w,h);
    //   TBox* box_yz = new TBox(position_along_wire.Z() - h/2., position_along_wire.Y() - w/2., position_along_wire.Z() + h/2., position_along_wire.Y() + w/2.);
    //   TBox* box_xz = new TBox(position_along_wire.Z() - h/2., position_along_wire.X() - w/2., position_along_wire.Z() + h/2., position_along_wire.X() + w/2.);

    //   canvas_digitization->cd(1);
    //   el_yz->SetFillStyle(0);
    //   el_yz->Draw();
      
    //   box_yz->SetFillStyle(0);
    //   box_yz->SetLineWidth(1);
    //   canvas_digitization->cd(1);
    //   box_yz->Draw();

    //   canvas_digitization->cd(2);
    //   el_xz->SetFillStyle(0);
    //   el_xz->Draw();
      
    //   box_xz->SetFillStyle(0);
    //   box_xz->SetLineWidth(1);
    //   box_xz->Draw();

    //   for (auto& hi:digit.hindex) {
    //     const TG4HitSegment& hseg = ev->SegmentDetectors[digit.det].at(hi);
    //     TLine* l_yz = new TLine(hseg.Start.Z(), hseg.Start.Y(), hseg.Stop.Z(), hseg.Stop.Y());
    //     canvas_digitization->cd(1);
    //     l_yz->Draw();
    //     TLine* l_xz = new TLine(hseg.Start.Z(), hseg.Start.X(), hseg.Stop.Z(), hseg.Stop.X());
    //     canvas_digitization->cd(2);
    //     l_xz->Draw();
    //   }

    //   TMarker* mark_yz = new TMarker(position_along_wire.Z(), position_along_wire.Y(), 5);
    //   mark_yz->SetMarkerColor(1);
    //   mark_yz->SetMarkerSize(0.5);
    //   canvas_digitization->cd(1);
    //   mark_yz->Draw();

    //   TMarker* mark_xz = new TMarker(position_along_wire.Z(), position_along_wire.X(), 5);
    //   mark_xz->SetMarkerColor(1);
    //   mark_xz->SetMarkerSize(0.5);
    //   canvas_digitization->cd(2);
    //   mark_xz->Draw();
    // }
    // canvas_digitization->SaveAs("./c2D.png");
    // canvas_digitization->SaveAs("./c2D.C");



  }
  h_res->Write();

}