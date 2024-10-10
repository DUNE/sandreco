#include <TVector3.h>
#include <TTree.h>
#include <TFile.h>
#include <TMarker.h>
#include <TArrow.h>

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
#include "utils.h"

int main(int argc, char* argv[])
{
  TFile f(argv[2], "READ");
  TGeoManager* geo = 0;
  geo = (TGeoManager*)f.Get("EDepSimGeometry");


  // MC info tree
  TTree* t_h = (TTree*)f.Get("EDepSimEvents");
  TG4Event* ev = new TG4Event;
  t_h->SetBranchAddress("Event", &ev);
  t_h->GetEntry(std::stoi(argv[1]));

  SANDGeoManager sand_geo;
  sand_geo.init(geo);
  
  TFile f_d(argv[3], "READ");
  TTree* t = (TTree*)f_d.Get("tDigit");

  auto& planes = sand_geo.get_planes();

  // Event
  std::vector<dg_wire>* digits = 0;
  t->SetBranchAddress("dg_wire", &digits);
  
  TFile* h_out = new TFile("h_out.root", "RECREATE");
  TH1D*  h_res = new TH1D("h_res", "h_res", 1000,-100,100);
  TH1D*  h_minima1000 = new TH1D("minima1000", "minima1000", 1000,0,100000);
  TH1D*  h_minima_100 = new TH1D("minima100", "minima100", 1000,0,100);
  TH1D*  h_minima_0_1 = new TH1D("minima0.1", "minima0.1", 1000,0,0.1);
  TH1D*  h_minima_0_0001 = new TH1D("minima0.0001", "minima0.0001", 1000,0,0.0001);

  for (int i = 1; i < 2; i++) {
    t->GetEntry(i);

    gStyle->SetOptStat(0);
    int p[9] = {100, -2000, 2000, 100, -3200, -2300, 100, 23800, 26000};

    SANDTrackerDigitCollection::FillMap(digits);
    SANDTrackerClusterCollection clusters(&sand_geo, SANDTrackerDigitCollection::GetDigits(), SANDTrackerClusterCollection::ClusteringMethod::kCellAdjacency);
    auto digit_map =  SANDTrackerDigitCollection::GetDigits();
    
    TrackletFinder traklet_finder;
    traklet_finder.SetVolumeParameters(p);
    traklet_finder.SetSigmaPosition(0.2);
    traklet_finder.SetSigmaAngle(0.2);
    

    TCanvas* canvas_cluster = new TCanvas("canvas_cluster","canvas_cluster",20000,10000);
    canvas_cluster->cd();
    
    TH2D* h_cluster_yz = new TH2D("h","h", p[6],p[7], p[8], p[3],p[4], p[5]);
    h_cluster_yz->Draw();

    canvas_cluster->Print("clu.pdf(","pdf");

    std::map<double, std::vector<TVectorD>> z_to_tracklets;

    int color = 2;
    for (const auto& container:clusters.GetContainers()) {
      int gg = 0;
      for (const auto& cluster_in_container:container->GetClusters()) {
        std::cout << (double)gg / container->GetClusters().size() * 100 << std::endl;
        gg++;
        if(color > 9) color = 2;
        

        traklet_finder.SetCells(cluster_in_container);
        auto minima = traklet_finder.FindTracklets();
        // Draw tracklets
        if (minima.size() != 0) {
          canvas_cluster->cd();
          std::sort(minima.begin(), minima.end(),
                    [](TVectorD v1, TVectorD v2){ return v1[4] < v2[4];});
          double z_start = cluster_in_container.GetZ();
          for (int trk = 0; trk < minima.size(); trk++) {
            std::cout << minima[trk][4] << std::endl;
            h_minima1000->Fill(minima[trk][4]);
            h_minima_100->Fill(minima[trk][4]);
            h_minima_0_1->Fill(minima[trk][4]);
            h_minima_0_0001->Fill(minima[trk][4]);
            
            if (minima[trk][4] < 1E-4) {
              z_to_tracklets[cluster_in_container.GetZ()].push_back(minima[trk]);

              TVector2 start_tracklet(minima[trk][1], z_start);
              double z_end = z_start + 5 * cos(minima[trk][3]);
              double y_end = minima[trk][1] + 5 * sin(minima[trk][3]);
              TVector2 end_tracklet(y_end, z_end);
              
              TLine* line_yz_tracklet = new TLine(start_tracklet.Y(), start_tracklet.X(), end_tracklet.Y(), end_tracklet.X());
              line_yz_tracklet->SetLineColor(color);
              line_yz_tracklet->SetLineWidth(1);
              line_yz_tracklet->Draw();
            }
          }
        }

        auto digitId_to_drift_time = traklet_finder.GetDigitToDriftTimeMap();
        
        // Draw cells of all digits
        // for (auto digit:digit_map) {
        //   auto cell = sand_geo.get_cell_info(SANDTrackerCellID(digit.did));
        //   double h,w;
        //   cell->second.size(w,h);
        //   TBox* box_yz = new TBox(cell->second.wire().center().Z() - h/2., cell->second.wire().center().Y() - w/2., cell->second.wire().center().Z() + h/2., cell->second.wire().center().Y() + w/2.);
        //   box_yz->SetFillStyle(0);
        //   box_yz->SetLineColor(1);
        //   box_yz->SetLineWidth(1);
        //   canvas_cluster->cd();
        //   box_yz->Draw();
        // }

        std::vector<SANDTrackerDigitID> digits_cluster = cluster_in_container.GetDigits();
        for (int d = 0; d < digits_cluster.size(); d++) {
          canvas_cluster->cd();

          auto digit = SANDTrackerDigitCollection::GetDigit(digits_cluster[d]);
          auto cell = sand_geo.get_cell_info(SANDTrackerCellID(digit.did));

          // Draw lines connecting cells in cluster
          // if (d < digits_cluster.size() - 1) {
          //   auto next_digit = SANDTrackerDigitCollection::GetDigit(digits_cluster[d+1]);
          //   auto next_cell  = sand_geo.get_cell_info(SANDTrackerCellID(next_digit.did));
          //   TLine* line_yz1 = new TLine(cell->second.wire().center().Z(), cell->second.wire().center().Y(), next_cell->second.wire().center().Z(), next_cell->second.wire().center().Y());
          //   canvas_cluster->cd();
          //   line_yz1->SetLineWidth(1);
          //   line_yz1->SetLineColor(color);
          //   line_yz1->Draw();
          // }
          
          // Draw cells of cluster
          double h,w;
          cell->second.size(w,h);
          TBox* box_yz = new TBox(cell->second.wire().center().Z() - h/2., cell->second.wire().center().Y() - w/2., cell->second.wire().center().Z() + h/2., cell->second.wire().center().Y() + w/2.);
          box_yz->SetFillStyle(0);
          box_yz->SetLineColor(color);
          box_yz->SetLineWidth(1);
          box_yz->Draw();
          

          // Draw reco drift time of digits in cluster
          TEllipse* el_yz_comp = new TEllipse(cell->second.wire().center().Z(), cell->second.wire().center().Y(), 
                                sand_reco::stt::wire_radius + cell->second.driftVelocity() * digitId_to_drift_time[digits_cluster[d]]);
          el_yz_comp->SetFillStyle(0);
          el_yz_comp->SetLineColor(color);
          el_yz_comp->SetLineWidth(1);
          el_yz_comp->Draw();
          
          // Draw true drift time of digits in cluster
          TEllipse* el_yz = new TEllipse(cell->second.wire().center().Z(), cell->second.wire().center().Y(), 
                                sand_reco::stt::wire_radius + cell->second.driftVelocity() * digit.drift_time);
          el_yz->SetFillStyle(0);
          el_yz->SetLineWidth(1);
          el_yz->SetLineColor(1);
          el_yz->Draw();

          // Draw hit segments for the cluster
          for (auto& kk:digit.hindex) {
            const TG4HitSegment& hseg = ev->SegmentDetectors[digit.det].at(kk);
            TLine* l_yz = new TLine(hseg.Start.Z(), hseg.Start.Y(), hseg.Stop.Z(), hseg.Stop.Y());
            l_yz->SetLineColor(1);
            l_yz->Draw();
          }
          
          // h_res->Fill(digitId_to_drift_time[d] - digit.drift_time);
        }
        color++;
        canvas_cluster->Write();
        canvas_cluster->Print("clu.pdf","pdf");
        canvas_cluster->Clear();
        h_cluster_yz->Draw();
        traklet_finder.Clear();

      }
    }
    canvas_cluster->Print("clu.pdf)","pdf");
    h_minima1000->Write();
    h_minima_100->Write();
    h_minima_0_1->Write();
    h_minima_0_0001->Write();

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
      auto cell = sand_geo.get_cell_info(SANDTrackerCellID(digit.did));
          
      TEllipse* el_yz = new TEllipse(cell->second.wire().center().Z(), cell->second.wire().center().Y(), sand_reco::stt::wire_radius + cell->second.driftVelocity() * digit.drift_time);
      TEllipse* el_xz = new TEllipse(cell->second.wire().center().Z(), cell->second.wire().center().X(), sand_reco::stt::wire_radius + cell->second.driftVelocity() * digit.drift_time);
      double h,w;
      cell->second.size(w,h);
      TBox* box_yz = new TBox(cell->second.wire().center().Z() - h/2., cell->second.wire().center().Y() - w/2., cell->second.wire().center().Z() + h/2., cell->second.wire().center().Y() + w/2.);
      TBox* box_xz = new TBox(cell->second.wire().center().Z() - h/2., cell->second.wire().center().X() - w/2., cell->second.wire().center().Z() + h/2., cell->second.wire().center().X() + w/2.);

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

      for (auto& i:digit.hindex) {
        const TG4HitSegment& hseg = ev->SegmentDetectors[digit.det].at(i);
        TLine* l_yz = new TLine(hseg.Start.Z(), hseg.Start.Y(), hseg.Stop.Z(), hseg.Stop.Y());
        canvas_digitization->cd(1);
        l_yz->Draw();
        TLine* l_xz = new TLine(hseg.Start.Z(), hseg.Start.X(), hseg.Stop.Z(), hseg.Stop.X());
        canvas_digitization->cd(2);
        l_xz->Draw();
      }

      TMarker* mark_yz = new TMarker(cell->second.wire().center().Z(), cell->second.wire().center().Y(), 5);
      mark_yz->SetMarkerColor(1);
      mark_yz->SetMarkerSize(0.5);
      canvas_digitization->cd(1);
      mark_yz->Draw();

      TMarker* mark_xz = new TMarker(cell->second.wire().center().Z(), cell->second.wire().center().X(), 5);
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