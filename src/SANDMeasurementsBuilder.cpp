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
  t->GetEntry(std::stoi(argv[1]));

  gStyle->SetOptStat(0);
  int p[9] = {100, -2000, 2000, 100, -3200, -2300, 100, 23800, 26000};

  SANDTrackerDigitCollection::FillMap(digits);
  SANDTrackerClusterCollection clusters(&sand_geo, SANDTrackerDigitCollection::GetDigits(), SANDTrackerClusterCollection::ClusteringMethod::kCellAdjacency);
  auto digit_map =  SANDTrackerDigitCollection::GetDigits();
  
  
  TCanvas* canvas_cluster = new TCanvas("canvas_cluster","canvas_cluster",2000,1000);
  
  TH2D* h_cluster_yz = new TH2D("h","h", p[6],p[7], p[8], p[3],p[4], p[5]);
  h_cluster_yz->Draw();

  canvas_cluster->Print("clu.pdf(","pdf");

  int color = 2;
  for (const auto& container:clusters.GetContainers()) {
    for (const auto& cluster_in_container:container->GetClusters()) {
      for (auto dg:digit_map) {
        auto cell = sand_geo.get_cell_info(SANDTrackerCellID(dg.did));
        double h,w;
        cell->second.size(w,h);
        TBox* box_yz = new TBox(cell->second.wire().center().Z() - h/2., cell->second.wire().center().Y() - w/2., cell->second.wire().center().Z() + h/2., cell->second.wire().center().Y() + w/2.);
        box_yz->SetFillStyle(0);
        box_yz->SetLineColor(1);
        box_yz->SetLineWidth(1);
        canvas_cluster->cd();
        box_yz->Draw();  
      }

      if(color > 9) color = 2;
      std::vector<SANDTrackerDigitID> digits_cluster = cluster_in_container.GetDigits();
      auto d1 = SANDTrackerDigitCollection::GetDigit(digits_cluster[0]);
      auto c1 = sand_geo.get_cell_info(SANDTrackerCellID(d1.did));
      auto d2 = SANDTrackerDigitCollection::GetDigit(digits_cluster[1]);
      auto c2 = sand_geo.get_cell_info(SANDTrackerCellID(d2.did));
      auto d3 = SANDTrackerDigitCollection::GetDigit(digits_cluster[2]);
      auto c3 = sand_geo.get_cell_info(SANDTrackerCellID(d3.did));
      TLine line_yz1(c1->second.wire().center().Z(), c1->second.wire().center().Y(), c2->second.wire().center().Z(), c2->second.wire().center().Y());
      TLine line_yz2(c2->second.wire().center().Z(), c2->second.wire().center().Y(), c3->second.wire().center().Z(), c3->second.wire().center().Y());
      canvas_cluster->cd();
      line_yz1.SetLineColor(color);
      line_yz2.SetLineColor(color);
      line_yz1.Draw();
      line_yz2.Draw();
      
      for (const auto& d:digits_cluster) {
        auto digit = SANDTrackerDigitCollection::GetDigit(d);
        auto cell = sand_geo.get_cell_info(SANDTrackerCellID(digit.did));
        
        double h,w;
        cell->second.size(w,h);
        TBox* box_yz = new TBox(cell->second.wire().center().Z() - h/2., cell->second.wire().center().Y() - w/2., cell->second.wire().center().Z() + h/2., cell->second.wire().center().Y() + w/2.);
        box_yz->SetFillStyle(0);
        box_yz->SetLineColor(color);
        box_yz->SetLineWidth(1);
        box_yz->Draw();
      }
      color++;
      canvas_cluster->Print("clu.pdf","pdf");
      canvas_cluster->Clear();
      h_cluster_yz->Draw();
    }
  }

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
    canvas_cluster->cd();
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