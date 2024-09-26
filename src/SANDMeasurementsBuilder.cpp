#include <TVector3.h>
#include <TTree.h>
#include <TFile.h>
#include <TMarker.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <random>

#include "SANDGeoManager.h"
#include "SANDTrackletFinder.h"
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

  auto modules_map = sand_geo.get_modules_map();

  // Event
  std::vector<dg_wire>* digits = 0;
  t->SetBranchAddress("dg_wire", &digits);
  t->GetEntry(std::stoi(argv[1]));

  // TrackletFinder tracklet_finder;
  // tracklet_finder.SetSigmaPosition(0.2); // mm
  // tracklet_finder.SetSigmaAngle(0.2);    // rad

  std::vector<std::pair<dg_wire, SANDTrackerCell>> digit_cell_pairs;
  for (int i = 0; i < (int)digits->size(); i++) {
    for (const auto& m : modules_map) {
      auto planes = m.second.planes();
      for (auto& p:planes) {
        auto cells = p.second.getIdToCellMap();
        for (const auto& c:cells) {
          if (c.second.wire().id() == digits->at(i).did) {
            digit_cell_pairs.push_back({digits->at(i),c.second});
          }
        }
      }
    }
  }


  gStyle->SetOptStat(0);
  int p[9] = {100, -2000, 2000, 100, -4000, 0, 100, 21000, 28000};

  TCanvas* c = new TCanvas("c2DMinimization","c2DMinimization",2000,1000);
  c->Divide(2,1);

  TH2D* h_yz = new TH2D("h","h", p[6],p[7], p[8],
                               p[3],p[4], p[5]);
  TH2D* h_xz = new TH2D("h","h", p[6],p[7], p[8],
                               p[0],p[1], p[2]);
  c->cd(1);
  h_yz->Draw();
  c->cd(2);
  h_xz->Draw();


  // Setting color id based on moduli id. Ugly but it works
  int ccc = 1;
  std::map<long, int> id_to_color;
  for (const auto& dc:digit_cell_pairs) {
    if (ccc > 3) {
      ccc = 1;
    }
    long module_unique_id, plane_global_id, plane_local_id, plane_type, wire_local_id;

    SANDGeoManager::decode_wire_id(dc.second.id(), plane_global_id, wire_local_id);
    SANDGeoManager::decode_plane_id(plane_global_id, module_unique_id, 
                                    plane_local_id, plane_type);
    id_to_color[module_unique_id] = ccc;
    ccc++;
  }
  ccc = 1;
  for (auto& m:id_to_color) {
    if (ccc > 4) {
      ccc = 1;
    }
    m.second = ccc;
    ccc++; 
  }

  // Used to scale the alpha value based on adc value
  std::vector<double> adcs;
  for (auto& dc:digit_cell_pairs) {
    adcs.push_back(log(dc.first.adc));
  }

  double min = *std::min_element(adcs.begin(), adcs.end());
  double max = *std::max_element(adcs.begin(), adcs.end());



  for (const auto& dc:digit_cell_pairs) {

    long module_unique_id, plane_global_id, plane_local_id, plane_type, wire_local_id;
    SANDGeoManager::decode_wire_id(dc.second.id(), plane_global_id, wire_local_id);
    SANDGeoManager::decode_plane_id(plane_global_id, module_unique_id, 
                                    plane_local_id, plane_type);


    std::cout << dc.second.wire().center().Z() << " " << dc.second.wire().center().Y() << " " << dc.second.wire().center().X() << " " << sand_reco::stt::wire_radius << " " << dc.second.driftVelocity() << " " << dc.first.drift_time << std::endl;
    // std::cout << dc.first.Z() << " " << dc.first.Y() << " " << dc.first.X() << " " << sand_reco::stt::wire_radius + dc.second.driftVelocity() * dc.first.drift_time << std::endl;
    TEllipse* el_yz = new TEllipse(dc.second.wire().center().Z(), dc.second.wire().center().Y(), sand_reco::stt::wire_radius + dc.second.driftVelocity() * dc.first.drift_time);
    TEllipse* el_xz = new TEllipse(dc.second.wire().center().Z(), dc.second.wire().center().X(), sand_reco::stt::wire_radius + dc.second.driftVelocity() * dc.first.drift_time);
    double h,w;
    dc.second.size(w,h);
    TBox* box_yz = new TBox(dc.second.wire().center().Z() - h/2., dc.second.wire().center().Y() - w/2., dc.second.wire().center().Z() + h/2., dc.second.wire().center().Y() + w/2.);
    TBox* box_xz = new TBox(dc.second.wire().center().Z() - h/2., dc.second.wire().center().X() - w/2., dc.second.wire().center().Z() + h/2., dc.second.wire().center().X() + w/2.);

    c->cd(1);
    el_yz->SetFillStyle(0);
    el_yz->SetLineColor(id_to_color[module_unique_id]);
    el_yz->Draw();
    
    box_yz->SetLineColor(id_to_color[module_unique_id]);
    box_yz->SetFillColorAlpha(id_to_color[module_unique_id], (log(dc.first.adc) - min) / (max - min) / 2.);
    box_yz->SetLineWidth(1);
    box_yz->Draw();

    c->cd(2);
    el_xz->SetFillStyle(0);
    el_xz->SetLineColor(id_to_color[module_unique_id]);
    el_xz->Draw();
    
    box_xz->SetLineColor(id_to_color[module_unique_id]);
    box_xz->SetFillColorAlpha(id_to_color[module_unique_id], (log(dc.first.adc) - min) / (max - min) / 2.);
    box_xz->SetLineWidth(1);
    box_xz->Draw();


    for (auto& h:dc.first.hindex) {
      const TG4HitSegment& hseg = ev->SegmentDetectors[dc.first.det].at(h);
      TLine* l_yz = new TLine(hseg.Start.Z(), hseg.Start.Y(), hseg.Stop.Z(), hseg.Stop.Y());
      l_yz->SetLineColor(id_to_color[module_unique_id]);
      c->cd(1);
      l_yz->Draw();
      TLine* l_xz = new TLine(hseg.Start.Z(), hseg.Start.X(), hseg.Stop.Z(), hseg.Stop.X());
      l_xz->SetLineColor(id_to_color[module_unique_id]);
      c->cd(2);
      l_xz->Draw();
    }

    TMarker* mark_yz = new TMarker(dc.second.wire().center().Z(), dc.second.wire().center().Y(), 5);
    mark_yz->SetMarkerColor(id_to_color[module_unique_id]);
    mark_yz->SetMarkerSize(0.5);
    c->cd(1);
    mark_yz->Draw();

    TMarker* mark_xz = new TMarker(dc.second.wire().center().Z(), dc.second.wire().center().Y(), 5);
    mark_xz->SetMarkerColor(id_to_color[module_unique_id]);
    mark_xz->SetMarkerSize(0.5);
    c->cd(2);
    mark_xz->Draw();
  }


  




  c->SaveAs("./c2D.C");
  c->SaveAs("./c2D.png");







































  // int p[9] = {100, -2000, 2000, 100, -4000, 0, 100, 21000, 28000};
  // tracklet_finder.SetVolumeParameters(p);

  // tracklet_finder.SetCells(&cell_digit_map);
  // // tracklet_finder.FindTracklets();

  // tracklet_finder.Draw3DWires();
  // // tracklet_finder.Draw3D();
  // // tracklet_finder.Draw2DDistance();
  // // tracklet_finder.Draw2DWires();
  // tracklet_finder.Draw2DDigits();
  // tracklet_finder.Clear();
}