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

  TCanvas* canvas = new TCanvas("c2DMinimization","c2DMinimization",2000,1000);
  // canvas->Divide(2,1);

  TH2D* h_yz = new TH2D("h","h", p[6],p[7], p[8],
                               p[3],p[4], p[5]);
  TH2D* h_xz = new TH2D("h2","h2", p[6],p[7], p[8],
                               p[0],p[1], p[2]);
  canvas->cd(1);
  h_yz->Draw();
  // canvas->cd(2);
  // h_xz->Draw();

  canvas->Print("clu.pdf(","pdf");

  std::cout << "SANDTrackerDigitCollection::GetDigits().size(): " << digits->size() << std::endl;
  SANDTrackerDigitCollection::FillMap(digits);
  std::cout << "SANDTrackerDigitCollection::GetDigits().size(): " << SANDTrackerDigitCollection::GetDigits().size() << std::endl;
  SANDTrackerClusterCollection clusters(&sand_geo, SANDTrackerDigitCollection::GetDigits(), SANDTrackerClusterCollection::ClusteringMethod::kCellAdjacency);
  std::cout << "SANDTrackerDigitCollection::GetDigits().size(): " << SANDTrackerDigitCollection::GetDigits().size() << std::endl;
  auto digit_map =  SANDTrackerDigitCollection::GetDigits();
  std::cout << "SANDTrackerDigitCollection::GetDigits().size(): " << SANDTrackerDigitCollection::GetDigits().size() << std::endl;
  clusters.GetFirstDownstreamCluster();


  int color = 2;
  int sum = 0;
  for (const auto& clusters_in_plane:clusters.GetPlanes()) {

    std::cout << " size: " << clusters_in_plane->GetClusters().size() << std::endl;
    for (const auto& single_cluster_in_plane:clusters_in_plane->GetClusters()) {
      // if (sum == 100) {
      //   canvas->Print("clu.pdf)","pdf");
      //   break;
      // }
      for (auto dg:digit_map) {
        auto cell = sand_geo.get_cell_info(SANDTrackerCellID(dg.did));
        double h,w;
        cell->second.size(w,h);
        TBox* box_yz = new TBox(cell->second.wire().center().Z() - h/2., cell->second.wire().center().Y() - w/2., cell->second.wire().center().Z() + h/2., cell->second.wire().center().Y() + w/2.);
        // std::cout << cell->second.wire().center().Z() - h/2. << std::endl;
        box_yz->SetFillStyle(0);
        box_yz->SetLineColor(1);
        box_yz->SetLineWidth(1);
        canvas->cd(1);
        box_yz->Draw("same");  
        // delete box_yz;  
      }
      if(color > 9) color = 2;
      std::vector<SANDTrackerDigitID> digits_cluster = single_cluster_in_plane.GetDigits();
      std::cout << "     " << digits_cluster[0]() << std::endl;
      sum += 1;
      auto d1 = SANDTrackerDigitCollection::GetDigit(digits_cluster[0]);
      auto c1 = sand_geo.get_cell_info(SANDTrackerCellID(d1.did));
      auto d2 = SANDTrackerDigitCollection::GetDigit(digits_cluster[1]);
      auto c2 = sand_geo.get_cell_info(SANDTrackerCellID(d2.did));
      auto d3 = SANDTrackerDigitCollection::GetDigit(digits_cluster[2]);
      auto c3 = sand_geo.get_cell_info(SANDTrackerCellID(d3.did));
      std::cout << c1->second.wire().center().Z()  << " " << c1->second.wire().center().Y() << " " << c2->second.wire().center().Z() << " " << c2->second.wire().center().Y() << std::endl;
      TLine line_yz1(c1->second.wire().center().Z(), c1->second.wire().center().Y(), c2->second.wire().center().Z(), c2->second.wire().center().Y());
      TLine line_yz2(c2->second.wire().center().Z(), c2->second.wire().center().Y(), c3->second.wire().center().Z(), c3->second.wire().center().Y());
      canvas->cd(1);
      line_yz1.SetLineColor(color);
      line_yz2.SetLineColor(color);
      line_yz1.Draw("same");
      line_yz2.Draw("same");
      
      for (const auto& d:digits_cluster) {
        std::cout << "     " << d() << std::endl;
        auto digit = SANDTrackerDigitCollection::GetDigit(d);
        std::cout << "     " << digit.did << std::endl;
        
        auto cell = sand_geo.get_cell_info(SANDTrackerCellID(digit.did));
        
        std::cout << "TMP SUM: " << sum << std::endl;
        
        // TEllipse* el_yz = new TEllipse(cell->second.wire().center().Z(), cell->second.wire().center().Y(), sand_reco::stt::wire_radius + cell->second.driftVelocity() * digit.drift_time);
        // TEllipse* el_xz = new TEllipse(cell->second.wire().center().Z(), cell->second.wire().center().X(), sand_reco::stt::wire_radius + cell->second.driftVelocity() * digit.drift_time);
        double h,w;
        cell->second.size(w,h);
        TBox* box_yz = new TBox(cell->second.wire().center().Z() - h/2., cell->second.wire().center().Y() - w/2., cell->second.wire().center().Z() + h/2., cell->second.wire().center().Y() + w/2.);
        // TBox* box_xz = new TBox(cell->second.wire().center().Z() - h/2., cell->second.wire().center().X() - w/2., cell->second.wire().center().Z() + h/2., cell->second.wire().center().X() + w/2.);

        canvas->cd(1);
        // el_yz->SetFillStyle(0);
        // el_yz->SetLineColor(color);
        // el_yz->Draw();
        
        box_yz->SetFillStyle(0);
        box_yz->SetLineColor(color);
        // box_yz->SetFillColorAlpha(color, (log(digit.adc) - min) / (max - min) / 2.);
        box_yz->SetLineWidth(1);
        box_yz->Draw();

        // canvas->cd(2);
        // el_xz->SetFillStyle(0);
        // el_xz->SetLineColor(color);
        // el_xz->Draw();
        
        // box_xz->SetFillStyle(0);
        // box_xz->SetLineColor(color);
        // box_xz->SetFillColorAlpha(color, (log(digit.adc) - min) / (max - min) / 2.);
        // box_xz->SetLineWidth(1);
        // box_xz->Draw();


        // for (auto& i:digit.hindex) {
        //   const TG4HitSegment& hseg = ev->SegmentDetectors[digit.det].at(i);
        //   TLine* l_yz = new TLine(hseg.Start.Z(), hseg.Start.Y(), hseg.Stop.Z(), hseg.Stop.Y());
        //   l_yz->SetLineColor(color);
        //   canvas->cd(1);
        //   l_yz->Draw();
        //   TLine* l_xz = new TLine(hseg.Start.Z(), hseg.Start.X(), hseg.Stop.Z(), hseg.Stop.X());
        //   l_xz->SetLineColor(color);
        //   canvas->cd(2);
        //   l_xz->Draw();
        // }

        // TMarker* mark_yz = new TMarker(cell->second.wire().center().Z(), cell->second.wire().center().Y(), 5);
        // mark_yz->SetMarkerColor(color);
        // mark_yz->SetMarkerSize(0.5);
        // canvas->cd(1);
        // mark_yz->Draw();

        // TMarker* mark_xz = new TMarker(cell->second.wire().center().Z(), cell->second.wire().center().X(), 5);
        // mark_xz->SetMarkerColor(color);
        // mark_xz->SetMarkerSize(0.5);
        // canvas->cd(2);
        // mark_xz->Draw();
      }
      color++;
      canvas->Print("clu.pdf","pdf");
      canvas->Clear();
      // canvas->Divide(2,1);
      canvas->cd(1);
      h_yz->Draw();
      // canvas->cd(2);
      // h_xz->Draw();





      digits_cluster = single_cluster_in_plane.GetExtendedDigits();
      std::cout << "EXTENDED SIZE: " << digits_cluster.size() << std::endl; 
      for (const auto& d:digits_cluster) {
        auto cell = sand_geo.get_cell_info(SANDTrackerCellID(d()));
        
        double h,w;
        cell->second.size(w,h);
        TBox* box_yz = new TBox(cell->second.wire().center().Z() - h/2., cell->second.wire().center().Y() - w/2., cell->second.wire().center().Z() + h/2., cell->second.wire().center().Y() + w/2.);
        TBox* box_xz = new TBox(cell->second.wire().center().Z() - h/2., cell->second.wire().center().X() - w/2., cell->second.wire().center().Z() + h/2., cell->second.wire().center().X() + w/2.);

        canvas->cd(1);
        box_yz->SetFillStyle(0);
        box_yz->SetLineColor(1);
        box_yz->SetLineWidth(1);
        box_yz->Draw();

        canvas->cd(2);
        box_xz->SetFillStyle(0);
        box_xz->SetLineColor(1);
        box_xz->SetLineWidth(1);
        box_xz->Draw();

        TMarker* mark_yz = new TMarker(cell->second.wire().center().Z(), cell->second.wire().center().Y(), 5);
        mark_yz->SetMarkerColor(1);
        mark_yz->SetMarkerSize(0.5);
        canvas->cd(1);
        mark_yz->Draw();

        TMarker* mark_xz = new TMarker(cell->second.wire().center().Z(), cell->second.wire().center().X(), 5);
        mark_xz->SetMarkerColor(1);
        mark_xz->SetMarkerSize(0.5);
        canvas->cd(2);
        mark_xz->Draw();
      }
    }
  }
  canvas->Print("clu.pdf)","pdf");
  std::cout << sum << std::endl;

  canvas->SaveAs("./c2D.png");
  canvas->SaveAs("./c2D.C");


  TCanvas* canvas2 = new TCanvas("adjacent_cells","adjacent_cells",20000,10000);
  canvas2->Divide(2,1);

  TH2D* h_yz2 = new TH2D("h","h", p[6],p[7], p[8],
                               p[3],p[4], p[5]);
  TH2D* h_xz2 = new TH2D("h2","h2", p[6],p[7], p[8],
                               p[0],p[1], p[2]);
  canvas2->cd(1);
  h_yz2->Draw();
  canvas2->cd(2);
  h_xz2->Draw();

  int c = 0;
  for (auto& plane:planes) {
    if (c > 5) {
      break;
    }
    for(auto& cell_it:plane.getIdToCellMap()) {
      auto& cell = cell_it.second;
      double h,w;
      cell.size(w,h);
      TBox* box_yz2 = new TBox(cell.wire().center().Z() - h/2., cell.wire().center().Y() - w/2., cell.wire().center().Z() + h/2., cell.wire().center().Y() + w/2.);
      // TBox* box_xz2 = new TBox(cell.wire().center().Z() - h/2., cell.wire().center().X() - w/2., cell.wire().center().Z() + h/2., cell.wire().center().X() + w/2.);

      canvas2->cd(1);
      box_yz2->SetFillStyle(0);
      box_yz2->SetLineColor(1);
      box_yz2->SetLineWidth(1);
      box_yz2->Draw();

      // canvas2->cd(2);
      // box_xz2->SetFillStyle(0);
      // box_xz2->SetLineColor(1);
      // box_xz2->SetLineWidth(1);
      // box_xz2->Draw();

      for (auto& adj_cell:cell.getAdjacentCell()) {
        TLine* arrow_yz = new TLine(cell.wire().center().Z(), cell.wire().center().Y(), adj_cell->wire().center().Z(), adj_cell->wire().center().Y());
        // TLine* arrow_xz = new TLine(cell.wire().center().Z(), cell.wire().center().X(), adj_cell->wire().center().Z(), adj_cell->wire().center().X());
        
        canvas2->cd(1);
        arrow_yz->SetLineColor(3);
        arrow_yz->Draw();
        // canvas2->cd(2);
        // arrow_xz->SetLineColor(3);
        // arrow_xz->Draw();
      }
    }
    c++;
  }

  canvas2->SaveAs("./adj.png");
  // canvas2->SaveAs("./adj.C");

}