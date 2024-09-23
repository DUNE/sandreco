#include <TVector3.h>
#include <TTree.h>
#include <TFile.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <random>

#include "SANDGeoManager.h"
#include "SANDTrackletFinder.h"
#include "utils.h"

int main(int argc, char* argv[])
{
  TFile f(argv[2], "READ");
  TGeoManager* geo = 0;
  geo = (TGeoManager*)f.Get("EDepSimGeometry");

  SANDGeoManager sand_geo;
  sand_geo.init(geo);
  
  TFile f_d(argv[3], "READ");
  TTree* t = (TTree*)f_d.Get("tDigit");

  auto modules_map = sand_geo.get_modules_map();

  // Event
  std::vector<dg_wire>* digits = 0;
  t->SetBranchAddress("dg_wire", &digits);
  t->GetEntry(std::stoi(argv[1]));
  std::cout << digits->size() << std::endl;

  TrackletFinder tracklet_finder;
  tracklet_finder.SetSigmaPosition(0.2); // mm
  tracklet_finder.SetSigmaAngle(0.2);    // rad

  std::map<double, SANDTrackerCell> cell_digit_map;
  // int i = 0;
  // int j = 0;
  for (int i = 0; i < digits->size(); i++) {
    for (const auto& m : modules_map) {
      auto planes = m.second.planes();
      for (const auto& p:planes) {
        auto cells = p.second.getIdToCellMap();
        for (const auto& c:cells) {
          if (c.second.wire().id() == digits->at(i).did) {
            cell_digit_map.insert({i, c.second});
            std::cout << digits->at(i).did << std::endl;

            // i++;
          }
        }
      }
      // if (j == 3) break;
      // j++;
    }
  }
  std::cout << cell_digit_map.size() << std::endl;
  int p[9] = {100, -2000, 2000, 100, -4000, 0, 100, 21000, 28000};
  tracklet_finder.SetVolumeParameters(p);

  tracklet_finder.SetCells(&cell_digit_map);
  // tracklet_finder.FindTracklets();

  tracklet_finder.Draw3DWires();
  // tracklet_finder.Draw3D();
  // tracklet_finder.Draw2DDistance();
  // tracklet_finder.Draw2DWires();
  tracklet_finder.Clear();
}