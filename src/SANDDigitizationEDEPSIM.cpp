#include "SANDDigitizationEDEPSIM.h"
#include "SANDDigitization.h"

#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "utils.h"

using namespace sand_reco;

namespace digitization
{

namespace edep_sim
{

namespace ecal
{
// convert deposited energy into
// into the mean number of pe
double energy_to_photo_electrons(double E)
{
  if (debug)
    std::cout << "E = " << E
              << " -> p.e. = " << sand_reco::ecal::photo_sensor::e2pe * E
              << std::endl;

  if (debug && flukatype == 1)
    std::cout << "E = " << E
              << " -> p.e. = " << sand_reco::ecal::fluka::e2p2_fluka * E
              << std::endl;

  return sand_reco::ecal::photo_sensor::e2pe * E;
}

// process calorimeter hits
// return:
//  - module id
//  - plane id
//  - cell id
//  - d1: distance from pmt1
//  - d2: distance from pmt2
//  - t: time of the hit
//  - de: energy deposit
bool process_hit(TGeoManager* g, const TG4HitSegment& hit, int& detID,
                 int& modID, int& planeID, int& cellID, double& d1, double& d2,
                 double& t, double& de)
{
  if (debug) {
    std::cout << "ProcessHit" << std::endl;
  }

  detID = -999;
  modID = -999;
  planeID = -999;
  cellID = -999;
  d1 = -999;
  d2 = -999;

  double x = 0.5 * (hit.Start.X() + hit.Stop.X());
  double y = 0.5 * (hit.Start.Y() + hit.Stop.Y());
  double z = 0.5 * (hit.Start.Z() + hit.Stop.Z());

  t = 0.5 * (hit.Start.T() + hit.Stop.T());
  de = hit.EnergyDeposit;

  /////
  TGeoNode* node = g->FindNode(x, y, z);

  if (node == 0) return false;

  TString str = node->GetName();
  TString str2 = g->GetPath();

  if (debug) {
    std::cout << "node name: " << str.Data() << std::endl;
  }

  if (sand_reco::ecal::geometry::CheckAndProcessPath(str2) == false)
    return false;
  //////

  // barrel modules
  if (sand_reco::ecal::geometry::isBarrel(str)) {

    sand_reco::ecal::geometry::BarrelModuleAndLayer(str, str2, detID, modID,
                                                    planeID);

    sand_reco::ecal::geometry::BarrelCell(x, y, z, g, node, cellID, d1, d2);

    if (debug) {
      std::cout << "hit: " << str.Data() << std::endl;
      std::cout << "\t[x,y,z]                " << x << " " << y << " " << z
                << std::endl;
      std::cout << "\t[detID,modID,planeID,cellID] " << detID << " " << modID
                << " " << planeID << " " << cellID << std::endl;
      std::cout << "\t[d1,d2,t,de]           " << d1 << " " << d2 << " " << t
                << " " << de << std::endl;
    }

    return true;
  }
  // end cap modules
  else if (sand_reco::ecal::geometry::isEndCap(str)) {

    sand_reco::ecal::geometry::EndCapModuleAndLayer(str, str2, detID, modID,
                                                    planeID);

    sand_reco::ecal::geometry::EndCapCell(x, y, z, g, node, cellID, d1, d2);

    if (debug) {
      std::cout << "hit: " << str.Data() << std::endl;
      std::cout << "\t[x,y,z]                " << x << " " << y << " " << z
                << std::endl;
      std::cout << "\t[detID,modID,planeID,cellID] " << detID << " " << modID
                << " " << planeID << " " << cellID << std::endl;
    }
    return true;
  } else {
    return false;
  }
}

void simulate_photo_electrons(TG4Event* ev, TGeoManager* g,
                              std::map<int, std::vector<pe> >& photo_el,
                              std::map<int, double>& L)
{
  int detID, modID, planeID, cellID, uniqID;
  double d1, d2, t0, de;

  for (std::map<std::string, std::vector<TG4HitSegment> >::iterator it =
           ev->SegmentDetectors.begin();
       it != ev->SegmentDetectors.end(); ++it) {
    if (it->first == "EMCalSci") {
      for (unsigned int j = 0; j < it->second.size(); j++) {

        if (g != NULL && (digitization::edep_sim::ecal::process_hit(
                              g, it->second[j], detID, modID, planeID, cellID,
                              d1, d2, t0, de) == true)) {
          double en1 =
              de * sand_reco::ecal::attenuation::AttenuationFactor(d1, planeID);
          double en2 =
              de * sand_reco::ecal::attenuation::AttenuationFactor(d2, planeID);

          double ave_pe1 =
              digitization::edep_sim::ecal::energy_to_photo_electrons(en1);
          double ave_pe2 =
              digitization::edep_sim::ecal::energy_to_photo_electrons(en2);

          int pe1 = digitization::rand.Poisson(ave_pe1);
          int pe2 = digitization::rand.Poisson(ave_pe2);

          uniqID =
              sand_reco::ecal::decoder::EncodeID(detID, modID, planeID, cellID);

          if (debug) {
            std::cout << "cell ID: " << uniqID << std::endl;
            std::cout << "\t" << de << " " << en1 << " " << en2 << std::endl;
            std::cout << "\t" << ave_pe1 << " " << ave_pe2 << std::endl;
            std::cout << "\t" << pe1 << " " << pe2 << std::endl;
          }

          // cellend 1 -> x < 0 -> ID > 0 -> left
          // cellend 2 -> x > 0 -> ID < 0 -> right

          for (int i = 0; i < pe1; i++) {
            pe this_pe;
            this_pe.time =
                digitization::ecal::photo_electron_time_to_pmt_arrival_time(t0,
                                                                            d1);
            this_pe.h_index = j;
            photo_el[uniqID].push_back(this_pe);
            L[uniqID] = d1 + d2;
          }

          for (int i = 0; i < pe2; i++) {
            pe this_pe;
            this_pe.time =
                digitization::ecal::photo_electron_time_to_pmt_arrival_time(t0,
                                                                            d2);
            this_pe.h_index = j;
            photo_el[-1 * uniqID].push_back(this_pe);
            L[-1 * uniqID] = d1 + d2;
          }
        }
      }
    }
  }
}

// simulate calorimeter responce for whole event
void digitize_ecal(TG4Event* ev, TGeoManager* geo,
                   std::vector<dg_cell>& vec_cell)
{
  std::map<int, std::vector<pe> > photo_el;
  std::map<int, std::vector<dg_ps> > ps;
  std::map<int, double> L;

  vec_cell.clear();

  if (debug) {
    std::cout << "SimulatePE" << std::endl;
  }

  digitization::edep_sim::ecal::simulate_photo_electrons(ev, geo, photo_el, L);
  if (debug) {
    std::cout << "TimeAndSignal" << std::endl;
  }
  digitization::ecal::eval_adc_and_tdc_from_photo_electrons(photo_el, ps);
  if (debug) {
    std::cout << "CollectSignal" << std::endl;
  }
  digitization::ecal::group_pmts_in_cells(geo, ps, L, vec_cell);
}

}  // namespace ecal

namespace stt
{
// Group hits into tube
void group_hits_by_tube(TG4Event* ev, TGeoManager* geo,
                        std::map<int, std::vector<hit> >& hits2Tube)
{
  hits2Tube.clear();

  for (unsigned int j = 0; j < ev->SegmentDetectors["Straw"].size(); j++) {
    const TG4HitSegment& hseg = ev->SegmentDetectors["Straw"].at(j);

    double x = 0.5 * (hseg.Start.X() + hseg.Stop.X());
    double y = 0.5 * (hseg.Start.Y() + hseg.Stop.Y());
    double z = 0.5 * (hseg.Start.Z() + hseg.Stop.Z());

    std::string sttname = "NULL";
    int stid = -999;  // should be implemented for FLUKA

    sttname = geo->FindNode(x, y, z)->GetName();

    stid = sand_reco::stt::getSTUniqID(geo, x, y, z);
    if (stid == -999) continue;

    hit h;
    h.det = sttname;
    h.did = stid;
    h.x1 = hseg.Start.X();
    h.y1 = hseg.Start.Y();
    h.z1 = hseg.Start.Z();
    h.t1 = hseg.Start.T();
    h.x2 = hseg.Stop.X();
    h.y2 = hseg.Stop.Y();
    h.z2 = hseg.Stop.Z();
    h.t2 = hseg.Stop.T();
    h.de = hseg.EnergyDeposit;
    h.pid = hseg.PrimaryId;
    h.index = j;

    hits2Tube[stid].push_back(h);
  }
}

// simulate stt responce for whole event
void digitize_stt(TG4Event* ev, TGeoManager* geo,
                  std::vector<dg_tube>& digit_vec)
{
  std::map<int, std::vector<hit> > hits2Tube;
  digit_vec.clear();

  group_hits_by_tube(ev, geo, hits2Tube);
  digitization::stt::create_digits_from_hits(hits2Tube, digit_vec);
}
}  // namespace stt

// digitize event
void digitize(const char* finname, const char* foutname)
{
  TFile f(finname, "READ");

  std::cout << "This is a standard Geant4-edepsim SIMULATION" << std::endl;

  // MC info tree
  TTree* t = (TTree*)f.Get("EDepSimEvents");

  // Event
  TG4Event* ev = new TG4Event;
  t->SetBranchAddress("Event", &ev);

  // TGeoManager pointer used in case of geant4-based simulation
  TGeoManager* geo = 0;

  // Get TGeoManager or additional Tree depending on the simulation chain
  geo = (TGeoManager*)f.Get("EDepSimGeometry");

  if (debug) std::cout << "Inizializzo la geometria" << std::endl;

  // Initialization of detector-geometry-related
  // usefull variables defined in utils.h
  // Specifically:
  // - calo cell center and size
  //   - double czlay[nLay]: "radial" position of the cells
  //   - double cxlay[nLay][nCel]: transverse position of the cells
  //   - double ec_r: ECAL endcap radius
  //   - double ec_dz: ECAL endcap thickness
  // - stt center and length
  //   - double stt_center[3]: SAND internal volume center
  //   - std::map<int, double> sand_reco::stL
  //       - KEY  : id of the tube
  //       - VALUE: length of the tubes
  //   - std::map<int, std::map<double, int>> sand_reco::stX
  //       - KEY  : id of the plane
  //       - VALUE: map with
  //         - KEY  : beam transversal coordinate of the tubes
  //         - VALUE: id of the tube
  //   - std::map<int, std::map<int, TVector2>> sand_reco::stPos
  //       - KEY  : id of the plane
  //       - VALUE: map with
  //         - KEY  : local id of the tube
  //         - VALUE: position of the center of the tube
  //                  (X in TVector2 is Z in the global system)
  //                  (Y in TVector2 is X or Y in the global system)
  //   - std::map<int, TVector2> sand_reco::tubePos
  //       - KEY  : id of the tube
  //       - VALUE: map with
  sand_reco::init(geo);

  // vector of ECAL and STT digits
  std::vector<dg_tube> digit_vec;
  std::vector<dg_cell> vec_cell;

  // output
  TFile fout(foutname, "RECREATE");
  TTree tout("tDigit", "Digitization");
  tout.Branch("dg_cell", "std::vector<dg_cell>", &vec_cell);
  tout.Branch("dg_tube", "std::vector<dg_tube>", &digit_vec);

  // number of events
  const int nev = t->GetEntries();

  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;

  // loop on all input events
  for (int i = 0; i < nev; i++) {
    t->GetEntry(i);

    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i) / nev * 100)
              << "%]" << std::flush;

    // define the T0 for this event
    // for each straw tubs:
    // std::map<int, double> sand_reco::t0
    sand_reco::stt::initT0(ev);

    // digitize ECAL and STT
    digitization::edep_sim::ecal::digitize_ecal(ev, geo, vec_cell);
    digitization::edep_sim::stt::digitize_stt(ev, geo, digit_vec);

    tout.Fill();
  }
  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;

  // write output
  fout.cd();
  tout.Write();
  geo->Write();
  fout.Close();

  f.Close();

  // cleaning
  sand_reco::stt::stL.clear();
  sand_reco::stt::stX.clear();
  sand_reco::stt::stPos.clear();
  sand_reco::stt::t0.clear();
  sand_reco::stt::tubePos.clear();
}

}  // namespace edep_sim

}  // namespace digitization