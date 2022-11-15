#include "SANDDigitizationFLUKA.h"
#include "SANDDigitization.h"

#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "utils.h"
#include "transf.h"

using namespace sand_reco;

namespace digitization
{

namespace fluka
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

  return sand_reco::ecal::fluka::e2p2_fluka * E;
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
bool process_hit(const TG4HitSegment& hit, int& detID, int& modID, int& planeID,
                 int& cellID, double& d1, double& d2, double& t, double& de)
{
  if (debug) {
    std::cout << "ProcessHit FLUKA" << std::endl;
  }

  detID = -999;
  modID = -999;
  planeID = -999;
  cellID = -999;
  d1 = -999;
  d2 = -999;
  t = -999;

  double x = 0.5 * (hit.Start.X() + hit.Stop.X());
  double y = 0.5 * (hit.Start.Y() + hit.Stop.Y());
  double z = 0.5 * (hit.Start.Z() + hit.Stop.Z());

  t = 0.5 * (hit.Start.T() + hit.Stop.T());
  de = hit.EnergyDeposit;

  // Global to local coordinates
  //
  TLorentzVector globalPos(x, y, z, 0);
  TLorentzVector localPos = GlobalToLocalCoordinates(globalPos);
  x = localPos.X();
  y = localPos.Y();
  z = localPos.Z();
  double radius = sqrt(y * y + z * z);
  if (debug)
    std::cout << "coord locali x: " << x << "\ty: " << y << "\tz: " << z
              << "\tr: " << radius << std::endl;

  // hitAngle, cellAngle, modAngle
  //
  double modDeltaAngle = 2.0 * TMath::Pi() / 24;  // 24 modules in a ring
  double cellDeltaAngle =
      2.0 * TMath::Pi() /
      (24 * 12);  // 24 modules * 12 cells/module = number of cells in a ring
  double hitAngle = 9999;
  // This is the angle w.r.t. the y-axis. In Ideal2RealCal I used the angle
  // w.r.t. the z-axis!
  if (z != 0) {
    if (z < 0) hitAngle = 2 * atan(-z / (y + sqrt(y * y + z * z)));
    if (z > 0)
      hitAngle = 2 * atan(-z / (y + sqrt(y * y + z * z))) + 2 * TMath::Pi();
  } else if (z == 0) {
    if (y < 0) hitAngle = TMath::Pi();
    if (y > 0) hitAngle = 0;
    if (y == 0) return false;
  }
  double cellAngle =
      int(hitAngle / cellDeltaAngle) * cellDeltaAngle + cellDeltaAngle / 2;
  double modAngle =
      int((hitAngle + 0.5 * modDeltaAngle) / modDeltaAngle) * modDeltaAngle;

  // Coordinates rotation and volume finding
  TString str = "";
  double rotated_z = z * cos(-modAngle) - y * sin(-modAngle);
  double rotated_y = z * sin(-modAngle) + y * cos(-modAngle);
  if ((rotated_y > sand_reco::ecal::fluka::kloe_int_R_f) &&
      (rotated_y < sand_reco::ecal::fluka::kloe_int_R_f +
                       2 * sand_reco::ecal::fluka::ec_dzf) &&
      (abs(x) < sand_reco::ecal::barrel::lCalBarrel / 2) &&
      (abs(rotated_z) < abs(rotated_y * tan(modDeltaAngle / 2))))
    str = "volECAL";  // ECAL barrel
  else if ((rotated_y < sand_reco::ecal::fluka::ec_rf) &&
           (abs(x) > sand_reco::ecal::fluka::kloe_int_dx_f) &&
           (abs(x) < sand_reco::ecal::fluka::kloe_int_dx_f +
                         2 * sand_reco::ecal::fluka::ec_dzf))
    str = "endvolECAL";  // ECAL endcaps
  else if ((rotated_y < sand_reco::ecal::fluka::ec_rf) &&
           (abs(x) < sand_reco::ecal::fluka::kloe_int_dx_f))
    str = "tracker";  // tracker
  else
    str = "outside";  // outside

  if (debug) std::cout << "\tVol: " << str;

  // modID, planeID, cellID, d1, d2
  //
  double cellD = 0;
  if (str == "volECAL") {
    // detID
    detID = 2;
    // modID
    modID = int((hitAngle + 0.5 * modDeltaAngle) / modDeltaAngle) % 24;
    // planeID
    planeID = int((rotated_y - sand_reco::ecal::fluka::kloe_int_R_f) / 44);
    if (planeID > 4) planeID = 4;
    // cellID
    cellID = int((hitAngle + 0.5 * modDeltaAngle) / cellDeltaAngle) %
             12;  // dal punto centrale in alto in senso antiorario        // d1
                  // distance from right end (x>0)
    d1 = sand_reco::ecal::barrel::lCalBarrel / 2 - x;
    // d2 distance from left end (x<0)
    d2 = sand_reco::ecal::barrel::lCalBarrel / 2 + x;
    // cellCoord
    cellD =
        sand_reco::ecal::fluka::kloe_int_R_f + sand_reco::ecal::dzlay[0] / 2;
    for (int planeindex = 1; planeindex < planeID + 1; planeindex++)
      cellD += sand_reco::ecal::dzlay[planeindex - 1] / 2 +
               sand_reco::ecal::dzlay[planeindex] / 2;
    sand_reco::ecal::fluka::cellCoordBarrel[modID][planeID][cellID][0] = 0;
    sand_reco::ecal::fluka::cellCoordBarrel[modID][planeID][cellID][2] =
        +cellD * sin(-modAngle) -
        cellD * tan(cellAngle - modAngle) * cos(-modAngle);
    sand_reco::ecal::fluka::cellCoordBarrel[modID][planeID][cellID][1] =
        +cellD * cos(-modAngle) +
        cellD * tan(cellAngle - modAngle) * sin(-modAngle);
  } else if (str == "endvolECAL") {

    if (debug)
      std::cout << "coord ENDCAP locali x: " << x << "\ty: " << y
                << "\tz: " << z << "\tr: " << radius;
    // modID
    if (x < 0) {
      detID = 1;
      modID = 40;
    } else if (x > 0) {
      detID = 3;
      modID = 30;
    }
    // planeID
    planeID = int((abs(x) - sand_reco::ecal::fluka::kloe_int_dx_f) / 44);
    if (planeID > 4) planeID = 4;
    // cellID
    if (modID == 40)
      cellID = int((z + sand_reco::ecal::fluka::ec_rf) /
                   44);  // crescono all'aumentare di z
    else
      cellID = int((sand_reco::ecal::fluka::ec_rf - z) /
                   44);  // decrescosno all'aumentare di z
    // d1 distance from top (y>0)
    d1 = sqrt(sand_reco::ecal::fluka::ec_rf * sand_reco::ecal::fluka::ec_rf -
              z * z) -
         y;
    // d2 distance from bottom (y<0)
    d2 = sqrt(sand_reco::ecal::fluka::ec_rf * sand_reco::ecal::fluka::ec_rf -
              z * z) +
         y;
    // cellCoord
    cellD = TMath::Sign(1.0, x) * (sand_reco::ecal::fluka::kloe_int_dx_f +
                                   sand_reco::ecal::dzlay[0] / 2);
    for (int planeindex = 1; planeindex < planeID + 1; planeindex++)
      cellD +=
          TMath::Sign(1.0, x) * (sand_reco::ecal::dzlay[planeindex - 1] / 2 +
                                 sand_reco::ecal::dzlay[planeindex] / 2);
    sand_reco::ecal::fluka::cellCoordEndcap[int(modID / 10)][planeID][cellID]
                                           [0] = cellD;
    sand_reco::ecal::fluka::cellCoordEndcap[int(modID / 10)][planeID][cellID]
                                           [1] = 0;
    if (modID == 40)
      sand_reco::ecal::fluka::cellCoordEndcap[int(modID /
                                                  10)][planeID][cellID][2] =
          44 / 2 + cellID * 44 -
          sand_reco::ecal::fluka::ec_rf;  // crescono all'aumentare di cellID
    else
      sand_reco::ecal::fluka::cellCoordEndcap[int(modID /
                                                  10)][planeID][cellID][2] =
          44 / 2 - (cellID)*44 +
          sand_reco::ecal::fluka::ec_rf;  // crescono al diminuire di cellID

  } else if (str == "tracker" || str == "outside") {
    if (debug) std::cout << std::endl;
    return false;
  }
  return true;
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

        if (g == NULL && (digitization::fluka::ecal::process_hit(
                              it->second[j], detID, modID, planeID, cellID, d1,
                              d2, t0, de) == true)) {
          double en1 =
              de * sand_reco::ecal::attenuation::AttenuationFactor(d1, planeID);
          double en2 =
              de * sand_reco::ecal::attenuation::AttenuationFactor(d2, planeID);

          double ave_pe1 =
              digitization::fluka::ecal::energy_to_photo_electrons(en1);
          double ave_pe2 =
              digitization::fluka::ecal::energy_to_photo_electrons(en2);

          int pe1 = rand.Poisson(ave_pe1);
          int pe2 = rand.Poisson(ave_pe2);

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

  digitization::fluka::ecal::simulate_photo_electrons(ev, geo, photo_el, L);
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
void group_hits_by_tube(TG4Event* ev, TGeoManager* geo, int NHits,
                        Int_t DetType[10000], Float_t xPos[10000],
                        Float_t yPos[10000], Float_t zPos[10000],
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

    bool found = false;
    for (int k = 0; k < NHits; k++) {
      if (abs(x - xPos[k]) < 1 && abs(y - yPos[k]) < 1 &&
          abs(z - zPos[k]) < 1) {  // 1 mm
        if (DetType[k] == 0)
          std::cout << "ERROR: this is not a point of stt " << std::endl;
        else if (DetType[k] == 1)
          sttname = "Horizontal";
        else if (DetType[k] == 2)
          sttname = "Vertical";
        else if (DetType[k] == 3)
          sttname = "Vertical";
        else
          std::cout
              << "ERROR: this point is not in standard detector!! DetType "
              << DetType[k] << std::endl;
        found = true;
        // STT module id put to k for fluka STT digit
        // ST id put to k for fluka STT digit
        int planeid = sand_reco::stt::encodePlaneID(k, 0, DetType[k]);
        stid = sand_reco::stt::encodeSTID(planeid, k);
        break;
      }
    }

    if (found == false) {
      std::cout << "ERROR: Point not FOUND!! " << std::endl;
      exit(1);
    }

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
void digitize_stt(TG4Event* ev, TGeoManager* geo, int NHits,
                  Int_t DetType[10000], Float_t xPos[10000],
                  Float_t yPos[10000], Float_t zPos[10000],
                  std::vector<dg_tube>& digit_vec)
{
  std::map<int, std::vector<hit> > hits2Tube;
  digit_vec.clear();

  group_hits_by_tube(ev, geo, NHits, DetType, xPos, yPos, zPos, hits2Tube);
  digitization::stt::create_digits_from_hits(hits2Tube, digit_vec);
}
}  // namespace stt

// digitize event
void digitize(const char* finname, const char* foutname)
{
  TFile f(finname, "READ");

  std::cout << "This is a FLUKA SIMULATION" << std::endl;

  // MC info tree
  TTree* t = (TTree*)f.Get("EDepSimEvents");

  // Event
  TG4Event* ev = new TG4Event;
  t->SetBranchAddress("Event", &ev);

  // TGeoManager pointer used in case of geant4-based simulation
  TGeoManager* geo = 0;

  // Tree pointer and variables used in case of FLUKA simulation chain
  TTree* MapTree;

  Int_t EvtNum;
  Int_t NHits;
  Float_t xHits[100000];
  Float_t yHits[100000];
  Float_t zHits[100000];
  Int_t DetType[100000];

  // Get TGeoManager or additional Tree depending on the simulation chain

  MapTree = (TTree*)f.Get("MapTree");
  if (MapTree->GetEntries() < 1) {
    std::cout << "MapTree Empty" << std::endl;
  }

  MapTree->SetBranchAddress("EvtNum", &EvtNum);
  MapTree->SetBranchAddress("NHits", &NHits);
  MapTree->SetBranchAddress("xHits", xHits);
  MapTree->SetBranchAddress("yHits", yHits);
  MapTree->SetBranchAddress("zHits", zHits);
  MapTree->SetBranchAddress("DetType", DetType);
  t->AddFriend(MapTree);

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
    digitization::fluka::ecal::digitize_ecal(ev, geo, vec_cell);
    digitization::fluka::stt::digitize_stt(ev, geo, NHits, DetType, xHits,
                                           yHits, zHits, digit_vec);

    tout.Fill();
  }
  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;

  // write output
  fout.cd();
  tout.Write();
  fout.Close();

  f.Close();

  // cleaning
  sand_reco::stt::stL.clear();
  sand_reco::stt::stX.clear();
  sand_reco::stt::stPos.clear();
  sand_reco::stt::t0.clear();
  sand_reco::stt::tubePos.clear();
}

}  // namespace fluka

}  // namespace digitization