#include <TGeoManager.h>
#include <TString.h>
#include <TGeoNode.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TGeoTrd2.h>
#include <TGeoTube.h>
#include <TChain.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMath.h>

#include "TG4Event.h"
#include "TG4HitSegment.h"

#include <vector>
#include <map>
#include <iostream>
#include <iomanip>

#include "struct.h"
#include "utils.h"
#include "transf.h"

// Energy MeV
// Distance mm
// Time ns

using namespace kloe_simu;

TRandom3 r(0);

/*
// get fiber attenuation factor.
// It depends on distance from pmt (d)
// and planeID (planes have different fibers)
double Attenuation(double d, int planeID)
{
  double atl2 = 0.0;

  switch (planeID) {
    case 0:
    case 1:
      atl2 = atl2_01;
      break;

    case 2:
      atl2 = atl2_2;
      break;

    case 3:
    case 4:
      atl2 = atl2_34;
      break;

    default:
      atl2 = -999.0;
      break;
  }

  if (debug) {
    std::cout << "planeID = " << planeID << std::endl;
    std::cout << "\tp1   = " << p1 << std::endl;
    std::cout << "\tatl1 = " << atl1 << std::endl;
    std::cout << "\tatl2 = " << atl2 << std::endl;
    std::cout << "\tatt  = "
              << p1* TMath::Exp(-d / atl1) + (1. - p1) * TMath::Exp(-d / atl2)
              << std::endl;
  }

  return p1 * TMath::Exp(-d / atl1) + (1. - p1) * TMath::Exp(-d / atl2);
}
*/
// convert deposited energy into
// into the mean number of pe
double E2PE(double E)
{
  if (debug) std::cout << "E = " << E << " -> p.e. = " << e2p2* E << std::endl;
 
  if (debug && flukatype==1) std::cout << "E = " << E << " -> p.e. = " << e2p2_fluka * E << std::endl;
  
  if(flukatype==1) return e2p2_fluka * E;
  else return e2p2 * E;
}

// simulate pe arrival time to pmt
double petime(double t0, double d)
{
  /*
     - For each photoelectron: Time for TDC simulation obtained from

C  PHOTOELECTRON TIME :  Particle TIME in the cell
C                      + SCINTILLATION DECAY TIME +
C                      + signal propagation to the cell
C                      + 1ns  uncertainty

               TPHE = Part_time+TSDEC+DPM1*VLFB+Gauss(1ns)

      VLFB = 5.85 ns/m
!!!! Input-TDC Scintillation time -
               TSDEC = TSCIN*(1./RNDMPH(1)-1)**TSCEX  (ns)

      TSCIN  3.08  ns
      TSCEX  0.588
  */

  double tdec = tscin * TMath::Power(1. / r.Uniform() - 1., tscex);

  double time = t0 + tdec + vlfb * d * mm_to_m + r.Gaus();

  if (debug) {
    std::cout << "time : " << time << std::endl;
    std::cout << "t0   : " << t0 << std::endl;
    std::cout << "scint: " << tdec << std::endl;
    std::cout << "prop : " << vlfb* d* mm_to_m << std::endl;
  }

  return time;
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
bool ProcessHit(TGeoManager* g, const TG4HitSegment& hit, int& modID,
                int& planeID, int& cellID, double& d1, double& d2, double& t,
                double& de)
{
  if (debug) {
    std::cout << "ProcessHit" << std::endl;
  }

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

  if (CheckAndProcessPath(str2) == false) return false;
  //////

  // barrel modules
  if (isBarrel(str)) {

    BarrelModuleAndLayer(str, str2, modID, planeID);

    BarrelCell(x, y, z, g, node, cellID, d1, d2);

    if (debug) {
      std::cout << "hit: " << str.Data() << std::endl;
      std::cout << "\t[x,y,z]                " << x << " " << y << " " << z
                << std::endl;
      std::cout << "\t[modID,planeID,cellID] " << modID << " " << planeID << " "
                << cellID << std::endl;
      std::cout << "\t[d1,d2,t,de]           " << d1 << " " << d2 << " " << t
                << " " << de << std::endl;
    }

    return true;
  }
  // end cap modules
  else if (isEndCap(str)) {

    if (debug) {
      TLorentzVector gPos(x, y, z, 0);
      TLorentzVector lPos = GlobalToLocalCoordinates(gPos);

      std::cout << "coord locali " << lPos.X() << " " << lPos.Y() << " "
                << lPos.Z() << std::endl;
      std::cout << "coord globali" << gPos.X() << " " << gPos.Y() << " "
                << gPos.Z() << std::endl;
    }

    EndCapModuleAndLayer(str, str2, modID, planeID);

    EndCapCell(x, y, z, g, node, cellID, d1, d2);

    if (debug) {
      std::cout << "hit: " << str.Data() << std::endl;
      std::cout << "\t[x,y,z]                " << x << " " << y << " " << z
                << std::endl;
      std::cout << "\t[modID,planeID,cellID] " << modID << " " << planeID << " "
                << cellID << std::endl;
    }
    return true;
  } else {
    return false;
  }
}

bool ProcessHitFluka(const TG4HitSegment& hit, int& modID, int& planeID,
                     int& cellID, double& d1, double& d2, double& t, double& de)
{
  if (debug) {
    std::cout << "ProcessHit FLUKA" << std::endl;
  }

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
  if ((rotated_y > kloe_int_R_f) && (rotated_y < kloe_int_R_f + 2 * ec_dzf) &&
      (abs(x) < lCalBarrel / 2) &&
      (abs(rotated_z) < abs(rotated_y * tan(modDeltaAngle / 2))))
    str = "volECAL";  // ECAL barrel
  else if ((rotated_y < ec_rf) && (abs(x) > kloe_int_dx_f) &&
           (abs(x) < kloe_int_dx_f + 2 * ec_dzf))
    str = "endvolECAL";  // ECAL endcaps
  else if ((rotated_y < ec_rf) && (abs(x) < kloe_int_dx_f))
    str = "tracker";  // tracker
  else
    str = "outside";  // outside

  if (debug) std::cout << "\tVol: " << str;

  // modID, planeID, cellID, d1, d2
  //
  double cellD = 0;
  if (str == "volECAL") {
    // modID
    modID = int((hitAngle + 0.5 * modDeltaAngle) / modDeltaAngle) % 24;
    // planeID
    planeID = int((rotated_y - kloe_int_R_f) / 44);
    if (planeID > 4) planeID = 4;
    // cellID
    cellID = int((hitAngle + 0.5 * modDeltaAngle) / cellDeltaAngle) %
             12;  // dal punto centrale in alto in senso antiorario        // d1
                  // distance from right end (x>0)
    d1 = lCalBarrel / 2 - x;
    // d2 distance from left end (x<0)
    d2 = lCalBarrel / 2 + x;
    // cellCoord
    cellD = kloe_int_R_f + dzlay[0] / 2;
    for (int planeindex = 1; planeindex < planeID + 1; planeindex++)
      cellD += dzlay[planeindex - 1] / 2 + dzlay[planeindex] / 2;
    cellCoordBarrel[modID][planeID][cellID][0] = 0;
    cellCoordBarrel[modID][planeID][cellID][2] =
        +cellD * sin(-modAngle) -
        cellD * tan(cellAngle - modAngle) * cos(-modAngle);
    cellCoordBarrel[modID][planeID][cellID][1] =
        +cellD * cos(-modAngle) +
        cellD * tan(cellAngle - modAngle) * sin(-modAngle);
  } else if (str == "endvolECAL") {

    if (debug)
      std::cout << "coord ENDCAP locali x: " << x << "\ty: " << y
                << "\tz: " << z << "\tr: " << radius;
    // modID
    if (x < 0)
      modID = 40;
    else if (x > 0)
      modID = 30;
    // planeID
    planeID = int((abs(x) - kloe_int_dx_f) / 44);
    if (planeID > 4) planeID = 4;
    // cellID
    if (modID == 40)
      cellID = int((z + ec_rf) / 44);  // crescono all'aumentare di z
    else
      cellID = int((ec_rf - z) / 44);  // decrescosno all'aumentare di z
    // d1 distance from top (y>0)
    d1 = sqrt(ec_rf * ec_rf - z * z) - y;
    // d2 distance from bottom (y<0)
    d2 = sqrt(ec_rf * ec_rf - z * z) + y;
    // cellCoord
    cellD = TMath::Sign(1.0, x) * (kloe_int_dx_f + dzlay[0] / 2);
    for (int planeindex = 1; planeindex < planeID + 1; planeindex++)
      cellD += TMath::Sign(1.0, x) *
               (dzlay[planeindex - 1] / 2 + dzlay[planeindex] / 2);
    cellCoordEndcap[int(modID / 10)][planeID][cellID][0] = cellD;
    cellCoordEndcap[int(modID / 10)][planeID][cellID][1] = 0;
    if (modID == 40)
      cellCoordEndcap[int(modID / 10)][planeID][cellID][2] =
          44 / 2 + cellID * 44 - ec_rf;  // crescono all'aumentare di cellID
    else
      cellCoordEndcap[int(modID / 10)][planeID][cellID][2] =
          44 / 2 - (cellID) * 44 + ec_rf;  // crescono al diminuire di cellID

  } else if (str == "tracker" || str == "outside") {
    if (debug) std::cout << std::endl;
    return false;
  }
  return true;
}

void SimulatePE(TG4Event* ev, TGeoManager* g,
                std::map<int, std::vector<pe> >& photo_el,
                std::map<int, double>& L)
{
  int modID, planeID, cellID, id;
  double d1, d2, t0, de;

  for (std::map<std::string, std::vector<TG4HitSegment> >::iterator it =
           ev->SegmentDetectors.begin();
       it != ev->SegmentDetectors.end(); ++it) {
    if (it->first == "EMCalSci") {
      for (unsigned int j = 0; j < it->second.size(); j++) {

        if ((g != NULL && (ProcessHit(g, it->second[j], modID, planeID, cellID,
                                      d1, d2, t0, de) == true)) ||
            (g == NULL && (ProcessHitFluka(it->second[j], modID, planeID,
                                           cellID, d1, d2, t0, de) == true))) {
          double en1 = de * AttenuationFactor(d1, planeID);
          double en2 = de * AttenuationFactor(d2, planeID);

          double ave_pe1 = E2PE(en1);
          double ave_pe2 = E2PE(en2);

          int pe1 = r.Poisson(ave_pe1);
          int pe2 = r.Poisson(ave_pe2);

          id = EncodeID(modID, planeID, cellID);

          if (debug) {
            std::cout << "cell ID: " << id << std::endl;
            std::cout << "\t" << de << " " << en1 << " " << en2 << std::endl;
            std::cout << "\t" << ave_pe1 << " " << ave_pe2 << std::endl;
            std::cout << "\t" << pe1 << " " << pe2 << std::endl;
          }

          // cellend 1 -> x < 0 -> ID > 0 -> left
          // cellend 2 -> x > 0 -> ID < 0 -> right

          for (int i = 0; i < pe1; i++) {
            pe this_pe;
            this_pe.time = petime(t0, d1);
            this_pe.h_index = j;
            photo_el[id].push_back(this_pe);
            L[id] = d1 + d2;
          }

          for (int i = 0; i < pe2; i++) {
            pe this_pe;
            this_pe.time = petime(t0, d2);
            this_pe.h_index = j;
            photo_el[-1 * id].push_back(this_pe);
            L[-1 * id] = d1 + d2;
          }
        }
      }
    }
  }
}

// from simulated pe produce adc e tdc of calo cell
void TimeAndSignal(std::map<int, std::vector<pe> >& photo_el,
                   std::map<int, std::vector<dg_ps> >& map_pmt)
{
  /*
    -  ADC - Proportional to NPHE
    -  TDC - Constant fraction - simulated
             TPHE(1...NPHE) in increasing time order
             IND_SEL= 0.15*NPHE
             TDC_cell = TPHE(IND_SEL)
  */

  // https://www-sciencedirect-com.ezproxy.cern.ch/science/article/pii/S0168900297013491

  double int_start;
  int pe_count;
  int start_index;
  int index;

  std::vector<pe> photo_el_digit;

  for (std::map<int, std::vector<pe> >::iterator it = photo_el.begin();
       it != photo_el.end(); ++it) {
    // order by arrival time
    std::sort(it->second.begin(), it->second.end(), kloe_simu::isPeBefore);

    auto side = 2 * (it->first > 0) - 1;

    photo_el_digit.clear();

    int_start = it->second.front().time;
    pe_count = 0;
    start_index = 0;
    index = 0;

    for (std::vector<pe>::iterator this_pe = it->second.begin();
         this_pe != it->second.end(); ++this_pe) {
      // integrate for int_time
      if (this_pe->time < int_start + int_time) {
        pe_count++;
        photo_el_digit.push_back(*this_pe);
      } else if (this_pe->time > int_start + int_time + dead_time) {
        // above threshold -> digit
        if (pe_count > pe_threshold) {
          dg_ps signal;
          signal.side = side;
          signal.adc = pe2ADC * pe_count;
          index = int(costant_fraction * pe_count) + start_index;
          signal.tdc = it->second[index].time;
          signal.photo_el = photo_el_digit;
          map_pmt[it->first].push_back(signal);
        }
        // get ready for next digiit
        pe_count = 1;
        photo_el_digit.clear();
        int_start = this_pe->time;
        start_index = this_pe - it->second.begin();
      }
    }

    if (pe_count > pe_threshold) {
      dg_ps signal;
      signal.side = side;
      signal.adc = pe2ADC * pe_count;
      index = int(costant_fraction * pe_count) + start_index;
      signal.tdc = it->second[index].time;
      signal.photo_el = photo_el_digit;
      map_pmt[it->first].push_back(signal);
    }
  }
}

// construct calo digit and collect them in a vector
void CollectSignal(TGeoManager* geo,
                   std::map<int, std::vector<dg_ps> >& ps,
                   std::map<int, double>& L,
                   std::vector<dg_cell>& vec_cell)
{
  std::map<int, dg_cell> map_cell;
  dg_cell* c;

  for (std::map<int, std::vector<dg_ps> >::iterator it = ps.begin(); it != ps.end();
       ++it) {
    int id = abs(it->first);

    c = &(map_cell[id]);

    c->id = id;
    DecodeID(c->id, c->mod, c->lay, c->cel);
    c->l = L[it->first];

    if (it->first >= 0) {
      c->ps1 = it->second;
    } else {
      c->ps2 = it->second;
    }
    CellPosition(geo, c->mod, c->lay, c->cel, c->x, c->y,
                 c->z);  // ok per fluka e geant4
  }

  for (std::map<int, dg_cell>::iterator it = map_cell.begin();
       it != map_cell.end(); ++it) {
    vec_cell.push_back(it->second);
  }
}

// simulate calorimeter responce for whole event
void DigitizeCal(TG4Event* ev, TGeoManager* geo, std::vector<dg_cell>& vec_cell)
{
  std::map<int, std::vector<pe> > photo_el;
  std::map<int, std::vector<dg_ps> > ps;
  std::map<int, double> L;

  vec_cell.clear();

  if (debug) {
    std::cout << "SimulatePE" << std::endl;
  }

  SimulatePE(ev, geo, photo_el, L);
  if (debug) {
    std::cout << "TimeAndSignal" << std::endl;
  }
  TimeAndSignal(photo_el, ps);
  if (debug) {
    std::cout << "CollectSignal" << std::endl;
  }
  CollectSignal(geo, ps, L, vec_cell);
}

// Group hits into tube
void CollectHits(TG4Event* ev, TGeoManager* geo, int NHits,
                 Int_t DetType[10000], Float_t xPos[10000], Float_t yPos[10000],
                 Float_t zPos[10000],
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

    if (flukatype == false) {
      sttname = geo->FindNode(x, y, z)->GetName();

      stid = getSTUniqID(geo, x, y, z);
      if (stid == -999) continue;
    } else {
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
          int planeid = kloe_simu::encodePlaneID(k, DetType[k]);
          stid = kloe_simu::encodeSTID(planeid, k);
          break;
        }
      }

      if (found == false) {
        std::cout << "ERROR: Point not FOUND!! " << std::endl;
        exit(1);
      }
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

// for each tube simulate tdc and adc
// tdc is the time of closest point to wire + drift time
// adc is the sum of energy deposit within integration time window
void Hits2Digit(std::map<int, std::vector<hit> >& hits2Tube,
                std::vector<dg_tube>& digit_vec)
{
  digit_vec.clear();

  for (std::map<int, std::vector<hit> >::iterator it = hits2Tube.begin();
       it != hits2Tube.end(); ++it) {
    double min_time_tub = 1E9;  // mm
    int did = it->first;

    int mod, tub, type, pla;
    double dwire = 0.;

    decodeSTID(did, pla, tub);
    decodePlaneID(pla, mod, type);

    TVector2 wire = tubePos[did];

    dg_tube d;
    d.det = it->second[0].det;
    d.did = did;
    d.de = 0;
    d.hor = (type % 2 == 0);
    d.t0 = kloe_simu::t0[pla];

    if (d.hor == true) {
      d.x = kloe_simu::stt_center[0];
      d.y = wire.Y();
      d.z = wire.X();
      dwire = d.x - 0.5 * kloe_simu::stL[did];
    } else {
      d.x = wire.Y();
      d.y = kloe_simu::stt_center[1];
      d.z = wire.X();
      dwire = d.y - 0.5 * kloe_simu::stL[did];
    }

    for (unsigned int i = 0; i < it->second.size(); i++) {
      double x1 = it->second[i].z1;
      double x2 = it->second[i].z2;
      double t1 = it->second[i].t1;
      double t2 = it->second[i].t2;

      double y1, y2;

      if (type == 2) {
        y1 = it->second[i].y1;
        y2 = it->second[i].y2;
      } else {
        y1 = it->second[i].x1;
        y2 = it->second[i].x2;
      }

      double l = getT(y1, y2, wire.Y(), x1, x2, wire.X());
      double x = x1 + (x2 - x1) * l;
      double y = y1 + (y2 - y1) * l;
      double t = t1 + (t2 - t1) * l;

      TVector2 min_dist_point(x, y);
      double min_dist_hit = (min_dist_point - wire).Mod();
      double min_time_hit =
          t + (min_dist_hit - kloe_simu::wire_radius) / kloe_simu::v_drift +
          dwire / kloe_simu::v_signal_inwire;

      if (min_time_hit < min_time_tub) min_time_tub = min_time_hit;

      if (t - d.t0 < kloe_simu::stt_int_time) d.de += it->second[i].de;

      d.hindex.push_back(it->second[i].index);
    }

    d.tdc = min_time_tub + r.Gaus(0, kloe_simu::tm_stt_smearing);
    d.adc = d.de;

    digit_vec.push_back(d);
  }
}

// simulate stt responce for whole event
void DigitizeStt(TG4Event* ev, TGeoManager* geo, int NHits,
                 Int_t DetType[10000], Float_t xPos[10000], Float_t yPos[10000],
                 Float_t zPos[10000], std::vector<dg_tube>& digit_vec)
{
  std::map<int, std::vector<hit> > hits2Tube;
  digit_vec.clear();

  CollectHits(ev, geo, NHits, DetType, xPos, yPos, zPos, hits2Tube);
  Hits2Digit(hits2Tube, digit_vec);
}
/*
void DigitizeFlukaStt(TG4Event* ev, int NHits, Int_t DetType[10000], Float_t
xPos[10000], Float_t yPos[10000], Float_t zPos[10000],
                std::vector<digit>& digit_vec)
{
    std::map<std::string, std::vector<hit> > cluster_map;
    digit_vec.clear();

    ClusterFluka(ev, NHits, DetType, xPos, yPos, zPos, cluster_map);
    Cluster2Digit(cluster_map, digit_vec);
}
*/

// digitize event
void Digitize(const char* finname, const char* foutname)
{
  TFile f(finname, "READ");

  if (TString(finname).Contains("fluka2edep") == true) {
    flukatype = true;  // dobbiamo leggere GeneratorName ..cambiare quando fatto
  }
  if (flukatype == true)
    std::cout << "This is a FLUKA SIMULATION" << std::endl;
  else
    std::cout << "This is a standard Geant4-edepsim SIMULATION" << std::endl;

  TTree* t = (TTree*)f.Get("EDepSimEvents");

  TG4Event* ev = new TG4Event;
  t->SetBranchAddress("Event", &ev);

  TGeoManager* geo = 0;

  TTree* MapTree;

  Int_t EvtNum;
  Int_t NHits;
  Float_t xHits[100000];
  Float_t yHits[100000];
  Float_t zHits[100000];
  Int_t DetType[100000];

  if (flukatype == false) {
    geo = (TGeoManager*)f.Get("EDepSimGeometry");
  } else {
    MapTree = (TTree*)f.Get("MapTree");  // geometry loaded for fluka file
    if (MapTree->GetEntries() < 1) {
      std::cout << "MapTree Empty" << std::endl;
    }

    // vector<float> xPos;
    // vector<float> yPos;
    MapTree->SetBranchAddress("EvtNum", &EvtNum);
    MapTree->SetBranchAddress("NHits", &NHits);
    MapTree->SetBranchAddress("xHits", xHits);
    MapTree->SetBranchAddress("yHits", yHits);
    MapTree->SetBranchAddress("zHits", zHits);
    MapTree->SetBranchAddress("DetType", DetType);
    t->AddFriend(MapTree);
  }
  if (debug) std::cout << "Inizializzo la geometria" << std::endl;
  init(geo);  // vale sia per geant che per fluka

  std::vector<dg_tube> digit_vec;
  std::vector<dg_cell> vec_cell;

  TFile fout(foutname, "RECREATE");
  TTree tout("tDigit", "Digitization");
  tout.Branch("dg_cell", "std::vector<dg_cell>", &vec_cell);
  tout.Branch("dg_tube", "std::vector<dg_tube>", &digit_vec);

  const int nev = t->GetEntries();

  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;

  for (int i = 0; i < nev; i++) {
    t->GetEntry(i);

    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i) / nev * 100)
              << "%]" << std::flush;

    initT0(ev);

    DigitizeCal(ev, geo, vec_cell);
    DigitizeStt(ev, geo, NHits, DetType, xHits, yHits, zHits, digit_vec);

    tout.Fill();
  }
  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;

  fout.cd();
  tout.Write();
  if (flukatype == false) geo->Write();
  fout.Close();

  f.Close();

  kloe_simu::stL.clear();
  kloe_simu::stX.clear();
  kloe_simu::stPos.clear();
  kloe_simu::t0.clear();
  kloe_simu::tubePos.clear();

}

void help_digit()
{
  std::cout << "Digitize <MC file> <digit file>" << std::endl;
  std::cout << "MC file name could contain wild card" << std::endl;
}

int main(int argc, char* argv[])
{
  if (argc != 3)
    help_digit();
  else
    Digitize(argv[1], argv[2]);
}
