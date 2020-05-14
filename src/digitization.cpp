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

#include "TG4Event.h"
#include "TG4HitSegment.h"

#include <vector>
#include <map>
#include <iostream>
#include <iomanip>

#include "struct.h"
#include "utils.h"

// Energy MeV
// Distance mm
// Time ns

using namespace kloe_simu;

TRandom3 r(0);

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
      // std::cout << "planeID out if range" << std::endl;
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

double E2PE(double E)
{
  if (debug) std::cout << "E = " << E << " -> p.e. = " << e2p2* E << std::endl;

  return e2p2 * E;
}

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

  TGeoNode* node = g->FindNode(x, y, z);

  if (node == 0) return false;

  TString str = node->GetName();
  TString str2 = g->GetPath();

  if (debug) {
    std::cout << "node name: " << str.Data() << std::endl;
  }

  if (CheckAndProcessPath(str2) == false) return false;

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

void SimulatePE(TG4Event* ev, TGeoManager* g,
                std::map<int, std::vector<double> >& time_pe,
                std::map<int, std::vector<int> >& id_hit,
                std::map<int, double>& L)
{
  int modID, planeID, cellID, id;
  double d1, d2, t0, de;

  for (std::map<std::string, std::vector<TG4HitSegment> >::iterator it =
           ev->SegmentDetectors.begin();
       it != ev->SegmentDetectors.end(); ++it) {
    if (it->first == "EMCalSci") {
      for (unsigned int j = 0; j < it->second.size(); j++) {
        if (ProcessHit(g, it->second[j], modID, planeID, cellID, d1, d2, t0,
                       de) == true) {
          double en1 = de * Attenuation(d1, planeID);
          double en2 = de * Attenuation(d2, planeID);

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
            time_pe[id].push_back(petime(t0, d1));
            id_hit[id].push_back(j);
            L[id] = d1 + d2;
          }

          for (int i = 0; i < pe2; i++) {
            time_pe[-1 * id].push_back(petime(t0, d2));
            id_hit[-1 * id].push_back(j);
            L[-1 * id] = d1 + d2;
          }
        }
      }
    }
  }
}

void TimeAndSignal(std::map<int, std::vector<double> >& time_pe,
                   std::map<int, double>& adc, std::map<int, double>& tdc)
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

  for (std::map<int, std::vector<double> >::iterator it = time_pe.begin();
       it != time_pe.end(); ++it) {
    // order by arrival time
    std::sort(it->second.begin(), it->second.end());

    int_start = it->second.front();
    pe_count = 0;
    start_index = 0;
    int index = 0;

    for (std::vector<double>::iterator pe_time = it->second.begin();
         pe_time != it->second.end(); ++pe_time) {
      // integrate for int_time
      if (*pe_time - int_start <= int_time) {
        pe_count++;
      }
      // below threshold -> reset
      else if (pe_count < pe_threshold) {
        pe_count = 1;
        int_start = *pe_time;
        start_index = pe_time - it->second.begin();
      }
      // above threshold -> stop integration and acquire
      else {
        break;
      }
    }

    if (pe_count >= pe_threshold) {
      adc[it->first] = pe2ADC * pe_count;
      index = int(costant_fraction * pe_count) + start_index;
      tdc[it->first] = it->second[index];
    }
  }
}

void CollectSignal(TGeoManager* geo,
                   std::map<int, std::vector<double> >& time_pe,
                   std::map<int, double>& adc, std::map<int, double>& tdc,
                   std::map<int, double>& L,
                   std::map<int, std::vector<int> >& id_hit,
                   std::vector<cell>& vec_cell)
{
  std::map<int, cell> map_cell;
  cell* c;

  for (std::map<int, double>::iterator it = adc.begin(); it != adc.end();
       ++it) {
    int id = abs(it->first);

    c = &(map_cell[id]);

    c->id = id;
    DecodeID(c->id, c->mod, c->lay, c->cel);
    c->l = L[c->id];

    if (it->first >= 0) {
      c->adc1 = adc[it->first];
      c->tdc1 = tdc[it->first];
      c->pe_time1 = time_pe[it->first];
      c->hindex1 = id_hit[it->first];
    } else {
      c->adc2 = adc[it->first];
      c->tdc2 = tdc[it->first];
      c->pe_time2 = time_pe[it->first];
      c->hindex2 = id_hit[it->first];
    }
    CellPosition(geo, c->mod, c->lay, c->cel, c->x, c->y, c->z);
  }

  for (std::map<int, cell>::iterator it = map_cell.begin();
       it != map_cell.end(); ++it) {
    vec_cell.push_back(it->second);
  }
}

void DigitizeCal(TG4Event* ev, TGeoManager* geo, std::vector<cell>& vec_cell)
{
  std::map<int, std::vector<double> > time_pe;
  std::map<int, std::vector<int> > id_hit;
  std::map<int, double> adc;
  std::map<int, double> tdc;
  std::map<int, double> L;

  vec_cell.clear();

  if (debug) {
    std::cout << "SimulatePE" << std::endl;
  }
  SimulatePE(ev, geo, time_pe, id_hit, L);
  if (debug) {
    std::cout << "TimeAndSignal" << std::endl;
  }
  TimeAndSignal(time_pe, adc, tdc);
  if (debug) {
    std::cout << "CollectSignal" << std::endl;
  }
  CollectSignal(geo, time_pe, adc, tdc, L, id_hit, vec_cell);
}

void Cluster(TG4Event* ev, TGeoManager* geo,
             std::map<std::string, std::vector<hit> >& cluster_map)
{
  cluster_map.clear();

  for (unsigned int j = 0; j < ev->SegmentDetectors["Straw"].size(); j++) {
    const TG4HitSegment& hseg = ev->SegmentDetectors["Straw"].at(j);

    double x = 0.5 * (hseg.Start.X() + hseg.Stop.X());
    double y = 0.5 * (hseg.Start.Y() + hseg.Stop.Y());
    double z = 0.5 * (hseg.Start.Z() + hseg.Stop.Z());

    std::string sttname = geo->FindNode(x, y, z)->GetName();

    hit h;
    h.det = sttname;
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

    std::string cluster_name(sttname);
    cluster_name += "_" + std::to_string(hseg.PrimaryId);

    cluster_map[cluster_name].push_back(h);
  }
}

void Cluster2Digit(std::map<std::string, std::vector<hit> >& cluster_map,
                   std::vector<digit>& digit_vec)
{
  for (std::map<std::string, std::vector<hit> >::iterator it =
           cluster_map.begin();
       it != cluster_map.end(); ++it) {
    digit d;
    d.de = 0;
    d.det = it->second[0].det;

    for (unsigned int k = 0; k < it->second.size(); k++) {
      d.hindex.push_back(it->second[k].index);
      d.de += it->second[k].de;
    }

    if (d.de < e_threshold) continue;

    d.hor = (d.det.find("hor") != std::string::npos) ? false : true;

    std::sort(it->second.begin(), it->second.end(), isHitBefore);
    
    if(d.hor)
    {
      d.x = 0.0;
      d.y = 0.5 * (it->second.front().y1 + it->second.back().y2) +
            r.Gaus(0., res_x);
    }
    else
    {
      d.x = 0.5 * (it->second.front().x1 + it->second.back().x2) +
            r.Gaus(0., res_x);
      d.y = 0.0;
    }

    d.t = 0.5 * (it->second.front().t1 + it->second.back().t2) +
          r.Gaus(0., res_t);
    d.z = 0.5 * (it->second.front().z1 + it->second.back().z2);

    digit_vec.push_back(d);
  }
}

void DigitizeStt(TG4Event* ev, TGeoManager* geo, std::vector<digit>& digit_vec)
{
  std::map<std::string, std::vector<hit> > cluster_map;
  digit_vec.clear();

  Cluster(ev, geo, cluster_map);
  Cluster2Digit(cluster_map, digit_vec);
}

void Digitize(const char* finname, const char* foutname)
{
  // TChain* t = new TChain("EDepSimEvents","EDepSimEvents");
  // t->Add(finname);
  // TFile f(t->GetListOfFiles()->At(0)->GetTitle());
  TFile f(finname, "READ");
  TTree* t = (TTree*)f.Get("EDepSimEvents");
  TGeoManager* geo = (TGeoManager*)f.Get("EDepSimGeometry");
  TTree* gRooTracker = (TTree*)f.Get("DetSimPassThru/gRooTracker");
  TTree* InputKinem = (TTree*)f.Get("DetSimPassThru/InputKinem");
  TTree* InputFiles = (TTree*)f.Get("DetSimPassThru/InputFiles");

  init(geo);

  TG4Event* ev = new TG4Event;
  t->SetBranchAddress("Event", &ev);

  std::vector<digit> digit_vec;
  std::vector<cell> vec_cell;

  TFile fout(foutname, "RECREATE");
  TTree tout("tDigit", "Digitization");
  tout.Branch("cell", "std::vector<cell>", &vec_cell);
  tout.Branch("Stt", "std::vector<digit>", &digit_vec);

  const int nev = t->GetEntries();

  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;

  for (int i = 0; i < nev; i++) {
    t->GetEntry(i);

    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i) / nev * 100)
              << "%]" << std::flush;

    DigitizeCal(ev, geo, vec_cell);
    DigitizeStt(ev, geo, digit_vec);

    tout.Fill();
  }
  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;

  fout.cd();
  tout.Write();
  geo->Write();
  t->CloneTree()->Write();
  if (gRooTracker) gRooTracker->CloneTree()->Write();
  if (InputKinem) InputKinem->CloneTree()->Write();
  if (InputFiles) InputFiles->CloneTree()->Write();
  fout.Close();

  f.Close();
}

void help_digit()
{
  std::cout << "Digitize <input file> <output file>" << std::endl;
  std::cout << "input file name could contain wild card" << std::endl;
}

int main(int argc, char* argv[])
{
  if (argc != 3)
    help_digit();
  else
    Digitize(argv[1], argv[2]);
}
