// File loader.C 
#include <vector> 
#include <map>

#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4Event.h"
#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4HitSegment.h"

#include <TRandom3.h>
#include <TChain.h>
#include <TGeoManager.h>
#include <TCanvas.h>

#ifndef LOADER_C
#define LOADER_C

struct cell {
  int id;
  double z;
  double y;
  double adc1;
  double tdc1;
  double adc2;
  double tdc2;
  int mod;
  int lay;
  int cel;
  std::vector<double> pe_time1;
  std::vector<int> hindex1;
  std::vector<double> pe_time2;
  std::vector<int> hindex2;
};

struct hit {
  std::string det;
  double x1;
  double y1;
  double z1;
  double t1;
  double x2;
  double y2;
  double z2;
  double t2;
  double de;
  int pid;
  int index;
};

struct digit {
  std::string det;
  double x;
  double y;
  double z;
  double t;
  double de;
  bool hor;
  std::vector<int> hindex;
};

struct track {
  int tid;
  double yc;
  double zc;
  double r;
  double a;
  double b;
  double h;
  double x0;
  double y0;
  double z0;
  double t0;
  std::vector<digit> digits;
};

class gcell {
  public:
    int id;
    double Z[4];
    double Y[4];
    double adc;
    double tdc;
};

bool isHitBefore(hit h1, hit h2)
{
  return h1.t1 < h2.t1;
}

bool isDigBefore(digit d1, digit d2)
{
  return d1.t < d2.t;
}

namespace ns_Digit {
  const bool debug = false;
  
  static const int nMod = 24;
  static const int nLay = 5;
  static const int nCel = 12;
    
  double dzlay[ns_Digit::nLay+1] = {115, 115-22, 115-22-22, 115-22-22-22, 115-22-22-22-22, 115-22-22-22-22-27};
  double czlay[ns_Digit::nLay];
  double cxlay[ns_Digit::nLay][ns_Digit::nCel];
  
  const char* path_template = "volWorld_PV/volDetEnclosure_PV_0/volKLOEFULLECALSENSITIVE_EXTTRK_NEWGAP_PV_0/KLOEBarrelECAL_%d_volume_PV_0";
 
  TRandom3 r;
}

namespace ns_draw {
  const bool debug = false;
  
  static const int nMod = 24;
  static const int nLay = 5;
  static const int nCel = 12;
  
  static const int nTotCells = nMod * nLay * nCel; 
  static const int nCellModule = nLay * nCel;
  
  double centerKLOE[3];
  double CellLocalX[nCellModule][4];
  double CellLocalZ[nCellModule][4];
  
  int palette = 87;
  
  bool initialized = false;
  
  double dwx = 2500.;
  double dwy = 2500.;
  double dwz = 2500.;
  
  TChain* t = 0;
  TG4Event* ev = new TG4Event;
  TGeoManager* geo = 0;
  TCanvas* cev = 0;

  std::vector<cell>* vec_cell;
  std::vector<digit>* vec_digi;
  std::vector<track>* vec_tr;
  std::map<int, gcell> calocell;
}

#endif

#ifdef __MAKECINT__ 
#pragma link C++ class std::map<int,std::vector<double> >+; 
#pragma link C++ class std::map<int,std::vector<int> >+;
#pragma link C++ class std::map<int,double>+;
#pragma link C++ class std::vector<cell>+; 
#pragma link C++ class std::map<std::string,std::vector<hit> >+; 
#pragma link C++ class std::vector<digit>+; 
#pragma link C++ class std::vector<track>+;
#endif
