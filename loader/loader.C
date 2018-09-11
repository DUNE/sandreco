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
  double x;
  double l;
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

struct cluster {
  int tid;
  double x;
  double z;
  double y;
  double t;
  double e;
  std::vector<cell> cells;
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
  int ret_ln;
  double chi2_ln;
  int ret_cr;
  double chi2_cr;
  std::vector<digit> digits;
};

struct particle {
  int primary;
  int pdg;
  int tid;
  int parent_tid;
  double charge;
  double mass;
  double pxtrue;
  double pytrue;
  double pztrue;
  double Etrue;
  double xtrue;
  double ytrue;
  double ztrue;
  double ttrue;
  
  double pxreco;
  double pyreco;
  double pzreco;
  double Ereco;
  double xreoc;
  double yreco;
  double zreco;
  double treco;
  
  bool has_track;
  double charge_reco;
  double px_tr;
  double py_tr;
  double pz_tr;
  double E_tr;
  double x_tr;
  double y_tr;
  double z_tr;
  double t_tr;
  track tr;
  
  bool has_cluster;
  double px_cl;
  double py_cl;
  double pz_cl;
  double E_cl;
  double x_cl;
  double y_cl;
  double z_cl;
  double t_cl;
  cluster cl;
  double x_cl_var;
  double y_cl_var;
  double z_cl_var;
};

struct event {
  double Enu;
  double pxnu;
  double pynu;
  double pznu;
  double Enureco;
  double pxnureco;
  double pynureco;
  double pznureco;
  std::vector<particle> particles;
};

struct gcell {
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
  
  double ec_r;
  double ec_dz;
  
  const char* path_barrel_template = "volWorld_PV/volDetEnclosure_PV_0/volKLOEFULLECALSENSITIVE_EXTTRK_NEWGAP_PV_0/KLOEBarrelECAL_%d_volume_PV_0";
  const char* path_endcapL_template = "volWorld_PV/volDetEnclosure_PV_0/volKLOEFULLECALSENSITIVE_EXTTRK_NEWGAP_PV_0/KLOEEndcapECALL_volume_PV_0";
  const char* path_endcapR_template = "volWorld_PV/volDetEnclosure_PV_0/volKLOEFULLECALSENSITIVE_EXTTRK_NEWGAP_PV_0/KLOEEndcapECALR_volume_PV_0";
 
  TRandom3 r;
  
  const double tscin = 3.08;
  const double tscex = 0.588;
  const double vlfb = 5.85;
  
  const double lCalBarrel = 4.3; // meter
}

#ifdef __MAKECINT__ 
#pragma link C++ class std::map<int,std::vector<double> >+; 
#pragma link C++ class std::map<int,std::vector<int> >+;
#pragma link C++ class std::map<int,double>+;
#pragma link C++ class std::vector<cell>+; 
#pragma link C++ class std::map<std::string,std::vector<hit> >+; 
#pragma link C++ class std::vector<digit>+; 
#pragma link C++ class std::vector<track>+;
#pragma link C++ class std::vector<cluster>+;
#pragma link C++ class std::vector<particle>+;
#pragma link C++ class digit+;
#pragma link C++ class cell+;
#pragma link C++ class cluster+;
#pragma link C++ class track+;
#pragma link C++ class particle+;
#pragma link C++ class event+;
#endif
#endif
