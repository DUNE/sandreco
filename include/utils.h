#include <TString.h>
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TPRegexp.h>
#include <TVector2.h>

#include "struct.h"

#include <TG4Event.h>

#ifndef UTILS_H
#define UTILS_H

namespace kloe_simu
{
const bool debug = false;

bool flukatype = false;  // for FLUKA

const double mm_to_m = 1E-3;
const double m_to_mm = 1000.;

/*
     dE/dx attenuation - Ea=p1*exp(-d/atl1)+(1.-p1)*exp(-d/atl2)
       d    distance from photocatode - 2 cells/cell; d1 and d2
      atl1  50. cm
      atl2  430 cm planes 1-2    innermost plane is 1
            380 cm plane 3
            330 cm planes 4-5
       p1   0.35
*/
const double p1 = 0.35;
const double atl1 = 500.;
const double atl2_01 = 4300.0;
const double atl2_2 = 3800.0;
const double atl2_34 = 3300.0;

// Average number of photoelectrons = 25*Ea(MeV)
// corrected to 18.5 to have mean number of pe of 40
// for mip crossing in the middle of barrel module
const double e2p2 = 18.5;

const int nMod = 24;
const int nLay = 5;
const int nCel = 12;
const int nCel_ec = 90;

// ecal dimension for fluka
const double xmin_f = 262.55;
const double xmax_f = 292.85;
const double dz_f = 115.0;

const double ec_rf = 2000.0;  // ad essere precisi nella realtà è 1980
const double ec_dzf = 115.0;
const double lCalBarrel = 4300;

const double kloe_int_R_f = 2000.;
const double kloe_int_dx_f = 1690.;

// coordinates of the cells for FLUKA
extern double cellCoordBarrel[nMod][nLay][nCel][3];
extern double cellCoordEndcap[5][nLay][90][3];

// thickness of the layers in mm
const double dzlay[nLay] = {44., 44., 44., 44., 54.};
extern double czlay[nLay];
extern double cxlay[nLay][nCel];

extern double ec_r;
extern double ec_dz;

/*
////////////////////////////////////////////////////////////////////////
// geometry v0
const char* const path_barrel_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/"
    "kloe_calo_volume_PV_0/ECAL_lv_PV_%d";
const char* const path_endcapL_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/"
    "kloe_calo_volume_PV_0/ECAL_end_lv_PV_0";
const char* const path_endcapR_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/"
    "kloe_calo_volume_PV_0/ECAL_end_lv_PV_1";
//////////////////////////////////////////////////////////////////////////
*/

////////////////////////////////////////////////////////////////////////
// geometry v1
const char* const path_barrel_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/"
    "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_lv_PV_%d";
const char* const path_endcapL_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/"
    "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_end_lv_PV_0";
const char* const path_endcapR_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/"
    "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_end_lv_PV_1";
//////////////////////////////////////////////////////////////////////////

const char* const path_internal_volume =
    "volWorld_PV/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/"
    "MagIntVol_volume_PV_0/volSTTLAR_PV_0/";
const char* const name_internal_volume = "volSTTLAR_PV";

const char* const rST_string =
    "_(C3H6|C|Tr)Mod_([0-9]+)_(ST_|)(hor|ver|hor2)_ST_stGas_(Xe|Ar)19_vol_PV_(["
    "0-9]+)";
const char* const rSTplane_string =
    "_(C3H6|C|Tr)Mod_([0-9]+)_(ST_|)(hor|ver|hor2)_vol_PV_0";

extern TPRegexp* rST;
extern TPRegexp* rSTplane;

const double tscin = 3.08;
const double tscex = 0.588;
const double vlfb = 5.85;

/*
// da qui in poi non ci sono più nella master
//


namespace ns_Draw
{
const bool debug = false;

static const int nMod = 24;
static const int nLay = 5;
static const int nCel = 12;
static const int nLay_ec = 5;
static const int nCel_ec = 90;

static const int nTotCells = nMod * nLay * nCel;
static const int nCellModule = nLay * nCel;

static const double dt = 500;

double centerKLOE[3];
double CellLocalX[nCellModule][4];
double CellLocalZ[nCellModule][4];

int palette = 87;

bool initialized = false;

double dwx = 2500.;
double dwy = 2500.;
double dwz = 2500.;

double kloe_int_R = 2000.;
double kloe_int_dx = 1690.;

TChain* t = 0;
TG4Event* ev = new TG4Event;
TGeoManager* geo = 0;
TCanvas* cev = 0;
TCanvas* cpr = 0;

std::vector<cell>* vec_cell;
std::vector<digit>* vec_digi;
std::vector<track>* vec_tr;
std::vector<cluster>* vec_cl;
std::map<int, gcell> calocell;
=======
*/

// photoelectron/counts = 0.25
const double pe2ADC = 1 / .25;
// ADC integration time = 30 ns
const double int_time = 30.;
// dead time
const double dead_time = 50.;

// https://www.sciencedirect.com/science/article/pii/S0168900201015029
// threshold 3-4 p.e. at 2 m distance
const double pe_threshold = 2.5;

// costant fraction 15%
const double costant_fraction = 0.15;

// stt resolution and threshold
const double res_x = 0.2;            // 0.2 mm
const double res_t = 0.;             // 1 ns
const double e_threshold = 0.25E-3;  // 0.25E-3 MeV

// ADC to MeV
const double adc2MeV = 1. / 10.;

const double k = 0.299792458;
const double B = 0.6;
const double GeV_to_MeV = 1000.;
const double c = k * 1E3;  // mm/ns
const double emk = 1.;
const double hadk = 1.;

// STT
extern double stt_center[3];
const double stt_int_time = 400.;     // ns
const double bucket_rms = 1.;         // ns
const double wire_radius = 0.02;      // mm
const double v_drift = 0.05;          // mm/ns
const double v_signal_inwire = 200.;  // mm/ns
const double tm_stt_smearing = 3.5;   // ns

extern std::map<int, std::map<double, int> > stX;
extern std::map<int, double> stL;
extern std::map<int, std::map<int, TVector2> > stPos;
extern std::map<int, TVector2> tubePos;
extern std::map<int, double> t0;

bool isPeBefore(const pe& p1, const pe& p2);
bool isCluBigger(const std::vector<dg_tube>& v1,
                 const std::vector<dg_tube>& v2);
bool isDigUpstream(const dg_tube& d1, const dg_tube& d2);
bool isHitBefore(hit h1, hit h2);
bool isDigBefore(dg_tube d1, dg_tube d2);
bool isCellBefore(dg_cell c1, dg_cell c2);
bool isAfter(particle p1, particle p2);
bool isBarrel(TString& str);
bool isEndCap(TString& str);
void BarrelModuleAndLayer(TString& str, TString& str2, int& modID,
                          int& planeID);
void EndCapModuleAndLayer(TString& str, TString& str2, int& modID,
                          int& planeID);
void BarrelCell(double x, double y, double z, TGeoManager* g, TGeoNode* node,
                int& cellID, double& d1, double& d2);
void EndCapCell(double x, double y, double z, TGeoManager* g, TGeoNode* node,
                int& cellID, double& d1, double& d2);
bool CheckAndProcessPath(TString& str2);
void CellPosition(TGeoManager* geo, int mod, int lay, int cel, double& x,
                  double& y, double& z);
void init(TGeoManager* geo);
int EncodeID(int mod, int lay, int cel);
void DecodeID(int id, int& mod, int& lay, int& cel);
double mindist(double s1x, double s1y, double s1z, double s2x, double s2y,
               double s2z, double px, double py, double pz);
double angle(double x1, double y1, double z1, double x2, double y2, double z2);
double AttenuationFactor(double d, int planeID);
double TfromTDC(double t1, double t2, double L);
double XfromTDC(double t1, double t2);
double EfromADC(double adc1, double adc2, double d1, double d2, int planeID);
void CellXYZTE(dg_cell c, double& x, double& y, double& z, double& t,
               double& e);

bool isST(TString name);
bool isSTPlane(TString name);
int getSTId(TString name);
int getPlaneID(TString name);
void getSTinfo(TGeoNode* nod, TGeoHMatrix mat, int pid,
               std::map<double, int>& stX, std::map<int, double>& stL,
               std::map<int, TVector2>& stPos);
void getSTPlaneinfo(TGeoNode* nod, TGeoHMatrix mat,
                    std::map<int, std::map<double, int> >& stX,
                    std::map<int, double>& stL,
                    std::map<int, std::map<int, TVector2> >& stPos);
int getSTUniqID(TGeoManager* g, double x, double y, double z);
int encodeSTID(int planeid, int tubeid);
void decodeSTID(int id, int& planeid, int& tubeid);
int encodePlaneID(int moduleid, int type);
void decodePlaneID(int id, int& moduleid, int& type);
double getT(double y1, double y2, double y, double z1, double z2, double z);
void initT0(TG4Event* ev);
}

#endif
