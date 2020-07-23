#include <TString.h>
#include <TGeoManager.h>
#include <TGeoNode.h>

#include "struct.h"

#ifndef UTILS_H
#define UTILS_H

namespace kloe_simu
{
const bool debug = false;

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

// thickness of the layers in mm
double dzlay[nLay] = {44., 44., 44., 44., 54.};
double czlay[nLay];
double cxlay[nLay][nCel];

double ec_r;
double ec_dz;

/*
////////////////////////////////////////////////////////////////////////
// geometry v0
const char* path_barrel_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/"
    "kloe_calo_volume_PV_0/ECAL_lv_PV_%d";
const char* path_endcapL_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/"
    "kloe_calo_volume_PV_0/ECAL_end_lv_PV_0";
const char* path_endcapR_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/"
    "kloe_calo_volume_PV_0/ECAL_end_lv_PV_1";
//////////////////////////////////////////////////////////////////////////
*/

////////////////////////////////////////////////////////////////////////
// geometry v1
const char* path_barrel_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/"
    "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_lv_PV_%d";
const char* path_endcapL_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/"
    "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_end_lv_PV_0";
const char* path_endcapR_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/"
    "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_end_lv_PV_1";
//////////////////////////////////////////////////////////////////////////

const char* path_internal_volume =
    "volWorld_PV/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/"
    "MagIntVol_volume_PV_0/volSTTFULL_PV_0/";
const char* name_internal_volume = "volSTTFULL_PV";

const double tscin = 3.08;
const double tscex = 0.588;
const double vlfb = 5.85;

// photoelectron/counts = 0.25
const double pe2ADC = 1 / .25;
// ADC integration time = 400 ns
const double int_time = 400.;

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

bool isCluBigger(const std::vector<digit>& v1, const std::vector<digit>& v2);
bool isDigUpstream(const digit& d1, const digit& d2);
bool isHitBefore(hit h1, hit h2);
bool isDigBefore(digit d1, digit d2);
bool isCellBefore(cell c1, cell c2);
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
void CellXYZTE(cell c, double& x, double& y, double& z, double& t, double& e);
}

#endif
