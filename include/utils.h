#include "struct.h"

#include "TG4Event.h"
#include "TG4HitSegment.h"

#include <TChain.h>
#include <TGeoManager.h>
#include <TCanvas.h>
#include <TGeoTrd2.h>
#include <TGeoTube.h>

#include <iostream>

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

static const int nMod = 24;
static const int nLay = 5;
static const int nCel = 12;
static const int nCel_ec = 90;

// thickness of the layers in mm
double dzlay[nLay] = {44., 44., 44., 44., 54.};
double czlay[nLay];
double cxlay[nLay][nCel];

double ec_r;
double ec_dz;

const char* path_barrel_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/"
    "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_lv_PV_%d";
const char* path_endcapL_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/"
    "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_end_lv_PV_0";
const char* path_endcapR_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/"
    "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_end_lv_PV_1";

const double tscin = 3.08;
const double tscex = 0.588;
const double vlfb = 5.85;

// photoelectron/counts = 0.25
const double pe2ADC = 1 / .25;
// ADC integration time = 400 ns
const double int_time = 400.;

// https://www.sciencedirect.com/science/article/pii/S0168900201015029
// threshold 3-4 p.e.
const double pe_threshold = 4;

// costant fraction 15%
const double costant_fraction = 0.15;

// stt resolution and threshold
const double res_x = 0.;        // 0.2 mm
const double res_t = 0.;        // 1 ns
const double e_threshold = 0.;  // 0.2E-3 MeV

// ADC to MeV
const double adc2MeV = 1. / 10.;

const double k = 0.299792458;
const double B = 0.6;
const double GeV_to_MeV = 1000.;
const double c = k * 1E3;  // mm/ns
const double emk = 1.;
const double hadk = 1.;
}

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
}

bool isHitBefore(hit h1, hit h2)
{
  return h1.t1 < h2.t1;
}

bool isDigBefore(digit d1, digit d2)
{
  return d1.t < d2.t;
}

bool isCellBefore(cell c1, cell c2)
{
  if (c1.adc1 == 0 || c1.adc2 == 0)
    return false;
  else if (c2.adc1 == 0 || c2.adc2 == 0)
    return true;
  else
    return ((c1.tdc1 + c1.tdc2) < (c2.tdc1 + c2.tdc2));
}

bool isAfter(particle p1, particle p2)
{
  return p1.tid > p2.tid;
}

bool isBarrel(TString& str)
{
  // something like: volECALActiveSlab_21_PV_0
  return str.Contains("volECAL") == true && str.Contains("Active") == true &&
         str.Contains("end") == false;
}

bool isEndCap(TString& str)
{
  // something like: endvolECALActiveSlab_0_PV_0
  return str.Contains("endvolECAL") == true && str.Contains("Active") == true;
}

void BarrelModuleAndLayer(TString& str, TString& str2, int& modID, int& planeID)
{
  TObjArray* obja = str.Tokenize("_");    // BARERL => volECALActiveSlab_21_PV_0
  TObjArray* obja2 = str2.Tokenize("_");  // BARREL => ECAL_lv_PV_18

  int slabID;
  // top module => modID == 0
  // increasing modID counterclockwise as seen from positive x
  //(i.e. z(modID==1) < z(modID==0) & z(modID==0) < z(modID==23))
  modID = ((TObjString*)obja2->At(3))->GetString().Atoi();
  slabID = ((TObjString*)obja->At(1))->GetString().Atoi();

  delete obja;
  delete obja2;

  // planeID==0 -> smallest slab -> internal
  // planeID==208 -> biggest slab -> external
  planeID = slabID / 40;

  if (planeID > 4) planeID = 4;
}

void EndCapModuleAndLayer(TString& str, TString& str2, int& modID, int& planeID)
{
  TObjArray* obja = str.Tokenize("_");  // ENDCAP => endvolECALActiveSlab_0_PV_0
  TObjArray* obja2 = str2.Tokenize("_");  // ENDCAP => ECAL_end_lv_PV_0

  int slabID;
  modID = ((TObjString*)obja2->At(4))->GetString().Atoi();
  slabID = ((TObjString*)obja->At(1))->GetString().Atoi();

  // mod == 40 -> left
  // mod == 30 -> right
  if (modID == 0)
    modID = 40;
  else if (modID == 1)
    modID = 30;

  delete obja;
  delete obja2;

  // planeID==0 -> internal
  // planeID==208 -> external
  planeID = slabID / 40;

  if (planeID > 4) planeID = 4;
}

void BarrelCell(double x, double y, double z, TGeoManager* g, TGeoNode* node,
                int& cellID, double& d1, double& d2)
{
  double Pmaster[3];
  double Plocal[3];
  Pmaster[0] = x;
  Pmaster[1] = y;
  Pmaster[2] = z;

  g->GetCurrentNavigator()->MasterToLocal(Pmaster, Plocal);

  TGeoTrd2* trd = (TGeoTrd2*)node->GetVolume()->GetShape();

  if (kloe_simu::debug) {
    std::cout << "pointer: " << trd << std::endl;
  }

  double dx1 = trd->GetDx1();
  double dx2 = trd->GetDx2();
  double dz = trd->GetDz();
  double dy1 = trd->GetDy1();
  double dy2 = trd->GetDy2();

  // d1 distanza da estremo left (x<0)
  // d2 distanza da estremo right (x>0)
  d1 = dy1 + Plocal[1];
  d2 = dy1 - Plocal[1];

  // http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/html/Detector/Geometry/geomSolids.html
  // if z = -dz -> dx = 2*dx1
  // if z =  dz -> dx = 2*dx2
  // semilarghezza della slab di scintillatore alla quota Plocal[2]
  double dx = 0.5 * Plocal[2] / dz * (dx2 - dx1) + 0.5 * (dx2 + dx1);

  // Cell width at z = Plocal[2]
  double cellw = 2. * dx / kloe_simu::nCel;

  // cellID = distanza dall'estremo diviso larghezza cella
  cellID = (Plocal[0] + dx) / cellw;
}

void EndCapCell(double x, double y, double z, TGeoManager* g, TGeoNode* node,
                int& cellID, double& d1, double& d2)
{
  double Pmaster[3];
  double Plocal[3];
  Pmaster[0] = x;
  Pmaster[1] = y;
  Pmaster[2] = z;

  g->GetCurrentNavigator()->MasterToLocal(Pmaster, Plocal);

  TGeoTube* tub = (TGeoTube*)node->GetVolume()->GetShape();

  if (kloe_simu::debug) {
    std::cout << "pointer: " << tub << std::endl;
  }

  double rmin = tub->GetRmin();
  double rmax = tub->GetRmax();
  double dz = tub->GetDz();

  // d1 distanza da estremo up (y>0)
  // d2 distanza da estremo down (y<0)
  d1 = rmax * TMath::Sin(TMath::ACos(Plocal[0] / rmax)) - Plocal[1];
  d2 = rmax * TMath::Sin(TMath::ACos(Plocal[0] / rmax)) + Plocal[1];

  cellID = int((Plocal[0] / rmax + 1.) * kloe_simu::nCel_ec * 0.5);
}

bool CheckAndProcessPath(TString& str2)
{
  // ENDCAP ==> something like:
  // "/volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_lv_PV_18/volECALActiveSlab_21_PV_0"
  // BARREL ==> ìsomething like:
  // "/volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_end_lv_PV_0/endvolECALActiveSlab_0_PV_0"
  TObjArray* obj = str2.Tokenize("/");

  int size = obj->GetEntries();
  if (size < 8) {
    return false;
  };

  // BARREL => ECAL_lv_PV_18
  // ENDCAP => ECAL_end_lv_PV_0
  str2 = ((TObjString*)obj->At(6))->GetString();
  delete obj;

  return true;
}

void CellPosition(TGeoManager* geo, int mod, int lay, int cel, double& x,
                  double& y, double& z)
{
  x = 0;
  y = 0;
  z = 0;

  double dummyLoc[3];
  double dummyMas[3];

  if (mod < 24) {
    dummyLoc[0] = kloe_simu::cxlay[lay][cel];
    dummyLoc[1] = 0.;
    dummyLoc[2] = kloe_simu::czlay[lay];

    geo->cd(TString::Format(kloe_simu::path_barrel_template, mod).Data());

  } else if (mod == 30 || mod == 40)
      // right x > 0 : c->mod = 30
      // left  x < 0 : c->mod = 40
  {
    dummyLoc[0] =
        kloe_simu::ec_r / kloe_simu::nCel_ec * (0.5 + cel) - kloe_simu::ec_r;
    dummyLoc[1] = 0.;
    dummyLoc[2] = kloe_simu::czlay[lay];

    if (mod == 30) {
      geo->cd(kloe_simu::path_endcapR_template);
    } else if (mod == 40) {
      geo->cd(kloe_simu::path_endcapL_template);
    }
  }

  geo->LocalToMaster(dummyLoc, dummyMas);

  x = dummyMas[0];
  y = dummyMas[1];
  z = dummyMas[2];
}

void init(TGeoManager* geo)
{
  TGeoTrd2* mod = (TGeoTrd2*)geo->FindVolumeFast("ECAL_lv_PV")->GetShape();

  // https://root.cern.ch/root/htmldoc/guides/users-guide/Geometry.html#shapes
  // GetDx1() half length in x at -Dz
  // GetDx2() half length in x at +Dz
  // Dx1 < Dx2 => -Dz corresponds to minor width => internal side
  double xmin = mod->GetDx1();
  double xmax = mod->GetDx2();
  double dz = mod->GetDz();

  double m = 0.5 * (xmax - xmin) / dz;
  double q = 0.5 * (xmax + xmin);

  // z edge of the cells
  double zlevel[kloe_simu::nLay + 1];
  zlevel[0] = -dz;

  for (int i = 0; i < kloe_simu::nLay; i++) {
    zlevel[i + 1] = zlevel[i] + kloe_simu::dzlay[i];
  }

  // z position of the center of the cells
  for (int i = 0; i < kloe_simu::nLay; i++) {
    kloe_simu::czlay[i] = 0.5 * (zlevel[i] + zlevel[i + 1]);

    // total module width at the z position of the center of the cell
    double xwidth = 2 * (m * kloe_simu::czlay[i] + q);

    // cell width at the z position of the center of the cell
    double dx = xwidth / kloe_simu::nCel;

    // x position of the center of the cells
    for (int j = 0; j < kloe_simu::nCel; j++) {
      kloe_simu::cxlay[i][j] = dx * (j + 0.5) - xwidth * 0.5;
    }
  }

  TGeoTube* ec = (TGeoTube*)geo->FindVolumeFast("ECAL_end_lv_PV")->GetShape();

  kloe_simu::ec_r = ec->GetRmax();
  kloe_simu::ec_dz = ec->GetDz();
}

int EncodeID(int mod, int lay, int cel)
{
  return cel + 100 * lay + 1000 * mod;
}

void DecodeID(int id, int& mod, int& lay, int& cel)
{
  mod = id / 1000;
  lay = (id - mod * 1000) / 100;
  cel = id - mod * 1000 - lay * 100;
}

double mindist(double s1x, double s1y, double s1z, double s2x, double s2y,
               double s2z, double px, double py, double pz)
{
  double segmod = (s1x - s2x) * (s1x - s2x) + (s1y - s2y) * (s1y - s2y) +
                  (s1z - s2z) * (s1z - s2z);

  double prod = (px - s1x) * (s2x - s1x) + (py - s1y) * (s2y - s1y) +
                (pz - s1z) * (s2z - s1z);

  double t = std::min(std::max(prod / segmod, 0.), 1.);

  double s3x = s1x + (s2x - s1x) * t;
  double s3y = s1y + (s2y - s1y) * t;
  double s3z = s1z + (s2z - s1z) * t;

  return sqrt((px - s3x) * (px - s3x) + (py - s3y) * (py - s3y) +
              (pz - s3z) * (pz - s3z));
}

double angle(double x1, double y1, double z1, double x2, double y2, double z2)
{
  double prod = x1 * x2 + y1 * y2 + z1 * z2;
  double mag1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
  double mag2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);

  return TMath::ACos(prod / (mag1 * mag2));
}

double AttenuationFactor(double d, int planeID)
{
  /*
       dE/dx attenuation - Ea=p1*exp(-d/atl1)+(1.-p1)*exp(-d/atl2)
         d    distance from photocatode - 2 cells/cell; d1 and d2
        atl1  50. cm
        atl2  430 cm planes 1-2    innermost plane is 1
              380 cm plane 3
              330 cm planes 4-5
         p1   0.35
  */
  double atl2 = 0.0;

  switch (planeID) {
    case 0:
    case 1:
      atl2 = kloe_simu::atl2_01;
      break;

    case 2:
      atl2 = kloe_simu::atl2_2;
      break;

    case 3:
    case 4:
      atl2 = kloe_simu::atl2_34;
      break;

    default:
      // std::cout << "planeID out if range" << std::endl;
      atl2 = -999.0;
      break;
  }

  if (kloe_simu::debug) {
    std::cout << "planeID = " << planeID << std::endl;
    std::cout << "\tp1   = " << kloe_simu::p1 << std::endl;
    std::cout << "\talt1 = " << kloe_simu::atl1 << std::endl;
    std::cout << "\talt2 = " << atl2 << std::endl;
    std::cout << "\tatt  = "
              << kloe_simu::p1* TMath::Exp(-d / kloe_simu::atl1) +
                     (1. - kloe_simu::p1) * TMath::Exp(-d / atl2) << std::endl;
  }

  return kloe_simu::p1 * TMath::Exp(-d / kloe_simu::atl1) +
         (1. - kloe_simu::p1) * TMath::Exp(-d / atl2);
}

double TfromTDC(double t1, double t2, double L)
{
  return 0.5 * (t1 + t2 - kloe_simu::vlfb * L / kloe_simu::m_to_mm);
}

double XfromTDC(double t1, double t2)
{
  return 0.5 * (t1 - t2) / kloe_simu::vlfb * kloe_simu::m_to_mm;
}

double EfromADC(double adc1, double adc2, double d1, double d2, int planeID)
{
  double f1 = AttenuationFactor(d1, planeID);
  double f2 = AttenuationFactor(d2, planeID);

  return 0.5 * (adc1 / f1 + adc2 / f2) * kloe_simu::adc2MeV;
}

void CellXYZTE(cell c, double& x, double& y, double& z, double& t, double& e)
{
  if (c.id < 25000)  // Barrel
  {
    x = c.x + XfromTDC(c.tdc1, c.tdc2);
    y = c.y;
  } else {
    x = c.x;
    y = c.y - XfromTDC(c.tdc1, c.tdc2);
  }
  double d1 = 0.5 * c.l + XfromTDC(c.tdc1, c.tdc2);
  double d2 = 0.5 * c.l - XfromTDC(c.tdc1, c.tdc2);
  z = c.z;
  t = TfromTDC(c.tdc1, c.tdc2, c.l);
  e = EfromADC(c.adc1, c.adc2, d1, d2, c.lay);
}

#endif