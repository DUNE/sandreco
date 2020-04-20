#include "struct.h"

#include "TG4Event.h"
#include "TG4HitSegment.h"

#include <TChain.h>
#include <TGeoManager.h>
#include <TCanvas.h>

#ifndef UTILS_H
#define UTILS_H

bool isHitBefore(hit h1, hit h2)
{
  return h1.t1 < h2.t1;
}

bool isDigBefore(digit d1, digit d2)
{
  return d1.t < d2.t;
}

namespace ns_Digit
{
const bool debug = false;

static const int nMod = 24;
static const int nLay = 5;
static const int nCel = 12;

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

const double lCalBarrel = 4.3;  // meter
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

#endif