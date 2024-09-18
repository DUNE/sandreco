#include <TApplication.h>
#include <TArrow.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TFile.h>
#include <TGeoEltu.h>
#include <TGeoManager.h>
#include <TGeoTrd2.h>
#include <TGeoTube.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TLine.h>
#include <TMarker.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <iomanip>
#include <iostream>
#include <map>
#include <stdlib.h>

#include "../include/struct.h"
#include "../include/utils.h"

#include "TG4Event.h"
#include "TG4HitSegment.h"

namespace display
{
static const int nMod = 24;
static const int nLay = 5;
static const int nCel = 12;
static const int nLay_ec = 5;
static const int nCel_ec = 90;

static const int nTotCells = nMod * nLay * nCel;
static const int nCellModule = nLay * nCel;

static const double dt = 500;

const bool debug = false;

bool initialized = false;

bool showtrj      = false;
bool showtrj2reco = false;
bool showede      = false;
bool showede2reco = false;
bool showdig      = false;
bool showrec      = false;
bool compare      = false;

bool is_ev_number_set = false;
bool is_mc_file_set = false;
bool is_out_file_set = false;
bool is_batch_mode_set = false;

int palette = 87;

double centerKLOE[3];
double CellLocalX[nCellModule][4];
double CellLocalZ[nCellModule][4];
double centerGRAIN[3];

double GRAIN_dx;
double GRAIN_dy;
double GRAIN_dz;

double dwx = 2500.;
double dwy = 2500.;
double dwz = 2500.;

double kloe_int_R = 2000.;
double kloe_int_dx = 1690.;

//files
TFile* f = 0;
TFile* fmc = 0;
std::vector<TFile*> vf;
TString fout;

TTree* t = 0;
TTree* tEvent = 0;
TTree* tReco = 0;
TTree* tDigit = 0;

TG4Event* ev = new TG4Event;
event* evt = new event;
TGeoManager* geo = 0;
TCanvas* cev = 0;
TCanvas* cpr = 0;

std::vector<dg_cell>* vec_cell = new std::vector<dg_cell>;
std::vector<dg_tube>* vec_tube = new std::vector<dg_tube>;
std::vector<dg_wire>* vec_wire = new std::vector<dg_wire>;
std::vector<track>* vec_tr = new std::vector<track>;
std::vector<cluster>* vec_cl = new std::vector<cluster>;
std::map<int, gcell> calocell;

const char* path_intreg =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
    "MagIntVol_volume_PV_0/sand_inner_volume_PV_0";

const char* path_barrel_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
    "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_lv_PV_%d";

const char* path_endcapR_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
    "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_end_lv_PV_1";

const char* path_endcapL_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
    "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_end_lv_PV_0";

const char* path_GRAIN =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
    "MagIntVol_volume_PV_0/sand_inner_volume_PV_0/GRAIN_lv_PV_0/"
    "GRAIN_Ext_vessel_outer_layer_lv_PV_0/"
    "GRAIN_Honeycomb_layer_lv_PV_0/GRAIN_Ext_vessel_inner_layer_lv_PV_0/"
    "GRAIN_gap_between_vessels_lv_PV_0/"
    "GRAIN_inner_vessel_lv_PV_0/GRAIN_LAr_lv_PV_0";

const char* path_GRIAN =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
    "MagIntVol_volume_PV_0/sand_inner_volume_PV_0/GRAIN_lv_PV_0/"
    "GRAIN_Ext_vessel_outer_layer_lv_PV_0/"
    "GRAIN_Honeycomb_layer_lv_PV_0/GRAIN_Ext_vessel_inner_layer_lv_PV_0/"
    "GRAIN_gap_between_vessels_lv_PV_0/"
    "GRAIN_inner_vessel_lv_PV_0/GRIAN_LAr_lv_PV_0";

const char* barrel_mod_vol_name = "ECAL_lv_PV";
const char* endcap_mod_vol_name = "ECAL_end_lv_PV";
const char* GRAIN_vol_name = "GRAIN_LAr_lv_PV";
const char* GRIAN_vol_name = "GRIAN_LAr_lv_PV";

}  // namespace display

using namespace display;

void init()
{
  gStyle->SetPalette(palette);

  TTree* tEdep = reinterpret_cast<TTree*>(fmc->Get("EDepSimEvents"));
  if (!tEdep) return;
  tEdep->SetBranchAddress("Event", &ev);

  TTree* tGenie =
      reinterpret_cast<TTree*>(fmc->Get("DetSimPassThru/gRooTracker"));
  if (tGenie) tEdep->AddFriend(tGenie);

  geo = reinterpret_cast<TGeoManager*>(fmc->Get("EDepSimGeometry"));
  if (!geo) return;

  for (auto f : vf) {
    TTree* tt = nullptr;
    tt = reinterpret_cast<TTree*>(f->Get("tEvent"));
    if (tt) tEvent = tt;
    tt = reinterpret_cast<TTree*>(f->Get("tReco"));
    if (tt) tReco = tt;
    tt = reinterpret_cast<TTree*>(f->Get("tDigit"));
    if (tt) tDigit = tt;
  }

  if (tEvent) tEdep->AddFriend(tEvent);
  if (tEvent) tEvent->SetBranchAddress("event", &evt);

  if (tDigit) tEdep->AddFriend(tDigit);
  //if (tDigit) tDigit->SetBranchAddress("dg_cell", &vec_cell);
  if(geo->FindVolumeFast("STTtracker_PV")){
    //stt based digitization
    if (tDigit) tDigit->SetBranchAddress("dg_tube", &vec_tube);
  }else{
    //drift based digitization
    if (tDigit) tDigit->SetBranchAddress("dg_wire", &vec_wire);
  }

  if (tReco) tEdep->AddFriend(tReco);
  if (tReco) tReco->SetBranchAddress("track", &vec_tr);
  if (tReco) tReco->SetBranchAddress("cluster", &vec_cl);

  t = tEdep;

  double dummyLoc[3];
  double dummyMas[3];

  geo->cd(path_intreg);

  dummyLoc[0] = 0.;
  dummyLoc[1] = 0.;
  dummyLoc[2] = 0.;
  geo->LocalToMaster(dummyLoc, centerKLOE);

  double dzlay[nLay] = {44., 44., 44., 44., 54.};
  double dx1[nLay];
  double dx2[nLay];

  TGeoTrd2* mod =
      (TGeoTrd2*)geo->FindVolumeFast(barrel_mod_vol_name)->GetShape();

  double xmin = mod->GetDx1();
  double xmax = mod->GetDx2();
  double dz = mod->GetDz();

  if (debug) {
    std::cout << dz << " " << xmax << " " << xmin << std::endl;
  }

  double m = 0.5 * (xmax - xmin) / dz;
  double q = 0.5 * (xmax + xmin);

  // z edge of the cells
  double zlevel[nLay + 1];
  zlevel[0] = -dz;

  for (int i = 0; i < nLay; i++) {
    zlevel[i + 1] = zlevel[i] + dzlay[i];
  }

  for (int i = 0; i < nLay; i++) {
    dx1[i] = 2 * (m * zlevel[i] + q);
    dx2[i] = 2 * (m * zlevel[i + 1] + q);
  }

  if (debug) {
    for (int i = 0; i < nLay; i++) {
      std::cout << dx1[i] << " " << dx2[i] << " " << zlevel[i] << " "
                << zlevel[i + 1] << std::endl;
    }
  }

  for (int i = 0; i < nLay; i++) {
    for (int j = 0; j < nCel; j++) {
      // from bottom-left to top-right

      CellLocalX[i * nCel + j][0] = dx1[i] * (j / 12. - 0.5);
      CellLocalX[i * nCel + j][1] = dx1[i] * ((j + 1) / 12. - 0.5);
      CellLocalX[i * nCel + j][2] = dx2[i] * ((j + 1) / 12. - 0.5);
      CellLocalX[i * nCel + j][3] = dx2[i] * (j / 12. - 0.5);

      CellLocalZ[i * nCel + j][0] = zlevel[i];
      CellLocalZ[i * nCel + j][1] = zlevel[i];
      CellLocalZ[i * nCel + j][2] = zlevel[i + 1];
      CellLocalZ[i * nCel + j][3] = zlevel[i + 1];

      if (debug)
        std::cout << CellLocalZ[i * nCel + j][0] << " "
                  << CellLocalX[i * nCel + j][0] << " "
                  << CellLocalZ[i * nCel + j][1] << " "
                  << CellLocalX[i * nCel + j][1] << " "
                  << CellLocalZ[i * nCel + j][2] << " "
                  << CellLocalX[i * nCel + j][2] << " "
                  << CellLocalZ[i * nCel + j][3] << " "
                  << CellLocalX[i * nCel + j][3] << std::endl;
    }
  }

  double CellMasterY[nCellModule][4];
  double CellMasterZ[nCellModule][4];

  if (debug) {
    TCanvas* cd1 = new TCanvas("cd1", "", 700, 700);
    cd1->DrawFrame(centerKLOE[2] - 2500, centerKLOE[1] - 2500,
                   centerKLOE[2] + 2500, centerKLOE[1] + 2500);
  }

  for (int i = 0; i < nMod; i++) {
    geo->cd(TString::Format(path_barrel_template, i).Data());

    if (debug)
      std::cout << "node: " << i << " " << geo->GetCurrentNode() << " "
                << geo->GetCurrentNode()->GetName() << " "
                << TString::Format(path_barrel_template, i).Data() << std::endl;

    for (int j = 0; j < nLay; j++) {
      for (int k = 0; k < nCel; k++) {

        int index = i * (nLay * nCel) + j * (nCel) + k;
        int id = k + 100 * j + 1000 * i;

        int local_index = j * nCel + k;

        calocell[id].id = id;

        if (debug)
          std::cout << i << " " << j << " " << k << " " << index << " " << id
                    << " " << local_index << " " << nMod << " " << nLay << " "
                    << nCel << std::endl;

        for (int mm = 0; mm < 4; mm++) {
          dummyLoc[0] = CellLocalX[local_index][mm];
          dummyLoc[1] = 0.;
          dummyLoc[2] = CellLocalZ[local_index][mm];

          geo->LocalToMaster(dummyLoc, dummyMas);

          if (debug) {
            std::cout << "local : " << dummyLoc[0] << " " << dummyLoc[1] << " "
                      << dummyLoc[2] << std::endl;
            std::cout << "master: " << dummyMas[0] << " " << dummyMas[1] << " "
                      << dummyMas[2] << std::endl;
          }

          calocell[id].Y[mm] = dummyMas[1];
          calocell[id].Z[mm] = dummyMas[2];

          CellMasterY[local_index][mm] = dummyMas[1];
          CellMasterZ[local_index][mm] = dummyMas[2];
        }

        if (debug) {
          TGraph* gr1 =
              new TGraph(4, CellMasterZ[local_index], CellMasterY[local_index]);
          gr1->Draw("f");
        }
      }
    }
  }

  if (debug) {
    TCanvas* cd2 = new TCanvas("cd2", "", 700, 700);
    cd2->DrawFrame(centerKLOE[2] - 2500, centerKLOE[0] - 2500,
                   centerKLOE[2] + 2500, centerKLOE[0] + 2500);
  }

  TGeoTube* ec =
      (TGeoTube*)geo->FindVolumeFast(endcap_mod_vol_name)->GetShape();

  double rmax = ec->GetRmax();
  double rmin = ec->GetRmin();
  double dz_ec = ec->GetDz();

  double dummyLoc_ec[4][3];

  geo->cd(path_endcapR_template);

  for (int j = 0; j < nLay_ec; j++) {
    for (int k = 0; k < nCel_ec; k++) {
      int id = k + 100 * j + 1000 * 30;

      calocell[id].id = id;

      dummyLoc_ec[0][0] = rmax / 45. * k - rmax;
      dummyLoc_ec[0][1] = 0.;
      dummyLoc_ec[0][2] = zlevel[j];

      dummyLoc_ec[1][0] = rmax / 45. * k - rmax;
      dummyLoc_ec[1][1] = 0.;
      dummyLoc_ec[1][2] = zlevel[j + 1];

      dummyLoc_ec[2][0] = rmax / 45. * (k + 1) - rmax;
      dummyLoc_ec[2][1] = 0.;
      dummyLoc_ec[2][2] = zlevel[j + 1];

      dummyLoc_ec[3][0] = rmax / 45. * (k + 1) - rmax;
      dummyLoc_ec[3][1] = 0.;
      dummyLoc_ec[3][2] = zlevel[j];

      for (int mm = 0; mm < 4; mm++) {
        geo->LocalToMaster(dummyLoc_ec[mm], dummyMas);

        if (debug) {
          std::cout << "local : " << dummyLoc_ec[mm][0] << " "
                    << dummyLoc_ec[mm][1] << " " << dummyLoc_ec[mm][2]
                    << std::endl;
          std::cout << "master: " << dummyMas[0] << " " << dummyMas[1] << " "
                    << dummyMas[2] << std::endl;
        }

        calocell[id].Y[mm] = dummyMas[0];
        calocell[id].Z[mm] = dummyMas[2];
      }
      if (debug) {
        TGraph* gr1 = new TGraph(4, calocell[id].Z, calocell[id].Y);
        gr1->Draw("f");
      }
    }
  }
  geo->cd(path_endcapL_template);

  for (int j = 0; j < nLay_ec; j++) {
    for (int k = 0; k < nCel_ec; k++) {
      int id = k + 100 * j + 1000 * 40;

      calocell[id].id = id;

      dummyLoc_ec[0][0] = rmax / 45. * k - rmax;
      dummyLoc_ec[0][1] = 0.;
      dummyLoc_ec[0][2] = zlevel[j];

      dummyLoc_ec[1][0] = rmax / 45. * k - rmax;
      dummyLoc_ec[1][1] = 0.;
      dummyLoc_ec[1][2] = zlevel[j + 1];

      dummyLoc_ec[2][0] = rmax / 45. * (k + 1) - rmax;
      dummyLoc_ec[2][1] = 0.;
      dummyLoc_ec[2][2] = zlevel[j + 1];

      dummyLoc_ec[3][0] = rmax / 45. * (k + 1) - rmax;
      dummyLoc_ec[3][1] = 0.;
      dummyLoc_ec[3][2] = zlevel[j];

      for (int mm = 0; mm < 4; mm++) {
        geo->LocalToMaster(dummyLoc_ec[mm], dummyMas);

        if (debug) {
          std::cout << "local : " << dummyLoc_ec[mm][0] << " "
                    << dummyLoc_ec[mm][1] << " " << dummyLoc_ec[mm][2]
                    << std::endl;
          std::cout << "master: " << dummyMas[0] << " " << dummyMas[1] << " "
                    << dummyMas[2] << std::endl;
        }

        calocell[id].Y[mm] = dummyMas[0];
        calocell[id].Z[mm] = dummyMas[2];
      }
      if (debug) {
        TGraph* gr1 = new TGraph(4, calocell[id].Z, calocell[id].Y);
        gr1->Draw("f");
      }
    }
  }

  TGeoVolume* GRAIN_vol = geo->FindVolumeFast(GRAIN_vol_name);
  if(!GRAIN_vol) GRAIN_vol = geo->FindVolumeFast(GRIAN_vol_name);

  TGeoEltu* grain = (TGeoEltu*)GRAIN_vol->GetShape();

  GRAIN_dz = grain->GetRmin();
  GRAIN_dy = grain->GetRmax();
  GRAIN_dx = grain->GetDz();

  if(geo->cd(path_GRAIN) == kFALSE) geo->cd(path_GRIAN);

  dummyLoc[0] = 0.;
  dummyLoc[1] = 0.;
  dummyLoc[2] = 0.;

  geo->LocalToMaster(dummyLoc, centerGRAIN);

  initialized = true;
}

//Add trajectories to a graph
void showTrj(int index, int frame1, int frame2, bool over)
{
  std::string same = "l";
  if (over)
    same = "same l";
  t->GetEntry(index);

  for (unsigned int i = 0; i < ev->Trajectories.size(); i++) {
    TGraph* tr_zy = new TGraph(ev->Trajectories[i].Points.size());
    TGraph* tr_zx = new TGraph(ev->Trajectories[i].Points.size());

    for (unsigned int j = 0; j < ev->Trajectories[i].Points.size(); j++) {
      tr_zy->SetPoint(j, ev->Trajectories[i].Points[j].GetPosition().Z(),
                      ev->Trajectories[i].Points[j].GetPosition().Y());
      tr_zx->SetPoint(j, ev->Trajectories[i].Points[j].GetPosition().Z(),
                      ev->Trajectories[i].Points[j].GetPosition().X());
    }

    switch (ev->Trajectories[i].GetPDGCode()) {
      // photons
      case 22:
        tr_zy->SetLineStyle(7);
        tr_zx->SetLineStyle(7);
      // e+/e-
      case 11:
      case -11:
        tr_zy->SetLineColor(kRed);
        tr_zx->SetLineColor(kRed);
        break;

      // mu+/mu-
      case 13:
      case -13:
        tr_zy->SetLineColor(kBlue);
        tr_zx->SetLineColor(kBlue);
        break;

      // proton
      case 2212:
        tr_zy->SetLineColor(kBlack);
        tr_zx->SetLineColor(kBlack);
        break;

      // neutron
      case 2112:
        tr_zy->SetLineStyle(7);
        tr_zx->SetLineStyle(7);
        tr_zy->SetLineColor(kGray);
        tr_zx->SetLineColor(kGray);
        break;

      // pion0
      case 111:
        tr_zy->SetLineStyle(7);
        tr_zx->SetLineStyle(7);
        tr_zy->SetLineColor(kMagenta);
        tr_zx->SetLineColor(kMagenta);
        break;

      // pion+/pion-
      case 211:
      case -211:;
        tr_zy->SetLineColor(kCyan);
        tr_zx->SetLineColor(kCyan);
        break;

      default:
        tr_zy->SetLineColor(8);
        tr_zx->SetLineColor(8);
        break;
    }

    cev->cd(frame1);
    tr_zy->Draw(same.c_str());
    cev->cd(frame2);
    tr_zx->Draw(same.c_str());
  }
}

//Add ede to a graph
void showEde(int index, int frame1, int frame2, bool over)
{
  std::string same = "";
  if (over)
    same = "same";
  t->GetEntry(index);

  //for (auto det : {"Straw", "EMCalSci", "LArHit","DriftVolume"}) {
  for (auto det : {"EMCalSci", "LArHit","DriftVolume"}) {
    for (auto& h : ev->SegmentDetectors[det]) {
      TLine* lzx =
          new TLine(h.Start.Z(), h.Start.X(), h.Stop.Z(), h.Stop.X());
      TLine* lzy =
          new TLine(h.Start.Z(), h.Start.Y(), h.Stop.Z(), h.Stop.Y());
      cev->cd(frame1);
      lzy->Draw(same.c_str());
      cev->cd(frame2);
      lzx->Draw(same.c_str());
    }
  }
}

void show(int index)
{
  if (!initialized) {
    std::cout << "not initialized" << std::endl;
    return;
  }

  if (cev == 0) {
    cev = new TCanvas("cev", TString::Format("Event: %d", index).Data(), 1200,
                      1100);
    cev->Divide(2, 2);
  } else {
    cev->SetTitle(TString::Format("Event: %d", index).Data());
  }

  auto hframe = cev->cd(1)->DrawFrame(centerKLOE[2] - dwz, centerKLOE[1] - dwy,
                                      centerKLOE[2] + dwz, centerKLOE[1] + dwy,
                                      "ZY (side);[mm]; [mm]");
  hframe->GetXaxis()->SetNdivisions(505);

  hframe = cev->cd(2)->DrawFrame(centerKLOE[2] - dwz, centerKLOE[0] - dwx,
                                 centerKLOE[2] + dwz, centerKLOE[0] + dwx,
                                 "XZ (top); [mm]; [mm]");
  hframe->GetXaxis()->SetNdivisions(505);

  cev->cd(2);
  TBox* kloe_int_xz =
      new TBox(centerKLOE[2] - kloe_int_R, centerKLOE[0] - kloe_int_dx,
               centerKLOE[2] + kloe_int_R, centerKLOE[0] + kloe_int_dx);
  kloe_int_xz->SetFillStyle(0);
  kloe_int_xz->Draw();

  TBox* grain_xz =
      new TBox(centerGRAIN[2] - GRAIN_dz, centerGRAIN[0] - GRAIN_dx,
               centerGRAIN[2] + GRAIN_dz, centerGRAIN[0] + GRAIN_dx);
  grain_xz->SetFillStyle(0);
  grain_xz->Draw();

  cev->cd(1);
  TEllipse* grain_zy =
      new TEllipse(centerGRAIN[2], centerGRAIN[1], GRAIN_dz, GRAIN_dy);
  grain_zy->SetFillStyle(0);
  grain_zy->Draw();

  t->GetEntry(index);

  for (std::map<int, gcell>::iterator it = calocell.begin();
       it != calocell.end(); ++it) {
    it->second.adc = 0.;
    it->second.tdc = 0.;
  }

  for (std::map<int, gcell>::iterator it = calocell.begin();
       it != calocell.end(); ++it) {
    if (it->first < 0) continue;

    TGraph* gr = new TGraph(4, it->second.Z, it->second.Y);

    gr->SetFillColor(17);
    if (it->first < 25000)
      cev->cd(1);
    else
      cev->cd(2);
    gr->Draw("f");
  }

  if (showtrj) {
    showTrj(index, 1, 2, false);
  }

  if (showede) {
    showEde(index, 1, 2, false);
  }

  if (showdig) {
    for (unsigned int i = 0; i < vec_tube->size(); i++) {
      if (vec_tube->at(i).hor) {
        TMarker* m = new TMarker(vec_tube->at(i).z, vec_tube->at(i).y, 6);
        cev->cd(1);
        m->Draw();
      } else {
        TMarker* m = new TMarker(vec_tube->at(i).z, vec_tube->at(i).x, 6);
        cev->cd(2);
        m->Draw();
      }
    }

    for (unsigned int i = 0; i < vec_wire->size(); i++) {
      if (vec_wire->at(i).hor) {
        TMarker* m = new TMarker(vec_wire->at(i).z, vec_wire->at(i).y, 6);
        cev->cd(1);
        m->Draw();
      } else {
        TMarker* m = new TMarker(vec_wire->at(i).z, vec_wire->at(i).x, 6);
        cev->cd(2);
        m->Draw();
      }
    }

    for (unsigned int j = 0; j < vec_cell->size(); j++) {
      int id = vec_cell->at(j).id;

      TGraph* gr = new TGraph(4, calocell[id].Z, calocell[id].Y);

      gr->SetFillColor(kBlack);
      if (id < 25000)
        cev->cd(1);
      else
        cev->cd(2);
      gr->Draw("f");
    }
  }

  if (showrec) {
    for (unsigned int i = 0; i < vec_tr->size(); i++) {
      if (vec_tr->at(i).ret_cr == 0 && vec_tr->at(i).ret_ln == 0) {
        cev->cd(1);

        double minz = std::max(vec_tr->at(i).zc - vec_tr->at(i).r,
                               vec_tr->at(i).clY.front().z);
        double maxz = std::min(vec_tr->at(i).zc + vec_tr->at(i).r,
                               vec_tr->at(i).clY.back().z);

        TF1* ffy = new TF1("", "[0]+[1]*TMath::Sqrt([2]*[2]-(x-[3])*(x-[3]))",
                           minz, maxz);
        ffy->SetParameter(0, vec_tr->at(i).yc);
        ffy->SetParameter(1, vec_tr->at(i).ysig);
        ffy->SetParameter(2, vec_tr->at(i).r);
        ffy->SetParameter(3, vec_tr->at(i).zc);
        ffy->Draw("same");

        /*
        TEllipse* e =
            new TEllipse(vec_tr->at(i).zc, vec_tr->at(i).yc, vec_tr->at(i).r);
        e->SetFillStyle(0);
        e->Draw();*/

        cev->cd(2);

        double x0 = vec_tr->at(i).clX.front().x;
        double z0 = vec_tr->at(i).clX.front().z;

        double radq;
        if (abs(z0 - vec_tr->at(i).zc) > vec_tr->at(i).r)
          radq = 0;
        else
          radq = TMath::Sqrt(vec_tr->at(i).r * vec_tr->at(i).r -
                             (z0 - vec_tr->at(i).zc) * (z0 - vec_tr->at(i).zc));

        double yexp = vec_tr->at(i).yc + vec_tr->at(i).ysig * radq;

        double phi0 =
            TMath::ATan2(yexp - vec_tr->at(i).yc, z0 - vec_tr->at(i).zc);

        minz = std::max(vec_tr->at(i).zc - vec_tr->at(i).r,
                        vec_tr->at(i).clX.front().z);
        maxz = std::min(vec_tr->at(i).zc + vec_tr->at(i).r,
                        vec_tr->at(i).clX.back().z);

        TF1* ffx = new TF1("",
                           "[0] - "
                           "[1]/"
                           "[2]*(TMath::ATan2(([4]+[7]*TMath::Sqrt([1]*[1]-(x-["
                           "5])*(x-[5]))) - [4],x - [5]) - [6])*[3]",
                           minz, maxz);

        ffx->SetParameter(0, x0);
        ffx->SetParameter(1, vec_tr->at(i).r);
        ffx->SetParameter(2, vec_tr->at(i).h);
        ffx->SetParameter(3, vec_tr->at(i).b);
        ffx->SetParameter(4, vec_tr->at(i).yc);
        ffx->SetParameter(5, vec_tr->at(i).zc);
        ffx->SetParameter(6, phi0);
        ffx->SetParameter(7, vec_tr->at(i).ysig);
        ffx->Draw("same");

        /*
        TLine* l = new TLine(
            vec_tr->at(i).z0, vec_tr->at(i).x0, centerKLOE[2] + dwz,
            vec_tr->at(i).x0 +
                vec_tr->at(i).b * (centerKLOE[2] + dwz - vec_tr->at(i).z0));
        l->Draw();*/
      }
    }

    for (unsigned int i = 0; i < vec_cl->size(); i++) {

      auto lw = 1. + 9. * (vec_cl->at(i).e - 100) / 1000.;
      if (lw < 1.)
        lw = 1.;
      else if (lw > 10.)
        lw = 10.;

      if (isnan(vec_cl->at(i).sz)) {
        TMarker* m1 = new TMarker(vec_cl->at(i).z, vec_cl->at(i).y, 20);
        m1->SetMarkerColor(kBlue);
        m1->SetMarkerSize(lw);

        cev->cd(1);
        m1->Draw();

        TMarker* m2 = new TMarker(vec_cl->at(i).z, vec_cl->at(i).x, 20);
        m2->SetMarkerColor(kBlue);
        m2->SetMarkerSize(lw);

        cev->cd(2);
        m2->Draw();
      } else {
        TArrow* arr1 = new TArrow(vec_cl->at(i).z - vec_cl->at(i).sz * 0.5 * dt,
                                  vec_cl->at(i).y - vec_cl->at(i).sy * 0.5 * dt,
                                  vec_cl->at(i).z + vec_cl->at(i).sz * 0.5 * dt,
                                  vec_cl->at(i).y + vec_cl->at(i).sy * 0.5 * dt,
                                  0.01, ">");

        arr1->SetLineWidth(lw);
        arr1->SetLineColor(kBlue);
        arr1->SetFillColor(kBlue);
        cev->cd(1);
        arr1->Draw();

        TArrow* arr2 = new TArrow(vec_cl->at(i).z - vec_cl->at(i).sz * 0.5 * dt,
                                  vec_cl->at(i).x - vec_cl->at(i).sx * 0.5 * dt,
                                  vec_cl->at(i).z + vec_cl->at(i).sz * 0.5 * dt,
                                  vec_cl->at(i).x + vec_cl->at(i).sx * 0.5 * dt,
                                  0.01, ">");

        arr2->SetLineWidth(lw);
        arr2->SetLineColor(kBlue);
        arr2->SetFillColor(kBlue);
        cev->cd(2);
        arr2->Draw();
      }
    }
  }
}

void showPri(int index)
{
  if (!initialized) {
    std::cout << "not initialized" << std::endl;
    return;
  }

  if (cpr == 0) {
    cpr = new TCanvas("cpr", TString::Format("Event: %d", index).Data(), 1400,
                      1000);
    cpr->Divide(3, 2);
  } else {
    cpr->SetTitle(TString::Format("Event: %d", index).Data());
  }

  t->GetEntry(index);

  TLorentzVector vpos = ev->Primaries.at(0).GetPosition();

  cpr->cd(1)->DrawFrame(
      -1 + vpos.X(), -1 + vpos.Y(), 1 + vpos.X(), 1 + vpos.Y(),
      TString::Format("XY (front) [Event: %d]", index).Data());
  cpr->cd(2)->DrawFrame(-1 + vpos.Z(), -1 + vpos.Y(), 1 + vpos.Z(),
                        1 + vpos.Y(),
                        TString::Format("ZY (side) [Event: %d]", index).Data());
  cpr->cd(3)->DrawFrame(-1 + vpos.Z(), -1 + vpos.X(), 1 + vpos.Z(),
                        1 + vpos.X(),
                        TString::Format("ZX (top) [Event: %d]", index).Data());

  double maxpXY = 0;
  double maxpXZ = 0;
  double maxpYZ = 0;

  std::vector<TVector3> pmom;

  std::cout << "TRUE "
               "==============================================================="
               "===="
            << std::endl;

  std::cout << std::setw(10) << "PDG"
            << " |" << std::setw(10) << "ID"
            << " |" << std::setw(10) << "PX"
            << " |" << std::setw(10) << "PY"
            << " |" << std::setw(10) << "PZ"
            << " |" << std::setw(10) << "E"
            << " |" << std::endl;

  std::cout << "==============================================================="
               "========="
            << std::endl;

  for (unsigned int i = 0; i < ev->Primaries.at(0).Particles.size(); i++) {
    TVector3 mom;
    mom.SetX(ev->Primaries.at(0).Particles.at(i).GetMomentum().X());
    mom.SetY(ev->Primaries.at(0).Particles.at(i).GetMomentum().Y());
    mom.SetZ(ev->Primaries.at(0).Particles.at(i).GetMomentum().Z());

    double pXY = TMath::Sqrt(mom.X() * mom.X() + mom.Y() * mom.Y());
    double pXZ = TMath::Sqrt(mom.X() * mom.X() + mom.Z() * mom.Z());
    double pYZ = TMath::Sqrt(mom.Y() * mom.Y() + mom.Z() * mom.Z());

    if (pXY > maxpXY) maxpXY = pXY;
    if (pXZ > maxpXZ) maxpXZ = pXZ;
    if (pYZ > maxpYZ) maxpYZ = pYZ;

    pmom.push_back(mom);

    std::cout << std::setw(10)
              << ev->Primaries.at(0).Particles.at(i).GetPDGCode() << " |"
              << std::setw(10)
              << ev->Primaries.at(0).Particles.at(i).GetTrackId() << " |"
              << std::setw(10)
              << ev->Primaries.at(0).Particles.at(i).GetMomentum().X() << " |"
              << std::setw(10)
              << ev->Primaries.at(0).Particles.at(i).GetMomentum().Y() << " |"
              << std::setw(10)
              << ev->Primaries.at(0).Particles.at(i).GetMomentum().Z() << " |"
              << std::setw(10)
              << ev->Primaries.at(0).Particles.at(i).GetMomentum().T() << " |"
              << std::endl;
  }
  std::cout << "==============================================================="
               "========="
            << std::endl;

  cpr->cd(1);
  for (unsigned int i = 0; i < pmom.size(); i++) {
    TArrow* par =
        new TArrow(vpos.X(), vpos.Y(), vpos.X() + pmom.at(i).X() / maxpXY,
                   vpos.Y() + pmom.at(i).Y() / maxpXY, 0.01, "|>");

    switch (ev->Primaries.at(0).Particles.at(i).GetPDGCode()) {
      // photons
      case 22:
        par->SetLineStyle(7);
      // e+/e-
      case 11:
      case -11:
        par->SetLineColor(kRed);
        par->SetFillColor(kRed);
        break;

      // mu+/mu-
      case 13:
      case -13:
        par->SetLineColor(kBlue);
        par->SetFillColor(kBlue);
        break;

      // proton
      case 2212:
        par->SetLineColor(kBlack);
        par->SetFillColor(kBlack);
        break;

      // neutron
      case 2112:
        par->SetLineStyle(7);
        par->SetLineColor(kGray);
        par->SetFillColor(kGray);
        break;

      // pion0
      case 111:
        par->SetLineStyle(7);
        par->SetLineColor(kMagenta);
        par->SetFillColor(kMagenta);
        break;

      // pion+/pion-
      case 211:
      case -211:;
        par->SetLineColor(kCyan);
        par->SetFillColor(kCyan);
        break;

      default:
        par->SetLineColor(8);
        par->SetFillColor(8);
        break;
    }

    par->SetLineWidth(2);
    par->Draw();
  }

  cpr->cd(2);
  for (unsigned int i = 0; i < pmom.size(); i++) {
    TArrow* par =
        new TArrow(vpos.Z(), vpos.Y(), vpos.Z() + pmom.at(i).Z() / maxpYZ,
                   vpos.Y() + pmom.at(i).Y() / maxpYZ, 0.01, "|>");

    switch (ev->Primaries.at(0).Particles.at(i).GetPDGCode()) {
      // photons
      case 22:
        par->SetLineStyle(7);
      // e+/e-
      case 11:
      case -11:
        par->SetLineColor(kRed);
        par->SetFillColor(kRed);
        break;

      // mu+/mu-
      case 13:
      case -13:
        par->SetLineColor(kBlue);
        par->SetFillColor(kBlue);
        break;

      // proton
      case 2212:
        par->SetLineColor(kBlack);
        par->SetFillColor(kBlack);
        break;

      // neutron
      case 2112:
        par->SetLineStyle(7);
        par->SetLineColor(kGray);
        par->SetFillColor(kGray);
        break;

      // pion0
      case 111:
        par->SetLineStyle(7);
        par->SetLineColor(kMagenta);
        par->SetFillColor(kMagenta);
        break;

      // pion+/pion-
      case 211:
      case -211:;
        par->SetLineColor(kCyan);
        par->SetFillColor(kCyan);
        break;

      default:
        par->SetLineColor(8);
        par->SetFillColor(8);
        break;
    }

    par->SetLineWidth(2);
    par->Draw();
  }

  cpr->cd(3);
  for (unsigned int i = 0; i < pmom.size(); i++) {
    TArrow* par =
        new TArrow(vpos.Z(), vpos.X(), vpos.Z() + pmom.at(i).Z() / maxpXZ,
                   vpos.X() + pmom.at(i).X() / maxpXZ, 0.01, "|>");

    switch (ev->Primaries.at(0).Particles.at(i).GetPDGCode()) {
      // photons
      case 22:
        par->SetLineStyle(7);
      // e+/e-
      case 11:
      case -11:
        par->SetLineColor(kRed);
        par->SetFillColor(kRed);
        break;

      // mu+/mu-
      case 13:
      case -13:
        par->SetLineColor(kBlue);
        par->SetFillColor(kBlue);
        break;

      // proton
      case 2212:
        par->SetLineColor(kBlack);
        par->SetFillColor(kBlack);
        break;

      // neutron
      case 2112:
        par->SetLineStyle(7);
        par->SetLineColor(kGray);
        par->SetFillColor(kGray);
        break;

      // pion0
      case 111:
        par->SetLineStyle(7);
        par->SetLineColor(kMagenta);
        par->SetFillColor(kMagenta);
        break;

      // pion+/pion-
      case 211:
      case -211:;
        par->SetLineColor(kCyan);
        par->SetFillColor(kCyan);
        break;

      default:
        par->SetLineColor(8);
        par->SetFillColor(8);
        break;
    }

    par->SetLineWidth(2);
    par->Draw();
  }

  int Nreco = 0;

  for (unsigned int i = 0; i < evt->particles.size(); i++) {
    if (evt->particles.at(i).primary == 1 &&
        evt->particles.at(i).tr.ret_ln == 0 &&
        evt->particles.at(i).tr.ret_cr == 0 &&
        !(evt->particles.at(i).pxreco == 0 &&
          evt->particles.at(i).pyreco == 0 &&
          evt->particles.at(i).pzreco == 0)) {
      Nreco++;
    }
  }

  if (Nreco > 0) {
    double maxpXY_reco = 0;
    double maxpXZ_reco = 0;
    double maxpYZ_reco = 0;

    double xmax = -1E8;
    double xmin = 1E8;
    double ymax = -1E8;
    double ymin = 1E8;
    double zmax = -1E8;
    double zmin = 1E8;

    std::vector<TVector3> pmom_reco;
    std::vector<TVector3> ppos_reco;
    std::vector<int> ppdg_reco;

    std::cout << "RECO "
                 "============================================================="
                 "======"
              << std::endl;

    std::cout << std::setw(10) << "PDG"
              << " |" << std::setw(10) << "ID"
              << " |" << std::setw(10) << "PX"
              << " |" << std::setw(10) << "PY"
              << " |" << std::setw(10) << "PZ"
              << " |" << std::setw(10) << "E"
              << " |" << std::endl;

    std::cout << "============================================================="
                 "==========="
              << std::endl;

    for (unsigned int i = 0; i < evt->particles.size(); i++) {
      if (evt->particles.at(i).primary == 1 &&
          evt->particles.at(i).tr.ret_ln == 0 &&
          evt->particles.at(i).tr.ret_cr == 0 &&
          !(evt->particles.at(i).pxreco == 0 &&
            evt->particles.at(i).pyreco == 0 &&
            evt->particles.at(i).pzreco == 0)) {
        TVector3 mom;
        mom.SetX(evt->particles.at(i).pxreco);
        mom.SetY(evt->particles.at(i).pyreco);
        mom.SetZ(evt->particles.at(i).pzreco);

        TVector3 pos;
        pos.SetX(evt->particles.at(i).xreco);
        pos.SetY(evt->particles.at(i).yreco);
        pos.SetZ(evt->particles.at(i).zreco);

        double pXY = TMath::Sqrt(mom.X() * mom.X() + mom.Y() * mom.Y());
        double pXZ = TMath::Sqrt(mom.X() * mom.X() + mom.Z() * mom.Z());
        double pYZ = TMath::Sqrt(mom.Y() * mom.Y() + mom.Z() * mom.Z());

        if (pXY > maxpXY_reco) maxpXY_reco = pXY;
        if (pXZ > maxpXZ_reco) maxpXZ_reco = pXZ;
        if (pYZ > maxpYZ_reco) maxpYZ_reco = pYZ;

        if (xmax < evt->particles.at(i).xreco)
          xmax = evt->particles.at(i).xreco;
        if (xmin > evt->particles.at(i).xreco)
          xmin = evt->particles.at(i).xreco;
        if (ymax < evt->particles.at(i).yreco)
          ymax = evt->particles.at(i).yreco;
        if (ymin > evt->particles.at(i).yreco)
          ymin = evt->particles.at(i).yreco;
        if (zmax < evt->particles.at(i).zreco)
          zmax = evt->particles.at(i).zreco;
        if (zmin > evt->particles.at(i).zreco)
          zmin = evt->particles.at(i).zreco;

        pmom_reco.push_back(mom);
        ppos_reco.push_back(pos);
        ppdg_reco.push_back(evt->particles.at(i).pdg);

        std::cout << std::setw(10) << evt->particles.at(i).pdg << " |"
                  << std::setw(10) << evt->particles.at(i).tid << " |"
                  << std::setw(10) << evt->particles.at(i).pxreco << " |"
                  << std::setw(10) << evt->particles.at(i).pyreco << " |"
                  << std::setw(10) << evt->particles.at(i).pzreco << " |"
                  << std::setw(10) << evt->particles.at(i).Ereco << " |"
                  << std::endl;
      }
    }
    std::cout << "============================================================="
                 "==========="
              << std::endl;

    double frac = 0.5;
    double dx = std::max<double>(10., xmax - xmin);
    double dy = std::max<double>(10., ymax - ymin);
    double dz = std::max<double>(10., zmax - zmin);

    double dXYmax = std::max<double>(dx, dy);
    double dXZmax = std::max<double>(dx, dz);
    double dYZmax = std::max<double>(dy, dz);

    double xc = 0.5 * (xmax + xmin);
    double yc = 0.5 * (ymax + ymin);
    double zc = 0.5 * (zmax + zmin);

    cpr->cd(4)->DrawFrame(
        xc - (1. + frac) * dXYmax, yc - (1. + frac) * dXYmax,
        xc + (1. + frac) * dXYmax, yc + (1. + frac) * dXYmax,
        TString::Format("XY (front) [Event: %d]", index).Data());

    cpr->cd(5)->DrawFrame(
        zc - (1. + frac) * dYZmax, yc - (1. + frac) * dYZmax,
        zc + (1. + frac) * dYZmax, yc + (1. + frac) * dYZmax,
        TString::Format("ZY (side) [Event: %d]", index).Data());

    cpr->cd(6)->DrawFrame(
        zc - (1. + frac) * dXZmax, xc - (1. + frac) * dXZmax,
        zc + (1. + frac) * dXZmax, xc + (1. + frac) * dXZmax,
        TString::Format("ZX (top) [Event: %d]", index).Data());

    cpr->cd(4);
    for (unsigned int i = 0; i < pmom_reco.size(); i++) {
      TArrow* par =
          new TArrow(ppos_reco.at(i).X(), ppos_reco.at(i).Y(),
                     ppos_reco.at(i).X() +
                         pmom_reco.at(i).X() / maxpXY_reco * dXYmax * frac,
                     ppos_reco.at(i).Y() +
                         pmom_reco.at(i).Y() / maxpXY_reco * dXYmax * frac,
                     0.01, "|>");

      switch (ppdg_reco.at(i)) {
        // photons
        case 22:
          par->SetLineStyle(7);
        // e+/e-
        case 11:
        case -11:
          par->SetLineColor(kRed);
          par->SetFillColor(kRed);
          break;

        // mu+/mu-
        case 13:
        case -13:
          par->SetLineColor(kBlue);
          par->SetFillColor(kBlue);
          break;

        // proton
        case 2212:
          par->SetLineColor(kBlack);
          par->SetFillColor(kBlack);
          break;

        // neutron
        case 2112:
          par->SetLineStyle(7);
          par->SetLineColor(kGray);
          par->SetFillColor(kGray);
          break;

        // pion0
        case 111:
          par->SetLineStyle(7);
          par->SetLineColor(kMagenta);
          par->SetFillColor(kMagenta);
          break;

        // pion+/pion-
        case 211:
        case -211:;
          par->SetLineColor(kCyan);
          par->SetFillColor(kCyan);
          break;

        default:
          par->SetLineColor(8);
          par->SetFillColor(8);
          break;
      }

      par->SetLineWidth(2);
      par->Draw();
    }

    cpr->cd(5);
    for (unsigned int i = 0; i < pmom_reco.size(); i++) {
      TArrow* par =
          new TArrow(ppos_reco.at(i).Z(), ppos_reco.at(i).Y(),
                     ppos_reco.at(i).Z() +
                         pmom_reco.at(i).Z() / maxpYZ_reco * dYZmax * frac,
                     ppos_reco.at(i).Y() +
                         pmom_reco.at(i).Y() / maxpYZ_reco * dYZmax * frac,
                     0.01, "|>");

      switch (ppdg_reco.at(i)) {
        // photons
        case 22:
          par->SetLineStyle(7);
        // e+/e-
        case 11:
        case -11:
          par->SetLineColor(kRed);
          par->SetFillColor(kRed);
          break;

        // mu+/mu-
        case 13:
        case -13:
          par->SetLineColor(kBlue);
          par->SetFillColor(kBlue);
          break;

        // proton
        case 2212:
          par->SetLineColor(kBlack);
          par->SetFillColor(kBlack);
          break;

        // neutron
        case 2112:
          par->SetLineStyle(7);
          par->SetLineColor(kGray);
          par->SetFillColor(kGray);
          break;

        // pion0
        case 111:
          par->SetLineStyle(7);
          par->SetLineColor(kMagenta);
          par->SetFillColor(kMagenta);
          break;

        // pion+/pion-
        case 211:
        case -211:;
          par->SetLineColor(kCyan);
          par->SetFillColor(kCyan);
          break;

        default:
          par->SetLineColor(8);
          par->SetFillColor(8);
          break;
      }

      par->SetLineWidth(2);
      par->Draw();
    }

    cpr->cd(6);
    for (unsigned int i = 0; i < pmom_reco.size(); i++) {
      TArrow* par =
          new TArrow(ppos_reco.at(i).Z(), ppos_reco.at(i).X(),
                     ppos_reco.at(i).Z() +
                         pmom_reco.at(i).Z() / maxpXZ_reco * dXZmax * frac,
                     ppos_reco.at(i).X() +
                         pmom_reco.at(i).X() / maxpXZ_reco * dXZmax * frac,
                     0.01, "|>");

      switch (ppdg_reco.at(i)) {
        // photons
        case 22:
          par->SetLineStyle(7);
        // e+/e-
        case 11:
        case -11:
          par->SetLineColor(kRed);
          par->SetFillColor(kRed);
          break;

        // mu+/mu-
        case 13:
        case -13:
          par->SetLineColor(kBlue);
          par->SetFillColor(kBlue);
          break;

        // proton
        case 2212:
          par->SetLineColor(kBlack);
          par->SetFillColor(kBlack);
          break;

        // neutron
        case 2112:
          par->SetLineStyle(7);
          par->SetLineColor(kGray);
          par->SetFillColor(kGray);
          break;

        // pion0
        case 111:
          par->SetLineStyle(7);
          par->SetLineColor(kMagenta);
          par->SetFillColor(kMagenta);
          break;

        // pion+/pion-
        case 211:
        case -211:;
          par->SetLineColor(kCyan);
          par->SetFillColor(kCyan);
          break;

        default:
          par->SetLineColor(8);
          par->SetFillColor(8);
          break;
      }

      par->SetLineWidth(2);
      par->Draw();
    }
  }
}

void DumpPri(int nev, int* ids)
{
  gROOT->SetBatch(true);
  for (int i = 0; i < nev; i++) {
    showPri(ids[i]);

    if (i == 0)
      cpr->SaveAs("display.pdf(");
    else if (i == nev - 1)
      cpr->SaveAs("display.pdf)");
    else
      cpr->SaveAs("display.pdf");
  }
  gROOT->SetBatch(false);
}

void DumpPri(int nev = 100, int istart = 0)
{
  int* ids = new int[nev];
  for (int i = 0; i < nev; i++) ids[i] = istart + i;

  DumpPri(nev, ids);
}

//To compare edepsim to reco
void comparison()
{
  if (!tDigit || !tReco) {
    std::cerr << "tDigit or tReco TTree non found" << std::endl;
    for (auto f : vf)
      f->Close();
    return;
  }

  //Loop on tDigit entries
  Long64_t nEntries = tDigit->GetEntries();

  TH1* hist1;
  TH1* hist2;

  for (Long64_t i=0; i<nEntries; i++)
  {
    show(i);

    //3rd plot
    cev->cd(3);
    tDigit->Draw(Form("fired_wires.y:fired_wires.z >> htemp%lld", i), Form("Entry$ == %lld && fired_wires.hor==1", i), "P");
    hist1 = gDirectory->Get<TH1>(Form("htemp%lld", i));
    hist1->SetStats(kFALSE);
    hist1->SetMarkerStyle(20);
    hist1->SetMarkerSize(0.5);
    hist1->Draw();

    //4th plot
    cev->cd(4);
    tDigit->Draw(Form("fired_wires.x:fired_wires.z >> htemp%lld_a", i),Form("Entry$ == %lld && fired_wires.hor==0", i), "P");
    hist2 = gDirectory->Get<TH1>(Form("htemp%lld_a", i));
    hist2->SetStats(kFALSE);
    hist2->SetMarkerStyle(20);
    hist2->SetMarkerSize(0.5);
    hist2->Draw();

    if (showtrj2reco)
      showTrj(i, 3, 4, true);
    if (showede2reco)
      showEde(i, 3, 4, true);

    if (is_out_file_set == true) {
      if (i == 0)
          cev->Print(fout + "(",Form("Title: Event N: %lld", i));
      else if (i == nEntries-1)
          cev->Print(fout + ")",Form("Title: Event N: %lld", i));
      else
      {
          cev->Print(fout,Form("Title: Event N: %lld", i));
      }

      std::cout << " writing entry N: "<< i << std::endl;
    }

    delete hist1;
    delete hist2;
    delete gDirectory->Get<TH1>(Form("htemp%lld", i));
    delete gDirectory->Get<TH1>(Form("htemp%lld_a", i));
  }
}

void help()
{
  std::cout << "Display -e <event number> -mc <MC file>"
               "[-f <input file1> -f <input file2> ... ] [-o <output file>] "
               "[--batch] [options]\n\n"
               "--compare      -- to wirte a file over all events\n"
               "--trj          -- to show trajectories\n"
               "--ede          -- to show energy deposits\n"
               "--dgt          -- to show digits\n"
               "--rec          -- to show reco objects\n"
            << std::endl;
}

int main(int argc, char* argv[])
{
  TApplication* myapp = new TApplication("myapp", 0, 0);

  int evid = 0;
  int index = 1;

  while (index < argc) {
    TString opt = argv[index];
    if (opt.CompareTo("-e") == 0) {
      try {
        evid = atoi(argv[++index]);
        is_ev_number_set = true;
      } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
      }
    } else if (opt.CompareTo("-mc") == 0) {
      try {
        fmc = new TFile(argv[++index]);
        is_mc_file_set = true;
      } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
      }
    } else if (opt.CompareTo("-f") == 0) {
      try {
        vf.push_back(new TFile(argv[++index]));
      } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
      }
    } else if (opt.CompareTo("-o") == 0) {
      try {
        fout = argv[++index];
        is_out_file_set = true;
      } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
      }
    } else if (opt.CompareTo("--batch") == 0) {
      try {
        is_batch_mode_set = true;
      } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
      }
    }
    else if (opt.CompareTo("--trj") == 0)
      showtrj = true;
    else if (opt.CompareTo("--trj2reco") == 0)
      showtrj2reco = true;
    else if (opt.CompareTo("--ede") == 0)
      showede = true;
    else if (opt.CompareTo("--ede2reco") == 0)
      showede2reco = true;
    else if (opt.CompareTo("--dgt") == 0)
      showdig = true;
    else if (opt.CompareTo("--rec") == 0)
      showrec = true;
    else if (opt.CompareTo("--compare") == 0)
      compare = true;
    index++;
  }

  if ((is_ev_number_set == false && compare == false) || is_mc_file_set == false) {
    help();
    return 1;
  }

  if (is_batch_mode_set == true) {
    gROOT->SetBatch();
  }

  init();

  if (compare == true) {
    comparison();
  }
  else
  {
    show(evid);

    if (is_out_file_set == true) {
      cev->SaveAs(fout.Data());
    }

    if (is_batch_mode_set == false) {
      myapp->Run();
    }
  }

  return 0;
}