#include <TApplication.h>
#include <TArrow.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoTrd2.h>
#include <TGeoTube.h>
#include <TGraph.h>
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

TFile* f = 0;
TFile* fmc = 0;
TTree* t = 0;
TG4Event* ev = new TG4Event;
event* evt = new event;
TGeoManager* geo = 0;
TCanvas* cev = 0;
TCanvas* cpr = 0;

std::vector<dg_cell>* vec_cell = new std::vector<dg_cell>;
std::vector<dg_tube>* vec_digi = new std::vector<dg_tube>;
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

const char* barrel_mod_vol_name = "ECAL_lv_PV";
const char* endcap_mod_vol_name = "ECAL_end_lv_PV";
}  // namespace display

using namespace display;

void init(const char* mcfile, const char* ifile)
{
  gStyle->SetPalette(palette);

  fmc = new TFile(mcfile);
  f = new TFile(ifile);
  TTree* tEvent = reinterpret_cast<TTree*>(f->Get("tEvent"));
  TTree* tReco = reinterpret_cast<TTree*>(f->Get("tReco"));
  TTree* tDigit = reinterpret_cast<TTree*>(f->Get("tDigit"));
  TTree* tEdep = reinterpret_cast<TTree*>(fmc->Get("EDepSimEvents"));
  TTree* tGenie =
      reinterpret_cast<TTree*>(fmc->Get("DetSimPassThru/gRooTracker"));

  if (!tEdep) return;

  if (tReco) tEdep->AddFriend(tReco);
  if (tDigit) tEdep->AddFriend(tDigit);
  if (tEvent) tEdep->AddFriend(tEvent);
  if (tGenie) tEdep->AddFriend(tGenie);

  tEdep->SetBranchAddress("Event", &ev);
  if (tDigit) tDigit->SetBranchAddress("dg_cell", &vec_cell);
  if (tDigit) tDigit->SetBranchAddress("dg_tube", &vec_digi);
  if (tReco) tReco->SetBranchAddress("track", &vec_tr);
  if (tReco) tReco->SetBranchAddress("cluster", &vec_cl);
  if (tEvent) tEvent->SetBranchAddress("event", &evt);

  t = tEdep;

  geo = reinterpret_cast<TGeoManager*>(fmc->Get("EDepSimGeometry"));

  if (!geo) return;

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

  initialized = true;
}

void show(int index, bool showtrj = true, bool showfit = true,
          bool showdig = true)
{
  if (!initialized) {
    std::cout << "not initialized" << std::endl;
    return;
  }

  if (cev == 0) {
    cev = new TCanvas("cev", TString::Format("Event: %d", index).Data(), 1200,
                      600);
    cev->Divide(2, 1);
  } else {
    cev->SetTitle(TString::Format("Event: %d", index).Data());
  }

  cev->cd(1)->DrawFrame(centerKLOE[2] - dwz, centerKLOE[1] - dwy,
                        centerKLOE[2] + dwz, centerKLOE[1] + dwy,
                        "ZY (side); (mm); (mm)");

  cev->cd(2)->DrawFrame(centerKLOE[2] - dwz, centerKLOE[0] - dwx,
                        centerKLOE[2] + dwz, centerKLOE[0] + dwx,
                        "XZ (top); (mm); (mm)");

  cev->cd(2);
  TBox* kloe_int_xz =
      new TBox(centerKLOE[2] - kloe_int_R, centerKLOE[0] - kloe_int_dx,
               centerKLOE[2] + kloe_int_R, centerKLOE[0] + kloe_int_dx);
  kloe_int_xz->SetFillStyle(0);
  kloe_int_xz->Draw();

  t->GetEntry(index);

  for (std::map<int, gcell>::iterator it = calocell.begin();
       it != calocell.end(); ++it) {
    it->second.adc = 0.;
    it->second.tdc = 0.;
  }
  /*
  double max = 0.0;

  for(unsigned int i = 0; i < vec_cell->size(); i++)
  {
    calocell[vec_cell->at(i).id].adc = vec_cell->at(i).adc1;
    calocell[vec_cell->at(i).id].tdc = vec_cell->at(i).tdc1;
    calocell[-1*vec_cell->at(i).id].adc = vec_cell->at(i).adc2;
    calocell[-1*vec_cell->at(i).id].tdc = vec_cell->at(i).tdc2;

    if((vec_cell->at(i).adc1 + vec_cell->at(i).adc2) > max)
      max = (vec_cell->at(i).adc1 + vec_cell->at(i).adc2);
  }

  for(std::map<int, gcell>::iterator it=calocell.begin(); it != calocell.end();
  ++it)
  {
    if(it->first < 0)
      continue;

    TGraph* gr = new TGraph(4, it->second.Z, it->second.Y);

    int ndiv = max > 255 ? 255 : max;
    double scale = ndiv / max;
    int color = 0.01 + (it->second.adc + calocell[-1*it->first].adc) * scale;
    int ncolors = gStyle->GetNumberOfColors();
    int theColor = (color + 0.99) * double(ncolors) / double(ndiv);
    int thiscolor = gStyle->GetColorPalette(theColor);

    gr->SetFillColor(thiscolor);

    cev->cd(1);

    if((it->second.adc + calocell[-1*it->first].adc) == 0.)
    {
      gr->SetFillColor(19);
    }
    else
    {

      //std::cout << "ID: " << it->first  << "\tADC1: " << setw(3) <<
  it->second.adc
      //                                  << "\tTDC1: " << setw(3) <<
  calocell[it->first].tdc
      //                                  << "\tADC2: " << setw(3) <<
  calocell[-1*it->first].adc
      //                                  << "\tTDC2: " << setw(3) <<
  calocell[-1*it->first].tdc << std::endl;
    }

    gr->Draw("f");
  }
  */

  for (std::map<int, gcell>::iterator it = calocell.begin();
       it != calocell.end(); ++it) {
    if (it->first < 0) continue;

    TGraph* gr = new TGraph(4, it->second.Z, it->second.Y);

    gr->SetFillColor(19);
    if (it->first < 25000)
      cev->cd(1);
    else
      cev->cd(2);
    gr->Draw("f");
  }

  if (showtrj) {
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

      cev->cd(1);
      tr_zy->Draw("l");
      cev->cd(2);
      tr_zx->Draw("l");
    }
  }

  for (unsigned int j = 0; j < vec_cl->size(); j++) {
    for (unsigned int i = 0; i < vec_cl->at(j).cells.size(); i++) {
      int id = vec_cl->at(j).cells.at(i).id;

      calocell[id].adc = vec_cl->at(j).cells.at(i).ps1.at(0).adc;
      calocell[id].tdc = vec_cl->at(j).cells.at(i).ps1.at(0).tdc;
      calocell[-id].adc = vec_cl->at(j).cells.at(i).ps2.at(0).adc;
      calocell[-id].tdc = vec_cl->at(j).cells.at(i).ps2.at(0).tdc;

      if (showdig) {
        TGraph* gr = new TGraph(4, calocell[id].Z, calocell[id].Y);
        int color = (vec_cl->at(j).tid == 0) ? 632 : vec_cl->at(j).tid;
        gr->SetFillColor(color);
        if (id < 25000)
          cev->cd(1);
        else
          cev->cd(2);
        gr->Draw("f");
      }
    }
  }

  if (showdig) {
    for (unsigned int i = 0; i < vec_digi->size(); i++) {
      if (vec_digi->at(i).hor) {
        TMarker* m = new TMarker(vec_digi->at(i).z, vec_digi->at(i).y, 6);
        cev->cd(1);
        m->Draw();
      } else {
        TMarker* m = new TMarker(vec_digi->at(i).z, vec_digi->at(i).x, 6);
        cev->cd(2);
        m->Draw();
      }
    }
  }

  if (showfit) {
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
      int color = (vec_cl->at(i).tid == 0) ? 632 : vec_cl->at(i).tid;

      TMarker* m1 = new TMarker(vec_cl->at(i).z, vec_cl->at(i).y, 34);
      // m1->SetMarkerColor(color);
      cev->cd(1);
      m1->Draw();

      TMarker* m2 = new TMarker(vec_cl->at(i).z, vec_cl->at(i).x, 34);
      // m2->SetMarkerColor(color);
      cev->cd(2);
      m2->Draw();

      TArrow* arr1 =
          new TArrow(vec_cl->at(i).z - vec_cl->at(i).sz * 0.5 * dt,
                     vec_cl->at(i).y - vec_cl->at(i).sy * 0.5 * dt,
                     vec_cl->at(i).z + vec_cl->at(i).sz * 0.5 * dt,
                     vec_cl->at(i).y + vec_cl->at(i).sy * 0.5 * dt, 0.01, ">");
      cev->cd(1);
      arr1->Draw();

      TArrow* arr2 =
          new TArrow(vec_cl->at(i).z - vec_cl->at(i).sz * 0.5 * dt,
                     vec_cl->at(i).x - vec_cl->at(i).sx * 0.5 * dt,
                     vec_cl->at(i).z + vec_cl->at(i).sz * 0.5 * dt,
                     vec_cl->at(i).x + vec_cl->at(i).sx * 0.5 * dt, 0.01, ">");
      cev->cd(2);
      arr2->Draw();
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

int main(int argc, char* argv[])
{
  TApplication* myapp = new TApplication("myapp", 0, 0);

  bool showtrj = true;
  bool showfit = true;
  bool showdig = true;

  int evid = 0;
  TString fname;
  TString fmc;
  TString tmp;

  if (argc < 4) {
    std::cout
        << "Display <event number> <MC file> <input file> [show trajectories] "
           "[show fits] [show digits]"
        << std::endl;
    return 1;
  } else {
    evid = atoi(argv[1]);
    fmc = argv[2];
    fname = argv[3];

    if (argc > 4) {
      tmp = argv[4];
      if (tmp.CompareTo("false") == 0) showtrj = false;
    }

    if (argc > 5) {
      tmp = argv[5];
      if (tmp.CompareTo("false") == 0) showfit = false;
    }

    if (argc > 6) {
      tmp = argv[6];
      if (tmp.CompareTo("false") == 0) showdig = false;
    }
  }

  init(fmc.Data(), fname.Data());

  show(evid, showtrj, showfit, showdig);

  myapp->Run();

  return 0;
}
