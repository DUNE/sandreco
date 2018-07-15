#include <TFile.h>
#include <TTree.h>
#include <TGeoManager.h>
#include <TGeoTrd2.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TChain.h>
#include <TMarker.h>
#include <TEllipse.h>

#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4Event.h"
#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4HitSegment.h"

#include <iostream>
#include <map>

#include "/mnt/nas01/users/mtenti/wd/analysis/KLOEcal/loader/loader.C"

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
}

std::vector<cell>* vec_cell;
std::vector<digit>* vec_digi;
std::vector<track>* vec_tr;
std::map<int, gcell> calocell;

TChain* t = 0;
TG4Event* ev = new TG4Event;
TGeoManager* geo = 0;
TCanvas* cev = 0;

void init(const char* fTrueMC, const char* fCalDigi, const char* fSttDigi, const char* fTrack)
{
  gStyle->SetPalette(ns_draw::palette);
  
  TChain* tTrueMC = new TChain("EDepSimEvents","EDepSimEvents");
  tTrueMC->Add(fTrueMC);
  TChain* tCalDigi = new TChain("tCalDigi","tCalDigi");
  tCalDigi->Add(fCalDigi);
  TChain* tSttDigi = new TChain("tSttDigi","tSttDigi");
  tSttDigi->Add(fSttDigi);
  TChain* tTrack = new TChain("tTrack","tTrack");
  tTrack->Add(fTrack);
  tTrueMC->SetBranchAddress("Event",&ev);
  tCalDigi->SetBranchAddress("cell",&vec_cell);
  tSttDigi->SetBranchAddress("Stt",&vec_digi);
  tTrack->SetBranchAddress("track",&vec_tr);
  
  t = tTrueMC;
  t->AddFriend(tCalDigi);
  t->AddFriend(tSttDigi);
  t->AddFriend(tTrack);
  
  TFile f(tTrueMC->GetListOfFiles()->At(0)->GetTitle());
  geo = (TGeoManager*) f.Get("EDepSimGeometry");
  
  double dummyLoc[3];
  double dummyMas[3];
  
  geo->cd("volWorld_PV/volDetEnclosure_PV_0/volKLOEFULLECALSENSITIVE_EXTTRK_NEWGAP_PV_0/KLOETrackingRegion_volume_PV_0");
  
  dummyLoc[0] = 0.;
  dummyLoc[1] = 0.;
  dummyLoc[2] = 0.;
  geo->LocalToMaster(dummyLoc, ns_draw::centerKLOE);
  
  double dzlay[ns_draw::nLay+1] = {115, 115-22, 115-22-22, 115-22-22-22, 115-22-22-22-22, 115-22-22-22-22-27};
  double dx1[ns_draw::nLay];
  double dx2[ns_draw::nLay];
  
  TGeoTrd2* mod = (TGeoTrd2*) geo->FindVolumeFast("KLOEBarrelECAL_0_volume_PV")->GetShape();
  
  double xmax = mod->GetDx1();
  double xmin = mod->GetDx2();
  double dz = mod->GetDz();
  
  if(ns_draw::debug)
  {
    std::cout << dz << " " << xmax << " " << xmin << std::endl;
  }
  
  for(int i = 0; i < ns_draw::nLay; i++)
  {    
    dx1[i] = xmax - (xmax - xmin)/dz * dzlay[i];
    dx2[i] = xmax - (xmax - xmin)/dz * dzlay[i+1];
  }
  
  if(ns_draw::debug)
  {
    for(int i = 0; i < ns_draw::nLay; i++)
    {  
      std::cout << dx1[i] << " " << dx2[i] << " " << dzlay[i] << " " << dzlay[i+1] << std::endl;
    }
  }
  
  for(int i = 0; i < ns_draw::nLay; i++)
  {
    for(int j = 0; j < ns_draw::nCel; j++)
    {
      // from bottom-left to top-right 
    
      ns_draw::CellLocalX[i*ns_draw::nCel+j][0] =  -dx1[i] + 2 * dx1[i]/12. * j;
      ns_draw::CellLocalX[i*ns_draw::nCel+j][1] =  -dx1[i] + 2 * dx1[i]/12. * (j+1);
      ns_draw::CellLocalX[i*ns_draw::nCel+j][2] =  -dx2[i] + 2 * dx2[i]/12. * (j+1);
      ns_draw::CellLocalX[i*ns_draw::nCel+j][3] =  -dx2[i] + 2 * dx2[i]/12. * j;
            
      ns_draw::CellLocalZ[i*ns_draw::nCel+j][0] =  -dz + 2 * dzlay[i];
      ns_draw::CellLocalZ[i*ns_draw::nCel+j][1] =  -dz + 2 * dzlay[i];
      ns_draw::CellLocalZ[i*ns_draw::nCel+j][2] =  -dz + 2 * dzlay[i+1];
      ns_draw::CellLocalZ[i*ns_draw::nCel+j][3] =  -dz + 2 * dzlay[i+1];
      
      if(ns_draw::debug)
        std::cout << ns_draw::CellLocalZ[i*ns_draw::nCel+j][0] << " " << ns_draw::CellLocalX[i*ns_draw::nCel+j][0] << " " << 
                     ns_draw::CellLocalZ[i*ns_draw::nCel+j][1] << " " << ns_draw::CellLocalX[i*ns_draw::nCel+j][1] << " " << 
                     ns_draw::CellLocalZ[i*ns_draw::nCel+j][2] << " " << ns_draw::CellLocalX[i*ns_draw::nCel+j][2] << " " << 
                     ns_draw::CellLocalZ[i*ns_draw::nCel+j][3] << " " << ns_draw::CellLocalX[i*ns_draw::nCel+j][3] << std::endl;
    }
  }
  
  const char* path_template = "volWorld_PV/volDetEnclosure_PV_0/volKLOEFULLECALSENSITIVE_EXTTRK_NEWGAP_PV_0/KLOEBarrelECAL_%d_volume_PV_0";
  
  double CellMasterY[ns_draw::nCellModule][4];
  double CellMasterZ[ns_draw::nCellModule][4];
  
  if(ns_draw::debug)
  {
    cev = new TCanvas("cev","", 700, 700);
    cev->DrawFrame(ns_draw::centerKLOE[2] - 2500,ns_draw::centerKLOE[1] - 2500,ns_draw::centerKLOE[2] + 2500,ns_draw::centerKLOE[1] + 2500);
  }
  
  for(int i = 0; i < ns_draw::nMod; i++)
  {
    geo->cd(TString::Format(path_template,i).Data());
    
    if(ns_draw::debug)
      std::cout << "node: " << i << " " << geo->GetCurrentNode() << " " << geo->GetCurrentNode()->GetName() << " " << TString::Format(path_template,i).Data() << std::endl;
    
    for(int j = 0; j < ns_draw::nLay; j++)
    {
      for(int k = 0; k < ns_draw::nCel; k++)
      {
        int index = i * (ns_draw::nLay * ns_draw::nCel) + j * (ns_draw::nCel) + k;
        int id = k + 100 * j + 1000 * i;
        
        int local_index = j*ns_draw::nCel+k;
        
        calocell[id].id = id;
        
        if(ns_draw::debug)
          std::cout << i << " " << j << " " << k << " " << index << " " << id << " " << local_index << " " << ns_draw::nMod << " " << ns_draw::nLay << " " << ns_draw::nCel << std::endl;
        
        for(int m = 0; m < 4; m++)
        {
          dummyLoc[0] = ns_draw::CellLocalX[local_index][m];
          dummyLoc[1] = 0.;
          dummyLoc[2] = ns_draw::CellLocalZ[local_index][m];
          
          geo->LocalToMaster(dummyLoc, dummyMas);
          
          if(ns_draw::debug)
          {
            std::cout << "local : " << dummyLoc[0] << " " << dummyLoc[1] << " " << dummyLoc[2] << std::endl;
            std::cout << "master: " << dummyMas[0] << " " << dummyMas[1] << " " << dummyMas[2] << std::endl;
          } 
          
          calocell[id].Y[m] = dummyMas[1];
          calocell[id].Z[m] = dummyMas[2];
          
          CellMasterY[local_index][m] = dummyMas[1];
          CellMasterZ[local_index][m] = dummyMas[2];
        }                    
        
        if(ns_draw::debug)
        {
          TGraph* gr1 = new TGraph(4, CellMasterZ[local_index], CellMasterY[local_index]);
          gr1->Draw("f");
        }
      }
    }
  }
  
  ns_draw::initialized = true;
}

void show(int index)
{
  if(!ns_draw::initialized)
  {
    std::cout << "not initialized" << std::endl;
    return;
  }


  if(cev == 0 )
  {
    cev = new TCanvas("cev",TString::Format("event: %d",index).Data(), 1200, 600);
    cev->Divide(2,1);
  }
  
  cev->cd(1)->DrawFrame(ns_draw::centerKLOE[2] - ns_draw::dwz,
                 ns_draw::centerKLOE[1] - ns_draw::dwy,
                 ns_draw::centerKLOE[2] + ns_draw::dwz,
                 ns_draw::centerKLOE[1] + ns_draw::dwy);
  
  cev->cd(2)->DrawFrame(ns_draw::centerKLOE[2] - ns_draw::dwz,
                 ns_draw::centerKLOE[0] - ns_draw::dwx,
                 ns_draw::centerKLOE[2] + ns_draw::dwz,
                 ns_draw::centerKLOE[0] + ns_draw::dwx);
  
  t->GetEntry(index);
  
  for(std::map<int, gcell>::iterator it=calocell.begin(); it != calocell.end(); ++it)
  {
    it->second.adc = 0.;
    it->second.tdc = 0.;
  }
  
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
  
  for(std::map<int, gcell>::iterator it=calocell.begin(); it != calocell.end(); ++it)
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
      /*
      std::cout << "ID: " << it->first  << "\tADC1: " << setw(3) << it->second.adc 
                                        << "\tTDC1: " << setw(3) << calocell[it->first].tdc 
                                        << "\tADC2: " << setw(3) << calocell[-1*it->first].adc 
                                        << "\tTDC2: " << setw(3) << calocell[-1*it->first].tdc << std::endl;
      */
    }
    
    gr->Draw("f");
  }
  
  for(unsigned int i = 0; i < vec_digi->size(); i++)
  {
    if(vec_digi->at(i).hor)
    {
      TMarker* m = new TMarker(vec_digi->at(i).z,vec_digi->at(i).y,6);
      cev->cd(1);
      m->Draw();
    }
    else
    {
      TMarker* m = new TMarker(vec_digi->at(i).z,vec_digi->at(i).x,6);
      cev->cd(2);
      m->Draw();
    }
  }
  
  for(unsigned int i = 0; i < vec_tr->size(); i++)
  {
    cev->cd(1);
    TEllipse* e = new TEllipse(vec_tr->at(i).zc, vec_tr->at(i).yc, vec_tr->at(i).r);
    e->SetFillStyle(0);
    e->Draw();
    
    cev->cd(2);
    TLine* l = new TLine(vec_tr->at(i).z0, vec_tr->at(i).x0, 
                         ns_draw::centerKLOE[2] + ns_draw::dwz, 
                         vec_tr->at(i).x0 + vec_tr->at(i).b * (ns_draw::centerKLOE[2] + ns_draw::dwz - vec_tr->at(i).z0));
    l->Draw();
  }
  
  for(unsigned int i = 0; i < ev->Trajectories.size(); i++)
  {
    TGraph* tr_zy = new TGraph(ev->Trajectories[i].Points.size());
    TGraph* tr_zx = new TGraph(ev->Trajectories[i].Points.size());
    
    for(unsigned int j = 0; j < ev->Trajectories[i].Points.size(); j++)
    {
      tr_zy->SetPoint(j, ev->Trajectories[i].Points[j].Position.Z(),ev->Trajectories[i].Points[j].Position.Y());
      tr_zx->SetPoint(j, ev->Trajectories[i].Points[j].Position.Z(),ev->Trajectories[i].Points[j].Position.X());
    }
    
    switch(ev->Trajectories[i].PDGCode)
    {
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