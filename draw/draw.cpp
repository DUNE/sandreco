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
#include <TBox.h>
#include <TGeoTube.h>
#include <TArrow.h>

#include <iostream>
#include <map>

#include "/wd/dune-it/enurec/analysis/kloe-simu/loader/loader.C"

#include "/wd/sw/EDEPSIM/edep-sim.binary/include/EDepSim/TG4Event.h"
#include "/wd/sw/EDEPSIM/edep-sim.binary/include/EDepSim/TG4HitSegment.h"

namespace ns_draw {
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
  TTree* t = 0;
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

using namespace ns_draw;

void init(const char* ifile)
{
  gStyle->SetPalette(palette);

  f = new TFile(ifile);
  TTree* tEvent = reinterpret_cast<TTree*>(f->Get("tEvent"));
  TTree* tReco = reinterpret_cast<TTree*>(f->Get("tReco"));
  TTree* tDigit = reinterpret_cast<TTree*>(f->Get("tDigit"));
  TTree* tEdep = reinterpret_cast<TTree*>(f->Get("EDepSimEvents"));
  TTree* tGenie = reinterpret_cast<TTree*>(f->Get("gRooTracker"));

  tEvent->AddFriend(tReco);
  tEvent->AddFriend(tDigit);
  tEvent->AddFriend(tEdep);
  tEvent->AddFriend(tGenie);
  
  tEdep->SetBranchAddress("Event",&ev);
  tDigit->SetBranchAddress("cell",&vec_cell);
  tDigit->SetBranchAddress("Stt",&vec_digi);
  tReco->SetBranchAddress("track",&vec_tr);
  tReco->SetBranchAddress("cluster",&vec_cl);

  t = tEvent;

  geo = reinterpret_cast<TGeoManager*>(f->Get("EDepSimGeometry"));
  
  double dummyLoc[3];
  double dummyMas[3];
  
  geo->cd("volWorld_PV/volDetEnclosure_PV_0/volKLOEFULLECALSENSITIVE_EXTTRK_NEWGAP_PV_0/KLOETrackingRegion_volume_PV_0");
  
  dummyLoc[0] = 0.;
  dummyLoc[1] = 0.;
  dummyLoc[2] = 0.;
  geo->LocalToMaster(dummyLoc, centerKLOE);
  
  double dzlay[nLay+1] = {115, 115-22, 115-22-22, 115-22-22-22, 115-22-22-22-22, 115-22-22-22-22-27};
  double dx1[nLay];
  double dx2[nLay];
  
  TGeoTrd2* mod = (TGeoTrd2*) geo->FindVolumeFast("KLOEBarrelECAL_0_volume_PV")->GetShape();
  
  double xmax = mod->GetDx1();
  double xmin = mod->GetDx2();
  double dz = mod->GetDz();
  
  if(debug)
  {
    std::cout << dz << " " << xmax << " " << xmin << std::endl;
  }
  
  for(int i = 0; i < nLay; i++)
  {    
    dx1[i] = xmax - (xmax - xmin)/dz * dzlay[i];
    dx2[i] = xmax - (xmax - xmin)/dz * dzlay[i+1];
  }
  
  if(debug)
  {
    for(int i = 0; i < nLay; i++)
    {  
      std::cout << dx1[i] << " " << dx2[i] << " " << dzlay[i] << " " << dzlay[i+1] << std::endl;
    }
  }
  
  for(int i = 0; i < nLay; i++)
  {
    for(int j = 0; j < nCel; j++)
    {
      // from bottom-left to top-right 
    
      CellLocalX[i*nCel+j][0] =  -dx1[i] + 2 * dx1[i]/12. * j;
      CellLocalX[i*nCel+j][1] =  -dx1[i] + 2 * dx1[i]/12. * (j+1);
      CellLocalX[i*nCel+j][2] =  -dx2[i] + 2 * dx2[i]/12. * (j+1);
      CellLocalX[i*nCel+j][3] =  -dx2[i] + 2 * dx2[i]/12. * j;
            
      CellLocalZ[i*nCel+j][0] =  -dz + 2 * dzlay[i];
      CellLocalZ[i*nCel+j][1] =  -dz + 2 * dzlay[i];
      CellLocalZ[i*nCel+j][2] =  -dz + 2 * dzlay[i+1];
      CellLocalZ[i*nCel+j][3] =  -dz + 2 * dzlay[i+1];
      
      if(debug)
        std::cout << CellLocalZ[i*nCel+j][0] << " " << CellLocalX[i*nCel+j][0] << " " << 
                     CellLocalZ[i*nCel+j][1] << " " << CellLocalX[i*nCel+j][1] << " " << 
                     CellLocalZ[i*nCel+j][2] << " " << CellLocalX[i*nCel+j][2] << " " << 
                     CellLocalZ[i*nCel+j][3] << " " << CellLocalX[i*nCel+j][3] << std::endl;
    }
  }
  
  const char* path_template = "volWorld_PV/volDetEnclosure_PV_0/volKLOEFULLECALSENSITIVE_EXTTRK_NEWGAP_PV_0/KLOEBarrelECAL_%d_volume_PV_0";
  
  double CellMasterY[nCellModule][4];
  double CellMasterZ[nCellModule][4];
  
  if(debug)
  {
    TCanvas* cd1 = new TCanvas("cd1","", 700, 700);
    cd1->DrawFrame(centerKLOE[2] - 2500,
                            centerKLOE[1] - 2500,
                            centerKLOE[2] + 2500,
                            centerKLOE[1] + 2500);
  }
  
  for(int i = 0; i < nMod; i++)
  {
    geo->cd(TString::Format(path_template,i).Data());
    
    if(debug)
      std::cout << "node: " << i << " " << geo->GetCurrentNode() 
                                 << " " << geo->GetCurrentNode()->GetName() 
                                 << " " << TString::Format(path_template,i).Data() << std::endl;
    
    for(int j = 0; j < nLay; j++)
    {
      for(int k = 0; k < nCel; k++)
      {
        
        int index = i * (nLay * nCel) + j * (nCel) + k;
        int id = k + 100 * j + 1000 * i;
        
        int local_index = j*nCel+k;
        
        calocell[id].id = id;
        
        if(debug)
          std::cout << i << " " 
                    << j << " " 
                    << k << " " 
                    << index << " " 
                    << id << " " 
                    << local_index << " " 
                    << nMod << " " 
                    << nLay << " " 
                    << nCel << std::endl;
        
        for(int m = 0; m < 4; m++)
        {
          dummyLoc[0] = CellLocalX[local_index][m];
          dummyLoc[1] = 0.;
          dummyLoc[2] = CellLocalZ[local_index][m];
          
          geo->LocalToMaster(dummyLoc, dummyMas);
          
          if(debug)
          {
            std::cout << "local : " << dummyLoc[0] << " " << dummyLoc[1] << " " << dummyLoc[2] << std::endl;
            std::cout << "master: " << dummyMas[0] << " " << dummyMas[1] << " " << dummyMas[2] << std::endl;
          } 
          
          calocell[id].Y[m] = dummyMas[1];
          calocell[id].Z[m] = dummyMas[2];
          
          CellMasterY[local_index][m] = dummyMas[1];
          CellMasterZ[local_index][m] = dummyMas[2];
        }                 
        
        if(debug)
        {
          TGraph* gr1 = new TGraph(4, CellMasterZ[local_index], CellMasterY[local_index]);
          gr1->Draw("f");
        }
      }
    }
  }
  
  if(debug)
  {
    TCanvas* cd2 = new TCanvas("cd2","", 700, 700);
    cd2->DrawFrame(centerKLOE[2] - 2500,
                            centerKLOE[0] - 2500,
                            centerKLOE[2] + 2500,
                            centerKLOE[0] + 2500);
  }
  
  TGeoTube* ec = (TGeoTube*) geo->FindVolumeFast("KLOEEndcapECALL_volume_PV")->GetShape();
  
  double rmax = ec->GetRmax();
  double rmin = ec->GetRmin();
  double dz_ec = ec->GetDz();
  
  double dummyLoc_ec[4][3];
  
  const char* path_endcapR_template = "volWorld_PV/volDetEnclosure_PV_0/volKLOEFULLECALSENSITIVE_EXTTRK_NEWGAP_PV_0/KLOEEndcapECALR_volume_PV_0";
  
  geo->cd(path_endcapR_template);
  
  for(int j = 0; j < nLay_ec; j++)
  {
    for(int k = 0; k < nCel_ec; k++)
    {
      int id = k + 100 * j + 1000 * 30;
      
      calocell[id].id = id;
      
      dummyLoc_ec[0][0] = rmax / 45. * k - rmax;
      dummyLoc_ec[0][1] = 0.;
      dummyLoc_ec[0][2] = dz - 2 * dzlay[j];
      
      dummyLoc_ec[1][0] = rmax / 45. * k - rmax;
      dummyLoc_ec[1][1] = 0.;
      dummyLoc_ec[1][2] = dz - 2 * dzlay[j+1];
      
      dummyLoc_ec[2][0] = rmax / 45. * (k + 1) - rmax;
      dummyLoc_ec[2][1] = 0.;
      dummyLoc_ec[2][2] = dz - 2 * dzlay[j+1];
      
      dummyLoc_ec[3][0] = rmax / 45. * (k + 1) - rmax;
      dummyLoc_ec[3][1] = 0.;
      dummyLoc_ec[3][2] = dz - 2 * dzlay[j];
      
      for(int m = 0; m < 4; m++)
      {
        geo->LocalToMaster(dummyLoc_ec[m], dummyMas);
          
        if(debug)
        {
          std::cout << "local : " << dummyLoc_ec[m][0] << " " << dummyLoc_ec[m][1] << " " << dummyLoc_ec[m][2] << std::endl;
          std::cout << "master: " << dummyMas[0] << " " << dummyMas[1] << " " << dummyMas[2] << std::endl;
        } 
        
        calocell[id].Y[m] = dummyMas[0];
        calocell[id].Z[m] = dummyMas[2];
      }
      if(debug)
      {
        TGraph* gr1 = new TGraph(4, calocell[id].Z, calocell[id].Y);
        gr1->Draw("f");
      }
    }
  }
  
  const char* path_endcapL_template = "volWorld_PV/volDetEnclosure_PV_0/volKLOEFULLECALSENSITIVE_EXTTRK_NEWGAP_PV_0/KLOEEndcapECALL_volume_PV_0";
  geo->cd(path_endcapL_template);
  
  for(int j = 0; j < nLay_ec; j++)
  {
    for(int k = 0; k < nCel_ec; k++)
    {
      int id = k + 100 * j + 1000 * 40;
      
      calocell[id].id = id;
      
      dummyLoc_ec[0][0] = rmax / 45. * k - rmax;
      dummyLoc_ec[0][1] = 0.;
      dummyLoc_ec[0][2] = -dz + 2 * dzlay[j];
      
      dummyLoc_ec[1][0] = rmax / 45. * k - rmax;
      dummyLoc_ec[1][1] = 0.;
      dummyLoc_ec[1][2] = -dz + 2 * dzlay[j+1];
      
      dummyLoc_ec[2][0] = rmax / 45. * (k + 1) - rmax;
      dummyLoc_ec[2][1] = 0.;
      dummyLoc_ec[2][2] = -dz + 2 * dzlay[j+1];
      
      dummyLoc_ec[3][0] = rmax / 45. * (k + 1) - rmax;
      dummyLoc_ec[3][1] = 0.;
      dummyLoc_ec[3][2] = -dz + 2 * dzlay[j];
      
      for(int m = 0; m < 4; m++)
      {        
        geo->LocalToMaster(dummyLoc_ec[m], dummyMas);
          
        if(debug)
        {
          std::cout << "local : " << dummyLoc_ec[m][0] << " " << dummyLoc_ec[m][1] << " " << dummyLoc_ec[m][2] << std::endl;
          std::cout << "master: " << dummyMas[0] << " " << dummyMas[1] << " " << dummyMas[2] << std::endl;
        }
        
        calocell[id].Y[m] = dummyMas[0];
        calocell[id].Z[m] = dummyMas[2];
      }
      if(debug)
      {
        TGraph* gr1 = new TGraph(4, calocell[id].Z, calocell[id].Y);
        gr1->Draw("f");
      }
    }
  }
    
  initialized = true;
}

void show(int index, bool showtrj = true, bool showfit = true, bool showdig = true)
{
  if(!initialized)
  {
    std::cout << "not initialized" << std::endl;
    return;
  }


  if(cev == 0 )
  {
    cev = new TCanvas("cev",TString::Format("Event: %d",index).Data(), 1200, 600);
    cev->Divide(2,1);
  }
  else
  {
    cev->SetTitle(TString::Format("Event: %d",index).Data());
  }
  
  cev->cd(1)->DrawFrame(centerKLOE[2] - dwz,
                 centerKLOE[1] - dwy,
                 centerKLOE[2] + dwz,
                 centerKLOE[1] + dwy,
                 "ZY (side); (mm); (mm)");
  
  cev->cd(2)->DrawFrame(centerKLOE[2] - dwz,
                 centerKLOE[0] - dwx,
                 centerKLOE[2] + dwz,
                 centerKLOE[0] + dwx,
                 "XZ (top); (mm); (mm)");
                 
  
  cev->cd(2);
  TBox* kloe_int_xz = new TBox(centerKLOE[2] - kloe_int_R, centerKLOE[0] - kloe_int_dx, centerKLOE[2] + kloe_int_R, centerKLOE[0] + kloe_int_dx);
  kloe_int_xz->SetFillStyle(0);
  kloe_int_xz->Draw();
  
  t->GetEntry(index);
  
  for(std::map<int, gcell>::iterator it=calocell.begin(); it != calocell.end(); ++it)
  {
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
      
      //std::cout << "ID: " << it->first  << "\tADC1: " << setw(3) << it->second.adc 
      //                                  << "\tTDC1: " << setw(3) << calocell[it->first].tdc 
      //                                  << "\tADC2: " << setw(3) << calocell[-1*it->first].adc 
      //                                  << "\tTDC2: " << setw(3) << calocell[-1*it->first].tdc << std::endl;
    }
    
    gr->Draw("f");
  }
  */
  
  for(std::map<int, gcell>::iterator it=calocell.begin(); it != calocell.end(); ++it)
  {
    if(it->first < 0)
      continue;
  
    TGraph* gr = new TGraph(4, it->second.Z, it->second.Y);  
    
    gr->SetFillColor(19);
    if(it->first < 25000) 
      cev->cd(1);
    else
      cev->cd(2);
    gr->Draw("f");
  }
  
  if(showtrj)
  {
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
  
  for(unsigned int j = 0; j < vec_cl->size(); j++)
  {
    for(unsigned int i = 0; i < vec_cl->at(j).cells.size(); i++)
    {
      int id = vec_cl->at(j).cells.at(i).id;
      
      calocell[id].adc = vec_cl->at(j).cells.at(i).adc1;
      calocell[id].tdc = vec_cl->at(j).cells.at(i).tdc1;
      calocell[-id].adc = vec_cl->at(j).cells.at(i).adc2;
      calocell[-id].tdc = vec_cl->at(j).cells.at(i).tdc2;
      
      if(showdig)
      {
        TGraph* gr = new TGraph(4, calocell[id].Z, calocell[id].Y); 
        int color = (vec_cl->at(j).tid == 0) ? 632 : vec_cl->at(j).tid;
        gr->SetFillColor(color);
        if(id < 25000) 
          cev->cd(1);
        else
          cev->cd(2);
        gr->Draw("f");
      }
    }
  }
  
  if(showdig)
  {
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
  }
  
  if(showfit)
  {
    for(unsigned int i = 0; i < vec_tr->size(); i++)
    {
      if(vec_tr->at(i).ret_cr == 0 && vec_tr->at(i).ret_ln == 0)
      {
        cev->cd(1);
        TEllipse* e = new TEllipse(vec_tr->at(i).zc, vec_tr->at(i).yc, vec_tr->at(i).r);
        e->SetFillStyle(0);
        e->Draw();
        
        cev->cd(2);
        TLine* l = new TLine(vec_tr->at(i).z0, vec_tr->at(i).x0, 
                             centerKLOE[2] + dwz, 
                             vec_tr->at(i).x0 + vec_tr->at(i).b * (centerKLOE[2] + dwz - vec_tr->at(i).z0));
        l->Draw();
      }
    }
    
    for(unsigned int i = 0; i < vec_cl->size(); i++)
    {
      int color = (vec_cl->at(i).tid == 0) ? 632 : vec_cl->at(i).tid;
      
      TMarker* m1 = new TMarker(vec_cl->at(i).z,vec_cl->at(i).y,34);
      //m1->SetMarkerColor(color);
      cev->cd(1);
      m1->Draw();
      
      TMarker* m2 = new TMarker(vec_cl->at(i).z,vec_cl->at(i).x,34);
      //m2->SetMarkerColor(color);
      cev->cd(2);
      m2->Draw();
      
      TArrow* arr1 = new TArrow(vec_cl->at(i).z - vec_cl->at(i).sz * 0.5 * dt, 
         vec_cl->at(i).y - vec_cl->at(i).sy * 0.5 * dt,
         vec_cl->at(i).z + vec_cl->at(i).sz * 0.5 * dt, 
         vec_cl->at(i).y + vec_cl->at(i).sy * 0.5 * dt, 0.01, ">");
      cev->cd(1);
      arr1->Draw();
      
     	TArrow* arr2 = new TArrow(vec_cl->at(i).z - vec_cl->at(i).sz * 0.5 * dt, 
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
  if(!initialized)
  {
    std::cout << "not initialized" << std::endl;
    return;
  }


  if(cpr == 0 )
  {
    cpr = new TCanvas("cpr",TString::Format("Event: %d",index).Data(), 1400, 500);
    cpr->Divide(3,1);
  }
  else
  {
    cpr->SetTitle(TString::Format("Event: %d",index).Data());
  }
  
  cpr->cd(1)->DrawFrame(-1,-1,1,1,"XY (front)");
  cpr->cd(2)->DrawFrame(-1,-1,1,1,"ZY (side)");
  cpr->cd(3)->DrawFrame(-1,-1,1,1,"ZX (top)");
  
  t->GetEntry(index);
  
  double maxpXY = 0;
  double maxpXZ = 0;
  double maxpYZ = 0;
  
  std::vector<TVector3> pmom;
    
  std::cout << "=============================================================" << std::endl;
  
  std::cout << std::setw(10) << "PDG" << " |" <<
    std::setw(10) << "PX" << " |" <<
    std::setw(10) << "PY" << " |" <<
    std::setw(10) << "PZ" << " |" <<
    std::setw(10) << "E"  << " |" << std::endl;
    
  std::cout << "=============================================================" << std::endl;
  
  for(unsigned int i = 0; i < ev->Primaries.at(0).Particles.size(); i++)
  {
    TVector3 mom;
    mom.SetX(ev->Primaries.at(0).Particles.at(i).Momentum.X());
    mom.SetY(ev->Primaries.at(0).Particles.at(i).Momentum.Y());
    mom.SetZ(ev->Primaries.at(0).Particles.at(i).Momentum.Z());
    
    double pXY = TMath::Sqrt(mom.X()*mom.X()+mom.Y()*mom.Y());
    double pXZ = TMath::Sqrt(mom.X()*mom.X()+mom.Z()*mom.Z());
    double pYZ = TMath::Sqrt(mom.Y()*mom.Y()+mom.Z()*mom.Z());
    
    if(pXY > maxpXY)
      maxpXY = pXY;
    if(pXZ > maxpXZ)
      maxpXZ = pXZ;
    if(pYZ > maxpYZ)
      maxpYZ = pYZ;
      
    pmom.push_back(mom);
    
    std::cout << std::setw(10) << ev->Primaries.at(0).Particles.at(i).PDGCode << " |" <<
      std::setw(10) << ev->Primaries.at(0).Particles.at(i).Momentum.X() << " |" <<
      std::setw(10) << ev->Primaries.at(0).Particles.at(i).Momentum.Y() << " |" <<
      std::setw(10) << ev->Primaries.at(0).Particles.at(i).Momentum.Z() << " |" <<
      std::setw(10) << ev->Primaries.at(0).Particles.at(i).Momentum.T() << " |" << std::endl;
  }
  std::cout << "=============================================================" << std::endl;
  
  cpr->cd(1);
  for(unsigned int i = 0; i < pmom.size(); i++)
  {
    TArrow* par = new TArrow(0.,0.,pmom.at(i).X()/maxpXY,pmom.at(i).Y()/maxpXY,0.01,"|>");
    
    switch(ev->Primaries.at(0).Particles.at(i).PDGCode)
    {
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
  for(unsigned int i = 0; i < pmom.size(); i++)
  {
    TArrow* par = new TArrow(0.,0.,pmom.at(i).Z()/maxpYZ,pmom.at(i).Y()/maxpYZ,0.01,"|>");
    
    switch(ev->Primaries.at(0).Particles.at(i).PDGCode)
    {
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
  for(unsigned int i = 0; i < pmom.size(); i++)
  {
    TArrow* par = new TArrow(0.,0.,pmom.at(i).Z()/maxpXZ,pmom.at(i).X()/maxpXZ,0.01,"|>");
    
    switch(ev->Primaries.at(0).Particles.at(i).PDGCode)
    {
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
