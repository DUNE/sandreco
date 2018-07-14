#include <TFile.h>
#include <TTree.h>
#include <TGeoManager.h>
#include <TGeoTrd2.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>

#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4Event.h"
#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4HitSegment.h"

#include <iostream>
#include <map>

const bool debug = false;

class cell {
  public:
    int id;
    double Z[4];
    double Y[4];
    double adc;
    double tdc;
};

static const int nMod = 24;
static const int nLay = 5;
static const int nCel = 12;

static const int nTotCells = nMod * nLay * nCel; 
static const int nCellModule = nLay * nCel;

std::map<int, cell> calocell;

TFile* fev;
TTree* tev;
TG4Event* ev;
TGeoManager* geo;
TCanvas* cev;
      
std::map<int, std::vector<double> >* list_pe;
std::map<int, double>* adc;
std::map<int, double>* tdc;

double centerKLOE[3];
double CellLocalX[nCellModule][4];
double CellLocalZ[nCellModule][4];

int palette = 87;

void init(const char* fname)
{
  gStyle->SetPalette(palette);

  fev = new TFile(fname);
  tev = (TTree*) fev->Get("Events");
  tev->SetBranchAddress("Event",&ev);
  
  tev->SetBranchAddress("cellPE",&list_pe);
  tev->SetBranchAddress("cellADC",&adc);
  tev->SetBranchAddress("cellTDC",&tdc);
  
  double dummyLoc[3];
  double dummyMas[3];
  
  geo = (TGeoManager*) fev->Get("EDepSimGeometry");
  
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
    cev = new TCanvas("cev","", 700, 700);
    cev->DrawFrame(centerKLOE[2] - 2500,centerKLOE[1] - 2500,centerKLOE[2] + 2500,centerKLOE[1] + 2500);
  }
  
  for(int i = 0; i < nMod; i++)
  {
    geo->cd(TString::Format(path_template,i).Data());
    
    if(debug)
      std::cout << "node: " << i << " " << geo->GetCurrentNode() << " " << geo->GetCurrentNode()->GetName() << " " << TString::Format(path_template,i).Data() << std::endl;
    
    for(int j = 0; j < nLay; j++)
    {
      for(int k = 0; k < nCel; k++)
      {
        int index = i * (nLay * nCel) + j * (nCel) + k;
        int id = k + 100 * j + 1000 * i;
        
        int local_index = j*nCel+k;
        
        calocell[id].id = id;
        
        if(debug)
          std::cout << i << " " << j << " " << k << " " << index << " " << id << " " << local_index << " " << nMod << " " << nLay << " " << nCel << std::endl;
        
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
}

void show(int index)
{
  cev = new TCanvas("cev",TString::Format("event: %d",index).Data(), 700, 700);
  cev->DrawFrame(centerKLOE[2] - 2500,centerKLOE[1] - 2500,centerKLOE[2] + 2500,centerKLOE[1] + 2500);
  
  tev->GetEntry(index);
  
  for(std::map<int, cell>::iterator it=calocell.begin(); it != calocell.end(); ++it)
  {
    it->second.adc = 0.;
    it->second.tdc = 0.;
  }
  
  double max = 0.0;
  
  for(std::map<int, double>::iterator it=adc->begin(); it != adc->end(); ++it)
  {
    calocell[it->first].adc = it->second;
    if(it->second > max) max = it->second;
  }
  
  for(std::map<int, double>::iterator it=tdc->begin(); it != tdc->end(); ++it)
  {
    calocell[it->first].tdc = it->second;
  }
  
  for(std::map<int, cell>::iterator it=calocell.begin(); it != calocell.end(); ++it)
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
    
    if((it->second.adc + calocell[-1*it->first].adc) == 0.)
    {
      gr->SetFillColor(19);
    }
    else
    {
      std::cout << "ID: " << it->first  << "\tADC1: " << setw(3) << it->second.adc 
                                        << "\tTDC1: " << setw(3) << calocell[it->first].tdc 
                                        << "\tADC2: " << setw(3) << calocell[-1*it->first].adc 
                                        << "\tTDC2: " << setw(3) << calocell[-1*it->first].tdc << std::endl;
    }
    
    gr->Draw("f");
  }
  
  for(unsigned int i = 0; i < ev->Trajectories.size(); i++)
  {
    TGraph* tr = new TGraph(ev->Trajectories[i].Points.size());
    for(unsigned int j = 0; j < ev->Trajectories[i].Points.size(); j++)
    {
      tr->SetPoint(j, ev->Trajectories[i].Points[j].Position.Z(),ev->Trajectories[i].Points[j].Position.Y());
    }
    
    switch(ev->Trajectories[i].PDGCode)
    {
      // photons
      case 22:
        tr->SetLineStyle(7);
      // e+/e-
      case 11:
      case -11:
        tr->SetLineColor(kRed);
      break;
      
      // mu+/mu-
      case 13:
      case -13:
        tr->SetLineColor(kBlue);
      break;
      
      // proton
      case 2212:
        tr->SetLineColor(kBlack);
      break;
      
      // neutron
      case 2112:
        tr->SetLineStyle(7);
        tr->SetLineColor(kGray);
      break;
      
      // pion0
      case 111:
        tr->SetLineStyle(7);
        tr->SetLineColor(kMagenta);
      break;
      
      // pion+/pion- 
      case 211:
      case -211:;
        tr->SetLineColor(kCyan);
      break;
      
      default:
        tr->SetLineColor(8);
      break;        
    }
    
    tr->Draw("l");
  }
}