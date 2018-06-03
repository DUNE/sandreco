#include <TFile.h>
#include <TTree.h>
#include <TGeoManager.h>
#include <TGeoTrd2.h>
#include <TGraph.h>
#include <TCanvas.h>

#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4Event.h"
#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4HitSegment.h"

#include <iostream>

const bool debug = true;

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

cell calocell[nTotCells];

TFile* fev;
TTree* tev;
TG4Event* ev;
TGeoManager* geo;

double CellLocalX[nCellModule][5];
double CellLocalZ[nCellModule][5];

void init(const char* fname)
{
  fev = new TFile(fname);
  tev = (TTree*) fev->Get("Events");
  tev->SetBranchAddress("Event",&ev);
  
  geo = (TGeoManager*) fev->Get("EDepSimGeometry");
  
  int index[nLay+1] = {0, 40, 80, 120, 160, 209};
  double dx1[nLay];
  double dx2[nLay];
  double dz1[nLay];
  double dz2[nLay];
  
  TGeoTrd2* scin = (TGeoTrd2*) geo->FindVolumeFast("KLOEBarrelECAL_0_sci_slab_0_volume_PV")->GetShape();
  TGeoTrd2* lead = (TGeoTrd2*) geo->FindVolumeFast("KLOEBarrelECAL_0_lead_slab_0_volume_PV")->GetShape();
  
  TGeoTrd2* mod = (TGeoTrd2*) geo->FindVolumeFast("KLOEBarrelECAL_0_volume_PV")->GetShape();
  double dzmod = mod->GetDz();
    
  double dzs = scin->GetDz();
  double dzl = lead->GetDz();
  
  for(int i = 0; i < nLay; i++)
  {
    TGeoTrd2* shape1 = (TGeoTrd2*) geo->FindVolumeFast(TString::Format("KLOEBarrelECAL_0_sci_slab_%i_volume_PV", index[i]).Data())->GetShape();
    TGeoTrd2* shape2 = (TGeoTrd2*) geo->FindVolumeFast(TString::Format("KLOEBarrelECAL_0_lead_slab_%i_volume_PV", index[i+1]-1).Data())->GetShape();
    
    if(debug)
      std::cout << index[i] << " " << shape1 << " " << index[i+1]-1 << " " << shape2 << std::endl;
    
    dx1[i] = shape1->GetDx1();
    dx2[i] = shape2->GetDx2();    
    
    dz1[i] = index[i] * (dzs + dzl);
    dz2[i] = index[i+1] * (dzs + dzl);
  }
  
  dz2[nLay-1] = dzmod;
  
  if(debug)
  {
    for(int i = 0; i < nLay; i++)
    {  
      std::cout << dx1[i] << " " << dx2[i] << " " << dz1[i] << " " << dz2[i] << std::endl;
    }
  }
  
  TCanvas* c1 = new TCanvas();
  c1->DrawFrame(-300,-150,300,150);
  
  for(int i = 0; i < nLay; i++)
  {
    for(int j = 0; j < nCel; j++)
    {
      CellLocalX[i*nCel+j][0] =  -dx1[i] + 2 * dx1[i]/12. * j;
      CellLocalX[i*nCel+j][1] =  -dx1[i] + 2 * dx1[i]/12. * (j+1);
      CellLocalX[i*nCel+j][2] =  -dx2[i] + 2 * dx2[i]/12. * (j+1);
      CellLocalX[i*nCel+j][3] =  -dx2[i] + 2 * dx2[i]/12. * j;
      CellLocalX[i*nCel+j][4] =  -dx1[i] + 2 * dx1[i]/12. * j;
            
      CellLocalZ[i*nCel+j][0] =  -dzmod + 2 * dz1[i];
      CellLocalZ[i*nCel+j][1] =  -dzmod + 2 * dz1[i];
      CellLocalZ[i*nCel+j][2] =  -dzmod + 2 * dz2[i];
      CellLocalZ[i*nCel+j][3] =  -dzmod + 2 * dz2[i];
      CellLocalZ[i*nCel+j][4] =  -dzmod + 2 * dz1[i];
      
      if(debug)
        std::cout << CellLocalZ[i*nCel+j][0] << " " << CellLocalX[i*nCel+j][0] << " " << 
                    CellLocalZ[i*nCel+j][1] << " " << CellLocalX[i*nCel+j][1] << " " << 
                    CellLocalZ[i*nCel+j][2] << " " << CellLocalX[i*nCel+j][2] << " " << 
                    CellLocalZ[i*nCel+j][3] << " " << CellLocalX[i*nCel+j][3] << std::endl;
                    
    
      TGraph* gr1 = new TGraph(5, CellLocalX[i*nCel+j], CellLocalZ[i*nCel+j]);
      gr1->Draw("l");
      //gr1->Draw("f");
    }
  }
  
} 