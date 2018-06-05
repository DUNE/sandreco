#include <TGeoManager.h>
#include <TString.h>
#include <TGeoNode.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TGeoTrd2.h>
#include <TChain.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TRandom3.h>

#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4Event.h"
#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4HitSegment.h"

#include <vector>
#include <map>
#include <iostream>

// Energy MeV
// Distance mm
// Time ns
 
TRandom3 r;

const bool debug = false;

double Attenuation(double d, int planeID)
{
/*
     dE/dx attenuation - Ea=p1*exp(-d/atl1)+(1.-p1)*exp(-d/atl2)
       d    distance from photocatode - 2 PMTs/cell; d1 and d2 
      atl1  50. cm
      atl2  430 cm planes 1-2    innermost plane is 1 
            380 cm plane 3
            330 cm planes 4-5
       p1   0.35
*/
  const double p1 = 0.35;
  const double alt1 = 50.;
  double alt2 = 0.0;
  
  switch (planeID)
  {
    case 0: 
    case 1: 
      alt2 = 430.0;
    break;
     
    case 2: 
      alt2 = 380.0;
    break;
     
    case 3: 
    case 4: 
      alt2 = 330.0;
    break;
     
    default: 
      //std::cout << "planeID out if range" << std::endl;
      alt2 = -999.0;
    break;
  }
  
  if(debug)
  {
    std::cout << "planeID = " << planeID << std::endl;
    std::cout << "\tp1   = " << p1 << std::endl;
    std::cout << "\talt1 = " << alt1 << std::endl;
    std::cout << "\talt2 = " << alt2 << std::endl;
    std::cout << "\tatt  = " << p1 * TMath::Exp(-d/alt1) + (1.-p1) * TMath::Exp(-d/alt2) << std::endl;
  }
  
  return p1 * TMath::Exp(-d/alt1) + (1.-p1) * TMath::Exp(-d/alt2);
}

double E2PE(double E)
{
  // Average number of photoelectrons = 25*Ea(MeV)
  const double e2p2 = 25.;
  
  if(debug)
    std::cout << "E = " << E << " -> p.e. = " << e2p2*E << std::endl;
  
  return e2p2*E;
}

double petime(double t0, double d)
{
  /*
     - For each photoelectron: Time for TDC simulation obtained from 

C  PHOTOELECTRON TIME :  Particle TIME in the cell  
C                      + SCINTILLATION DECAY TIME + 
C                      + signal propagation to the PMT
C                      + 1ns  uncertainty
     
               TPHE = Part_time+TSDEC+DPM1*VLFB+Gauss(1ns)

      VLFB = 5.85 ns/m 
!!!! Input-TDC Scintillation time - 
               TSDEC = TSCIN*(1./RNDMPH(1)-1)**TSCEX  (ns)

      TSCIN  3.08  ns
      TSCEX  0.588
  */
  
  double tscin = 3.08;
  double tscex = 0.588;
  double vlfb = 5.85;
  
  double mm_to_m = 1E-3;
  
  double tdec = tscin * TMath::Power(1./r.Uniform()-1.,tscex);
  
  double time = t0 + tdec + vlfb * d * mm_to_m + r.Gaus();
  
  if(debug)
    std::cout << "time: " << time << std::endl;
    
  
  return time;
}

bool ProcessHit(TGeoManager* g, const TG4HitSegment& hit, int& modID, int& planeID, int& cellID, double& d1, double& d2, double& t, double& de)
{
  modID = -999;
  planeID = -999;
  cellID = -999;
  d1 = -999;
  d2 = -999;
  t = -999;
  
  double x = 0.5*(hit.Start.X()+hit.Stop.X());
  double y = 0.5*(hit.Start.Y()+hit.Stop.Y());
  double z = 0.5*(hit.Start.Z()+hit.Stop.Z());
  
  t = 0.5*(hit.Start.T()+hit.Stop.T());
  de = hit.EnergyDeposit;
  
  TGeoNode* node = g->FindNode(x,y,z);
  
  if(node == 0) return false;
  
  TString str = node->GetName();
  
  if(str.Contains("KLOEBarrelECAL") == false || str.Contains("lead_slab") == true) return false;
  
  TObjArray* obja = str.Tokenize("_");
  
  int slabID;
  modID  = ((TObjString*) obja->At(1))->GetString().Atoi();
  slabID = ((TObjString*) obja->At(4))->GetString().Atoi();
  
  delete obja;
  
  // planeID==0 -> smallest slab
  // planeID==208 -> biggest slab
  planeID = slabID/40;
  
  if (planeID > 4) planeID = 4;
  
  double Pmaster[3];
  double Plocal[3];
  Pmaster[0] = x;
  Pmaster[1] = y;
  Pmaster[2] = z;
  
  g->GetCurrentNavigator()->MasterToLocal(Pmaster,Plocal);
  
  TGeoTrd2* trd = (TGeoTrd2*) node->GetVolume()->GetShape();
  
  double dx1 = trd->GetDx1();
  double dx2 = trd->GetDx2();
  double dz  = trd->GetDz(); 
  double dy1 = trd->GetDy1(); 
  double dy2 = trd->GetDy2(); 
    
  d1 = dy1 + Plocal[1];
  d2 = dy1 - Plocal[1];
  
  // http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/html/Detector/Geometry/geomSolids.html
  // if z = -dz -> dx = 2*dx1
  // if z =  dz -> dx = 2*dx2
  double dx = - (dx1 - dx2) / dz * Plocal[2];
  dx =+ dx1 + dx2;
  
  // Cell width at z
  double cellw = dx / 12.;
  
  cellID = (Plocal[0] + dx*0.5) / cellw;
  
  if(debug)
  {
    std::cout << "hit" << std::endl;
    std::cout << "\t[x,y,z]                " << x << " " << y << " " << z << std::endl;
    std::cout << "\t[modID,planeID,cellID] " << modID << " " << planeID << " " << cellID << std::endl;
    std::cout << "\t[d1,d2,t,de]           " << d1 << " " << d2 << " " << t << " " << de << std::endl;
  }
  
  return true;
}

void SimulatePE(TG4Event* ev, TGeoManager* g, int index, std::map<int, std::vector<double> >& list_pe)
{    
    int modID, planeID, cellID, id;
    double d1, d2, t0, de;

    for (std::map<std::string,std::vector<TG4HitSegment> >::iterator it=ev->SegmentDetectors.begin();
            it!=ev->SegmentDetectors.end(); ++it)
    {
      if(it->first == "EMCalSci")
      {
        for(unsigned int j = 0; j < it->second.size(); j++)
        {          
          if(ProcessHit(g, it->second[j], modID, planeID, cellID, d1, d2, t0, de) == true)
          {
            double en1 = de * Attenuation(d1, planeID);
            double en2 = de * Attenuation(d2, planeID);
            
            double ave_pe1 = E2PE(en1);
            double ave_pe2 = E2PE(en2);
            
            int pe1 = r.Poisson(ave_pe1);
            int pe2 = r.Poisson(ave_pe2);
            
            id = cellID + 100 * planeID + 1000 * modID;
            
            if(debug)
            {
              std::cout << "cell ID: " << id << std::endl;
              std::cout << "\t" << de << " " << en1 << " " << en2 << std::endl;
              std::cout << "\t" << ave_pe1 << " " << ave_pe2 << std::endl;
              std::cout << "\t" << pe1 << " " << pe2 << std::endl;
            }
            
            //PMT 1 -> x < 0 -> ID > 0
            //PMT 2 -> x > 0 -> ID < 0
            
            for(int i = 0; i < pe1; i++)
            {
              list_pe[id].push_back(petime(de, d1));
            }
            
            for(int i = 0; i < pe2; i++)
            {
              list_pe[-1*id].push_back(petime(de, d2));
            }            
          }
        }
      }
    }
}

void TestPoisson(double l, double nlayer, int nexp)
{  
  TH1I* h1 = new TH1I("h1","",100,0,100);
  TH1I* h2 = new TH1I("h2","",100,0,100);
  
  for(int i = 0; i < nexp; i++)
  {
    double s1 = 0;
    double s2 = r.Poisson(l*nlayer);
  
    for(int j = 0; j < nlayer; j++)
    {
      s1 += r.Poisson(l);
    }
    
    h1->Fill(s1);
    h2->Fill(s2);
  }
  
  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlack);
  
  TCanvas* c = new TCanvas();
  h1->Draw("E0");
  h2->Draw("E0same");
}

void Digitize(std::map<int, std::vector<double> >& list_pe, std::map<int, double>& adc, std::map<int, double>& tdc)
{
  /*
    -  ADC - Proportional to NPHE 
    -  TDC - Constant fraction - simulated 
             TPHE(1...NPHE) in increasing time order 
             IND_SEL= 0.15*NPHE            
             TDC_cell = TPHE(IND_SEL)
  */
  
  const double k = 1.0;
  
  for(std::map<int, std::vector<double> >::iterator it=list_pe.begin(); it != list_pe.end(); ++it)
  {
    adc[it->first] = k * it->second.size();
    std::sort(it->second.begin(), it->second.end());
    int index = 0.15 * it->second.size();
    tdc[it->first] = it->second[index];
  }
}

void analize(const char* finname, const char* foutname)
{
	  TChain* t = new TChain("EDepSimEvents","EDepSimEvents");
	  t->Add(finname);
	  TFile f(t->GetListOfFiles()->At(0)->GetTitle());
    TGeoManager* geo = (TGeoManager*) f.Get("EDepSimGeometry");
    
    TG4Event* ev = 0;
    t->SetBranchAddress("Event",&ev);
      
    std::map<int, std::vector<double> > list_pe;
    std::map<int, double> adc;
    std::map<int, double> tdc;
    
    TFile fout(foutname,"RECREATE");
    TTree tev("Events","Events");
    tev.Branch("Event","TG4Event",&ev);
    tev.Branch("cellPE","std::map<int, std::vector<double> >",&list_pe);
    tev.Branch("cellADC","std::map<int, double>",&adc);
    tev.Branch("cellTDC","std::map<int, double>",&tdc);
    
    const int nev = t->GetEntries();
    //const int nev = 100;
    
    std::cout << "Events: " << nev << " [";
    std::cout << std::setw(3) << int(0) << "%]" << std::flush;

    for(int i = 0; i < nev; i++)
    {
      t->GetEntry(i);
    
      std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;
      
      list_pe.clear();
      adc.clear();
      tdc.clear();
      
      SimulatePE(ev, geo, i, list_pe); 
      Digitize(list_pe,adc,tdc);
         
      tev.Fill();
    }
    std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
    std::cout << std::endl;
    
    fout.cd();
    tev.Write();
    geo->Write();
    fout.Close();
    
    delete t;
    f.Close();
}

void checkfast(const char* fname)
{
  TFile f(fname);
  TTree* t = (TTree*) f.Get("Events");
      
  std::map<int, std::vector<double> >* list_pe = new std::map<int, std::vector<double> >;
  std::map<int, double>* adc = new std::map<int, double>;
  std::map<int, double>* tdc = new std::map<int, double>;
  
  t->SetBranchAddress("cellPE",  &list_pe);
  t->SetBranchAddress("cellADC", &adc);
  t->SetBranchAddress("cellTDC", &tdc);
  
  TH1I* h_id  = new TH1I("h_id","Cell ID",50000,-25000,25000);
  TH1I* h_pe  = new TH1I("h_pe","pe",100,0,100);
  TH1I* h_adc = new TH1I("h_adc","adc",100,0,100);
  TH1D* h_tdc = new TH1D("h_tdc","tdc",200,0,100);
  TH1D* h_tdif = new TH1D("h_tdif","diff",500,-100,100);
  TH1D* h_tsum = new TH1D("h_tsum","sum",200,0,100);
  TH2D* h1 = new TH2D("h1","",100,0,100,100,0,100);
  
  h_id ->SetDirectory(0);
  h_pe ->SetDirectory(0);
  h_adc->SetDirectory(0);
  h_tdc->SetDirectory(0);
  
  h_tdif->SetDirectory(0);
  h_tsum->SetDirectory(0);
  h1->SetDirectory(0);
  
  //const int nev = t->GetEntries();
  const int nev = 1000;
    
  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;
  
  for(int i = 0; i < nev; i++)
  {
    t->GetEntry(i);
    
    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;
    
    for(std::map<int, std::vector<double> >::iterator it=list_pe->begin(); it != list_pe->end(); ++it)
    {
      h_id->Fill(it->first);
      h_pe->Fill(it->second.size());
    }
    
    for(std::map<int, double>::iterator it=adc->begin(); it != adc->end(); ++it)
    {
      h_adc->Fill(it->second);
      h1->Fill(it->second,(*tdc)[it->first]);
    }
    
    for(std::map<int, double>::iterator it=tdc->begin(); it != tdc->end(); ++it)
    {
      h_tdc->Fill(it->second);
      if(it->first > 0)
      {
        if(it->second > 0 && (*tdc)[-1*it->first] > 0)
        {
          h_tsum->Fill(it->second + (*tdc)[-1*it->first]);
          h_tdif->Fill(it->second - (*tdc)[-1*it->first]);
        }
      }
    }
  }
  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;
  
  TCanvas* c1 = new TCanvas();
  h_id->Draw();
  TCanvas* c2 = new TCanvas();
  h_pe->Draw();
  TCanvas* c3 = new TCanvas();
  h_adc->Draw();
  TCanvas* c4 = new TCanvas();
  h_tdc->Draw();
  TCanvas* c5 = new TCanvas();
  h_tsum->Draw();
  TCanvas* c6 = new TCanvas();
  h_tdif->Draw();
  TCanvas* c7 = new TCanvas();
  h1->Draw("colz");
}
