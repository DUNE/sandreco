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
#include <TStyle.h>

#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4Event.h"
#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4HitSegment.h"

#include <vector>
#include <map>
#include <iostream>

// Energy MeV
// Distance mm
// Time ns
 
TRandom3 r;

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
  TH1D* h_tdc = new TH1D("h_tdc","tdc",200,10,40);
  
  TH1D* h_tdif = new TH1D("h_tdif","diff",100,-20,20);
  TH1D* h_tsum = new TH1D("h_tsum","sum",100,30,60);
  
  TH2D* h1 = new TH2D("h1","",100,0,100,100,0,100);
  
  TH1D* hlaysum1 = new TH1D("hlaysum1","",100,0,20);
  TH1D* hlaysum2 = new TH1D("hlaysum2","",100,0,20);
  TH1D* hlaysum3 = new TH1D("hlaysum3","",100,0,20);
  TH1D* hlaysum4 = new TH1D("hlaysum4","",100,0,20);
  TH1D* hlaysum5 = new TH1D("hlaysum5","",100,0,20);
  
  TH1D* hlaysumADC1 = new TH1D("hlaysumADC1","",100,0,5000);
  TH1D* hlaysumADC2 = new TH1D("hlaysumADC2","",100,0,5000);
  TH1D* hlaysumADC3 = new TH1D("hlaysumADC3","",100,0,5000);
  TH1D* hlaysumADC4 = new TH1D("hlaysumADC4","",100,0,5000);
  TH1D* hlaysumADC5 = new TH1D("hlaysumADC5","",100,0,5000);
  
  TH1D* htres = new TH1D("htres","",100,4,14);
  
  TH1D* hADCtotSum = new TH1D("hADCtotSum","",100,1500,5000);
  TH1D* hADCtotSumContained = new TH1D("hADCtotSumContained","",100,1500,5000);
  
  h_id ->SetDirectory(0);
  h_pe ->SetDirectory(0);
  h_adc->SetDirectory(0);
  h_tdc->SetDirectory(0);
  
  h_tdif->SetDirectory(0);
  h_tsum->SetDirectory(0);
  h1->SetDirectory(0);
  
  hlaysum1->SetDirectory(0);
  hlaysum2->SetDirectory(0);
  hlaysum3->SetDirectory(0);
  hlaysum4->SetDirectory(0);
  hlaysum5->SetDirectory(0);
  
  hlaysumADC1->SetDirectory(0);
  hlaysumADC2->SetDirectory(0);
  hlaysumADC3->SetDirectory(0);
  hlaysumADC4->SetDirectory(0);
  hlaysumADC5->SetDirectory(0);
  
  htres->SetDirectory(0);
  
  hADCtotSum->SetDirectory(0);
  hADCtotSumContained->SetDirectory(0);
  
  const int nev = t->GetEntries();
    
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
    
    double ADClay1 = 0;
    double ADClay2 = 0;
    double ADClay3 = 0;
    double ADClay4 = 0;
    double ADClay5 = 0;
    
    for(std::map<int, double>::iterator it=adc->begin(); it != adc->end(); ++it)
    {
      h_adc->Fill(it->second);
      h1->Fill(it->second,(*tdc)[it->first]);
      
      int id_lay = TMath::Abs(int(it->first/100)%10);
          
      switch(id_lay)
      {
         case 0:
           ADClay1 += it->second;
         break;
         case 1:
           ADClay2 += it->second;
         break;
         case 2:
           ADClay3 += it->second;
         break;
         case 3:
           ADClay4 += it->second;
         break;
         case 4:
           ADClay5 += it->second;
         break;
         default:
           std::cout << "error: layer id: " << id_lay << std::endl;
         break; 
      }
    }
    
    hlaysumADC1->Fill(ADClay1);
    hlaysumADC2->Fill(ADClay2);
    hlaysumADC3->Fill(ADClay3);
    hlaysumADC4->Fill(ADClay4);
    hlaysumADC5->Fill(ADClay5);
    
    if(ADClay5 < 175)
      hADCtotSumContained->Fill(ADClay1+ADClay2+ADClay3+ADClay4+ADClay5);
    
    hADCtotSum->Fill(ADClay1+ADClay2+ADClay3+ADClay4+ADClay5);
    
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
      
      if(it->first > 0)
      {
        if(it->second > 0 && (*tdc)[-1*it->first] > 0)
        {
          int id_lay = int(it->first/100)%10;
          
          switch(id_lay)
          {
             case 0:
               hlaysum1->Fill(it->second + (*tdc)[-1*it->first]);
               htres->Fill(0.5*(it->second + (*tdc)[-1*it->first]) - 0.5*4.3*5.85);
             break;
             case 1:
               hlaysum2->Fill(it->second + (*tdc)[-1*it->first]);
             break;
             case 2:
               hlaysum3->Fill(it->second + (*tdc)[-1*it->first]);
             break;
             case 3:
               hlaysum4->Fill(it->second + (*tdc)[-1*it->first]);
             break;
             case 4:
               hlaysum5->Fill(it->second + (*tdc)[-1*it->first]);
             break;
             default:
               std::cout << "error: layer id: " << id_lay << std::endl;
             break; 
          }
        }
      }
    }
  }
  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;
  
  gStyle->SetOptFit(1);
  
  TCanvas* c1 = new TCanvas();
  h_id->Draw();
  TCanvas* c2 = new TCanvas();
  h_pe->Draw();
  TCanvas* c3 = new TCanvas();
  h_adc->Draw();
  TCanvas* c4 = new TCanvas();
  h_tdc->Draw();
  TCanvas* c5 = new TCanvas();
  h_tsum->Fit("gaus");
  h_tsum->Draw("E0");
  TCanvas* c6 = new TCanvas();
  h_tdif->Fit("gaus");
  h_tdif->Draw("E0");
  TCanvas* c7 = new TCanvas();
  h1->Draw("colz");
  TCanvas* c8 = new TCanvas();
  hlaysum1->Fit("gaus");
  hlaysum1->Draw("E0");
  TCanvas* c9 = new TCanvas();
  hlaysum2->Fit("gaus");
  hlaysum2->Draw("E0");
  TCanvas* c10 = new TCanvas();
  hlaysum3->Fit("gaus");
  hlaysum3->Draw("E0");
  TCanvas* c11 = new TCanvas();
  hlaysum4->Fit("gaus");
  hlaysum4->Draw("E0");
  TCanvas* c12 = new TCanvas();
  hlaysum5->Fit("gaus");
  hlaysum5->Draw("E0");
  
  TCanvas* c13 = new TCanvas();
  //hlaysumADC1->Fit("landau");
  hlaysumADC1->Draw("E0");
  TCanvas* c14 = new TCanvas();
  //hlaysumADC2->Fit("landau");
  hlaysumADC2->Draw("E0");
  TCanvas* c15 = new TCanvas();
  //hlaysumADC3->Fit("landau");
  hlaysumADC3->Draw("E0");
  TCanvas* c16 = new TCanvas();
  //hlaysumADC4->Fit("landau");
  hlaysumADC4->Draw("E0");
  TCanvas* c17 = new TCanvas();
  //hlaysumADC5->Fit("landau");
  hlaysumADC5->Draw("E0");
  TCanvas* c18 = new TCanvas();
  hADCtotSum->Fit("gaus");
  hADCtotSum->Draw("E0");
  TCanvas* c19 = new TCanvas();
  hADCtotSumContained->Fit("gaus");
  hADCtotSumContained->Draw("E0");
  
  TCanvas* c20 = new TCanvas();
  htres->Fit("gaus");
  htres->Draw("E0");
}

// time resolution with muons
void TimeRes()
{
  TFile f("calo_cal_muons_1k.01.root");
  TTree* t = (TTree*) f.Get("Events");
      
  std::map<int, std::vector<double> >* list_pe = new std::map<int, std::vector<double> >;
  std::map<int, double>* adc = new std::map<int, double>;
  std::map<int, double>* tdc = new std::map<int, double>;
  
  t->SetBranchAddress("cellPE",  &list_pe);
  t->SetBranchAddress("cellADC", &adc);
  t->SetBranchAddress("cellTDC", &tdc);
  
  TH1D* h_tdc = new TH1D("h_tdc","tdc",200,10,40);
  
  TH1D* h_tdif = new TH1D("h_tdif","diff",100,-10,10);
  TH1D* h_tsum = new TH1D("h_tsum","sum",100,30,60);
  
  TH1D* hlaysum1 = new TH1D("hlaysum1","",100,35,55);
  TH1D* hlaysum2 = new TH1D("hlaysum2","",100,35,55);
  TH1D* hlaysum3 = new TH1D("hlaysum3","",100,35,55);
  TH1D* hlaysum4 = new TH1D("hlaysum4","",100,35,55);
  TH1D* hlaysum5 = new TH1D("hlaysum5","",100,35,55);
  
  TH1D* htres = new TH1D("htres","",100,4,14);
  
  h_tdc->SetDirectory(0);
  h_tdif->SetDirectory(0);
  h_tsum->SetDirectory(0);
  
  hlaysum1->SetDirectory(0);
  hlaysum2->SetDirectory(0);
  hlaysum3->SetDirectory(0);
  hlaysum4->SetDirectory(0);
  hlaysum5->SetDirectory(0);
  
  htres->SetDirectory(0);
  
  const int nev = t->GetEntries();
    
  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;
  
  for(int i = 0; i < nev; i++)
  {
    t->GetEntry(i);
    
    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;
    
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
      
      if(it->first > 0)
      {
        if(it->second > 0 && (*tdc)[-1*it->first] > 0)
        {
          int id_lay = int(it->first/100)%10;
          
          switch(id_lay)
          {
             case 0:
               hlaysum1->Fill(it->second + (*tdc)[-1*it->first]);
               htres->Fill(0.5*(it->second + (*tdc)[-1*it->first]) - 0.5*4.3*5.85);
             break;
             case 1:
               hlaysum2->Fill(it->second + (*tdc)[-1*it->first]);
             break;
             case 2:
               hlaysum3->Fill(it->second + (*tdc)[-1*it->first]);
             break;
             case 3:
               hlaysum4->Fill(it->second + (*tdc)[-1*it->first]);
             break;
             case 4:
               hlaysum5->Fill(it->second + (*tdc)[-1*it->first]);
             break;
             default:
               std::cout << "error: layer id: " << id_lay << std::endl;
             break; 
          }
        }
      }
    }
  }
  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;
  
  gStyle->SetOptFit(1);
  
  TCanvas* c4 = new TCanvas();
  h_tdc->Draw();
  TCanvas* c5 = new TCanvas();
  h_tsum->Fit("gaus");
  h_tsum->Draw("E0");
  TCanvas* c6 = new TCanvas();
  h_tdif->Fit("gaus");
  h_tdif->Draw("E0");
  TCanvas* c8 = new TCanvas();
  hlaysum1->Fit("gaus");
  hlaysum1->Draw("E0");
  TCanvas* c9 = new TCanvas();
  hlaysum2->Fit("gaus");
  hlaysum2->Draw("E0");
  TCanvas* c10 = new TCanvas();
  hlaysum3->Fit("gaus");
  hlaysum3->Draw("E0");
  TCanvas* c11 = new TCanvas();
  hlaysum4->Fit("gaus");
  hlaysum4->Draw("E0");
  TCanvas* c12 = new TCanvas();
  hlaysum5->Fit("gaus");
  hlaysum5->Draw("E0");
  
  TCanvas* c20 = new TCanvas();
  htres->Fit("gaus");
  htres->Draw("E0");
}

void EnergyRes()
{
  TFile f("photons_1GeV_10k.root");
  TTree* t = (TTree*) f.Get("Events");
      
  std::map<int, std::vector<double> >* list_pe = new std::map<int, std::vector<double> >;
  std::map<int, double>* adc = new std::map<int, double>;
  std::map<int, double>* tdc = new std::map<int, double>;
  
  t->SetBranchAddress("cellPE",  &list_pe);
  t->SetBranchAddress("cellADC", &adc);
  t->SetBranchAddress("cellTDC", &tdc);
  
  TH1I* h_pe  = new TH1I("h_pe","pe",100,0,100);
  TH1I* h_adc = new TH1I("h_adc","adc",100,0,100);
  
  TH1D* hlaysumADC1 = new TH1D("hlaysumADC1","",100,0,5000);
  TH1D* hlaysumADC2 = new TH1D("hlaysumADC2","",100,0,5000);
  TH1D* hlaysumADC3 = new TH1D("hlaysumADC3","",100,0,5000);
  TH1D* hlaysumADC4 = new TH1D("hlaysumADC4","",100,0,5000);
  TH1D* hlaysumADC5 = new TH1D("hlaysumADC5","",100,0,5000);
  
  TH1D* hADCtotSum = new TH1D("hADCtotSum","",100,2000,5000);
  TH1D* hADCtotSumContained = new TH1D("hADCtotSumContained","",100,2000,5000);
  
  h_pe ->SetDirectory(0);
  h_adc->SetDirectory(0);
  
  hlaysumADC1->SetDirectory(0);
  hlaysumADC2->SetDirectory(0);
  hlaysumADC3->SetDirectory(0);
  hlaysumADC4->SetDirectory(0);
  hlaysumADC5->SetDirectory(0);
  
  hADCtotSum->SetDirectory(0);
  hADCtotSumContained->SetDirectory(0);
  
  const int nev = t->GetEntries();
    
  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;
  
  for(int i = 0; i < nev; i++)
  {
    t->GetEntry(i);
    
    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;
    
    for(std::map<int, std::vector<double> >::iterator it=list_pe->begin(); it != list_pe->end(); ++it)
    {
      h_pe->Fill(it->second.size());
    }
    
    double ADClay1 = 0;
    double ADClay2 = 0;
    double ADClay3 = 0;
    double ADClay4 = 0;
    double ADClay5 = 0;
    
    for(std::map<int, double>::iterator it=adc->begin(); it != adc->end(); ++it)
    {
      h_adc->Fill(it->second);
      
      int id_lay = TMath::Abs(int(it->first/100)%10);
          
      switch(id_lay)
      {
         case 0:
           ADClay1 += it->second;
         break;
         case 1:
           ADClay2 += it->second;
         break;
         case 2:
           ADClay3 += it->second;
         break;
         case 3:
           ADClay4 += it->second;
         break;
         case 4:
           ADClay5 += it->second;
         break;
         default:
           std::cout << "error: layer id: " << id_lay << std::endl;
         break; 
      }
    }
    
    hlaysumADC1->Fill(ADClay1);
    hlaysumADC2->Fill(ADClay2);
    hlaysumADC3->Fill(ADClay3);
    hlaysumADC4->Fill(ADClay4);
    hlaysumADC5->Fill(ADClay5);
    
    if(ADClay5 < 175)
      hADCtotSumContained->Fill(ADClay1+ADClay2+ADClay3+ADClay4+ADClay5);
    
    hADCtotSum->Fill(ADClay1+ADClay2+ADClay3+ADClay4+ADClay5);
    
  }
  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;
  
  gStyle->SetOptFit(1);
  
  TCanvas* c2 = new TCanvas();
  h_pe->Draw();
  TCanvas* c3 = new TCanvas();
  h_adc->Draw();
  
  TCanvas* c13 = new TCanvas();
  hlaysumADC1->Draw("E0");
  TCanvas* c14 = new TCanvas();
  hlaysumADC2->Draw("E0");
  TCanvas* c15 = new TCanvas();
  hlaysumADC3->Draw("E0");
  TCanvas* c16 = new TCanvas();
  hlaysumADC4->Draw("E0");
  TCanvas* c17 = new TCanvas();
  hlaysumADC5->Draw("E0");
  TCanvas* c18 = new TCanvas();
  hADCtotSum->Fit("gaus");
  hADCtotSum->Draw("E0");
  TCanvas* c19 = new TCanvas();
  hADCtotSumContained->Fit("gaus");
  hADCtotSumContained->Draw("E0");
}
