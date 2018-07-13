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
#include <algorithm>

// Energy MeV
// Distance mm
// Time ns
 
void partCaloDE(const char* fname)
{
  TFile f(fname);
  TTree* t = (TTree*) f.Get("Events");
      
  std::map<int, std::vector<double> >* list_pe = new std::map<int, std::vector<double> >;
  std::map<int, std::vector<int> >* list_hit = new std::map<int, std::vector<int> >;
  std::map<int, double>* adc = new std::map<int, double>;
  std::map<int, double>* tdc = new std::map<int, double>;
  
  t->SetBranchAddress("cellPE",  &list_pe);
  t->SetBranchAddress("hitIndex",  &list_hit);
  t->SetBranchAddress("cellADC", &adc);
  t->SetBranchAddress("cellTDC", &tdc);
  
  const int nev = t->GetEntries();

  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;
 
  TH1D* hNpartCell = new TH1D("hNpartCell","",1000,0,1000);
  
  hNpartCell->SetDirectory(0);
  
  std::vector<int> idhits;
  
  for(int i = 0; i < nev; i++)
  {
    t->GetEntry(i);
    
    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;
    
	/*
    for(std::map<int, double>::iterator it=tdc->begin(); it != tdc->end(); ++it)
    {
    }
	*/

  	for(std::map<int, double>::iterator it1=adc->begin(); it1 != adc->end(); ++it1)
  	{
       int cellID = it1->first;
       
       if(cellID < 0)
         continue;
       
       std::map<int, double>::iterator it2 = adc->find(-cellID);
       
       if(it2 == adc->end())
         continue;
       
       int layer = (cellID - (cellID % 100)) % 1000;
       
       if(layer != 0)
         continue;
       
  		 if(it1->second > 0 && it2->second > 0)
       {
         double sig = it1->second > 0 && it2->second > 0;
         
         idhits = list_hit->at(cellID);
         idhits.insert(idhits.end(), list_hit->at(-cellID).begin(), list_hit->at(-cellID).end());
         
         std::sort(idhits.begin(),idhits.end());
         std::unique(idhits.begin(),idhits.end());
         
         hNpartCell->Fill(idhits.size());
       }
  	}
  }

  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;
  
  gStyle->SetOptFit(1);
  
  TCanvas* c1 = new TCanvas();
  hNpartCell->Draw();
}

