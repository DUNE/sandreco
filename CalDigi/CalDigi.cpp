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

#include "/mnt/nas01/users/mtenti/wd/analysis/KLOEcal/loader/loader.C"

// Energy MeV
// Distance mm
// Time ns

namespace ns_CalDigi {
  const bool debug = false;
  
  static const int nMod = 24;
  static const int nLay = 5;
  static const int nCel = 12;
 
  TRandom3 r;
}

double Attenuation(double d, int planeID)
{
/*
     dE/dx attenuation - Ea=p1*exp(-d/atl1)+(1.-p1)*exp(-d/atl2)
       d    distance from photocatode - 2 cells/cell; d1 and d2 
      atl1  50. cm
      atl2  430 cm planes 1-2    innermost plane is 1 
            380 cm plane 3
            330 cm planes 4-5
       p1   0.35
*/
  const double p1 = 0.35;
  const double alt1 = 500.;
  double alt2 = 0.0;
  
  switch (planeID)
  {
    case 0: 
    case 1: 
      alt2 = 4300.0;
    break;
     
    case 2: 
      alt2 = 3800.0;
    break;
     
    case 3: 
    case 4: 
      alt2 = 3300.0;
    break;
     
    default: 
      //std::cout << "planeID out if range" << std::endl;
      alt2 = -999.0;
    break;
  }
  
  if(ns_CalDigi::debug)
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
  
  if(ns_CalDigi::debug)
    std::cout << "E = " << E << " -> p.e. = " << e2p2*E << std::endl;
  
  return e2p2*E;
}

double petime(double t0, double d)
{
  /*
     - For each photoelectron: Time for TDC simulation obtained from 

C  PHOTOELECTRON TIME :  Particle TIME in the cell  
C                      + SCINTILLATION DECAY TIME + 
C                      + signal propagation to the cell
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
  
  double tdec = tscin * TMath::Power(1./ns_CalDigi::r.Uniform()-1.,tscex);
  
  double time = t0 + tdec + vlfb * d * mm_to_m + ns_CalDigi::r.Gaus();
  
  if(ns_CalDigi::debug)
  {
    std::cout << "time : " << time << std::endl;
    std::cout << "t0   : " << t0 << std::endl;
    std::cout << "scint: " << tdec << std::endl;
    std::cout << "prop : " << vlfb * d * mm_to_m << std::endl;
  }
    
  
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
  
  if(ns_CalDigi::debug)
  {
    std::cout << "hit" << std::endl;
    std::cout << "\t[x,y,z]                " << x << " " << y << " " << z << std::endl;
    std::cout << "\t[modID,planeID,cellID] " << modID << " " << planeID << " " << cellID << std::endl;
    std::cout << "\t[d1,d2,t,de]           " << d1 << " " << d2 << " " << t << " " << de << std::endl;
  }
  
  return true;
}

void SimulatePE(TG4Event* ev, TGeoManager* g, int index, std::map<int, std::vector<double> >& time_pe, std::map<int, std::vector<int> >& id_hit)
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
            
            int pe1 = ns_CalDigi::r.Poisson(ave_pe1);
            int pe2 = ns_CalDigi::r.Poisson(ave_pe2);
            
            id = cellID + 100 * planeID + 1000 * modID;
            
            if(ns_CalDigi::debug)
            {
              std::cout << "cell ID: " << id << std::endl;
              std::cout << "\t" << de << " " << en1 << " " << en2 << std::endl;
              std::cout << "\t" << ave_pe1 << " " << ave_pe2 << std::endl;
              std::cout << "\t" << pe1 << " " << pe2 << std::endl;
            }
            
            //cell 1 -> x < 0 -> ID > 0
            //cell 2 -> x > 0 -> ID < 0
            
            for(int i = 0; i < pe1; i++)
            {
              time_pe[id].push_back(petime(t0, d1));
              id_hit[id].push_back(j);
            }
            
            for(int i = 0; i < pe2; i++)
            {
              time_pe[-1*id].push_back(petime(t0, d2));
              id_hit[-1*id].push_back(j);
            }            
          }
        }
      }
    }
}

void Digitize(std::map<int, std::vector<double> >& time_pe, std::map<int, double>& adc, std::map<int, double>& tdc)
{
  /*
    -  ADC - Proportional to NPHE 
    -  TDC - Constant fraction - simulated 
             TPHE(1...NPHE) in increasing time order 
             IND_SEL= 0.15*NPHE            
             TDC_cell = TPHE(IND_SEL)
  */
  
  const double pe2ADC = 1.0;
  
  for(std::map<int, std::vector<double> >::iterator it=time_pe.begin(); it != time_pe.end(); ++it)
  {
    adc[it->first] = pe2ADC * it->second.size();
    std::sort(it->second.begin(), it->second.end());
    int index = 0.15 * it->second.size();
    tdc[it->first] = it->second[index];
  }
}

void DigitizeCal(const char* finname, const char* foutname)
{
	  TChain* t = new TChain("EDepSimEvents","EDepSimEvents");
	  t->Add(finname);
	  TFile f(t->GetListOfFiles()->At(0)->GetTitle());
    TGeoManager* geo = (TGeoManager*) f.Get("EDepSimGeometry");
  
    const char* path_template = "volWorld_PV/volDetEnclosure_PV_0/volKLOEFULLECALSENSITIVE_EXTTRK_NEWGAP_PV_0/KLOEBarrelECAL_%d_volume_PV_0";
  
    TGeoTrd2* mod = (TGeoTrd2*) geo->FindVolumeFast("KLOEBarrelECAL_0_volume_PV")->GetShape();
    
    double xmax = mod->GetDx1();
    double xmin = mod->GetDx2();
    double dz = mod->GetDz();
    
    double dzlay[ns_CalDigi::nLay+1] = {115, 115-22, 115-22-22, 115-22-22-22, 115-22-22-22-22, 115-22-22-22-22-27};
    double czlay[ns_CalDigi::nLay];
    double cxlay[ns_CalDigi::nLay][ns_CalDigi::nCel];
    
    for(int i = 0; i < ns_CalDigi::nLay; i++)
    {
      czlay[i] = 0.5 * (dzlay[i] + dzlay[i+1]);
      
      double dx = xmax - (xmax - xmin)/dz * czlay[i];
      
      for(int j = 0; j < ns_CalDigi::nCel; j++)
      {
        cxlay[i][j] = 2*dx/ns_CalDigi::nCel * (j - ns_CalDigi::nCel/2. + 0.5);
      }
    }    
    
    TG4Event* ev = new TG4Event;
    t->SetBranchAddress("Event",&ev);
      
    std::map<int, std::vector<double> > time_pe;
    std::map<int, std::vector<int> > id_hit;
    std::map<int, double> adc;
    std::map<int, double> tdc;
    
    std::vector<cell> vec_cell;
    
    TFile fout(foutname,"RECREATE");
    TTree tcal("tCalDigi","KLOE Calorimeter Digitization");
    tcal.Branch("cell","std::vector<cell>",&vec_cell);
    
    const int nev = t->GetEntries();
    
    std::cout << "Events: " << nev << " [";
    std::cout << std::setw(3) << int(0) << "%]" << std::flush;

    for(int i = 0; i < nev; i++)
    {
      t->GetEntry(i);
    
      std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;
      
      time_pe.clear();
      id_hit.clear();
      adc.clear();
      tdc.clear();
      vec_cell.clear();
      
      SimulatePE(ev, geo, i, time_pe, id_hit); 
      Digitize(time_pe,adc,tdc);
      
      for(std::map<int, std::vector<double> >::iterator it = time_pe.begin(); it != time_pe.end(); ++it)
      {
        if(it->first < 0)
          continue;
        
        cell c;
        c.id = it->first;
        c.adc1 = adc[it->first];
        c.tdc1 = tdc[it->first];
        c.adc2 = adc[-1*it->first];
        c.tdc2 = tdc[-1*it->first];
        c.pe_time1 = time_pe[it->first];
        c.pe_time2 = time_pe[-1*it->first];
        c.hindex1 = id_hit[it->first];
        c.hindex2 = id_hit[-1*it->first];
            
        c.mod = c.id / 1000;
        c.lay = (c.id - c.mod * 1000) / 100;
        c.cel = c.id - c.mod * 1000 - c.lay *100;
                
        double dummyLoc[3];
        double dummyMas[3];
        
        dummyLoc[0] = cxlay[c.lay][c.cel];
        dummyLoc[1] = 0.;
        dummyLoc[2] = czlay[c.lay];
        
        geo->cd(TString::Format(path_template,c.mod).Data());
        
        geo->LocalToMaster(dummyLoc, dummyMas);
        
        c.y = dummyMas[1];
        c.z = dummyMas[2];
        
        vec_cell.push_back(c);
      }
         
      tcal.Fill();
    }
    std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
    std::cout << std::endl;
    
    fout.cd();
    tcal.Write();
    fout.Close();
    
    delete t;
    f.Close();
}

void CalDigi()
{
  std::cout << "DigitizeCal(const char* finname, const char* foutname)" << std::endl;
  std::cout << "finname could contain wild card" << std::endl;
}