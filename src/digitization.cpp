#include <TGeoManager.h>
#include <TString.h>
#include <TGeoNode.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TGeoTrd2.h>
#include <TGeoTube.h>
#include <TChain.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TSystem.h>

#include "TG4Event.h"
#include "TG4HitSegment.h"

#include <vector>
#include <map>
#include <iostream>

#include "struct.h"
#include "utils.h"

// Energy MeV
// Distance mm
// Time ns

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
  
  if(ns_Digit::debug)
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
  
  if(ns_Digit::debug)
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

  TRandom3 r(0);
  
  double mm_to_m = 1E-3;
  
  double tdec = ns_Digit::tscin * TMath::Power(1./r.Uniform()-1.,ns_Digit::tscex);
  
  double time = t0 + tdec + ns_Digit::vlfb * d * mm_to_m + r.Gaus();
  
  if(ns_Digit::debug)
  {
    std::cout << "time : " << time << std::endl;
    std::cout << "t0   : " << t0 << std::endl;
    std::cout << "scint: " << tdec << std::endl;
    std::cout << "prop : " << ns_Digit::vlfb * d * mm_to_m << std::endl;
  }
  
  return time;
}

bool ProcessHit(TGeoManager* g, const TG4HitSegment& hit, int& modID, int& planeID, int& cellID, double& d1, double& d2, double& t, double& de)
{
  if(ns_Digit::debug)
  {
    std::cout << "ProcessHit" << std::endl;
  }

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
    
  if(ns_Digit::debug)
  {
    std::cout << "node name: " << str.Data() << std::endl;
  }
  if(str.Contains("KLOEEndcapECALL_volume_PV_0") == true || str.Contains("KLOEEndcapECALR_volume_PV_0") == true)
  {
    return false;
  }
  else if(str.Contains("KLOEBarrelECAL") == true && str.Contains("lead_slab") == false)
  {
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
    
    if(ns_Digit::debug)
    {
      std::cout << "pointer: " << trd << std::endl;
    }
    
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
    
    if(ns_Digit::debug)
    {
      std::cout << "hit: " << str.Data() << std::endl;
      std::cout << "\t[x,y,z]                " << x << " " << y << " " << z << std::endl;
      std::cout << "\t[modID,planeID,cellID] " << modID << " " << planeID << " " << cellID << std::endl;
      std::cout << "\t[d1,d2,t,de]           " << d1 << " " << d2 << " " << t << " " << de << std::endl;
    }
    
    return true;
  }
  else if(str.Contains("KLOEEndcapECALR") == true && str.Contains("lead_slab") == false)
  {
    TObjArray* obja = str.Tokenize("_");
    
    int slabID;
    modID  = 30;
    slabID = ((TObjString*) obja->At(4))->GetString().Atoi();
    
    delete obja;
    
    // planeID==0 -> smallest slab
    // planeID==208 -> biggest slab
    planeID = (208 - slabID)/40;
    
    if (planeID > 4) planeID = 4;
    
    double Pmaster[3];
    double Plocal[3];
    Pmaster[0] = x;
    Pmaster[1] = y;
    Pmaster[2] = z;
    
    g->GetCurrentNavigator()->MasterToLocal(Pmaster,Plocal);
    
    TGeoTube* tub = (TGeoTube*) node->GetVolume()->GetShape();
    
    if(ns_Digit::debug)
    {
      std::cout << "pointer: " << tub << std::endl;
    }
    
    double rmin = tub->GetRmin();
    double rmax = tub->GetRmax();
    double dz  = tub->GetDz();
    
    d1 = rmax * TMath::Sin(TMath::ACos(Plocal[0]/rmax)) - Plocal[1];
    d2 = rmax * TMath::Sin(TMath::ACos(Plocal[0]/rmax)) + Plocal[1];
    
    cellID = int((Plocal[0]/rmax + 1.) * 45);
    
    if(ns_Digit::debug)
    {
      std::cout << "hit: " << str.Data() << std::endl;
      std::cout << "\t[x,y,z]                " << x << " " << y << " " << z << std::endl;
      std::cout << "\t[modID,planeID,cellID] " << modID << " " << planeID << " " << cellID << std::endl;
    }
    
    return true;
  }
  else if(str.Contains("KLOEEndcapECALL") == true && str.Contains("lead_slab") == false)
  {
    TObjArray* obja = str.Tokenize("_");
    
    int slabID;
    modID  = 40;
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
    
    TGeoTube* tub = (TGeoTube*) node->GetVolume()->GetShape();
    
    if(ns_Digit::debug)
    {
      std::cout << "pointer: " << tub << std::endl;
    }
    
    double rmin = tub->GetRmin();
    double rmax = tub->GetRmax();
    double dz  = tub->GetDz();
    
    d1 = rmax * TMath::Sin(TMath::ACos(Plocal[0]/rmax)) - Plocal[1];
    d2 = rmax * TMath::Sin(TMath::ACos(Plocal[0]/rmax)) + Plocal[1];
    
    cellID = int((Plocal[0]/rmax + 1.) * 45);
    
    if(ns_Digit::debug)
    {
      std::cout << "hit: " << str.Data() << std::endl;
      std::cout << "\t[x,y,z]                " << x << " " << y << " " << z << std::endl;
      std::cout << "\t[modID,planeID,cellID] " << modID << " " << planeID << " " << cellID << std::endl;
    }
    
    return true;
  }
  else
  {
    return false;
  }
}

void SimulatePE(TG4Event* ev, TGeoManager* g, std::map<int, std::vector<double> >& time_pe, std::map<int, std::vector<int> >& id_hit, std::map<int, double>& L)
{    
    int modID, planeID, cellID, id;
    double d1, d2, t0, de;

    TRandom3 r(0);

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
            
            if(ns_Digit::debug)
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
              L[id] = d1 + d2;
            }
            
            for(int i = 0; i < pe2; i++)
            {
              time_pe[-1*id].push_back(petime(t0, d2));
              id_hit[-1*id].push_back(j);
              L[-1*id] = d1 + d2;
            }            
          }
        }
      }
    }
}

void TimeAndSignal(std::map<int, std::vector<double> >& time_pe, std::map<int, double>& adc, std::map<int, double>& tdc)
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

void CollectSignal(TGeoManager* geo, 
    std::map<int, std::vector<double> >& time_pe, 
    std::map<int, double>& adc, 
    std::map<int, double>& tdc, 
    std::map<int, double>& L, 
    std::map<int, std::vector<int> >& id_hit,
    std::vector<cell>& vec_cell)
{
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
      c.l = L[it->first];
      c.x = 0;
      c.y = 0;
      c.z = 0;
          
      c.mod = c.id / 1000;
      c.lay = (c.id - c.mod * 1000) / 100;
      c.cel = c.id - c.mod * 1000 - c.lay *100;
      
      double dummyLoc[3];
      double dummyMas[3];
      
      if(c.mod < 24)
      {        
        dummyLoc[0] = ns_Digit::cxlay[c.lay][c.cel];
        dummyLoc[1] = 0.;
        dummyLoc[2] = ns_Digit::czlay[c.lay];
        
        geo->cd(TString::Format(ns_Digit::path_barrel_template,c.mod).Data());
        
        geo->LocalToMaster(dummyLoc, dummyMas);
        
        c.x = dummyMas[0];
        c.y = dummyMas[1];
        c.z = dummyMas[2];
      }
      else if(c.mod == 30) // right x < 0
      {        
        dummyLoc[0] = ns_Digit::ec_r/45. * (0.5 + c.cel) - ns_Digit::ec_r;
        dummyLoc[1] = 0.;
        dummyLoc[2] = -ns_Digit::czlay[c.lay];
        
        geo->cd(ns_Digit::path_endcapR_template);
        
        geo->LocalToMaster(dummyLoc, dummyMas);
        
        c.x = dummyMas[0];
        c.y = dummyMas[1];
        c.z = dummyMas[2];
      }
      else if(c.mod == 40) // left x > 0
      {        
        dummyLoc[0] = ns_Digit::ec_r/45. * (0.5 + c.cel) - ns_Digit::ec_r;
        dummyLoc[1] = 0.;
        dummyLoc[2] = ns_Digit::czlay[c.lay];
        
        geo->cd(ns_Digit::path_endcapL_template);
        
        geo->LocalToMaster(dummyLoc, dummyMas);
        
        c.x = dummyMas[0];
        c.y = dummyMas[1];
        c.z = dummyMas[2];
      }
      
      vec_cell.push_back(c);
    }
}

void init(TGeoManager* geo)
{
    TGeoTrd2* mod = (TGeoTrd2*) geo->FindVolumeFast("KLOEBarrelECAL_0_volume_PV")->GetShape();
    
    double xmax = mod->GetDx1();
    double xmin = mod->GetDx2();
    double dz = mod->GetDz();
    
    for(int i = 0; i < ns_Digit::nLay; i++)
    {
      ns_Digit::czlay[i] = (ns_Digit::dzlay[i] + ns_Digit::dzlay[i+1]) - ns_Digit::dzlay[0];
      
      double dx = xmax - (xmax - xmin)/dz * (0.5 * (ns_Digit::dzlay[i] + ns_Digit::dzlay[i+1]));
      
      for(int j = 0; j < ns_Digit::nCel; j++)
      {
        ns_Digit::cxlay[i][j] = 2*dx/ns_Digit::nCel * (j - ns_Digit::nCel/2. + 0.5);
      }
    } 
      
    TGeoTube* ec = (TGeoTube*) geo->FindVolumeFast("KLOEEndcapECALL_volume_PV")->GetShape();
    
    ns_Digit::ec_r = ec->GetRmax();
    ns_Digit::ec_dz = ec->GetDz();
}

void DigitizeCal(TG4Event* ev, TGeoManager* geo, std::vector<cell>& vec_cell)
{     
    std::map<int, std::vector<double> > time_pe;
    std::map<int, std::vector<int> > id_hit;
    std::map<int, double> adc;
    std::map<int, double> tdc;
    std::map<int, double> L;
    
    vec_cell.clear();
    
    
    if(ns_Digit::debug)
    {
      std::cout << "SimulatePE" << std::endl;
    }    
    SimulatePE(ev, geo, time_pe, id_hit, L);  
    if(ns_Digit::debug)
    {
      std::cout << "TimeAndSignal" << std::endl;
    }
    TimeAndSignal(time_pe,adc,tdc); 
    if(ns_Digit::debug)
    {
      std::cout << "CollectSignal" << std::endl;
    }
    CollectSignal(geo, time_pe, adc, tdc, L, id_hit, vec_cell);
}

void Cluster(TG4Event* ev, TGeoManager* geo, std::map<std::string,std::vector<hit> >& cluster_map)
{
  cluster_map.clear();
    
  for(unsigned int j = 0; j < ev->SegmentDetectors["StrawTracker"].size(); j++)
  {
    const TG4HitSegment& hseg = ev->SegmentDetectors["StrawTracker"].at(j);
    
    double x = 0.5 * (hseg.Start.X() + hseg.Stop.X());
    double y = 0.5 * (hseg.Start.Y() + hseg.Stop.Y());
    double z = 0.5 * (hseg.Start.Z() + hseg.Stop.Z());
    
    std::string sttname = geo->FindNode(x,y,z)->GetName();
    
    hit h;
    h.det = sttname;
    h.x1 = hseg.Start.X();
    h.y1 = hseg.Start.Y();
    h.z1 = hseg.Start.Z();
    h.t1 = hseg.Start.T();
    h.x2 = hseg.Stop.X();
    h.y2 = hseg.Stop.Y();
    h.z2 = hseg.Stop.Z();
    h.t2 = hseg.Stop.T();
    h.de = hseg.EnergyDeposit;
    h.pid = hseg.PrimaryId;
    h.index = j;
    
    std::string cluster_name(sttname);
    cluster_name += "_" + std::to_string(hseg.PrimaryId);
    
    cluster_map[cluster_name].push_back(h);
  }
}

void Cluster2Digit(std::map<std::string,std::vector<hit> >& cluster_map, std::vector<digit>& digit_vec)
{    
  for(std::map<std::string,std::vector<hit> >::iterator it = cluster_map.begin(); it != cluster_map.end(); ++it)
  {      
    digit d;
    d.de = 0;    
    d.det = it->second[0].det;  
  
    for(unsigned int k = 0; k < it->second.size(); k++)
    {
      d.hindex.push_back(it->second[k].index);
      d.de += it->second[k].de;
    }
    
    d.hor = (d.det.find("STTPlane1FULL") != std::string::npos) ? false : true;
    
    std::sort(it->second.begin(), it->second.end(), isHitBefore);
    
    d.t = it->second.at(0).t1;
    d.x = 0.5 * (it->second.front().x1 + it->second.back().x2);
    d.y = 0.5 * (it->second.front().y1 + it->second.back().y2);
    d.z = 0.5 * (it->second.front().z1 + it->second.back().z2);
    
    digit_vec.push_back(d);
  }
}

void DigitizeStt(TG4Event* ev, TGeoManager* geo, std::vector<digit>& digit_vec)
{
  std::map<std::string,std::vector<hit> > cluster_map;
  digit_vec.clear();
  
  Cluster(ev, geo, cluster_map);
  Cluster2Digit(cluster_map, digit_vec);
}

void Digitize(const char* finname, const char* foutname)
{
	  //TChain* t = new TChain("EDepSimEvents","EDepSimEvents");
	  //t->Add(finname);
	  //TFile f(t->GetListOfFiles()->At(0)->GetTitle());
    TFile f(finname,"READ");
    TTree* t = (TTree*) f.Get("EDepSimEvents");
    TGeoManager* geo = (TGeoManager*) f.Get("EDepSimGeometry"); 
    TTree* gRooTracker = (TTree*) f.Get("DetSimPassThru/gRooTracker");
    TTree* InputKinem = (TTree*) f.Get("DetSimPassThru/InputKinem");
    TTree* InputFiles = (TTree*) f.Get("DetSimPassThru/InputFiles");
    
    init(geo);
    
    TG4Event* ev = new TG4Event;
    t->SetBranchAddress("Event",&ev);
  
    std::vector<digit> digit_vec;    
    std::vector<cell> vec_cell;
    
    TFile fout(foutname,"RECREATE");
    TTree tout("tDigit","Digitization");
    tout.Branch("cell","std::vector<cell>",&vec_cell);
    tout.Branch("Stt","std::vector<digit>",&digit_vec);
    
    //const int nev = t->GetEntries();
    const int nev = 5;
    
    std::cout << "Events: " << nev << " [";
    std::cout << std::setw(3) << int(0) << "%]" << std::flush;

    for(int i = 0; i < nev; i++)
    {
      t->GetEntry(i);
    
      std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;
      
      DigitizeCal(ev, geo, vec_cell);
      DigitizeStt(ev, geo, digit_vec);
         
      tout.Fill();
    }
    std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
    std::cout << std::endl;
    
    fout.cd();
    tout.Write();
    geo->Write();
    t->CloneTree()->Write();
    gRooTracker->CloneTree()->Write();
    InputKinem->CloneTree()->Write();;
    InputFiles->CloneTree()->Write();;
    fout.Close();
    
    f.Close();
}

void help_digit()
{
  std::cout << "Digitize <input file> <output file>" << std::endl;
  std::cout << "input file name could contain wild card" << std::endl;
} 

int main(int argc, char* argv[])
{
  if(argc != 3)
    help_digit();
  else
    Digitize(argv[1], argv[2]);
}
