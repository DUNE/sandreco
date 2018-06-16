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

#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4Event.h"
#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4HitSegment.h"

#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
#include <string>
#include <limits>

#ifdef __MAKECINT__ 
#pragma link C++ class map<int,vector<double> >+; 
#pragma link C++ class map<int,vector<int> >+;
#pragma link C++ class map<int,double>+; 
#endif

// Energy MeV
// Distance mm
// Time ns

const bool debug = false;
const int neutPDG = 2112;
const std::string hitCalName("EMCalSci");
const std::string hitSttName("StrawTracker");
const int nMaxPrim = 200;
const double rgap = 150.;
const double vel = 5.85 /*ns/m*/;
const double tconst = 4.3 /*m*/ * vel * 0.5;
const double c = 299.792458; //mm/ns

const char* path_template = "volWorld_PV/volDetEnclosure_PV_0/volKLOEFULLECALSENSITIVE_EXTTRK_NEWGAP_PV_0/KLOEBarrelECAL_%d_volume_PV_0";
 
void primNeut(const char* finname,const char* foutname)
{
  TFile f(finname);
  TTree* t = (TTree*) f.Get("Events");
  TTree* InputKinem = (TTree*) f.Get("InputKinem");
  TTree* InputFiles = (TTree*) f.Get("InputFiles");
  TTree* gRooTracker = (TTree*) f.Get("gRooTracker");
  
  TGeoManager* geo = (TGeoManager*) f.Get("EDepSimGeometry");
      
  std::map<int, std::vector<double> >* list_pe = new std::map<int, std::vector<double> >;
  std::map<int, std::vector<int> >* list_hit = new std::map<int, std::vector<int> >;
  std::map<int, double>* adc = new std::map<int, double>;
  std::map<int, double>* tdc = new std::map<int, double>;
  
  TG4Event* ev = new TG4Event();
  
  t->SetBranchAddress("cellPE",  &list_pe);
  t->SetBranchAddress("hitIndex",  &list_hit);
  t->SetBranchAddress("cellADC", &adc);
  t->SetBranchAddress("cellTDC", &tdc);
  t->SetBranchAddress("Event", &ev);
  
  int inputFileNum, inputEntryNum;
  
  InputKinem->SetBranchAddress("inputFileNum",&inputFileNum);
  InputKinem->SetBranchAddress("inputEntryNum",&inputEntryNum);
  
  int fileEntries;
  
  InputFiles->SetBranchAddress("fileEntries",&fileEntries);
  
  double StdHepP4[nMaxPrim][4];
  double EvtVtx[4];
  
  gRooTracker->SetBranchAddress("StdHepP4",StdHepP4);
  gRooTracker->SetBranchAddress("EvtVtx",EvtVtx);
  
  const int nev = t->GetEntries();
  //const int nev = 1000;

  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;
  
  TFile fout(foutname,"RECREATE");
  TTree tPrimaryNeut("tPrimaryNeut","Primary Neutrons");
  
  int eventID, runID;
  double pnu[4];
  double vtx[4];
  
  double p[4];
  double ptot;
  double mass;
  double ek;
  double ekCalo;
  
  int isSTTScat;
  int isScat;
  std::vector<double>* thSTTScat = new std::vector<double>();
  std::vector<double>* deSTTScat = new std::vector<double>();
  std::vector<double>* mtSTTScat = new std::vector<double>();
  std::vector<int>* sgSTTScat = new std::vector<int>();
  
  std::vector<double>* thScat = new std::vector<double>();
  std::vector<double>* deScat = new std::vector<double>();
  std::vector<double>* mtScat = new std::vector<double>();
  
  int sgCalo;
  int ekCaloOK;
  double hitCalo[4];
  double pathCalo;
  double tdcCellA;
  double tdcCellB;
  double adcCellA;
  double adcCellB;
  double posCell[4];
  
  double pathReco;
  double betaReco;
  double ekReco;
  double ptotReco;
  double pReco[4];
  
  tPrimaryNeut.Branch("eventID",&eventID,"eventID/I");
  tPrimaryNeut.Branch("runID",&runID,"runID/I");
  tPrimaryNeut.Branch("pnu",pnu,"pnu[4]/D");;
  tPrimaryNeut.Branch("vtx",vtx,"vtx[4]/D");;
  tPrimaryNeut.Branch("p",p,"p[4]/D");
  tPrimaryNeut.Branch("ptot",&ptot,"ptot/D");
  tPrimaryNeut.Branch("mass",&mass,"mass/D");
  tPrimaryNeut.Branch("ek",&ek,"ek/D");
  tPrimaryNeut.Branch("ekCalo",&ekCalo,"ekCalo/D");
  tPrimaryNeut.Branch("isSTTScat",&isSTTScat,"isSTTScat/I");
  tPrimaryNeut.Branch("thSTTScat","std::vector<double>",thSTTScat);
  tPrimaryNeut.Branch("deSTTScat","std::vector<double>",deSTTScat);
  tPrimaryNeut.Branch("mtSTTScat","std::vector<double>",mtSTTScat);
  tPrimaryNeut.Branch("sgSTTScat","std::vector<int>",sgSTTScat);
  tPrimaryNeut.Branch("isScat",&isScat,"isScat/I");
  tPrimaryNeut.Branch("thScat","std::vector<double>",thScat);
  tPrimaryNeut.Branch("deScat","std::vector<double>",deScat);
  tPrimaryNeut.Branch("mtScat","std::vector<double>",mtScat);
  tPrimaryNeut.Branch("sgCalo",&sgCalo,"sgCalo/I");
  tPrimaryNeut.Branch("ekCaloOK",&ekCaloOK,"ekCaloOK/I");
  tPrimaryNeut.Branch("hitCalo",hitCalo,"hitCalo[4]/D");
  tPrimaryNeut.Branch("pathCalo",&pathCalo,"pathCalo/D");
  tPrimaryNeut.Branch("tdcCellA",&tdcCellA,"tdcCellA/D");
  tPrimaryNeut.Branch("tdcCellB",&tdcCellB,"tdcCellB/D");
  tPrimaryNeut.Branch("adcCellA",&adcCellA,"adcCellA/D");
  tPrimaryNeut.Branch("adcCellB",&adcCellB,"adcCellB/D");
  tPrimaryNeut.Branch("posCell",posCell,"posCell[4]/D");
  tPrimaryNeut.Branch("pathReco",&pathReco,"pathReco/D");
  tPrimaryNeut.Branch("betaReco",&betaReco,"betaReco/D");
  tPrimaryNeut.Branch("ekReco",&ekReco,"ekReco/D");
  tPrimaryNeut.Branch("ptotReco",&ptotReco,"ptotReco/D");
  tPrimaryNeut.Branch("pReco",pReco,"pReco[4]/D");
  
  geo->cd("volWorld_PV/volDetEnclosure_PV_0/volKLOEFULLECALSENSITIVE_EXTTRK_NEWGAP_PV_0/KLOETrackingRegion_volume_PV_0");
  
  double dummyLoc[3];
  double centerKLOE[3];
   
  dummyLoc[0] = 0.;
  dummyLoc[1] = 0.;
  dummyLoc[2] = 0.;
  geo->LocalToMaster(dummyLoc, centerKLOE);
  
  TGeoTube* tub = (TGeoTube*) geo->GetVolume("KLOETrackingRegion_volume_PV")->GetShape();
  TGeoTrd2* trd = (TGeoTrd2*) geo->GetVolume("KLOEBarrelECAL_0_volume_PV")->GetShape();
  
  const double rcalo = tub->GetRmax();
  const double rtrak = rcalo - rgap;
  const double xtrak = trd->GetDy1();
  
  const int nFile = InputFiles->GetEntries();
  
  int* nFileEvents = new int[nFile+1];
  nFileEvents[0] = 0;
    
  for(int j = 1; j < nFile+1; j++)
  {
    InputFiles->GetEntry(j-1);
    nFileEvents[j] = nFileEvents[j-1] + fileEntries;
  }
  
  for(int i = 0; i < nev; i++)
  {
    t->GetEntry(i);
    InputKinem->GetEntry(i);  
    
    int index = 1;
    
    while(i >= nFileEvents[index])
    {
      index++;
    }
    
    index--;
    
    gRooTracker->GetEntry(nFileEvents[index] + inputEntryNum);
    
    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;
    
    eventID = ev->EventId;
    runID = ev->RunId;
    pnu[0] = StdHepP4[0][0];
    pnu[1] = StdHepP4[0][1];
    pnu[2] = StdHepP4[0][2];
    pnu[3] = StdHepP4[0][3];
    
    // only 1 primary vertex
    TG4PrimaryVertex* pV = &(ev->Primaries[0]);
    
    vtx[0] = pV->Position.X();
    vtx[1] = pV->Position.Y();
    vtx[2] = pV->Position.Z();
    vtx[3] = pV->Position.T();
    
    if(debug)
    {
      std::cout << vtx[0] << " " << vtx[1] << " " << vtx[2] << " " << EvtVtx[0] << " " << EvtVtx[1] << " " << EvtVtx[2] << std::endl;
    }
    
    // check
    if(vtx[0] - EvtVtx[0]*1000)
      std::cout << "error: " << vtx[0] - EvtVtx[0]*1000 << " " << inputFileNum << " " 
                << inputEntryNum << " " << index << " " << nFileEvents[index] + inputEntryNum << std::endl;
    
    int id;
    
    // Loop over primary particles
    for(unsigned int j = 0; j < pV->Particles.size(); j++)
    {
      isSTTScat = 0;
      isScat = 0;
      sgCalo = 0;
      ekCaloOK = 0;
      
      thSTTScat->clear();
      deSTTScat->clear();
      mtSTTScat->clear();
      sgSTTScat->clear();
      
      thScat->clear();
      deScat->clear();
      mtScat->clear();
      
      ekCalo = -999999999.;;
      
      hitCalo[0] = -999999999.;
      hitCalo[1] = -999999999.;
      hitCalo[2] = -999999999.;
      hitCalo[3] = -999999999.;
      
      pathCalo = -999999999.;
      tdcCellA = -999999999.;
      tdcCellB = -999999999.;
      adcCellA = -999999999.;
      adcCellB = -999999999.;
      
      posCell[0] = -999999999.;
      posCell[1] = -999999999.;
      posCell[2] = -999999999.;
      posCell[3] = -999999999.;
      
      pathReco = -999999999.;
      betaReco = -999999999.;
      
      ekReco = -999999999.;
      ptotReco = -999999999.;
      pReco[0] = -999999999.;
      pReco[1] = -999999999.;
      pReco[2] = -999999999.;
      pReco[3] = -999999999.;
    
      // selcet only neutrons by PDG code 
      if(pV->Particles[j].PDGCode == neutPDG)
      {
        TG4PrimaryParticle* part = &(pV->Particles[j]);
        p[0] = part->Momentum.X();
        p[1] = part->Momentum.Y();
        p[2] = part->Momentum.Z();
        p[3] = part->Momentum.T();
        ptot = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
        mass = TMath::Sqrt(p[3]*p[3] - ptot*ptot);
        ek = p[3] - mass;
        
        id = part->TrackId;
      
        // Loop over trajectories
        for(unsigned k = 0; k < ev->Trajectories.size(); k++)
        {
          // select neutrons trajectory by TrackId
          if(ev->Trajectories[k].TrackId == id)
          {
            TG4Trajectory* tr = &(ev->Trajectories[k]);
            
            // Loop over neutron trajectory
            for(unsigned m = 1; m < tr->Points.size(); m++)
            {
              TVector3* p1 = &(tr->Points[m-1].Momentum);
              TVector3* p2 = &(tr->Points[m].Momentum);
              
              double dl1 = p1->Mag(); 
              double dl2 = p2->Mag();
              double prod = p1->Dot(*p2);
              
              if(dl1 == 0. || dl2 == 0.)
              {
                if(m < tr->Points.size() - 1)
                {
                  double dl1_x = tr->Points[m-1].Position.X() - tr->Points[m].Position.X();
                  double dl1_y = tr->Points[m-1].Position.Y() - tr->Points[m].Position.Y();
                  double dl1_z = tr->Points[m-1].Position.Z() - tr->Points[m].Position.Z();
                  
                  double dl2_x = tr->Points[m].Position.X() - tr->Points[m+1].Position.X();
                  double dl2_y = tr->Points[m].Position.Y() - tr->Points[m+1].Position.Y();
                  double dl2_z = tr->Points[m].Position.Z() - tr->Points[m+1].Position.Z();
                  
                  dl1 = TMath::Sqrt(dl1_x*dl1_x+dl1_y*dl1_y+dl1_z*dl1_z);
                  dl1 = TMath::Sqrt(dl2_x*dl2_x+dl2_y*dl2_y+dl1_z*dl2_z);
                  
                  prod = dl1_x*dl2_x + dl1_y*dl2_y + dl1_z*dl2_z;
                }
                else
                {
                  dl1 = 1.;
                  dl1 = 1.;
                  
                  prod = 0;
                }
              }
              
              double angle = TMath::ACos(prod/(dl1*dl2));
              
              // check if neutrons scatters by angle between following directions
              if(TMath::Abs(angle) > 0.000001)
              {
                isScat = 1;
                
                double e1 = TMath::Sqrt(dl1*dl1 + mass*mass);
                double e2 = TMath::Sqrt(dl2*dl2 + mass*mass);
                
                /*
                Energy   conservation: E1i + M2 = E1f + E2f -> E2f = E1i + M2 - E1f
                Momentum conservation: p1i = p1f + p2f      -> p2f = p1i - p2f
                
                M2**2 = E1i**2 + M2**2 + E1f**2 + 2 * M2 * E1i - 2 * E1i * E1f - 2 * M2 * E1f - p1i**2 - p1f**2 + 2 * p1i * p1f * cos
                    0 = 2 * M1**2 + 2 * M2 * E1i - 2 * E1i * E1f - 2 * M2 * E1f + 2 * p1i * p1f * cos  
                
                --------> M2 = (E1i * E1f - M1**2 - p1i * p1f * cos ) / ( E1i - E1f ) 
                */
        
                double mt = (e1 * e2 - mass * mass - prod) / (e1 - e2);
                
                thScat->push_back(angle);
                deScat->push_back(e1 - e2);
                mtScat->push_back(mt);
                
              
                double x = tr->Points[m].Position.X();
                double y = tr->Points[m].Position.Y();
                double z = tr->Points[m].Position.Z();
                
                double r = TMath::Sqrt((x - centerKLOE[0])*(x - centerKLOE[0]) + 
                                       (y - centerKLOE[1])*(y - centerKLOE[1]) +
                                       (z - centerKLOE[2])*(z - centerKLOE[2]));
                  
                // check if neutrons scatters in STT 
                if(r <= rtrak && TMath::Abs(x) < xtrak)
                {                  
                  isSTTScat = 1;
                  thSTTScat->push_back(angle);
                  deSTTScat->push_back(e1 - e2);
                  mtSTTScat->push_back(mt);
                  sgSTTScat->push_back(0);
                  
                  // Loop over trcjectories 
                  for(unsigned int l = 0; l < ev->Trajectories.size(); l++)
                  {
                    // Search for neutrons daughter
                    if(ev->Trajectories[l].ParentId == id)
                    {
                      // check daughter origin corresponds to scatter point
                      if(x == ev->Trajectories[l].Points[0].Position.X() &&
                         y == ev->Trajectories[l].Points[0].Position.Y() &&
                         z == ev->Trajectories[l].Points[0].Position.Z())
                      {
                        if(debug)
                        {
                          std::cout << "[" << i << "]: " << y << " - " << ev->Trajectories[l].Points[0].Position.Y() << "  " 
                                                         << z << " - " << ev->Trajectories[l].Points[0].Position.Z() << "  "
                                                         << ev->Trajectories[l].PDGCode << "  "
                                                         << ev->Trajectories[l].TrackId << "  "
                                                         << ev->Trajectories[l].Points[0].Momentum.Mag() << std::endl; 
                          
                          for(unsigned int kk = 0; kk < ev->Trajectories[l].Points.size(); kk++)
                          {
                            std::cout << "[" << kk << "]: " << ev->Trajectories[l].Points[kk].Position.X() << "  " 
                                                            << ev->Trajectories[l].Points[kk].Position.Y() << "  "
                                                            << ev->Trajectories[l].Points[kk].Position.Z() << std::endl;
                          }
                        }
                        
                        // Loop over hit
                        for(unsigned int ii = 0; ii < ev->SegmentDetectors[hitSttName].size(); ii++)
                        {
                          // check if neutron daughter leaves hit on STT
                          if(ev->SegmentDetectors[hitSttName][ii].PrimaryId == ev->Trajectories[l].TrackId)
                          {
                            sgSTTScat->back() = 1;
                            
                            if(debug)
                              std::cout << "found: " << ev->SegmentDetectors[hitSttName][ii].Start.X() << " " 
                                                     << ev->SegmentDetectors[hitSttName][ii].Start.Y() << " " 
                                                     << ev->SegmentDetectors[hitSttName][ii].Start.Z() << std::endl; 
                          }                            
                        }
                      }
                    }
                  }
                }
              }
              
              // evaluate r1 and r2 of the neutron step
              double dx1 = tr->Points[m-1].Position.X() - centerKLOE[0];
              double dy1 = tr->Points[m-1].Position.Y() - centerKLOE[1];
              double dz1 = tr->Points[m-1].Position.Z() - centerKLOE[2];
              
              double dx2 = tr->Points[m].Position.X() - centerKLOE[0];
              double dy2 = tr->Points[m].Position.Y() - centerKLOE[1];
              double dz2 = tr->Points[m].Position.Z() - centerKLOE[2]; 
              
              double r1 = TMath::Sqrt(dy1*dy1+dz1*dz1);
              double r2 = TMath::Sqrt(dy2*dy2+dz2*dz2);
              
              double xx = (dx2 - dx1) / (r2 - r1) * (rcalo - r1) + dx1;
              
              if(debug)
              {
                std::cout << "index: " << m << " of " << tr->Points.size() << std::endl;
                std::cout << "1) " << tr->Points[m-1].Position.X() << " " << tr->Points[m-1].Position.Y() << " " << tr->Points[m-1].Position.Z() << std::endl;
                std::cout << "2) " << tr->Points[m].Position.X() << " " << tr->Points[m].Position.Y() << " " << tr->Points[m].Position.Z() << std::endl;
                std::cout << "r1: " << r1 << "\tr2: " << r2 << "\t" << rcalo << "\t" << xtrak << std::endl;
                std::cout << "mom: " << tr->Points[m].Momentum.X() << " " << tr->Points[m].Momentum.Y() << " " << tr->Points[m].Momentum.Z() << std::endl;
              }
              
              // get neutron momentum when leaves tracking reagion (before entering Calo)
              if(r1 < rcalo && r2 >= rcalo && TMath::Abs(xx) < xtrak)
              {
                double mom = tr->Points[m-1].Momentum.Mag2();
                ekCalo = TMath::Sqrt(mom + mass*mass) - mass;
                ekCaloOK = 1;
              }
            }
          }
        }
        
        // First hit of neutron in Calo
        TG4HitSegment cellFirstHit;
        cellFirstHit.Start.SetT(DBL_MAX); 
        
        // First hit Calo cell
        int cellFirst = -99999999;
        
        // Loop over cells hit
        std::vector<TG4HitSegment>* vhit = &(ev->SegmentDetectors[hitCalName]);
        for(std::map<int, std::vector<int> >::iterator it = list_hit->begin(); it != list_hit->end(); ++it)
        {
          for(unsigned int k = 0; k < it->second.size(); k++)
          {
            if(vhit->at(it->second.at(k)).PrimaryId == id)
            {
              if(vhit->at(it->second.at(k)).Start.T() < cellFirstHit.Start.T())
              {
                if(adc->find(it->first) != adc->end() && adc->find(-1*it->first) != adc->end())
                {
                  cellFirstHit = vhit->at(it->second.at(k));
                  cellFirst = it->first;
                }
              }
            } 
          }
        }
        
        sgCalo = 0;
        if(cellFirst != -99999999)
        {
          sgCalo = 1;
          
          hitCalo[0] = 0.5 * (cellFirstHit.Start.X() + cellFirstHit.Stop.X());
          hitCalo[1] = 0.5 * (cellFirstHit.Start.Y() + cellFirstHit.Stop.Y());
          hitCalo[2] = 0.5 * (cellFirstHit.Start.Z() + cellFirstHit.Stop.Z());
          hitCalo[3] = 0.5 * (cellFirstHit.Start.T() + cellFirstHit.Stop.T());
          
          double pathX = hitCalo[0] - vtx[0];
          double pathY = hitCalo[1] - vtx[1];
          double pathZ = hitCalo[2] - vtx[2];
          
          pathCalo = TMath::Sqrt(pathX*pathX+pathY*pathY+pathZ*pathZ);
          
          int cid = cellFirst >= 0 ? cellFirst : -cellFirst;
          
          tdcCellA = tdc->at(cid);
          tdcCellB = tdc->at(-cid);
          adcCellA = adc->at(cid);
          adcCellB = adc->at(-cid);
          
          int celID = cid % 100;
          int layID = int(cid/100) % 10;
          int modID = int(cid/1000);
          
          double ylay[] = {93., 49., 5., -39., -88.};
          double dy = 230;
          double bmax = 262.55;
          double bmin = 292.85;
          
          double local[3];
          double master[3];
          
          local[2] = ylay[layID];
          local[1] = 0.;
          local[0] = ((bmax - bmin) / dy * (local[2] + dy*0.5) + bmin) * ((celID + 0.5) / 12 - 0.5 );
          
          geo->cd(TString::Format(path_template,modID).Data());
          
          geo->LocalToMaster(local, master);
          
          posCell[1] = master[1];
          posCell[2] = master[2];
          
          /*
          t1 = t + (L/2 - x)/v
          t2 = t + (L/2 + x)/v
          
          t1 - t2 = t + (L/2 - x)/v - t - (L/2 + x)/v = (L/2 - x)/v - (L/2 + x)/v = (L/2 - x - L/2 - x)/v = -2x/v
          */
          
          posCell[3] = (tdcCellA + tdcCellB) * 0.5 - tconst;
          posCell[0] = (tdcCellA - tdcCellB) * 0.5 / vel * 1000.;
          pathReco = TMath::Sqrt(TMath::Power(posCell[0] - vtx[0],2) + 
                    TMath::Power(posCell[1] - vtx[1],2) +
                    TMath::Power(posCell[2] - vtx[2],2) );
                    
          betaReco = pathReco/(posCell[3] - vtx[3])/c;
          
          double gamma = 1./TMath::Sqrt(1. - betaReco*betaReco);
          
          ekReco = mass * (gamma - 1.) ;
          
          ptotReco = mass * gamma * betaReco;
          
          pReco[0] = ptotReco/pathReco * (posCell[0] - vtx[0]);
          pReco[1] = ptotReco/pathReco * (posCell[1] - vtx[1]);
          pReco[2] = ptotReco/pathReco * (posCell[2] - vtx[2]);
          pReco[3] = mass * gamma;  
        }
        
        tPrimaryNeut.Fill();
      }      
    }
  }

  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;
  
  fout.cd();
  tPrimaryNeut.Write();
  fout.Close();
  
}

