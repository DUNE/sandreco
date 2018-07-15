#include <map>
#include <iostream>
#include <vector>

#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4Event.h"
#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4HitSegment.h"

#include <TFile.h>
#include <TTree.h>
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoTrd2.h>
#include <TChain.h>

#include "/mnt/nas01/users/mtenti/wd/analysis/KLOEcal/loader/loader.C"

void STTDigi()
{
  std::cout << "DigitizeSTT(const char* finname, const char* foutname)" << std::endl;
  std::cout << "finname could contain wild card" << std::endl;
}

double mindist(double s1x, double s1y, double s1z,
               double s2x, double s2y, double s2z,
               double px, double py, double pz)
{
    double segmod = (s1x - s2x)*(s1x - s2x)
                    +(s1y - s2y)*(s1y - s2y)
                    +(s1z - s2z)*(s1z - s2z);

    double prod = (px - s1x)*(s2x - s1x)
                  +(py - s1y)*(s2y - s1y)
                  +(pz - s1z)*(s2z - s1z);

    double t = std::min(std::max(prod/segmod,0.),1.);

    double s3x = s1x + (s2x - s1x) * t;
    double s3y = s1y + (s2y - s1y) * t;
    double s3z = s1z + (s2z - s1z) * t;
    
    return sqrt((px - s3x)*(px - s3x)
                +(py - s3y)*(py - s3y)
                +(pz - s3z)*(pz - s3z));
}

double angle(double x1, double y1, double z1,
             double x2, double y2, double z2)
{
    double prod = x1*x2+y1*y2+z1*z2;
    double mag1 = sqrt(x1*x1+y1*y1+z1*z1);
    double mag2 = sqrt(x2*x2+y2*y2+z2*z2);

    return TMath::ACos(prod/(mag1*mag2));
}

bool ishitok(TG4Event* ev, TG4HitSegment hit,
             double postol = 5., double angtol = 0.3)
{
    double x = 0.5*(hit.Start.X()+hit.Stop.X());
    double y = 0.5*(hit.Start.Y()+hit.Stop.Y());
    double z = 0.5*(hit.Start.Z()+hit.Stop.Z());

    std::vector<double> dpos;
    std::vector<double> dang;
    for(unsigned int jj = 0; jj < ev->Trajectories.size(); jj++)
    {
        if(ev->Trajectories[jj].TrackId == hit.PrimaryId)
        {
            for(unsigned int kk = 0; kk < ev->Trajectories[jj].Points.size()-1; kk++)
            {
                dpos.push_back(mindist(ev->Trajectories[jj].Points[kk].Position.X(),
                                       ev->Trajectories[jj].Points[kk].Position.Y(),
                                       ev->Trajectories[jj].Points[kk].Position.Z(),
                                       ev->Trajectories[jj].Points[kk+1].Position.X(),
                                       ev->Trajectories[jj].Points[kk+1].Position.Y(),
                                       ev->Trajectories[jj].Points[kk+1].Position.Z(),
                                       x,y,z));
                dang.push_back(angle(ev->Trajectories[jj].Points[kk+1].Position.X()-ev->Trajectories[jj].Points[kk].Position.X(),
                                     ev->Trajectories[jj].Points[kk+1].Position.Y()-ev->Trajectories[jj].Points[kk].Position.Y(),
                                     ev->Trajectories[jj].Points[kk+1].Position.Z()-ev->Trajectories[jj].Points[kk].Position.Z(),
                                     hit.Stop.X()-hit.Start.X(),
                                     hit.Stop.Y()-hit.Start.Y(),
                                     hit.Stop.Z()-hit.Start.Z()));
            }
        }
    }
    int index = std::distance(dpos.begin(), std::min_element(dpos.begin(),dpos.end()));

    //std::cout << dpos[index] << " " << dang[index] << std::endl;

    if(dpos[index] > postol || dang[index] > angtol)
        return false;
    else
        return true;
}

void DigitizeSTT(const char* finname, const char* foutname)
{
  
  TChain* t = new TChain("EDepSimEvents","EDepSimEvents");
  t->Add(finname);
  TFile f(t->GetListOfFiles()->At(0)->GetTitle());
  TGeoManager* geo = (TGeoManager*) f.Get("EDepSimGeometry");
  
  TG4Event* ev = new TG4Event;
  
  t->SetBranchAddress("Event",&ev);
  
  std::map<std::string,std::vector<hit> > cluster_map;
  std::vector<digit> digit_vec;
  
  TFile fout(foutname,"RECREATE");
  TTree tstt("tSttDigi","KLOE Stt Digitization");
  tstt.Branch("Stt","std::vector<digit>",&digit_vec);
    
  const int nev = t->GetEntries();
  
  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;

  for(int i = 0; i < nev; i++)
  {
    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;
      
    t->GetEntry(i);
    
    digit_vec.clear();
    cluster_map.clear();
    
    for(unsigned int j = 0; j < ev->SegmentDetectors["StrawTracker"].size(); j++)
    {
      const TG4HitSegment& hseg = ev->SegmentDetectors["StrawTracker"].at(j);
      
      if(!ishitok(ev, hseg))
        continue;
      
      double x = 0.5 * (hseg.Start.X() + hseg.Stop.X());
      double y = 0.5 * (hseg.Start.Y() + hseg.Stop.Y());
      double z = 0.5 * (hseg.Start.Z() + hseg.Stop.Z());
      
      std:string sttname = geo->FindNode(x,y,z)->GetName();
      
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
    
    tstt.Fill();
  }
  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;
    
  fout.cd();
  tstt.Write();
  fout.Close();
  
  delete t;
  f.Close();
}