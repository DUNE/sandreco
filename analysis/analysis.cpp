#include <TChain.h>
#include <TFile.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include "/mnt/nas01/users/mtenti/wd/analysis/KLOEcal/loader/loader.C"

#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4Event.h"
#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4HitSegment.h"

#include <iostream>

  
const double k = 0.299792458;
const double B = 0.6;

TDatabasePDG db; 

void reset(particle& p)
{
  p.primary = false;
  p.pdg = -999;
  p.tid = -999;
  p.mass = -999.;
  p.pxtrue = -999.;
  p.pytrue = -999.;
  p.pztrue = -999.;
  p.Etrue = -999.;
  p.pxreco = -999.;
  p.pyreco = -999.;
  p.pzreco = -999.;
  p.Ereco = -999.;
  p.x0dig = -999.;
  p.y0dig = -999.;
  p.z0dig = -999.;
  p.t0dig = -999.;
  p.x0trj = -999.;
  p.y0trj = -999.;
  p.z0trj = -999.;
  p.t0trj = -999.;
}

bool IsPrimary(TG4Event* ev, int tid)
{
  for(unsigned int i = 0; i < ev->Primaries.at(0).Particles.size(); i++)
  {
    if(ev->Primaries.at(0).Particles.at(i).TrackId == tid)
      return true;
  }
  return false;
}

void FillParticleInfo(TG4Event* ev, particle& p)
{
  p.primary = IsPrimary(ev, p.tid) == true ? 1 : 0;  
  
  for(unsigned int i = 0; i < ev->Trajectories.size(); i++)
  {
    if(ev->Trajectories.at(i).TrackId == p.tid)
    {
      p.pdg = ev->Trajectories.at(i).PDGCode;
      TParticlePDG* part = db.GetParticle(p.pdg);
      if(part != 0) p.mass = part->Mass();
      
      p.pxtrue = ev->Trajectories.at(i).InitialMomentum.X();
      p.pytrue = ev->Trajectories.at(i).InitialMomentum.Y();
      p.pztrue = ev->Trajectories.at(i).InitialMomentum.Z();
      p.Etrue = ev->Trajectories.at(i).InitialMomentum.T();
      
      p.x0trj = ev->Trajectories.at(i).Points.at(0).Position.X();
      p.y0trj = ev->Trajectories.at(i).Points.at(0).Position.Y();
      p.z0trj = ev->Trajectories.at(i).Points.at(0).Position.Z();
      p.t0trj = ev->Trajectories.at(i).Points.at(0).Position.T();
    }
  }
}

void FillTrackInfo(const track& tr, particle& p)
{
  double mom_yz = k * tr.r * B;
  double dz = (tr.zc - tr.z0);
  double dy = (tr.yc - tr.y0);
  double ang_yz = TMath::ATan2(dz, -dy); 
  double ang_x = 0.5 * TMath::Pi() - TMath::ATan(1./tr.b);
  
  p.pxreco = mom_yz * TMath::Tan(ang_x);
  p.pyreco = mom_yz * TMath::Sin(ang_yz);
  p.pzreco = mom_yz * TMath::Cos(ang_yz);
  p.Ereco = TMath::Sqrt(p.pxreco*p.pxreco + p.pyreco*p.pyreco + p.pzreco*p.pzreco + p.mass*p.mass);
  p.x0dig = tr.x0;
  p.y0dig = tr.y0;
  p.z0dig = tr.z0;
  p.t0dig = tr.t0;
}

void analysis(const char* fReco, const char* fTrueMC, const char* fOut)
{
  TChain tReco("tReco");
  TChain tTrueMC("EDepSimEvents");
  tReco.Add(fReco);
  tTrueMC.Add(fTrueMC);
  
  tReco.AddFriend(&tTrueMC);
  
  TChain* t = &tReco;
  
  std::vector<track>* vec_tr = new std::vector<track>;
  std::vector<cluster>* vec_cl = new std::vector<cluster>;
  
  TG4Event* ev = new TG4Event;
  t->SetBranchAddress("Event",&ev);
  t->SetBranchAddress("track",&vec_tr);
  t->SetBranchAddress("cluster",&vec_cl);
  
  std::vector<particle> vec_part;
    
  TFile fout(fOut,"RECREATE");
  TTree tout("tParticle","tParticle");
  tout.Branch("particle","std::vector<particle>",&vec_part);
    
  const int nev = t->GetEntries();
  
  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;
  
  for(int i = 0; i < nev; i++)
  {
    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;
    
    t->GetEntry(i);
    vec_part.clear();
    
    for(unsigned int j = 0; j < vec_tr->size(); j++)
    {
      if(vec_tr->at(j).ret_ln == 0 && vec_tr->at(j).ret_cr == 0)
      {
        particle p;
        reset(p);
        
        p.tid = vec_tr->at(j).tid;
        
        FillParticleInfo(ev, p);
        FillTrackInfo(vec_tr->at(j), p);
        
        vec_part.push_back(p);
      }
    }    
    tout.Fill();
  }
  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;
  
  fout.cd();
  tout.Write();
  fout.Close();
}