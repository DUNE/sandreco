#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4Event.h"
#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4HitSegment.h"

#include <iostream>
#include <sstream>

#include <TFile.h>
#include <TTree.h>

TFile* f;
TTree* t;
TG4Event* ev;

void init(const char* fname, const char* tname, const char* bname)
{
  f = new TFile(fname);
  t = (TTree*) f->Get(tname);
  t->SetBranchAddress(bname,&ev);
}

void dumpPrt(int entry)
{
  for(unsigned int i = 0; i < ev->Primaries[entry].Particles.size(); i++)
  {
    std::cout << entry << "." << i <<  ")" << std::endl; 
    std::cout << "TID: " << ev->Primaries[entry].Particles[i].TrackId << 
                 " NAM: " << ev->Primaries[entry].Particles[i].Name << 
                 " PDG: " << ev->Primaries[entry].Particles[i].PDGCode << std::endl;
    std::cout << "MOM: " << ev->Primaries[entry].Particles[i].Momentum.X() << " " 
                         << ev->Primaries[entry].Particles[i].Momentum.Y() << " " 
                         << ev->Primaries[entry].Particles[i].Momentum.Z() << " " 
                         << ev->Primaries[entry].Particles[i].Momentum.T() << std::endl;
  }
}

void dumpVtx()
{
  for(unsigned int i = 0; i < ev->Primaries.size(); i++)
  {
    std::cout << i << ")" << std::endl;
    std::cout << "GEN: " << ev->Primaries[i].GeneratorName.c_str() << std::endl;
    std::cout << "REA: " << ev->Primaries[i].Reaction.c_str() << std::endl;
    std::cout << "FIL: " << ev->Primaries[i].Filename.c_str() << std::endl;
    std::cout << "POS: " << ev->Primaries[i].Position.X() << " " << ev->Primaries[i].Position.Y() << " " << ev->Primaries[i].Position.Z() << std::endl;
    std::cout << "INT: " << ev->Primaries[i].InteractionNumber 
              << " XSC: " << ev->Primaries[i].CrossSection 
              << " DFX: " << ev->Primaries[i].DiffCrossSection << std::endl;
    std::cout << "WGT: " << ev->Primaries[i].Weight << " PRB: " << ev->Primaries[i].Probability << " PRT: " << ev->Primaries[i].Particles.size() << std::endl;
    
    dumpPrt(i);
  }
}

void dumpPnt(int entry)
{
  for(unsigned int i = 0; i < ev->Trajectories[entry].Points.size(); i++)
  {
    std::cout << entry << "." << i << ")" << std::endl;
    std::cout << "POS: " << ev->Trajectories[entry].Points[i].Position.X() << " " 
                         << ev->Trajectories[entry].Points[i].Position.Y() << " "
                         << ev->Trajectories[entry].Points[i].Position.Z() << " "
                         << ev->Trajectories[entry].Points[i].Position.T() << std::endl;
    std::cout << "MOM: " << ev->Trajectories[entry].Points[i].Momentum.X() << " " 
                         << ev->Trajectories[entry].Points[i].Momentum.Y() << " " 
                         << ev->Trajectories[entry].Points[i].Momentum.Z() << std::endl;
    std::cout << "PRC: " << ev->Trajectories[entry].Points[i].Process << " "
                 " SPR: " << ev->Trajectories[entry].Points[i].Subprocess << std::endl;
  }
}

void dumpTrj()
{
  for(unsigned int i = 0; i < ev->Trajectories.size(); i++)
  {
    std::cout << i << ")" << std::endl;
    std::cout << "TID: " << ev->Trajectories[i].TrackId << 
                 " PID: " << ev->Trajectories[i].ParentId <<
                 " NAM: " << ev->Trajectories[i].Name.c_str() << 
                 " PDG: " << ev->Trajectories[i].PDGCode << 
                 " NPT: " << ev->Trajectories[i].Points.size() << std::endl;
    std::cout << "MOM: " << ev->Trajectories[i].InitialMomentum.X() << " " 
                         << ev->Trajectories[i].InitialMomentum.Y() << " " 
                         << ev->Trajectories[i].InitialMomentum.Z() << " " 
                         << ev->Trajectories[i].InitialMomentum.T() << std::endl;
    dumpPnt(i);
  }
}

void dumpHit()
{
  for(map<string,vector<TG4HitSegment> >::iterator it = ev->SegmentDetectors.begin(); it != ev->SegmentDetectors.end(); ++it)
  {
    std::cout << "----------------------------------------------------------------------------" << std::endl;  
    std::cout << "HIT: " << it->first.c_str() << "\t[ " << it->second.size() << " ]" << std::endl;
    std::cout << "----------------------------------------------------------------------------" << std::endl;  
    
    for(unsigned int j = 0; j < it->second.size(); j++)
    {
      std::cout << j << ")" << std::endl;
      std::cout << "STR: " << it->second.at(j).Start.X() << " " << it->second.at(j).Start.Y() << " " << it->second.at(j).Start.Z() << std::endl;
      std::cout << "STP: " << it->second.at(j).Stop.X() << " " << it->second.at(j).Stop.Y() << " " << it->second.at(j).Stop.Z() << std::endl;
      std::cout << "PRI: " << it->second.at(j).PrimaryId << 
                  " DE : " << it->second.at(j).EnergyDeposit << 
                  " DE2: " << it->second.at(j).SecondaryDeposit << 
                  " TRL: " << it->second.at(j).TrackLength << std::endl;
      std::ostringstream str;
      str << "CRB: ";
      for(unsigned int k = 0; k < it->second.at(j).Contrib.size(); k++)
      {
        str << it->second.at(j).Contrib[k] << ", ";
      }
      
      std::cout << str.str() << std::endl;
    }
  }
}

void dumpEvent()
{
  std::cout << "==========================================================================" << std::endl;
  std::cout << "run: " << ev->RunId << " event: " << ev->EventId << " vtx: " << ev->Primaries.size() << " trj: " << ev->Trajectories.size() << std::endl;
  std::cout << "==========================================================================" << std::endl;
  
  dumpVtx();
  dumpTrj();
  dumpHit();
}

void dump(int entry)
{
  t->GetEntry(entry);
  dumpEvent();
}