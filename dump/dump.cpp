#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4Event.h"
#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4HitSegment.h"

#include <iostream>
#include <sstream>

#include <TFile.h>
#include <TTree.h>

namespace ns_dump {
  TFile* f;
  TTree* t;
  TG4Event* ev;
}

using namespace ns_dump;

void init(const char* fname, const char* tname = "EDepSimEvents", const char* bname = "Event")
{
  f = new TFile(fname);
  t = (TTree*) f->Get(tname);
  t->SetBranchAddress(bname,&ev);
}

void dumpPrt()
{
  std::cout << std::setw(7)  << std::right << "vertex"  
            << std::setw(7)  << std::right << "track"   
            << std::setw(7)  << std::right << "id"    
            << std::setw(20) << std::right << "name"   
            << std::setw(10) << std::right << "pdg"    
            << std::setw(10) << std::right << "px"    
            << std::setw(10) << std::right << "py"    
            << std::setw(10) << std::right << "pz"    
            << std::setw(10) << std::right << "E"      
            << std::endl;
                              
  for(unsigned int j = 0; j < ev->Primaries.size(); j++)
  {
    for(unsigned int i = 0; i < ev->Primaries[j].Particles.size(); i++)
    {
      std::cout << std::setw(7)  << std::right << j 
                << std::setw(7)  << std::right << i 
                << std::setw(7)  << std::right << ev->Primaries[j].Particles[i].TrackId 
                << std::setw(20) << std::right << ev->Primaries[j].Particles[i].Name 
                << std::setw(10) << std::right << ev->Primaries[j].Particles[i].PDGCode 
                << std::setw(10) << std::right << ev->Primaries[j].Particles[i].Momentum.X()
                << std::setw(10) << std::right << ev->Primaries[j].Particles[i].Momentum.Y()
                << std::setw(10) << std::right << ev->Primaries[j].Particles[i].Momentum.Z()
                << std::setw(10) << std::right << ev->Primaries[j].Particles[i].Momentum.T() 
                << std::endl;
    }
  }
  std::cout << "==========================================================================" << std::endl;
}

void dumpVtx()
{
  std::cout << std::setw(7)   << std::right << "entry" 
            << std::setw(30)  << std::right << "generator" 
            << std::setw(100) << std::right << "reaction" 
            << std::setw(100) << std::right << "file" 
            << std::setw(10)  << std::right << "X"  
            << std::setw(10)  << std::right << "Y"  
            << std::setw(10)  << std::right << "Z" 
            << std::setw(10)  << std::right << "T"  
            << std::setw(10)  << std::right << "inter"
            << std::setw(20)  << std::right << "XS" 
            << std::setw(20)  << std::right << "diffXS" 
            << std::setw(20)  << std::right << "weigth" 
            << std::setw(20)  << std::right << "prob" 
            << std::setw(10)  << std::right << "npart"
            << std::endl;

  for(unsigned int i = 0; i < ev->Primaries.size(); i++)
  {
    std::cout << std::setw(7)   << std::right << i
              << std::setw(30)  << std::right << ev->Primaries[i].GeneratorName.c_str() 
              << std::setw(100) << std::right << ev->Primaries[i].Reaction.c_str() 
              << std::setw(100) << std::right << ev->Primaries[i].Filename.c_str() 
              << std::setw(10)  << std::right << ev->Primaries[i].Position.X() 
              << std::setw(10)  << std::right << ev->Primaries[i].Position.Y() 
              << std::setw(10)  << std::right << ev->Primaries[i].Position.Z() 
              << std::setw(10)  << std::right << ev->Primaries[i].Position.T() 
              << std::setw(10)  << std::right << ev->Primaries[i].InteractionNumber 
              << std::setw(20)  << std::right << ev->Primaries[i].CrossSection 
              << std::setw(20)  << std::right << ev->Primaries[i].DiffCrossSection 
              << std::setw(20)  << std::right << ev->Primaries[i].Weight 
              << std::setw(20)  << std::right << ev->Primaries[i].Probability 
              << std::setw(10)  << std::right << ev->Primaries[i].Particles.size() 
              << std::endl;
  }
  std::cout << "==========================================================================" << std::endl;
}

void dumpTrj()
{
  std::cout << std::setw(7)  << std::right << "entry" 
            << std::setw(7)  << std::right << "id" 
            << std::setw(7)  << std::right << "pid"
            << std::setw(10) << std::right << "name" 
            << std::setw(10) << std::right << "pdg" 
            << std::setw(10) << std::right << "npoints" 
            << std::setw(10) << std::right << "ipx" 
            << std::setw(10) << std::right << "ipy" 
            << std::setw(10) << std::right << "ipz" 
            << std::setw(10) << std::right << "iE" 
            << std::setw(7)  << std::right << "point" 
            << std::setw(10) << std::right << "X" 
            << std::setw(10) << std::right << "Y" 
            << std::setw(10) << std::right << "Z" 
            << std::setw(10) << std::right << "T" 
            << std::setw(10) << std::right << "px" 
            << std::setw(10) << std::right << "py" 
            << std::setw(10) << std::right << "pz" 
            << std::setw(5)  << std::right << "proc" 
            << std::setw(5)  << std::right << "subp"
            << std::endl;
                
  for(unsigned int i = 0; i < ev->Trajectories.size(); i++)
  {
    for(unsigned int j = 0; j < ev->Trajectories[i].Points.size(); j++)
    {
      std::cout << std::setw(7) << i 
                << std::setw(7)  << std::right << ev->Trajectories[i].TrackId 
                << std::setw(7)  << std::right << ev->Trajectories[i].ParentId 
                << std::setw(10) << std::right << ev->Trajectories[i].Name.c_str() 
                << std::setw(10) << std::right << ev->Trajectories[i].PDGCode 
                << std::setw(10) << std::right << ev->Trajectories[i].Points.size() 
                << std::setw(10) << std::right << ev->Trajectories[i].InitialMomentum.X() 
                << std::setw(10) << std::right << ev->Trajectories[i].InitialMomentum.Y() 
                << std::setw(10) << std::right << ev->Trajectories[i].InitialMomentum.Z() 
                << std::setw(10) << std::right << ev->Trajectories[i].InitialMomentum.T() 
                << std::setw(7)  << std::right << j 
                << std::setw(10) << std::right << ev->Trajectories[i].Points[j].Position.X() 
                << std::setw(10) << std::right << ev->Trajectories[i].Points[j].Position.Y() 
                << std::setw(10) << std::right << ev->Trajectories[i].Points[j].Position.Z() 
                << std::setw(10) << std::right << ev->Trajectories[i].Points[j].Position.T() 
                << std::setw(10) << std::right << ev->Trajectories[i].Points[j].Momentum.X() 
                << std::setw(10) << std::right << ev->Trajectories[i].Points[j].Momentum.Y() 
                << std::setw(10) << std::right << ev->Trajectories[i].Points[j].Momentum.Z() 
                << std::setw(5) <<  std::right << ev->Trajectories[i].Points[j].Process 
                << std::setw(5) <<  std::right << ev->Trajectories[i].Points[j].Subprocess 
                << std::endl;
    }
    std::cout << "----------------------------------------------------------------------------" << std::endl;
  }
  std::cout << "==========================================================================" << std::endl;
}

void dumpHit()
{
  std::cout << std::setw(7)  << std::right << "index"
            << std::setw(20) << std::right << "type" 
            << std::setw(10) << std::right << "X1"   
            << std::setw(10) << std::right << "Y1"   
            << std::setw(10) << std::right << "Z1"   
            << std::setw(10) << std::right << "T1"  
            << std::setw(10) << std::right << "X2"   
            << std::setw(10) << std::right << "Y2"   
            << std::setw(10) << std::right << "Z2"   
            << std::setw(10) << std::right << "T2"    
            << std::setw(7)  << std::right << "pid"   
            << std::setw(15) << std::right << "dE"    
            << std::setw(15) << std::right << "dE2"   
            << std::setw(15) << std::right << "L" 
            << std::setw(20) << std::right << "Contrib"
            << std::endl;
                  
  for(map<string,vector<TG4HitSegment> >::iterator it = ev->SegmentDetectors.begin(); it != ev->SegmentDetectors.end(); ++it)
  {
    for(unsigned int j = 0; j < it->second.size(); j++)
    {
      std::cout << std::setw(7)  << std::right << j
                << std::setw(20) << std::right << it->first 
                << std::setw(10) << std::right << it->second.at(j).Start.X()
                << std::setw(10) << std::right << it->second.at(j).Start.Y()
                << std::setw(10) << std::right << it->second.at(j).Start.Z()
                << std::setw(10) << std::right << it->second.at(j).Start.T() 
                << std::setw(10) << std::right << it->second.at(j).Stop.X()  
                << std::setw(10) << std::right << it->second.at(j).Stop.Y()  
                << std::setw(10) << std::right << it->second.at(j).Stop.Z()  
                << std::setw(10) << std::right << it->second.at(j).Stop.T()
                << std::setw(7)  << std::right << it->second.at(j).PrimaryId   
                << std::setw(15) << std::right << it->second.at(j).EnergyDeposit    
                << std::setw(15) << std::right << it->second.at(j).SecondaryDeposit   
                << std::setw(15) << std::right << it->second.at(j).TrackLength
                << "\t";
                
                for(unsigned int k = 0; k < it->second.at(j).Contrib.size(); k++)
                {
                  std::cout << it->second.at(j).Contrib[k] << ", ";
                }
                
                std::cout << std::endl;
    }
    std::cout << "----------------------------------------------------------------------------" << std::endl;
  }
  std::cout << "==========================================================================" << std::endl;
}

void dumpEvent()
{
  std::cout << "==========================================================================" << std::endl;
  std::cout << "run: " << std::setw(6) << ev->RunId << 
            " event: " << std::setw(6) << ev->EventId << 
         " vertices: " << std::setw(6) << ev->Primaries.size() << 
           " tracks: " << std::setw(6) << ev->Trajectories.size() << std::endl;
  std::cout << "==========================================================================" << std::endl;
  
  dumpVtx();
  dumpPrt();
  dumpTrj();
  dumpHit();
}

void dump(int entry)
{
  t->GetEntry(entry);
  dumpEvent();
}