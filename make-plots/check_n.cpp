#include <iostream>

#include <TTree.h>
#include <TFile.h>
#include <TGeoManager.h>

#include "/wd/sw/EDEPSIM/edep-sim.binary/include/EDepSim/TG4Event.h"
#include "/wd/sw/EDEPSIM/edep-sim.binary/include/EDepSim/TG4HitSegment.h"

#include "../include/struct.h"

TGeoManager* geo = 0;
TTree* t = 0;
TFile* f = 0;

double m_to_mm = 1000.;
double vlfb = 5.85;
const double k = 0.299792458;
const double vc = k * 1E3; // mm/ns

void CellXYZTE(cell c, double& x, double& y, double& z, double& t, double& e)
{      
  if(c.id < 25000) //Barrel
  {
    x = 0.5 * (c.tdc1 - c.tdc2)/vlfb * m_to_mm + c.x;
    y = c.y;
  }
  else
  {
    x = c.x;
    y = -0.5 * (c.tdc1 - c.tdc2)/vlfb * m_to_mm + c.y;
  }
  z = c.z;
  t = 0.5 * (c.tdc1 + c.tdc2 - vlfb * c.l / m_to_mm );
  e = c.adc1 + c.adc2;
}

double min_dist(double x0, double y0, double z0,
                double x1, double y1, double z1,
                double x2, double y2, double z2)
{
  double k = ((x0-x1)*(x2-x1)+(y0-y1)*(y2-y1)+(z0-z1)*(z2-z1))/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
  
  if(k < 0.)
    k = 0.;
  else if(k > 1.)
    k = 1.;
  
  double xv = x1 + k * (x2 - x1);
  double yv = y1 + k * (y2 - y1);
  double zv = z1 + k * (z2 - z1);
  
  double d = sqrt((x0 - xv)*(x0 - xv)+(y0 - yv)*(y0 - yv)+(z0 - zv)*(z0 - zv));
  
  return d;
}

void check_n(const char* ifile)
{
  f = new TFile(ifile);
  TTree* tEvent = reinterpret_cast<TTree*>(f->Get("tEvent"));
  TTree* tReco = reinterpret_cast<TTree*>(f->Get("tReco"));
  TTree* tDigit = reinterpret_cast<TTree*>(f->Get("tDigit"));
  TTree* tEdep = reinterpret_cast<TTree*>(f->Get("EDepSimEvents"));
  TTree* tGenie = reinterpret_cast<TTree*>(f->Get("gRooTracker"));

  tEvent->AddFriend(tReco);
  tEvent->AddFriend(tDigit);
  tEvent->AddFriend(tEdep);
  tEvent->AddFriend(tGenie);

  t = tEvent;

  geo = reinterpret_cast<TGeoManager*>(f->Get("EDepSimGeometry"));
  
  event* ev = new event;
  particle* p = 0;
  cell* c = 0;
  TG4HitSegment* seg = 0;
  TG4Trajectory* trj = 0;

  int hidx = -1;
  
  TG4Event* g4e = new TG4Event;
  
  double x,y,z,tm,e,beta,dx,dy,dz,dr,dt,bt,br;
  
  t->SetBranchAddress("event",&ev);
  t->SetBranchAddress("Event",&g4e);
  
  int nev = t->GetEntries();

  int err1 = 0;
  int err2 = 0;
  int err3 = 0;
  int err4 = 0;
  int ok = 0;
  int good = 0;
  
  for(int i = 0; i < nev; i++)
  {
    t->GetEntry(i);
    
    unsigned int np = ev->particles.size();
    
    for(unsigned int j = 0; j < np; j++)
    {
      p = &ev->particles.at(j);
    
      if(p->pdg == 2112 && p->primary == 1 && p->Etrue < p->Ereco)
      {
        unsigned int nc = p->cl.cells.size();
        
        trj = &g4e->Trajectories[p->tid];

        if(trj->GetTrackId() != p->tid)
        {
          std::cout << "Error: ID(trj) != ID(part)" << std::endl;
          err1++;
          continue;
        }

        std::vector<double> times(nc);
        double xdummy, ydummy, zdummy, edummy;
        
        for(unsigned int k = 0; k < nc; k++)
        {
          CellXYZTE(p->cl.cells.at(k), xdummy, ydummy, zdummy, times.at(k), edummy);
        }

        int idx_min = std::distance(times.begin(), std::min_element(times.begin(), times.end()));
        
        //for(unsigned int k = 0; k < nc; k++)
        //{
        c = &p->cl.cells.at(idx_min);
        
        CellXYZTE(*c, x, y, z, tm, e);
        
        dx = x - ev->x;
        dy = y - ev->y;
        dz = z - ev->z;
        dt = tm - ev->t;
        
        dr = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
        
        beta = dr / dt / vc;
        
        bt = TMath::Sqrt(p->pxtrue*p->pxtrue+p->pytrue*p->pytrue+p->pztrue*p->pztrue)/p->Etrue;
        br = TMath::Sqrt(p->pxreco*p->pxreco+p->pyreco*p->pyreco+p->pzreco*p->pzreco)/p->Ereco;
        
        if(abs(1-beta/br) > 0.0001)
        {
          std::cout << "Error: beta not correctly reconstructed" << std::endl;
          err2++;
          continue;
        }
        
        //if(abs(1-beta/br) < 0.0001)
        //{
      
        std::cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" <<std::endl;
                
        std::cout << std::setw(3) << "id" << "|" <<
                      std::setprecision(5) << std::setw(7) << "vx" << "|" <<
                      std::setprecision(5) << std::setw(6) << "vy" << "|" <<
                      std::setprecision(5) << std::setw(6) << "vz" << "|" <<
                      std::setw(2) << "vt" << "|" <<
                      std::setw(3) << "pid" << "|" <<
                      std::setw(5) << "cid" << "|" <<
                      std::setw(5) << "cx" << "|" <<
                      std::setw(6) << "cy" << "|" << 
                      std::setw(6) << "cz" << "|" <<
                      std::setw(6) << "cl" << "|" <<
                      std::setw(4) << "adc1" << "|" <<
                      std::setw(4) << "adc2" << "|" <<
                      std::setw(6) << "tdc1" << "|" <<
                      std::setw(6) << "tdc2" << "|" << 
                      std::setprecision(5) << std::setw(7) << "x" << "|" <<  
                      std::setprecision(5) << std::setw(6) << "y" << "|" <<  
                      std::setprecision(5) << std::setw(6) << "z" << "|" <<  
                      std::setprecision(3) << std::setw(4) << "t" << "|" << 
                      std::setw(3) << "e" << "|" <<  
                      std::setprecision(5) << std::setw(7) << "dx" << "|" <<  
                      std::setprecision(5) << std::setw(10) << "dy" << "|" <<  
                      std::setprecision(5) << std::setw(7) << "dz" << "|" <<   
                      std::setprecision(5) << std::setw(7) << "dt" << "|" <<   
                      std::setprecision(5) << std::setw(7) << "dr" << "|" <<  
                      std::setprecision(5) << std::setw(8) << "beta" << "|" <<
                      std::setprecision(5) << std::setw(8) << "br" << "|" <<
                      std::setprecision(5) << std::setw(8) << "bt" << "|" << std::endl;
                      
        std::cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" <<std::endl;

        std::cout << std::setw(3) << i << "|" <<
                      std::setprecision(5) << std::setw(7) << ev->x << "|" <<
                      std::setprecision(5) << std::setw(6) << ev->y << "|" <<
                      std::setprecision(5) << std::setw(6) << ev->z << "|" <<
                      std::setw(2) << ev->t << "|" <<
                      std::setw(3) << p->tid << "|" <<
                      std::setw(5) << c->id << "|" <<
                      std::setw(5) << c->x << "|" <<
                      std::setw(6) << c->y << "|" << 
                      std::setw(6) << c->z << "|" <<
                      std::setw(6) << c->l << "|" <<
                      std::setw(4) << c->adc1 << "|" <<
                      std::setw(4) << c->adc2 << "|" <<
                      std::setw(6) << c->tdc1 << "|" <<
                      std::setw(6) << c->tdc2 << "|" << 
                      std::setprecision(5) << std::setw(7) << x << "|" <<  
                      std::setprecision(5) << std::setw(6) << y << "|" <<  
                      std::setprecision(5) << std::setw(6) << z << "|" <<  
                      std::setprecision(3) << std::setw(4) << tm << "|" << 
                      std::setw(3) << e << "|" <<  
                      std::setprecision(5) << std::setw(7) << dx << "|" <<  
                      std::setprecision(5) << std::setw(10) << dy << "|" <<  
                      std::setprecision(5) << std::setw(7) << dz << "|" <<   
                      std::setprecision(5) << std::setw(7) << dt << "|" <<   
                      std::setprecision(5) << std::setw(7) << dr << "|" <<  
                      std::setprecision(5) << std::setw(8) << beta << "|" <<
                      std::setprecision(5) << std::setw(8) << br << "|" <<
                      std::setprecision(5) << std::setw(8) << bt << "|" << std::endl;
          
          std::cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" <<std::endl;
          
          bool found = false;

          hidx = -999;

          std::sort(c->hindex1.begin(), c->hindex1.end());
          c->hindex1.resize( std::distance(c->hindex1.begin(),std::unique (c->hindex1.begin(), c->hindex1.end())) );

          for(unsigned int l = 0; l < c->hindex1.size(); l++)
          {
            seg = &g4e->SegmentDetectors["EMCalSci"].at(c->hindex1.at(l));

            if(seg->GetPrimaryId() != p->tid)
            {
              continue;
            }

            found = true;
            
            double xstep = 0.5 * (seg->GetStart().X() + seg->GetStop().X());
            double ystep = 0.5 * (seg->GetStart().Y() + seg->GetStop().Y());
            double zstep = 0.5 * (seg->GetStart().Z() + seg->GetStop().Z());
            double tstep = 0.5 * (seg->GetStart().T() + seg->GetStop().T());
            
            std::cout << std::setw(3) << "idx" << " " 
              << std::setw(6) << "hid" << " " 
              << std::setw(4) << "np" << " " 
              << std::setw(6) << "pid" << " " 
              << std::setw(10) << "de" << " " 
              << std::setw(10) << "de2" << " " 
              << std::setw(11) << "l" << " " 
              << std::setw(11) << "x1" << " " 
              << std::setw(11) << "y1" << " " 
              << std::setw(11) << "z1" << " " 
              << std::setw(11) << "t1" << " " 
              << std::setw(11) << "x2"  << " " 
              << std::setw(11) << "y2"  << " " 
              << std::setw(11) << "z2"  << " " 
              << std::setw(11) << "t2"  << " " 
              << std::setw(11) << "xhit"  << " " 
              << std::setw(11) << "yhit"  << " " 
              << std::setw(11) << "zhit"  << " " 
              << std::setw(11) << "thit"  << " " 
              << std::endl; 
            
            std::cout << std::setw(3) << l << " " 
              << std::setw(6) << c->hindex1.at(l) << " " 
              << std::setw(4) << seg->Contrib.size() << " " 
              << std::setw(6) << seg->GetPrimaryId() << " " 
              << std::setw(10) << seg->GetEnergyDeposit() << " " 
              << std::setw(10) << seg->GetSecondaryDeposit() << " " 
              << std::setw(11) << seg->GetTrackLength() << " " 
              << std::setw(11) << seg->GetStart().X() << " " 
              << std::setw(11) << seg->GetStart().Y() << " " 
              << std::setw(11) << seg->GetStart().Z() << " " 
              << std::setw(11) << seg->GetStart().T() << " " 
              << std::setw(11) << seg->GetStop().X()  << " " 
              << std::setw(11) << seg->GetStop().Y()  << " " 
              << std::setw(11) << seg->GetStop().Z()  << " " 
              << std::setw(11) << seg->GetStop().T()  << " " 
              << std::setw(11) << xstep  << " " 
              << std::setw(11) << ystep  << " " 
              << std::setw(11) << zstep  << " " 
              << std::setw(11) << tstep  << " " 
              << std::endl; 
            
            unsigned int istep = 1;
            
            while(trj->Points[istep].GetPosition().T() < tstep && istep < trj->Points.size()){istep++;}

            if(istep == trj->Points.size())
            {
              std::cout << "Error: trajectory segment not found" << std::endl;
              std::cout << "--------------------------------------------------------------------------------------------------------------------------------------------------------" <<std::endl;
              err4++;
              continue;
            }

            double mind = min_dist(xstep, ystep, zstep, 
                trj->Points[istep-1].Position.X(), trj->Points[istep-1].Position.Y(), trj->Points[istep-1].Position.Z(), 
                trj->Points[istep].Position.X(), trj->Points[istep].Position.Y(), trj->Points[istep].Position.Z());
            
            std::cout << std::setw(3) << "ist" << " "
              << std::setw(6) << "pid" << " " 
              << std::setw(4) << "pdg" << " "  
              << std::setw(6) << "trid" << " " 
              << std::setw(10) << "ke" << " " 
              << std::setw(10) << "np" << " "
              << std::setw(11) << "tstep" << " " 
              << std::setw(11) << "x1" << " " 
              << std::setw(11) << "y1" << " " 
              << std::setw(11) << "z1" << " " 
              << std::setw(11) << "t1" << " " 
              << std::setw(11) << "x2" << " " 
              << std::setw(11) << "y2" << " " 
              << std::setw(11) << "z2" << " " 
              << std::setw(11) << "t2" << " "  
              << std::setw(11) << "mean dist." << " "
              << std::endl;
            
            std::cout << std::setw(3) << istep << " "
              << std::setw(6) << trj->GetParentId() << " " 
              << std::setw(4) << trj->GetPDGCode() << " "  
              << std::setw(6) << trj->GetTrackId() << " " 
              << std::setw(10) << trj->GetInitialMomentum().T() - trj->GetInitialMomentum().Mag() << " " 
              << std::setw(10) << trj->Points.size() << " "
              << std::setw(11) << tstep << " " 
              << std::setw(11) << trj->Points[istep-1].Position.X() << " " 
              << std::setw(11) << trj->Points[istep-1].Position.Y() << " " 
              << std::setw(11) << trj->Points[istep-1].Position.Z() << " " 
              << std::setw(11) << trj->Points[istep-1].Position.T() << " " 
              << std::setw(11) << trj->Points[istep].Position.X() << " " 
              << std::setw(11) << trj->Points[istep].Position.Y() << " " 
              << std::setw(11) << trj->Points[istep].Position.Z() << " " 
              << std::setw(11) << trj->Points[istep].Position.T() << " "  
              << std::setw(11) << mind << " "
              << std::endl;
          
            std::cout << "--------------------------------------------------------------------------------------------------------------------------------------------------------" <<std::endl;

            ok++;
            if(mind < 1.)
              good++;
          }

          if(found == false)
          {
            std::cout << "Error: particle hit not found" << std::endl;
            err3++;
          }

          found = false;

          hidx = -999;

          std::sort(c->hindex2.begin(), c->hindex2.end());
          c->hindex2.resize( std::distance(c->hindex2.begin(),std::unique (c->hindex2.begin(), c->hindex2.end())) );

          for(unsigned int l = 0; l < c->hindex2.size(); l++)
          {
            seg = &g4e->SegmentDetectors["EMCalSci"].at(c->hindex2.at(l));

            if(seg->GetPrimaryId() != p->tid)
            {
              continue;
            }

            found = true;
            
            double xstep = 0.5 * (seg->GetStart().X() + seg->GetStop().X());
            double ystep = 0.5 * (seg->GetStart().Y() + seg->GetStop().Y());
            double zstep = 0.5 * (seg->GetStart().Z() + seg->GetStop().Z());
            double tstep = 0.5 * (seg->GetStart().T() + seg->GetStop().T());

            std::cout << std::setw(3) << "idx" << " " 
              << std::setw(6) << "hid" << " " 
              << std::setw(4) << "np" << " " 
              << std::setw(6) << "pid" << " " 
              << std::setw(10) << "de" << " " 
              << std::setw(10) << "de2" << " " 
              << std::setw(11) << "l" << " " 
              << std::setw(11) << "x1" << " " 
              << std::setw(11) << "y1" << " " 
              << std::setw(11) << "z1" << " " 
              << std::setw(11) << "t1" << " " 
              << std::setw(11) << "x2"  << " " 
              << std::setw(11) << "y2"  << " " 
              << std::setw(11) << "z2"  << " " 
              << std::setw(11) << "t2"  << " " 
              << std::setw(11) << "xhit"  << " " 
              << std::setw(11) << "yhit"  << " " 
              << std::setw(11) << "zhit"  << " " 
              << std::setw(11) << "thit"  << " " 
              << std::endl; 
            
            std::cout << std::setw(3) << l << " " 
              << std::setw(6) << c->hindex2.at(l) << " " 
              << std::setw(4) << seg->Contrib.size() << " " 
              << std::setw(6) << seg->GetPrimaryId() << " " 
              << std::setw(10) << seg->GetEnergyDeposit() << " " 
              << std::setw(10) << seg->GetSecondaryDeposit() << " " 
              << std::setw(11) << seg->GetTrackLength() << " " 
              << std::setw(11) << seg->GetStart().X() << " " 
              << std::setw(11) << seg->GetStart().Y() << " " 
              << std::setw(11) << seg->GetStart().Z() << " " 
              << std::setw(11) << seg->GetStart().T() << " " 
              << std::setw(11) << seg->GetStop().X()  << " " 
              << std::setw(11) << seg->GetStop().Y()  << " " 
              << std::setw(11) << seg->GetStop().Z()  << " " 
              << std::setw(11) << seg->GetStop().T()  << " " 
              << std::setw(11) << xstep  << " " 
              << std::setw(11) << ystep  << " " 
              << std::setw(11) << zstep  << " " 
              << std::setw(11) << tstep  << " " 
              << std::endl; 
            
            unsigned int istep = 1;
            
            while(trj->Points[istep].GetPosition().T() < tstep && istep < trj->Points.size()){istep++;}

            if(istep == trj->Points.size())
            {
              std::cout << "Error: trajectory segment not found" << std::endl;
              std::cout << "--------------------------------------------------------------------------------------------------------------------------------------------------------" <<std::endl;
              err4++;
              continue;
            }

            double mind = min_dist(xstep, ystep, zstep, 
                trj->Points[istep-1].Position.X(), trj->Points[istep-1].Position.Y(), trj->Points[istep-1].Position.Z(), 
                trj->Points[istep].Position.X(), trj->Points[istep].Position.Y(), trj->Points[istep].Position.Z());
            
            std::cout << std::setw(3) << "ist" << " "
              << std::setw(6) << "pid" << " " 
              << std::setw(4) << "pdg" << " "  
              << std::setw(6) << "trid" << " " 
              << std::setw(10) << "ke" << " " 
              << std::setw(10) << "np" << " "
              << std::setw(11) << "tstep" << " " 
              << std::setw(11) << "x1" << " " 
              << std::setw(11) << "y1" << " " 
              << std::setw(11) << "z1" << " " 
              << std::setw(11) << "t1" << " " 
              << std::setw(11) << "x2" << " " 
              << std::setw(11) << "y2" << " " 
              << std::setw(11) << "z2" << " " 
              << std::setw(11) << "t2" << " "  
              << std::setw(11) << "mean dist." << " "
              << std::endl;
            
            std::cout << std::setw(3) << istep << " "
              << std::setw(6) << trj->GetParentId() << " " 
              << std::setw(4) << trj->GetPDGCode() << " "  
              << std::setw(6) << trj->GetTrackId() << " " 
              << std::setw(10) << trj->GetInitialMomentum().T() - trj->GetInitialMomentum().Mag() << " " 
              << std::setw(10) << trj->Points.size() << " "
              << std::setw(11) << tstep << " " 
              << std::setw(11) << trj->Points[istep-1].Position.X() << " " 
              << std::setw(11) << trj->Points[istep-1].Position.Y() << " " 
              << std::setw(11) << trj->Points[istep-1].Position.Z() << " " 
              << std::setw(11) << trj->Points[istep-1].Position.T() << " " 
              << std::setw(11) << trj->Points[istep].Position.X() << " " 
              << std::setw(11) << trj->Points[istep].Position.Y() << " " 
              << std::setw(11) << trj->Points[istep].Position.Z() << " " 
              << std::setw(11) << trj->Points[istep].Position.T() << " "  
              << std::setw(11) << mind << " "
              << std::endl;
          
            std::cout << "--------------------------------------------------------------------------------------------------------------------------------------------------------" <<std::endl;

            ok++;
            if(mind < 1.)
              good++;
          } 
        //}
        //}

          if(found == false)
          {
            std::cout << "Error: particle hit not found" << std::endl;
            err3++;
          }
      }
    }
  }

  std::cout << "Summary: " << std::endl;
  std::cout << "- " << std::setw(5) << err1 << " errors : [Error: ID(trj) != ID(part)]" << std::endl;
  std::cout << "- " << std::setw(5) << err2 << " errors : [Error: beta not correctly reconstructed]" << std::endl;
  std::cout << "- " << std::setw(5) << err3 << " errors : [Error: particle hit not found]" << std::endl;
  std::cout << "- " << std::setw(5) << err4 << " errors : [Error: trajectory segment not found]" << std::endl;
  std::cout << "- " << std::setw(5) << ok << " cells with trajectory" << std::endl;
  std::cout << "- " << std::setw(5) << good << " cells with good trajectory [min distance < 1 mm]" << std::endl;
}

