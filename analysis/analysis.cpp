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
const double m_to_mm = 1000.;
const double GeV_to_MeV = 1000.;
const double c = k * 1E3; // mm/ns
const double emk = 2.76455e-01;

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
  p.xtrue = -999.;
  p.ytrue = -999.;
  p.ztrue = -999.;
  p.ttrue = -999.;
  p.has_track = false;
  p.px_tr = -999.;
  p.py_tr = -999.;
  p.pz_tr = -999.;
  p.E_tr = -999.;
  p.x_tr = -999.;
  p.y_tr = -999.;
  p.z_tr = -999.;
  p.t_tr = -999.;
  p.has_cluster = false;
  p.px_cl = -999.;
  p.py_cl = -999.;
  p.pz_cl = -999.;
  p.E_cl = -999.;
  p.x_cl = -999.;
  p.y_cl = -999.;
  p.z_cl = -999.;
  p.t_cl = -999.;
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
      if(part != 0) p.mass = part->Mass() * GeV_to_MeV;
      
      p.pxtrue = ev->Trajectories.at(i).InitialMomentum.X();
      p.pytrue = ev->Trajectories.at(i).InitialMomentum.Y();
      p.pztrue = ev->Trajectories.at(i).InitialMomentum.Z();
      p.Etrue = ev->Trajectories.at(i).InitialMomentum.T();
      
      p.xtrue = ev->Trajectories.at(i).Points.at(0).Position.X();
      p.ytrue = ev->Trajectories.at(i).Points.at(0).Position.Y();
      p.ztrue = ev->Trajectories.at(i).Points.at(0).Position.Z();
      p.ttrue = ev->Trajectories.at(i).Points.at(0).Position.T();
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
  
  p.has_track = true;
  p.px_tr = mom_yz * TMath::Tan(ang_x);
  p.py_tr = mom_yz * TMath::Sin(ang_yz);
  p.pz_tr = mom_yz * TMath::Cos(ang_yz);
  p.E_tr = TMath::Sqrt(p.px_tr*p.px_tr + p.py_tr*p.py_tr + p.pz_tr*p.pz_tr + p.mass*p.mass);
  p.x_tr = tr.x0;
  p.y_tr = tr.y0;
  p.z_tr = tr.z0;
  p.t_tr = tr.t0;
}

bool isCellBefore(cell c1, cell c2) 
{
  if(c1.adc1 == 0 || c1.adc2 == 0)
    return false;
  else if(c2.adc1 == 0 || c2.adc2 == 0)
    return true;
  else
    return ((c1.tdc1 + c1.tdc2) < (c2.tdc1 + c2.tdc2));
}

double TfromTDC(double t1, double t2)
{
  return 0.5 * (t1 + t2 - ns_Digit::vlfb * ns_Digit::lCalBarrel);
}

double XfromTDC(double t1, double t2)
{
  return 0.5 * (t1 - t2)/ns_Digit::vlfb * m_to_mm;
}

void FillClusterInfo(TG4Event* ev, const cluster& cl, particle& p)
{ 
  p.has_cluster = true;
  
  switch(p.pdg)
  {
    case 2212: // proton
    case 2112: // neutron
    {
      // evaluate neutron velocity using earlier cell 
      int idx_min = std::distance(cl.cells.begin(), std::min_element(cl.cells.begin(), cl.cells.end(), isCellBefore));
      
      double tmin = TfromTDC(cl.cells.at(idx_min).tdc1, cl.cells.at(idx_min).tdc2);
      double xmin = XfromTDC(cl.cells.at(idx_min).tdc1, cl.cells.at(idx_min).tdc2);
      double ymin = cl.cells.at(idx_min).y;
      double zmin = cl.cells.at(idx_min).z;
      
      double dt = tmin - ev->Primaries.at(0).Position.T();
      double dx = xmin - ev->Primaries.at(0).Position.X();
      double dy = ymin - ev->Primaries.at(0).Position.Y();
      double dz = zmin - ev->Primaries.at(0).Position.Z();
      double dr = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
      double beta = dr / dt / c;
      double gamma = 1. / TMath::Sqrt(1 - beta*beta);
      
      p.E_cl = p.mass * gamma;
      p.px_cl = p.mass * gamma * beta * dx/dr;
      p.py_cl = p.mass * gamma * beta * dy/dr;
      p.pz_cl = p.mass * gamma * beta * dz/dr;
      
      //std::cout << idx_min << " " << tmin << " " << xmin << " " << ymin << " " << zmin << " " << dt << " " << dx << " " << dy << " " << dz << " " << dr << " " << beta << " " << gamma << " " << p.mass << std::endl;
      break;
    }
    case 211: // pi+
    case -211: // pi-
    {
      std::vector<double> x(5);
      std::vector<double> y(5);
      std::vector<double> z(5);
      std::vector<double> elayer(5);
      
      double etot = 0.0;
      
      for(unsigned int i = 0; i < cl.cells.size(); i++)
      {
        double e = cl.cells.at(i).adc1 + cl.cells.at(i).adc2;
        
        etot += e;
        
        elayer[cl.cells.at(i).lay] += e; 
        
        x[cl.cells.at(i).lay] += e * XfromTDC(cl.cells.at(i).tdc1, cl.cells.at(i).tdc2);
        y[cl.cells.at(i).lay] += e * cl.cells.at(i).y;
        z[cl.cells.at(i).lay] += e * cl.cells.at(i).z;
      }
      
      double sx = 0;
      double sy = 0;
      double sz = 0;
      double sxz = 0;
      double syz = 0;
      double sx2 = 0;
      double sy2 = 0;
      double sz2 = 0;
      int nlayer = 0;
      
      for(int i = 0; i < 5; i++)
      {
        if(elayer[i] != 0.)
        {
          x[i] /= elayer[i];
          y[i] /= elayer[i];
          z[i] /= elayer[i];
          
          sx += x[i];
          sy += y[i];
          sz += z[i];
          sxz += x[i]*z[i];
          syz += y[i]*z[i];
          sx2 += x[i]*x[i];
          sy2 += y[i]*y[i];
          sz2 += z[i]*z[i];
          
          nlayer++;
        }
      }
      
      double ax = (nlayer * sxz - sz*sx)/(nlayer*sz2 - sz*sz);
      double bx = (sx*sz2 - sz*sxz)/(nlayer*sz2 - sz*sz);
      double ay = (nlayer * syz - sz*sy)/(nlayer*sz2 - sz*sz);
      double by = (sy*sz2 - sz*syz)/(nlayer*sz2 - sz*sz);

      p.E_cl = etot;
      p.px_cl = etot * ax / TMath::Sqrt(1+ax*ax+ay*ay);
      p.py_cl = etot * ay / TMath::Sqrt(1+ax*ax+ay*ay);
      p.pz_cl = etot / TMath::Sqrt(1+ax*ax+ay*ay);
      
      break;
    }
    case 22: // gamma
    {
      std::vector<double> x(5);
      std::vector<double> y(5);
      std::vector<double> z(5);
      std::vector<double> elayer(5);
      
      double etot = 0.0;
      
      for(unsigned int i = 0; i < cl.cells.size(); i++)
      {
        double e = cl.cells.at(i).adc1 + cl.cells.at(i).adc2;
        
        etot += e;
        
        elayer[cl.cells.at(i).lay] += e; 
        
        x[cl.cells.at(i).lay] += e * XfromTDC(cl.cells.at(i).tdc1, cl.cells.at(i).tdc2);
        y[cl.cells.at(i).lay] += e * cl.cells.at(i).y;
        z[cl.cells.at(i).lay] += e * cl.cells.at(i).z;
      }
      
      double sx = 0;
      double sy = 0;
      double sz = 0;
      double sxz = 0;
      double syz = 0;
      double sx2 = 0;
      double sy2 = 0;
      double sz2 = 0;
      int nlayer = 0;
      double dir = -999.;
      bool ok = false;
      
      for(int i = 0; i < 5; i++)
      {
        if(elayer[i] != 0.)
        {        
          x[i] /= elayer[i];
          y[i] /= elayer[i];
          z[i] /= elayer[i];
          
          if(dir == -999.)
          {
            dir = z[i];
          }
          else if(!ok)
          {
            dir -= z[i];
            dir /= TMath::Abs(dir);
            ok = true;
          }
          
          sx += x[i];
          sy += y[i];
          sz += z[i];
          sxz += x[i]*z[i];
          syz += y[i]*z[i];
          sx2 += x[i]*x[i];
          sy2 += y[i]*y[i];
          sz2 += z[i]*z[i];
          
          nlayer++;
        }
      }
      
      double ax = (nlayer * sxz - sz*sx)/(nlayer*sz2 - sz*sz);
      double bx = (sx*sz2 - sz*sxz)/(nlayer*sz2 - sz*sz);
      double ay = (nlayer * syz - sz*sy)/(nlayer*sz2 - sz*sz);
      double by = (sy*sz2 - sz*syz)/(nlayer*sz2 - sz*sz);

      p.E_cl = etot * emk;
      p.px_cl = etot * emk * dir * ax / TMath::Sqrt(1+ax*ax+ay*ay);
      p.py_cl = etot * emk * dir * ay / TMath::Sqrt(1+ax*ax+ay*ay);
      p.pz_cl = etot * emk * dir / TMath::Sqrt(1+ax*ax+ay*ay);
      
      break;
    }
    case 111: // pi zero
    {
      std::vector<double> x(5);
      std::vector<double> y(5);
      std::vector<double> z(5);
      std::vector<double> elayer(5);
      
      double etot = 0.0;
      
      for(unsigned int i = 0; i < cl.cells.size(); i++)
      {
        double e = cl.cells.at(i).adc1 + cl.cells.at(i).adc2;
        
        etot += e;
        
        elayer[cl.cells.at(i).lay] += e; 
        
        x[cl.cells.at(i).lay] += e * XfromTDC(cl.cells.at(i).tdc1, cl.cells.at(i).tdc2);
        y[cl.cells.at(i).lay] += e * cl.cells.at(i).y;
        z[cl.cells.at(i).lay] += e * cl.cells.at(i).z;
      }
      
      double sx = 0;
      double sy = 0;
      double sz = 0;
      double sxz = 0;
      double syz = 0;
      double sx2 = 0;
      double sy2 = 0;
      double sz2 = 0;
      int nlayer = 0;
      double dir = -999.;
      bool ok = false;
      
      for(int i = 0; i < 5; i++)
      {
        if(elayer[i] != 0.)
        {        
          x[i] /= elayer[i];
          y[i] /= elayer[i];
          z[i] /= elayer[i];
          
          if(dir == -999.)
          {
            dir = z[i];
          }
          else if(!ok)
          {
            dir -= z[i];
            dir /= TMath::Abs(dir);
            ok = true;
          }
          
          sx += x[i];
          sy += y[i];
          sz += z[i];
          sxz += x[i]*z[i];
          syz += y[i]*z[i];
          sx2 += x[i]*x[i];
          sy2 += y[i]*y[i];
          sz2 += z[i]*z[i];
          
          nlayer++;
        }
      }
      
      double ax = (nlayer * sxz - sz*sx)/(nlayer*sz2 - sz*sz);
      double bx = (sx*sz2 - sz*sxz)/(nlayer*sz2 - sz*sz);
      double ay = (nlayer * syz - sz*sy)/(nlayer*sz2 - sz*sz);
      double by = (sy*sz2 - sz*syz)/(nlayer*sz2 - sz*sz);

      p.E_cl = etot * emk;
      p.px_cl = etot * emk * dir * ax / TMath::Sqrt(1+ax*ax+ay*ay);
      p.py_cl = etot * emk * dir * ay / TMath::Sqrt(1+ax*ax+ay*ay);
      p.pz_cl = etot * emk * dir / TMath::Sqrt(1+ax*ax+ay*ay);
      
      break;
    }
    default:
    {
      break;
    }
  }
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
  std::map<int, particle> map_part;
    
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
    map_part.clear();
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
        
        map_part[p.tid] = p;
      }
    } 
    
    for(unsigned int j = 0; j < vec_cl->size(); j++)
    {
      std::map<int, particle>::iterator it = map_part.find(vec_cl->at(j).tid);
    
      if(it != map_part.end())
      {
        FillClusterInfo(ev, vec_cl->at(j), it->second);
      }
      else
      {
        particle p;
        reset(p);
        
        p.tid = vec_cl->at(j).tid;
        
        FillParticleInfo(ev, p);
        FillClusterInfo(ev, vec_cl->at(j), p);
        
        map_part[p.tid] = p;
      }
    }
    
    for(std::map<int, particle>::iterator it = map_part.begin(); it != map_part.end(); ++it) 
    {
        vec_part.push_back(it->second);
    }
       
    tout.Fill();
  }
  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;
  
  fout.cd();
  tout.Write();
  fout.Close();
}