#include <TChain.h>
#include <TFile.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TDirectoryFile.h>

#include "TG4Event.h"
#include "TG4HitSegment.h"

#include "struct.h"
#include "utils.h"

#include <iostream>
#include <string>

const double k = 0.299792458;
const double B = 0.6;
const double m_to_mm = 1000.;
const double GeV_to_MeV = 1000.;
const double c = k * 1E3; // mm/ns
const double emk = 2.76455e-01;
const double hadk = 3.15234e-01;

TDatabasePDG db; 

void reset(particle& p)
{
  p.primary = false;
  p.pdg = 0;
  p.tid = 0;
  p.mass = 0.;
  p.charge = 0.;
  p.pxtrue = 0.;
  p.pytrue = 0.;
  p.pztrue = 0.;
  p.Etrue = 0.;
  p.xtrue = 0.;
  p.ytrue = 0.;
  p.ztrue = 0.;
  p.ttrue = 0.;
  p.has_track = false;
  p.charge_reco = 0.;
  p.pxreco = 0.;
  p.pyreco = 0.;
  p.pzreco = 0.;
  p.Ereco = 0.;
  p.xreco = 0.;
  p.yreco = 0.;
  p.zreco = 0.;
  p.treco = 0.;
  p.has_cluster = false;
  p.has_daughter = false;
}

bool isDigitBefore(digit d1, digit d2)
{
  return (d1.t < d2.t);
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

void FillParticleInfo(TG4Event* ev, std::map<int, particle>& map_part)
{
  for(unsigned int j = 0; j < ev->Trajectories.size(); j++)
  {
    particle p;
    reset(p);
    
    p.tid = ev->Trajectories.at(j).TrackId;
    p.pdg = ev->Trajectories.at(j).PDGCode;
    p.parent_tid = ev->Trajectories.at(j).ParentId;
    
    TParticlePDG* part = db.GetParticle(p.pdg);
    if(part != 0) 
    {
      p.mass = part->Mass() * GeV_to_MeV;
      p.charge = part->Charge()/3.;
    }
    
    p.primary = IsPrimary(ev, p.tid) == true ? 1 : 0;

    p.pxtrue = ev->Trajectories.at(j).InitialMomentum.X();
    p.pytrue = ev->Trajectories.at(j).InitialMomentum.Y();
    p.pztrue = ev->Trajectories.at(j).InitialMomentum.Z();
    p.Etrue = ev->Trajectories.at(j).InitialMomentum.T();
    
    p.xtrue = ev->Trajectories.at(j).Points.at(0).Position.X();
    p.ytrue = ev->Trajectories.at(j).Points.at(0).Position.Y();
    p.ztrue = ev->Trajectories.at(j).Points.at(0).Position.Z();
    p.ttrue = ev->Trajectories.at(j).Points.at(0).Position.T();
    
    map_part[p.tid] = p;
  }
}

/*
void FillTrackInfo(track& tr, particle& p)
{
  p.has_track = (tr.ret_ln == 0 && tr.ret_cr == 0) ? true : false;
  
  if(p.has_track)
  {
    std::sort(tr.digits.begin(), tr.digits.end(), isDigBefore);
    
    double yc = tr.yc;
    double zc = tr.zc;
    
    double l = 0.;
    
    for(unsigned int i = 0; i < tr.digits.size() - 1; i++)
    {
      double y0 = tr.digits.at(i).y;
      double z0 = tr.digits.at(i).z;
      double y1 = tr.digits.at(i+1).y;
      double z1 = tr.digits.at(i+1).z;
      
      double rz = z0 - zc;
      double ry = y0 - yc;
      double vz = z1 - z0;
      double vy = y1 - y0;
      
      l += ry * vz - rz * vy;
    }
    
    l /= TMath::Abs(l);
    
    double r0z = tr.z0 - tr.zc;
    double r0y = tr.y0 - tr.yc;
    
    double mom_yz = k * tr.r * B;
    double ang_yz = TMath::ATan2(r0z, -r0y); 
    double ang_x = 0.5 * TMath::Pi() - TMath::ATan(1./tr.b);
    
    p.charge_reco = -l;
    p.px_tr = mom_yz * TMath::Tan(ang_x);
    p.py_tr = p.charge_reco * mom_yz * TMath::Sin(ang_yz);
    p.pz_tr = p.charge_reco * mom_yz * TMath::Cos(ang_yz);
    p.E_tr = TMath::Sqrt(p.px_tr*p.px_tr + p.py_tr*p.py_tr + p.pz_tr*p.pz_tr + p.mass*p.mass);
    p.x_tr = tr.x0;
    p.y_tr = tr.y0;
    p.z_tr = tr.z0;
    p.t_tr = tr.t0;
    
    p.tr = tr;
  }
}
*/

bool isCellBefore(cell c1, cell c2) 
{
  if(c1.adc1 == 0 || c1.adc2 == 0)
    return false;
  else if(c2.adc1 == 0 || c2.adc2 == 0)
    return true;
  else
    return ((c1.tdc1 + c1.tdc2) < (c2.tdc1 + c2.tdc2));
}

double TfromTDC(double t1, double t2, double L)
{
  return 0.5 * (t1 + t2 - ns_Digit::vlfb * L / m_to_mm );
}

double XfromTDC(double t1, double t2, double x0)
{
  return 0.5 * (t1 - t2)/ns_Digit::vlfb * m_to_mm + x0;
}

void CellXYZTE(cell c, double& x, double& y, double& z, double& t, double& e)
{      
  if(c.id < 25000) //Barrel
  {
    x = 0.5 * (c.tdc1 - c.tdc2)/ns_Digit::vlfb * m_to_mm + c.x;
    y = c.y;
  }
  else
  {
    x = c.x;
    y = -0.5 * (c.tdc1 - c.tdc2)/ns_Digit::vlfb * m_to_mm + c.y;
  }
  z = c.z;
  t = 0.5 * (c.tdc1 + c.tdc2 - ns_Digit::vlfb * c.l / m_to_mm );
  e = c.adc1 + c.adc2;
}

void RecoFromTrack(particle& p)
{
  if(p.tr.ret_ln == 0 && p.tr.ret_cr == 0)
  {
    std::sort(p.tr.digits.begin(), p.tr.digits.end(), isDigBefore);
    
    double yc = p.tr.yc;
    double zc = p.tr.zc;
    
    double l = 0.;
    
    for(unsigned int i = 0; i < p.tr.digits.size() - 1; i++)
    {
      double y0 = p.tr.digits.at(i).y;
      double z0 = p.tr.digits.at(i).z;
      double y1 = p.tr.digits.at(i+1).y;
      double z1 = p.tr.digits.at(i+1).z;
      
      double rz = z0 - zc;
      double ry = y0 - yc;
      double vz = z1 - z0;
      double vy = y1 - y0;
      
      l += ry * vz - rz * vy;
    }
    
    l /= TMath::Abs(l);
    
    double r0z = p.tr.z0 - p.tr.zc;
    double r0y = p.tr.y0 - p.tr.yc;
    
    double mom_yz = k * p.tr.r * B;
    double ang_yz = TMath::ATan2(r0z, -r0y); 
    double ang_x = 0.5 * TMath::Pi() - TMath::ATan(1./p.tr.b);
    
    p.charge_reco = -l;
    p.pxreco = mom_yz * TMath::Tan(ang_x);
    p.pyreco = p.charge_reco * mom_yz * TMath::Sin(ang_yz);
    p.pzreco = p.charge_reco * mom_yz * TMath::Cos(ang_yz);
    p.Ereco = TMath::Sqrt(p.pxreco*p.pxreco + p.pyreco*p.pyreco + p.pzreco*p.pzreco + p.mass*p.mass);
    p.xreco = p.tr.x0;
    p.yreco = p.tr.y0;
    p.zreco = p.tr.z0;
    p.treco = p.tr.t0;
  }
  else
  {
    p.charge_reco = 0.;
    p.pxreco = 0.;
    p.pyreco = 0.;
    p.pzreco = 0.;
    p.Ereco = 0.;
    p.xreco = 0.;
    p.yreco = 0.;
    p.zreco = 0.;
    p.treco = 0.;
  }
}

void RecoFromBeta(particle& p, double x0, double y0, double z0, double t0)
{      
  // evaluate neutron velocity using earlier cell 
  std::vector<double> cell_t(p.cl.cells.size());
  std::vector<double> cell_x(p.cl.cells.size());
  std::vector<double> cell_y(p.cl.cells.size());
  std::vector<double> cell_z(p.cl.cells.size());
  
  double e;
  
  for(unsigned int i = 0; i < p.cl.cells.size(); i++)
  {
    CellXYZTE(p.cl.cells.at(i), cell_x.at(i), cell_y.at(i), cell_z.at(i), cell_t.at(i), e);
  }
  
  int idx_min = std::distance(cell_t.begin(), std::min_element(cell_t.begin(), cell_t.end()));
  
  double dt = cell_t.at(idx_min) - t0;
  double dx = cell_x.at(idx_min) - x0;
  double dy = cell_y.at(idx_min) - y0;
  double dz = cell_z.at(idx_min) - z0;
  double dr = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
  double beta = dr / dt / c;
  double gamma = 1. / TMath::Sqrt(1 - beta*beta);
  
  if(beta < 1.)
  {
    p.Ereco = p.mass * gamma;
    p.pxreco = p.mass * gamma * beta * dx/dr;
    p.pyreco = p.mass * gamma * beta * dy/dr;
    p.pzreco = p.mass * gamma * beta * dz/dr;
  }
  else
  {
    p.Ereco = -0.;
    p.pxreco = -0.;
    p.pyreco = -0.;
    p.pzreco = -0.;
  }
}

void RecoFromEMShower(particle& p)
{
  p.Ereco = p.cl.e * emk;
  
  double mom = TMath::Sqrt(p.Ereco*p.Ereco - p.mass*p.mass);
  
  p.pxreco = mom * emk * p.cl.sx;
  p.pyreco = mom * emk * p.cl.sy;
  p.pzreco = mom * emk * p.cl.sz;
}

void RecoFromHadShower(particle& p)
{
  p.Ereco = p.cl.e * hadk;
  
  double mom = TMath::Sqrt(p.Ereco*p.Ereco - p.mass*p.mass);
  
  p.pxreco = hadk * emk * p.cl.sx;
  p.pyreco = hadk * emk * p.cl.sy;
  p.pzreco = hadk * emk * p.cl.sz;
}

void RecoFromDaugthers(particle& p)
{
  for(unsigned int i = 0; i < p.daughters.size(); i++)
  {      
    p.Ereco += p.daughters.at(i).Ereco;
    p.pxreco += p.daughters.at(i).pxreco;
    p.pyreco += p.daughters.at(i).pyreco;
    p.pzreco += p.daughters.at(i).pzreco;
  }
}

void RecoGamma(particle& p)
{
  bool ok = false;
  if(p.has_daughter == 1)
  {
    for(unsigned int i = 0; i < p.daughters.size(); i++)
    {
      if(p.daughters.at(i).has_track == 1 && p.daughters.at(i).tr.ret_ln == 0 && p.daughters.at(i).tr.ret_cr == 0)
      {
        ok = true;
      }
    }
  }  
  
  if(ok)
    RecoFromDaugthers(p);
  else if(p.has_cluster == 1)
    RecoFromEMShower(p);
}

void RecoPi0(particle& p)
{
  for(unsigned int i = 0; i < p.daughters.size(); i++)
  {
    if(p.daughters.at(i).pdg == 22 && p.daughters.at(i).has_daughter == 1)
    {
      RecoGamma(p.daughters.at(i));
    }
  }
  
  RecoFromDaugthers(p);
}

void FindGammaConversion(event& ev, particle& p)
{     
  for(unsigned int j = 0; j < ev.particles.size(); j++)
  {
    if(ev.particles.at(j).parent_tid == p.tid)
    {
      p.has_daughter = 1;
      p.daughters.push_back(ev.particles.at(j));
    }
  }
}

void FindPriGammaConversion(event& ev)
{
  for(unsigned int i = 0; i < ev.particles.size(); i++)
  {
    if(ev.particles.at(i).primary == 1 && ev.particles.at(i).pdg == 22)
    {
      FindGammaConversion(ev, ev.particles.at(i));
    }
  }
}

void FindPi0Decay(event& ev, particle& p)
{
  for(unsigned int j = 0; j < ev.particles.size(); j++)
  {
    if(ev.particles.at(j).parent_tid == p.tid)
    {
      if(ev.particles.at(j).pdg == 22)
      {
        FindGammaConversion(ev, ev.particles.at(j));
        
        p.has_daughter = 1;
        p.daughters.push_back(ev.particles.at(j));
        
      }
    }
  }
}

void FindPriPi0Decay(event& ev)
{
  for(unsigned int i = 0; i < ev.particles.size(); i++)
  {
    if(ev.particles.at(i).primary == 1 && ev.particles.at(i).pdg == 111)
    { 
      FindPi0Decay(ev, ev.particles.at(i));
    }
  }
}

void ProcessParticle(event& evt, int index)
{
  particle& p = evt.particles.at(index);

  switch(p.pdg)
  {
    case -2212: // antiproton
    case  2212: // proton
    {
      if(p.has_track == 1 && p.tr.ret_ln == 0 && p.tr.ret_cr == 0)
        RecoFromTrack(p);
      else if(p.has_cluster == 1)
        RecoFromBeta(p,evt.x,evt.y,evt.z,evt.t);
      break;
    }
    case  -211: // antipion
    case   211: // pion
    {
      if(p.has_track == 1 && p.tr.ret_ln == 0 && p.tr.ret_cr == 0)
        RecoFromTrack(p);
      else if(p.has_cluster == 1)
        RecoFromHadShower(p);
      break;
    }
    case    11: // electron
    case   -11: // positron
    {
      if(p.has_track == 1 && p.tr.ret_ln == 0 && p.tr.ret_cr == 0)
        RecoFromTrack(p);
      else if(p.has_cluster == 1)
        RecoFromEMShower(p);
      break;
    }
    case  2112: // neutron
    case -2112: // antineutron
    {
      if(p.has_cluster == 1)
        RecoFromBeta(p,evt.x,evt.y,evt.z,evt.t);
      break;
    }
    case 22: // gamma
    {
      if(p.primary==1)
      {
        FindGammaConversion(evt, p);
      }
      RecoGamma(p);
      break;
    }
    case 111: // pi zero
    {
      if(p.primary==1)
      {
        FindPi0Decay(evt, p);
      }
      RecoPi0(p);
      break;
    }
    default: // other (assuming hadron)
    {
      if(p.has_track == 1 && p.tr.ret_ln == 0 && p.tr.ret_cr == 0)
        RecoFromTrack(p);
      else if(p.has_cluster == 1)
        RecoFromHadShower(p);
      break;
    }
  }
}

void ProcessParticles(event& evt)
{
  for(unsigned int i = 0; i < evt.particles.size(); i++)
  { 
    ProcessParticle(evt, i);
  }
}

/*
void FillClusterInfo(TG4Event* ev, const cluster& cl, particle& p)
{ 
  p.has_cluster = true;
  
  switch(p.pdg)
  {
    case 2212: // proton
    case 2112: // neutron
    case -2212: // antiproton
    case -2112: // antineutron
    {
      // evaluate neutron velocity using earlier cell 
      std::vector<double> cell_time(cl.cells.size());
      
      for(unsigned int i = 0; i < cl.cells.size(); i++)
      {
        cell_time[i] = TfromTDC(cl.cells.at(i).tdc1, cl.cells.at(i).tdc2, cl.cells.at(i).l);
        
        //std::cout << cell_time[i] << " " << cl.cells.at(i).tdc1 << " " << cl.cells.at(i).tdc2 << " " << cl.cells.at(i).l << std::endl;
      }
      
      int idx_min = std::distance(cell_time.begin(), std::min_element(cell_time.begin(), cell_time.end()));
      
      double xmin, ymin;
      
      if(cl.cells.at(idx_min).id < 25000)
      {
        xmin = XfromTDC(cl.cells.at(idx_min).tdc1, cl.cells.at(idx_min).tdc2, cl.cells.at(idx_min).x);
        ymin = cl.cells.at(idx_min).y;
      }
      else
      {
        ymin = XfromTDC(cl.cells.at(idx_min).tdc1, cl.cells.at(idx_min).tdc2, cl.cells.at(idx_min).y);
        xmin = cl.cells.at(idx_min).x;
      }
      
      double tmin = cell_time.at(idx_min);
      double zmin = cl.cells.at(idx_min).z;
      
      double dt = tmin - ev->Primaries.at(0).Position.T();
      double dx = xmin - ev->Primaries.at(0).Position.X();
      double dy = ymin - ev->Primaries.at(0).Position.Y();
      double dz = zmin - ev->Primaries.at(0).Position.Z();
      double dr = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
      double beta = dr / dt / c;
      double gamma = 1. / TMath::Sqrt(1 - beta*beta);
      
      //std::cout << "-> " << idx_min << " " << cell_time[idx_min] << " " 
      //  << cl.cells.at(idx_min).tdc1 << " " << cl.cells.at(idx_min).tdc2 << " " << cl.cells.at(idx_min).l << " "  
      //  << xmin << " " << ymin << " " << zmin << " " 
      //  << ev->Primaries.at(0).Position.X() << " " << ev->Primaries.at(0).Position.Y() << " " << ev->Primaries.at(0).Position.Z() << " "
      //  << ev->Primaries.at(0).Position.T() << " " << dt << " " << dr << " " << beta << std::endl;
      
      if(beta < 1.)
      {
        p.E_cl = p.mass * gamma;
        p.px_cl = p.mass * gamma * beta * dx/dr;
        p.py_cl = p.mass * gamma * beta * dy/dr;
        p.pz_cl = p.mass * gamma * beta * dz/dr;
      }
      else
      {
        p.E_cl = -999.;
        p.px_cl = -999.;
        p.py_cl = -999.;
        p.pz_cl = -999.;
      }
      
      p.cl = cl; 
      
      //std::cout << idx_min << " " << tmin << " " << xmin << " " << ymin << " " << zmin << " " << dt << " " << dx << " " << dy << " " << dz << " " << dr << " " << beta << " " << gamma << " " << p.mass << std::endl;
      break;
    }
    case 22: // gamma
    case 11: // electron
    case -11: // positron
    case 111: // pi zero
    {
      std::vector<double> x(5);
      std::vector<double> y(5);
      std::vector<double> z(5);
      std::vector<double> elayer(5);
      
      double x2_cl_mean = 0.;
      double y2_cl_mean = 0.;
      double z2_cl_mean = 0.;
      
      double x_cl_mean = 0.;
      double y_cl_mean = 0.;
      double z_cl_mean = 0.;
      
      double etot = 0.0;
        
      double xv,yv,zv;
      
      for(unsigned int i = 0; i < cl.cells.size(); i++)
      {
        double e = cl.cells.at(i).adc1 + cl.cells.at(i).adc2;
        
        etot += e;
        
        elayer[cl.cells.at(i).lay] += e; 
      
        if(cl.cells.at(i).id < 25000)
        {
          xv = XfromTDC(cl.cells.at(i).tdc1, cl.cells.at(i).tdc2, cl.cells.at(i).x);
          yv = cl.cells.at(i).y;
          zv = cl.cells.at(i).z;
        }
        else
        {
          xv = cl.cells.at(i).x;
          yv = XfromTDC(cl.cells.at(i).tdc1, cl.cells.at(i).tdc2, cl.cells.at(i).y);
          zv = cl.cells.at(i).z;
        }
          
        x[cl.cells.at(i).lay] += e * xv;
        y[cl.cells.at(i).lay] += e * yv;
        z[cl.cells.at(i).lay] += e * zv;
        
        x2_cl_mean += xv*xv;
        x_cl_mean += xv;
        y2_cl_mean += yv*yv;
        y_cl_mean += yv;
        z2_cl_mean += zv*zv;
        z_cl_mean += zv;
      }
      
      x2_cl_mean /= cl.cells.size();
      x_cl_mean /= cl.cells.size();
      y2_cl_mean /= cl.cells.size();
      y_cl_mean /= cl.cells.size();
      z2_cl_mean /= cl.cells.size();
      z_cl_mean /= cl.cells.size();
      
      p.x_cl_var = x2_cl_mean - x_cl_mean*x_cl_mean;
      p.y_cl_var = y2_cl_mean - y_cl_mean*y_cl_mean;
      p.z_cl_var = z2_cl_mean - z_cl_mean*z_cl_mean;
      
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
      
      p.cl = cl;
      
      break;
    }
    default: // other (assuming hadron)
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
      
        if(cl.cells.at(i).id < 25000)
        {
          x[cl.cells.at(i).lay] += e * XfromTDC(cl.cells.at(i).tdc1, cl.cells.at(i).tdc2, cl.cells.at(i).x);
          y[cl.cells.at(i).lay] += e * cl.cells.at(i).y;
          z[cl.cells.at(i).lay] += e * cl.cells.at(i).z;
        }
        else
        {
          x[cl.cells.at(i).lay] += e * cl.cells.at(i).x;
          y[cl.cells.at(i).lay] += e * XfromTDC(cl.cells.at(i).tdc1, cl.cells.at(i).tdc2, cl.cells.at(i).y);
          z[cl.cells.at(i).lay] += e * cl.cells.at(i).z;
        }
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

      p.E_cl = hadk * etot;
      p.px_cl = hadk * etot * ax / TMath::Sqrt(1+ax*ax+ay*ay);
      p.py_cl = hadk * etot * ay / TMath::Sqrt(1+ax*ax+ay*ay);
      p.pz_cl = hadk * etot / TMath::Sqrt(1+ax*ax+ay*ay);
      
      p.cl = cl;
      
      break;
    }
  }
}
*/

bool isAfter(particle p1, particle p2)
{
  return p1.tid > p2.tid;
}

void EvalNuEnergy(event& ev)
{
  ev.Enureco = 0.;
  ev.pxnureco = 0.;
  ev.pynureco = 0.;
  ev.pznureco = 0.;
  
  for(unsigned int i = 0; i < ev.particles.size(); i++)
  {
    if(ev.particles.at(i).primary == 1)
    {
      if((ev.particles.at(i).pdg == 2212 || ev.particles.at(i).pdg == 2112) && ev.particles.at(i).Ereco > 0) //proton
      {
        ev.Enureco  += ev.particles.at(i).Ereco - ev.particles.at(i).mass;
        ev.pxnureco += ev.particles.at(i).pxreco;
        ev.pynureco += ev.particles.at(i).pyreco;
        ev.pznureco += ev.particles.at(i).pzreco;
      }
      else
      {
        ev.Enureco  += ev.particles.at(i).Ereco;
        ev.pxnureco += ev.particles.at(i).pxreco;
        ev.pynureco += ev.particles.at(i).pyreco;
        ev.pznureco += ev.particles.at(i).pzreco;
      }
    }
  }
}

void Analyze(const char* fIn)
{
  TFile f(fIn,"UPDATE");
  TTree* tReco = (TTree*) f.Get("tReco");
  TTree* tTrueMC = (TTree*) f.Get("EDepSimEvents");
  TTree* gRooTracker = (TTree*) f.Get("gRooTracker");
  TGeoManager* geo = (TGeoManager*) f.Get("EDepSimGeometry");
  
  tReco->AddFriend(tTrueMC);
  tReco->AddFriend(gRooTracker);
  
  TTree* t = tReco;
  
  std::vector<track>* vec_tr = new std::vector<track>;
  std::vector<cluster>* vec_cl = new std::vector<cluster>;

  const int kMaxStdHepN = 200;

  double part_mom[kMaxStdHepN][4];
  int part_pdg[kMaxStdHepN];
  
  TG4Event* ev = new TG4Event;
  t->SetBranchAddress("Event",&ev);
  t->SetBranchAddress("track",&vec_tr);
  t->SetBranchAddress("cluster",&vec_cl);
  t->SetBranchAddress("StdHepP4",part_mom);
  t->SetBranchAddress("StdHepPdg",part_pdg);
  
  std::map<int, particle> map_part;
  
  event evt;

  TTree tout("tEvent","tEvent");
  tout.Branch("event","event",&evt);
    
  const int nev = t->GetEntries();
  
  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;
  
  for(int i = 0; i < nev; i++)
  {
    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;
    
    t->GetEntry(i);
    map_part.clear();
    evt.particles.clear();

    evt.x = ev->Primaries.at(0).Position.X();
    evt.y = ev->Primaries.at(0).Position.Y();
    evt.z = ev->Primaries.at(0).Position.Z();
    evt.t = ev->Primaries.at(0).Position.T();
    
    evt.pxnu = part_mom[0][0]*GeV_to_MeV;
    evt.pynu = part_mom[0][1]*GeV_to_MeV;
    evt.pznu = part_mom[0][2]*GeV_to_MeV;
    evt.Enu = part_mom[0][3]*GeV_to_MeV;

    //std::string volname = geo->FindNode(evt.x, evt.y, evt.z)->GetName();
    //strcpy(evt.vol, geo->FindNode(evt.x, evt.y, evt.z)->GetName());
    //strcpy(evt.intType, ev->Primaries.at(0).Reaction.c_str());
    //evt.vol = volname;
    //evt.pdgnu = 1;//part_pdg[0];
    //evt.isCC = false;
    //evt.intType = ev->Primaries.at(0).Reaction;
    //evt.isCC = (strstr(evt.intType,"CC") != 0);

    //std::cout << evt.vol << " " << evt.intType << std::endl;
    
    FillParticleInfo(ev, map_part);
    
    for(unsigned int j = 0; j < vec_tr->size(); j++)
    {
      std::map<int, particle>::iterator it = map_part.find(vec_tr->at(j).tid);
      //FillTrackInfo(vec_tr->at(j), it->second);
      it->second.has_track = true;
      it->second.tr = vec_tr->at(j);
    } 
    
    for(unsigned int j = 0; j < vec_cl->size(); j++)
    {
      std::map<int, particle>::iterator it = map_part.find(vec_cl->at(j).tid);
      //FillClusterInfo(ev, vec_cl->at(j), it->second);
      it->second.has_cluster = true;
      it->second.cl = vec_cl->at(j);
    }
    
    for(std::map<int, particle>::iterator it = map_part.begin(); it != map_part.end(); ++it) 
    {
        evt.particles.push_back(it->second);
    }
    
    std::sort(evt.particles.begin(), evt.particles.end(), isAfter);
    
    //FindPriGammaConversion(evt);
    //FindPriPi0Decay(evt);
    
    ProcessParticles(evt);
    
    EvalNuEnergy(evt);
    
    tout.Fill();
  }
  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;
  
  f.cd();
  tout.Write("",TObject::kOverwrite);
  f.Close();

  vec_tr->clear();
  vec_cl->clear();
  delete vec_tr;
  delete vec_cl;
}

void help_ana()
{
  std::cout << "Analyze <input file>" << std::endl;
} 

int main(int argc, char* argv[])
{
  if(argc != 2)
    help_ana();
  else
    Analyze(argv[1]);
}


