#include <TChain.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TDirectoryFile.h>

#include "TG4Event.h"
#include "TG4HitSegment.h"

#include <iostream>
#include <algorithm>

#include "struct.h"
#include "utils.h"

const double m_to_mm = 1000.;

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

bool ishitok(TG4Event* ev, int trackid, TG4HitSegment hit,
             double postol = 5., double angtol = 0.3)
{
    double x = 0.5*(hit.Start.X()+hit.Stop.X());
    double y = 0.5*(hit.Start.Y()+hit.Stop.Y());
    double z = 0.5*(hit.Start.Z()+hit.Stop.Z());
    
    for(unsigned int jj = 0; jj < ev->Trajectories.size(); jj++)
    {
        //if(ev->Trajectories[jj].TrackId == hit.PrimaryId)
        if(ev->Trajectories[jj].TrackId == trackid)
        {
            for(unsigned int kk = 0; kk < ev->Trajectories[jj].Points.size()-1; kk++)
            {
                double dpos = mindist(ev->Trajectories[jj].Points[kk].Position.X(),
                                       ev->Trajectories[jj].Points[kk].Position.Y(),
                                       ev->Trajectories[jj].Points[kk].Position.Z(),
                                       ev->Trajectories[jj].Points[kk+1].Position.X(),
                                       ev->Trajectories[jj].Points[kk+1].Position.Y(),
                                       ev->Trajectories[jj].Points[kk+1].Position.Z(),
                                       x,y,z);
                double dang = angle(ev->Trajectories[jj].Points[kk+1].Position.X()-ev->Trajectories[jj].Points[kk].Position.X(),
                                     ev->Trajectories[jj].Points[kk+1].Position.Y()-ev->Trajectories[jj].Points[kk].Position.Y(),
                                     ev->Trajectories[jj].Points[kk+1].Position.Z()-ev->Trajectories[jj].Points[kk].Position.Z(),
                                     hit.Stop.X()-hit.Start.X(),
                                     hit.Stop.Y()-hit.Start.Y(),
                                     hit.Stop.Z()-hit.Start.Z());
                /*
                std::cout << ev->Trajectories[jj].Points[kk].Position.X() << " " <<
                             ev->Trajectories[jj].Points[kk].Position.Y() << " " <<
                             ev->Trajectories[jj].Points[kk].Position.Z() << " " <<
                             ev->Trajectories[jj].Points[kk+1].Position.X() << " " <<
                             ev->Trajectories[jj].Points[kk+1].Position.Y() << " " <<
                             ev->Trajectories[jj].Points[kk+1].Position.Z() << " " <<
                             x << " " << y << " " << z << " " << dpos.back() << " " <<
                             ev->Trajectories[jj].TrackId << " " << ev->Trajectories[jj].PDGCode << std::endl;
                */
                
                if(dpos < postol && dang < angtol)
                  return true;
            }
        }
    }
    return false;
}

int fitCircle(int n, const std::vector<double>& x, const std::vector<double>& y, double &xc, double &yc, double &r, double &errr, double &chi2)
{
    xc = -999;
    yc = -999;
    r = -999;
    errr = -999;
    chi2 = -999;
    
    if(x.size() != y.size())
      return 1;
    
    if(n < 3)
      return 2;
    
    double sumx = 0, sumy = 0;                            // linear    terms
    double sumx2 = 0, sumy2 = 0, sumxy = 0;               // quadratic terms
    double sumxy2 = 0, sumx2y = 0, sumx3 = 0, sumy3 = 0;  // cubic     terms
    
    for (int i = 0; i < n; i++)
    {
            double xp = x.at(i);
            double yp = y.at(i);
            sumx   += xp;       sumy   += yp;
            sumx2  += xp*xp;    sumy2  += yp*yp;    sumxy += xp*yp;
            sumxy2 += xp*yp*yp; sumx2y += xp*xp*yp; sumx3 += xp*xp*xp; sumy3 += yp*yp*yp;
    }

    double a = n*sumx2 - sumx*sumx;
    double b = n*sumxy - sumx*sumy;
    double c = n*sumy2 - sumy*sumy;
    double d = 0.5*(n*sumxy2 - sumx*sumy2 + n*sumx3 - sumx*sumx2);
    double e = 0.5*(n*sumx2y - sumy*sumx2 + n*sumy3 - sumy*sumy2);
    
    if(a*c - b*b == 0.)
      return 3;

    xc = (d*c - b*e) / (a*c - b*b);
    yc = (a*e - b*d) / (a*c - b*b);

    double rMean = 0;
    double rrms = 0;

    for (int i = 0; i < n; i++)
    {
        double xp = x.at(i);
        double yp = y.at(i);
        double r2 = (xp - xc)*(xp - xc) + (yp - yc)*(yp - yc);

        rMean += sqrt(r2);
        rrms += r2;
    }

    rMean /= n;
    rrms /= n;
    r = rMean;
      
    errr = sqrt(rrms - rMean*rMean);
    
    chi2 = 0.0;
    
    for(int i = 0; i < n; i++)
    {
        chi2 += TMath::Abs((y.at(i) - yc)*(y.at(i) - yc) + (x.at(i) - xc) * (x.at(i) - xc) - r*r);
    }
    
    chi2 /= n;
    
    /*
    std::cout << "==== ZY =====" << std::endl;
    
    for(int i = 0; i < n; i++)
    {
      std::cout << x.at(i) << " " << y.at(i) << std::endl;
    }
    
    std::cout << "---> " << xc << " " << yc << " " << r << " " << chi2 << std::endl;
    */
    
    return 0;
}

int fitLinear(int n, const std::vector<double>& x, const std::vector<double>& y, double &a, double &b, double cov[][2], double &chi2)
{
    a = -999;
    b = -999;
    chi2 = -999;
    
    if(x.size() != y.size())
      return 1;
    
    if(n < 2)
      return 2;
    
    cov[0][0] = -999;;
    cov[0][1] = -999;
    cov[1][0] = -999;
    cov[1][1] = -999;
    
    double S1  = 0;
    double SX  = 0;
    double SY  = 0;
    double SYX = 0;
    double SXX = 0;
    
    for(int i = 0; i < n; i++)
    {
        //std::cout << "XY: " <<  x.at(i) << " " <<  y.at(i) << std::endl;
    
        S1 += 1;
        SX += x.at(i);
        SY += y.at(i);
        SYX += x.at(i)*y.at(i);
        SXX += x.at(i)*x.at(i);
    }
    
    double D = S1*SXX - SX*SX;
    
    //std::cout << "D: " << D << " " << S1 << " " << SX << " " << SY << " " << SXX << " " << SYX << std::endl;
    
    a = (SY*SXX - SX*SYX)/D;  // i.p. at x = 0
    b = (S1*SYX - SX*SY)/D;   // tg(theta)
    
    //std::cout << a << " " << b << std::endl;
    
    cov[0][0] = SXX/D;
    cov[0][1] = -SX/D;
    cov[1][0] = -SX/D;
    cov[1][1] = S1/D;
    
    chi2 = 0.0;
    
    for(int i = 0; i < n; i++)
    {
        chi2 += (y.at(i) - a - b * x.at(i)) * (y.at(i) - a - b * x.at(i));
    }
    
    chi2 /= n;
    
    /*
    std::cout << "==== ZX =====" << std::endl;
    
    for(int i = 0; i < n; i++)
    {
      std::cout << x.at(i) << " " << y.at(i) << std::endl;
    }
    
    std::cout << "---> " << a << " " << b << " " << chi2 << std::endl;
    */
    
    return 0;
}

void TrackFind(TG4Event* ev, std::vector<digit>* vec_digi, std::vector<track>& vec_tr)
{
  vec_tr.clear();    
          
  //for(unsigned int j = 0; j < ev->Primaries[0].Particles.size(); j++)
  for(unsigned int j = 0; j < ev->Trajectories.size(); j++)
  {
    track tr;
          
    tr.tid = ev->Trajectories.at(j).TrackId;
    
    for(unsigned int k = 0; k < vec_digi->size(); k++)
    {
      for(unsigned int m = 0; m < vec_digi->at(k).hindex.size(); m++)
      {
        const TG4HitSegment& hseg = ev->SegmentDetectors["StrawTracker"].at(vec_digi->at(k).hindex.at(m)); 
               
        //if(hseg.PrimaryId == tr.tid)
        //{
          //if(ishitok(ev, hseg->PrimaryId, hseg))
          if(ishitok(ev, tr.tid, hseg))
          {
            tr.digits.push_back(vec_digi->at(k));
            break;
          }
        //}
      }
    }
    
    std::sort(tr.digits.begin(), tr.digits.end(), isDigBefore);
    
    vec_tr.push_back(tr);
  }
}

void TrackFit(std::vector<track>& vec_tr)
{
    
  std::vector<double> y_h;
  std::vector<double> z_h;
  std::vector<double> x_v;
  std::vector<double> y_v;
  std::vector<double> z_v;
  std::vector<double> rho;
  
  for(unsigned int j = 0; j < vec_tr.size(); j++)
  {
    y_h.clear();
    z_h.clear();
    x_v.clear();
    y_v.clear();
    z_v.clear();
    rho.clear();
    
    for(unsigned int k = 0; k < vec_tr[j].digits.size(); k++)
    {    
      if(vec_tr[j].digits.at(k).hor)
      {
        y_h.push_back(vec_tr[j].digits.at(k).y);
        z_h.push_back(vec_tr[j].digits.at(k).z);
      }
      else
      {
        x_v.push_back(vec_tr[j].digits.at(k).x);
        z_v.push_back(vec_tr[j].digits.at(k).z);
      }
    }
    
    double chi2_cir, chi2_lin;
    double errr;
    
    double cov[2][2];
    
    int n_h = z_h.size();
    int n_v = z_v.size();
    
    //std::cout << n_h << " " << n_v << std::endl;
    
    if(n_v <= 2 || n_h <= 2)
      continue;
    
    n_v = 2;
    
    while((n_v < int(z_v.size())) && (z_v[n_v-1]-z_v[n_v-2])*(z_v[1]-z_v[0]) > 0.)
      n_v++;
    
    int ret1 = fitCircle(n_h, z_h, y_h, vec_tr[j].zc, vec_tr[j].yc, vec_tr[j].r, errr, chi2_cir);
      
    if( (z_h[1] - z_h[0]) * (vec_tr[j].yc - y_h[0]) > 0 )
      vec_tr[j].h = -1;
    else
      vec_tr[j].h = 1;
    
    int n_v_fit = 0;
      
    for(int k = 0; k < n_v; k++)
    {
      if(vec_tr[j].r > z_v[k] - vec_tr[j].zc)
      {
        y_v.push_back(vec_tr[j].yc + vec_tr[j].h * sqrt(vec_tr[j].r*vec_tr[j].r - (z_v[k] - vec_tr[j].zc)*(z_v[k] - vec_tr[j].zc)));
        n_v_fit++;
      }
      else
      {
        break;
      }
    }

	double sin = 0.0;
	double cos = 0.0;

	if(n_v_fit > 0)
	{
  
      cos =  vec_tr[j].h * (y_v[0] - vec_tr[j].yc)/vec_tr[j].r;
      sin = -vec_tr[j].h * (z_v[0] - vec_tr[j].zc)/vec_tr[j].r;
	}
    
    for(int k = 0; k < n_v_fit; k++)
    {
        rho.push_back(z_v[k] * cos + y_v[k] * sin);
    }
    
    x_v.resize(n_v_fit);
    
    int ret2 = fitLinear(n_v_fit, rho, x_v, vec_tr[j].a, vec_tr[j].b, cov, chi2_lin);  
    
    //std::cout << vec_tr[j].tid << " " << vec_tr[j].a << " " << vec_tr[j].b << std::endl;
    
    vec_tr[j].z0 = vec_tr[j].digits.front().z;
    vec_tr[j].y0 = vec_tr[j].yc + vec_tr[j].h * TMath::Sqrt(vec_tr[j].r*vec_tr[j].r - (vec_tr[j].z0 - vec_tr[j].zc)*(vec_tr[j].z0 - vec_tr[j].zc));
    vec_tr[j].x0 = vec_tr[j].a + vec_tr[j].b * (vec_tr[j].z0 * cos + vec_tr[j].y0 * sin);
    vec_tr[j].t0 = vec_tr[j].digits.front().t;
    
    vec_tr[j].ret_ln = ret2;
    vec_tr[j].chi2_ln = chi2_lin;
    vec_tr[j].ret_cr = ret1;
    vec_tr[j].chi2_cr = chi2_cir;
  }
}

bool IsContiguous(const cell& c1, const cell& c2)
{
  if(c1.mod == c2.mod)
  {
    if(TMath::Abs(c1.lay - c2.lay) <= 1)
    {
      if(TMath::Abs(c1.cel - c2.cel) <= 1)
      {
        return true;
      }
    }
  }
  else if(TMath::Abs(c1.mod - c2.mod) == 1)
  {
    if(TMath::Abs(c1.lay - c2.lay) <= 1)
    {
      if((c1.cel == 0 && c2.cel == 11) || (c2.cel == 0 && c1.cel == 11))
      {
        return true;
      }
    }
  }
  return false;
}

bool IsContiguous(cluster cl, const cell& c)
{
  for(unsigned int i = 0; i < cl.cells.size(); i++)
  {
    if(IsContiguous(cl.cells.at(i), c))
    {
      return true;
    }
  }
  return false;
}

void PreCluster(std::vector<cell>*  vec_cell, std::vector<cluster>& vec_precl)
{
  vec_precl.clear();

  std::vector<cell> vec_tmpcell(*vec_cell);
  
  while(vec_tmpcell.size() != 0)
  {
    cell c = vec_tmpcell.front();
    
    //std::cout << vec_tmpcell.size() << std::endl;
    
    bool found = false;
    
    for(unsigned int i = 0; i < vec_precl.size(); i++)
    {
      if(IsContiguous(vec_precl.at(i), c))
      {
        vec_precl.at(i).cells.push_back(c);
        found = true;
        break;
      }
    }
    
    if(found == false)
    {
      cluster cl;
      cl.cells.push_back(c);
      vec_precl.push_back(cl);
    }
    vec_tmpcell.erase(vec_tmpcell.begin());
  }
  
  //std::cout << vec_precl.size() << std::endl;
}

void Filter(std::vector<cluster>& vec_cl)
{
  for(unsigned int i = 0; i < vec_cl.size(); i++)
  {
    for(unsigned int k = 0; k < vec_cl.at(i).cells.size(); k++)
    {
      if(vec_cl.at(i).cells.at(k).adc1 == 0. || vec_cl.at(i).cells.at(k).adc2 == 0.)
      {
        vec_cl.at(i).cells.erase(vec_cl.at(i).cells.begin() + k);
      }
    }
  }
}

double TfromTDC(double t1, double t2, double L)
{
  return 0.5 * (t1 + t2 - ns_Digit::vlfb * L / m_to_mm );
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

void Merge(std::vector<cluster>& vec_cl)
{
  for(unsigned int i = 0; i < vec_cl.size(); i++)
  {
    vec_cl.at(i).e = 0.;
    vec_cl.at(i).x = 0.;
    vec_cl.at(i).y = 0.;
    vec_cl.at(i).z = 0.;
    vec_cl.at(i).t = 0.;
    
    std::vector<double> xv(5);
    std::vector<double> yv(5);
    std::vector<double> zv(5);
    std::vector<double> tv(5);
    std::vector<double> elayer(5);
    
    double x2_cl_mean = 0.;
    double y2_cl_mean = 0.;
    double z2_cl_mean = 0.;    
    
    for(unsigned int k = 0; k < vec_cl.at(i).cells.size(); k++)
    {      
      double x,y,z,t,e;
      CellXYZTE(vec_cl.at(i).cells.at(k), x, y, z, t, e);
      
      // time reference at center of first layer
      // t -= ( (ns_Digit::dzlay[0] + ns_Digit::dzlay[1]) - 
      //       (ns_Digit::dzlay[vec_cl.at(i).cells.at(k).lay] + ns_Digit::dzlay[vec_cl.at(i).cells.at(k).lay + 1]) ) 
      //       * 3.335640951981520495756E-3;
      
      vec_cl.at(i).e += e;
      vec_cl.at(i).x += e * x;
      vec_cl.at(i).y += e * y;
      vec_cl.at(i).z += e * z;
      vec_cl.at(i).t += e * t;
      
      elayer[vec_cl.at(i).cells.at(k).lay] += e; 
        
      xv[vec_cl.at(i).cells.at(k).lay] += e * x;
      yv[vec_cl.at(i).cells.at(k).lay] += e * y;
      zv[vec_cl.at(i).cells.at(k).lay] += e * z;
      tv[vec_cl.at(i).cells.at(k).lay] += e * t;
      
      x2_cl_mean += x*x;
      y2_cl_mean += y*y;
      z2_cl_mean += z*z;
    }
    
    vec_cl.at(i).x /= vec_cl.at(i).e;
    vec_cl.at(i).y /= vec_cl.at(i).e;
    vec_cl.at(i).z /= vec_cl.at(i).e;
    vec_cl.at(i).t /= vec_cl.at(i).e;
    
    x2_cl_mean /= vec_cl.at(i).cells.size();
    y2_cl_mean /= vec_cl.at(i).cells.size();
    z2_cl_mean /= vec_cl.at(i).cells.size();
    
    vec_cl.at(i).varx = x2_cl_mean - vec_cl.at(i).x*vec_cl.at(i).x;
    vec_cl.at(i).vary = y2_cl_mean - vec_cl.at(i).y*vec_cl.at(i).y;
    vec_cl.at(i).varz = z2_cl_mean - vec_cl.at(i).z*vec_cl.at(i).z;
    
    double sx = 0;
    double sy = 0;
    double sz = 0;
    double sxz = 0;
    double syz = 0;
    double sx2 = 0;
    double sy2 = 0;
    double sz2 = 0;
    int nlayer = 0;
    
    double tinner = -999.;
    double touter = -999.;
    double zinner = -999.;
    double zouter = -999.;
    double dir = 0.;
    
    for(int i = 0; i < 5; i++)
    {
      if(elayer[i] != 0.)
      {        
        xv[i] /= elayer[i];
        yv[i] /= elayer[i];
        zv[i] /= elayer[i];
        tv[i] /= elayer[i];
        
        if(tinner < 0)
        {
          tinner = tv[i]; 
          zinner = zv[i]; 
        }
        
        touter = tv[i];
        zouter = zv[i];
        
        sx += xv[i];
        sy += yv[i];
        sz += zv[i];
        sxz += xv[i]*zv[i];
        syz += yv[i]*zv[i];
        sx2 += xv[i]*xv[i];
        sy2 += yv[i]*yv[i];
        sz2 += zv[i]*zv[i];
        
        nlayer++;
      }
    }
    
    if(tinner != touter)
      dir = (zouter - zinner)*(touter - tinner)/TMath::Abs((touter - tinner)*(zouter - zinner));
    else
      dir = 0.;
    
    double ax = (nlayer * sxz - sz*sx)/(nlayer*sz2 - sz*sz);
    double bx = (sx*sz2 - sz*sxz)/(nlayer*sz2 - sz*sz);
    double ay = (nlayer * syz - sz*sy)/(nlayer*sz2 - sz*sz);
    double by = (sy*sz2 - sz*syz)/(nlayer*sz2 - sz*sz);
    
    double mod = TMath::Sqrt(1+ax*ax+ay*ay);
    
    vec_cl.at(i).sx = dir * ax/mod;
    vec_cl.at(i).sy = dir * ay/mod;
    vec_cl.at(i).sz = dir * 1./mod;
  }
}

bool value_comparer(std::map<int, int>::value_type &i1, std::map<int, int>::value_type &i2)
{
  return i1.second < i2.second;
}

void PidBasedClustering(TG4Event* ev, std::vector<cell>* vec_cell, std::vector<cluster>& vec_cl)
{
  std::vector<int> pid(vec_cell->size());
  std::map<int, int> hit_pid;
  
  for(unsigned int i = 0; i < vec_cell->size(); i++)
  {
    // find particle corresponding to more p.e.
    hit_pid.clear();
    
    for(unsigned int j = 0; j < vec_cell->at(i).hindex1.size(); j++)
    {
      hit_pid[ev->SegmentDetectors["EMCalSci"].at(vec_cell->at(i).hindex1.at(j)).PrimaryId]++;
    }
    
    for(unsigned int j = 0; j < vec_cell->at(i).hindex2.size(); j++)
    {
      hit_pid[ev->SegmentDetectors["EMCalSci"].at(vec_cell->at(i).hindex2.at(j)).PrimaryId]++;
    }
    
    pid[i] = std::max_element(hit_pid.begin(), hit_pid.end(), value_comparer)->first;
  }
  
  std::vector<int> unique_pid = pid;
  std::sort(unique_pid.begin(), unique_pid.end());
  std::vector<int>::iterator last = std::unique(unique_pid.begin(), unique_pid.end());
  unique_pid.erase(last, unique_pid.end());
  
  for(unsigned int i = 0; i < unique_pid.size(); i++)
  {
    cluster cl;
    cl.tid = unique_pid[i];
    
    for(unsigned int j = 0; j < pid.size(); j++)
    {
      if(pid[j] == unique_pid[i])
      {
        if(vec_cell->at(j).adc1 == 0 || vec_cell->at(j).adc2 == 0)
          continue;
      
        cl.cells.push_back(vec_cell->at(j));
      }
    }
    if(cl.cells.size() != 0) vec_cl.push_back(cl);
  }
}

void Reconstruct(const char* fDigit, const char* fTrueMC, const char* fOut)
{
  TChain tDigit("tDigit");
  TChain tTrueMC("EDepSimEvents");
  tDigit.Add(fDigit);
  tTrueMC.Add(fTrueMC);
  
  tDigit.AddFriend(&tTrueMC);
  
  TChain* t = &tDigit;
  TFile f(t->GetListOfFiles()->At(0)->GetTitle());
  TGeoManager* geo = (TGeoManager*) f.Get("EDepSimGeometry"); 
  TDirectoryFile* dirfile = (TDirectoryFile*) f.Get("DetSimPassThru");
  
  TG4Event* ev = new TG4Event;
  t->SetBranchAddress("Event",&ev);
  
  std::vector<digit>* vec_digi = new std::vector<digit>;
  std::vector<cell>*  vec_cell = new std::vector<cell>;
  
  t->SetBranchAddress("Stt",&vec_digi);
  t->SetBranchAddress("cell",&vec_cell);
  
  std::vector<track> vec_tr;
  std::vector<cluster> vec_cl;
    
  TFile fout(fOut,"RECREATE");
  TTree tout("tReco","tReco");
  tout.Branch("track","std::vector<track>",&vec_tr);
  tout.Branch("cluster","std::vector<cluster>",&vec_cl);
    
  const int nev = t->GetEntries();
  
  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;
  
  for(int i = 0; i < nev; i++)
  {
    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;
    
    t->GetEntry(i);
    
    vec_tr.clear();
    vec_cl.clear();
    
    TrackFind(ev, vec_digi, vec_tr);
    TrackFit(vec_tr);
    //PreCluster(vec_cell, vec_cl);
    //Filter(vec_cl);
    PidBasedClustering(ev, vec_cell, vec_cl);
    Merge(vec_cl);
    tout.Fill();
  }
  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;
  
  vec_tr.clear();
  vec_cl.clear();
  
  vec_digi->clear();
  vec_cell->clear();
  
  delete vec_digi;
  delete vec_cell;
  
  fout.cd();
  tout.Write();
  geo->Write();
  fout.Close();
}

void help_reco()
{
  std::cout << "Reconstruct(const char* fDigit, const char* fTrueMC, const char* fOut)" << std::endl;
  std::cout << "input file names could contain wild card" << std::endl;
} 

int main(int argc, char* argv[])
{
  if(argc != 4)
    help_reco();
  else
    Reconstruct(argv[1], argv[2], argv[3]);
}
