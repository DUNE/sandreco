#include <TChain.h>
#include <TFile.h>

#include "/mnt/nas01/users/mtenti/wd/analysis/KLOEcal/loader/loader.C"

#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4Event.h"
#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4HitSegment.h"

#include <iostream>
#include <algorithm>

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
        S1 += 1;
        SX += x.at(i);
        SY += y.at(i);
        SYX += x.at(i)*y.at(i);
        SXX += x.at(i)*x.at(i);
    }
    
    double D = S1*SXX - SX*SX;
    
    a = (SY*SXX - SX*SYX)/D;  // i.p. at x = 0
    b = (S1*SYX - SX*SY)/D;   // tg(theta)
    
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

void RecoSTTTrack(const char* fSttDigi, const char* fTrueMC, const char* fOut)
{
  TChain tSttDigi("tSttDigi");
  TChain tTrueMC("EDepSimEvents");
  tSttDigi.Add(fSttDigi);
  tTrueMC.Add(fTrueMC);
  
  tSttDigi.AddFriend(&tTrueMC);
  
  TChain* t = &tSttDigi;
  
  TG4Event* ev = new TG4Event;
  t->SetBranchAddress("Event",&ev);
  
  std::vector<digit>* digit_vec = new std::vector<digit>;
  t->SetBranchAddress("Stt",&digit_vec);
  
  std::vector<track> vec_tr;
    
  TFile fout(fOut,"RECREATE");
  TTree ttr("tTrack","Track");
  ttr.Branch("track","std::vector<track>",&vec_tr);
    
  const int nev = t->GetEntries();
  
  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;
  
  for(int i = 0; i < nev; i++)
  {
    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;
    
    t->GetEntry(i);
    
    /*
    std::cout << "===================================" << std::endl;
    std::cout << i << std::endl;
    std::cout << "===================================" << std::endl;
    */
    
    vec_tr.clear();
    
    //for(unsigned int j = 0; j < ev->Primaries[0].Particles.size(); j++)
    for(unsigned int j = 0; j < ev->Trajectories.size(); j++)
    {
      track tr;
            
      tr.tid = ev->Trajectories.at(j).TrackId;
      
      std::vector<double> y_h;
      std::vector<double> z_h;
      std::vector<double> x_v;
      std::vector<double> y_v;
      std::vector<double> z_v;
      std::vector<double> rho;
      
      for(unsigned int k = 0; k < digit_vec->size(); k++)
      {
        for(unsigned int m = 0; m < digit_vec->at(k).hindex.size(); m++)
        {
          const TG4HitSegment& hseg = ev->SegmentDetectors["StrawTracker"].at(digit_vec->at(k).hindex.at(m)); 
                 
          //if(hseg.PrimaryId == tr.tid)
          //{
            //if(ishitok(ev, hseg->PrimaryId, hseg))
            if(ishitok(ev, tr.tid, hseg))
            {
              tr.digits.push_back(digit_vec->at(k));
              break;
            }
          //}
        }
      }
      
      std::sort(tr.digits.begin(), tr.digits.end(), isDigBefore);
      
      for(unsigned int k = 0; k < tr.digits.size(); k++)
      {
        if(tr.digits.at(k).hor)
        {
          y_h.push_back(tr.digits.at(k).y);
          z_h.push_back(tr.digits.at(k).z);
        }
        else
        {
          x_v.push_back(tr.digits.at(k).x);
          z_v.push_back(tr.digits.at(k).z);
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
      
      int ret1 = fitCircle(n_h, z_h, y_h, tr.zc, tr.yc, tr.r, errr, chi2_cir);
        
      if( (z_h[1] - z_h[0]) * (tr.yc - y_h[0]) > 0 )
        tr.h = -1;
      else
        tr.h = 1;
        
      for(int k = 0; k < n_v; k++)
      {
        y_v.push_back(tr.yc + tr.h * sqrt(tr.r*tr.r - (z_v[k] - tr.zc)*(z_v[k] - tr.zc)));
      }
    
      double cos =  tr.h * (y_v[0] - tr.yc)/tr.r;
      double sin = -tr.h * (z_v[0] - tr.zc)/tr.r;
      
      for(int k = 0; k < n_v; k++)
      {
          rho.push_back(z_v[k] * cos + y_v[i] * sin);
      }
      
      x_v.resize(n_v);
      
      int ret2 = fitLinear(n_v, rho, x_v, tr.a, tr.b, cov, chi2_lin);  
      
      tr.z0 = tr.digits.front().z;
      tr.y0 = tr.yc + tr.h * TMath::Sqrt(tr.r*tr.r - (tr.z0 - tr.zc)*(tr.z0 - tr.zc));
      tr.x0 = tr.a + tr.b * (tr.z0 * cos + tr.y0 * sin);
      tr.t0 = tr.digits.front().t;
      
      //std::cout << "ret1: " << ret1 << " ret2: " << ret2 << std::endl;
           
      if(ret1 == 0 && ret2 == 0)
      {
        vec_tr.push_back(tr);
      }
    }
    ttr.Fill();
  }
  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;
  
  fout.cd();
  ttr.Write();
  fout.Close();
}

void Reco()
{
}