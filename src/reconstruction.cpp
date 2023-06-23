#include <TChain.h>
#include <TDirectoryFile.h>
#include <TFile.h>
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TH1D.h>
#include <TRandom3.h>

#include "TG4Event.h"
#include "TG4HitSegment.h"

#include <algorithm>
#include <iostream>

#include "struct.h"
#include "utils.h"
#include <iomanip>

using namespace sand_reco;

void reset(track& tr)
{
  tr.tid = -1;
  tr.yc = 0.;
  tr.zc = 0.;
  tr.r = 0.;
  tr.a = 0.;
  tr.b = 0.;
  tr.h = 0.;
  tr.ysig = 0.;
  tr.x0 = 0.;
  tr.y0 = 0.;
  tr.z0 = 0.;
  tr.t0 = 0.;
  tr.ret_ln = -1;
  tr.chi2_ln = 0.;
  tr.ret_cr = -1;
  tr.chi2_cr = 0.;
  tr.clX.clear();
  tr.clY.clear();
}

void evalUV(double& u, double& v, double zv, double yv, double z, double y)
{
  u = (z - zv);
  v = (yv - y);
  double d = (u * u + v * v);

  u /= d;
  v /= d;
}

void evalPhi(double& phi, double u, double v) { phi = TMath::ATan2(v, u); }

void findNearDigPhi(int& idx, double exp_phi, std::vector<dg_tube>& vd,
                    double zv, double yv)
{
  double dphi = 1E3;

  double phi, u, v;

  for (unsigned int k = 0; k < vd.size(); k++) {
    evalUV(u, v, zv, yv, vd.at(k).z, vd.at(k).y);
    evalPhi(phi, u, v);
    if (std::abs(exp_phi - phi) < dphi) {
      dphi = std::abs(exp_phi - phi);
      idx = k;
    }
  }
}

bool findNearDigX(int& idx, double exp_x, std::vector<dg_tube>& vd, double xvtx)
{
  double dx = 10000.;

  double x;

  idx = -1;

  bool found = false;

  const double xtol_on_vtx = 50.;

  for (unsigned int k = 0; k < vd.size(); k++) {
    x = vd.at(k).x;
    if (std::abs(exp_x - x) < dx &&
        (std::abs(x - xvtx) < std::abs(exp_x - xvtx) ||
         std::abs(x - xvtx) < xtol_on_vtx)) {
      dx = std::abs(exp_x - x);
      idx = k;
      found = true;
    }
  }

  return found;
}

bool ishitok(TG4Event* ev, int trackid, TG4HitSegment hit, double postol = 5.,
             double angtol = 0.3)
{
  double x = 0.5 * (hit.Start.X() + hit.Stop.X());
  double y = 0.5 * (hit.Start.Y() + hit.Stop.Y());
  double z = 0.5 * (hit.Start.Z() + hit.Stop.Z());

  for (unsigned int jj = 0; jj < ev->Trajectories.size(); jj++) {
    // if(ev->Trajectories[jj].TrackId == hit.PrimaryId)
    if (ev->Trajectories[jj].TrackId == trackid) {
      for (unsigned int kk = 0; kk < ev->Trajectories[jj].Points.size() - 1;
           kk++) {
        double dpos =
            mindist(ev->Trajectories[jj].Points[kk].Position.X(),
                    ev->Trajectories[jj].Points[kk].Position.Y(),
                    ev->Trajectories[jj].Points[kk].Position.Z(),
                    ev->Trajectories[jj].Points[kk + 1].Position.X(),
                    ev->Trajectories[jj].Points[kk + 1].Position.Y(),
                    ev->Trajectories[jj].Points[kk + 1].Position.Z(), x, y, z);
        double dang =
            angle(ev->Trajectories[jj].Points[kk + 1].Position.X() -
                      ev->Trajectories[jj].Points[kk].Position.X(),
                  ev->Trajectories[jj].Points[kk + 1].Position.Y() -
                      ev->Trajectories[jj].Points[kk].Position.Y(),
                  ev->Trajectories[jj].Points[kk + 1].Position.Z() -
                      ev->Trajectories[jj].Points[kk].Position.Z(),
                  hit.Stop.X() - hit.Start.X(), hit.Stop.Y() - hit.Start.Y(),
                  hit.Stop.Z() - hit.Start.Z());

        if (dpos < postol && dang < angtol) return true;
      }
    }
  }
  return false;
}

int EvalDirection(const std::vector<double>& z, const std::vector<double>& y,
                  double& zc, double& yc)
{
  double cross_prod = 0.;
  double vz;
  double vy;
  double rz;
  double ry;

  for (unsigned int i = 0; i < z.size() - 1; i++) {
    rz = z[i] - zc;
    ry = y[i] - yc;

    vz = z[i + 1] - z[i];
    vy = y[i + 1] - y[i];

    cross_prod += ry * vz - rz * vy;
  }

  // clockwise direction (as seen from positive x) is positive
  return cross_prod >= 0 ? 1 : -1;
}
/*
int getSignY(const std::vector<double>& z_h, std::vector<double>& y_h,
             const track& tr)
{
  int forward = z_h[1] - z_h[0] > 0 ? 1 : -1;

  double dy = (y_h[0] - tr.yc) * abs(y_h[0] - tr.yc) / tr.r;

  for (unsigned int i = 1; i < z_h.size(); i++) {
    if ((z_h[i] - z_h[i - 1]) * forward > 0.) {
      dy += (y_h[i] - tr.yc) * abs(y_h[i] - tr.yc) / tr.r;
    } else {
      break;
    }
  }

  return dy > 0. ? +1 : -1;
}
*/
void getVertCoord(const std::vector<double>& z_v, std::vector<double>& y_v,
                  int sign, const track& tr)
{
  y_v.clear();

  int forward = z_v[1] - z_v[0] > 0 ? 1 : -1;

  double dy_sq, dy;

  dy_sq = tr.r * tr.r - (z_v[0] - tr.zc) * (z_v[0] - tr.zc);

  if (dy_sq < 0.) dy_sq = 0.;

  dy = TMath::Sqrt(dy_sq);
  y_v.push_back(tr.yc + sign * dy);

  for (unsigned int i = 1; i < z_v.size(); i++) {
    if ((z_v[i] - z_v[i - 1]) * forward >= 0.) {
      dy_sq = tr.r * tr.r - (z_v[i] - tr.zc) * (z_v[i] - tr.zc);

      if (dy_sq < 0.) dy_sq = 0.;

      dy = TMath::Sqrt(dy_sq);

      y_v.push_back(tr.yc + sign * dy);
    } else {
      break;
    }
  }
}

void GetRho(const std::vector<double>& z_v, const std::vector<double>& y_v,
            std::vector<double>& rho, const track& tr, double cos, double sin)
{
  rho.clear();

  for (unsigned int k = 0; k < y_v.size(); k++) {
    rho.push_back(z_v[k] * cos + y_v[k] * sin);
    // if(isnan(z_v[k] * cos + y_v[k] * sin))
    //{
    //  std::cout << z_v[k] << " " << cos << " " << y_v[k] << " " << sin <<
    // std::endl;
    //}
  }
}

void FillPositionInfo(track& tr, int signy, double cos, double sin)
{
  dg_tube d = (sand_reco::stt::isDigBefore(tr.clX.front(), tr.clY.front())
                   ? tr.clX.front()
                   : tr.clY.front());
  tr.z0 = d.z;
  tr.t0 = d.tdc;
  if (d.hor == true) {
    tr.y0 = d.y;
    tr.x0 = tr.a + tr.b * (tr.z0 * cos + tr.y0 * sin);
  } else {
    tr.x0 = d.x;
    double dy2 = tr.r * tr.r - (tr.z0 - tr.zc) * (tr.z0 - tr.zc);
    double dy = dy2 < 0. ? 0. : TMath::Sqrt(dy2);
    tr.y0 = tr.yc + signy * dy;
  }
}

int evalYSign(const std::vector<double>& ys, double yc)
{
  auto yres = 0.;

  for (const auto y : ys) yres += y - yc;

  if (yres > 0.)
    return +1;
  else
    return -1;
}

int evalYSign(const track& tr)
{
  std::vector<double> ys;

  for (const auto d : tr.clY) ys.push_back(d.y);

  return evalYSign(ys, tr.yc);
}

int fitCircle(int n, const std::vector<double>& x, const std::vector<double>& y,
              double& xc, double& yc, double& r, double& errr, double& chi2)
{
  xc = -999;
  yc = -999;
  r = -999;
  errr = -999;
  chi2 = -999;

  if (x.size() != y.size()) return 1;

  double sumx = 0, sumy = 0;                            // linear    terms
  double sumx2 = 0, sumy2 = 0, sumxy = 0;               // quadratic terms
  double sumxy2 = 0, sumx2y = 0, sumx3 = 0, sumy3 = 0;  // cubic     terms

  for (int i = 0; i < n; i++) {
    double xp = x.at(i);
    double yp = y.at(i);
    sumx += xp;
    sumy += yp;
    sumx2 += xp * xp;
    sumy2 += yp * yp;
    sumxy += xp * yp;
    sumxy2 += xp * yp * yp;
    sumx2y += xp * xp * yp;
    sumx3 += xp * xp * xp;
    sumy3 += yp * yp * yp;
  }

  double a = n * sumx2 - sumx * sumx;
  double b = n * sumxy - sumx * sumy;
  double c = n * sumy2 - sumy * sumy;
  double d = 0.5 * (n * sumxy2 - sumx * sumy2 + n * sumx3 - sumx * sumx2);
  double e = 0.5 * (n * sumx2y - sumy * sumx2 + n * sumy3 - sumy * sumy2);

  if (a * c - b * b == 0.) return 2;

  xc = (d * c - b * e) / (a * c - b * b);
  yc = (a * e - b * d) / (a * c - b * b);

  double rMean = 0;
  double rrms = 0;

  for (int i = 0; i < n; i++) {
    double xp = x.at(i);
    double yp = y.at(i);
    double r2 = (xp - xc) * (xp - xc) + (yp - yc) * (yp - yc);

    rMean += sqrt(r2);
    rrms += r2;
  }

  rMean /= n;
  rrms /= n;
  r = rMean;

  errr = sqrt(rrms - rMean * rMean);

  chi2 = 0.0;

  for (int i = 0; i < n; i++) {
    chi2 += TMath::Abs((y.at(i) - yc) * (y.at(i) - yc) +
                       (x.at(i) - xc) * (x.at(i) - xc) - r * r);
  }

  chi2 /= n;

  /*
  std::cout << "==== ZY =====" << std::endl;

  for(int i = 0; i < n; i++)
  {
    std::cout << x.at(i) << " " << y.at(i) << std::endl;
  }

  std::cout << "---> " << xc << " " << yc << " " << r << " " << chi2 <<
  std::endl;
  */

  return 0;
}

int fitLinear(int n, const std::vector<double>& x, const std::vector<double>& y,
              double& a, double& b, double cov[][2], double& chi2)
{
  a = -999;
  b = -999;
  chi2 = -999;

  // if (x.size() != y.size()) return 1;

  cov[0][0] = -999;
  ;
  cov[0][1] = -999;
  cov[1][0] = -999;
  cov[1][1] = -999;

  double S1 = 0.;
  double SX = 0.;
  double SY = 0.;
  double SYX = 0.;
  double SXX = 0.;

  for (int i = 0; i < n; i++) {
    // std::cout << "XY: " <<  x.at(i) << " " <<  y.at(i) << std::endl;

    // if(isnan(x.at(i)) || isnan(y.at(i)))
    //  std::cout << x.at(i) << " " << y.at(i) << std::endl;

    S1 += 1.;
    SX += x.at(i);
    SY += y.at(i);
    SYX += x.at(i) * y.at(i);
    SXX += x.at(i) * x.at(i);
  }

  double D = S1 * SXX - SX * SX;

  if (D == 0.) return 2;

  // std::cout << "D: " << D << " " << S1 << " " << SX << " " << SY << " " <<
  // SXX << " " << SYX << std::endl;

  a = (SY * SXX - SX * SYX) / D;  // i.p. at x = 0
  b = (S1 * SYX - SX * SY) / D;   // tg(theta)

  // if(isnan(a) || isnan(b))
  //{
  //  std::cout << a << " " << b << " " << SY << " " << SXX << " " << SX << " "
  // << SYX << " " << S1 << " " << D << std::endl;
  //}

  // std::cout << a << " " << b << std::endl;

  cov[0][0] = SXX / D;
  cov[0][1] = -SX / D;
  cov[1][0] = -SX / D;
  cov[1][1] = S1 / D;

  chi2 = 0.0;

  for (int i = 0; i < n; i++) {
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

void fillInfoCircFit(int n, const std::vector<double>& z,
                     const std::vector<double>& y, track& tr)
{
  double errr;

  tr.ret_cr = fitCircle(n, z, y, tr.zc, tr.yc, tr.r, errr, tr.chi2_cr);

  if (tr.ret_cr != 0) return;

  tr.h = EvalDirection(z, y, tr.zc, tr.yc);
  tr.ysig = evalYSign(y, tr.yc);
}

enum class TrackFilter { all_tracks, only_primaries };

void TrackFind(TG4Event* ev, std::vector<dg_tube>* vec_digi,
               std::vector<track>& vec_tr,
               TrackFilter const track_filter = TrackFilter::all_tracks)
{
  vec_tr.clear();

  // for(unsigned int j = 0; j < ev->Primaries[0].Particles.size(); j++)
  for (unsigned int j = 0; j < ev->Trajectories.size(); j++) {
    // exclude not primaries particles
    if (track_filter == TrackFilter::only_primaries &&
        ev->Trajectories.at(j).ParentId != -1) {
      continue;
    }

    track tr;

    reset(tr);

    tr.tid = ev->Trajectories.at(j).TrackId;

    TRandom3 rand(0);

    std::map<double, dg_tube> time_ordered_XZdigit;
    std::map<double, dg_tube> time_ordered_YZdigit;

    for (unsigned int k = 0; k < vec_digi->size(); k++) {

      std::vector<TG4HitSegment> vhits;

      for (unsigned int m = 0; m < vec_digi->at(k).hindex.size(); m++) {
        const TG4HitSegment& hseg =
            ev->SegmentDetectors["Straw"].at(vec_digi->at(k).hindex.at(m));

        if (hseg.PrimaryId == tr.tid)
          if (ishitok(ev, tr.tid, hseg)) vhits.push_back(hseg);
      }

      if (vhits.size() > 0u) {
        std::sort(vhits.begin(), vhits.end(),
                  [](const TG4HitSegment& h1, const TG4HitSegment& h2) {
                    return (h1.GetStart().T() + h1.GetStop().T()) <
                           (h2.GetStart().T() + h2.GetStop().T());
                  });
        auto& fdig = vhits.front().GetStart();
        auto& ldig = vhits.back().GetStop();

        auto t = 0.5 * (ldig.T() + fdig.T());
        if (vec_digi->at(k).hor == 0) {
          vec_digi->at(k).x = fdig.X() +
                              (ldig.X() - fdig.X()) / (ldig.Z() - fdig.Z()) *
                                  (vec_digi->at(k).z - fdig.Z()) +
                              rand.Gaus(0, 0.2);
          time_ordered_XZdigit[t] = vec_digi->at(k);
        } else {
          vec_digi->at(k).y = fdig.Y() +
                              (ldig.Y() - fdig.Y()) / (ldig.Z() - fdig.Z()) *
                                  (vec_digi->at(k).z - fdig.Z()) +
                              rand.Gaus(0, 0.2);
          time_ordered_YZdigit[t] = vec_digi->at(k);
        }
      }
    }

    for (const auto& p : time_ordered_XZdigit) tr.clX.push_back(p.second);
    for (const auto& p : time_ordered_YZdigit) tr.clY.push_back(p.second);

    vec_tr.push_back(tr);
  }
}

void fillLayers(std::map<int, std::vector<dg_tube> >& m,
                std::vector<dg_tube>& d, TH1D& hdummy, int hor)
{
  for (unsigned int k = 0; k < d.size(); k++) {
    if (d.at(k).hor == hor)
      m[-1 * hdummy.FindBin(d.at(k).z)].push_back(d.at(k));
  }
}

void findTracksY(std::map<int, std::vector<dg_tube> >& mdY,
                 std::vector<std::vector<dg_tube> >& clustersY, TH1D& hdummy,
                 double xvtx_reco, double yvtx_reco, double zvtx_reco,
                 double phi_tol = 0.1, double dlay_tol = 5)
{
  int idx = 0;
  int prev_mod;
  double prev_phi, exp_phi;
  double prev_z, slope;

  double u, v, phi, dphi;

  dg_tube current_digit;

  // loop on modules
  while (mdY.size() != 0) {
    // get most downstream module
    std::map<int, std::vector<dg_tube> >::iterator ite = mdY.begin();
    std::vector<dg_tube>* layer = &(ite->second);

    // loop on digits in the module
    while (layer->size() != 0) {
      // get forst digit
      current_digit = layer->front();

      // create cluster and insert first digit
      std::vector<dg_tube> clY;
      clY.push_back(std::move(current_digit));

      // remove the current digit from the module
      layer->erase(layer->begin());

      // get module id
      prev_mod = ite->first;

      // get iterator to the next module
      std::map<int, std::vector<dg_tube> >::iterator nite = std::next(ite);

      // reset parameters
      evalUV(u, v, zvtx_reco, yvtx_reco, current_digit.z, current_digit.y);
      evalPhi(phi, u, v);
      slope = 0.;
      prev_phi = phi;
      prev_z = hdummy.GetBinCenter(-1 * ite->first);
      prev_mod = ite->first;

      // loop on upstream modules
      while (nite != mdY.end()) {
        // check if module are downstream the reco vtx
        if (hdummy.GetXaxis()->GetBinUpEdge(-1 * nite->first) > zvtx_reco) {
          // get next layer
          std::vector<dg_tube>* nlayer = &(nite->second);

          // evaluate the distance (in number of modules) between this module
          // and the previous one
          int dlayer = nite->first - prev_mod;

          // check if the distance is within the tolerance
          if (dlayer <= dlay_tol) {
            // eval prediction
            exp_phi = slope * (prev_z - hdummy.GetBinCenter(-1 * nite->first)) +
                      prev_phi;

            // find nearest digit
            findNearDigPhi(idx, exp_phi, *nlayer, zvtx_reco, yvtx_reco);

            // eval phi
            evalUV(u, v, zvtx_reco, yvtx_reco, nlayer->at(idx).z,
                   nlayer->at(idx).y);
            evalPhi(phi, u, v);

            // eval residual
            dphi = exp_phi - phi;

            // check if distance is within tolerance
            if (std::abs(dphi) <= phi_tol) {
              // get selected digit
              current_digit = nlayer->at(idx);

              // add to the cluster
              clY.push_back(std::move(current_digit));

              // update parameter
              prev_mod = nite->first;
              slope = (prev_phi - phi) /
                      (prev_z - hdummy.GetBinCenter(-1 * nite->first));
              prev_phi = phi;
              prev_z = hdummy.GetBinCenter(-1 * nite->first);

              // remove from the module
              nlayer->erase(nlayer->begin() + idx);

              // remove module if it is empty
              if (nlayer->size() == 0) {
                std::map<int, std::vector<dg_tube> >::iterator dummy =
                    std::next(nite);
                mdY.erase(nite);
                nite = dummy;
                continue;
              }
            }
          }
        }
        // go to upstream module
        nite = std::next(nite);
      }
      // save cluster
      clustersY.push_back(clY);
    }
    // remove module from map
    mdY.erase(ite);
  }
}

void findTracksX(std::map<int, std::vector<dg_tube> >& mdX,
                 std::vector<std::vector<dg_tube> >& clustersX, TH1D& hdummy,
                 double xvtx_reco, double yvtx_reco, double zvtx_reco,
                 double x_tol = 100, double dlay_tol = 5)
{
  int idx = 0;
  int prev_mod;

  double prev_x, dx;

  dg_tube current_digit;

  // loop on modules from downstream to upstream
  while (mdX.size() != 0) {
    // get most downstream module
    std::map<int, std::vector<dg_tube> >::iterator ite = mdX.begin();
    std::vector<dg_tube>* layer = &(ite->second);

    // loop on digits in the module
    while (layer->size() != 0) {
      // get forst digit
      current_digit = layer->front();

      // create cluster and insert first digit
      std::vector<dg_tube> clX;
      clX.push_back(std::move(current_digit));

      // remove the current digit from the module
      layer->erase(layer->begin());

      // get module id
      prev_mod = ite->first;

      // get iterator to the next module
      std::map<int, std::vector<dg_tube> >::iterator nite = std::next(ite);

      // reset parameters
      prev_x = current_digit.x;
      prev_mod = ite->first;

      // loop on upstream modules
      while (nite != mdX.end()) {
        // check if module are downstream the reco vtx
        if (hdummy.GetXaxis()->GetBinUpEdge(-1 * nite->first) > zvtx_reco) {
          // get next layer
          std::vector<dg_tube>* nlayer = &(nite->second);

          // evaluate the distance (in number of modules) between this module
          // and the previous one
          int dlayer = nite->first - prev_mod;

          // check if the distance is within the tolerance
          if (dlayer <= dlay_tol) {
            // find nearest digit
            if (findNearDigX(idx, prev_x, *nlayer, xvtx_reco) == true) {
              // eval residual
              dx = prev_x - nlayer->at(idx).x;

              // check if distance is within tolerance
              if (std::abs(dx) <= x_tol) {
                // get selected digit
                current_digit = nlayer->at(idx);

                // add to the cluster
                clX.push_back(std::move(current_digit));

                // update parameter
                prev_mod = nite->first;
                prev_x = nlayer->at(idx).x;

                // remove from the module
                nlayer->erase(nlayer->begin() + idx);

                // remove module if it is empty
                if (nlayer->size() == 0) {
                  std::map<int, std::vector<dg_tube> >::iterator dummy =
                      std::next(nite);
                  mdX.erase(nite);
                  nite = dummy;
                  continue;
                }
              }
            }
          }
        }
        // go to upstream module
        nite = std::next(nite);
      }
      // save cluster
      clustersX.push_back(clX);
    }
    // remove module from map
    mdX.erase(ite);
  }
}

void mergeXYTracks(std::vector<std::vector<dg_tube> >& clustersX,
                   std::vector<std::vector<dg_tube> >& clustersY,
                   std::vector<track>& tr3D, double dn_tol, double dz_tol)
{
  double dzend;
  int index = 0;

  // loop on Y clusters
  for (unsigned int jj = 0; jj < clustersY.size(); jj++) {
    // resent downstream z residual
    dzend = 10000.;

    // loop on X clusters
    for (unsigned int kk = 0; kk < clustersX.size(); kk++) {
      // evaluate the percentege difference in number of digits: abs(dn)/mn
      double dn = int(clustersY.at(jj).size()) - int(clustersX.at(kk).size());
      double mn = 0.5 * (clustersY.at(jj).size() + clustersX.at(kk).size());

      // search cluster X with abs(dn)/mn less than tolerance and
      // residual in the most downstream digit z less than tolerance
      if (std::abs(dn) / mn < dn_tol &&
          std::abs(clustersY.at(jj).back().z - clustersX.at(kk).back().z) <
              dz_tol &&
          std::abs(clustersY.at(jj).front().z - clustersX.at(kk).front().z) <
              dz_tol) {
        // find cluster X best matching Y cluster in term of
        // residual in the most downstream and upstream digit z
        if (std::abs(clustersY.at(jj).back().z - clustersX.at(kk).back().z) +
                std::abs(clustersY.at(jj).front().z -
                         clustersX.at(kk).front().z) <
            dzend) {
          dzend =
              std::abs(clustersY.at(jj).back().z - clustersX.at(kk).back().z) +
              std::abs(clustersY.at(jj).front().z - clustersX.at(kk).front().z);
          track tr;
          tr.tid = index++;
          tr.clX = std::move(clustersX.at(kk));
          tr.clY = std::move(clustersY.at(jj));
          clustersY.erase(clustersY.begin() + jj--);
          clustersX.erase(clustersX.begin() + kk--);

          tr3D.push_back(std::move(tr));
          break;
        }
      }
    }
  }
}

void TrackFind(std::vector<track>& tracks, std::vector<dg_tube> digits,
               std::vector<double>& binning, double xvtx_reco, double yvtx_reco,
               double zvtx_reco, double tol_phi, double tol_x, int tol_mod,
               unsigned int mindigtr, const double dn_tol, const double dz_tol)
{
  TH1D hdummy("hdummy", "hdummy;Z (mm); multipliciy", binning.size() - 1,
              binning.data());
  std::vector<std::vector<dg_tube> > clustersY;
  std::vector<std::vector<dg_tube> > clustersX;

  // track finding with clustering in arctg(v/u) VS z
  std::map<int, std::vector<dg_tube> > mdY;
  fillLayers(mdY, digits, hdummy, 1);
  clustersY.clear();

  // find track on XZ view
  std::map<int, std::vector<dg_tube> > mdX;
  fillLayers(mdX, digits, hdummy, 0);
  clustersX.clear();

  // tolerance in phi: 0.1 rad
  // tolerance in module: 5
  findTracksY(mdY, clustersY, hdummy, xvtx_reco, yvtx_reco, zvtx_reco, tol_phi,
              tol_mod);

  // tolerance in x: 10 cm
  // tolerance in module: 5
  findTracksX(mdX, clustersX, hdummy, xvtx_reco, yvtx_reco, zvtx_reco, tol_x,
              tol_mod);

  // remove clusters with less than 3 digits
  for (unsigned int nn = 0; nn < clustersY.size(); nn++) {
    if (clustersY.at(nn).size() <= mindigtr) {
      clustersY.erase(clustersY.begin() + nn);
      nn--;
    }
  }

  // remove clusters with less than 3 digits
  for (unsigned int nn = 0; nn < clustersX.size(); nn++) {
    if (clustersX.at(nn).size() <= mindigtr) {
      clustersX.erase(clustersX.begin() + nn);
      nn--;
    }
  }

  // order tracks by length
  std::sort(clustersY.begin(), clustersY.end(), sand_reco::ecal::isCluBigger);

  // order tracks by length
  std::sort(clustersX.begin(), clustersX.end(), sand_reco::ecal::isCluBigger);

  // order digits: upstream first
  for (unsigned int jj = 0; jj < clustersX.size(); jj++) {
    std::sort(clustersX.at(jj).begin(), clustersX.at(jj).end(),
              sand_reco::stt::isDigUpstream);
  }

  for (unsigned int jj = 0; jj < clustersY.size(); jj++) {
    std::sort(clustersY.at(jj).begin(), clustersY.at(jj).end(),
              sand_reco::stt::isDigUpstream);
  }

  // merge XZ and YZ clusters
  // tolerance in dn/n: 0.7
  // tolerance dz: 20 cm
  mergeXYTracks(clustersX, clustersY, tracks, dn_tol, dz_tol);
}

void TrackFit(std::vector<track>& vec_tr)
{

  std::vector<double> y_h;
  std::vector<double> z_h;
  std::vector<double> x_v;
  std::vector<double> y_v;
  std::vector<double> z_v;
  std::vector<double> rho;

  for (auto& tr : vec_tr) {

    y_h.clear();
    z_h.clear();
    x_v.clear();
    y_v.clear();
    z_v.clear();
    rho.clear();

    for (auto const& d : tr.clX) {
      x_v.push_back(d.x);
      z_v.push_back(d.z);
    }

    for (auto const& d : tr.clY) {
      y_h.push_back(d.y);
      z_h.push_back(d.z);
    }

    double chi2_cir, chi2_lin;
    double errr;

    double cov[2][2];

    int n_h = z_h.size();
    int n_v = z_v.size();

    // std::cout << n_h << " " << n_v << std::endl;

    if (n_h <= 2) {
      tr.ret_cr = -2;
      continue;
    }

    fillInfoCircFit(n_h, z_h, y_h, tr);

    if (tr.ret_cr != 0) continue;

    if (n_v <= 2) {
      tr.ret_ln = -2;
      continue;
    }

    int quadrant;

    if (z_h[0] - tr.zc >= 0)
      quadrant = 1;
    else
      quadrant = 2;

    if (y_h[0] - tr.yc < 0) quadrant *= -1;

    int signy = quadrant > 0 ? 1 : -1;

    // int signy = getSignY(z_h, y_h, tr);

    getVertCoord(z_v, y_v, signy, tr);

    if (y_v.size() <= 2) {
      tr.ret_ln = -3;
      continue;
    }

    double cos = tr.h * (y_v[0] - tr.yc) / tr.r;
    double sin = -tr.h * (z_v[0] - tr.zc) / tr.r;

    GetRho(z_v, y_v, rho, tr, cos, sin);

    tr.ret_ln = fitLinear(y_v.size(), rho, x_v, tr.a, tr.b, cov, tr.chi2_ln);

    if (tr.ret_ln != 0) continue;

    // std::cout << tr.tid << " " << tr.a << " " << tr.b <<
    // std::endl;

    FillPositionInfo(tr, signy, cos, sin);

    // if(isnan(tr.x0))
    //{
    //  std::cout << tr.a  << " " << tr.b << " " << tr.z0
    // << " " << cos << " " << tr.y0 << " " <<  sin << std::endl;
    //  std::cout << tr.yc << " " << tr.h << " " << tr.r <<
    // " " << tr.zc << " " << dz << " " << y_v[0] << std::endl;
    //}
  }
}

void fitCircle(std::vector<track>& tr3D, TH1D& hdummy)
{

  double errr;

  std::vector<double> z;
  std::vector<double> y;

  // loop on clusters
  for (unsigned int nn = 0; nn < tr3D.size(); nn++) {
    // reset z, y digit position vector
    z.clear();
    y.clear();

    // loop on the digits of clusters
    for (unsigned int tt = 0; tt < tr3D.at(nn).clY.size(); tt++) {
      // fill z, y digit position vector
      z.push_back(tr3D.at(nn).clY.at(tt).z);
      y.push_back(tr3D.at(nn).clY.at(tt).y);
    }
    // fit and save results
    tr3D.at(nn).ret_cr =
        fitCircle(tr3D.at(nn).clY.size(), z, y, tr3D.at(nn).zc, tr3D.at(nn).yc,
                  tr3D.at(nn).r, errr, tr3D.at(nn).chi2_cr);
  }

  // evaluate if the track is up or down its y center and the direction in the
  // circle
  // evaluate particle direction in the circle
  double prodVect = 0.;
  double yres = 0.;
  double ry, rz;
  double vely, velz;

  for (unsigned int nn = 0; nn < tr3D.size(); nn++) {
    prodVect = 0.;
    yres = 0.;

    for (unsigned int tt = 0; tt < tr3D.at(nn).clY.size() - 1; tt++) {
      yres += (tr3D.at(nn).clY.at(tt).y - tr3D.at(nn).yc);

      if (tr3D.at(nn).clY.at(tt + 1).z < tr3D.at(nn).clY.at(tt).z)
        std::cout << "problem: " << tr3D.at(nn).clY.at(tt).z << " "
                  << tr3D.at(nn).clY.at(tt + 1).z << std::endl;

      ry = tr3D.at(nn).clY.at(tt).y - tr3D.at(nn).yc;
      rz = tr3D.at(nn).clY.at(tt).z - tr3D.at(nn).zc;
      vely = tr3D.at(nn).clY.at(tt + 1).y - tr3D.at(nn).clY.at(tt).y;
      velz = tr3D.at(nn).clY.at(tt + 1).z - tr3D.at(nn).clY.at(tt).z;

      prodVect += ry * velz - rz * vely;
    }
    if (prodVect != 0.) prodVect /= std::abs(prodVect);
    tr3D.at(nn).h = int(prodVect);

    tr3D.at(nn).ysig = evalYSign(tr3D.at(nn));
  }
}

void fitLine(std::vector<track>& tr3D, double xvtx_reco, double yvtx_reco,
             double zvtx_reco)
{
  std::vector<double> vrho;
  std::vector<double> vx;
  double rho, yexp, x, phi0, sin, cos;

  double cov[2][2];

  for (unsigned int i = 0; i < tr3D.size(); i++) {
    vrho.clear();
    vx.clear();

    // evaluate phi0
    phi0 = TMath::ATan2(yvtx_reco - tr3D.at(i).yc, zvtx_reco - tr3D.at(i).zc);

    // evalute sin and cos
    // http://www2.fisica.unimi.it/andreazz/AA_TrackingSystems.pdf
    sin = -1 * tr3D.at(i).h * TMath::Cos(phi0);
    cos = tr3D.at(i).h * TMath::Sin(phi0);

    // evaluate rho and x
    for (unsigned int kk = 0; kk < tr3D.at(i).clX.size(); kk++) {
      double radq;
      if (std::abs(tr3D.at(i).clX.at(kk).z - tr3D.at(i).zc) > tr3D.at(i).r)
        radq = 0;
      else
        radq = TMath::Sqrt(tr3D.at(i).r * tr3D.at(i).r -
                           (tr3D.at(i).clX.at(kk).z - tr3D.at(i).zc) *
                               (tr3D.at(i).clX.at(kk).z - tr3D.at(i).zc));

      yexp = tr3D.at(i).yc + tr3D.at(i).ysig * radq;

      rho = tr3D.at(i).clX.at(kk).z * cos + yexp * sin;

      vrho.push_back(rho);
      vx.push_back(tr3D.at(i).clX.at(kk).x);
    }

    // linear fit
    tr3D.at(i).ret_ln = fitLinear(vrho.size(), vrho, vx, tr3D.at(i).a,
                                  tr3D.at(i).b, cov, tr3D.at(i).chi2_ln);
  }
}

void fillPosAndTime(std::vector<track>& tracks)
{
  for (unsigned int i = 0; i < tracks.size(); i++) {

    // std::cout << tracks.at(i).clY.front().z << " " <<
    // tracks.at(i).clX.front().z << std::endl;

    if (tracks.at(i).clY.front().z < tracks.at(i).clX.front().z) {
      tracks.at(i).y0 = tracks.at(i).clY.front().y;
      tracks.at(i).z0 = tracks.at(i).clY.front().z;
      tracks.at(i).t0 = tracks.at(i).clY.front().tdc;

      double cos =
          tracks.at(i).h * (tracks.at(i).y0 - tracks.at(i).yc) / tracks.at(i).r;
      double sin = -tracks.at(i).h * (tracks.at(i).z0 - tracks.at(i).zc) /
                   tracks.at(i).r;

      tracks.at(i).x0 =
          tracks.at(i).a +
          tracks.at(i).b * (tracks.at(i).z0 * cos + tracks.at(i).y0 * sin);

      // std::cout << tracks.at(i).x0 << " " << tracks.at(i).a << " " <<
      // tracks.at(i).b << " " << tracks.at(i).z0 << " " << cos << " " <<
      // tracks.at(i).y0 << " " << sin << std::endl;
    } else {
      tracks.at(i).x0 = tracks.at(i).clX.front().x;
      tracks.at(i).z0 = tracks.at(i).clX.front().z;
      tracks.at(i).t0 = tracks.at(i).clX.front().tdc;

      double dy2 = tracks.at(i).r * tracks.at(i).r -
                   (tracks.at(i).z0 - tracks.at(i).zc) *
                       (tracks.at(i).z0 - tracks.at(i).zc);
      double dy = dy2 < 0. ? 0. : TMath::Sqrt(dy2);

      tracks.at(i).y0 = tracks.at(i).yc + tracks.at(i).ysig * dy;
    }
  }
}

void TrackFit(std::vector<track>& tracks, std::vector<double>& binning,
              double xvtx_reco, double yvtx_reco, double zvtx_reco)
{
  TH1D hdummy("hdummy", "hdummy;Z (mm); multipliciy", binning.size() - 1,
              binning.data());

  // circular fit
  fitCircle(tracks, hdummy);

  // linear fit
  fitLine(tracks, xvtx_reco, yvtx_reco, zvtx_reco);

  // track position
  fillPosAndTime(tracks);
}

bool IsContiguous(const dg_cell& c1, const dg_cell& c2)
{
  if (c1.mod == c2.mod) {
    if (TMath::Abs(c1.lay - c2.lay) <= 1) {
      if (TMath::Abs(c1.cel - c2.cel) <= 1) {
        return true;
      }
    }
  } else if (TMath::Abs(c1.mod - c2.mod) == 1) {
    if (TMath::Abs(c1.lay - c2.lay) <= 1) {
      if ((c1.cel == 0 && c2.cel == 11) || (c2.cel == 0 && c1.cel == 11)) {
        return true;
      }
    }
  }
  return false;
}

bool IsContiguous(cluster cl, const dg_cell& c)
{
  for (unsigned int i = 0; i < cl.cells.size(); i++) {
    if (IsContiguous(cl.cells.at(i), c)) {
      return true;
    }
  }
  return false;
}

void PreCluster(std::vector<dg_cell>* vec_cell, std::vector<cluster>& vec_precl)
{
  vec_precl.clear();

  std::vector<dg_cell> vec_tmpcell(*vec_cell);

  while (vec_tmpcell.size() != 0) {
    dg_cell c = vec_tmpcell.front();

    // std::cout << vec_tmpcell.size() << std::endl;

    bool found = false;

    for (unsigned int i = 0; i < vec_precl.size(); i++) {
      if (IsContiguous(vec_precl.at(i), c)) {
        vec_precl.at(i).cells.push_back(c);
        found = true;
        break;
      }
    }

    if (found == false) {
      cluster cl;
      cl.cells.push_back(c);
      vec_precl.push_back(cl);
    }
    vec_tmpcell.erase(vec_tmpcell.begin());
  }

  // std::cout << vec_precl.size() << std::endl;
}

void Filter(std::vector<cluster>& vec_cl)
{
  for (unsigned int i = 0; i < vec_cl.size(); i++) {
    for (unsigned int k = 0; k < vec_cl.at(i).cells.size(); k++) {
      if (vec_cl.at(i).cells.at(k).ps1.at(0).adc == 0. ||
          vec_cl.at(i).cells.at(k).ps2.at(0).adc == 0.) {
        vec_cl.at(i).cells.erase(vec_cl.at(i).cells.begin() + k);
      }
    }
  }
}

void Merge(std::vector<cluster>& vec_cl)
{
  for (unsigned int i = 0; i < vec_cl.size(); i++) {
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

    for (unsigned int k = 0; k < vec_cl.at(i).cells.size(); k++) {
      double x, y, z, t, e;
      sand_reco::ecal::reco::CellXYZTE(vec_cl.at(i).cells.at(k), x, y, z, t, e);

      // time reference at center of first layer
      // t -= ( (dzlay[0] + dzlay[1]) -
      //       (dzlay[vec_cl.at(i).cells.at(k).lay] +
      // dzlay[vec_cl.at(i).cells.at(k).lay + 1]) )
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

      x2_cl_mean += x * x * e;
      y2_cl_mean += y * y * e;
      z2_cl_mean += z * z * e;
    }

    vec_cl.at(i).x /= vec_cl.at(i).e;
    vec_cl.at(i).y /= vec_cl.at(i).e;
    vec_cl.at(i).z /= vec_cl.at(i).e;
    vec_cl.at(i).t /= vec_cl.at(i).e;

    x2_cl_mean /= vec_cl.at(i).e;
    y2_cl_mean /= vec_cl.at(i).e;
    z2_cl_mean /= vec_cl.at(i).e;

    vec_cl.at(i).varx = x2_cl_mean - vec_cl.at(i).x * vec_cl.at(i).x;
    vec_cl.at(i).vary = y2_cl_mean - vec_cl.at(i).y * vec_cl.at(i).y;
    vec_cl.at(i).varz = z2_cl_mean - vec_cl.at(i).z * vec_cl.at(i).z;

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

    for (int ii = 0; ii < 5; ii++) {
      if (elayer[ii] != 0.) {
        xv[ii] /= elayer[ii];
        yv[ii] /= elayer[ii];
        zv[ii] /= elayer[ii];
        tv[ii] /= elayer[ii];

        if (tinner < 0) {
          tinner = tv[ii];
          zinner = zv[ii];
        }

        touter = tv[ii];
        zouter = zv[ii];

        sx += xv[ii];
        sy += yv[ii];
        sz += zv[ii];
        sxz += xv[ii] * zv[ii];
        syz += yv[ii] * zv[ii];
        sx2 += xv[ii] * xv[ii];
        sy2 += yv[ii] * yv[ii];
        sz2 += zv[ii] * zv[ii];

        nlayer++;
      }
    }

    if (tinner != touter)
      dir = (zouter - zinner) * (touter - tinner) /
            TMath::Abs((touter - tinner) * (zouter - zinner));
    else
      dir = 0.;

    double ax = (nlayer * sxz - sz * sx) / (nlayer * sz2 - sz * sz);
    double bx = (sx * sz2 - sz * sxz) / (nlayer * sz2 - sz * sz);
    double ay = (nlayer * syz - sz * sy) / (nlayer * sz2 - sz * sz);
    double by = (sy * sz2 - sz * syz) / (nlayer * sz2 - sz * sz);

    double mod = TMath::Sqrt(1 + ax * ax + ay * ay);

    vec_cl.at(i).sx = dir * ax / mod;
    vec_cl.at(i).sy = dir * ay / mod;
    vec_cl.at(i).sz = dir * 1. / mod;
  }
}

bool value_comparer(std::map<int, int>::value_type& i1,
                    std::map<int, int>::value_type& i2)
{
  return i1.second < i2.second;
}

void PidBasedClustering(TG4Event* ev, std::vector<dg_cell>* vec_cell,
                        std::vector<cluster>& vec_cl)
{
  const double cell_max_dt = 30.;  // ns -> dt > 30. ns is unphysical

  std::vector<int> pid(vec_cell->size());
  std::map<int, int> hit_pid;

  for (unsigned int i = 0; i < vec_cell->size(); i++) {
    // find particle corresponding to more p.e.
    hit_pid.clear();

    if (vec_cell->at(i).ps1.size() > 0)
      for (unsigned int j = 0; j < vec_cell->at(i).ps1.at(0).photo_el.size();
           j++) {
        hit_pid[ev->SegmentDetectors["EMCalSci"]
                    .at(vec_cell->at(i).ps1.at(0).photo_el.at(j).h_index)
                    .PrimaryId]++;
      }

    if (vec_cell->at(i).ps2.size() > 0)
      for (unsigned int j = 0; j < vec_cell->at(i).ps2.at(0).photo_el.size();
           j++) {
        hit_pid[ev->SegmentDetectors["EMCalSci"]
                    .at(vec_cell->at(i).ps2.at(0).photo_el.at(j).h_index)
                    .PrimaryId]++;
      }

    pid[i] =
        std::max_element(hit_pid.begin(), hit_pid.end(), value_comparer)->first;
  }

  std::vector<int> unique_pid = pid;
  std::sort(unique_pid.begin(), unique_pid.end());
  std::vector<int>::iterator last =
      std::unique(unique_pid.begin(), unique_pid.end());
  unique_pid.erase(last, unique_pid.end());

  for (unsigned int i = 0; i < unique_pid.size(); i++) {
    cluster cl;
    cl.tid = unique_pid[i];

    for (unsigned int j = 0; j < pid.size(); j++) {
      if (pid[j] == unique_pid[i]) {
        // good cell should have signal on both side and a tdc different less
        // than 30 ns (5.85 ns/m * 4 m)
        if (vec_cell->at(j).ps1.size() == 0 ||
            vec_cell->at(j).ps2.size() == 0 ||
            std::abs(vec_cell->at(j).ps1.at(0).tdc -
                     vec_cell->at(j).ps2.at(0).tdc) > cell_max_dt)
          continue;

        cl.cells.push_back(vec_cell->at(j));
      }
    }
    if (cl.cells.size() != 0) vec_cl.push_back(cl);
  }
}

void MeanAndRMS(std::vector<dg_tube>& digits, TH1D& hmeanX, TH1D& hrmsX,
                TH1I& hnX, TH1D& hmeanY, TH1D& hrmsY, TH1I& hnY)
{

  double mean, mean2, var, rms;
  int n;

  for (unsigned int j = 0; j < digits.size(); j++) {
    if (digits.at(j).hor == 1) {
      hrmsY.Fill(digits.at(j).z, digits.at(j).y * digits.at(j).y);
      hmeanY.Fill(digits.at(j).z, digits.at(j).y);
      hnY.Fill(digits.at(j).z);
    } else {
      hrmsX.Fill(digits.at(j).z, digits.at(j).x * digits.at(j).x);
      hmeanX.Fill(digits.at(j).z, digits.at(j).x);
      hnX.Fill(digits.at(j).z);
    }
  }

  for (int j = 0; j < hmeanX.GetNbinsX(); j++) {
    mean = -1;
    rms = -1;
    n = hnX.GetBinContent(j + 1);

    if (n != 0) {
      mean = hmeanX.GetBinContent(j + 1) / n;
      mean2 = hrmsX.GetBinContent(j + 1) / n;
      var = mean2 - mean * mean;
      rms = sqrt(var);
    }
    if (n > 1) n = 2;

    hrmsX.SetBinContent(j + 1, rms);
    hmeanX.SetBinContent(j + 1, mean);

    mean = -1;
    rms = -1;
    n = hnY.GetBinContent(j + 1);
    if (n != 0) {
      mean = hmeanY.GetBinContent(j + 1) / n;
      mean2 = hrmsY.GetBinContent(j + 1) / n;
      var = mean2 - mean * mean;
      rms = sqrt(var);
    }
    if (n > 1) n = 2;

    hrmsY.SetBinContent(j + 1, rms);
    hmeanY.SetBinContent(j + 1, mean);
  }
}

void filterDigitsModule(std::map<int, double>& yd, std::vector<int>& toremove,
                        const double epsilon)
{
  double mean = 0.;
  double var = 0.;
  double rms = 0.;
  double rms_thr = 0.;

  for (std::map<int, double>::iterator it = yd.begin(); it != yd.end(); ++it) {
    mean += it->second;
    var += it->second * it->second;
  }

  mean /= yd.size();
  var /= yd.size();
  var -= mean * mean;
  rms = sqrt(var) * (1. + epsilon);

  double maxd = 0.;
  std::map<int, double>::iterator idx;

  for (std::map<int, double>::iterator it = yd.begin(); it != yd.end(); ++it) {
    if (std::abs(it->second - mean) > maxd) {
      maxd = std::abs(it->second - mean);
      idx = it;
    }
  }

  if (maxd > rms && rms > rms_thr) {
    toremove.push_back(idx->first);
    yd.erase(idx);
    filterDigitsModule(yd, toremove, epsilon);
  }
}

void filterDigits(std::vector<dg_tube>& digits, TH1D& hdummy,
                  const double epsilon)
{
  std::map<int, std::map<int, double> > dy;
  std::map<int, std::map<int, double> > dx;
  std::vector<int> toremove;
  for (unsigned int i = 0; i < digits.size(); i++) {
    if (digits.at(i).hor == 1) {
      dy[hdummy.FindBin(digits.at(i).z)][i] = digits.at(i).y;
    } else {
      dx[hdummy.FindBin(digits.at(i).z)][i] = digits.at(i).x;
    }
  }

  for (std::map<int, std::map<int, double> >::iterator it = dy.begin();
       it != dy.end(); ++it) {
    filterDigitsModule(it->second, toremove, epsilon);
  }

  for (std::map<int, std::map<int, double> >::iterator it = dx.begin();
       it != dx.end(); ++it) {
    filterDigitsModule(it->second, toremove, epsilon);
  }

  std::sort(toremove.begin(), toremove.end());

  for (int i = toremove.size() - 1; i >= 0; i--) {
    digits.erase(digits.begin() + toremove[i]);
  }
}

void vtxfinding(double& xvtx_reco, double& yvtx_reco, double& zvtx_reco,
                int& VtxType, TH1D& hmeanX, TH1D& hmeanY, TH1D& hrmsX,
                TH1D& hrmsY, TH1I& hnX, TH1I& hnY)
{
  int idx = 0;
  double rms = 10000.;
  double rmsX, rmsY;
  VtxType = -1;

  for (int j = 0; j < hrmsX.GetNbinsX(); j++) {
    if (hnX.GetBinContent(j + 1) >= 2 && hnY.GetBinContent(j + 1) >= 2) {
      rmsX = hrmsX.GetBinContent(j + 1);
      rmsY = hrmsY.GetBinContent(j + 1);

      if (rmsX * rmsX + rmsY * rmsY < rms) {
        rms = rmsX * rmsX + rmsY * rmsY;
        idx = j + 1;
        VtxType = 2;
      }
    }
  }

  if (VtxType == -1) {
    for (int j = 0; j < hrmsX.GetNbinsX(); j++) {
      if (hnX.GetBinContent(j + 1) > 0 && hnY.GetBinContent(j + 1) > 0) {
        idx = j + 1;
        VtxType = 1;
        break;
      }
    }
  }

  if (VtxType != -1) {
    xvtx_reco = hmeanX.GetBinContent(idx);
    yvtx_reco = hmeanY.GetBinContent(idx);
    zvtx_reco = hrmsX.GetXaxis()->GetBinLowEdge(idx);
  } else {
    xvtx_reco = -1E10;
    yvtx_reco = -1E10;
    zvtx_reco = -1E10;
  }
}

void VertexFind(double& xvtx_reco, double& yvtx_reco, double& zvtx_reco,
                int& VtxType, std::vector<dg_tube> digits,
                std::vector<double>& binning, double epsilon)
{
  TH1D hrmsX("hrmsX", "rmsX;Z (mm); rmsX (mm)", binning.size() - 1,
             binning.data());
  TH1D hrmsY("hrmsY", "rmsY;Z (mm); rmsY (mm)", binning.size() - 1,
             binning.data());

  TH1D hmeanX("hmeanX", "meanX;Z (mm); meanX (mm)", binning.size() - 1,
              binning.data());
  TH1D hmeanY("hmeanY", "meanY;Z (mm); meanY (mm)", binning.size() - 1,
              binning.data());

  TH1I hnX("hnX", "nX;Z (mm); nX", binning.size() - 1, binning.data());
  TH1I hnY("hnY", "nY;Z (mm); nY", binning.size() - 1, binning.data());

  TH1D hdummy("hdummy", "hdummy;Z (mm); multipliciy", binning.size() - 1,
              binning.data());

  // filter digits
  // digits distant more than rms from mean are removed
  filterDigits(digits, hdummy, epsilon);

  // eveluate rms, mean and multiplicity in Y and X as a function of the STT
  // module
  MeanAndRMS(digits, hmeanX, hrmsX, hnX, hmeanY, hrmsY, hnY);

  // find vertex as modules with less spreas between tracks
  vtxfinding(xvtx_reco, yvtx_reco, zvtx_reco, VtxType, hmeanX, hmeanY, hrmsX,
             hrmsY, hnX, hnY);
}

///////// to be reimplemented with SANDGeoManager
void DetermineModulesPosition(TGeoManager* g, std::vector<double>& binning)
{
  TString path_prefix(sand_geometry::path_internal_volume);
  TGeoVolume* v = g->FindVolumeFast(sand_geometry::name_internal_volume);

  double origin[3];
  double master[3];
  double last_z = 0.;
  double dz;
  bool is_first = true;

  binning.clear();

  origin[0] = 0.;
  origin[1] = 0.;
  origin[2] = 0.;

  // assuming they are order by Z position
  for (int i = 0; i < v->GetNdaughters(); i++) {
    TString name = v->GetNode(i)->GetName();

    if (name.Contains("TrMod") || name.Contains("C3H6Mod") ||
        name.Contains("CMod")) {
      TString path = path_prefix + name;
      g->cd(path.Data());
      g->LocalToMaster(origin, master);
      TGeoBBox* b = (TGeoBBox*)v->GetNode(i)->GetVolume()->GetShape();
      dz = b->GetDX();

      if (is_first) {
        is_first = false;
        binning.push_back(master[2] - dz);
      } else {
        binning.push_back(0.5 * (last_z + master[2] - dz));
      }
      last_z = master[2] + dz;
    }
  }
  binning.push_back(last_z);
}

enum class STT_Mode { fast_only_primaries, fast, full };
enum class ECAL_Mode { fast };

void Reconstruct(std::string const& fname_hits, std::string const& fname_digits,
                 std::string const& fname_out, STT_Mode stt_mode,
                 ECAL_Mode ecal_mode)
{
  std::cout << "Reconstruct\ninput hits: " << fname_hits
            << "\ninput digits: " << fname_digits
            << "\noutput (update): " << fname_out << '\n';

  TFile f_hits(fname_hits.data(), "READ");
  TFile f_digits(fname_digits.data(), "READ");
  TFile f_out(fname_out.data(), "UPDATE");

  if (f_hits.IsZombie() || f_digits.IsZombie() || f_out.IsZombie()) {
    std::cout << "Error in opening file\n";
    exit(1);
  }

  TTree* tTrueMC = (TTree*)f_hits.Get("EDepSimEvents");
  TGeoManager* geo = (TGeoManager*)f_hits.Get("EDepSimGeometry");
  TTree* tDigit = (TTree*)f_digits.Get("tDigit");

  if (tTrueMC == nullptr || geo == nullptr || tDigit == nullptr) {
    std::cout << "Error in retrieving objects from root file: "
              << (tTrueMC == nullptr ? "EDepSimEvents " : "")
              << (geo == nullptr ? "EDepSimGeometry " : "")
              << (tDigit == nullptr ? "tDigit " : "") << '\n';
    exit(-1);
  }

  std::vector<double> sampling;

  DetermineModulesPosition(geo, sampling);

  tDigit->AddFriend(tTrueMC);

  TTree* t = tDigit;

  TG4Event* ev = new TG4Event;
  t->SetBranchAddress("Event", &ev);

  std::vector<dg_tube>* vec_digi = new std::vector<dg_tube>;
  std::vector<dg_cell>* vec_cell = new std::vector<dg_cell>;

  t->SetBranchAddress("dg_tube", &vec_digi);
  t->SetBranchAddress("dg_cell", &vec_cell);

  std::vector<track> vec_tr;
  std::vector<cluster> vec_cl;

  TTree tout("tReco", "tReco");
  tout.Branch("track", "std::vector<track>", &vec_tr);
  tout.Branch("cluster", "std::vector<cluster>", &vec_cl);

  const int nev = t->GetEntries();
  const double epsilon = 0.5;
  const double tol_phi = 0.1;
  const double tol_x = 100.;
  const int tol_mod = 4;
  const int mindigtr = 3;
  const double dn_tol = 1.E7;
  const double dz_tol = 1.E7;

  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;

  for (int i = 0; i < nev; i++) {
    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i) / nev * 100)
              << "%]" << std::flush;

    t->GetEntry(i);

    vec_tr.clear();
    vec_cl.clear();

    double xvtx_reco, yvtx_reco, zvtx_reco;
    int VtxType;

    std::vector<dg_tube> clustersY;
    std::vector<dg_tube> clustersX;

    switch (stt_mode) {
      case STT_Mode::fast_only_primaries:
        TrackFind(ev, vec_digi, vec_tr, TrackFilter::only_primaries);
        TrackFit(vec_tr);
        break;
      case STT_Mode::fast:
        TrackFind(ev, vec_digi, vec_tr);
        TrackFit(vec_tr);
        break;
      case STT_Mode::full:
        VertexFind(xvtx_reco, yvtx_reco, zvtx_reco, VtxType, *vec_digi,
                   sampling, epsilon);
        TrackFind(vec_tr, *vec_digi, sampling, xvtx_reco, yvtx_reco, zvtx_reco,
                  tol_phi, tol_x, tol_mod, mindigtr, dn_tol, dz_tol);
        TrackFit(vec_tr, sampling, xvtx_reco, yvtx_reco, zvtx_reco);
        break;
    }

    switch (ecal_mode) {
      case ECAL_Mode::fast:
        // PreCluster(vec_cell, vec_cl);
        // Filter(vec_cl);
        PidBasedClustering(ev, vec_cell, vec_cl);
        Merge(vec_cl);
        break;
    }
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

  f_out.cd();
  tout.Write("", TObject::kOverwrite);
  f_out.Close();
}

void help_reco()
{
  std::cout
      << "usage: Reconstruct hit_file digit_file output_file [stt_mode]\n";
  std::cout << "    - stt_mode: 'stt_mode::fast_only_primaries' (default) \n";
  std::cout << "                'stt_mode::fast' \n";
  std::cout << "                'stt_mode::full' \n";
}

int main(int argc, char* argv[])
{
  // boost::program_options wuold be great here....

  if (argc < 4 || argc > 6) {
    help_reco();
    return -1;
  }

  auto stt_mode = STT_Mode::fast_only_primaries;
  if (argc > 4 && strcmp(argv[4], "stt_mode::full") == 0) {
    stt_mode = STT_Mode::full;
    std::cout << "STT_Mode: full\n";
  } else if (argc > 4 && strcmp(argv[4], "stt_mode::fast") == 0) {
    stt_mode = STT_Mode::fast;
    std::cout << "STT_Mode: fast\n";
  } else {
    std::cout << "STT_Mode: fast_only_primaries\n";
  }

  Reconstruct(argv[1], argv[2], argv[3], stt_mode, ECAL_Mode::fast);
  return 0;
}
