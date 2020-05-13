#include <TGeoManager.h>
#include <TGeoBBox.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TParameter.h>

#include <iostream>
#include <vector>
#include <map>

#include "../include/track3D.h"

#include "/wd/dune-it/ext_bkg/kloe-simu/src/display.cpp"

bool isCluBigger(const std::vector<digit>& v1, const std::vector<digit>& v2)
{
  return v1.size() > v2.size();
}

bool isDigUpstream(const digit& d1, const digit& d2)
{
  return d1.z < d2.z;
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

void fillLayers(std::map<int, std::vector<digit> >& m, std::vector<digit>* d,
                TH1D& hdummy, int hor)
{
  for (unsigned int k = 0; k < d->size(); k++) {
    if (d->at(k).hor == hor)
      m[-1 * hdummy.FindBin(d->at(k).z)].push_back(d->at(k));
  }
}

void evalUV(double& u, double& v, double zv, double yv, double z, double y)
{
  u = (z - zv);
  v = (yv - y);
  double d = (u * u + v * v);

  u /= d;
  v /= d;
}

void evalPhi(double& phi, double u, double v)
{
  phi = TMath::ATan2(v, u);
}

void findNearDigPhi(int& idx, double exp_phi, std::vector<digit>& vd, double zv,
                    double yv)
{
  double dphi = 1E3;

  double phi, u, v;

  for (unsigned int k = 0; k < vd.size(); k++) {
    evalUV(u, v, zv, yv, vd.at(k).z, vd.at(k).y);
    evalPhi(phi, u, v);
    if (abs(exp_phi - phi) < dphi) {
      dphi = abs(exp_phi - phi);
      idx = k;
    }
  }
}

bool findNearDigX(int& idx, double exp_x, std::vector<digit>& vd, double xvtx)
{
  double dx = 10000.;

  double x;

  idx = -1;

  bool found = false;

  const double xtol_on_vtx = 50.;

  for (unsigned int k = 0; k < vd.size(); k++) {
    x = vd.at(k).x;
    if (abs(exp_x - x) < dx &&
        (abs(x - xvtx) < abs(exp_x - xvtx) || abs(x - xvtx) < xtol_on_vtx)) {
      dx = abs(exp_x - x);
      idx = k;
      found = true;
    }
  }

  return found;
}

void MeanAndRMS(std::vector<digit>* digits, TH1D& hmeanX, TH1D& hrmsX,
                TH1I& hnX, TH1D& hmeanY, TH1D& hrmsY, TH1I& hnY)
{

  double mean, mean2, var, rms;
  int n;

  for (unsigned int j = 0; j < digits->size(); j++) {
    if (digits->at(j).hor == 1) {
      hrmsY.Fill(digits->at(j).z, digits->at(j).y * digits->at(j).y);
      hmeanY.Fill(digits->at(j).z, digits->at(j).y);
      hnY.Fill(digits->at(j).z);
    } else {
      hrmsX.Fill(digits->at(j).z, digits->at(j).x * digits->at(j).x);
      hmeanX.Fill(digits->at(j).z, digits->at(j).x);
      hnX.Fill(digits->at(j).z);
    }
  }

  for (unsigned int j = 0; j < hmeanX.GetNbinsX(); j++) {
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
    if (abs(it->second - mean) > maxd) {
      maxd = abs(it->second - mean);
      idx = it;
    }
  }

  if (maxd > rms && rms > rms_thr) {
    toremove.push_back(idx->first);
    yd.erase(idx);
    filterDigitsModule(yd, toremove, epsilon);
  }
}

void filterDigits(std::vector<digit>* digits, TH1D& hdummy,
                  const double epsilon)
{
  std::map<int, std::map<int, double> > dy;
  std::map<int, std::map<int, double> > dx;
  std::vector<int> toremove;
  for (unsigned int i = 0; i < digits->size(); i++) {
    if (digits->at(i).hor == 1) {
      dy[hdummy.FindBin(digits->at(i).z)][i] = digits->at(i).y;
    } else {
      dx[hdummy.FindBin(digits->at(i).z)][i] = digits->at(i).x;
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
    digits->erase(digits->begin() + toremove[i]);
  }
}

void vtxFinding(double& xvtx_reco, double& yvtx_reco, double& zvtx_reco,
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

void findTracksY(std::map<int, std::vector<digit> >& mdY,
                 std::vector<std::vector<digit> >& clustersY, TH1D& hdummy,
                 double xvtx_reco, double yvtx_reco, double zvtx_reco,
                 double phi_tol = 0.1, double dlay_tol = 5)
{
  int idx = 0;
  int prev_mod;
  double prev_phi, exp_phi;
  double prev_z, slope;

  double u, v, phi, dphi;

  digit current_digit;

  // loop on modules
  while (mdY.size() != 0) {
    // get most downstream module
    std::map<int, std::vector<digit> >::iterator ite = mdY.begin();
    std::vector<digit>* layer = &(ite->second);

    // loop on digits in the module
    while (layer->size() != 0) {
      // get forst digit
      current_digit = layer->front();

      // create cluster and insert first digit
      std::vector<digit> clY;
      clY.push_back(std::move(current_digit));

      // remove the current digit from the module
      layer->erase(layer->begin());

      // get module id
      prev_mod = ite->first;

      // get iterator to the next module
      std::map<int, std::vector<digit> >::iterator nite = std::next(ite);

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
          std::vector<digit>* nlayer = &(nite->second);

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
            if (abs(dphi) <= phi_tol) {
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
                std::map<int, std::vector<digit> >::iterator dummy =
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

void findTracksX(std::map<int, std::vector<digit> >& mdX,
                 std::vector<std::vector<digit> >& clustersX, TH1D& hdummy,
                 double xvtx_reco, double yvtx_reco, double zvtx_reco,
                 double x_tol = 100, double dlay_tol = 5)
{
  int idx = 0;
  int prev_mod;

  double prev_x, dx;

  digit current_digit;

  // loop on modules from downstream to upstream
  while (mdX.size() != 0) {
    // get most downstream module
    std::map<int, std::vector<digit> >::iterator ite = mdX.begin();
    std::vector<digit>* layer = &(ite->second);

    // loop on digits in the module
    while (layer->size() != 0) {
      // get forst digit
      current_digit = layer->front();

      // create cluster and insert first digit
      std::vector<digit> clX;
      clX.push_back(std::move(current_digit));

      // remove the current digit from the module
      layer->erase(layer->begin());

      // get module id
      prev_mod = ite->first;

      // get iterator to the next module
      std::map<int, std::vector<digit> >::iterator nite = std::next(ite);

      // reset parameters
      prev_x = current_digit.x;
      prev_mod = ite->first;

      // loop on upstream modules
      while (nite != mdX.end()) {
        // check if module are downstream the reco vtx
        if (hdummy.GetXaxis()->GetBinUpEdge(-1 * nite->first) > zvtx_reco) {
          // get next layer
          std::vector<digit>* nlayer = &(nite->second);

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
              if (abs(dx) <= x_tol) {
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
                  std::map<int, std::vector<digit> >::iterator dummy =
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

void mergeXYTracks(std::vector<std::vector<digit> >& clustersX,
                   std::vector<std::vector<digit> >& clustersY,
                   std::vector<track3D>& tr3D, double dn_tol, double dz_tol,
                   std::map<int, int>& trMCIDY, std::map<int, int>& trMCIDX,
                   int& ntr3DOK, int& ntr3DNO, TG4Event* ev)
{
  double dzend;
  int index = 0;

  // loop on Y clusters
  for (int jj = 0; jj < clustersY.size(); jj++) {
    // resent downstream z residual
    dzend = 10000.;

    // loop on X clusters
    for (int kk = 0; kk < clustersX.size(); kk++) {
      // evaluate the percentege difference in number of digits: abs(dn)/mn
      double dn = int(clustersY.at(jj).size()) - int(clustersX.at(kk).size());
      double mn = 0.5 * (clustersY.at(jj).size() + clustersX.at(kk).size());

      // search cluster X with abs(dn)/mn less than tolerance and
      // residual in the most downstream digit z less than tolerance
      if (abs(dn) / mn < dn_tol &&
          abs(clustersY.at(jj).back().z - clustersX.at(kk).back().z) < dz_tol &&
          abs(clustersY.at(jj).front().z - clustersX.at(kk).front().z) <
              dz_tol) {
        // find cluster X best matching Y cluster in term of
        // residual in the most downstream and upstream digit z
        if (abs(clustersY.at(jj).back().z - clustersX.at(kk).back().z) +
                abs(clustersY.at(jj).front().z - clustersX.at(kk).front().z) <
            dzend) {
          dzend = abs(clustersY.at(jj).back().z - clustersX.at(kk).back().z) +
                  abs(clustersY.at(jj).front().z - clustersX.at(kk).front().z);
          track3D tr;
          tr.tid = index++;
          tr.clX = std::move(clustersX.at(kk));
          tr.clY = std::move(clustersY.at(jj));
          clustersY.erase(clustersY.begin() + jj--);
          clustersX.erase(clustersX.begin() + kk--);

          if (trMCIDY[jj] == trMCIDX[kk]) {
            ntr3DOK++;

            TG4Trajectory pMC = ev->Trajectories.at(trMCIDY[jj]);
            if (pMC.GetTrackId() != trMCIDY[jj])
              std::cout << "Error: track MC does not macth" << std::endl;
            tr.px = pMC.GetInitialMomentum().X();
            tr.py = pMC.GetInitialMomentum().Y();
            tr.pz = pMC.GetInitialMomentum().Z();
          } else {
            ntr3DNO++;
            tr.px = -999;
            tr.py = -999;
            tr.pz = -999;
          }

          tr3D.push_back(std::move(tr));
          break;
        }
      }
    }
  }
}

void fitCircle(std::vector<track3D>& tr3D, TH1D& hdummy)
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
    if (prodVect != 0.) prodVect /= abs(prodVect);
    tr3D.at(nn).h = int(prodVect);

    if (yres > 0)
      tr3D.at(nn).ysig = +1;
    else
      tr3D.at(nn).ysig = -1;
  }
}

void fitLine(std::vector<track3D>& tr3D, double xvtx_reco, double yvtx_reco,
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
      if (abs(tr3D.at(i).clX.at(kk).z - tr3D.at(i).zc) > tr3D.at(i).r)
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

void processEvents(const double tol_phi, const double tol_x,
                   const double tol_mod, const int mindigtr,
                   const double frac = 1.0, double dn_tol = 0.7,
                   double dz_tol = 200.)
{
  gROOT->SetBatch();

  const int minreg = 1;
  const double epsilon = 0.5;
  const bool wideMultReg = false;

  bool display = false;
  bool save2root = true && display;
  bool save2pdf = true && display;

  // root [0] TGeoManager::Import("../geo/nd_hall_kloe_empty.gdml")
  // (TGeoManager *) 0x3037200
  // root [12]
  // gGeoManager->cd("volWorld/rockBox_lv_0/volDetEnclosure_0/volKLOE_0")
  // (bool) true
  // root [16] double vorigin[3]
  // (double [3]) { 0.0000000, 0.0000000, 0.0000000 }
  // root [17] double vmaster[3]
  // (double [3]) { 0.0000000, 0.0000000, 0.0000000 }
  // root [18] gGeoManager->LocalToMaster(vorigin,vmaster)
  // root [19] vmaster
  // (double [3]) { 0.0000000, -238.47300, 2391.0000 }
  // root [7] v = gGeoManager->GetVolume("volSTTFULL")
  // (TGeoVolume *) 0x17fc5b10
  // root [9] TGeoTube* tb = (TGeoTube*) v->GetShape()
  // root [11] tb->GetRmin()
  // (double) 0.0000000
  // root [12] tb->GetRmax()
  // (double) 200.00000
  // root [13] tb->GetDz()
  // (double) 169.00000

  double kloe_center[] = {0.0000000, -2384.7300, 23910.000};
  double kloe_size[] = {2 * 1690.00000, 2 * 2000.00000, 2 * 2000.00000};

  TFile f("../files/reco/numu_internal_10k.0.reco.root");

  TTree* tDigit = (TTree*)f.Get("tDigit");
  TTree* tMC = (TTree*)f.Get("EDepSimEvents");
  TGeoManager* g = (TGeoManager*)f.Get("EDepSimGeometry");

  TString path_prefix =
      "volWorld_PV/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/"
      "MagIntVol_volume_PV_0/volSTTFULL_PV_0/";
  TGeoVolume* v = g->FindVolumeFast("volSTTFULL_PV");

  double origin[3];
  double master[3];
  double last_z = 0.;
  double dz;
  std::vector<double> binning;
  bool is_first = true;

  origin[0] = 0.;
  origin[1] = 0.;
  origin[2] = 0.;

  // assuming they are order by Z position
  for (int i = 0; i < v->GetNdaughters(); i++) {
    TString name = v->GetNode(i)->GetName();

    if (name.Contains("sttmod") || name.Contains("volfrontST")) {
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

  const double z_sampl = 44.4;
  const double xy_sampl = 5;

  TDatabasePDG* pdg = new TDatabasePDG();
  int nPriPart = 0;
  int nPriPartXFound = 0;
  int nPriPartYFound = 0;
  int nPriPartRecoOK = 0;

  int nall = 0;
  int nok = 0;

  double avePurityY = 0;
  double avePurityX = 0;
  double aveNY = 0;
  double aveNX = 0;
  int countClY = 0;
  int countClX = 0;

  double vrmsX = 0.;
  double vrmsY = 0.;
  double vrmsZ = 0.;
  double vrms3D = 0.;
  double vmeanX = 0.;
  double vmeanY = 0.;
  double vmeanZ = 0.;
  double vmean3D = 0.;
  double norm_dist = 0.;

  std::cout << std::setw(20) << "Z sampling (mm)" << std::setw(20)
            << "XY sampling (mm)" << std::endl;
  std::cout << std::setw(20) << z_sampl << std::setw(20) << xy_sampl
            << std::endl;

  /*
  for(unsigned int i = 0; i < binning.size()-1; i++)
  {
    std::cout << i << " " << binning.data()[i] << " " << binning.data()[i+1] <<
  " " << binning.data()[i+1] - binning.data()[i] << std::endl;
  }
  */

  gSystem->Load("lib/libTrack3D.so");

  std::vector<digit>* digits = new std::vector<digit>;
  TG4Event* ev = new TG4Event;
  std::vector<std::vector<digit> > clustersY;
  std::vector<std::vector<digit> > clustersX;

  std::map<int, int> trMCIDY;
  std::map<int, int> trMCIDX;

  int ntr3DOK = 0;
  int ntr3DNO = 0;
  int n100MeV = 0;
  int n500MeV = 0;
  int n1GeV = 0;

  tDigit->SetBranchAddress("Stt", &digits);
  tMC->SetBranchAddress("Event", &ev);

  TFile fout("vtx.root", "RECREATE");
  TTree tv("tv", "tv");

  double xvtx_true, yvtx_true, zvtx_true;
  double xvtx_reco, yvtx_reco, zvtx_reco;

  int isCC;
  int VtxMulti;
  int VtxMono;
  int VtxType;

  std::vector<track3D> tracks;

  tv.Branch("xvtx_true", &xvtx_true, "xvtx_true/D");
  tv.Branch("yvtx_true", &yvtx_true, "yvtx_true/D");
  tv.Branch("zvtx_true", &zvtx_true, "zvtx_true/D");

  tv.Branch("xvtx_reco", &xvtx_reco, "xvtx_reco/D");
  tv.Branch("yvtx_reco", &yvtx_reco, "yvtx_reco/D");
  tv.Branch("zvtx_reco", &zvtx_reco, "zvtx_reco/D");

  tv.Branch("VtxMulti", &VtxMulti, "VtxMulti/I");
  tv.Branch("VtxMono", &VtxMono, "VtxMono/I");
  tv.Branch("isCC", &isCC, "isCC/I");

  tv.Branch("tracks", "std::vector<track3D>", &tracks);

  TH1D hrmsX("hrmsX", "rmsX;Z (mm); rmsX (mm)", binning.size() - 1,
             binning.data());
  TH1D hrmsY("hrmsY", "rmsY;Z (mm); rmsY (mm)", binning.size() - 1,
             binning.data());
  hrmsY.SetLineColor(kRed);

  TH1D hmeanX("hmeanX", "meanX;Z (mm); meanX (mm)", binning.size() - 1,
              binning.data());
  TH1D hmeanY("hmeanY", "meanY;Z (mm); meanY (mm)", binning.size() - 1,
              binning.data());
  hmeanY.SetLineColor(kRed);

  TH1I hnX("hnX", "nX;Z (mm); nX", binning.size() - 1, binning.data());
  TH1I hnY("hnY", "nY;Z (mm); nY", binning.size() - 1, binning.data());
  hnY.SetLineColor(kRed);

  TH1D hdummy("hdummy", "hdummy;Z (mm); multipliciy", binning.size() - 1,
              binning.data());

  hrmsX.SetStats(false);
  hrmsY.SetStats(false);
  hmeanX.SetStats(false);
  hmeanY.SetStats(false);
  hnX.SetStats(false);
  hnY.SetStats(false);

  TParameter<double> xv("xv", 0);
  TParameter<double> yv("yv", 0);
  TParameter<double> zv("zv", 0);

  init("../files/reco/numu_internal_10k.0.reco.root");

  int first = 0;
  int last = 9999; /*tDigit->GetEntries()-1;*/
  int nev = last - first + 1;

  TCanvas c("c", "revo event", 1200, 600);
  c.SaveAs("rms.pdf(");

  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;

  for (int i = first; i < last + 1; i++) {
    tracks.clear();

    tDigit->GetEntry(i);
    tMC->GetEntry(i);

    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i) / nev * 100)
              << "%]" << std::flush;

    TString reaction = ev->Primaries.at(0).GetReaction();
    VtxMulti = 0;
    VtxMono = 0;
    VtxType = 0;

    if (reaction.Contains("CC") == false)
      isCC = 0;
    else
      isCC = 1;

    xv.SetVal(ev->Primaries.at(0).GetPosition().X());
    yv.SetVal(ev->Primaries.at(0).GetPosition().Y());
    zv.SetVal(ev->Primaries.at(0).GetPosition().Z());

    xvtx_true = xv.GetVal();
    yvtx_true = yv.GetVal();
    zvtx_true = zv.GetVal();

    TDirectoryFile* fd;

    if (save2root)
      fd = new TDirectoryFile(TString::Format("ev_%d", i).Data(),
                              TString::Format("ev_%d", i).Data(), "", &fout);

    hrmsX.Reset("ICESM");
    hrmsY.Reset("ICESM");
    hmeanX.Reset("ICESM");
    hmeanY.Reset("ICESM");
    hnX.Reset("ICESM");
    hnY.Reset("ICESM");

    hrmsX.SetTitle(TString::Format("event: %d", i).Data());
    hrmsY.SetTitle(TString::Format("event: %d", i).Data());
    hmeanX.SetTitle(TString::Format("event: %d", i).Data());
    hmeanY.SetTitle(TString::Format("event: %d", i).Data());
    hnX.SetTitle(TString::Format("event: %d", i).Data());
    hnY.SetTitle(TString::Format("event: %d", i).Data());

    // filter digits
    // digits distant more than rms from mean are removed
    filterDigits(digits, hdummy, epsilon);

    // eveluate rms, mean and multiplicity in Y and X as a function of the STT
    // module
    MeanAndRMS(digits, hmeanX, hrmsX, hnX, hmeanY, hrmsY, hnY);

    // find vertex as modules with less spreas between tracks
    vtxFinding(xvtx_reco, yvtx_reco, zvtx_reco, VtxType, hmeanX, hmeanY, hrmsX,
               hrmsY, hnX, hnY);

    if (VtxType == 2)
      VtxMulti = 1;
    else if (VtxType == 1)
      VtxMono = 1;

    tDigit->GetEntry(i);

    // track finding with clustering in arctg(v/u) VS z
    std::map<int, std::vector<digit> > mdY;
    fillLayers(mdY, digits, hdummy, 1);
    clustersY.clear();

    // find track on XZ view
    std::map<int, std::vector<digit> > mdX;
    fillLayers(mdX, digits, hdummy, 0);
    clustersX.clear();

    // tolerance in phi: 0.1 rad
    // tolerance in module: 5
    findTracksY(mdY, clustersY, hdummy, xvtx_reco, yvtx_reco, zvtx_reco,
                tol_phi, tol_mod);

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
    std::sort(clustersY.begin(), clustersY.end(), isCluBigger);

    // order tracks by length
    std::sort(clustersX.begin(), clustersX.end(), isCluBigger);

    // order digits: upstream first
    for (int jj = 0; jj < clustersX.size(); jj++) {
      std::sort(clustersX.at(jj).begin(), clustersX.at(jj).end(),
                isDigUpstream);
    }

    for (int jj = 0; jj < clustersY.size(); jj++) {
      std::sort(clustersY.at(jj).begin(), clustersY.at(jj).end(),
                isDigUpstream);
    }

    //////////////////////////////////////////////////

    tDigit->GetEntry(i);
    std::map<int, int> primaryDigX;
    std::map<int, int> primaryDigY;

    vector<TG4HitSegment> hits = ev->SegmentDetectors["Straw"];

    for (unsigned int ll = 0; ll < digits->size(); ll++) {
      if (digits->at(ll).hor == 1)
        primaryDigY[hits[digits->at(ll).hindex.front()].GetPrimaryId()]++;
      else
        primaryDigX[hits[digits->at(ll).hindex.front()].GetPrimaryId()]++;
    }

    std::vector<int> primaryidY;
    int nhitYOK, nhitXOK;

    trMCIDY.clear();
    trMCIDX.clear();

    for (unsigned int ii = 0; ii < ev->Primaries.at(0).Particles.size(); ii++) {
      TParticlePDG* partPDG =
          pdg->GetParticle(ev->Primaries.at(0).Particles.at(ii).GetPDGCode());
      if (partPDG != 0) {
        if (partPDG->Charge() != 0.) {
          if (ev->Primaries.at(0).Particles.at(ii).GetMomentum().Z() > 0) {
            if (primaryDigY[ev->Primaries.at(0).Particles.at(ii).GetTrackId()] >
                    3 &&
                primaryDigX[ev->Primaries.at(0).Particles.at(ii).GetTrackId()] >
                    3) {
              primaryidY.push_back(
                  ev->Primaries.at(0).Particles.at(ii).GetTrackId());
              nPriPart++;

              double KE =
                  ev->Primaries.at(0).Particles.at(ii).GetMomentum().E() -
                  1E3 * partPDG->Mass();

              if (KE > 100) n100MeV++;
              if (KE > 500) n500MeV++;
              if (KE > 1000) n1GeV++;
            }
          }
        }
      } /*
       else
       {
         std::cout <<  ev->Primaries.at(0).Particles.at(ii).GetPDGCode() <<
       std::endl;
       }*/
    }

    std::vector<int> primaryidX = primaryidY;

    for (unsigned int ii = 0; ii < clustersY.size(); ii++) {
      std::map<int, int> pr;
      int tot = 0;

      for (unsigned int jj = 0; jj < clustersY.at(ii).size(); jj++) {
        for (unsigned int kk = 0; kk < clustersY.at(ii).at(jj).hindex.size();
             kk++) {
          pr[hits[clustersY.at(ii).at(jj).hindex.at(kk)].GetPrimaryId()]++;
          tot++;
        }
      }

      nhitYOK = 0;

      for (std::map<int, int>::iterator it = pr.begin(); it != pr.end(); it++) {
        if (it->second > nhitYOK) {
          trMCIDY[ii] = it->first;
          nhitYOK = it->second;
        }
      }

      for (std::map<int, int>::iterator it = pr.begin(); it != pr.end(); it++) {
        if (it->second >= frac * tot) {
          std::vector<int>::iterator iter =
              std::find(primaryidY.begin(), primaryidY.end(), it->first);
          if (iter != primaryidY.end()) {
            primaryidY.erase(iter);
            nPriPartYFound++;
            break;
          }
        }
      }

      avePurityY += double(nhitYOK) / tot;
      aveNY += double(nhitYOK);
      countClY++;
    }

    //////////////////////////////////////////////////

    //////////////////////////////////////////////////

    for (unsigned int ii = 0; ii < clustersX.size(); ii++) {
      std::map<int, int> pr;
      int tot = 0;

      for (unsigned int jj = 0; jj < clustersX.at(ii).size(); jj++) {
        for (unsigned int kk = 0; kk < clustersX.at(ii).at(jj).hindex.size();
             kk++) {
          pr[hits[clustersX.at(ii).at(jj).hindex.at(kk)].GetPrimaryId()]++;
          tot++;
        }
      }

      nhitXOK = 0;

      for (std::map<int, int>::iterator it = pr.begin(); it != pr.end(); it++) {
        if (it->second > nhitXOK) {
          trMCIDX[ii] = it->first;
          nhitXOK = it->second;
        }
      }

      for (std::map<int, int>::iterator it = pr.begin(); it != pr.end(); it++) {
        if (it->second >= frac * tot) {
          std::vector<int>::iterator iter =
              std::find(primaryidX.begin(), primaryidX.end(), it->first);
          if (iter != primaryidX.end()) {
            primaryidX.erase(iter);
            nPriPartXFound++;
            break;
          }
        }
      }

      avePurityX += double(nhitXOK) / tot;
      aveNX += double(nhitXOK);
      countClX++;
    }
    //////////////////////////////////////////////////

    // merge XZ and YZ clusters
    // tolerance in dn/n: 0.7
    // tolerance dz: 20 cm
    mergeXYTracks(clustersX, clustersY, tracks, dn_tol, dz_tol, trMCIDY,
                  trMCIDX, ntr3DOK, ntr3DNO, ev);

    // circular fit
    fitCircle(tracks, hdummy);

    // linear fit
    fitLine(tracks, xvtx_reco, yvtx_reco, zvtx_reco);

    const double ptol = 0.1;
    const int tol = 3;

    double phi_true, pt_true;
    double phi_reco, pt_reco;
    double px_true, py_true, pz_true;

    for (unsigned int kk = 0; kk < ev->Primaries.at(0).Particles.size(); kk++) {
      px_true = ev->Primaries.at(0).Particles.at(kk).GetMomentum().X();
      py_true = ev->Primaries.at(0).Particles.at(kk).GetMomentum().Y();
      pz_true = ev->Primaries.at(0).Particles.at(kk).GetMomentum().Z();

      pt_true = TMath::Sqrt(py_true * py_true + pz_true * pz_true);

      phi_true = TMath::ATan2(px_true, pt_true);

      for (unsigned int jj = 0; jj < tracks.size(); jj++) {
        phi_reco = TMath::ATan(tracks.at(jj).b);
        pt_reco = 0.299792458 * tracks.at(jj).r * 0.6;

        if (abs(pt_reco - pt_true) / pt_true < 0.2 &&
            abs(phi_reco - phi_true) < 0.02) {
          nPriPartRecoOK++;
          break;
        }
      }
    }

    if (VtxMulti == 1) {
      nall++;
      if (abs(xvtx_reco - xvtx_true) < tol * xy_sampl &&
          abs(yvtx_reco - yvtx_true) < tol * xy_sampl &&
          abs(zvtx_reco - zvtx_true) < tol * z_sampl)
        nok++;

      norm_dist = sqrt(pow((xvtx_reco - xvtx_true) / xy_sampl, 2) +
                       pow((yvtx_reco - yvtx_true) / xy_sampl, 2) +
                       pow((zvtx_reco - zvtx_true) / z_sampl, 2));

      vmeanX += xvtx_reco - xvtx_true;
      vmeanY += yvtx_reco - yvtx_true;
      vmeanZ += zvtx_reco - zvtx_true;
      vmean3D += norm_dist;

      vrmsX += (xvtx_reco - xvtx_true) * (xvtx_reco - xvtx_true);
      vrmsY += (yvtx_reco - yvtx_true) * (yvtx_reco - yvtx_true);
      vrmsZ += (zvtx_reco - zvtx_true) * (zvtx_reco - zvtx_true);
      vrms3D += norm_dist * norm_dist;
    }

    double minz, maxz;

    if (display) {

      // display of the event
      show(i, true, false, true);

      TCanvas* cev = (TCanvas*)gROOT->FindObject("cev");

      TGraph* gr;
      TMarker* m;
      TEllipse* el;
      TF1* ffx;
      TF1* ffy;
      double yexp;
      double phi0;

      // save to pdf
      if (save2pdf) {
        if (cev) {
          cev->SaveAs("rms.pdf");

          c.Clear();
          c.Divide(2, 1);

          c.cd(1)->DrawFrame(21500, -5000, 26500, 0,
                             TString::Format("Event: %d (ZY)", i).Data());
          c.cd(2)->DrawFrame(21500, -2500, 26500, 2500,
                             TString::Format("Event: %d (ZX)", i).Data());

          std::vector<double> dXx;
          std::vector<double> dXz;
          std::vector<double> dYy;
          std::vector<double> dYz;

          TGraph* grY;
          TGraph* grX;

          for (unsigned int j = 0; j < tracks.size(); j++) {
            dXx.clear();
            dXz.clear();
            dYy.clear();
            dYz.clear();

            c.cd(1);

            for (unsigned int kk = 0; kk < tracks.at(j).clY.size(); kk++) {
              dYy.push_back(tracks.at(j).clY.at(kk).y);
              dYz.push_back(tracks.at(j).clY.at(kk).z);
            }

            grY = new TGraph(dYz.size(), dYz.data(), dYy.data());
            grY->SetMarkerColor(j + 1);
            grY->SetMarkerStyle(7);
            grY->Draw("samep");

            minz = std::max(tracks.at(j).zc - tracks.at(j).r,
                            tracks.at(j).clY.front().z);
            maxz = std::min(tracks.at(j).zc + tracks.at(j).r,
                            tracks.at(j).clY.back().z);

            ffy = new TF1("", "[0]+[1]*TMath::Sqrt([2]*[2]-(x-[3])*(x-[3]))",
                          minz, maxz);
            ffy->SetParameter(0, tracks.at(j).yc);
            ffy->SetParameter(1, tracks.at(j).ysig);
            ffy->SetParameter(2, tracks.at(j).r);
            ffy->SetParameter(3, tracks.at(j).zc);
            ffy->SetLineColor(j + 1);
            ffy->Draw("same");

            c.cd(2);

            for (unsigned int kk = 0; kk < tracks.at(j).clX.size(); kk++) {
              dXx.push_back(tracks.at(j).clX.at(kk).x);
              dXz.push_back(tracks.at(j).clX.at(kk).z);
            }

            grX = new TGraph(dXz.size(), dXz.data(), dXx.data());
            grX->SetMarkerColor(j + 1);
            grX->SetMarkerStyle(7);
            grX->Draw("samep");

            double x0 = tracks.at(j).clX.front().x;
            double z0 = tracks.at(j).clX.front().z;

            double radq;
            if (abs(z0 - tracks.at(j).zc) > tracks.at(j).r)
              radq = 0;
            else
              radq =
                  TMath::Sqrt(tracks.at(j).r * tracks.at(j).r -
                              (z0 - tracks.at(j).zc) * (z0 - tracks.at(j).zc));

            yexp = tracks.at(j).yc + tracks.at(j).ysig * radq;

            phi0 = TMath::ATan2(yexp - tracks.at(j).yc, z0 - tracks.at(j).zc);

            minz = std::max(tracks.at(j).zc - tracks.at(j).r,
                            tracks.at(j).clX.front().z);
            maxz = std::min(tracks.at(j).zc + tracks.at(j).r,
                            tracks.at(j).clX.back().z);

            ffx = new TF1("",
                          "[0] - "
                          "[1]/"
                          "[2]*(TMath::ATan2(([4]+[7]*TMath::Sqrt([1]*[1]-(x-["
                          "5])*(x-[5]))) - [4],x - [5]) - [6])*[3]",
                          minz, maxz);

            ffx->SetParameter(0, x0);
            ffx->SetParameter(1, tracks.at(j).r);
            ffx->SetParameter(2, tracks.at(j).h);
            ffx->SetParameter(3, tracks.at(j).b);
            ffx->SetParameter(4, tracks.at(j).yc);
            ffx->SetParameter(5, tracks.at(j).zc);
            ffx->SetParameter(6, phi0);
            ffx->SetParameter(7, tracks.at(j).ysig);
            ffx->SetLineColor(j + 1);
            ffx->Draw("same");
          }
          c.SaveAs("rms.pdf");
        }
      }
    }

    if (save2root) {
      fd->Add(&hmeanX);
      fd->Add(&hnX);
      fd->Add(&hrmsX);

      fd->Add(&hmeanY);
      fd->Add(&hnY);
      fd->Add(&hrmsY);

      fd->Add(&xv);
      fd->Add(&yv);
      fd->Add(&zv);

      fout.cd();
      fd->Write();
    }
    tv.Fill();
  }

  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;

  vmeanX /= nall;
  vmeanY /= nall;
  vmeanZ /= nall;
  vmean3D /= nall;

  vrmsX /= nall;
  vrmsY /= nall;
  vrmsZ /= nall;
  vrms3D /= nall;

  vrmsX -= vmeanX * vmeanX;
  vrmsY -= vmeanY * vmeanY;
  vrmsZ -= vmeanZ * vmeanZ;
  vrms3D -= vmean3D * vmean3D;

  vrmsX = sqrt(vrmsX);
  vrmsY = sqrt(vrmsY);
  vrmsZ = sqrt(vrmsZ);
  vrms3D = sqrt(vrms3D);

  std::cout << std::setw(10) << "minreg" << std::setw(10) << "epsilon"
            << std::setw(15) << "wideMultReg" << std::setw(15) << "gefficiency"
            << std::setw(15) << "efficiency" << std::setw(10) << "rmsX"
            << std::setw(10) << "rmsY" << std::setw(10) << "rmsZ"
            << std::setw(10) << "rms3D" << std::endl;
  std::cout << std::setw(10) << minreg << std::setw(10) << epsilon
            << std::setw(15) << wideMultReg << std::setw(15)
            << double(nall) / nev << std::setw(15) << double(nok) / nall
            << std::setw(10) << vrmsX << std::setw(10) << vrmsY << std::setw(10)
            << vrmsZ << std::setw(10) << vrms3D << std::endl;

  avePurityY /= countClY;
  avePurityX /= countClX;
  aveNY /= countClY;
  aveNX /= countClX;

  std::cout << std::setw(10) << "fraction" << std::setw(10) << "tol_phi"
            << std::setw(10) << "tol_mod" << std::setw(10) << "mindigtr"
            << std::setw(10) << "dz_tol" << std::setw(10) << "dn_tol"
            << std::setw(10) << "total" << std::setw(10) << "fracY"
            << std::setw(10) << "avePY" << std::setw(10) << "aveNY"
            << std::setw(10) << "3DtrackOK" << std::setw(10) << "3DtrackNO"
            << std::endl;
  std::cout << std::setw(10) << frac << std::setw(10) << tol_phi
            << std::setw(10) << tol_mod << std::setw(10) << mindigtr
            << std::setw(10) << dz_tol << std::setw(10) << dn_tol
            << std::setw(10) << nPriPart << std::setw(10)
            << double(nPriPartYFound) / nPriPart << std::setw(10) << avePurityY
            << std::setw(10) << aveNY << std::setw(10) << ntr3DOK
            << std::setw(10) << ntr3DNO << std::endl;

  std::cout << std::setw(10) << "fraction" << std::setw(10) << "tol_x"
            << std::setw(10) << "tol_mod" << std::setw(10) << "mindigtr"
            << std::setw(10) << "dz_tol" << std::setw(10) << "dn_tol"
            << std::setw(10) << "total" << std::setw(10) << "fracX"
            << std::setw(10) << "avePX" << std::setw(10) << "aveNX"
            << std::setw(10) << "3DtrackOK" << std::setw(10) << "3DtrackNO"
            << std::endl;
  std::cout << std::setw(10) << frac << std::setw(10) << tol_x << std::setw(10)
            << tol_mod << std::setw(10) << mindigtr << std::setw(10) << dz_tol
            << std::setw(10) << dn_tol << std::setw(10) << nPriPart
            << std::setw(10) << double(nPriPartXFound) / nPriPart
            << std::setw(10) << avePurityX << std::setw(10) << aveNX
            << std::setw(10) << ntr3DOK << std::setw(10) << ntr3DNO
            << std::endl;

  std::cout << std::setw(10) << "total" << std::setw(10) << "n100MeV"
            << std::setw(10) << "n500MeV" << std::setw(10) << "n1GeV"
            << std::setw(10) << "3DtrOK" << std::setw(10) << "3DtrNO"
            << std::setw(10) << "3DtrGOOD" << std::endl;
  std::cout << std::setw(10) << nPriPart << std::setw(10) << n100MeV
            << std::setw(10) << n500MeV << std::setw(10) << n1GeV
            << std::setw(10) << ntr3DOK << std::setw(10) << ntr3DNO
            << std::setw(10) << nPriPartRecoOK << std::endl;

  c.SaveAs("rms.pdf)");
  fout.cd();
  tv.Write();
  fout.Close();
}

void findvtx()
{
  // find best value of epsilon: 0.5
  // filter out small region of different multiplicity: no
  // start searching vtx from wider region: no
  /*
  processEvents(1, 0.0, false);
  processEvents(1, 0.00000001, false);
  processEvents(1, 0.0000001, false);
  processEvents(1, 0.000001, false);
  processEvents(1, 0.00001, false);
  processEvents(1, 0.0001, false);
  processEvents(1, 0.001, false);
  processEvents(1, 0.01, false);
  processEvents(1, 0.1, false);
  processEvents(1, 0.2, false);
  processEvents(1, 0.3, false);
  processEvents(1, 0.4, false);
  processEvents(1, 0.5, false);
  processEvents(1, 0.6, false);
  processEvents(1, 0.7, false);
  processEvents(1, 0.8, false);
  processEvents(1, 0.9, false);
  processEvents(1, 1., false);
  processEvents(1, 2., false);
  processEvents(1, 0.5, true);
  processEvents(2, 0.5, true);
  processEvents(3, 0.5, true);
  processEvents(4, 0.5, true);
  processEvents(5, 0.5, true);
  processEvents(6, 0.5, true);
  processEvents(7, 0.5, true);
  processEvents(8, 0.5, true);
  processEvents(9, 0.5, true);
  processEvents(10, 0.5, true);
  processEvents(2, 0.5, false);
  processEvents(3, 0.5, false);
  processEvents(4, 0.5, false);
  processEvents(5, 0.5, false);
  processEvents(6, 0.5, false);
  processEvents(7, 0.5, false);
  processEvents(8, 0.5, false);
  processEvents(9, 0.5, false);
  processEvents(10, 0.5, false);*/

  /*
    processEvents(1E7,1E7,1,0,1.,1E7,1E7);
    processEvents(1E7,1E7,1,1,1.,1E7,1E7);
    processEvents(1E7,1E7,1,2,1.,1E7,1E7);
    processEvents(1E7,1E7,1,3,1.,1E7,1E7);
    processEvents(1E7,1E7,1,4,1.,1E7,1E7);
    processEvents(1E7,1E7,1,5,1.,1E7,1E7);
    processEvents(1E7,1E7,1,10,1.,1E7,1E7);

  fraction tol_phi  tol_x tol_mod  mindigtr    dz_tol dn_tol total     fracY
  fracX     avePY     avePX     aveNY     aveNX    3DtrOK    3DtrNO
       1     1e+07  1e+07       1         0     1e+07  1e+07 31001  0.585529
  0.666688  0.951308  0.956523   27.1775   20.2663     83410         0
       1     1e+07  1e+07       1         1     1e+07  1e+07 31001  0.362956
  0.466727  0.858389  0.883189   68.5525   44.3859     29058         0
       1     1e+07  1e+07       1         2     1e+07  1e+07 31001  0.331118
  0.427244   0.85393  0.887535   75.1021   52.5307     25875         0
       1     1e+07  1e+07       1         3     1e+07  1e+07 31001   0.30928
  0.400181  0.851226  0.888905   79.2529   57.0791     23924         0
       1     1e+07  1e+07       1         4     1e+07  1e+07 31001  0.292249
  0.376439  0.849571  0.888555   82.7189   60.4982     22385         0
       1     1e+07  1e+07       1         5     1e+07  1e+07 31001  0.274411
  0.354892  0.847312  0.888689   86.2045    63.515     20962         0
       1     1e+07  1e+07       1        10     1e+07  1e+07 31001  0.218864
  0.271862  0.841057  0.888397   99.0147   75.9628     16139         0
  */
  /*
    processEvents(0.0001,1E7,1,3,1.,1E7,1E7);
    processEvents(0.001,1E7,1,3,1.,1E7,1E7);
    processEvents(0.01,1E7,1,3,1.,1E7,1E7);
    processEvents(0.1,1E7,1,3,1.,1E7,1E7);
    processEvents(0.2,1E7,1,3,1.,1E7,1E7);
    processEvents(0.5,1E7,1,3,1.,1E7,1E7);
    processEvents(1.,1E7,1,3,1.,1E7,1E7);
    processEvents(2.,1E7,1,3,1.,1E7,1E7);
    processEvents(10.,1E7,1,3,1.,1E7,1E7);

  fraction   tol_phi tol_x tol_mod  mindigtr dz_tol    dn_tol     total
  fracY     fracX     avePY     avePX     aveNY     aveNX    3DtrOK    3DtrNO
         1    0.0001 1e+07       1         3  1e+07     1e+07
  310010.00170962  0.400181         1  0.888905   10.4444   57.0791        72 0
         1     0.001 1e+07       1         3  1e+07     1e+07     31001
  0.0626109  0.400181  0.987784  0.888905    17.716   57.0791      3323 0
         1      0.01 1e+07       1         3  1e+07     1e+07     31001
  0.328989  0.400181  0.967039  0.888905   40.0435   57.0791     19067         0
         1       0.1 1e+07       1         3  1e+07     1e+07     31001
  0.421244  0.400181  0.924476  0.888905   63.6038   57.0791     27199         0
         1       0.2 1e+07       1         3  1e+07     1e+07     31001
  0.390149  0.400181  0.902731  0.888905   68.9421   57.0791     26616         0
         1       0.5 1e+07       1         3  1e+07     1e+07     31001
  0.346118  0.400181  0.874969  0.888905   75.4198   57.0791     25103         0
         1         1 1e+07       1         3  1e+07     1e+07     31001
  0.321119  0.400181  0.859481  0.888905   78.2855   57.0791     24226         0
         1         2 1e+07       1         3  1e+07     1e+07     31001
  0.309893  0.400181  0.851923  0.888905   79.2079   57.0791     23940         0
         1        10 1e+07       1         3  1e+07     1e+07     31001
  0.30928  0.400181  0.851226  0.888905   79.2529   57.0791     23924         0
  */

  /*
    processEvents(0.1,1,1,3,1.,1E7,1E7);
    processEvents(0.1,2,1,3,1.,1E7,1E7);
    processEvents(0.1,5,1,3,1.,1E7,1E7);
    processEvents(0.1,10,1,3,1.,1E7,1E7);
    processEvents(0.1,20,1,3,1.,1E7,1E7);
    processEvents(0.1,50,1,3,1.,1E7,1E7);
    processEvents(0.1,100,1,3,1.,1E7,1E7);
    processEvents(0.1,1000,1,3,1.,1E7,1E7);

  fraction   tol_phi tol_x   tol_mod  mindigtr dz_tol dn_tol total     fracY
  fracX     avePY     avePX     aveNY     aveNX 3DtrackOK 3DtrackNO
         1       0.1     1         1         3  1e+07  1e+07 31001  0.421244
  0.0393858  0.924476  0.983681   63.6038   23.0531      2821         0
         1       0.1     2         1         3  1e+07  1e+07 31001  0.421244
  0.0714816  0.924476  0.976192   63.6038   29.4961      5205         0
         1       0.1     5         1         3  1e+07  1e+07 31001  0.421244
  0.150931  0.924476  0.966015   63.6038   37.7544     10768         0
         1       0.1    10         1         3  1e+07  1e+07 31001  0.421244
  0.238896  0.924476   0.95911   63.6038    43.071     16509         0
         1       0.1    20         1         3  1e+07  1e+07 31001  0.421244
  0.334667  0.924476  0.954244   63.6038   48.8484     21732         0
         1       0.1    50         1         3  1e+07  1e+07 31001  0.421244
  0.420599  0.924476  0.945698   63.6038    54.109     25857         0
         1       0.1   100         1         3  1e+07  1e+07 31001  0.421244
  0.43605  0.924476  0.936404   63.6038   55.9506     26848         0
         1       0.1  1000         1         3  1e+07  1e+07 31001  0.421244
  0.407116  0.924476  0.895416   63.6038   57.2511     27194         0
  */

  /*
  processEvents(0.1,100,2,3,1.,1E7,1E7);
  processEvents(0.1,100,3,3,1.,1E7,1E7);
  processEvents(0.1,100,4,3,1.,1E7,1E7);
  processEvents(0.1,100,5,3,1.,1E7,1E7);
  processEvents(0.1,100,10,3,1.,1E7,1E7);
  processEvents(0.1,100,100,3,1.,1E7,1E7);


  fraction   tol_phi tol_x   tol_mod  mindigtr dz_tol dn_tol     total     fracY
  fracX     avePY     avePX     aveNY     aveNX 3DtrackOK 3DtrackNO
         1       0.1   100         2         3  1e+07  1e+07     31001  0.403019
  0.406051  0.915534  0.925119   66.9911   64.5242     25337         0
         1       0.1   100         3         3  1e+07  1e+07     31001   0.39689
  0.394116  0.910148   0.91927   67.5035   66.6706     24943         0
         1       0.1   100         4         3  1e+07  1e+07     31001  0.393503
  0.386858  0.906048  0.914574    67.637   67.5364     24749         0
         1       0.1   100         5         3  1e+07  1e+07     31001  0.390665
  0.38402  0.902626  0.911501   67.6032   67.5336     24746         0
         1       0.1   100        10         3  1e+07  1e+07     31001  0.381762
  0.37702  0.890486  0.898828    67.195    66.542     24996         0
         1       0.1   100       100         3  1e+07  1e+07     31001  0.366053
  0.363859   0.87025   0.86753    65.845   62.0067     25857         0

  */

  processEvents(0.1, 100, 4, 3, 1., 1E7, 1E7);

  gStyle->SetPalette(53);

  TFile f("vtx.root");
  TTree* tv = (TTree*)f.Get("tv");
  TCanvas c;
  TH1D* h;

  TCut vtxReco = "VtxMono == 1 || VtxMulti == 1";
  TCut monoprong = "VtxMono == 1";
  TCut multiprong = "VtxMulti == 1";
  TCut notLAr = "zvtx_true > 22100";

  // x VS z
  tv->Draw("xvtx_true:zvtx_true>>hxVSz", "", "");
  TH2D* hxVSz = (TH2D*)gROOT->FindObject("hxVSz");
  hxVSz->SetTitle(";Z (mm); X (mm); #Deltar_{3D} (mm)");
  hxVSz->Draw();
  c.SaveAs("vtx.pdf(");

  // y VS z
  tv->Draw("yvtx_true:zvtx_true>>hyVSz", "", "");
  TH2D* hyVSz = (TH2D*)gROOT->FindObject("hyVSz");
  hyVSz->SetTitle(";Z (mm); Y (mm); #Deltar_{3D} (mm)");
  hyVSz->Draw();
  c.SaveAs("vtx.pdf");

  // x VS z reco
  tv->Draw("xvtx_reco:zvtx_reco>>hxVSz_reco", vtxReco, "");
  TH2D* hxVSz_reco = (TH2D*)gROOT->FindObject("hxVSz_reco");
  hxVSz_reco->SetTitle(";Z (mm); X (mm); #Deltar_{3D} (mm)");
  hxVSz_reco->Draw();
  c.SaveAs("vtx.pdf");

  // y VS z reco
  tv->Draw("yvtx_reco:zvtx_reco>>hyVSz_reco", vtxReco, "");
  TH2D* hyVSz_reco = (TH2D*)gROOT->FindObject("hyVSz_reco");
  hyVSz_reco->SetTitle(";Z (mm); Y (mm); #Deltar_{3D} (mm)");
  hyVSz_reco->Draw();
  c.SaveAs("vtx.pdf");

  // x VS z VS mean dr
  tv->Draw("xvtx_reco:zvtx_reco>>hxVSz_n(25,21500,26500,25,-2000,2000)",
           "VtxMulti == 1", "colz");
  TH2D* hxVSz_n = new TH2D(*((TH2D*)gROOT->FindObject("hxVSz_n")));
  tv->Draw("xvtx_reco:zvtx_reco>>hxVSz_dv(25,21500,26500,25,-2000,2000)",
           "(VtxMulti == "
           "1)*sqrt(pow(xvtx_reco-xvtx_true,2)+pow(yvtx_reco-yvtx_true,2)+pow("
           "zvtx_reco-zvtx_true,2))",
           "colz");
  TH2D* hxVSz_dv = new TH2D(*((TH2D*)gROOT->FindObject("hxVSz_dv")));
  TH2D hxVSz_meandv = (*hxVSz_dv) / (*hxVSz_n);
  hxVSz_meandv.SetStats(false);
  hxVSz_meandv.SetTitle(";Z (mm); Y (mm); #Deltar_{3D} (mm)");
  hxVSz_meandv.Draw("colz");
  c.SaveAs("vtx.pdf");

  // y VS z VS mean dr
  tv->Draw("yvtx_reco:zvtx_reco>>hyVSz_n(25,21500,26500,25,-5000,0)",
           "VtxMulti == 1", "colz");
  TH2D* hyVSz_n = new TH2D(*((TH2D*)gROOT->FindObject("hyVSz_n")));
  tv->Draw("yvtx_reco:zvtx_reco>>hyVSz_dv(25,21500,26500,25,-5000,0)",
           "(VtxMulti == "
           "1)*sqrt(pow(xvtx_reco-xvtx_true,2)+pow(yvtx_reco-yvtx_true,2)+pow("
           "zvtx_reco-zvtx_true,2))",
           "colz");
  TH2D* hyVSz_dv = new TH2D(*((TH2D*)gROOT->FindObject("hyVSz_dv")));
  TH2D hyVSz_meandv = (*hyVSz_dv) / (*hyVSz_n);
  hyVSz_meandv.SetStats(false);
  hyVSz_meandv.SetTitle(";Z (mm); Y (mm); #Deltar_{3D} (mm)");
  hyVSz_meandv.Draw("colz");
  c.SaveAs("vtx.pdf");

  c.SetLogy(true);
  // dx
  tv->Draw("xvtx_reco-xvtx_true>>dvx_multi(100,-4000,4000)", multiprong, "");
  TH1D* dvx_multi = new TH1D(*((TH1D*)gROOT->FindObject("dvx_multi")));
  dvx_multi->SetFillColorAlpha(kBlue, 0.3);
  std::cout << "dx (multi): " << dvx_multi->GetMean() << " "
            << dvx_multi->GetRMS() << std::endl;

  tv->Draw("xvtx_reco-xvtx_true>>dvx_mono(100,-4000,4000)", monoprong, "");
  TH1D* dvx_mono = new TH1D(*((TH1D*)gROOT->FindObject("dvx_mono")));
  dvx_mono->SetFillColorAlpha(kRed, 0.3);
  std::cout << "dx (mono) : " << dvx_mono->GetMean() << " "
            << dvx_mono->GetRMS() << std::endl;

  THStack dvx_all("dvx_all", ";xreco-xtrue (mm)");
  dvx_all.Add(dvx_multi);
  dvx_all.Add(dvx_mono);

  TLegend leg(0.79, 0.79, 0.99, 0.99);
  leg.AddEntry(dvx_multi, "mtr", "fl");
  leg.AddEntry(dvx_mono, "1tr", "fl");

  dvx_all.Draw("nostack");
  leg.Draw();
  c.SaveAs("vtx.pdf");

  // dy
  tv->Draw("yvtx_reco-yvtx_true>>dvy_multi(100,-5000,5000)", multiprong, "");
  TH1D* dvy_multi = new TH1D(*((TH1D*)gROOT->FindObject("dvy_multi")));
  dvy_multi->SetFillColorAlpha(kBlue, 0.3);
  std::cout << "dy (multi): " << dvy_multi->GetMean() << " "
            << dvy_multi->GetRMS() << std::endl;

  tv->Draw("yvtx_reco-yvtx_true>>dvy_mono(100,-5000,5000)", monoprong, "");
  TH1D* dvy_mono = new TH1D(*((TH1D*)gROOT->FindObject("dvy_mono")));
  dvy_mono->SetFillColorAlpha(kRed, 0.3);
  std::cout << "dy (mono) : " << dvy_mono->GetMean() << " "
            << dvy_mono->GetRMS() << std::endl;

  THStack dvy_all("dvy_all", ";yreco-ytrue (mm)");
  dvy_all.Add(dvy_multi);
  dvy_all.Add(dvy_mono);

  dvy_all.Draw("nostack");
  leg.Draw();
  c.SaveAs("vtx.pdf");

  // dz
  tv->Draw("zvtx_reco-zvtx_true>>dvz_multi(100,-4000,4000)", multiprong, "");
  TH1D* dvz_multi = new TH1D(*((TH1D*)gROOT->FindObject("dvz_multi")));
  dvz_multi->SetFillColorAlpha(kBlue, 0.3);
  std::cout << "dz (multi): " << dvz_multi->GetMean() << " "
            << dvz_multi->GetRMS() << std::endl;

  tv->Draw("zvtx_reco-zvtx_true>>dvz_mono(100,-4000,4000)", monoprong, "");
  TH1D* dvz_mono = new TH1D(*((TH1D*)gROOT->FindObject("dvz_mono")));
  dvz_mono->SetFillColorAlpha(kRed, 0.3);
  std::cout << "dz (mono) : " << dvz_mono->GetMean() << " "
            << dvz_mono->GetRMS() << std::endl;

  THStack dvz_all("dvz_all", ";zreco-ztrue (mm)");
  dvz_all.Add(dvz_multi);
  dvz_all.Add(dvz_mono);

  dvz_all.Draw("nostack");
  leg.Draw();
  c.SaveAs("vtx.pdf");

  // dxy
  tv->Draw(
      "sqrt(pow(xvtx_reco-xvtx_true,2)+pow(yvtx_reco-yvtx_true,2))>>dvabsxy_"
      "multi(200,0,5000)",
      multiprong, "");
  TH1D* dvabsxy_multi = new TH1D(*((TH1D*)gROOT->FindObject("dvabsxy_multi")));
  dvabsxy_multi->SetFillColorAlpha(kBlue, 0.3);
  std::cout << "dr_xy (multi): " << dvabsxy_multi->GetMean() << " "
            << dvabsxy_multi->GetRMS() << std::endl;

  tv->Draw(
      "sqrt(pow(xvtx_reco-xvtx_true,2)+pow(yvtx_reco-yvtx_true,2))>>dvabsxy_"
      "mono(200,0,5000)",
      monoprong, "");
  TH1D* dvabsxy_mono = new TH1D(*((TH1D*)gROOT->FindObject("dvabsxy_mono")));
  dvabsxy_mono->SetFillColorAlpha(kRed, 0.3);
  std::cout << "dr_xy (mono) : " << dvabsxy_mono->GetMean() << " "
            << dvabsxy_mono->GetRMS() << std::endl;

  THStack dvabsxy_all("dvabsxy_all", ";#Deltar_{xy} (mm)");
  dvabsxy_all.Add(dvabsxy_multi);
  dvabsxy_all.Add(dvabsxy_mono);

  dvabsxy_all.Draw("nostack");
  leg.Draw();
  c.SaveAs("vtx.pdf");

  c.SetLogy(false);

  // dxy cumulative
  TH1D hintxy(*dvabsxy_multi);
  TH1D hintxy_mono(*dvabsxy_multi);
  TH1D hintxy_multi(*dvabsxy_multi);
  hintxy.SetFillColor(0);
  hintxy.SetStats(false);
  hintxy_multi.SetFillColor(0);
  hintxy_multi.SetLineColor(kBlue);
  hintxy_multi.SetLineStyle(2);
  hintxy_multi.SetStats(false);
  hintxy_mono.SetFillColor(0);
  hintxy_mono.SetLineColor(kRed);
  hintxy_mono.SetLineStyle(2);
  hintxy_mono.SetStats(false);
  for (int i = 0; i < hintxy.GetNbinsX(); i++) {
    hintxy.SetBinContent(
        i + 1,
        (dvabsxy_multi->Integral(1, i + 1) + dvabsxy_mono->Integral(1, i + 1)) /
            (dvabsxy_multi->Integral() + dvabsxy_mono->Integral()));
    hintxy_mono.SetBinContent(
        i + 1, dvabsxy_mono->Integral(1, i + 1) / dvabsxy_mono->Integral());
    hintxy_multi.SetBinContent(
        i + 1, dvabsxy_multi->Integral(1, i + 1) / dvabsxy_multi->Integral());
  }
  c.DrawFrame(0, 0, 1000, 1,
              "fraction of events with #Deltar_{xy} < "
              "#Deltar^{thr}_{xy};#Deltar^{thr}_{xy} (mm)");
  hintxy.Draw("csame");
  hintxy_mono.Draw("csame");
  hintxy_multi.Draw("csame");

  TLegend leg2(0.69, 0.11, 0.89, 0.31);
  leg2.AddEntry(&hintxy, "all vtx", "l");
  leg2.AddEntry(&hintxy_multi, "mtr", "l");
  leg2.AddEntry(&hintxy_mono, "1tr", "l");
  leg2.Draw();

  c.SaveAs("vtx.pdf");

  c.SetLogy(true);

  // dz
  tv->Draw("abs(zvtx_reco-zvtx_true)>>dvabsz_multi(100,0,4000)", multiprong,
           "");
  TH1D* dvabsz_multi = new TH1D(*((TH1D*)gROOT->FindObject("dvabsz_multi")));
  dvabsz_multi->SetFillColorAlpha(kBlue, 0.3);
  std::cout << "|dz| (multi): " << dvabsz_multi->GetMean() << " "
            << dvabsz_multi->GetRMS() << std::endl;

  tv->Draw("abs(zvtx_reco-zvtx_true)>>dvabsz_mono(100,0,4000)", monoprong, "");
  TH1D* dvabsz_mono = new TH1D(*((TH1D*)gROOT->FindObject("dvabsz_mono")));
  dvabsz_mono->SetFillColorAlpha(kRed, 0.3);
  std::cout << "|dz| (mono) : " << dvabsz_mono->GetMean() << " "
            << dvabsz_mono->GetRMS() << std::endl;

  THStack dvabsz_all("dvabsz_all", ";|zreco-ztrue| (mm)");
  dvabsz_all.Add(dvabsz_multi);
  dvabsz_all.Add(dvabsz_mono);

  dvabsz_all.Draw("nostack");
  leg.Draw();
  c.SaveAs("vtx.pdf");

  c.SetLogy(false);
  // dz cumulative
  TH1D hintz(*dvabsz_multi);
  TH1D hintz_mono(*dvabsz_multi);
  TH1D hintz_multi(*dvabsz_multi);
  hintz.SetFillColor(0);
  hintz.SetStats(false);
  hintz_mono.SetFillColor(0);
  hintz_mono.SetLineColor(kRed);
  hintz_mono.SetLineStyle(2);
  hintz_mono.SetStats(false);
  hintz_multi.SetFillColor(0);
  hintz_multi.SetLineColor(kBlue);
  hintz_multi.SetLineStyle(2);
  hintz_multi.SetStats(false);
  for (int i = 0; i < hintz.GetNbinsX(); i++) {
    hintz.SetBinContent(
        i + 1,
        (dvabsz_multi->Integral(1, i + 1) + dvabsz_mono->Integral(1, i + 1)) /
            (dvabsz_multi->Integral() + dvabsz_mono->Integral()));
    hintz_mono.SetBinContent(
        i + 1, dvabsz_mono->Integral(1, i + 1) / dvabsz_mono->Integral());
    hintz_multi.SetBinContent(
        i + 1, dvabsz_multi->Integral(1, i + 1) / dvabsz_multi->Integral());
  }
  c.DrawFrame(0, 0, 1000, 1,
              "fraction of events with |zreco-ztrue| < "
              "|zreco-ztrue|_{thr};|zreco-ztrue|_{thr} (mm)");
  hintz.Draw("csame");
  hintz_mono.Draw("csame");
  hintz_multi.Draw("csame");
  leg2.Draw();
  c.SaveAs("vtx.pdf");
  c.SetLogy(true);

  // 3D
  tv->Draw(
      "sqrt(pow(xvtx_reco-xvtx_true,2)+pow(yvtx_reco-yvtx_true,2)+pow(zvtx_"
      "reco-zvtx_true,2))>>dv3D_multi(250,0,5500)",
      multiprong, "");
  TH1D* dv3D_multi = new TH1D(*((TH1D*)gROOT->FindObject("dv3D_multi")));
  dv3D_multi->SetFillColorAlpha(kBlue, 0.3);
  std::cout << "|dr_3D| (multi): " << dv3D_multi->GetMean() << " "
            << dv3D_multi->GetRMS() << std::endl;

  tv->Draw(
      "sqrt(pow(xvtx_reco-xvtx_true,2)+pow(yvtx_reco-yvtx_true,2)+pow(zvtx_"
      "reco-zvtx_true,2))>>dv3D_mono(250,0,5500)",
      monoprong, "");
  TH1D* dv3D_mono = new TH1D(*((TH1D*)gROOT->FindObject("dv3D_mono")));
  dv3D_mono->SetFillColorAlpha(kRed, 0.3);
  std::cout << "|dr_3D| (mono) : " << dv3D_mono->GetMean() << " "
            << dv3D_mono->GetRMS() << std::endl;

  THStack dv3D_all("dv3D_all", ";#Deltar_{3D} (mm)");
  dv3D_all.Add(dv3D_multi);
  dv3D_all.Add(dv3D_mono);

  dv3D_all.Draw("nostack");
  leg.Draw();
  c.SaveAs("vtx.pdf");

  c.SetLogy(false);

  // 3D cumulative
  TH1D hint3D(*dv3D_multi);
  TH1D hint3D_mono(*dv3D_multi);
  TH1D hint3D_multi(*dv3D_multi);
  hint3D.SetFillColor(0);
  hint3D.SetStats(false);
  hint3D_mono.SetFillColor(0);
  hint3D_mono.SetLineColor(kRed);
  hint3D_mono.SetLineStyle(2);
  hint3D_mono.SetStats(false);
  hint3D_multi.SetFillColor(0);
  hint3D_multi.SetLineColor(kBlue);
  hint3D_multi.SetLineStyle(2);
  hint3D_multi.SetStats(false);
  for (int i = 0; i < hint3D.GetNbinsX(); i++) {
    hint3D.SetBinContent(
        i + 1,
        (dv3D_multi->Integral(1, i + 1) + dv3D_mono->Integral(1, i + 1)) /
            (dv3D_multi->Integral() + dv3D_mono->Integral()));
    hint3D_mono.SetBinContent(
        i + 1, dv3D_mono->Integral(1, i + 1) / dv3D_mono->Integral());
    hint3D_multi.SetBinContent(
        i + 1, dv3D_multi->Integral(1, i + 1) / dv3D_multi->Integral());
  }
  c.DrawFrame(0, 0, 1000, 1,
              "fraction of events with #Deltar_{3D} < "
              "#Deltar^{thr}_{3D};#Deltar^{thr}_{3D} (mm)");
  hint3D.Draw("csame");
  hint3D_mono.Draw("csame");
  hint3D_multi.Draw("csame");
  leg2.Draw();
  c.SaveAs("vtx.pdf");

  // xy
  tv->Draw(
      "yvtx_reco-yvtx_true:xvtx_reco-xvtx_true>>dvxy(500,-1000,1000,500,-1000,"
      "1000)",
      vtxReco, "colz");
  TH2D* h2D = (TH2D*)gROOT->FindObject("dvxy");
  h2D->SetTitle(";xreco-xtrue (mm);yreco-ytrue (mm)");
  h2D->Draw("colz");
  TPad p("p", "", 0.12, 0.55, 0.45, 0.88);
  p.cd();
  p.DrawFrame(-100, -100, 100, 100);
  h2D->SetStats(false);
  h2D->Draw("colsame");
  c.cd();
  p.Draw();
  c.SaveAs("vtx.pdf)");
}