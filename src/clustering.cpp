#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Linkdef.h"
#include "struct.h"
#include "utils.h"

using namespace std;

bool RepetitionCheck(std::vector<int>, int);

std::pair<vector<int>, vector<int>> FindNeighbours(int[], std::vector<int>, int,
                                                   int, std::vector<int>);

std::pair<int, int> Next(std::vector<int>, int[], int);
std::tuple<double, double, double, double> fit_ls(int, double[], double[],
                                                  double[]);
std::vector<cluster> TrackFit(std::vector<cluster>);
void Clust_info(cluster);
std::vector<cluster> Merge(std::vector<cluster>);
std::vector<cluster> Split(std::vector<cluster>);
std::vector<cluster> RecoverIncomplete(std::vector<cluster>,
                                       std::vector<dg_cell>);
cluster Calc_variables(std::vector<dg_cell>);

double TfromTDC(double t1, double t2, double L);
double AttenuationFactor(double d, int planeID);
double EfromADC(double adc1, double adc2, double d1, double d2, int planeID);
double EfromADCsingle(double adc, double f);
double DfromADC(double, double);

// PROMEMORIA PER FRANCESCO: RICORDATI DI CORREGGERE LE BROKEN CELLS
std::vector<cluster> Clusterize(std::vector<dg_cell>* vec_cellraw)
{
  int dg_size = vec_cellraw->size();
  int dg_idvec[dg_size], i_id = 0;
  std::vector<dg_cell> incomplete_cells, complete_cells, broken_cells,
      weird_cells;
  std::vector<cluster> vec_clust;
  // Create vector of complete and incomplete cells
  for (auto const& cell : *vec_cellraw) {
    // for (int i = 0; i < dg_size; i++) {
    if (cell.ps1.size() == 0 && cell.ps2.size() == 0) {
      continue;
    } else if (cell.ps1.size() == 0 || cell.ps2.size() == 0) {
      // cout << "broken" << endl;
      broken_cells.push_back(cell);
    } /*
    if ((cell.ps1.adc.at(0) == 0 || cell.ps2.adc.at(0) == 0 ||
    cell.ps1.tdc.at(0) == 0 || cell.ps2.tdc.at(0) == 0)) {
        cout << "incomplete" << endl;
        //Incomplete cell
        incomplete_cells.push_back(cell);
    }
    if (abs(cell.ps1.tdc.at(0) - cell.ps2.tdc.at(0)) > 30) {
        cout << "incomplete time" << endl;
        weird_cells.push_back(cell);
    }*/
    else {
      // Complete cell
      complete_cells.push_back(cell);
      dg_idvec[i_id] = cell.id;
      i_id++;
    }
  }
  int cluster[i_id], n_cluster = 0;
  std::vector<int> checked_array;
  std::vector<int> vec_cell;
  for (int j = 0; j < complete_cells.size(); j++) {
    if (j == 0) {
      cluster[j] = dg_idvec[0];
    } else {
      std::pair<int, int> clu_nclu =
          Next(checked_array, dg_idvec, complete_cells.size());
      cluster[j] = clu_nclu.first;
      if (cluster[j] == 0) {
        break;
      }
      n_cluster = clu_nclu.second;
    }
    vec_cell.push_back(n_cluster);
    checked_array.push_back(cluster[j]);
    std::pair<vector<int>, vector<int>> check_cell =
        FindNeighbours(dg_idvec, checked_array, i_id, cluster[j], vec_cell);
    checked_array = check_cell.first;
    vec_cell = check_cell.second;
    std::vector<int> vec_cellid;
    for (std::vector<int>::const_iterator it = vec_cell.begin();
         it != vec_cell.end(); it++) {
      vec_cellid.push_back(*it);
    }
    std::vector<dg_cell> cluster_cells;
    for (int its = 0; its < vec_cellid.size(); its++) {
      int cell = vec_cellid[its];
      cluster_cells.push_back(complete_cells.at(cell));
    }
    vec_cellid.clear();
    struct cluster clust;
    clust = Calc_variables(cluster_cells);
    vec_clust.push_back(clust);
    vec_cell.clear();
  }
  // SPLIT
  vec_clust = Split(vec_clust);
  // MERGE
  vec_clust = Merge(vec_clust);
  // Track Fit
  vec_clust = TrackFit(vec_clust);
  // std::vector<cluster> Og_clus = vec_clust;
  vec_clust = RecoverIncomplete(vec_clust, broken_cells);
  return vec_clust;
}

int Clustering(std::string const& input)
{
  gSystem->Load("libStruct.so");
  const char* finname = input.c_str();
  TFile f(finname, "READ");
  TTree* t = (TTree*)f.Get("tDigit");
  int nEvents = t->GetEntries();
  std::vector<dg_cell>* cell = new std::vector<dg_cell>;
  std::vector<cluster> f_clust, og_clust;
  TString output = input;
  output.ReplaceAll(".root", ".Clusters.root");
  TFile fout(output, "RECREATE");
  TTree tout("tCluster", "Clustering");
  tout.Branch("cluster", "std::vector<cluster>", &f_clust);
  t->SetBranchAddress("dg_cell", &cell);
  for (int i = 0; i < nEvents; i++) {
    cout << "--------------" << endl;
    cout << "Entry " << i << endl;
    t->GetEntry(i);
    std::vector<cluster> clust = Clusterize(std::move(cell));
     for (auto const& clu_info : clust) {
        Clust_info(clu_info);
    }
    f_clust = clust;
    tout.Fill();
    clust.clear();
    f_clust.clear();
    og_clust.clear();
  }
  fout.cd();
  tout.Write();
  fout.Close();
  delete cell;
  return 0;
}

std::vector<cluster> RecoverIncomplete(std::vector<cluster> clus,
                                       std::vector<dg_cell> incomplete_cells)
{
  for (auto const& brok_cells : incomplete_cells) {
    int isbarrel = 0;
    if (brok_cells.mod == 30) isbarrel = 1;
    if (brok_cells.mod == 40) isbarrel = 2;
    double cell_phi =
        atan((brok_cells.z - 23910.00) / (brok_cells.y + 2384.73)) * 180 /
        TMath::Pi();
    double cell_theta =
        atan((brok_cells.z - 23910.00) / (brok_cells.x)) * 180 /
        TMath::Pi();
    int minentry = 0;
    int found = 0;
    for (int j = 0; j < clus.size(); j++) {
      double rec_en = 0;
      double clus_phi =
          atan((clus.at(j).z - 23910.00) / (clus.at(j).y + 2384.73)) * 180 /
          TMath::Pi();
      double clus_theta =
          atan((clus.at(j).z - 23910.00) / (clus.at(j).x)) * 180 /
          TMath::Pi();
      double minphi = 999, mintheta = 999;
      int isbarrelc = 0;
      if (clus.at(j).cells[0].mod == 30) isbarrelc = 1;
      if (clus.at(j).cells[0].mod == 40) isbarrelc = 2;
      if (abs(cell_phi - clus_phi) < 3 && isbarrel == 0) {
        found = 1;
        if (abs(cell_phi - clus_phi) < minphi) {
          minphi = abs(cell_phi - clus_phi);
          minentry = j;
        }
        continue;
      }
      if (isbarrel == isbarrelc) {
        if (abs(cell_theta - clus_theta) < 3) {
          found = 1;
          mintheta = abs(cell_theta - clus_theta);
          minentry = j;
        }
        continue;
      }
    }
    if (found == 1 && isbarrel == 0) {
      double rec_en = 0;
      double DpmA = clus.at(minentry).x / 10 + 215;
      double DpmB = -clus.at(minentry).x / 10 + 215;
      if (brok_cells.ps1.size() != 0 && brok_cells.ps2.size() != 0) {
        double Ea = brok_cells.ps1.at(0).adc;
        double Eb = brok_cells.ps2.at(0).adc;
        rec_en = kloe_simu::EfromADC(Ea, Eb, DpmA, DpmB, brok_cells.lay);
        clus.at(minentry).e = clus.at(minentry).e + rec_en;
        clus.at(minentry).cells.push_back(brok_cells);
      } else if (brok_cells.ps1.size() != 0) {
        int laycell = brok_cells.lay;
        double f = kloe_simu::AttenuationFactor(DpmA, laycell);
        rec_en = EfromADCsingle(brok_cells.ps1.at(0).adc, f);
        clus.at(minentry).e = clus.at(minentry).e + rec_en;
        clus.at(minentry).cells.push_back(brok_cells);
      } else if (brok_cells.ps2.size() != 0) {
        int laycell = brok_cells.lay;
        double f = kloe_simu::AttenuationFactor(DpmB, laycell);
        rec_en = EfromADCsingle(brok_cells.ps2.at(0).adc, f);
        clus.at(minentry).e = clus.at(minentry).e + rec_en;
        clus.at(minentry).cells.push_back(brok_cells);
      }
    }
    if (found == 1 && isbarrel != 0) {
      double rec_en = 0;
      double ecl = brok_cells.l;
      double DpmA = clus.at(minentry).z / 10 + ecl / 20;
      double DpmB = -clus.at(minentry).z / 10 + ecl / 20;
      if (brok_cells.ps1.size() != 0 && brok_cells.ps2.size() != 0) {
        double Ea = brok_cells.ps1.at(0).adc;
        double Eb = brok_cells.ps2.at(0).adc;
        rec_en = kloe_simu::EfromADC(Ea, Eb, DpmA, DpmB, brok_cells.lay);
        clus.at(minentry).e = clus.at(minentry).e + rec_en;
        clus.at(minentry).cells.push_back(brok_cells);
      } else if (brok_cells.ps1.size() != 0) {
        int laycell = brok_cells.lay;
        double f = kloe_simu::AttenuationFactor(DpmA, laycell);
        rec_en = EfromADCsingle(brok_cells.ps1.at(0).adc, f);
        clus.at(minentry).e = clus.at(minentry).e + rec_en;
        clus.at(minentry).cells.push_back(brok_cells);
      } else if (brok_cells.ps2.size() != 0) {
        int laycell = brok_cells.lay;
        double f = kloe_simu::AttenuationFactor(DpmB, laycell);
        rec_en = EfromADCsingle(brok_cells.ps2.at(0).adc, f);
        clus.at(minentry).e = clus.at(minentry).e + rec_en;
        clus.at(minentry).cells.push_back(brok_cells);
      }
    }
  }
  return clus;
}

void Clust_info(cluster clus)
{
  cout << "=O=O=O=O=O=O=O=O=O=O=O=O=O" << endl;
  cout << "New cluster: Energy " << clus.e << " MeV" << endl;
  cout << "Coordinate centroide: " << clus.x << " [X] " << clus.y << " [Y] "
       << clus.z << " [z] e tempo di arrivo medio " << clus.t << " ns" << endl;
  cout << "Varianza: " << clus.varx << " [X] " << clus.vary << " [Y] "
       << clus.varz << " [z]" << endl;

  double phi = atan((clus.z - 23910.00) / (clus.y + 2384.73)) * 180 /
               3.14159265358979323846;
  cout << "Composto dalle seguenti celle: ";
  for (int i = 0; i < clus.cells.size(); i++) {
    cout << clus.cells.at(i).id << " ";
  }
  cout << endl;
}

std::vector<cluster> Split(std::vector<cluster> original_clu_vec)
{
  std::vector<cluster> clu_vec;
  std::vector<dg_cell> all_cells;
  for (auto const& clus : original_clu_vec) {
    int splitted = 0;
    double tA = 0, tB = 0, tA2 = 0, tB2 = 0;
    double EA, EB, EAtot = 0, EBtot = 0, EA2tot = 0, EB2tot = 0;
    double tRMS_A, tRMS_B, dist;
    all_cells = clus.cells;
    for (int j = 0; j < all_cells.size(); j++) {
      EA = all_cells.at(j).ps1.at(0).adc;
      EB = all_cells.at(j).ps2.at(0).adc;
      EAtot += EA;
      EA2tot += EA * EA;
      EBtot += EB;
      EB2tot += EB * EB;
      double d =
          DfromADC(all_cells[j].ps1.at(0).tdc, all_cells[j].ps2.at(0).adc);
      double d1, d2, d3;
      d1 = 0.5 * all_cells[j].l + d;
      d2 = 0.5 * all_cells[j].l - d;
      double cell_E = kloe_simu::EfromADC(all_cells[j].ps1.at(0).adc,
                                          all_cells[j].ps2.at(0).adc, d1, d2,
                                          all_cells[j].lay);
      tA += (all_cells.at(j).ps1.at(0).tdc -
             kloe_simu::vlfb * d1 / kloe_simu::m_to_mm) *
            EA;
      tA2 += std::pow(all_cells.at(j).ps1.at(0).tdc -
                          kloe_simu::vlfb * d1 / kloe_simu::m_to_mm,
                      2) *
             EA;
      tB += (all_cells.at(j).ps2.at(0).tdc -
             kloe_simu::vlfb * d2 / kloe_simu::m_to_mm) *
            EB;
      tB2 += std::pow(all_cells.at(j).ps2.at(0).adc -
                          kloe_simu::vlfb * d2 / kloe_simu::m_to_mm,
                      2) *
             EB;
    }
    tA = tA / EAtot;
    tA2 = tA2 / EAtot;
    tB = tB / EAtot;
    tB2 = tB2 / EBtot;
    tRMS_A = (tA2 - tA * tA) * (EA2tot - EAtot * EAtot) / EA2tot;
    tRMS_B = (tB2 - tB * tB) * (EB2tot - EBtot * EBtot) / EB2tot;
    dist = std::sqrt(tRMS_A * tRMS_A + tRMS_B * tRMS_B);
    if (dist > 5) {
      std::vector<dg_cell> q1_cells, q2_cells, q3_cells, q4_cells;
      for (auto const& a_cells : all_cells) {
        double t_difA = a_cells.ps1.at(0).tdc - tA;
        double t_difB = a_cells.ps2.at(0).tdc - tB;
        if (t_difA > 0) {
          if (t_difB > 0) {
            q1_cells.push_back(a_cells);
          } else {
            q2_cells.push_back(a_cells);
          }
        } else {
          if (t_difB > 0) {
            q3_cells.push_back(a_cells);
          } else {
            q4_cells.push_back(a_cells);
          }
        }
      }
      if (q1_cells.size() != 0) {
        cluster clus = Calc_variables(q1_cells);
        clu_vec.push_back(clus);
        splitted++;
      }
      if (q2_cells.size() != 0) {
        cluster clus = Calc_variables(q2_cells);
        clu_vec.push_back(clus);
        splitted++;
      }
      if (q3_cells.size() != 0) {
        cluster clus = Calc_variables(q3_cells);
        clu_vec.push_back(clus);
        splitted++;
      }
      if (q4_cells.size() != 0) {
        cluster clus = Calc_variables(q4_cells);
        clu_vec.push_back(clus);
        splitted++;
      }
      q1_cells.clear();
      q2_cells.clear();
      q3_cells.clear();
      q4_cells.clear();
    } else {
      clu_vec.push_back(clus);
    }
    all_cells.clear();
  }
  original_clu_vec.clear();
  return clu_vec;
}

std::vector<cluster> Merge(std::vector<cluster> Og_cluster)
{
  std::vector<cluster> mgd_cluster;
  std::vector<int> checked;
  for (int i = 0; i < Og_cluster.size(); i++) {
    double xi = Og_cluster.at(i).x;
    double yi = Og_cluster.at(i).y;
    double zi = Og_cluster.at(i).z;
    double ti = Og_cluster.at(i).t;
    double ei = Og_cluster.at(i).e;
    bool RepCheck = RepetitionCheck(checked, i);
    if (RepCheck == true) {
      continue;
    }
    checked.push_back(i);
    cluster clust;
    clust.x = Og_cluster.at(i).x;
    clust.y = Og_cluster.at(i).y;
    clust.z = Og_cluster.at(i).z;
    clust.t = Og_cluster.at(i).t;
    clust.e = Og_cluster.at(i).e;
    clust.varx = Og_cluster.at(i).varx;
    clust.vary = Og_cluster.at(i).vary;
    clust.varz = Og_cluster.at(i).varz;
    clust.cells = Og_cluster.at(i).cells;
    for (int j = i; j < Og_cluster.size(); j++) {
      RepCheck = RepetitionCheck(checked, j);
      if (RepCheck == true) {
        continue;
      }
      double xj = Og_cluster.at(j).x;
      double yj = Og_cluster.at(j).y;
      double zj = Og_cluster.at(j).z;
      double tj = Og_cluster.at(j).t;
      double ej = Og_cluster.at(j).e;
      double D = sqrt((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj) +
                      (zi - zj) * (zi - zj));
      double DT = abs(ti - tj);
      if (D < 40 && DT < 2.5) {
        bool endcap = false;
        if (clust.cells[0].id > 25000) {
          endcap = true;
        }
        if (endcap == true) {
          double Dz_ec = abs(yi - yj);
          D = sqrt((xi - xj) * (xi - xj) + (zi - zj) * (zi - zj));
          if (Dz_ec < 30 && D < 40) {
            std::vector<dg_cell> vec_cells_j = Og_cluster.at(j).cells;
            for (int k = 0; k < vec_cells_j.size(); k++) {
              clust.cells.push_back(vec_cells_j.at(k));
            }
            clust = Calc_variables(clust.cells);
            checked.push_back(j);
          }
        } else if (endcap == false) {
          double Dz_bar = abs(xi - xj);
          D = sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj));
          if (Dz_bar < 30 && D < 40) {
            std::vector<dg_cell> vec_cells_j = Og_cluster.at(j).cells;
            for (int k = 0; k < vec_cells_j.size(); k++) {
              clust.cells.push_back(vec_cells_j.at(k));
            }
            clust = Calc_variables(clust.cells);
            checked.push_back(j);
          }
        }
      }
    }
    mgd_cluster.push_back(clust);
  }
  return mgd_cluster;
}

std::vector<cluster> TrackFit(std::vector<cluster> clu_vec)
{
  const double xl[5] = {4.44, 4.44, 4.44, 4.44, 5.24};
  for (int i = 0; i < clu_vec.size(); i++) {
    double apx[3], eapx[3] = {0, 0, 0}, ctrk[3] = {0, 0, 0},
                   ectrk[3] = {0, 0, 0};
    std::vector<dg_cell> cell_vec_0, cell_vec_1, cell_vec_2, cell_vec_3,
        cell_vec_4;
    for (int j = 0; j < clu_vec.at(i).cells.size(); j++) {
      if (clu_vec.at(i).cells.at(j).lay == 0) {
        cell_vec_0.push_back(clu_vec.at(i).cells.at(j));
      } else if (clu_vec.at(i).cells.at(j).lay == 1) {
        cell_vec_1.push_back(clu_vec.at(i).cells.at(j));
      } else if (clu_vec.at(i).cells.at(j).lay == 2) {
        cell_vec_2.push_back(clu_vec.at(i).cells.at(j));
      } else if (clu_vec.at(i).cells.at(j).lay == 3) {
        cell_vec_3.push_back(clu_vec.at(i).cells.at(j));
      } else if (clu_vec.at(i).cells.at(j).lay == 4) {
        cell_vec_4.push_back(clu_vec.at(i).cells.at(j));
      }
    }
    cluster Lay0, Lay1, Lay2, Lay3, Lay4;
    Lay0 = Calc_variables(cell_vec_0);
    Lay1 = Calc_variables(cell_vec_1);
    Lay2 = Calc_variables(cell_vec_2);
    Lay3 = Calc_variables(cell_vec_3);
    Lay4 = Calc_variables(cell_vec_4);
    double LayE[5] = {Lay0.e, Lay1.e, Lay2.e, Lay3.e, Lay4.e};
    bool isBarrel = true;
    double yx[5] = {0, 0, 0, 0, 0}, yy[5] = {0, 0, 0, 0, 0},
           yz[5] = {0, 0, 0, 0, 0}, wx[5] = {0, 0, 0, 0, 0},
           wy[5] = {0, 0, 0, 0, 0}, wz[5] = {0, 0, 0, 0, 0};
    double X[5] = {0, 0, 0, 0, 0}, D;
    if (clu_vec.at(i).cells[0].id > 25000) {
      isBarrel = false;
    }
    int lay_cross = 0, first_lay = 0;
    if (Lay0.e > 0) {
      lay_cross++;
      if (lay_cross == 1) {
        first_lay = 1;
      }
      yx[lay_cross - 1] = Lay0.x;
      wx[lay_cross - 1] = 0.6;
      wy[lay_cross - 1] = 0.6;
      wz[lay_cross - 1] = 0.001 * Lay0.e;
      yy[lay_cross - 1] = Lay0.y;
      yz[lay_cross - 1] = Lay0.z;
      if (isBarrel == false) {
        wy[lay_cross - 1] = wz[lay_cross - 1];
        wz[lay_cross - 1] = 0.6;
      }
      D = D + xl[0];
      if (first_lay == 1) X[0] = 0.5 * xl[0];
    }
    if (Lay1.e > 0) {
      lay_cross++;
      yx[lay_cross - 1] = Lay1.x;
      yy[lay_cross - 1] = Lay1.y;
      yz[lay_cross - 1] = Lay1.z;
      wx[lay_cross - 1] = 0.001 * Lay1.e;
      wy[lay_cross - 1] = 0.6;
      wz[lay_cross - 1] = 0.6;
      if (lay_cross == 1) {
        first_lay = 2;
      }
      if (isBarrel == false) {
        wz[lay_cross - 1] = wx[lay_cross - 1];
        wx[lay_cross - 1] = 0.6;
      }
      D = D + xl[1];
      if (first_lay == 2)
        X[0] = 0.5 * xl[1];
      else if (first_lay == 1)
        X[1] = xl[0] + 0.5 * xl[1];
    }
    if (Lay2.e > 0) {
      lay_cross++;
      yx[lay_cross - 1] = Lay2.x;
      yy[lay_cross - 1] = Lay2.y;
      yz[lay_cross - 1] = Lay2.z;
      wx[lay_cross - 1] = 0.001 * Lay2.e;
      wy[lay_cross - 1] = 0.6;
      wz[lay_cross - 1] = 0.6;
      if (lay_cross == 1) {
        first_lay = 3;
      }
      if (isBarrel == false) {
        wz[lay_cross - 1] = wx[lay_cross - 1];
        wx[lay_cross - 1] = 0.6;
      }
      D = D + xl[2];
      if (first_lay == 3)
        X[0] = 0.5 * xl[2];
      else if (first_lay == 2)
        X[1] = xl[0] + 0.5 * xl[2];
      else if (first_lay == 1)
        X[2] = 2 * xl[0] + 0.5 * xl[3];
    }
    if (Lay3.e > 0) {
      lay_cross++;
      yx[lay_cross - 1] = Lay3.x;
      yy[lay_cross - 1] = Lay3.y;
      yz[lay_cross - 1] = Lay3.z;
      wx[lay_cross - 1] = 0.001 * Lay3.e;
      wy[lay_cross - 1] = 0.6;
      wz[lay_cross - 1] = 0.6;
      if (lay_cross == 1) {
        first_lay = 4;
      }
      if (isBarrel == false) {
        wz[lay_cross - 1] = wx[lay_cross - 1];
        wx[lay_cross - 1] = 0.6;
      }
      D = D + xl[3];
      if (first_lay == 4)
        X[0] = 0.5 * xl[3];
      else if (first_lay == 3)
        X[1] = xl[0] + 0.5 * xl[3];
      else if (first_lay == 2)
        X[2] = 2 * xl[0] + 0.5 * xl[3];
      else if (first_lay == 1)
        X[3] = 3 * xl[0] + 0.5 * xl[3];
    }
    if (Lay4.e > 0) {
      lay_cross++;
      yx[lay_cross - 1] = Lay4.x;
      yy[lay_cross - 1] = Lay4.y;
      yz[lay_cross - 1] = Lay4.z;
      wx[lay_cross - 1] = 0.001 * Lay4.e;
      wy[lay_cross - 1] = 0.6;
      wz[lay_cross - 1] = 0.6;
      if (lay_cross == 1) {
        first_lay = 5;
      }
      if (isBarrel == false) {
        wz[lay_cross - 1] = wx[lay_cross - 1];
        wx[lay_cross - 1] = 0.6;
      }
      D = D + xl[4];
      if (first_lay == 5)
        X[0] = 0.5 * xl[4];
      else if (first_lay == 4)
        X[1] = xl[3] + 0.5 * xl[4];
      else if (first_lay == 3)
        X[2] = 2 * xl[3] + 0.5 * xl[4];
      else if (first_lay == 2)
        X[3] = 3 * xl[3] + 0.5 * xl[4];
      else if (first_lay == 1)
        X[4] = 4 * xl[3] + 0.5 * xl[4];
    }
    if (lay_cross == 0) {
      continue;
    }
    D = D - 0.5 * xl[lay_cross];
    int Q = 0, L = 0;
    for (int k = first_lay; k <= 5; k++) {
      if (Q == 0 && (LayE[k - 1] >= 0.05 * clu_vec.at(i).e)) {
        Q = k;
      }
      L++;
    }
    double E1 = 0, E2 = 0;
    if (lay_cross > 1) {
      for (int k_i = 5; k_i >= Q; k_i--) {
        E2 = E1;
        E1 = E1 + LayE[k_i - 1];
      }
      double Rk = E2 / E1;
      double B = 3;
      if (clu_vec.at(i).e > 16.5) B = 1.5 / log(clu_vec.at(i).e / 10);
      double Zmin = 0;
      double Zapx = 0.5 * B * xl[Q - 1];
      double Zmax = B * xl[Q - 1];
      double R1;
      for (int j = 0; j < 4; j++) {
        R1 = exp(-Zapx) * (1 + Zapx);
        if (R1 > Rk) {
          Zmin = Zapx;
          Zapx = 0.5 * (Zmax + Zmin);
        } else if (R1 < Rk) {
          Zmax = Zapx;
          Zapx = 0.5 * (Zmax + Zmin);
        }
      }
      Zapx = -Zapx / B;
      for (int j = 0; j < Q; j++) {
        Zapx = Zapx + xl[j];
      }
      for (int j = 0; j < lay_cross; j++) {
        X[j] = X[j] - Zapx;
      }
      std::tuple<double, double, double, double> fit_varx =
          fit_ls(lay_cross, X, yx, wx);
      std::tuple<double, double, double, double> fit_vary =
          fit_ls(lay_cross, X, yy, wy);
      std::tuple<double, double, double, double> fit_varz =
          fit_ls(lay_cross, X, yz, wz);
      double trktot = sqrt(get<1>(fit_varx) * get<1>(fit_varx) +
                           get<1>(fit_vary) * get<1>(fit_vary) +
                           get<1>(fit_varz) * get<1>(fit_varz));
      ctrk[0] = get<1>(fit_varx) / trktot;
      ctrk[1] = get<1>(fit_vary) / trktot;
      ctrk[2] = get<1>(fit_varz) / trktot;
      ectrk[0] = get<3>(fit_varx) / trktot;
      ectrk[1] = get<3>(fit_vary) / trktot;
      ectrk[2] = get<3>(fit_varz) / trktot;
      apx[0] = get<0>(fit_varx);
      apx[1] = get<0>(fit_vary);
      apx[2] = get<0>(fit_varz);
      eapx[0] = get<2>(fit_varx);
      eapx[1] = get<2>(fit_vary);
      eapx[2] = get<2>(fit_varz);
      // cout << "Apx: [" << apx[0] << "; " << apx[1] << "; " << apx[2] << "] "
      // << " e direzione vettore: [" << ctrk[0] << "; " << ctrk[1] << "; " <<
      // ctrk[2] << "] " << endl;
    }
    if (lay_cross == 1) {
      // cout << "apx(x)=" << yx[0] << " apx(y)=" << yy[0] << " apx(z)=" <<
      // yz[0] << endl;
      // cout << "e_apx(x)=" << sqrt(1/wx[0]) << " e_apx(y)=" << sqrt(1 / wy[0])
      // << " e_apx(z)=" << sqrt(1 / wz[0]) << endl;
      apx[0] = yx[0];
      apx[1] = yy[0];
      apx[2] = yz[0];
    }
    clu_vec.at(i).x = apx[0];
    clu_vec.at(i).y = apx[1];
    clu_vec.at(i).z = apx[2];
    clu_vec.at(i).sx = ctrk[0];
    clu_vec.at(i).sy = ctrk[1];
    clu_vec.at(i).sz = ctrk[2];
  }
  return clu_vec;
}

std::tuple<double, double, double, double> fit_ls(int lay, double X[lay],
                                                  double Y[lay], double W[lay])
{
  double norm = 0, xa = 0, ya = 0, xya = 0, x2a = 0;
  double det, A, B, dA, dB;
  for (int i = 0; i < lay; i++) {
    norm = norm + W[i];
    xa = xa + W[i] * X[i];
    ya = ya + W[i] * Y[i];
    xya = xya + W[i] * Y[i] * X[i];
    x2a = x2a + W[i] * X[i] * X[i];
  }
  norm = norm / lay;
  xa = xa / lay;
  ya = ya / lay;
  xya = xya / lay;
  x2a = x2a / lay;
  det = x2a * norm - xa * xa;
  B = (norm * xya - xa * ya) / (norm * x2a - xa * xa);
  A = ya / norm - B * xa / norm;
  dB = 1 / sqrt(lay * det);
  dA = sqrt(x2a / (lay * det));
  return std::make_tuple(A, B, dA, dB);
}

cluster Calc_variables(std::vector<dg_cell> cells)
{
  double x_weighted = 0, y_weighted = 0, z_weighted = 0, t_weighted = 0,
         x2_weighted = 0, y2_weighted = 0, z2_weighted = 0, Etot = 0, E2tot,
         EvEtot = 0, EA, EAtot = 0, EB, EBtot = 0, TA = 0, TB = 0;
  for (auto const& cell : cells) {
    double d = DfromADC(cell.ps1.at(0).tdc, cell.ps2.at(0).tdc);
    double d1, d2, d3;
    d1 = 0.5 * cell.l + d;
    d2 = 0.5 * cell.l - d;
    double cell_E = kloe_simu::EfromADC(cell.ps1.at(0).adc, cell.ps2.at(0).adc,
                                        d1, d2, cell.lay);
    double cell_T =
        kloe_simu::TfromTDC(cell.ps1.at(0).tdc, cell.ps2.at(0).tdc, cell.l);
    if (cell.mod > 25) {
      d3 = cell.y - d;
      y_weighted = y_weighted + (d3 * cell_E);
      y2_weighted = y2_weighted + (d3 * d3 * cell_E);
      x_weighted = x_weighted + (cell.x * cell_E);
      x2_weighted = x2_weighted + (cell.x * cell.x * cell_E);
    } else {
      d3 = cell.x - d;
      x_weighted = x_weighted + (d3 * cell_E);
      x2_weighted = x2_weighted + (d3 * d3 * cell_E);
      y_weighted = y_weighted + (cell.y * cell_E);
      y2_weighted = y2_weighted + (cell.y * cell.y * cell_E);
    }
    t_weighted = t_weighted + cell_T * cell_E;
    z_weighted = z_weighted + (cell.z * cell_E);
    z2_weighted = z2_weighted + (cell.z * cell.z * cell_E);
    Etot = Etot + cell_E;
    E2tot = E2tot + cell_E * cell_E;
  }
  x_weighted = x_weighted / Etot;
  x2_weighted = x2_weighted / Etot;
  y_weighted = y_weighted / Etot;
  y2_weighted = y2_weighted / Etot;
  z_weighted = z_weighted / Etot;
  z2_weighted = z2_weighted / Etot;
  t_weighted = t_weighted / Etot;
  if (x_weighted > -0.000001 && x_weighted < 0.000001) x_weighted = 0;
  if (y_weighted > -0.000001 && y_weighted < 0.000001) y_weighted = 0;
  if (z_weighted > -0.000001 && z_weighted < 0.000001) z_weighted = 0;
  double dx, dy, dz;
  double neff = Etot * Etot / E2tot;
  double dum = neff / (neff - 1);
  if (cells.size() == 1) {
    dx = 0;
    dy = 0;
    dz = 0;
  } else {
    if (x2_weighted - x_weighted * x_weighted < 0) {
      dx = 0;
    } else {
      dx = sqrt(dum * (x2_weighted - x_weighted * x_weighted));
    }
    if (y2_weighted - y_weighted * y_weighted < 0) {
      dy = 0;
    } else {
      dy = sqrt(dum * (y2_weighted - y_weighted * y_weighted));
    }
    if (z2_weighted - z_weighted * z_weighted < 0) {
      dz = 0;
    } else {
      dz = sqrt(dum * (z2_weighted - z_weighted * z_weighted));
    }
  }
  cluster clust;
  clust.e = Etot;
  clust.x = x_weighted;
  clust.y = y_weighted;
  clust.z = z_weighted;
  clust.t = t_weighted;
  clust.varx = dx;
  clust.vary = dy;
  clust.varz = dz;
  clust.cells = cells;
  return clust;
}

bool RepetitionCheck(std::vector<int> v, int check)
{
  if (std::find(v.begin(), v.end(), check) != v.end()) {
    return true;
  } else {
    return false;
  }
}

std::pair<std::vector<int>, std::vector<int>> FindNeighbours(
    int digits[], std::vector<int> already_checked, int size, int check,
    std::vector<int> vec_cell)
{
  for (int i = 0; i < size; i++) {
    int id_check = digits[i];
    if (id_check == check) {
      continue;
    }
    // Middle of the module
    else if (id_check == check + 101 || id_check == check + 100 ||
             id_check == check + 99 || id_check == check + 1 ||
             id_check == check - 1 || id_check == check - 99 ||
             id_check == check - 100 || id_check == check - 101) {
      bool RepCheck = RepetitionCheck(already_checked, id_check);
      if (RepCheck == true) {
        continue;
      } else {
        vec_cell.push_back(i);
        already_checked.push_back(id_check);
        std::pair<vector<int>, vector<int>> check_cell =
            FindNeighbours(digits, already_checked, size, id_check, vec_cell);
        already_checked = check_cell.first;
        vec_cell = check_cell.second;
      }
    }
    // Right edge of the module
    else if (id_check == check - 889 || id_check == check - 989 ||
             id_check == check - 1089) {
      bool RepCheck = RepetitionCheck(already_checked, id_check);
      if (RepCheck == true) {
        continue;
      } else {
        vec_cell.push_back(i);
        already_checked.push_back(id_check);
        std::pair<vector<int>, vector<int>> check_cell =
            FindNeighbours(digits, already_checked, size, id_check, vec_cell);
        already_checked = check_cell.first;
        vec_cell = check_cell.second;
      }
    }
    // Right edge of 0 module
    else if (id_check == check + 23111 || id_check == check + 23011 ||
             id_check == check + 22911) {
      bool RepCheck = RepetitionCheck(already_checked, id_check);
      if (RepCheck == true) {
        continue;
      } else {
        vec_cell.push_back(i);
        already_checked.push_back(id_check);
        std::pair<vector<int>, vector<int>> check_cell =
            FindNeighbours(digits, already_checked, size, id_check, vec_cell);
        already_checked = check_cell.first;
        vec_cell = check_cell.second;
      }
    }
    // Left edge of the module
    else if (id_check == check + 1089 || id_check == check + 989 ||
             id_check == check + 889) {
      bool RepCheck = RepetitionCheck(already_checked, id_check);
      if (RepCheck == true) {
        continue;
      } else {
        vec_cell.push_back(i);
        already_checked.push_back(id_check);
        std::pair<vector<int>, vector<int>> check_cell =
            FindNeighbours(digits, already_checked, size, id_check, vec_cell);
        already_checked = check_cell.first;
        vec_cell = check_cell.second;
      }
    }
    // Left edge of the 23 module
    else if (id_check == check - 23011 || id_check == check - 22911 ||
             id_check == check - 23111) {
      bool RepCheck = RepetitionCheck(already_checked, id_check);
      if (RepCheck == true) {
        continue;
      } else {
        vec_cell.push_back(i);
        already_checked.push_back(id_check);
        std::pair<vector<int>, vector<int>> check_cell =
            FindNeighbours(digits, already_checked, size, id_check, vec_cell);
        already_checked = check_cell.first;
        vec_cell = check_cell.second;
      }
    }
  }
  return std::make_pair(already_checked, vec_cell);
}

std::pair<int, int> Next(std::vector<int> already_checked, int digits[],
                         int size)
{
  int next = 0;
  int entry = 0;
  for (int i = 0; i < size; i++) {
    bool RepCheck = RepetitionCheck(already_checked, digits[i]);
    if (RepCheck == true) {
      continue;
    } else {
      next = digits[i];
      entry = i;
      break;
    }
  }
  return std::make_pair(next, entry);
}

double kloe_simu::AttenuationFactor(double d, int planeID)
{
  /*
       dE/dx attenuation - Ea=p1*exp(-d/atl1)+(1.-p1)*exp(-d/atl2)
         d    distance from photocatode - 2 cells/cell; d1 and d2
        atl1  50. cm
        atl2  430 cm planes 1-2    innermost plane is 1
              380 cm plane 3
              330 cm planes 4-5
         p1   0.35
  */
  double atl2 = 0.0;

  switch (planeID) {
    case 0:
    case 1:
      atl2 = kloe_simu::atl2_01;
      break;

    case 2:
      atl2 = kloe_simu::atl2_2;
      break;

    case 3:
    case 4:
      atl2 = kloe_simu::atl2_34;
      break;

    default:
      // std::cout << "planeID out if range" << std::endl;
      atl2 = -999.0;
      break;
  }

  return kloe_simu::p1 * TMath::Exp(-d / kloe_simu::atl1) +
         (1. - kloe_simu::p1) * TMath::Exp(-d / atl2);
}

// reconstruct t of the hit from tdc1 and tdc2
double kloe_simu::TfromTDC(double t1, double t2, double L)
{
  return 0.5 * (t1 + t2 - kloe_simu::vlfb * L / kloe_simu::m_to_mm);
}

// energy deposit of the hit from adc1 and adc2 and
// reconstructed longidutinal coordinate
double kloe_simu::EfromADC(double adc1, double adc2, double d1, double d2,
                           int planeID)
{
  double f1 = AttenuationFactor(d1, planeID);
  double f2 = AttenuationFactor(d2, planeID);

  double const attpassratio = 0.187;
  return 0.5 * (adc1 / f1 + adc2 / f2) /
         (attpassratio * kloe_simu::pe2ADC * kloe_simu::e2p2);
}

double EfromADCsingle(double adc, double f)
{
  double const attpassratio = 0.187;
  return adc / (f * attpassratio * kloe_simu::pe2ADC * kloe_simu::e2p2);
}

double DfromADC(double ta, double tb)
{
  return 0.5 * (ta - tb) / kloe_simu::vlfb * kloe_simu::m_to_mm;
}
