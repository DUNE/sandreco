#include "SANDClustering.h"
#include "utils.h"

std::vector<cluster> Clusterize(std::vector<dg_cell>* vec_cellraw)
{
  std::vector<dg_cell> complete_cells, broken_cells, multicomplete_cells;
  std::vector<cluster> vec_clust;

  for (auto const& cell : *vec_cellraw) {
    if (cell.ps1.size() == 0 && cell.ps2.size() == 0) {

      continue;
    } else if (cell.ps1.size() == 0 || cell.ps2.size() == 0) {

      broken_cells.push_back(cell);
    } else if ((cell.ps1.size() != 0 && cell.ps2.size() != 0)) {
      // Complete cell

      complete_cells.push_back(cell);
    }
  }

  std::pair<std::vector<dg_cell>, std::vector<dg_cell>> processed_cells =
      ProcessMultiHits(complete_cells, broken_cells);
  multicomplete_cells = processed_cells.first;
  broken_cells = processed_cells.second;

  std::vector<int> checked_array;
  std::vector<int> vec_cell;

  std::vector<int> chck;

  for (int i = 0; i < multicomplete_cells.size(); i++) {

    std::vector<dg_cell> v_cell;

    if (RepetitionCheck(chck, i) == true) {
      continue;
    } else {
      chck.push_back(i);
    }

    v_cell.push_back(multicomplete_cells.at(i));

    std::pair<std::vector<dg_cell>, std::vector<int>> Neighbours =
        GetNeighbours(multicomplete_cells, i, chck, v_cell);
    v_cell = Neighbours.first;
    chck = Neighbours.second;
    struct cluster Clust;

    Clust = Calc_variables(v_cell);
    vec_clust.push_back(Clust);
  }
  // SPLIT

  vec_clust = Split(vec_clust);
  
  // MERGE

  vec_clust = Merge(vec_clust);
  
  // Track Fit
  vec_clust = TrackFit(vec_clust);
  
  // Recover Incomplete cells
  vec_clust = RecoverIncomplete(vec_clust, broken_cells);
  return vec_clust;
}

std::pair<std::vector<dg_cell>, std::vector<dg_cell>> ProcessMultiHits(
    std::vector<dg_cell> og_cell, std::vector<dg_cell> incomplete_cells)
{

  std::vector<dg_cell> complete_cells;
  for (auto const& cell : og_cell) {
    double delta = cell.l * sand_reco::ecal::scintillation::vlfb /
                   sand_reco::conversion::m_to_mm;

    for (int i = 0; i < cell.ps1.size(); i++) {
      int found = 0;
      for (int j = 0; j < cell.ps2.size(); j++) {
        if (abs(cell.ps1.at(i).tdc - cell.ps2.at(j).tdc) < delta) {
          dg_cell good_cell;
          good_cell.id = cell.id;
          good_cell.z = cell.z;
          good_cell.x = cell.x;
          good_cell.y = cell.y;
          good_cell.l = cell.l;
          good_cell.mod = cell.mod;
          good_cell.lay = cell.lay;
          good_cell.cel = cell.cel;
          good_cell.ps1.push_back(cell.ps1.at(i));
          good_cell.ps2.push_back(cell.ps2.at(j));
          complete_cells.push_back(good_cell);

          found++;
          break;
        }
      }
      if (found == 0) {

        dg_cell ps1bad_cell;
        ps1bad_cell.id = cell.id;
        ps1bad_cell.z = cell.z;
        ps1bad_cell.x = cell.x;
        ps1bad_cell.y = cell.y;
        ps1bad_cell.l = cell.l;
        ps1bad_cell.mod = cell.mod;
        ps1bad_cell.lay = cell.lay;
        ps1bad_cell.cel = cell.cel;
        ps1bad_cell.ps1.push_back(cell.ps1.at(i));
        incomplete_cells.push_back(ps1bad_cell);
      }
    }

    for (int k = 0; k < cell.ps2.size(); k++) {
      int found = 0;
      for (int l = 0; l < cell.ps1.size(); l++) {
        if (abs(cell.ps1.at(l).tdc - cell.ps2.at(k).tdc) < delta) {
          found++;
        }
      }
      if (found == 0) {

        dg_cell ps2bad_cell;
        ps2bad_cell.id = cell.id;
        ps2bad_cell.z = cell.z;
        ps2bad_cell.x = cell.x;
        ps2bad_cell.y = cell.y;
        ps2bad_cell.l = cell.l;
        ps2bad_cell.mod = cell.mod;
        ps2bad_cell.lay = cell.lay;
        ps2bad_cell.cel = cell.cel;
        ps2bad_cell.ps2.push_back(cell.ps2.at(k));
        incomplete_cells.push_back(ps2bad_cell);
      }
    }
  }
  return std::make_pair(complete_cells, incomplete_cells);
}

std::vector<cluster> RecoverIncomplete(std::vector<cluster> clus,
                                       std::vector<dg_cell> incomplete_cells)
{

  std::vector<double> tAvec;
  std::vector<double> tBvec;
  /*for (auto cl : clus) {
    int splitted = 0;
    double tA = 0, tB = 0, tA2 = 0, tB2 = 0;
    double EA, EB, EAtot = 0, EBtot = 0, EA2tot = 0, EB2tot = 0;
    double tRMS_A, tRMS_B, dist;
    std::vector<dg_cell> all_cells = cl.cells;
    for (int q = 0; q < all_cells.size(); q++) {
      EA = all_cells.at(q).ps1.at(0).adc;
      EB = all_cells.at(q).ps2.at(0).adc;
      EAtot += EA;
      EA2tot += EA * EA;
      EBtot += EB;
      EB2tot += EB * EB;
      tA += (all_cells.at(q).ps1.at(0).tdc) * EA;
      tA2 += std::pow(all_cells.at(q).ps1.at(0).tdc, 2) * EA;
      tB += (all_cells.at(q).ps2.at(0).tdc) * EB;
      tB2 += std::pow(all_cells.at(q).ps2.at(0).tdc, 2) * EB;
    }
    tA = tA / EAtot;
    tAvec.push_back(tA);
    tA2 = tA2 / EAtot;
    tB = tB / EBtot;
    tBvec.push_back(tB);
    tB2 = tB2 / EBtot;
  }*/

  for (auto const& brok_cells : incomplete_cells) {
    int isbarrel = 0;
    if (brok_cells.mod == 30) isbarrel = 1;
    if (brok_cells.mod == 40) isbarrel = 2;
    double cell_phi =
        atan((brok_cells.z - 23910.00) / (brok_cells.y + 2384.73)) * 180 /
        TMath::Pi();
    double cell_theta =
        atan((brok_cells.z - 23910.00) / (brok_cells.x)) * 180 / TMath::Pi();
    int minentry = 0;
    int found = 0;

    std::vector<int> clus_index;
    for (int j = 0; j < clus.size(); j++) {
      double rec_en = 0;
      bool hasNeigh = false;
      int nclusterC = 0;

      for (int i = 0; i < clus.at(j).cells.size(); i++) {

        if (isNeighbour(brok_cells.id, clus.at(j).cells.at(i).id)) {

          hasNeigh = true;
          nclusterC++;

          if (find(clus_index.begin(), clus_index.end(), j) ==
              clus_index.end()) {
            clus_index.push_back(j);
          }
        }
      }

      if (hasNeigh) {

        minentry = j;
        found = 1;

      } else if (found == 0) {
        double clus_phi =
            atan((clus.at(j).z - 23910.00) / (clus.at(j).y + 2384.73)) * 180 /
            TMath::Pi();
        // if (isbarrel != 0) cout << "Cluster Phi: " << clus_phi << std::endl;
        double clus_theta = atan((clus.at(j).z - 23910.00) / (clus.at(j).x)) *
                            180 / TMath::Pi();
        // if (isbarrel == 0) cout << "Cluster Theta: " << clus_theta <<
        // std::endl;
        double minphi = 999, mintheta = 999;
        int isbarrelc = 0;
        if (clus.at(j).cells[0].mod == 30) isbarrelc = 1;  
        if (clus.at(j).cells[0].mod == 40) isbarrelc = 2;
        if (abs(cell_phi - clus_phi) < 3 && (isbarrelc == isbarrel) &&
            isbarrelc == 0) {
          found = 1;
          double dist = sqrt(
              (brok_cells.z - clus.at(j).z) * (brok_cells.z - clus.at(j).z) +
              (brok_cells.y - clus.at(j).y) * (brok_cells.y - clus.at(j).y));
          if (abs(cell_phi - clus_phi) < minphi && dist < 2000) {
            // cout << "Barrel found 1: Clus Theta/Phi " << clus_theta << " / "
            // << clus_phi << " VS Cell Theta/Phi: " << cell_theta << " " <<
            // cell_phi
            // << std::endl; cout << " Cluster: " << clus.at(j).x << " " <<
            // clus.at(j).y
            // << " " << clus.at(j).z << " VS Cell: " << brok_cells.x << " " <<
            // brok_cells.y << " " << brok_cells.z << "  ->  Clus-cell: " <<
            // clus.at(j).cells.at(0).id << " and cell id: " << brok_cells.id <<
            // "
            // -> " << dist << std::endl;
            minphi = abs(cell_phi - clus_phi);
            minentry = j;
          }
          continue;
        }
        if (isbarrel == isbarrelc && isbarrel != 0) {
          if (abs(cell_theta - clus_theta) < 3) {
            // cout << "EndCap found 1: Clus Theta/Phi " << clus_theta << " / "
            // << clus_phi << " VS Cell Theta/Phi: " << cell_theta << " " <<
            // cell_phi
            // << std::endl; cout << " Cluster: " << clus.at(j).x << " " <<
            // clus.at(j).y
            // << " " << clus.at(j).z << " VS Cell: " << brok_cells.x << " " <<
            // brok_cells.y << " " << brok_cells.z << "  ->  Clus-cell: " <<
            // clus.at(j).cells.at(0).id << "-" << clus.at(j).cells.at(0).mod <<
            // " and cell id: " << brok_cells.id << std::endl;
            found = 1;
            mintheta = abs(cell_theta - clus_theta);
            minentry = j;
          }
          continue;
        }
      }
    }
    if (found == 1 && isbarrel == 0) {
      double rec_en = 0;
      double DpmA = clus.at(minentry).x / 10 + 215; 
      double DpmB = -clus.at(minentry).x / 10 + 215; 
      if (brok_cells.ps1.size() != 0 && brok_cells.ps2.size() != 0) {
        double Ea = brok_cells.ps1.at(0).adc;
        double Eb = brok_cells.ps2.at(0).adc;
        rec_en =
            sand_reco::ecal::reco::EfromADC(Ea, Eb, DpmA, DpmB, brok_cells.lay);
        clus.at(minentry).e = clus.at(minentry).e + rec_en;
        clus.at(minentry).cells.push_back(brok_cells);
        
      } else if (brok_cells.ps1.size() != 0) {
        int laycell = brok_cells.lay;
        double f =
            sand_reco::ecal::attenuation::AttenuationFactor(DpmA, laycell);
        rec_en = EfromADCsingle(brok_cells.ps1.at(0).adc, f);
        clus.at(minentry).e = clus.at(minentry).e + rec_en; 
        clus.at(minentry).cells.push_back(brok_cells);
      } else if (brok_cells.ps2.size() != 0) {
        int laycell = brok_cells.lay;
        double f =
            sand_reco::ecal::attenuation::AttenuationFactor(DpmB, laycell);
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
        rec_en =
            sand_reco::ecal::reco::EfromADC(Ea, Eb, DpmA, DpmB, brok_cells.lay);
        clus.at(minentry).e = clus.at(minentry).e + rec_en;
        clus.at(minentry).cells.push_back(brok_cells);
      } else if (brok_cells.ps1.size() != 0) {
        int laycell = brok_cells.lay;
        double f =
            sand_reco::ecal::attenuation::AttenuationFactor(DpmA, laycell);
        rec_en = EfromADCsingle(brok_cells.ps1.at(0).adc, f);
        clus.at(minentry).e = clus.at(minentry).e + rec_en;
        clus.at(minentry).cells.push_back(brok_cells);
      } else if (brok_cells.ps2.size() != 0) {
        int laycell = brok_cells.lay;
        double f =
            sand_reco::ecal::attenuation::AttenuationFactor(DpmB, laycell);
        rec_en = EfromADCsingle(brok_cells.ps2.at(0).adc, f);
        clus.at(minentry).e = clus.at(minentry).e + rec_en;
        clus.at(minentry).cells.push_back(brok_cells);
      }
    }
  }

  return clus;
}

std::vector<cluster> Split(std::vector<cluster> original_clu_vec)
{
  std::vector<cluster> clu_vec;
  std::vector<dg_cell> all_cells;
  int num_c = 0;
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
          DfromTDC(all_cells[j].ps1.at(0).tdc, all_cells[j].ps2.at(0).tdc);
      double d1, d2, d3;
      d1 = 0.5 * all_cells[j].l + d;
      d2 = 0.5 * all_cells[j].l - d;
      double cell_E = sand_reco::ecal::reco::EfromADC(
          all_cells[j].ps1.at(0).adc, all_cells[j].ps2.at(0).adc, d1, d2,
          all_cells[j].lay);
      tA += (all_cells.at(j).ps1.at(0).tdc -
             sand_reco::ecal::scintillation::vlfb * d1 /
                 sand_reco::conversion::m_to_mm) *
            EA;
      tA2 += std::pow(all_cells.at(j).ps1.at(0).tdc -
                          sand_reco::ecal::scintillation::vlfb * d1 /
                              sand_reco::conversion::m_to_mm,
                      2) *
             EA;
      tB += (all_cells.at(j).ps2.at(0).tdc -
             sand_reco::ecal::scintillation::vlfb * d2 /
                 sand_reco::conversion::m_to_mm) *
            EB;
      tB2 += std::pow(all_cells.at(j).ps2.at(0).tdc -
                          sand_reco::ecal::scintillation::vlfb * d2 /
                              sand_reco::conversion::m_to_mm,
                      2) *
             EB;
    }

    tA = tA / EAtot;

    tA2 = tA2 / EAtot;
    tB = tB / EBtot;

    tB2 = tB2 / EBtot;
    tRMS_A = sqrt(abs((tA2 - tA * tA) * (EA2tot - EAtot * EAtot) / EA2tot));

    tRMS_B = sqrt(abs((tB2 - tB * tB) * (EB2tot - EBtot * EBtot) / EB2tot));

    dist = std::sqrt(tRMS_A * tRMS_A + tRMS_B * tRMS_B);

    if (dist > 5) {

      std::vector<dg_cell> q1_cells, q2_cells, q3_cells, q4_cells;
      int n_cella = 0;
      for (auto const& a_cells : all_cells) {

        double d = DfromTDC(a_cells.ps1.at(0).tdc, a_cells.ps2.at(0).tdc);
        double d1, d2, d3;
        d1 = 0.5 * a_cells.l + d;
        d2 = 0.5 * a_cells.l - d;

        double t_difA = a_cells.ps1.at(0).tdc -
                        (sand_reco::ecal::scintillation::vlfb * d1 /
                         sand_reco::conversion::m_to_mm) -
                        tA;
        double t_difB = a_cells.ps2.at(0).tdc -
                        (sand_reco::ecal::scintillation::vlfb * d2 /
                         sand_reco::conversion::m_to_mm) -
                        tB;

        double t_difA_old = a_cells.ps1.at(0).tdc - tA;
        double t_difB_old = a_cells.ps2.at(0).tdc - tB;

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
        n_cella++;
      }
      std::vector<cluster> quadrant_cluster;
      if (q1_cells.size() != 0) {
        cluster clus = Calc_variables(q1_cells);
        quadrant_cluster.push_back(clus);

        splitted++;
      }
      if (q2_cells.size() != 0) {
        cluster clus = Calc_variables(q2_cells);
        quadrant_cluster.push_back(clus);

        splitted++;
      }
      if (q3_cells.size() != 0) {
        cluster clus = Calc_variables(q3_cells);
        quadrant_cluster.push_back(clus);

        splitted++;
      }
      if (q4_cells.size() != 0) {
        cluster clus = Calc_variables(q4_cells);
        quadrant_cluster.push_back(clus);

        splitted++;
      }
      q1_cells.clear();
      q2_cells.clear();
      q3_cells.clear();
      q4_cells.clear();
      auto new_vec_clu = Split(quadrant_cluster);
      for (auto cl : new_vec_clu) {
        clu_vec.push_back(cl);
      }
      quadrant_cluster.clear();
      new_vec_clu.clear();
    } else {
      clu_vec.push_back(clus);
    }
    all_cells.clear();
    num_c++;
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
    double varxi = Og_cluster.at(i).varx;
    double varyi = Og_cluster.at(i).vary;
    double varzi = Og_cluster.at(i).varz;
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
    clust.reco_cells = Og_cluster.at(i).reco_cells;

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
      double varxj = Og_cluster.at(j).varx;
      double varyj = Og_cluster.at(j).vary;
      double varzj = Og_cluster.at(j).varz;

      double D = sqrt((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj) +
                      (zi - zj) * (zi - zj));
      double DT = abs(ti - tj);
      double dx = sqrt(std::pow(xi - xj, 2));
      double dy = sqrt(std::pow(yi - yj, 2));
      double dz = sqrt(std::pow(zi - zj, 2));

      if (D < 400 && DT < 2.5) {

        bool endcap = false;
        if (clust.cells[0].id > 25000) {
          endcap = true;
        }
        if (endcap == true) {
          double Dz_ec = abs(yi - yj);
          D = sqrt((xi - xj) * (xi - xj) + (zi - zj) * (zi - zj));
          if (Dz_ec < 300 && D < 400) {
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
          if (Dz_bar < 300 && D < 400) {
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
    double apx[3] = {0, 0, 0}, eapx[3] = {0, 0, 0}, ctrk[3] = {0, 0, 0},
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
    double X[5] = {0, 0, 0, 0, 0}, D = 0;
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
      yy[lay_cross - 1] = Lay0.y;
      yz[lay_cross - 1] = Lay0.z;
      wz[lay_cross - 1] = 0.6;
      wy[lay_cross - 1] = 0.6;
      wx[lay_cross - 1] = 0.001 * Lay0.e;
      if (isBarrel == false) {
        wy[lay_cross - 1] = wx[lay_cross - 1];
        wx[lay_cross - 1] = 0.6;
      }
      // D = D + xl[0];
      // if (first_lay == 1) X[0] = 0.5 * xl[0];
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
        wy[lay_cross - 1] = wx[lay_cross - 1];
        wx[lay_cross - 1] = 0.6;
      }
      // D = D + xl[1];
      // if (first_lay == 2)
      //  X[0] = 0.5 * xl[1];
      // else if (first_lay == 1)
      //  X[1] = xl[0] + 0.5 * xl[1];
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
        wy[lay_cross - 1] = wx[lay_cross - 1];
        wx[lay_cross - 1] = 0.6;
      }
      // D = D + xl[2];
      // if (first_lay == 3)
      // X[0] = 0.5 * xl[2];
      // else if (first_lay == 2)
      //   X[1] = xl[1] + 0.5 * xl[2];
      // else if (first_lay == 1)
      //  X[2] = 2 * xl[0] + 0.5 * xl[2];
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
        wy[lay_cross - 1] = wx[lay_cross - 1];
        wx[lay_cross - 1] = 0.6;
      }
      // D = D + xl[3];
      // if (first_lay == 4)
      //  X[0] = 0.5 * xl[3];
      // else if (first_lay == 3)
      //  X[1] = xl[2] + 0.5 * xl[3];
      // else if (first_lay == 2)
      //  X[2] = 2 * xl[0] + 0.5 * xl[3];
      // else if (first_lay == 1)
      //  X[3] = 3 * xl[0] + 0.5 * xl[3];
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
        wy[lay_cross - 1] = wx[lay_cross - 1];
        wx[lay_cross - 1] = 0.6;
      }
      // D = D + xl[4];
      // if (first_lay == 5)
      //   X[0] = 0.5 * xl[4];
      // else if (first_lay == 4)
      //   X[1] = xl[3] + 0.5 * xl[4];
      // else if (first_lay == 3)
      //   X[2] = 2 * xl[3] + 0.5 * xl[4];
      // else if (first_lay == 2)
      //   X[3] = 3 * xl[3] + 0.5 * xl[4];
      // else if (first_lay == 1)
      //   X[4] = 4 * xl[3] + 0.5 * xl[4];
    }
    if (lay_cross == 0) {
      continue;
    }

    // D = D - 0.5 * xl[lay_cross];
    int Q = 0, L = 0;
    for (int k = first_lay; k <= 5; k++) {
      if (Q == 0 && (LayE[k - 1] >= 0.05 * clu_vec.at(i).e)) {  // 0.5?
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
      double R1 = 0;
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
      double XFix[5] = {0, 0, 0, 0, 0};
      if (first_lay == 1) {
        XFix[0] = 2.22;
        XFix[1] = 6.66;
        XFix[2] = 11.1;
        XFix[3] = 15.54;
        XFix[4] = 20.38;
      } else if (first_lay == 2) {
        XFix[0] = 6.66;
        XFix[1] = 11.1;
        XFix[2] = 15.54;
        XFix[3] = 20.38;
      } else if (first_lay == 3) {
        XFix[0] = 11.1;
        XFix[1] = 15.54;
        XFix[2] = 20.38;
      } else if (first_lay == 4) {
        XFix[0] = 15.54;
        XFix[1] = 20.38;
      } else if (first_lay == 5) {
        XFix[0] = 20.38;
      }
      for (int j = 0; j < lay_cross; j++) {
        X[j] = XFix[j] - Zapx;
      }

      std::tuple<double, double, double, double> fit_varx =
          fit_ls(lay_cross, X, yx, wx);

      std::tuple<double, double, double, double> fit_vary =
          fit_ls(lay_cross, X, yy, wy);

      std::tuple<double, double, double, double> fit_varz =
          fit_ls(lay_cross, X, yz, wz);
      double trktot = sqrt(std::get<1>(fit_varx) * std::get<1>(fit_varx) +
                           std::get<1>(fit_vary) * std::get<1>(fit_vary) +
                           std::get<1>(fit_varz) * std::get<1>(fit_varz));
      ctrk[0] = std::get<1>(fit_varx) / trktot;
      ctrk[1] = std::get<1>(fit_vary) / trktot;
      ctrk[2] = std::get<1>(fit_varz) / trktot;
      ectrk[0] = std::get<3>(fit_varx) / trktot;
      ectrk[1] = std::get<3>(fit_vary) / trktot;
      ectrk[2] = std::get<3>(fit_varz) / trktot;
      apx[0] = std::get<0>(fit_varx);
      apx[1] = std::get<0>(fit_vary);
      apx[2] = std::get<0>(fit_varz);
      eapx[0] = std::get<2>(fit_varx);
      eapx[1] = std::get<2>(fit_vary);
      eapx[2] = std::get<2>(fit_varz);
    }
    if (lay_cross == 1) {

      apx[0] = yx[0];
      apx[1] = yy[0];
      apx[2] = yz[0];
    }

    clu_vec.at(i).ax = apx[0];
    clu_vec.at(i).ay = apx[1];
    clu_vec.at(i).az = apx[2];

    clu_vec.at(i).sx = ctrk[0];
    clu_vec.at(i).sy = ctrk[1];
    clu_vec.at(i).sz = ctrk[2];
  }
  return clu_vec;
}

std::tuple<double, double, double, double> fit_ls(int lay, double* X, double* Y,
                                                  double* W)
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

  // std::vector<double> cells_e;
  std::vector<reco_cell> reconstructed_cells;
  
  for (auto& cell : cells) {
    
    reco_cell rec_cell;
    double d = DfromTDC(cell.ps1.at(0).tdc, cell.ps2.at(0).tdc);
    double d1, d2, d3;
    d1 = 0.5 * cell.l + d;

    d2 = 0.5 * cell.l - d;
    double cell_E = sand_reco::ecal::reco::EfromADC(
        cell.ps1.at(0).adc, cell.ps2.at(0).adc, d1, d2, cell.lay);
    // cells_e.push_back(cell_E);
    // cell.e = cell_E;

    double cell_T = sand_reco::ecal::reco::TfromTDC(cell.ps1.at(0).tdc,
                                                    cell.ps2.at(0).tdc, cell.l);

    rec_cell.id = cell.id;
    rec_cell.l = cell.l;
    rec_cell.mod = cell.mod;
    rec_cell.lay = cell.lay;
    rec_cell.e = cell_E;
    rec_cell.t = cell_T;
    

    if (cell.mod > 25) {  // endcap
      d3 = cell.y - d;

      rec_cell.y = d3;
      rec_cell.x = cell.x;
      rec_cell.z = cell.z;

      y_weighted = y_weighted + (d3 * cell_E);
      y2_weighted = y2_weighted + (d3 * d3 * cell_E);
      x_weighted = x_weighted + (cell.x * cell_E);
      x2_weighted = x2_weighted + (cell.x * cell.x * cell_E);
    } else {  // barrel
      d3 = cell.x - d;

      rec_cell.x = d3;
      rec_cell.y = cell.y;
      rec_cell.z = cell.z;

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

    reconstructed_cells.push_back(rec_cell);
    
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
  clust.reco_cells = reconstructed_cells;
  
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

bool isNeighbour(int id, int c_id)
{
  // Middle of the module
  if (c_id == id + 101 || c_id == id + 100 || c_id == id + 99 ||
      c_id == id + 1 || c_id == id - 1 || c_id == id - 99 || c_id == id - 100 ||
      c_id == id - 101) {
    return true;
  }
  // Right edge of the module
  else if (c_id == id - 889 || c_id == id - 989 || c_id == id - 1089) {
    return true;
  }
  // Right edge of 0 module
  else if ((c_id == id + 23111 || c_id == id + 23011 || c_id == id + 22911) &&
           id < 25000) {
    return true;
  }
  // Left edge of the module
  else if (c_id == id + 1089 || c_id == id + 989 || c_id == id + 889) {
    return true;
  }
  // Left edge of the 23 module
  else if ((c_id == id - 23011 || c_id == id - 22911 || c_id == id - 23111) &&
           id < 25000) {
    return true;
  }
  // Multiple hit on same cell
  else if (c_id == id) {
    return true;
  } else
    return false;
}

std::pair<std::vector<dg_cell>, std::vector<int>> GetNeighbours(
    std::vector<dg_cell> cells, int start, std::vector<int> checked,
    std::vector<dg_cell> neigh_chain)
{
  for (int i = 0; i < cells.size(); i++) {

    if (RepetitionCheck(checked, i) == true) continue;
    bool check = isNeighbour(cells.at(start).id, cells.at(i).id);

    if (check == true)

      if (check == true) {
        neigh_chain.push_back(cells.at(i));
        checked.push_back(i);

        std::pair<std::vector<dg_cell>, std::vector<int>> find_chain =
            GetNeighbours(cells, i, checked, neigh_chain);
        neigh_chain = find_chain.first;
        checked = find_chain.second;
      }
  }
  return std::make_pair(neigh_chain, checked);
}

double AttenuationFactor(double d, int planeID)
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
      atl2 = sand_reco::ecal::attenuation::atl2_01;
      break;
    case 1:
      atl2 = sand_reco::ecal::attenuation::atl2_01;
      break;

    case 2:
      atl2 = sand_reco::ecal::attenuation::atl2_2;
      break;

    case 3:
      atl2 = sand_reco::ecal::attenuation::atl2_34;
      break;

    case 4:
      atl2 = sand_reco::ecal::attenuation::atl2_34;
      break;

    default:

      atl2 = -999.0;
      break;
  }

  return sand_reco::ecal::attenuation::p1 *
             TMath::Exp(-d / sand_reco::ecal::attenuation::atl1) +
         (1. - sand_reco::ecal::attenuation::p1) * TMath::Exp(-d / atl2);
}

double TfromTDC(double t1, double t2, double L)
{
  return 0.5 * (t1 + t2 -
                sand_reco::ecal::scintillation::vlfb * L /
                    sand_reco::conversion::m_to_mm);
}

double EfromADC(double adc1, double adc2, double d1, double d2, int planeID)
{
  double f1 = AttenuationFactor(d1, planeID);
  double f2 = AttenuationFactor(d2, planeID);

  double const attpassratio = 0.187;
  return 0.5 * (adc1 / f1 + adc2 / f2) /
         (attpassratio * sand_reco::ecal::acquisition::pe2ADC *
          sand_reco::ecal::photo_sensor::e2pe);
}

double EfromADCsingle(double adc, double f)
{
  double const attpassratio = 0.187;
  return adc / (f * attpassratio * sand_reco::ecal::acquisition::pe2ADC *
                sand_reco::ecal::photo_sensor::e2pe);
}

double DfromTDC(double ta, double tb)
{
  return 0.5 * (ta - tb) / sand_reco::ecal::scintillation::vlfb *
         sand_reco::conversion::m_to_mm;
}

bool endsWith(const std::string& fullString, const std::string& ending)
{
  if (fullString.length() >= ending.length()) {
    return (0 == fullString.compare(fullString.length() - ending.length(),
                                    ending.length(), ending));
  } else {
    return false;
  }
}
