#include "SANDClustering.h"
#include "utils.h"
#include <numeric>

int evaluateClusterType(const cluster& clust) {
  // type 1 barrel
  // type 2 endcap
  // type 3 mixed
  int type = 0;
  bool barrel = false;
  bool endcap = false;
  for (const auto& c:clust.reco_cells) {
    if (c.id < 20000000 && !endcap) {
        type += 2;
        endcap = true;
    }
    if (c.id > 20000000 && !barrel) {
        type += 1;
        barrel = true;
    }
  }
  return type;
}

std::vector<cluster> clusterize(const SANDGeoManager* sand_geo, const std::vector<dg_cell>& cells)
{
  
  std::pair<std::vector<dg_cell>, std::vector<dg_cell>> processed_cells = processMultiHits(sand_geo, cells);
  std::vector<dg_cell> complete_cells   = processed_cells.first;
  std::vector<dg_cell> incomplete_cells = processed_cells.second;

  // std::cout << complete_cells.size() << " " << incomplete_cells.size() << std::endl;

  std::vector<int> chck;
  std::vector<cluster> vec_clust;
  for (uint i = 0; i < complete_cells.size(); i++) {

    std::vector<dg_cell> v_cell;

    if (repetitionCheck(chck, i) == true) {
      continue;
    } else {
      chck.push_back(i);
    }

    v_cell.push_back(complete_cells.at(i));

    std::pair<std::vector<dg_cell>, std::vector<int>> Neighbours =
        getNeighbours(complete_cells, i, chck, v_cell);
    v_cell = Neighbours.first;
    chck = Neighbours.second;

    // std::cout << v_cell.size() << std::endl;
    struct cluster Clust;

    Clust = createCluster(sand_geo, v_cell);

    if (!Clust.reco_cells.empty()) {
      vec_clust.push_back(Clust);
    }
  }

  // SPLIT
  int n_clu = 0;

  bool HasSplit = false;
  int iteration = 0;
  do {

    vec_clust = split(sand_geo, vec_clust, HasSplit);
    iteration++;
  } while (HasSplit);
    
  // MERGE
  vec_clust = merge(vec_clust);

  // Track Fit
  trackFit(vec_clust);

  // RecoverIncomplete(sand_geo, vec_clust, incomplete_cells);

  for (auto& c:vec_clust) {
    std::cout << c.e << std::endl;
    double total_pe = 0;
    for (auto cell:c.reco_cells) {
      total_pe += cell.ps1.adc + cell.ps2.adc;
    }
    std::cout << total_pe / 4 << std::endl;
    c.type = evaluateClusterType(c);
  }

  return vec_clust;
}

void clustInfo(cluster clus)
{
  std::cout << "Cluster Energy " << clus.e << " MeV" << std::endl;
  std::cout << "Centroid coordinates: " << clus.x << " [X] " << clus.y
            << " [Y] " << clus.z << " [z] , and mean arrival time: " << clus.t
            << " ns." << std::endl;
  std::cout << "Variance: " << clus.varx << " [X] " << clus.vary << " [Y] "
            << clus.varz << " [z]" << std::endl;
  std::cout << "Composed by the following cells: ";
  for (uint i = 0; i < clus.reco_cells.size(); i++) {
    std::cout << "Cell: " << clus.reco_cells.at(i).id
              << " X: " << clus.reco_cells.at(i).x
              << " Y: " << clus.reco_cells.at(i).y
              << " Z: " << clus.reco_cells.at(i).z
              << ", Energy: " << clus.reco_cells.at(i).e << std::endl;
  }
}

std::pair<std::vector<dg_cell>, std::vector<dg_cell>> processMultiHits(const SANDGeoManager* sand_geo,
    const std::vector<dg_cell>& cells)
{

  std::vector<dg_cell> complete_cells;
  std::vector<dg_cell> incomplete_cells;
  for (auto& cell : cells) {
    std::vector<int> found1;
    std::vector<int> found2;
    const auto& cell_info = sand_geo->get_ecal_cell_info(cell.id);
    double delta = cell_info.getLength() * sand_reco::ecal::scintillation::vlfb /
                   sand_reco::conversion::m_to_mm;

    for (uint i = 0; i < cell.ps1.size(); i++) {
      int found = 0;
      for (uint j = 0; j < cell.ps2.size(); j++) {
        if (fabs(cell.ps1.at(i).tdc - cell.ps2.at(j).tdc) < delta) {
          if (std::find(found1.begin(), found1.end(), i) == found1.end() &&
              std::find(found2.begin(), found2.end(), j) == found2.end()) {
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
            found1.push_back(i);
            found2.push_back(j);
            found++;
            break;
          }
        }
      }
      if (found == 0) {
        if (std::find(found1.begin(), found1.end(), i) == found1.end()) {
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
          found1.push_back(i);
        }
      }
    }

    for (uint k = 0; k < cell.ps2.size(); k++) {
      int found = 0;
      for (uint l = 0; l < cell.ps1.size(); l++) {
        if (fabs(cell.ps1.at(l).tdc - cell.ps2.at(k).tdc) < delta) {
          if (std::find(found1.begin(), found1.end(), l) == found1.end() &&
              std::find(found2.begin(), found2.end(), k) == found2.end()) {
              found++;
          }
        }
      }
      if (found == 0) {
        if (std::find(found2.begin(), found2.end(), k) == found2.end()) {
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
  }
  return std::make_pair(complete_cells, incomplete_cells);
}

void updateCluster(const dg_cell& incomplete_cell, double distance, 
                   int fired_pmt, sand_geometry::ecal::ECALCellInfo::ModuleType module_type, cluster& clus) {
  dg_ps ps1;
  dg_ps ps2;
  double adc;
  if (fired_pmt == 1) {
    ps1 = incomplete_cell.ps1.at(0);
    adc = incomplete_cell.ps1.at(0).adc;
    ps2.adc = 0;
  } else {
    ps2 = incomplete_cell.ps2.at(0);
    adc = incomplete_cell.ps2.at(0).adc;
    ps1.adc = 0;
  }
  
  int laycell = incomplete_cell.lay;
  double f = sand_reco::ecal::attenuation::AttenuationFactor(distance, laycell);
  double rec_en = sand_reco::ecal::reco::EfromADCsingle(adc, f);

  double x;
  double y;
  if (module_type == sand_geometry::ecal::ECALCellInfo::ModuleType::kBarrel) {
    x = -99999;
    y = incomplete_cell.y;
  } else {
    x = incomplete_cell.x;
    y = -99999;
  }

  reco_cell reco_cell_from_incomplete = {incomplete_cell.id, incomplete_cell.z, 
                                         y, x,
                                         incomplete_cell.l, incomplete_cell.mod, 
                                         laycell, rec_en, -99999, ps1, ps2, fired_pmt};

  clus.e = clus.e + rec_en;
  clus.reco_cells.push_back(reco_cell_from_incomplete);
}


void recoverIncomplete(const SANDGeoManager* sand_geo, std::vector<cluster>& clus,
                       const std::vector<dg_cell>& incomplete_cells)
{
  
  for (auto const& incomplete_cell : incomplete_cells) {
    auto module_type = sand_geo->get_ecal_cell_info(incomplete_cell.id).getModuleType();

    int closest_cluster_index;
    bool found = false;
    double spatial_range = 15;
    double time_range = 3;
    double min_distance = 10e5;
    double min_time = 10e5;
    for (uint j = 0; j < clus.size(); j++) {
      double rec_en = 0;
      bool isNeigh = false;
      
      double cell_time;
      if (incomplete_cell.ps1.size() != 0) {
        cell_time = incomplete_cell.ps1.at(0).tdc;
      } else {
        cell_time = incomplete_cell.ps2.at(0).tdc;
      }

      // Notice: one could compare the time of the closest 
      //         reco cell with the one of the incompete one
      double time_distance = fabs(clus.at(j).t - cell_time);
      
      double spatial_distance;
      if (module_type == sand_geometry::ecal::ECALCellInfo::ModuleType::kBarrel) {
        spatial_distance = sqrt(pow(clus[j].y - incomplete_cell.y, 2)+
                                pow(clus[j].z - incomplete_cell.z, 2));
      } else {
        spatial_distance = sqrt(pow(clus[j].x - incomplete_cell.x, 2)+
                                pow(clus[j].z - incomplete_cell.z, 2));
      }

      if (spatial_distance < spatial_range && time_distance < time_range) {
        // TODO: optimize the condition for the clostest cluster
        if (spatial_distance < min_distance && time_distance < min_time) {
          found = true;
          closest_cluster_index = j;
          min_distance = spatial_distance;
          min_time = time_distance;
        }
      }
    }

    reco_cell reco_cell_from_incomplete;
    if (found) {
      double DpmA;
      double DpmB;

      double inco_cell_lenght = sand_geo->get_ecal_cell_info(incomplete_cell.id).getLength();

      if (module_type == sand_geometry::ecal::ECALCellInfo::ModuleType::kBarrel) {
        DpmA =  clus.at(closest_cluster_index).x + inco_cell_lenght * 0.5;
        DpmB = -clus.at(closest_cluster_index).x + inco_cell_lenght * 0.5;
      } else {
        double shifted_y = clus.at(closest_cluster_index).y - (-2384.73);  // Shift relative to center
        DpmA =  shifted_y + inco_cell_lenght * 0.5;
        DpmB = -shifted_y + inco_cell_lenght * 0.5;
      }
      // std::cout << DpmA << " " << DpmB << std::endl;
      if (incomplete_cell.ps1.size() != 0) {
        updateCluster(incomplete_cell, DpmA, 1, module_type, clus.at(closest_cluster_index));
      } else {
        updateCluster(incomplete_cell, DpmB, 2, module_type, clus.at(closest_cluster_index));
      }
    } else {
        // TODO: think about it 
    }
    
  }
}

std::vector<cluster> split(const SANDGeoManager* sand_geo, const std::vector<cluster>& original_clu_vec,
                           bool& HasSplit)
{
  std::vector<cluster> clu_vec;
  std::vector<reco_cell> all_cells;

  HasSplit = false;
  for (auto const& clus : original_clu_vec) {

    int splitted = 0;
    double tA = 0, tB = 0, tA2 = 0, tB2 = 0;
    double EA, EB, EAtot = 0, EBtot = 0, EA2tot = 0, EB2tot = 0;
    double tRMS_A, tRMS_B, dist;

    all_cells = clus.reco_cells;
    for (uint j = 0; j < all_cells.size(); j++) {

      EA = all_cells.at(j).ps1.adc;

      EB = all_cells.at(j).ps2.adc;
      EAtot += EA;
      EA2tot += EA * EA;
      EBtot += EB;
      EB2tot += EB * EB;

      const auto& cell_info = sand_geo->get_ecal_cell_info(all_cells.at(j).id);
      double cell_lenght = cell_info.getLength();

      double d1, d2;
      d1 = sand_geo->compute_cell_d1(cell_lenght, all_cells.at(j).ps1.tdc, all_cells.at(j).ps2.tdc);
      d2 = sand_geo->compute_cell_d2(cell_lenght, all_cells.at(j).ps1.tdc, all_cells.at(j).ps2.tdc);

      double cell_E = sand_reco::ecal::reco::EfromADC(
          all_cells[j].ps1.adc, all_cells[j].ps2.adc, d1, d2, all_cells[j].lay);

      tA +=
          (all_cells.at(j).ps1.tdc - sand_reco::ecal::scintillation::vlfb * d1 /
                                         sand_reco::conversion::m_to_mm) *
          EA;

      tA2 += std::pow(all_cells.at(j).ps1.tdc -
                          sand_reco::ecal::scintillation::vlfb * d1 /
                              sand_reco::conversion::m_to_mm,
                      2) *
             EA;

      tB +=
          (all_cells.at(j).ps2.tdc - sand_reco::ecal::scintillation::vlfb * d2 /
                                         sand_reco::conversion::m_to_mm) *
          EB;

      tB2 += std::pow(all_cells.at(j).ps2.tdc -
                          sand_reco::ecal::scintillation::vlfb * d2 /
                              sand_reco::conversion::m_to_mm,
                      2) *
             EB;
    }
    tA = tA / EAtot;

    tA2 = tA2 / EAtot;
    tB = tB / EBtot;

    tB2 = tB2 / EBtot;
    tRMS_A = sqrt(fabs((tA2 - tA * tA) * (EA2tot - EAtot * EAtot) / EA2tot));

    tRMS_B = sqrt(fabs((tB2 - tB * tB) * (EB2tot - EBtot * EBtot) / EB2tot));

    dist = std::sqrt(tRMS_A * tRMS_A + tRMS_B * tRMS_B);

    if (dist > 5) {
      HasSplit = true;
      std::vector<reco_cell> q1_cells, q2_cells, q3_cells, q4_cells;

      int n_cella = 0;
      for (auto const& a_cells : all_cells) {

        const auto& cell_info = sand_geo->get_ecal_cell_info(a_cells.id);
        double cell_lenght = cell_info.getLength();
        double d1 = sand_geo->compute_cell_d1(cell_lenght, a_cells.ps1.tdc, a_cells.ps2.tdc);
        double d2 = sand_geo->compute_cell_d2(cell_lenght, a_cells.ps1.tdc, a_cells.ps2.tdc);

        double t_difA = a_cells.ps1.tdc -
                        (sand_reco::ecal::scintillation::vlfb * d1 /
                         sand_reco::conversion::m_to_mm) -
                        tA;

        double t_difB = a_cells.ps2.tdc -
                        (sand_reco::ecal::scintillation::vlfb * d2 /
                         sand_reco::conversion::m_to_mm) -
                        tB;

        double t_difA_old = a_cells.ps1.tdc - tA;

        double t_difB_old = a_cells.ps2.tdc - tB;

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
        cluster new_cluster = calcVariables(q1_cells);

        clu_vec.push_back(new_cluster);
        splitted++;
      }
      if (q2_cells.size() != 0) {
        cluster new_cluster = calcVariables(q2_cells);

        clu_vec.push_back(new_cluster);
        splitted++;
      }
      if (q3_cells.size() != 0) {
        cluster new_cluster = calcVariables(q3_cells);

        clu_vec.push_back(new_cluster);
        splitted++;
      }
      if (q4_cells.size() != 0) {
        cluster new_cluster = calcVariables(q4_cells);

        clu_vec.push_back(new_cluster);
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

  return clu_vec;
}

std::vector<cluster> merge(const std::vector<cluster>& Og_cluster)
{

  std::vector<cluster> mgd_cluster;
  std::vector<int> checked;

  for (uint i = 0; i < Og_cluster.size(); i++) {
    double xi = Og_cluster.at(i).x;
    double yi = Og_cluster.at(i).y;
    double zi = Og_cluster.at(i).z;
    double varxi = Og_cluster.at(i).varx;
    double varyi = Og_cluster.at(i).vary;
    double varzi = Og_cluster.at(i).varz;
    double ti = Og_cluster.at(i).t;
    double ei = Og_cluster.at(i).e;
    bool RepCheck = repetitionCheck(checked, i);
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
    clust.reco_cells = Og_cluster.at(i).reco_cells;

    for (uint j = i; j < Og_cluster.size(); j++) {
      RepCheck = repetitionCheck(checked, j);
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
      double DT = fabs(ti - tj);
      double dx = sqrt(std::pow(xi - xj, 2));
      double dy = sqrt(std::pow(yi - yj, 2));
      double dz = sqrt(std::pow(zi - zj, 2));

      if (D < 250 && DT < 2.5) {

        bool endcap = false;
        if (clust.reco_cells[0].id > 25000) {
          endcap = true;
        }
        if (endcap == true) {
          double Dz_ec = fabs(yi - yj);
          D = sqrt((xi - xj) * (xi - xj) + (zi - zj) * (zi - zj));
          if (Dz_ec < 250 && D < 250) {
            std::vector<reco_cell> vec_cells_j = Og_cluster.at(j).reco_cells;
            for (uint k = 0; k < vec_cells_j.size(); k++) {
              clust.reco_cells.push_back(vec_cells_j.at(k));
            }
            clust = calcVariables(clust.reco_cells);
            checked.push_back(j);
          }
        } else if (endcap == false) {
          double Dz_bar = fabs(xi - xj);
          D = sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj));
          if (Dz_bar < 250 && D < 250) {
            std::vector<reco_cell> vec_cells_j = Og_cluster.at(j).reco_cells;
            for (uint k = 0; k < vec_cells_j.size(); k++) {
              clust.reco_cells.push_back(vec_cells_j.at(k));
            }
            clust = calcVariables(clust.reco_cells);
            checked.push_back(j);
          }
        }
      }
    }
    mgd_cluster.push_back(clust);
  }

  return mgd_cluster;
}

void updateArrays(double* yx, double* yy, double* yz,
                  double* wx, double* wy, double* wz,
                  const cluster& lay, int lay_cross, bool isBarrel) 
{
  yx[lay_cross - 1] = lay.x;
  yy[lay_cross - 1] = lay.y;
  yz[lay_cross - 1] = lay.z;
  wx[lay_cross - 1] = 0.001 * lay.e;
  wy[lay_cross - 1] = 0.6;
  wz[lay_cross - 1] = 0.6;
  if (isBarrel == false) {
    wy[lay_cross - 1] = wx[lay_cross - 1];
    wx[lay_cross - 1] = 0.6;
  }
}

void trackFit(std::vector<cluster>& clu_vec)
{
  const double xl[5] = {4.44, 4.44, 4.44, 4.44, 5.24};
  for (uint i = 0; i < clu_vec.size(); i++) {
    double apx[3]   = {0, 0, 0};
    // double eapx[3]  = {0, 0, 0};
    double ctrk[3]  = {0, 0, 0};
    // double ectrk[3] = {0, 0, 0};
    std::vector<reco_cell> cell_vec[5];
        
    for (uint j = 0; j < clu_vec.at(i).reco_cells.size(); j++) {
      int layer_number = clu_vec.at(i).reco_cells.at(j).lay;
      cell_vec[layer_number].push_back(clu_vec.at(i).reco_cells.at(j));
    }

    cluster Lay[5];
    for (int j = 0; j < 5; j++) {
      Lay[j] = calcVariables(cell_vec[j]);
    }

    bool isBarrel = true;

    double yx[5] = {0, 0, 0, 0, 0};
    double yy[5] = {0, 0, 0, 0, 0};
    double yz[5] = {0, 0, 0, 0, 0};
    double wx[5] = {0, 0, 0, 0, 0};
    double wy[5] = {0, 0, 0, 0, 0};
    double wz[5] = {0, 0, 0, 0, 0};
    double X[5]  = {0, 0, 0, 0, 0};
    double D = 0;

    if (clu_vec.at(i).reco_cells[0].id > 25000) {
      isBarrel = false;
    }

    int lay_cross = 0, first_lay = 0;
    for (int layer_number = 0; layer_number < 5; layer_number++) {
      if (Lay[layer_number].e > 0) {
        lay_cross++;
        if (lay_cross == 1) {
          first_lay = layer_number + 1;
        }
        updateArrays(yx, yy, yz, wx, wy, wz, Lay[layer_number], lay_cross, isBarrel);
      }
    }
    
    if (lay_cross == 0) {
      continue;
    }

    int Q = 0, L = 0;
    double LayE[5] = {Lay[0].e, Lay[1].e, Lay[2].e, Lay[3].e, Lay[4].e};
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
      
      for (int j = 0; j < lay_cross; j++) {
        X[j] = sand_reco::ecal::geometry::XFix[first_lay + j - 1] - Zapx;
      }

      std::tuple<double, double, double, double> fit_varx =
          fitLs(lay_cross, X, yx, wx);

      std::tuple<double, double, double, double> fit_vary =
          fitLs(lay_cross, X, yy, wy);

      std::tuple<double, double, double, double> fit_varz =
          fitLs(lay_cross, X, yz, wz);
      double trktot = sqrt(std::get<1>(fit_varx) * std::get<1>(fit_varx) +
                           std::get<1>(fit_vary) * std::get<1>(fit_vary) +
                           std::get<1>(fit_varz) * std::get<1>(fit_varz));
      ctrk[0] = std::get<1>(fit_varx) / trktot;
      ctrk[1] = std::get<1>(fit_vary) / trktot;
      ctrk[2] = std::get<1>(fit_varz) / trktot;
      apx[0] = std::get<0>(fit_varx);
      apx[1] = std::get<0>(fit_vary);
      apx[2] = std::get<0>(fit_varz);
      // Notice: these are computed but not used. Why?
      // ectrk[0] = std::get<3>(fit_varx) / trktot;
      // ectrk[1] = std::get<3>(fit_vary) / trktot;
      // ectrk[2] = std::get<3>(fit_varz) / trktot;
      // eapx[0] = std::get<2>(fit_varx);
      // eapx[1] = std::get<2>(fit_vary);
      // eapx[2] = std::get<2>(fit_varz);
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
}

std::tuple<double, double, double, double> fitLs(int lay, double* X, double* Y,
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

cluster createCluster(const SANDGeoManager* sand_geo, const std::vector<dg_cell>& cells)
{

  double x_weighted = 0, y_weighted = 0, z_weighted = 0, t_weighted = 0,
         x2_weighted = 0, y2_weighted = 0, z2_weighted = 0, Etot = 0, E2tot = 0,
         EvEtot = 0, EA, EAtot = 0, EB, EBtot = 0, TA = 0, TB = 0;

  std::vector<reco_cell> reconstructed_cells;

  for (auto& cell : cells) {

    const auto& cell_length = sand_geo->get_ecal_cell_info(cell.id).getLength();

    reco_cell rec_cell;

    double d1, d2;
    d1 = sand_geo->compute_cell_d1(cell_length, cell.ps1.at(0).tdc, cell.ps2.at(0).tdc);
    d2 = sand_geo->compute_cell_d2(cell_length, cell.ps1.at(0).tdc, cell.ps2.at(0).tdc);
    // Notice: what do we want to do in these cases?
    if (d1 > cell_length || d1 < 0 || 
        d2 > cell_length || d2 < 0) {
      continue;
    }
    double cell_E = sand_reco::ecal::reco::EfromADC(
        cell.ps1.at(0).adc, cell.ps2.at(0).adc, d1, d2, cell.lay);

    double cell_T = sand_reco::ecal::reco::TfromTDC(cell.ps1.at(0).tdc,
                                                    cell.ps2.at(0).tdc, cell_length);

    rec_cell.id = cell.id;
    rec_cell.l = cell_length;
    rec_cell.mod = cell.mod;
    rec_cell.lay = cell.lay;
    rec_cell.e = cell_E;
    rec_cell.t = cell_T;
    rec_cell.ps1 = cell.ps1.at(0);
    rec_cell.ps2 = cell.ps2.at(0);
    rec_cell.fired_pmt = 3;

    double cell_x = -99999, cell_y = -99999, cell_z = -99999;
    sand_geo->get_reco_hit_pos(cell.id, cell_length, cell.ps1.at(0).tdc, cell.ps2.at(0).tdc, cell_x, cell_y, cell_z);
    
    if (cell_x == -99999 || cell_y == -99999 || cell_z == -99999) {
      continue;
    }
    rec_cell.x = cell_x;
    rec_cell.y = cell_y;
    rec_cell.z = cell_z;

    x_weighted = x_weighted + (rec_cell.x * cell_E);
    x2_weighted = x2_weighted + (rec_cell.x * rec_cell.x * cell_E);
    y_weighted = y_weighted + (rec_cell.y * cell_E);
    y2_weighted = y2_weighted + (rec_cell.y * rec_cell.y * cell_E);
  
    t_weighted = t_weighted + cell_T * cell_E;

    z_weighted = z_weighted + (rec_cell.z * cell_E);
    z2_weighted = z2_weighted + (rec_cell.z * rec_cell.z * cell_E);
    Etot = Etot + cell_E;
    E2tot = E2tot + cell_E * cell_E;

    reconstructed_cells.push_back(rec_cell);
  }

  if (reconstructed_cells.size() == 0) {
    return cluster();
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
  clust.reco_cells = reconstructed_cells;

  return clust;
}

cluster calcVariables(const std::vector<reco_cell>& cells)
{

  double x_weighted  = std::accumulate(cells.begin(), cells.end(), 0., [](double sum, const reco_cell& c) { return sum += c.x * c.e; });
  double x2_weighted = std::accumulate(cells.begin(), cells.end(), 0., [](double sum, const reco_cell& c) { return sum += c.x * c.x * c.e; });
  double y_weighted  = std::accumulate(cells.begin(), cells.end(), 0., [](double sum, const reco_cell& c) { return sum += c.y * c.e; });
  double y2_weighted = std::accumulate(cells.begin(), cells.end(), 0., [](double sum, const reco_cell& c) { return sum += c.y * c.y * c.e; });
  double z_weighted  = std::accumulate(cells.begin(), cells.end(), 0., [](double sum, const reco_cell& c) { return sum += c.z * c.e; });
  double z2_weighted = std::accumulate(cells.begin(), cells.end(), 0., [](double sum, const reco_cell& c) { return sum += c.z * c.z * c.e; });
  double t_weighted  = std::accumulate(cells.begin(), cells.end(), 0., [](double sum, const reco_cell& c) { return sum += c.t * c.e; });
  double Etot        = std::accumulate(cells.begin(), cells.end(), 0., [](double sum, const reco_cell& c) { return sum += c.e; });
  double E2tot       = std::accumulate(cells.begin(), cells.end(), 0., [](double sum, const reco_cell& c) { return sum += c.e * c.e; });

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
  clust.reco_cells = cells;

  return clust;
}

bool repetitionCheck(std::vector<int> v, int check)
{
  if (std::find(v.begin(), v.end(), check) != v.end()) {
    return true;
  } else {
    return false;
  }
}

bool isNeighbour(const dg_cell& cell, const dg_cell& check_cell)
{

  bool arebothendcap = false;
  bool arebothbarrel = false;
  
  if (fabs(cell.x) > 1500 && fabs(check_cell.x) > 1500) {
    arebothendcap = true;
  } else if (fabs(cell.x) < 1500 && fabs(check_cell.x) < 1500) {
    arebothbarrel = true;
  } else if ((fabs(cell.x) < 1500 && fabs(check_cell.x) > 1500) ||
             (fabs(cell.x) > 1500 && fabs(check_cell.x) < 1500)) {
    return false;
  }

  if (arebothendcap) {
    double distance = sqrt((cell.x - check_cell.x) * (cell.x - check_cell.x) +
                           (cell.z - check_cell.z) * (cell.z - check_cell.z));

    if (distance < 65.70) { //max distance between two cells in the endcap (layer 5-4, diagonal)
      return true;
    }

    else {
      return false;
    }
  } else if (arebothbarrel) {
    double distance = sqrt((cell.y - check_cell.y) * (cell.y - check_cell.y) +
                           (cell.z - check_cell.z) * (cell.z - check_cell.z));
    if (distance < 72.36) { //max distance between two cells in the barrel (layer5-4 digonal)
      return true;
    } else {
      return false;
    }
  }
  return false;
}

std::pair<std::vector<dg_cell>, std::vector<int>> getNeighbours(
    const std::vector<dg_cell>& cells, int start, std::vector<int> checked,
    std::vector<dg_cell> neigh_chain)
{
  for (uint i = 0; i < cells.size(); i++) {

    if (repetitionCheck(checked, i) == true) continue;

    bool check = isNeighbour(cells.at(start), cells.at(i));

      if (check == true) {
        neigh_chain.push_back(cells.at(i));
        checked.push_back(i);

        std::pair<std::vector<dg_cell>, std::vector<int>> find_chain =
            getNeighbours(cells, i, checked, neigh_chain);
        neigh_chain = find_chain.first;
        checked = find_chain.second;
      }
  }
  return std::make_pair(neigh_chain, checked);
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
