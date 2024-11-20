#include "SANDTrackerDriftCellMap.h"

#include <iostream>

// Default constructor
SANDTrackerDriftCellMap::SANDTrackerDriftCellMap() {}
// constructor with direct initialization of the parameters from the map path
SANDTrackerDriftCellMap::SANDTrackerDriftCellMap(std::string cell_map_path)
{
  // Initialize the configuration parameter attributes
  set_config_pars(cell_map_path);
  // initialize the pixelID map
  initialize_wf_map(cell_map_path);
}

// parameter and map initialization
void SANDTrackerDriftCellMap::init(std::string cell_map_path)
{
  // Initialize the configuration parameter attributes
  set_config_pars(cell_map_path);
  // initialize the pixelID map
  initialize_wf_map(cell_map_path);
}

int SANDTrackerDriftCellMap::encode_pix_id(const double &x,
                                           const double &y) const
{
  // y=0 is  the centre of the cell -> shift pixels up by a half thickness
  return std::round(x * 1e2) * 1e3 + std::round((y + cell_size_.at(1)) * 1e2);
}

void SANDTrackerDriftCellMap::set_config_pars(std::string c_map_path)
{
  std::unique_ptr<TFile> data_file(TFile::Open(
      Form("%s/track_drift_slice_ix0_iy0_.root", c_map_path.c_str())));
  std::unique_ptr<std::map<std::string, double>> c_map(
      data_file->Get<std::map<std::string, double>>("config_pars"));
  // initialize the cell size parameters
  cell_size_ = {c_map->at("c_hw"), c_map->at("c_ht")};
  sense_coords_ = {c_map->at("c_hw"), 0., 0.};
  time_window_ = {0., 1., c_map->at("n_t_bins")};
  cell_volt_ = {c_map->at("v_strip"), c_map->at("v_sens"),
                c_map->at("v_field")};
  d_sense_ = c_map->at("d_sense");
}

void SANDTrackerDriftCellMap::initialize_wf_map(std::string cell_map_path)
{

  auto map_rdf =
      ROOT::RDataFrame("out_tree", Form("%s/*.root", cell_map_path.c_str()));
  auto make_pair_lambda = [](int pixID, std::vector<double> wf) {
    return std::make_pair(pixID, wf);
  };
  auto encode_lambda = [this](const double &x0, const double &y0) {
    return this->encode_pix_id(x0, y0);
  };
  auto ID_wf_pair_lst =
      map_rdf.Filter("drift_st==-5")
          .Define("pixID", encode_lambda, {"x0", "y0"})
          .Define("ID_wf_pair", make_pair_lambda, {"pixID", "induced_wf"})
          .Take<std::pair<int, std::vector<double>>>("ID_wf_pair")
          .GetValue();
  // fill induced_wf_map with the pair vector entries
  for (const auto &pair_entry : ID_wf_pair_lst) {
    induced_wf_map_[pair_entry.first] = pair_entry.second;
  }
}

const std::vector<double> SANDTrackerDriftCellMap::get_cluster_waveform(
    const double &c_x, const double &c_y /*, const double& c_z */) const
{
  const auto c_ID = encode_pix_id(c_x, c_y);
  if (induced_wf_map_.count(c_ID) > 0)
    return induced_wf_map_.at(c_ID);
  else {
    // std::cout << "Found no map entry!\n";
    return std::vector<double>();
  }
}

std::vector<std::vector<double>>
    SANDTrackerDriftCellMap::generate_uniform_clusters(
        const double &hit_de, const std::array<double, 3> &x1,
        const std::array<double, 3> &x2, const double &t_0) const
{
  // initialize a track object
  std::vector<std::vector<double>> cluster_coords = {};
  // estimate the number of electron/ion pairs for hit_de
  const int n_clusters = (hit_de / cell_W_factor_);
  // compute the hit segment track length
  const double track_path_len = std::sqrt(std::pow(x1.at(0) - x2.at(0), 2) +
                                          std::pow(x1.at(1) - x2.at(1), 2) +
                                          std::pow(x1.at(2) - x2.at(2), 2));

  // distribute the electrons in single-clusters uniformly across the hit
  // segment
  for (int idx = 0; idx < n_clusters; idx++) {
    cluster_coords.push_back(std::vector<double>{
        x1.at(0) + (1.0 * idx / (n_clusters - 1)) * (x2.at(0) - x1.at(0)),
        x1.at(1) + (1.0 * idx / (n_clusters - 1)) * (x2.at(1) - x1.at(1)),
        x1.at(2) + (1.0 * idx / (n_clusters - 1)) * (x2.at(2) - x1.at(2)), t_0,
        1., 1.});
  }
  return cluster_coords;
}

std::vector<double> SANDTrackerDriftCellMap::build_induced_waveform(
    const std::array<double, 3> &hit_loc_start,
    const std::array<double, 3> &hit_loc_stop, const double &hit_de,
    const double &hit_t0) const
{
  // std::cout << "hit_x1: (" << hit_loc_start[0] << ", " << hit_loc_start[1]
  //           << ", " << hit_loc_start[2] << ")\n";
  // std::cout << "hit_x2: (" << hit_loc_stop[0] << ", " << hit_loc_stop[1] << ", "
  //           << hit_loc_stop[2] << ")\n";
  std::vector<double> track_waveform(time_window_.at(2));
  // loop over the cluster coordinate vectors (x,y,z,t,ne,ni)
  std::vector<std::vector<double>> avg_vecs;
  // std::cout << "Here 1\n";
  // generate uniformly distributed clusters along the segment from hit_de
  auto c_vec =
      generate_uniform_clusters(hit_de, hit_loc_start, hit_loc_stop, hit_t0);
  // std::cout << "Here 2\n";
  // find the waveforms for the active clusters and scale them by the number of
  // electrons
  for (const auto &cluster : c_vec) {
    if (std::abs(cluster.at(0) - sense_coords_.at(0)) > cell_size_.at(0) ||
        std::abs(cluster.at(1) - sense_coords_.at(1)) > cell_size_.at(1))
      continue;

    // std::cout << "cluster: (" << cluster.at(0) << ", " << cluster.at(1)<<")\n";

    std::vector<double> cluster_wf =
        get_cluster_waveform(cluster.at(0), cluster.at(1));

    if (!cluster_wf.size()) continue;

    // int cluster_t0 = cluster.at(3); // time shift for the single cluster is
    // not yet implemented
    auto cluster_ne = cluster.at(4);

    // NOTICE: the induction signals are directly set to positive!!
    std::transform(cluster_wf.begin(), cluster_wf.end(), cluster_wf.begin(),
                   [&cluster_ne](double &c) { return -1 * c * cluster_ne; });

    avg_vecs.push_back(cluster_wf);
  }
  // std::cout << "Here 3\n";
  // if (!avg_vecs.size())
  //   std::cout << ">Empty waveform\n";
  // else
  //   std::cout << ">Viable clusters: " << avg_vecs.size() << "\n";
  // std::cout << "\n";
  // add the cluster signals to the overall vector

  for (std::size_t idx = 0; idx < avg_vecs.size(); idx++) {
    std::transform(track_waveform.begin(), track_waveform.end(),
                   avg_vecs.at(idx).begin(), track_waveform.begin(),
                   std::plus<double>());
  }
  // std::cout << "Here 4\n";
  // shift the waveform by the hit t0 (+ propagation time)
  std::move(track_waveform.begin(), track_waveform.end() - hit_t0,
            track_waveform.begin() + hit_t0);
  std::fill(track_waveform.begin(), track_waveform.begin() + hit_t0, 0.);

  // for (int i = 0; i < 500; i++) std::cout << track_waveform[i] << ", ";
  // std::cout << "\n";
  return track_waveform;
}