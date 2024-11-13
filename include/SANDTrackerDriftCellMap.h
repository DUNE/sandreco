#include <ROOT/RDataFrame.hxx>
#include <TObject.h>
#include <TVector2.h>

#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "struct.h"

#ifndef SANDTrackerDriftCellMap_H
#define SANDTrackerDriftCellMap_H

class SANDTrackerDriftCellMap : public TObject
{
 private:
  // cell parameters
  std::array<double, 2> cell_size_;     // cell sense wire coordinates
  std::array<double, 3> sense_coords_;  // sensor time window
  std::vector<double> time_window_;
  std::array<double, 3> cell_volt_;  // cell voltages
  double cell_W_factor_ = 30;  // gas mix. energy per ion pair factor W [eV]
  double d_sense_;             // cell sense diameter (for now)
  // ROOT::RDataFrame cell_map_;  // Map RDataframe
  std::map<int, std::vector<double>> induced_wf_map_;  // std::map of waveform
                                                       // IDs from the pixel
                                                       // coordinates

  // compute an index for the map pixels/ track clusters given the pixel spacing
  int encode_pix_id(const double &x, const double &y) const;
  void set_config_pars(std::string cell_map_path);  // set the cell parameters
                                                      // from the map files
  void initialize_wf_map(std::string cell_map_path);
  // fill the induced waveform map single-e waveform generation given
  // coordinates
  const std::vector<double> get_cluster_waveform(const double &c_x,
                                                 const double &c_y) const;
  // generation of a cluster vector given the hit-info
  std::vector<std::vector<double>> generate_uniform_clusters(
      const double &hit_de, const std::array<double, 2> &x1,
      const std::array<double, 2> &x2, const double &t_0) const;

 public:
  SANDTrackerDriftCellMap();                         // Default constructor
  SANDTrackerDriftCellMap(std::string signal_map_path);  // constructor with direct
                                                     // initialization of the
                                                     // parameters from the map
                                                     // path
  void init(std::string cell_map_path);
  // getter methods
  const std::array<double, 2> cell_size() const { return cell_size_; };
  const std::array<double, 3> sense_coords() const { return sense_coords_; };
  const std::vector<double> time_window() const { return time_window_; };
  const std::array<double, 3> cell_voltages() { return cell_volt_; };
  // ROOT::RDataFrame GetRDF() const { return cell_map_; };
  // build the induced current waveform for a hit-segment
  std::vector<double> build_induced_waveform(
      const std::array<double, 2> &hit_loc_start,
      const std::array<double, 2> &hit_loc_stop, const double &hit_de,
      const double &hit_t0) const;

  ClassDef(SANDTrackerDriftCellMap, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class SANDTrackerDriftCellMap + ;
#endif

#endif