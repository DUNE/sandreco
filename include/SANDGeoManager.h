#include "SANDECALCellInfo.h"
#include "SANDENDCAPModInfo.h"
#include "SANDSTTTubeInfo.h"

#include <TGeoManager.h>
#include <TPRegexp.h>
#include <TVector3.h>

#include <map>

#ifndef SANDGEOMANAGER_H
#define SANDGEOMANAGER_H

namespace sand_geometry
{
namespace stt
{
const char* const path_internal_volume =
    "volWorld_PV/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
    "MagIntVol_volume_PV_0/sand_inner_volume_PV_0";
const char* const name_internal_volume = "sand_inner_volume_PV";

const char* const stt_tube_regex_string =
    "(C|C3H6|Trk)Mod_([0-9]+)_plane(XX|YY)_straw_PV_([0-9]+)(/|)";
// "(horizontalST_(Ar|Xe)|STT_([0-9]+)_(Trk|C3H6|C)Mod(_ST|)_vv_ST)_PV_([0-9]+"
// ")(/|)";
// "_(C3H6|C|Tr)Mod_([0-9]+)_(ST_|)(hor|ver|hor2)_ST_stGas_(Xe|Ar)19_vol_PV_(["
// "0-9]+)";
const char* const stt_two_tubes_regex_string =
    "STT_([0-9]+)_(Trk|C3H6|C)Mod_(ST_|)(hh|vv)_2straw_PV_([0-9]+)<(/|)>";
const char* const stt_plane_regex_string =
    // "STT_([0-9]+)_(Trk|C3H6|C)Mod(_ST|)_(hh|vv)_PV_([0-9]+)(/|)";
    // "_(C3H6|C|Tr)Mod_([0-9]+)_(ST_|)(hor|ver|hor2)_vol_PV_0";
    "(C|C3H6|Trk)Mod_([0-9]+)_plane(XX|YY)_PV_([0-9]+)(/|)";
const char* const stt_module_regex_string =
    // "STT_([0-9]+)_(Trk|C3H6|C)Mod_PV_([0-9]+)(/|)";
    "(C|C3H6|Trk)Mod_([0-9]+)_PV_([0-9]+)(/|)";
}  // namespace stt

namespace ecal
{
const char* const path_barrel_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
    "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_lv_PV_%d";
const char* const path_endcapL_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
    "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_end_lv_PV_0";
const char* const path_endcapR_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
    "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_end_lv_PV_1";
const char* const endcap_mod_regex_string =
    "ECAL_ec_mod_([0-9]+)_lv_PV_([0-9]+)(/|)";
const char* const endcap_mod_path_regex_string =
    "ECAL_endcap_lv_PV_([0-9]+)/ECAL_ec_mod_([0-9]+)_lv_PV_([0-9]+)(/|)";

const char* const barrel_module_name = "ECAL_lv_PV";
const char* const endcap_module_name = "ECAL_end_lv_PV";
const char* const barrel_last_passive_slab_name = "volECALPassiveSlab_208_PV";

const double endcap_cell_width = 4.44;
const int number_of_layers = 5;
const int number_of_cells_per_barrel_layer = 12;
const int number_of_barrel_modules = 24;
const int number_of_cells_per_endcap_layer = 6;

// thickness of the layers in mm (barrel)
const double layer_thickness[number_of_layers] = {44., 44., 44., 44., 54.};
// thickness of the cell layers in mm (endcap)
const double ec_layer_thickness[number_of_layers] = {44.4, 44.4, 44.4, 44.4,
                                                     52.4};

// endcap module id
const int endcap_module_ids[2] = {30, 40};

namespace fluka
{
// ecal dimension for fluka
const double barrel_module_xmin = 262.55;
const double barrel_module_xmax = 292.85;
const double barrel_module_thickness = 115.0;
const double barrel_module_length = 4300;

const double endcap_rmax = 2000.0;  // ad essere precisi nella realtà è 1980
const double endcap_thickness = 115.0;
}  // namespace fluka
}  // namespace ecal
}  // namespace sand_geometry

class SANDGeoManager : public TObject
{
 private:
  TGeoManager* geo_;  // TGeoManager pointer to ND site geometry
  std::map<int, SANDECALCellInfo> cellmap_;  // map of ecal cell (key: id,
                                             // value: info on cell)
  std::map<int, SANDSTTTubeInfo> sttmap_;    // map of stt tube (key: id, value:
                                             // info on tube)

  std::map<int, SANDENDCAPModInfo> endcapmap_;  // map of the endcap modules
                                                // (key: mod id, value: mod
                                                // info)

  mutable TPRegexp stt_tube_regex_{
      sand_geometry::stt::stt_tube_regex_string};  // regular expression
                                                   // to match relevant
                                                   // info about tube
                                                   // from volume path
  mutable TPRegexp stt_two_tubes_regex_{
      sand_geometry::stt::stt_two_tubes_regex_string};  // regular expression to
                                                        // match relevant info
                                                        // about tube couple
                                                        // from volume path
  mutable TPRegexp stt_plane_regex_{
      sand_geometry::stt::stt_plane_regex_string};  // regular expression to
                                                    // match relevant info about
                                                    // plane from volume path
  mutable TPRegexp stt_module_regex_{
      sand_geometry::stt::stt_module_regex_string};  // regular expression to
                                                     // match relevant info
                                                     // about module from volume
                                                     // path
  mutable TPRegexp endcap_mod_regex_{
      sand_geometry::ecal::endcap_mod_regex_string};  // regular expression to
                                                      // match relevant info
                                                      // about endcap module
                                                      // from volume path
  mutable TPRegexp endcap_mod_path_regex_{
      sand_geometry::ecal::endcap_mod_path_regex_string};  // regular expression
                                                           // to match relevant
                                                           // info about endcap
                                                           // module from volume
                                                           // path
  std::map<int, std::map<double, int> >
      stt_tube_tranverse_position_map_;  // map (key: plane id, value: map (key:
                                         // tube id, value: 2D position [i.e. x
                                         // = z, y = transversal coord]))

  // ECAL
  std::vector<double> get_levels_z(double half_module_height,
                                   const double (&layers_thickness)[5]) const;
  int encode_ecal_barrel_cell_local_id(int layer, int cell) const;
  int encode_ecal_endcap_cell_local_id(int layer, int cell) const;
  static int encode_endcap_mod_id(int module_id, int module_replica_id,
                                  int endcap_side_id);
  static void decode_endcap_mod_id(int endcap_mod_global_id, int& module_id,
                                   int& module_replica_id, int& endcap_side_id);
  std::pair<int, int> decode_ecal_barrel_cell_local_id(int id) const;
  std::pair<int, int> decode_ecal_endcap_cell_local_id(int id) const;
  std::map<int, TVector3> get_ecal_barrel_cell_center_local_position(
      const std::vector<double>& zlevels, double m, double q) const;
  std::map<int, TVector3> get_ecal_endcap_cell_center_local_position(
      const std::vector<double>& zlevels, double rmin, double rmax) const;
  // new (alternative version)
  std::map<int, TVector3> get_ec_cell_center_local_position(
      const std::vector<double>& zlevels,
      const SANDENDCAPModInfo& module) const;

  bool is_ecal_barrel(const TString& volume_name) const;
  bool is_ecal_endcap(const TString& volume_name) const;
  bool is_endcap_mod(const TString& volume_name) const;
  bool check_and_process_ecal_path(TString& volume_path) const;
  void get_ecal_barrel_module_and_layer(const TString& volume_name,
                                        const TString& volume_path,
                                        int& detector_id, int& module_id,
                                        int& plane_id) const;
  void get_ecal_endcap_module_and_layer(const TString& volume_name,
                                        const TString& volume_path,
                                        int& detector_id, int& module_id,
                                        int& plane_id) const;
  void get_ecal_barrel_cell_local_id(double x, double y, double z,
                                     const TGeoNode* const node,
                                     int& cell_local_id) const;
  //   void get_ecal_endcap_cell_local_id(double x, double y, double z,
  //                                      const TGeoNode* const node,
  //                                      int& cell_local_id) const;
  void get_ecal_endcap_cell_local_id(double x, double y, double z,
                                     const int& endcap_mod_id,
                                     int& cell_local_id) const;
  int get_barrel_path_len(const double& hx, const double& hy, const double& hz,
                          double& d1, double& d2) const;
  int get_endcap_path_len(const double& hx, const double& hy, const double& hz,
                          const int& endcap_mod_id, double& d1,
                          double& d2) const;
  // mod id for the new endcap modules
  int get_endcap_mod_id(const TString& volume_path) const;
  void set_ecal_info();
  void set_ecal_endcap_info(const TGeoHMatrix& matrix);
  void set_ecal_endcap_info();

  // STT
  bool is_stt_tube(const TString& volume_name) const;
  bool is_stt_plane(const TString& volume_name) const;
  int get_stt_plane_id(const TString& volume_path) const;
  void set_stt_tube_info(const TGeoNode* const node, const TGeoHMatrix& matrix,
                         int stt_plane_id);
  void set_stt_info(const TGeoHMatrix& matrix);
  void set_stt_info();

 public:
  SANDGeoManager()
      : cellmap_(),
        sttmap_(),
        stt_tube_regex_(sand_geometry::stt::stt_tube_regex_string),
        stt_two_tubes_regex_(sand_geometry::stt::stt_two_tubes_regex_string),
        stt_plane_regex_(sand_geometry::stt::stt_plane_regex_string),
        stt_module_regex_(sand_geometry::stt::stt_module_regex_string),
        stt_tube_tranverse_position_map_()
  {
  }
  void init(TGeoManager* const geo);
  int save_to_file(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)
  {
    geo_ = 0;
    return Write(name, option, bufsize);
  }
  const SANDECALCellInfo& get_ecal_cell_info(int ecal_cell_id) const
  {
    return cellmap_.at(ecal_cell_id);
  }
  const SANDSTTTubeInfo& get_stt_tube_info(int stt_tube_id) const
  {
    return sttmap_.at(stt_tube_id);
  }
  const std::map<int, SANDECALCellInfo>& get_ecal_cell_info() const
  {
    return cellmap_;
  }
  const std::map<int, SANDSTTTubeInfo>& get_stt_tube_info() const
  {
    return sttmap_;
  }
  int get_ecal_cell_id(double x, double y, double z) const;
  int get_stt_tube_id(double x, double y, double z) const;

  // ECAL
  static int encode_ecal_cell_id(int detector_id, int module_id, int layer_id,
                                 int cell_local_id);
  static void decode_ecal_cell_id(int cell_global_id, int& detector_id,
                                  int& module_id, int& layer_id,
                                  int& cell_local_id);
  int get_hit_path_len(const double& hx, const double& hy, const double& hz,
                       const int& global_cell_id, double& d1, double& d2) const;

  // STT
  static int encode_stt_tube_id(int stt_plane_global_id, int stt_tube_local_id);
  static void decode_stt_tube_id(int stt_tube_global_id,
                                 int& stt_plane_global_id,
                                 int& stt_tube_local_id);
  static int encode_stt_plane_id(int stt_module_id, int stt_plane_local_id,
                                 int stt_plane_type);
  static void decode_stt_plane_id(int stt_plane_global_id, int& stt_module_id,
                                  int& stt_plane_local_id, int& stt_plane_type);

  ClassDef(SANDGeoManager, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class SANDGeoManager + ;
#endif

#endif