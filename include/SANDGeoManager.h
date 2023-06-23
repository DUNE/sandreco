#include "SANDECALCellInfo.h"
#include "SANDSTTTubeInfo.h"

#include <TGeoManager.h>
#include <TPRegexp.h>
#include <TVector3.h>

#include <map>

#ifndef SANDGEOMANAGER_H
#define SANDGEOMANAGER_H

namespace sand_geometry
{
const char* const path_internal_volume =
    "volWorld_PV/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
    "MagIntVol_volume_PV_0/sand_inner_volume_PV_0";
const char* const name_internal_volume = "sand_inner_volume_PV";

namespace chamber
{
const char* const wire_regex_string =
   "(C|C3H6)DriftModule_([0-2]+)(_X0_|_X1_|_A_|_B_|_C_|_)(F|S)wire_PV_([0-9]+)(/|)";
const char* const drift_plane_regex_string = 
   "(C|C3H6)DriftModule_([0-2]+)(_X0_|_X1_|_A_|_B_|_C_|_)PV_0(/|)";
const char* const drift_chamber_regex_string = 
   "(C|C3H6)DriftChamber(_X0_|_X1_|_A_|_B_|_C_|_)PV_0(/|)";
const char* const module_regex_string =
   "(C|C3H6)Mod(_X0_|_X1_|_A_|_B_|_C_|_)PV_([0-9]+)(/|)";
const char* const supermodule_regex_string =
   "SuperMod(_X0_|_X1_|_A_|_B_|_C_|_)PV_([0-1]+)(/|)";
} // namespace chamber
    
namespace stt
{
const char* const stt_single_tube_regex_string =
    "(horizontalST_(Ar|Xe)|STT_([0-9]+)_(Trk|C3H6|C)Mod(_ST|)_vv_ST)_PV_([0-9]+"
    ")(/|)";
// "_(C3H6|C|Tr)Mod_([0-9]+)_(ST_|)(hor|ver|hor2)_ST_stGas_(Xe|Ar)19_vol_PV_(["
// "0-9]+)";
const char* const stt_two_tubes_regex_string =
    "STT_([0-9]+)_(Trk|C3H6|C)Mod_(ST_|)(hh|vv)_2straw_PV_([0-9]+)(/|)";
const char* const stt_plane_regex_string =
    "STT_([0-9]+)_(Trk|C3H6|C)Mod(_ST|)_(hh|vv)_PV_([0-9]+)(/|)";
// "_(C3H6|C|Tr)Mod_([0-9]+)_(ST_|)(hor|ver|hor2)_vol_PV_0";
const char* const stt_module_regex_string =
    "STT_([0-9]+)_(Trk|C3H6|C)Mod_PV_([0-9]+)(/|)";
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
const char* const barrel_module_name = "ECAL_lv_PV";
const char* const endcap_module_name = "ECAL_end_lv_PV";

const int number_of_layers = 5;
const int number_of_cells_per_barrel_layer = 12;
const int number_of_barrel_modules = 24;
const int number_of_cells_per_endcap_layer = 90;

// thickness of the layers in mm
const double layer_thickness[number_of_layers] = {44., 44., 44., 44., 54.};

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
  std::map<int, SANDWireInfo> sttmap_;    // map of stt tube (key: id, value:
                                             // info on tube)

  std::map<int, SANDWireInfo> wiremap_; // map of wire (key : id, value:
                                            // info on wire)
  mutable TPRegexp stt_single_tube_regex_{
      sand_geometry::stt::stt_single_tube_regex_string};  // regular expression
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
  std::map<int, std::map<double, int> >
      stt_tube_tranverse_position_map_;  // map (key: plane id, value: map (key:
                                         // tube id, value: 2D position [i.e. x
                                         // = z, y = transversal coord]))

  // DRIFT CHAMBER
    mutable TPRegexp wire_regex_{
      sand_geometry::chamber::wire_regex_string};
    mutable TPRegexp drift_plane_regex_{
      sand_geometry::chamber::drift_plane_regex_string};
    mutable TPRegexp drift_chamber_regex_{
      sand_geometry::chamber::drift_chamber_regex_string};
    mutable TPRegexp module_regex_{
      sand_geometry::chamber::module_regex_string};
    mutable TPRegexp supermodule_regex_{
      sand_geometry::chamber::supermodule_regex_string};
  }                                         
                                           

  // ECAL
  std::vector<double> get_levels_z(double half_module_height) const;
  int encode_ecal_barrel_cell_local_id(int layer, int cell) const;
  int encode_ecal_endcap_cell_local_id(int layer, int cell) const;
  std::pair<int, int> decode_ecal_barrel_cell_local_id(int id) const;
  std::pair<int, int> decode_ecal_endcap_cell_local_id(int id) const;
  std::map<int, TVector3> get_ecal_barrel_cell_center_local_position(
      const std::vector<double>& zlevels, double m, double q) const;
  std::map<int, TVector3> get_ecal_endcap_cell_center_local_position(
      const std::vector<double>& zlevels, double rmin, double rmax) const;
  bool is_ecal_barrel(const TString& volume_name) const;
  bool is_ecal_endcap(const TString& volume_name) const;
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
  void get_ecal_endcap_cell_local_id(double x, double y, double z,
                                     const TGeoNode* const node,
                                     int& cell_local_id) const;
  void set_ecal_info();

  void set_wire_info();
  // STT
  bool is_stt_tube(const TString& volume_name) const;
  bool is_stt_plane(const TString& volume_name) const;
  int get_stt_plane_id(const TString& volume_path) const;
  void set_stt_tube_info(const TGeoNode* const node, const TGeoHMatrix& matrix,
                         int stt_plane_id);
  void set_stt_info(const TGeoHMatrix& matrix);
  // void set_stt_info();

  // DRIFT CHAMEBER
  void set_wire_info(const TGeoHMatrix& matrix);

 public:
  SANDGeoManager()
      : cellmap_(),
        sttmap_(),
        // stt_single_tube_regex_(
        //     sand_geometry::stt::stt_single_tube_regex_string),
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
  const SANDWireInfo& get_stt_tube_info(int stt_tube_id) const
  {
    return sttmap_.at(stt_tube_id);
  }
  const SANDWireInfo& get_wire_info(int wire_id) const
  {
    return wiremap_.at(wire_id);
  }
  const std::map<int, SANDECALCellInfo>& get_ecal_cell_info() const
  {
    return cellmap_;
  }
  const std::map<int, SANDWireInfo>& get_stt_tube_info() const
  {
    return sttmap_;
  }
  const std::map<int, SANDWireInfo>& get_wire_info() const
  {
    return wiremap_;
  }
  int get_ecal_cell_id(double x, double y, double z) const;
  int get_stt_tube_id(double x, double y, double z) const;
  int get_wire_id(double x, double y, double z) const;

  // ECAL
  static int encode_ecal_cell_id(int detector_id, int module_id, int layer_id,
                                 int cell_local_id);
  static void decode_ecal_cell_id(int cell_global_id, int& detector_id,
                                  int& module_id, int& layer_id,
                                  int& cell_local_id);

  // STT
  static int encode_stt_tube_id(int stt_plane_global_id, int stt_tube_local_id);
  static void decode_stt_tube_id(int stt_tube_global_id,
                                 int& stt_plane_global_id,
                                 int& stt_tube_local_id);
  static int encode_stt_plane_id(int stt_module_id, int stt_plane_local_id,
                                 int stt_plane_type);
  static void decode_stt_plane_id(int stt_plane_global_id, int& stt_module_id,
                                  int& stt_plane_local_id, int& stt_plane_type);

    // DRIFT CHAMBER
  static int encode_wire_id(int drift_plane_global_id, int wire_local_id);
  static void decode_wire_id(int wire_global_id, int& drift_plane_global_id, int& wire_local_id);

  ClassDef(SANDGeoManager, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class SANDGeoManager + ;
#endif

#endif