#include "SANDECALCellInfo.h"
#include "SANDWireInfo.h"
#include "SANDTrackerModule.h"
#include "struct.h"

#include <TGeoManager.h>
#include <TPRegexp.h>
#include <TVector3.h>
#include <TG4HitSegment.h>

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
    "(C|C3H6)DriftModule_([0-2]+)(_X0_|_X1_|_A_|_B_|_C_|_)(F|S)wire_PV_([0-9]+)"
    "(/|)";
const char* const drift_plane_regex_string =
    "(C|C3H6)DriftModule_([0-2]+)(_X0_|_X1_|_A_|_B_|_C_|_)PV_([0-9])(/|)";
const char* const drift_chamber_regex_string =
    "(C|C3H6)DriftChamber(_X0_|_X1_|_A_|_B_|_C_|_)PV_0(/|)";
const char* const module_regex_string =
    "(C|C3H6)Mod(_X0_|_X1_|_A_|_B_|_C_|_)PV_([0-9]+)(/|)";
const char* const supermodule_regex_string =
    "(Trk|SuperMod)(_X0_|_X1_|_A_|_B_|_C_|_)PV_([0-1]+)(/|)";
}  // namespace chamber

namespace stt
{
const char* const path_internal_volume =
    "volWorld_PV/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
    "MagIntVol_volume_PV_0/sand_inner_volume_PV_0";
const char* const name_internal_volume = "sand_inner_volume_PV";

const char* const stt_single_tube_regex_string =
    "(C|C3H6|Trk)Mod_([0-9]+)_plane(XX|YY)_straw_PV_([0-9]+)(/|)";
// const char* const stt_two_tubes_regex_string =
//  "(Trk|C3H6|C)Mod_([0-9]+)_plane(XX|YY)_2straw_stGas_(Xe|Ar)19_PV_([0-9]+)(/|)";

const char* const stt_plane_regex_string =
    "(C|C3H6|Trk)Mod_([0-9]+)_plane(XX|YY)_PV_([0-9]+)(/|)";
const char* const stt_module_regex_string =
    "(C|C3H6|Trk)Mod_([0-9]+)_PV_([0-9]+)(/|)";

// const char* const stt_single_tube_regex_string =
//     "(horizontalST_(Ar|Xe)|STT_([0-9]+)_(Trk|C3H6|C)Mod(_ST|)_vv_ST)_PV_([0-9]+"
//     ")(/|)";
// //
// "_(C3H6|C|Tr)Mod_([0-9]+)_(ST_|)(hor|ver|hor2)_ST_stGas_(Xe|Ar)19_vol_PV_(["
// // "0-9]+)";
// const char* const stt_two_tubes_regex_string =
//     "STT_([0-9]+)_(Trk|C3H6|C)Mod_(ST_|)(hh|vv)_2straw_PV_([0-9]+)<(/|)>";
// const char* const stt_plane_regex_string =
//     // "STT_([0-9]+)_(Trk|C3H6|C)Mod(_ST|)_(hh|vv)_PV_([0-9]+)(/|)";
//     // "_(C3H6|C|Tr)Mod_([0-9]+)_(ST_|)(hor|ver|hor2)_vol_PV_0";
//     "(C|C3H6|Trk)Mod_([0-9]+)_plane(XX|YY)_PV_([0-9]+)(/|)";
// const char* const stt_module_regex_string =
//     "STT_([0-9]+)_(Trk|C3H6|C)Mod_PV_([0-9]+)(/|)";
const char* const stt_supermodule_regex_string =
    "(Trk|SuperMod)(_X0_|_X1_|_A_|_B_|_C_|_)PV_([0-1]+)(/|)";
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

class Counter
{
  // private:
 public:
  std::map<std::string, int> hit_counter_;
  void IncrementCounter(std::string k);
  void PrintCounter();
};

class SANDGeoManager : public TObject
{
 private:
  TGeoManager* geo_;  // TGeoManager pointer to ND site geometry
  std::map<int, SANDECALCellInfo> cellmap_;  // map of ecal cell (key: id,
                                             // value: info on cell)

  std::map<long, SANDWireInfo> wiremap_;  // map of wire (key : id, value:
                                          // info on wire)

  std::map<long, std::map<double, long>>
      wire_tranverse_position_map_;  // map (key: plane id, value: map (key:
                                     // tube id, value: 2D position [i.e. x
                                     // = z, y = transversal coord]))
                                     // map (key: plane id, value: map (key:
                                     // wire id, value: 2D position [i.e. x
                                     // = z, y = transversal coord]))
  std::map<long, SANDTrackerModule> _tracker_modules_map;

  mutable TPRegexp stt_tube_regex_{
      sand_geometry::stt::stt_single_tube_regex_string};  // regular expression
                                                          // to match relevant
                                                          // info about tube
                                                          // from volume path
  // mutable TPRegexp stt_two_tubes_regex_{
  // sand_geometry::stt::stt_two_tubes_regex_string};  // regular expression to
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
  mutable TPRegexp stt_supermodule_regex_{
      sand_geometry::stt::stt_supermodule_regex_string};  // regular expression
                                                          // to
                                                          // match relevant info
                                                          // about supermodule from volume
                                                          // path


  const TVector2 pointInRotatedSystem(TVector2 v, double angle) const;
  long GetClosestCellToHit(TVector3 hit_center, const SANDTrackerPlane& plane, bool checkCloseCells) const;
  double GetHitCellDistance(TVector2 rotated_local_yz_hit_position, 
                                        std::map<long, SANDTrackerCell>::const_iterator cell_it, 
                                        const SANDTrackerPlane& plane) const;
  bool getLineSegmentIntersection(TVector2 p, TVector2 dir, TVector2 A, TVector2 B, TVector3& intersection);
  void set_drift_plane_info(SANDTrackerPlane& plane, double angle);
  void PrintModulesInfo(int verbose = 1);
  void DrawModulesInfo();

  // DRIFT CHAMBER
  mutable TPRegexp wire_regex_{sand_geometry::chamber::wire_regex_string};
  mutable TPRegexp drift_plane_regex_{
      sand_geometry::chamber::drift_plane_regex_string};
  mutable TPRegexp drift_chamber_regex_{
      sand_geometry::chamber::drift_chamber_regex_string};
  mutable TPRegexp module_regex_{sand_geometry::chamber::module_regex_string};
  mutable TPRegexp supermodule_regex_{
      sand_geometry::chamber::supermodule_regex_string};

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
  long get_stt_module_id(const TString& volume_path) const;
  bool is_stt_tube(const TString& volume_name) const;
  bool is_stt_plane(const TString& volume_name) const;
  int get_stt_plane_id(const TString& volume_path, bool justLocal) const;
  void set_stt_wire_info(SANDTrackerPlane& plane, const TGeoNode* const node, const TGeoHMatrix& matrix);
  void set_stt_plane_info(const TGeoNode* const node, const TGeoHMatrix& matrix);

  // DRIFT CHAMEBER
  void set_wire_info(const TGeoHMatrix& matrix);
  void set_drift_plane_info(const TGeoNode* const node, const TGeoHMatrix& matrix);
  void set_drift_wire_info(SANDTrackerPlane& plane);
  long get_drift_plane_id(const TString& volume_path, bool JustLocalId) const;
  long get_drift_module_id(const TString& volume_path) const;
  int get_drift_supermodule_id(const TString& volume_path) const;
  long get_drift_module_replica_id(const TString& volume_path) const;
  long get_wire_id(const TString& volume_path) const;
  bool is_drift_plane(const TString& volume_name) const;
  bool isSwire(const TString& volume_path) const;
  void WriteMapOnFile(std::string fName,
                      const std::map<long, SANDWireInfo>& map);

 public:
  SANDGeoManager()
      : cellmap_(),
        wiremap_(),
        stt_tube_regex_(sand_geometry::stt::stt_single_tube_regex_string),
        // stt_two_tubes_regex_(sand_geometry::stt::stt_two_tubes_regex_string),
        stt_plane_regex_(sand_geometry::stt::stt_plane_regex_string),
        stt_module_regex_(sand_geometry::stt::stt_module_regex_string),
        stt_supermodule_regex_(
            sand_geometry::stt::stt_supermodule_regex_string),
        wire_tranverse_position_map_()
  {
  }
  void init(TGeoManager* const geo);
  void SetGeoCurrentPoint(double x, double y, double z) const;
  void SetGeoCurrentDirection(double x, double y, double z) const;
  void InitVolume(volume& v) const;
  void LOGVolumeInfo(volume& v) const;
  void PrintCounter();
  int save_to_file(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)
  {
    geo_ = 0;
    return Write(name, option, bufsize);
  }
  const SANDECALCellInfo& get_ecal_cell_info(int ecal_cell_id) const
  {
    return cellmap_.at(ecal_cell_id);
  }
  const SANDTrackerCell& get_cell_info(long wire_id) const;
  const std::map<int, SANDECALCellInfo>& get_ecal_cell_info() const
  {
    return cellmap_;
  }
  const std::map<long, SANDWireInfo>& get_wire_info() const
  {
    return wiremap_;
  }
  const std::map<long, std::map<double, long>>&
      get_wires_transverse_position_map() const
  {
    return wire_tranverse_position_map_;
  }
  int get_ecal_cell_id(double x, double y, double z) const;
  long get_stt_tube_id(double x, double y, double z) const;
  long print_stt_tube_id(double x, double y, double z) const;

  long get_wire_id(long drift_plane_id, double z,
                   double transverse_coord) const;
  std::vector<long> get_segment_ids(const TG4HitSegment& hseg) const;
  TVector3 FindClosestDrift(TVector3 point, double epsilon) const;
  TVector3 SmearPoint(TVector3 point, double epsilon) const;
  bool IsOnEdge(TVector3 point) const;

  // ECAL
  static int encode_ecal_cell_id(int detector_id, int module_id, int layer_id,
                                 int cell_local_id);
  static void decode_ecal_cell_id(int cell_global_id, int& detector_id,
                                  int& module_id, int& layer_id,
                                  int& cell_local_id);
  TString FindNextActiveLayer(const double* starting_point,
                              const double* direction) const;

  // STT
  static long encode_wire_id(long plane_global_id, long wire_local_id);
  static void decode_wire_id(long wire_global_id, long& plane_global_id,
                             long& wire_local_id);
  static long encode_plane_id(long unique_module_id,
                              long plane_local_id, long plane_type);
  static void decode_plane_id(long plane_global_id, long& unique_module_id, 
                              long& plane_local_id, long& plane_type);
  static long encode_module_id(long supermodule_id, 
                               long module_id, long module_replica_id);
  static void decode_module_id(long unique_module_id, long& supermodule_id, 
                               long& module_id, long& module_replica_id);
  // DRIFT CHAMBER
  // static void decode_chamber_plane_id(int wire_global_id,
  //                            int& drift_plane_global_id,
  //                            int& wire_local_id);

  ClassDef(SANDGeoManager, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class SANDGeoManager + ;
#endif

#endif
