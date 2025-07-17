#include "SANDECALCellInfo.h"
#include "SANDENDCAPModInfo.h"
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
const std::string path_internal_volume =
    "volWorld_PV/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
    "MagIntVol_volume_PV_0/sand_inner_volume_PV_0";
const std::string name_internal_volume = "sand_inner_volume_PV";

namespace tracker
{

namespace chamber
{
const std::string wire_regex_string =
    "(C|C3H6)DriftModule_([0-2]+)(_X[0-9]_|_[A-Z]_|_)(F|S)wire_PV_([0-9]+)"
    "(/|)";
const std::string drift_plane_regex_string =
    "(C|C3H6)DriftModule_([0-2]+)(_X[0-9]_|_[A-Z]_|_)PV_([0-9])(/|)";
const std::string drift_chamber_regex_string =
    "(C|C3H6)DriftChamber(_X[0-9]_|_[A-Z]_|_)PV_0(/|)";
const std::string module_regex_string =
    "(C|C3H6)Mod(_X[0-9]_|_[A-Z]_|_)PV_([0-9]+)(/|)";
const std::string supermodule_regex_string =
    "(Trk|SuperMod)(_X[0-9]_|_[A-Z]_|_)PV_([0-1]+)(/|)";
}  // namespace chamber

namespace stt
{
const std::string path_internal_volume =
    "volWorld_PV/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
    "MagIntVol_volume_PV_0/sand_inner_volume_PV_0";
const std::string name_internal_volume = "sand_inner_volume_PV";

const std::string stt_single_tube_regex_string =
    "(C|C3H6|Trk)Mod_([0-9]+)_plane(XX|YY)_straw_PV_([0-9]+)(/|)";
// const std::string stt_two_tubes_regex_string =
//  "(Trk|C3H6|C)Mod_([0-9]+)_plane(XX|YY)_2straw_stGas_(Xe|Ar)19_PV_([0-9]+)(/|)";

const std::string stt_plane_regex_string =
    "(C|C3H6|Trk)Mod_([0-9]+)_plane(XX|YY)_PV_([0-9]+)(/|)";
const std::string stt_module_regex_string =
    "(C|C3H6|Trk)Mod_([0-9]+)_PV_([0-9]+)(/|)";

const std::string stt_supermodule_regex_string =
    "(Trk|SuperMod)(_X0_|_X1_|_A_|_B_|_C_|_)PV_([0-1]+)(/|)";
}  // namespace stt
} // namespace tracker
namespace ecal
{
const std::string path_barrel_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
    "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_lv_PV_%d";
const std::string path_endcapL_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
    "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_end_lv_PV_0";
const std::string path_endcapR_template =
    "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
    "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_end_lv_PV_1";
const std::string endcap_mod_regex_string =
    "ECAL_ec_mod_([0-9]+)_lv_PV_([0-9]+)(/|)";
const std::string endcap_mod_path_regex_string =
    "ECAL_endcap_lv_PV_([0-9]+)/ECAL_ec_mod_([0-9]+)_lv_PV_([0-9]+)(/|)";

const std::string barrel_module_name = "ECAL_lv_PV";
const std::string endcap_module_name = "ECAL_end_lv_PV";
const std::string barrel_last_passive_slab_name = "volECALPassiveSlab_208_PV";

const double endcap_cell_width = 44.4;
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
  std::map<int,  sand_geometry::ecal::ECALCellInfo> cellmap_;  // map of ecal cell (key: id,
                                             // value: info on cell)

  std::vector<sand_geometry::tracker::Plane> planes_;
  std::map<sand_geometry::tracker::PlaneID, sand_geometry::tracker::plane_iterator> id_to_plane_;


  std::map<int,  sand_geometry::ecal::ENDCAPModInfo> endcapmap_;  // map of the endcap modules
                                                // (key: mod id, value: mod
                                                // info)

  mutable TPRegexp stt_tube_regex_ = TPRegexp(
      sand_geometry::tracker::stt::stt_single_tube_regex_string);  // regular expression
                                                          // to match relevant
                                                          // info about tube
                                                          // from volume path
  mutable TPRegexp stt_plane_regex_ = TPRegexp(
      sand_geometry::tracker::stt::stt_plane_regex_string);  // regular expression to
                                                    // match relevant info about
                                                    // plane from volume path
  mutable TPRegexp stt_module_regex_ = TPRegexp(
      sand_geometry::tracker::stt::stt_module_regex_string);  // regular expression to
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
  mutable TPRegexp stt_supermodule_regex_ = TPRegexp(
      sand_geometry::tracker::stt::stt_supermodule_regex_string);  // regular expression
                                                          // to
                                                          // match relevant info
                                                          // about supermodule from volume
                                                          // path


  bool getLineSegmentIntersection(TVector2 p, TVector2 dir, TVector2 A, TVector2 B, TVector2& intersection);
  void setDriftPlaneInfo(sand_geometry::tracker::Plane& plane, double angle);
  void printModulesInfo(int verbose = 1);
  void drawModulesInfo();

  // DRIFT CHAMBER
  mutable TPRegexp wire_regex_ = TPRegexp(sand_geometry::tracker::chamber::wire_regex_string);
  mutable TPRegexp drift_plane_regex_ = TPRegexp(
      sand_geometry::tracker::chamber::drift_plane_regex_string);
  mutable TPRegexp drift_chamber_regex_ = TPRegexp(
      sand_geometry::tracker::chamber::drift_chamber_regex_string);
  mutable TPRegexp module_regex_ = TPRegexp(sand_geometry::tracker::chamber::module_regex_string);
  mutable TPRegexp supermodule_regex_ = TPRegexp(
      sand_geometry::tracker::chamber::supermodule_regex_string);

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
  std::map<int, TVector3> get_ec_cell_center_local_position(
      const std::vector<double>& zlevels,
      const  sand_geometry::ecal::ENDCAPModInfo& module) const;

  bool is_ecal_barrel(const TString& volume_name, bool include_passive) const;
  bool is_ecal_endcap(const TString& volume_name, bool include_passive) const;
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
  int get_barrel_hit_pos(double d1, int global_cellID,
                         double& reco_x, double& reco_y, double& reco_z) const;
  int get_endcap_path_len(const double& hx, const double& hy, const double& hz,
                          const int& endcap_mod_id, double& d1,
                          double& d2) const;
  int get_endcap_hit_pos(const double& d1, const int& global_cellID,
                         const int& modID, double& reco_x, double& reco_y,
                         double& reco_z) const;
//   double compute_cell_d1(const double& cell_l, const double& tdc_1,
//                          const double& tdc_2) const;
//   double compute_cell_d2(const double& cell_l, const double& tdc_1,
//                          const double& tdc_2) const;
  // mod id for the new endcap modules
  int get_endcap_mod_id(const TString& volume_path) const;
  void set_ecal_info();
  void set_ecal_endcap_info(const TGeoHMatrix& matrix);
  void set_ecal_endcap_info();
  
  void setTrackerInfo();

  void rearrangePlanes();

  std::vector<TVector2> getLocalLinePlaneIntersections(const TVector2& local_2d_position,
                                                       const sand_geometry::tracker::Plane& plane);
  std::vector<TVector2> getGlobalLinePlaneIntersections(const TVector2& local_2d_position, 
                                                        const sand_geometry::tracker::Plane& plane);
  double getMinDistanceBetweenSegments(TVector3 a, TVector3 b,
                                       TVector3 c, TVector3 d);
  // STT
  sand_geometry::tracker::ModuleID getSttModuleId(const TString& volume_path) const;
  bool isSttTube(const TString& volume_name) const;
  bool isSttPlane(const TString& volume_name) const;
  sand_geometry::tracker::PlaneID getSttPlaneId(const TString& volume_path, bool justLocal) const;
  void setSttWireInfo(sand_geometry::tracker::Plane& plane, const TGeoNode* const node, const TGeoHMatrix& matrix);
  void setSttPlaneInfo(const TGeoNode* const node, const TGeoHMatrix& matrix);

  // DRIFT CHAMEBER
  void setPlaneInfo(const TGeoHMatrix& matrix);
  void setDriftPlaneInfo(const TGeoNode* const node, const TGeoHMatrix& matrix);
  void setDriftWireInfo(sand_geometry::tracker::Plane& plane);
  sand_geometry::tracker::PlaneID getDriftPlaneId(const TString& volume_path, bool JustLocalId) const;
  sand_geometry::tracker::ModuleID getDriftModuleId(const TString& volume_path) const;
  sand_geometry::tracker::ModuleID getDriftSupermoduleId(const TString& volume_path) const;
  sand_geometry::tracker::ModuleID getDriftModuleReplicaId(const TString& volume_path) const;
  sand_geometry::tracker::WireID getWireId(const TString& volume_path) const;
  bool isDriftPlane(const TString& volume_name) const;
  bool isSwire(const TString& volume_path) const;
  void writeMapOnFile(std::string fName,
                      const std::map<sand_geometry::tracker::WireID, sand_geometry::tracker::WireInfo>& map);

 public:
  SANDGeoManager()
      : cellmap_(),
        stt_tube_regex_(sand_geometry::tracker::stt::stt_single_tube_regex_string),
        // stt_two_tubes_regex_(sand_geometry::tracker::stt::stt_two_tubes_regex_string),
        stt_plane_regex_(sand_geometry::tracker::stt::stt_plane_regex_string),
        stt_module_regex_(sand_geometry::tracker::stt::stt_module_regex_string),
        stt_supermodule_regex_(
            sand_geometry::tracker::stt::stt_supermodule_regex_string)
  {
  }
  void init(TGeoManager* const geo);
  void setGeoCurrentPoint(double x, double y, double z) const;
  void setGeoCurrentDirection(double x, double y, double z) const;
  void initVolume(volume& v) const;
  const  sand_geometry::ecal::ECALCellInfo& get_ecal_cell_info(int ecal_cell_id) const
  {
    return cellmap_.at(ecal_cell_id);
  }
  void fillAdjacentCells(std::string geometry);
  std::map<sand_geometry::tracker::CellID, sand_geometry::tracker::Cell>::const_iterator getCellInfo(sand_geometry::tracker::CellID cell_id) const;
  sand_geometry::tracker::plane_iterator getPlaneInfo(sand_geometry::tracker::CellID cell_id) const;
  sand_geometry::tracker::plane_iterator getPlaneInfo(sand_geometry::tracker::PlaneID unique_plane_id) const;
  const std::map<int,  sand_geometry::ecal::ECALCellInfo>& get_ecal_cell_info() const
  {
    return cellmap_;
  }
  const std::vector<sand_geometry::tracker::Plane>&
      getPlanes() const
  {
    return planes_;
  }
  std::vector<sand_geometry::tracker::Plane>&
      getPlanes()
  {
    return planes_;
  }

  const sand_geometry::tracker::PlaneIndex getPlaneIndex(const sand_geometry::tracker::PlaneID& plane_uid) const
  {
    return sand_geometry::tracker::PlaneIndex(std::distance(planes_.cbegin(), id_to_plane_.at(plane_uid)));
  }
  int get_ecal_cell_id(double x, double y, double z, bool include_passive) const;
  sand_geometry::tracker::CellID getSttTubeId(double x, double y, double z) const;
  long printSttTubeId(double x, double y, double z) const;
  
  // Notice: is the non-const version needed?
  const TGeoManager* getTGeoManager() const {return geo_;};
  TGeoManager* getTGeoManager() {return geo_;};
  long get_wire_id(long drift_plane_id, double z,
                   double transverse_coord) const;
  std::vector<sand_geometry::tracker::CellID> getSegmentIds(const TG4HitSegment& hseg) const;
  TVector3 findClosestDrift(TVector3 point, double epsilon) const;
  TVector3 smearPoint(TVector3 point, double epsilon) const;
  
  const TVector2 pointInRotatedSystem(TVector2 v, double angle) const;
  const TVector2 globalToLocal(TVector2 global, const sand_geometry::tracker::Plane& plane) const;
  const TVector2 localToRotated(TVector2 local, const sand_geometry::tracker::Plane& plane) const;
  const TVector2 globalToRotated(TVector2 global, const sand_geometry::tracker::Plane& plane) const;
  const TVector2 rotatedToLocal(TVector2 rotated, const sand_geometry::tracker::Plane& plane) const;
  const TVector2 localToGlobal(TVector2 local, const sand_geometry::tracker::Plane& plane) const;
  const TVector2 rotatedToGlobal(TVector2 rotated, const sand_geometry::tracker::Plane& plane) const;
  sand_geometry::tracker::CellID getClosestCellToHit(TVector3 hit_center, const sand_geometry::tracker::Plane& plane, bool checkCloseCells) const;
  double getHitCellDistance(TVector2 rotated_local_yz_hit_position, 
                                        std::map<sand_geometry::tracker::CellID, sand_geometry::tracker::Cell>::const_iterator cell_it, 
                                        const sand_geometry::tracker::Plane& plane) const;

  double compute_cell_d1(const double& cell_l, const double& tdc_1,
                         const double& tdc_2) const;
  double compute_cell_d2(const double& cell_l, const double& tdc_1,
                         const double& tdc_2) const;

  // ECAL
  static int encode_ecal_cell_id(int detector_id, int module_id, int layer_id,
                                 int cell_local_id);
  static void decode_ecal_cell_id(int cell_global_id, int& detector_id,
                                  int& module_id, int& layer_id,
                                  int& cell_local_id);
  TString FindNextActiveLayer(const double* starting_point,
                              const double* direction) const;
  int get_hit_path_len(const double& hx, const double& hy, const double& hz,
                       const int& global_cell_id, double& d1, double& d2) const;
                       
  int get_reco_hit_pos(const int& cellID, const double& cell_l,
                       const double& tdc_1, const double& tdc_2, double& reco_x,
                       double& reco_y, double& reco_z) const;

  // STT
  static sand_geometry::tracker::CellID encodeCellId(sand_geometry::tracker::PlaneID plane_global_id, sand_geometry::tracker::CellID wire_local_id);
  static void decodeCellId(sand_geometry::tracker::CellID   cell_global_id, 
                             sand_geometry::tracker::PlaneID& plane_global_id,
                             sand_geometry::tracker::CellID&  cell_local_id);
  static sand_geometry::tracker::PlaneID encodePlaneId(sand_geometry::tracker::ModuleID unique_module_id,
                              sand_geometry::tracker::PlaneID plane_local_id, sand_geometry::tracker::PlaneID plane_type);
  static void decodePlaneId(sand_geometry::tracker::PlaneID plane_global_id, sand_geometry::tracker::ModuleID& unique_module_id, 
                              sand_geometry::tracker::PlaneID& plane_local_id, sand_geometry::tracker::PlaneID& plane_type);
  static sand_geometry::tracker::ModuleID encodeModuleId(sand_geometry::tracker::ModuleID supermodule_id, 
                               sand_geometry::tracker::ModuleID module_id, sand_geometry::tracker::ModuleID module_replica_id);
  static void decodeModuleId(sand_geometry::tracker::ModuleID unique_module_id, sand_geometry::tracker::ModuleID& supermodule_id, 
                               sand_geometry::tracker::ModuleID& module_id, sand_geometry::tracker::ModuleID& module_replica_id);

  ClassDef(SANDGeoManager, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class SANDGeoManager + ;
#endif

#endif
