#include "SANDGeoManager.h"
#include "SANDTrackerModuleConfig.h"
#include "utils.h"

#include <iostream>
#include <fstream>

#include <iomanip>


#include <TGeoTrd2.h>
#include <TGeoTube.h>
#include <TGeoBBox.h>
#include <TObjString.h>
#include <TRandom3.h>
#include <TH2D.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TMatrixD.h>
#include <TStyle.h>

//  #######################################################
//  ##                     ECAL_BARREL                   ##
//  #######################################################

int SANDGeoManager::encode_ecal_barrel_cell_local_id(int layer, int cell) const
{
  return cell * 100 + layer;
}

int SANDGeoManager::encode_ecal_endcap_cell_local_id(int layer, int cell) const
{
  return cell * 100 + layer;
}

std::pair<int, int> SANDGeoManager::decode_ecal_barrel_cell_local_id(int id)
    const
{
  int cell = id / 100;
  int layer = id % 100;
  return std::make_pair(layer, cell);
}

std::pair<int, int> SANDGeoManager::decode_ecal_endcap_cell_local_id(int id)
    const
{
  int cell = id / 100;
  int layer = id % 100;
  return std::make_pair(layer, cell);
}

std::vector<double> SANDGeoManager::get_levels_z(
    double half_module_height, const double (&layers_thickness)[5]) const
{
  // z edge of the cells
  std::vector<double> zlevel;
  zlevel.push_back(-half_module_height);

  for (int i = 0; i < sand_geometry::ecal::number_of_layers; i++) {
    zlevel.push_back(zlevel.back() + layers_thickness[i]);
  }
  return zlevel;
}

std::map<int, TVector3>
    SANDGeoManager::get_ecal_barrel_cell_center_local_position(
        const std::vector<double>& zlevels, double m, double q) const
{
  // z position of the center of the cells
  std::map<int, TVector3> ecal_barrel_cell_center_local_positions;
  for (auto i = 0u; i < zlevels.size() - 1u; i++) {
    auto z_this_layer = 0.5 * (zlevels.at(i) + zlevels.at(i + 1));

    // total module width at the z position of the center of the cell
    double x_module_width_at_z = 2 * (m * z_this_layer + q);

    // cell width at the z position of the center of the cell
    double x_cell_width = x_module_width_at_z /
                          sand_geometry::ecal::number_of_cells_per_barrel_layer;

    // position of the center of the cells
    for (int j = 0; j < sand_geometry::ecal::number_of_cells_per_barrel_layer;
         j++) {
      auto x = x_cell_width * (j + 0.5) - x_module_width_at_z * 0.5;
      auto y = 0.;
      auto z = z_this_layer;
      auto id = encode_ecal_barrel_cell_local_id(i, j);
      ecal_barrel_cell_center_local_positions[id] = TVector3(x, y, z);
    }
  }
  return ecal_barrel_cell_center_local_positions;
}

std::map<int, TVector3> SANDGeoManager::get_ec_cell_center_local_position(
    const std::vector<double>& zlevels, const sand_geometry::ecal::ENDCAPModInfo& module) const
{
  std::map<int, TVector3> ecal_endcap_cell_center_local_positions;
  for (auto i = 0u; i < zlevels.size() - 1u; i++) {
    auto z_this_layer = 0.5 * (zlevels.at(i) + zlevels.at(i + 1));
    // the cell width along x in the endcaps is fixed (44.4 mm)
    int n_cells = (module.getWidth() / sand_geometry::ecal::endcap_cell_width);
    for (int j = 0; j < n_cells; j++) {
      auto x = sand_geometry::ecal::endcap_cell_width * (j + 0.5) -
               0.5 * module.getWidth();
      auto z = z_this_layer;
      auto y = module.getCurvatureArcLength(z + module.getModDz()) + module.getLHor() +
               0.5 * module.getLVert() -
               0.5 * module.getCellTotalLength(z + module.getModDz());
      auto id = encode_ecal_endcap_cell_local_id(i, j);
      ecal_endcap_cell_center_local_positions[id] = TVector3(x, y, z);
    }
  }
  return ecal_endcap_cell_center_local_positions;
}

int SANDGeoManager::encode_ecal_cell_id(int detector_id, int module_id,
                                        int layer_id, int cell_local_id)
{
  return cell_local_id + 100 * layer_id + 1000 * module_id + detector_id * 1e7;
}

void SANDGeoManager::decode_ecal_cell_id(int cell_global_id, int& detector_id,
                                         int& module_id, int& layer_id,
                                         int& cell_local_id)
{
  detector_id = cell_global_id / 1e7;
  cell_global_id -= detector_id * 1e7;

  module_id = cell_global_id / 1000;
  cell_global_id -= module_id * 1000;

  layer_id = cell_global_id / 100;
  cell_global_id -= layer_id * 100;

  cell_local_id = cell_global_id;
}

TString SANDGeoManager::FindNextActiveLayer(const double* starting_point,
                                            const double* direction) const
{
  // step from current point haed till you find an Active layer
  int max_nof_steps = 3;

  int nof_steps = 0;

  geo_->SetCurrentPoint(starting_point[0], starting_point[1],
                        starting_point[2]);

  geo_->SetCurrentDirection(direction[0], direction[1], direction[2]);

  TString current_node = geo_->GetCurrentNode()->GetName();

  while (!current_node.Contains("Active")) {
    geo_->FindNextBoundaryAndStep();
    current_node = geo_->GetCurrentNode()->GetName();
    // std::cout << current_node << "\n";
    if (current_node.Contains("Passive"))
      current_node.ReplaceAll("Passive", "Active");
    if ((nof_steps > max_nof_steps) &&
        (!current_node.Contains(
              "Active"))) {  // invert direction to find Active volume
      double _direction[3] = {-direction[0], -direction[1], -direction[2]};
      FindNextActiveLayer(starting_point, _direction);
    }
    nof_steps++;
  }
  return current_node;
}

bool SANDGeoManager::is_ecal_barrel(const TString& volume_name, bool include_passive) const
{
  // something like: volECALActiveSlab_21_PV_0
  // To Do: add a more refined check for passive and active volumes
  if (include_passive) {
    return volume_name.Contains("volECAL") == true &&
           (volume_name.Contains("Active") == true || volume_name.Contains("Passive") == true) &&
           volume_name.Contains("end") == false;
  } else {
    return volume_name.Contains("volECAL") == true &&
           volume_name.Contains("Active") == true &&
           volume_name.Contains("end") == false;
  }
}

bool SANDGeoManager::is_ecal_endcap(const TString& volume_name, bool include_passive) const
{
  // something like: endvolECALActiveSlab_0_PV_0
  if (include_passive) {
    return volume_name.Contains("endvolECAL") == true &&
           (volume_name.Contains("Active") == true || volume_name.Contains("Passive") == true);
  } else {
    return volume_name.Contains("endvolECAL") == true &&
           volume_name.Contains("Active") == true;
  }
}

bool SANDGeoManager::is_endcap_mod(const TString& volume_name) const
{
  return volume_name.Contains(endcap_mod_regex_);
}

bool SANDGeoManager::check_and_process_ecal_path(TString& volume_path) const
{

  // BARREL ==> something like:
  // "/volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_lv_PV_18/volECALActiveSlab_21_PV_0"
  // ENDCAP ==> something like:
  // "/volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_end_lv_PV_0/ECAL_ec_mod_4_lv_PV_1/endvolECALActiveSlab_0_PV_0"
  TObjArray* obj = volume_path.Tokenize("/");

  int size = obj->GetEntries();
  if (size < 8) {
    return false;
  };

  // BARREL => ECAL_lv_PV_18
  // ENDCAP => ECAL_end_lv_PV_0/ECAL_ec_mod_4_lv_PV_1
  if (size == 8) {  // barrel module path
    volume_path = ((TObjString*)obj->At(6))->GetString();
  } else if (size == 10) {  // endcap module path
    volume_path = ((TObjString*)obj->At(6))->GetString() + "_" +
                  ((TObjString*)obj->At(7))->GetString() + "_" +
                  ((TObjString*)obj->At(8))->GetString();
  } else
    return false;
  delete obj;

  return true;
}

void SANDGeoManager::get_ecal_barrel_module_and_layer(
    const TString& volume_name, const TString& volume_path, int& detector_id,
    int& module_id, int& layer_id) const
{
  TObjArray* obja1 =
      volume_name.Tokenize("_");  // BARREL => volECALActiveSlab_21_PV_0
  TObjArray* obja2 = volume_path.Tokenize("_");  // BARREL => ECAL_lv_PV_18

  // top module => modID == 0
  // increasing modID counterclockwise as seen from positive x
  //(i.e. z(modID==1) < z(modID==0) & z(modID==0) < z(modID==23))
  detector_id = 2;
  module_id = ((TObjString*)obja2->At(3))->GetString().Atoi();
  int slab_id = ((TObjString*)obja1->At(1))->GetString().Atoi();  // 21

  delete obja1;
  delete obja2;

  // layer_id==0 -> smallest slab -> internal (slab_id==0 ?)
  // layer_id==208 -> biggest slab -> external (slab_id==208 ?)
  // why layer_id = slab_id / 40? Because const int number_of_layers = 5 =
  // 208//40
  layer_id = slab_id / 40;

  if (layer_id > 4) layer_id = 4;
}


// Updated version
void SANDGeoManager::get_ecal_endcap_module_and_layer(
    const TString& volume_name, const TString& volume_path, int& detector_id,
    int& module_id, int& layer_id) const
{
  TObjArray* obja1 =
      volume_name.Tokenize("_");  // ENDCAP => endvolECALActiveSlab_0_PV_0
  TObjArray* obja2 = volume_path.Tokenize(
      "_");  // ENDCAP =>
             // ECAL_end_lv_PV_0/ECAL_ec_mod_4_lv_PV_1/ECAL_ec_mod_vert_0_lv_PV

  // std::cout << "volume_name: " << volume_name << "\n";
  // for (int i = 0; i < obja1->GetEntries(); i++)
  //   std::cout << ((TObjString*)obja1->At(i))->GetString().Atoi() << "\n";

  int slab_id = 0;
  if (volume_name.Contains("curv"))
    slab_id = ((TObjString*)obja1->At(2))->GetString().Atoi();
  else
    slab_id = ((TObjString*)obja1->At(3))->GetString().Atoi();

  // std::cout << "> check side_id: " <<
  // ((TObjString*)obja2->At(4))->GetString()
  //           << "\n";
  detector_id = ((TObjString*)obja2->At(4))->GetString().Atoi();
  int mod_id = ((TObjString*)obja2->At(8))->GetString().Atoi();
  int replica_id = ((TObjString*)obja2->At(11))->GetString().Atoi();

  module_id = encode_endcap_mod_id(mod_id, replica_id, detector_id);

  // // mod == 40 -> left  -> detID = 1
  // // mod == 30 -> right -> detID = 3
  // // (see issue: https://baltig.infn.it/dune/sand-reco/-/issues/18)
  // if (module_id == 0) {
  //   detector_id = 1;
  //   module_id = 40;
  // } else if (module_id == 1) {
  //   detector_id = 3;
  //   module_id = 30;
  // }
  // delete obja1;
  // delete obja2;

  // layer_id==0 -> internal (slab_id==0 ?)
  // layer_id==208 -> external (slab_id==208 ?)
  // layer numbers grow in the opposite direction in vertical modules
  layer_id = slab_id / 40;

  if (layer_id > 4) layer_id = 4;

  delete obja1;
  delete obja2;
}

void SANDGeoManager::get_ecal_barrel_cell_local_id(double x, double y, double z,
                                                   const TGeoNode* const node,
                                                   int& cell_local_id) const
{
  double master[3];
  double local[3];
  master[0] = x;
  master[1] = y;
  master[2] = z;

  geo_->GetCurrentNavigator()->MasterToLocal(master, local);

  TString shape_name = node->GetVolume()->GetShape()->GetName();

  if (shape_name != "TGeoTrd2") {
    std::cout << __FILE__ << " " << __LINE__ << "\n";
    std::cout << "invalid shape : " << shape_name << "\n";
    throw "";
  }

  TGeoTrd2* trd = (TGeoTrd2*)node->GetVolume()->GetShape();

  double dx1 = trd->GetDx1();
  double dx2 = trd->GetDx2();
  double dz = trd->GetDz();

  // http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/html/Detector/Geometry/geomSolids.html
  // if z = -dz -> dx = 2*dx1
  // if z =  dz -> dx = 2*dx2
  // semilarghezza della slab di scintillatore alla quota Plocal[2]
  double dx = 0.5 * local[2] / dz * (dx2 - dx1) + 0.5 * (dx2 + dx1);

  // Cell width at z = Plocal[2]
  double cell_width =
      2. * dx / sand_geometry::ecal::number_of_cells_per_barrel_layer;

  // cellID = distanza dall'estremo diviso larghezza cella
  cell_local_id = (local[0] + dx) / cell_width;
  // set the cellID to 11 if the hit end-up at the outer boundary
  if(cell_local_id == sand_geometry::ecal::number_of_cells_per_barrel_layer)
    cell_local_id--;
}

int SANDGeoManager::get_barrel_path_len(const double& hx, const double& hy,
                                        const double& hz, double& d1,
                                        double& d2) const
{
  double master[3];
  double local[3];
  master[0] = hx;
  master[1] = hy;
  master[2] = hz;

  TGeoNode* layer_node = geo_->FindNode(hx, hy, hz);
  if (layer_node == 0) return -999;

  TGeoTrd2* trd = (TGeoTrd2*)layer_node->GetVolume()->GetShape();
  geo_->GetCurrentNavigator()->MasterToLocal(master, local);

  d1 = trd->GetDy1() - local[1];
  d2 = trd->GetDy1() + local[1];

  // std::cout << "Barrel layer: " << layer_node->GetName() << "\nd1: " << d1
  //           << ", d2: " << d2 << "\n";

  return 1;
}

int SANDGeoManager::get_barrel_hit_pos(double d1,
                                       int global_cellID, double& reco_x,
                                       double& reco_y, double& reco_z) const
{
  auto current_cell = cellmap_.at(global_cellID);
  double master[3];
  double local[3];
  master[0] = current_cell.getX();
  master[1] = current_cell.getY();
  master[2] = current_cell.getZ();

  TGeoNode* layer_node = geo_->FindNode(master[0], master[1], master[2]);
  if (layer_node == 0) return -999;

  TGeoTrd2* trd = (TGeoTrd2*)layer_node->GetVolume()->GetShape();
  geo_->GetCurrentNavigator()->MasterToLocal(master, local);

  local[1] = trd->GetDy1() - d1;

  geo_->GetCurrentNavigator()->LocalToMaster(local, master);

  reco_x = master[0];
  reco_y = current_cell.getY();
  reco_z = current_cell.getZ();
  return 1;
}

// NEW VERSION --> This needs to be reviewed
void SANDGeoManager::get_ecal_endcap_cell_local_id(double x, double y, double z,
                                                   const int& endcap_mod_id,
                                                   int& cell_local_id) const
{
  double master[3];
  double local[3];
  master[0] = x;
  master[1] = y;
  master[2] = z;

  // geo_->GetCurrentNavigator()->MasterToLocal(master, local);
  // TGeoTube* tub = (TGeoTube*)node->GetVolume()->GetShape();
  endcapmap_.at(endcap_mod_id).getModHMatrix().MasterToLocal(master, local);

  // double rmin = tub->GetRmin();
  // double rmax = tub->GetRmax();
  // double dz = tub->GetDz();

  // Cell width at z = Plocal[2]
  // double cell_width = endcapmap_.at(endcap_mod_id).width() /
  //                     sand_geometry::ecal::endcap_cell_width;

  // std::cout << "endcapmap_[" << endcap_mod_id
  //           << "].width(): " << endcapmap_.at(endcap_mod_id).width()
  //           << ", cell_width: " << cell_width << "\n";
  // cellID = distanza dall'estremo diviso larghezza cella
  cell_local_id = (local[0] + 0.5 * endcapmap_.at(endcap_mod_id).getWidth()) /
                  sand_geometry::ecal::endcap_cell_width;

  // std::cout << "local[0]: " << local[0]
  //           << ", Dx: " << local[0] + 0.5 *
  //           endcapmap_.at(endcap_mod_id).width()
  //           << "\n";
}

int SANDGeoManager::get_endcap_path_len(const double& hx, const double& hy,
                                        const double& hz,
                                        const int& endcap_mod_id, double& d1,
                                        double& d2) const
{
  double master[3];
  double local[3];
  master[0] = hx;
  master[1] = hy;
  master[2] = hz;

  // extract the module corresponding to the indxex from the map
  auto ec_mod = endcapmap_.at(endcap_mod_id);

  TGeoNode* layer_node = geo_->FindNode(hx, hy, hz);
  if (layer_node == 0) return -999;

  TString volume_name = layer_node->GetName();
  TString volume_path = geo_->GetPath();

  // check whether the layer is actually contained inside the
  // module

  if (!volume_path.Contains(ec_mod.getPath())) return -999;

  // convert to the section local coordinates (one level up)
  geo_->GetCurrentNavigator()->CdUp();
  geo_->MasterToLocal(master, local);

  // manage each section separately
  if (volume_path.Contains("vert")) {
    auto depth = local[2] + ec_mod.getModDz();
    d1 = 0.5 * ec_mod.getLVert() - local[1] + ec_mod.getCurvatureArcLength(depth) +
         ec_mod.getLHor();
    d2 = ec_mod.getCellTotalLength(depth) - d1;

  } else if (volume_path.Contains("hor") &&
             volume_path.Contains("lv_PV_0/endvol")) {
    auto depth = ec_mod.getModDz() - local[2];
    d1 = 0.5 * ec_mod.getLHor() + local[1];
    d2 = ec_mod.getCellTotalLength(depth) - d1;
  } else if (volume_path.Contains("hor") &&
             volume_path.Contains("lv_PV_1/endvol")) {
    auto depth = ec_mod.getModDz() - local[2];
    d2 = 0.5 * ec_mod.getLHor() - local[1];
    d1 = ec_mod.getCellTotalLength(depth) - d2;
  } else if (volume_path.Contains("curv") &&
             volume_path.Contains("lv_PV_0/endvol")) {
    auto depth = ec_mod.getRMax() -
                 std::sqrt(std::pow(local[0], 2) + std::pow(local[1], 2));
    d1 = std::sqrt(std::pow(local[0], 2) + std::pow(local[1], 2)) *
             std::atan(std::abs(local[0] / local[1])) +
         ec_mod.getLHor();
    d2 = ec_mod.getCellTotalLength(depth) - d1;

  } else if (volume_path.Contains("curv") &&
             volume_path.Contains("lv_PV_1/endvol")) {
    auto depth = ec_mod.getRMax() -
                 std::sqrt(std::pow(local[0], 2) + std::pow(local[1], 2));
    d2 = std::sqrt(std::pow(local[0], 2) + std::pow(local[1], 2)) *
             std::atan(std::abs(local[0] / local[1])) +
         ec_mod.getLHor();
    d1 = ec_mod.getCellTotalLength(depth) - d2;
  } else {
    return 0;
  }
  return 1;
}

int SANDGeoManager::get_endcap_hit_pos(const double& d1,
                                       const int& global_cellID,
                                       const int& modID, double& reco_x,
                                       double& reco_y, double& reco_z) const
{
  // extract the cell and module corresponding to the indexes from the
  // corresponding maps
  auto current_cell = cellmap_.at(global_cellID);
  auto ec_mod = endcapmap_.at(modID);

  double master[3];
  double local[3];
  master[0] = current_cell.getX();
  master[1] = current_cell.getY();
  master[2] = current_cell.getZ();

  TGeoNode* layer_node = geo_->FindNode(master[0], master[1], master[2]);
  if (layer_node == 0) return -999;

  TString volume_name = layer_node->GetName();
  TString volume_path = geo_->GetPath();

  // check whether the layer is actually contained inside the
  // module
  if (!volume_path.Contains(ec_mod.getPath())) return -999;

  // convert to the section local coordinates (one level up)
  geo_->GetCurrentNavigator()->CdUp();
  geo_->MasterToLocal(master, local);

  // the cell center will be in the vertical section
  auto depth = local[2] + ec_mod.getModDz();
  auto cell_rad = (ec_mod.getRMax() - depth);
  const double d_hor0 = ec_mod.getLHor(),
               d_curv0 = d_hor0 + ec_mod.getCurvatureArcLength(depth),
               d_vert = d_curv0 + ec_mod.getLVert(),
               d_curv1 = d_vert + ec_mod.getCurvatureArcLength(depth),
               d_hor1 = d_curv1 + ec_mod.getLHor();
  // the local coordinates will always refer to the vertical section 
  // find the right module section based on the d1 range
  if (d1 <= d_hor0) {
    local[1] = 0.5 * ec_mod.getLVert() + ec_mod.getRMax() - depth;
    local[2] = -ec_mod.getModDz() + ec_mod.getRMax() + (d_hor0 - d1);
  } else if (d1 > d_hor0 && d1 <= d_curv0) {
    const auto sec_angle = (d1 - d_hor0) / cell_rad;
    local[1] = 0.5 * ec_mod.getLVert() + cell_rad * std::cos(sec_angle);
    local[2] = -ec_mod.getModDz() + ec_mod.getRMax() - cell_rad * std::sin(sec_angle);
  } else if (d1 > d_curv0 && d1 <= d_vert) {
    local[1] = 0.5 * ec_mod.getLVert() - (d1 - ec_mod.getCurvatureArcLength(depth) - d_hor0);
  } else if (d1 > d_vert && d1 <= d_curv1) {
    const auto sec_angle = (d_curv1 - d1) / cell_rad;
    local[1] = -0.5 * ec_mod.getLVert() - cell_rad * std::cos(sec_angle);
    local[2] = -ec_mod.getModDz() + ec_mod.getRMax() - cell_rad * std::sin(sec_angle);
  } else if (d1 > d_curv1) {
    if (ec_mod.getNSections() == 5 && d1 <= d_hor1) {
      local[1] = -0.5 * ec_mod.getLVert() - ec_mod.getRMax() + depth;
      local[2] = -ec_mod.getModDz() + ec_mod.getRMax() + ec_mod.getLHor() - (d_hor1 - d1);
    } else {
      return -999;
    }
  }

  geo_->LocalToMaster(local, master);
  reco_x = master[0];
  reco_y = master[1];
  reco_z = master[2];
  return 1;
}

int SANDGeoManager::get_hit_path_len(const double& hx, const double& hy,
                                     const double& hz,
                                     const int& global_cell_id, double& d1,
                                     double& d2) const
{
  if (geo_ == 0) {
    std::cout << "ERROR: TGeoManager pointer not initialized" << std::endl;
    return -999;
  }

  // decode the global_cell_id to extract detID and modID
  int detID, modID, layerID, locID;
  decode_ecal_cell_id(global_cell_id, detID, modID, layerID, locID);

  int exit = 0;
  if (detID == 2) {  // barrel modules
    exit = get_barrel_path_len(hx, hy, hz, d1, d2);
  } else if (detID == 0 || detID == 1) {  // endcap modules
    exit = get_endcap_path_len(hx, hy, hz, modID, d1, d2);
  } else {
    std::cout << "> get_hit_path_len exiting with error:\n";
    return -999;
  }
  return exit;
}

double SANDGeoManager::compute_cell_d1(const double& cell_l,
                                       const double& tdc_1,
                                       const double& tdc_2) const
{
  return 0.5 *
         (cell_l + ((tdc_1 - tdc_2) / (sand_reco::ecal::scintillation::vlfb *
                                       sand_reco::conversion::mm_to_m)));
  // either this or 0.5*cell_l+sand_reco::ecal::reco::XfromTDC(double t1, double t2)
}

double SANDGeoManager::compute_cell_d2(const double& cell_l,
                                       const double& tdc_1,
                                       const double& tdc_2) const
{
  return 0.5 *
         (cell_l - ((tdc_1 - tdc_2) / (sand_reco::ecal::scintillation::vlfb *
                                       sand_reco::conversion::mm_to_m)));
  // either this or 0.5*cell_l-sand_reco::ecal::reco::XfromTDC(double t1, double t2)
}

int SANDGeoManager::get_reco_hit_pos(const int& cellID, const double& cell_l,
                                     const double& tdc_1, const double& tdc_2,
                                     double& reco_x, double& reco_y,
                                     double& reco_z) const
{

  if (geo_ == 0) {
    std::cout << "ERROR: TGeoManager pointer not initialized" << std::endl;
    return -999;
  }

  const double d1 = compute_cell_d1(cell_l, tdc_1, tdc_2);

  if (d1 < 0 || d1 > cell_l) {
    return -999;
  }

  // decode the global_cell_id to extract detID and modID
  int detID, modID, layerID, locID;
  decode_ecal_cell_id(cellID, detID, modID, layerID, locID);

  int exit = 0;
  if (detID == 2) {  // barrel modules
    exit = get_barrel_hit_pos(d1, cellID, reco_x, reco_y, reco_z);
  } else if (detID == 0 || detID == 1) {  // endcap modules
     exit = get_endcap_hit_pos(d1, cellID, modID, reco_x, reco_y, reco_z);
   } else {
    std::cout << "> get_hit_path_len exiting with error:\n";
    return -999;
  }
  return exit;
}

// -- -- -- -- -- -- -- -- -- -- --

void SANDGeoManager::set_ecal_info()
{
  // https://root.cern.ch/root/htmldoc/guides/users-guide/Geometry.html#shapes
  // GetDx1() half length in x at -Dz
  // GetDx2() half length in x at +Dz
  // Dx1 < Dx2 => -Dz corresponds to minor width => internal side
  //-- TGeoTrd2 methods, which defines a trapezoid with up and down faces
  // parallel to zy and oblique sides along xz

  // this is potentially outdated: barrel-module volumes contain both trapezoid
  // layers and the Al plate for dz the old code should be valid: 127.5 is the
  // same Dz as the endcaps (which have the same layers) I'm not sure about
  // xmax, though
  TGeoTrd2* mod =
      (TGeoTrd2*)geo_->FindVolumeFast(sand_geometry::ecal::barrel_module_name.c_str())
          ->GetShape();
  double ecal_barrel_xmin = mod->GetDx1();

  double ecal_barrel_dz = mod->GetDz();
  double ecal_barrel_dy = mod->GetDy1();
  TGeoTrd2* last_passive_slab =
      (TGeoTrd2*)geo_
          ->FindVolumeFast(sand_geometry::ecal::barrel_last_passive_slab_name.c_str())
          ->GetShape();
  double ecal_barrel_xmax = last_passive_slab->GetDx2();

  // This is outdated
  // TGeoTube* ec =
  //     (TGeoTube*)geo_->FindVolumeFast(sand_geometry::ecal::endcap_module_name)
  //         ->GetShape();
  // double ecal_endcap_rmax = ec->GetRmax();  // Maximum radius = 2000
  // double ecal_endcap_rmin = ec->GetRmin();

  // get z of the levels between the layers
  auto z_levels =
      get_levels_z(ecal_barrel_dz, sand_geometry::ecal::layer_thickness);

  // get slope of the edge of the barrel module
  //-- the height of the module is 2 * ecal_barrel_dz, hence 0.5
  auto ecal_barrel_edge_slope =
      0.5 * (ecal_barrel_xmax - ecal_barrel_xmin) / ecal_barrel_dz;
  auto ecal_barrel_edge_position = 0.5 * (ecal_barrel_xmax + ecal_barrel_xmin);

  // eval barrel cell (within each module) center local position (relative to
  // the module) this is an std::map<int,TVector3>
  auto ecal_barrel_cell_center_local_positions =
      get_ecal_barrel_cell_center_local_position(
          z_levels, ecal_barrel_edge_slope, ecal_barrel_edge_position);

  double local[3];
  double master[3];

  for (int module_id = 0;
       module_id < sand_geometry::ecal::number_of_barrel_modules; module_id++) {
    for (auto cell_position : ecal_barrel_cell_center_local_positions) {
      // first/second stands for the key/value in a map
      local[0] = cell_position.second.X();
      local[1] = cell_position.second.Y();
      local[2] = cell_position.second.Z();

      auto cell_and_layer_id =
          decode_ecal_barrel_cell_local_id(cell_position.first);
      auto layer_id = cell_and_layer_id.first;
      auto cell_local_id = cell_and_layer_id.second;

      // module_id is just an integer from o to n_modules
      // local and master are the coordinate vectors
      geo_->cd(
          TString::Format(sand_geometry::ecal::path_barrel_template.c_str(), module_id)
              .Data());

      geo_->LocalToMaster(local, master);

      // here we create new cellInfo
      int detector_id = 2;
      int cell_unique_id =
          encode_ecal_cell_id(detector_id, module_id, layer_id, cell_local_id);
      cellmap_[cell_unique_id] = sand_geometry::ecal::ECALCellInfo(
          cell_unique_id, master[0], master[1], master[2], 2 * ecal_barrel_dy,
          sand_geometry::ecal::ECALCellInfo::ModuleType::kBarrel);
    }
  }
  std::cout << "> Barrel cells info. set\n";

  // fill the endcap module info (unique_mod_id, width, n_sections)
  set_ecal_endcap_info();

  // module thickness along z is the same for all modules and cells (so use the
  // first item)
  z_levels = get_levels_z(endcapmap_.begin()->second.getModDz(),
                          sand_geometry::ecal::ec_layer_thickness);

  for (const auto& module : endcapmap_) {
    auto endcap_cell_center_local_positions =
        get_ec_cell_center_local_position(z_levels, module.second);

    for (auto cell_position : endcap_cell_center_local_positions) {
      local[0] = cell_position.second.X();
      local[1] = cell_position.second.Y();
      local[2] = cell_position.second.Z();

      double cell_length =
          module.second.getCellTotalLength(local[2] + module.second.getModDz());

      auto cell_and_layer_id =
          decode_ecal_endcap_cell_local_id(cell_position.first);
      auto layer_id = cell_and_layer_id.first;
      auto cell_local_id = cell_and_layer_id.second;

      int detector_id = 0;
      int replica_id = 0;
      int m_ID = 0;
      decode_endcap_mod_id(module.first, m_ID, replica_id, detector_id);

      // cd to the module path
      geo_->cd(module.second.getPath());
      geo_->LocalToMaster(local, master);

      // here we create new cellInfo
      int cell_unique_id = encode_ecal_cell_id(detector_id, module.first,
                                               layer_id, cell_local_id);

      cellmap_[cell_unique_id] =
          sand_geometry::ecal::ECALCellInfo(cell_unique_id, master[0], master[1], master[2],
                           cell_length, sand_geometry::ecal::ECALCellInfo::ModuleType::kEndcap);
    }
  }
  std::cout << "> Endcap cells info. set\n";
  std::cout << "> cellmap_ size: " << cellmap_.size() << "\n";
}

//  #######################################################
//  ##                     ECAL_ENDCAP                   ##
//  #######################################################
int SANDGeoManager::encode_endcap_mod_id(int module_id, int module_replica_id,
                                         int endcap_side_id)
{
  return module_id * 100 + module_replica_id * 10 + endcap_side_id;
}

void SANDGeoManager::decode_endcap_mod_id(int endcap_mod_global_id,
                                          int& module_id,
                                          int& module_replica_id,
                                          int& endcap_side_id)
{
  module_id = endcap_mod_global_id / 100;
  module_replica_id = (endcap_mod_global_id - module_id * 100) / 10;
  endcap_side_id = endcap_mod_global_id % 10;
}

int SANDGeoManager::get_endcap_mod_id(const TString& volume_path) const
{
  auto module_matches = endcap_mod_path_regex_.MatchS(volume_path);

  // endcap side matching
  int endcap_side_id = (reinterpret_cast<TObjString*>(module_matches->At(1)))
                           ->GetString()
                           .Atoi();
  int module_id = (reinterpret_cast<TObjString*>(module_matches->At(2)))
                      ->GetString()
                      .Atoi();
  int module_replica_id = (reinterpret_cast<TObjString*>(module_matches->At(3)))
                              ->GetString()
                              .Atoi();
  delete module_matches;

  return encode_endcap_mod_id(module_id, module_replica_id, endcap_side_id);
}
void SANDGeoManager::set_ecal_endcap_info(const TGeoHMatrix& matrix)
{
  TGeoNode* node = gGeoManager->GetCurrentNode();
  TString node_path = gGeoManager->GetPath();
  TString node_name = node->GetName();
  TGeoMatrix* node_matrix = node->GetMatrix();
  TGeoHMatrix node_hmatrix = matrix * (*node_matrix);

  if (is_endcap_mod(node_name)) {

    int mod_id = get_endcap_mod_id(node_path);
    // set the module info
    endcapmap_[mod_id] = sand_geometry::ecal::ENDCAPModInfo(mod_id, node, node_hmatrix);


  } else {
    for (int i = 0; i < node->GetNdaughters(); i++) {
      gGeoManager->CdDown(i);
      set_ecal_endcap_info(node_hmatrix);
      gGeoManager->CdUp();
    }
  }
}

void SANDGeoManager::set_ecal_endcap_info()
{
  geo_->CdTop();
  // std::cout << geo_->GetPath() << "\n";
  TGeoHMatrix matrix = *gGeoIdentity;
  std::cout << "> Checking endcap info\n";
  set_ecal_endcap_info(matrix);
}

//  #######################################################
//  ##                     STT ENCODING                  ##
//  #######################################################
sand_geometry::tracker::plane_iterator SANDGeoManager::getPlaneInfo(sand_geometry::tracker::PlaneID plane_global_id) const
{
  sand_geometry::tracker::ModuleID module_unique_id;
  sand_geometry::tracker::PlaneID  plane_local_id, plane_type;

  decodePlaneId(plane_global_id, module_unique_id, 
                  plane_local_id, plane_type);
  return id_to_plane_.at(plane_global_id);
}

sand_geometry::tracker::plane_iterator SANDGeoManager::getPlaneInfo(sand_geometry::tracker::CellID cell_global_id) const
{
  sand_geometry::tracker::PlaneID  plane_global_id;
  sand_geometry::tracker::CellID   cell_local_id;

  decodeCellId(cell_global_id, plane_global_id, cell_local_id);
  return id_to_plane_.at(plane_global_id);
}

std::map<sand_geometry::tracker::CellID, sand_geometry::tracker::Cell>::const_iterator SANDGeoManager::getCellInfo(sand_geometry::tracker::CellID cell_global_id) const
{
  sand_geometry::tracker::ModuleID module_unique_id;
  sand_geometry::tracker::PlaneID  plane_global_id, plane_local_id, plane_type;
  sand_geometry::tracker::CellID   cell_local_id;

  decodeCellId(cell_global_id, plane_global_id, cell_local_id);
  decodePlaneId(plane_global_id, module_unique_id, 
                  plane_local_id, plane_type);
  return id_to_plane_.at(plane_global_id)->getCell(cell_global_id);
}

sand_geometry::tracker::CellID SANDGeoManager::encodeCellId(sand_geometry::tracker::PlaneID plane_global_id, sand_geometry::tracker::CellID cell_local_id)
{
  return plane_global_id() * 10000 + cell_local_id();
}

void SANDGeoManager::decodeCellId(sand_geometry::tracker::CellID   cell_global_id, 
                                    sand_geometry::tracker::PlaneID& plane_global_id,
                                    sand_geometry::tracker::CellID&  cell_local_id)
{
  plane_global_id = cell_global_id() / 10000;  // global id
  cell_local_id = cell_global_id() % 10000;
}

sand_geometry::tracker::PlaneID SANDGeoManager::encodePlaneId(sand_geometry::tracker::ModuleID unique_module_id,
                                     sand_geometry::tracker::PlaneID plane_replica_id, sand_geometry::tracker::PlaneID plane_type)
{
  return unique_module_id() + (2 * plane_replica_id() + plane_type()) * 10 + plane_type();
}

void SANDGeoManager::decodePlaneId(sand_geometry::tracker::PlaneID plane_global_id, sand_geometry::tracker::ModuleID& unique_module_id, 
                                     sand_geometry::tracker::PlaneID& plane_replica_id, sand_geometry::tracker::PlaneID& plane_type)
{
  unique_module_id = plane_global_id() / 100 * 100;
  sand_geometry::tracker::PlaneID local_plane_id = plane_global_id() - unique_module_id();
  plane_type = local_plane_id() % 10;
  plane_replica_id = ((local_plane_id() / 10) - plane_type()) / 2;
}

sand_geometry::tracker::ModuleID SANDGeoManager::encodeModuleId(sand_geometry::tracker::ModuleID supermodule_id, sand_geometry::tracker::ModuleID module_id, sand_geometry::tracker::ModuleID module_replica_id)
{
  return supermodule_id() * 1E5 + (module_id() * 10 + module_replica_id()) * 100;
}

void SANDGeoManager::decodeModuleId(sand_geometry::tracker::ModuleID unique_module_id, sand_geometry::tracker::ModuleID& supermodule_id, 
                                      sand_geometry::tracker::ModuleID& module_id, sand_geometry::tracker::ModuleID& module_replica_id)
{
  supermodule_id = unique_module_id() / 1E5;
  sand_geometry::tracker::ModuleID local_module_id = (unique_module_id() - supermodule_id() * 1E5) / 100;
  module_id = local_module_id() / 10;
  module_replica_id = local_module_id() % 10;
}
bool SANDGeoManager::isSttTube(const TString& volume_name) const
{
  return volume_name.Contains(stt_tube_regex_);
}

bool SANDGeoManager::isSttPlane(const TString& volume_name) const
{
  return volume_name.Contains(stt_plane_regex_);
}
bool SANDGeoManager::isDriftPlane(const TString& volume_name) const
{
  return volume_name.Contains(drift_plane_regex_);
}

sand_geometry::tracker::ModuleID SANDGeoManager::getSttModuleId(const TString& volume_path) const
{
  auto supermodule_matches = stt_supermodule_regex_.MatchS(volume_path);

  long supermodule_id;
  if (supermodule_matches->GetEntries() == 0) {
    supermodule_id = 0;
  } else {
    // To Do: Currently there are no supermodules in the stt geometry
  }

  auto module_matches = stt_module_regex_.MatchS(volume_path);
  long module_id =
      (reinterpret_cast<TObjString*>(module_matches->At(2)))->GetString().Atoi();
  long module_replica_id = (reinterpret_cast<TObjString*>(module_matches->At(3)))
                            ->GetString()
                            .Atoi();
  return encodeModuleId(sand_geometry::tracker::ModuleID(supermodule_id), 
                          sand_geometry::tracker::ModuleID(module_id), 
                          sand_geometry::tracker::ModuleID(module_replica_id));
}

sand_geometry::tracker::PlaneID SANDGeoManager::getSttPlaneId(const TString& volume_path, bool justLocal = false) const
{

  auto plane_matches = stt_plane_regex_.MatchS(volume_path);

  if (plane_matches->GetEntries() < 5) {
    // Sometimes the volume path returned by the TGeoManager does not match
    // with the expected one for a tube...to be investigated!!!
    // std::cout << "Error: volume path for STT digit not expected!! returning
    // default value (0) for stt plane id" << std::endl;
    delete plane_matches;
    return -999;
  }

  int plane_replica_id =
      (reinterpret_cast<TObjString*>(plane_matches->At(4)))->GetString().Atoi();

  int plane_type = (reinterpret_cast<TObjString*>(plane_matches->At(3)))
                           ->GetString()
                           .EqualTo("XX")
                       ? 2
                       : 1;
delete plane_matches;
  if (justLocal) {
    return plane_type;
  } else {
    sand_geometry::tracker::ModuleID unique_module_id = getSttModuleId(volume_path);
    return encodePlaneId(unique_module_id, sand_geometry::tracker::PlaneID(plane_replica_id), sand_geometry::tracker::PlaneID(plane_type));
  }
}

sand_geometry::tracker::ModuleID SANDGeoManager::getDriftSupermoduleId(const TString& volume_path) const
{
  // upstram -> downstrea,
  // Trk, C1, B1, A1, A0, B0, C0, X0, X1
  //   0,  1,  2,  3,  4,  5,  6,  7,  8
  auto supermodule_matches = supermodule_regex_.MatchS(volume_path);
  TString supermodule_name =
      (reinterpret_cast<TObjString*>(supermodule_matches->At(2)))->GetString();
  int supermodule_replica =
      (reinterpret_cast<TObjString*>(supermodule_matches->At(3)))
          ->GetString()
          .Atoi();
  int supermodule_id;
  if (supermodule_name.Contains("X0")) {
    supermodule_id = 8;
  } else if (supermodule_name.Contains("X1")) {
    supermodule_id = 7;
  } else if (supermodule_name.Contains("C")) {
    supermodule_id = supermodule_replica ? 1 : 6;
  } else if (supermodule_name.Contains("B")) {
    supermodule_id = supermodule_replica ? 2 : 5;
  } else if (supermodule_name.Contains("A")) {
    supermodule_id = supermodule_replica ? 3 : 4;
  } else {
    supermodule_id = 0;
  }

  delete supermodule_matches;

  return sand_geometry::tracker::ModuleID(supermodule_id);
}

sand_geometry::tracker::ModuleID SANDGeoManager::getDriftModuleReplicaId(const TString& volume_path)
    const
{
  auto matches = module_regex_.MatchS(volume_path);
  auto type = (reinterpret_cast<TObjString*>(matches->At(1)))->GetString();
  long id = (reinterpret_cast<TObjString*>(matches->At(3)))->GetString().Atoi();
  if (type == "C") {
    id = 9;
  }
  return sand_geometry::tracker::ModuleID(id);
}

bool SANDGeoManager::isSwire(const TString& volume_path) const
{
  return (volume_path.Contains("Swire"));
}

sand_geometry::tracker::ModuleID SANDGeoManager::getDriftModuleId(const TString& volume_path) const
{
  sand_geometry::tracker::ModuleID supermodule_id(getDriftSupermoduleId(volume_path));
  sand_geometry::tracker::ModuleID module_id(0);
  sand_geometry::tracker::ModuleID module_replica_id(0);
  if (supermodule_id() != 0) {
    module_replica_id = getDriftModuleReplicaId(volume_path);
  }
  return encodeModuleId(supermodule_id, module_id, module_replica_id);
}

sand_geometry::tracker::PlaneID SANDGeoManager::getDriftPlaneId(const TString& volume_path,
                                       bool JustLocalId = false) const
{
  auto plane_matches = drift_plane_regex_.MatchS(volume_path);

  if (plane_matches->GetEntries() == 0) {
    delete plane_matches;
    return 0;
  }

  int plane_type =
      (reinterpret_cast<TObjString*>(plane_matches->At(2)))->GetString().Atoi();
  int plane_replica_id =
      (reinterpret_cast<TObjString*>(plane_matches->At(4)))->GetString().Atoi();

  if (JustLocalId) {
    return plane_type;
  } else {
    sand_geometry::tracker::ModuleID unique_module_id(getDriftModuleId(volume_path));
    return encodePlaneId(unique_module_id, plane_replica_id, plane_type);
  }
}

std::vector<TVector2> SANDGeoManager::getLocalLinePlaneIntersections(const TVector2& local_2d_position,
                                               const sand_geometry::tracker::Plane& plane)
{
  TVector2 local_plane_x_axis = rotatedToLocal(TVector2(1, 0), plane);
  std::vector<TVector2> vertices = plane.getPlaneVertices();

  std::vector<TVector2> intersection_points;

  TVector2 intersection(0, 0);
    if (getLineSegmentIntersection(local_2d_position, local_plane_x_axis, 
                               vertices[0], vertices[1], intersection)) {
      intersection_points.push_back(intersection);
    };
    if (getLineSegmentIntersection(local_2d_position, local_plane_x_axis, 
                               vertices[1], vertices[2], intersection)) {
      intersection_points.push_back(intersection);
    };
    if (getLineSegmentIntersection(local_2d_position, local_plane_x_axis, 
                               vertices[2], vertices[3], intersection)) {
      intersection_points.push_back(intersection);
    };
    if (getLineSegmentIntersection(local_2d_position, local_plane_x_axis, 
                               vertices[3], vertices[0], intersection)) {
      intersection_points.push_back(intersection);
    };
  return intersection_points;
}

std::vector<TVector2> SANDGeoManager::getGlobalLinePlaneIntersections(const TVector2& local_2d_position,
                                               const sand_geometry::tracker::Plane& plane)
{
  std::vector<TVector2> intersection_points = getLocalLinePlaneIntersections(local_2d_position, plane);

  for(auto& intersection:intersection_points) {
    intersection = localToGlobal(intersection, plane);
  }

  return intersection_points;
}

void SANDGeoManager::setSttPlaneInfo(const TGeoNode* const node,
                                        const TGeoHMatrix& matrix)
{
  TString node_path = gGeoManager->GetPath();
  sand_geometry::tracker::PlaneID stt_plane_unique_id = getSttPlaneId(node_path);
  sand_geometry::tracker::PlaneID stt_plane_local_id  = getSttPlaneId(node_path, true);

  planes_.push_back(sand_geometry::tracker::Plane(stt_plane_unique_id, stt_plane_local_id));
  id_to_plane_[planes_.back().uId()] = std::prev(planes_.end());

  auto& plane = planes_.back();
  double angle = tracker_module_configuration::stt::id_to_angle[std::to_string(stt_plane_local_id())];

  plane.setRotation(angle);
  
  TGeoMatrix* plane_matrix = node->GetMatrix();
  TGeoHMatrix plane_hmatrix = matrix * (*plane_matrix);
  TGeoBBox* plane_shape = (TGeoBBox*)node->GetVolume()->GetShape();
  TVector3 plane_dimension;
  // Notice: this is a workaround to the planes in the geometry being rotated sometimes
  if(stt_plane_local_id() == 2) {
    plane_dimension.SetX(2 * plane_shape->GetDZ());
    plane_dimension.SetY(2 * plane_shape->GetDY());
  } else {
    plane_dimension.SetY(2 * plane_shape->GetDZ());
    plane_dimension.SetX(2 * plane_shape->GetDY());
  }
  plane_dimension.SetZ(2 * plane_shape->GetDX());

  TVector3 plane_position;
  plane_position.SetX(matrix.GetTranslation()[0]);
  plane_position.SetY(matrix.GetTranslation()[1]);
  plane_position.SetZ(matrix.GetTranslation()[2]);

  plane.setPosition(plane_position);
  plane.setDimension(plane_dimension);
  
  plane.computePlaneVertices();
  plane.computeMaxTransversePosition();

  setSttWireInfo(plane, node, matrix);
}

void SANDGeoManager::setSttWireInfo(sand_geometry::tracker::Plane& plane,
                                       const TGeoNode* const node,
                                       const TGeoHMatrix& matrix)
{
  TVector2 local_plane_x_axis = rotatedToLocal(TVector2(1, 0), plane);
  std::vector<TVector2> vertices = plane.getPlaneVertices();

  for (int i = 0; i < node->GetNdaughters(); i++) {
    sand_geometry::tracker::WireInfo w;

    auto tube_node = node->GetDaughter(i);
    auto tube_matches = stt_tube_regex_.MatchS(tube_node->GetName());
    int tube_id = (reinterpret_cast<TObjString*>(tube_matches->At(4)))
                      ->GetString()
                      .Atoi();
    delete tube_matches;

    sand_geometry::tracker::CellID cell_unique_id = encodeCellId(plane.uId(), sand_geometry::tracker::CellID(tube_id));
    w.setId(sand_geometry::tracker::WireID(cell_unique_id()));
    w.setType(sand_geometry::tracker::WireInfo::Type::kSignal);

    TGeoMatrix* tube_matrix = tube_node->GetMatrix();
    TGeoHMatrix tube_hmatrix = matrix * (*tube_matrix);

    TVector3 tube_position;
    tube_position.SetX(tube_hmatrix.GetTranslation()[0]);
    tube_position.SetY(tube_hmatrix.GetTranslation()[1]);
    tube_position.SetZ(tube_hmatrix.GetTranslation()[2]);
    TVector2 local_2d_position = globalToLocal(TVector2(tube_position.X(), tube_position.Y()), plane);
    
    std::vector<TVector2> intersection_points = 
            getGlobalLinePlaneIntersections(local_2d_position, plane);
    for (const auto& point:intersection_points) {
      TVector3 intersection(point.X(), point.Y(), tube_position.Z());
      w.setPoint(intersection);
    }

    if (w.getPoints().size() == 2) {
      w.setCenter((w.getFirstPoint() + w.getSecondPoint()) * 0.5);
      w.setLength((w.getSecondPoint() - w.getFirstPoint()).Mag());

      if (fabs(w.getFirstPoint().Y() - plane.getPosition().Y() + vertices[0].Y()) < 1E-3) {
        w.setReadoutEnd(sand_geometry::tracker::WireInfo::ReadoutEnd::kFirst);
      }
      if (fabs(w.getFirstPoint().X() - plane.getPosition().X() + vertices[0].X()) < 1E-3) {
        w.setReadoutEnd(sand_geometry::tracker::WireInfo::ReadoutEnd::kFirst);
      }
      if (fabs(w.getSecondPoint().Y() - plane.getPosition().Y() + vertices[0].Y()) < 1E-3) {
        w.setReadoutEnd(sand_geometry::tracker::WireInfo::ReadoutEnd::kSecond);
      }
      if (fabs(w.getSecondPoint().X() - plane.getPosition().X() + vertices[0].X()) < 1E-3) {
        w.setReadoutEnd(sand_geometry::tracker::WireInfo::ReadoutEnd::kSecond);
      }
    }

    TGeoTubeSeg* tube_shape = (TGeoTubeSeg*)tube_node->GetVolume()->GetShape();
    TVector2 rotated_2d_position = localToRotated(local_2d_position, plane);
    plane.addCell(rotated_2d_position.Y(), 
                  sand_geometry::tracker::Cell(cell_unique_id, w, 2. * tube_shape->GetRmax(), 2. * tube_shape->GetRmax(), 
                  tracker_module_configuration::stt::id_to_velocity[std::to_string(plane.lId()())]));
  }
}

const TVector2 SANDGeoManager::pointInRotatedSystem(TVector2 v, double angle) const
{
  TVector2 rotated_v;
  rotated_v.SetX( v.X() * cos(angle) + v.Y() * sin(angle));
  rotated_v.SetY(-v.X() * sin(angle) + v.Y() * cos(angle));

  return rotated_v;
}
const TVector2 SANDGeoManager::globalToLocal(TVector2 global, const sand_geometry::tracker::Plane& plane) const
{
  return TVector2(global.X() - plane.getPosition().X(), global.Y() - plane.getPosition().Y());
}
const TVector2 SANDGeoManager::localToRotated(TVector2 local, const sand_geometry::tracker::Plane& plane) const
{
  return pointInRotatedSystem(local, plane.getRotation());
}
const TVector2 SANDGeoManager::globalToRotated(TVector2 global, const sand_geometry::tracker::Plane& plane) const
{
  return localToRotated(globalToLocal(global, plane), plane);
}
const TVector2 SANDGeoManager::rotatedToLocal(TVector2 rotated, const sand_geometry::tracker::Plane& plane) const
{
  return pointInRotatedSystem(rotated, -plane.getRotation());
}
const TVector2 SANDGeoManager::localToGlobal(TVector2 local, const sand_geometry::tracker::Plane& plane) const
{
  return TVector2(local.X() + plane.getPosition().X(), local.Y() + plane.getPosition().Y());
}
const TVector2 SANDGeoManager::rotatedToGlobal(TVector2 rotated, const sand_geometry::tracker::Plane& plane) const
{
  return localToGlobal(rotatedToLocal(rotated, plane), plane);
}


bool SANDGeoManager::getLineSegmentIntersection(TVector2 p, TVector2 dir, TVector2 A, TVector2 B, TVector2& intersection)
{
  double delta_x = A.X() - B.X();
  double delta_y = A.Y() - B.Y();
  double det = dir.X() * delta_y  - dir.Y() * delta_x;

  if (fabs(det) < 1E-9) {
    // std::cout << "Line and segment are parallel." << std::endl;
    return false;
  } else {
    
    double t = ((A.X() - p.X()) * delta_y - (A.Y() - p.Y()) * delta_x) / det;
    double s = ((p.X() - A.X()) * dir.Y() - (p.Y() - A.Y()) * dir.X()) / det;

    if (s >= 0 && s <= 1) {
      intersection.SetX(p.X() + t * dir.X());
      intersection.SetY(p.Y() + t * dir.Y());
      return true;
    }

    return false;
  }
}

void SANDGeoManager::setDriftPlaneInfo(const TGeoNode* const node,
                                          const TGeoHMatrix& matrix)
{
  // To Do:
  // Check rotation of modules in the geometry.
  // Currently the rotation of the wire is obtained by rotating 
  // the whole plane. This results is the x-y dimensions being swapped
  TString node_path = gGeoManager->GetPath();

  sand_geometry::tracker::PlaneID drift_plane_unique_id = getDriftPlaneId(node_path);
  sand_geometry::tracker::PlaneID drift_plane_local_id  = getDriftPlaneId(node_path, true);  // 0,1 or 2

  planes_.push_back(sand_geometry::tracker::Plane(drift_plane_unique_id, drift_plane_local_id));
  id_to_plane_[planes_.back().uId()] = std::prev(planes_.end());

  auto& plane = planes_.back();
  double angle = tracker_module_configuration::drift::id_to_angle[std::to_string(drift_plane_local_id())];

  plane.setRotation(angle);

  TGeoBBox* plane_shape = (TGeoBBox*)node->GetVolume()->GetShape();
  TVector3 plane_dimension;
  // Notice: this is a workaround to the planes in the geometry being rotated sometimes
  if(drift_plane_local_id() != 2) {
    plane_dimension.SetX(2 * plane_shape->GetDZ());
    plane_dimension.SetY(2 * plane_shape->GetDY());
  } else {
    plane_dimension.SetY(2 * plane_shape->GetDZ());
    plane_dimension.SetX(2 * plane_shape->GetDY());
  }
  plane_dimension.SetZ(2 * plane_shape->GetDX());

  TVector3 plane_position;
  plane_position.SetX(matrix.GetTranslation()[0]);
  plane_position.SetY(matrix.GetTranslation()[1]);
  plane_position.SetZ(matrix.GetTranslation()[2]);

  plane.setPosition(plane_position);
  plane.setDimension(plane_dimension);
  
  plane.computePlaneVertices();
  plane.computeMaxTransversePosition();

  setDriftWireInfo(plane);
}

void SANDGeoManager::setDriftWireInfo(sand_geometry::tracker::Plane& plane)
{

  TVector2 local_plane_x_axis = rotatedToLocal(TVector2(1, 0), plane);
  // local_plane_x_axis.Print();
  
  std::vector<TVector2> vertices = plane.getPlaneVertices();

  double transverse_position = plane.getMaxTransverseCoord() - tracker_module_configuration::drift::id_to_offset[std::to_string(plane.lId()())];
  long wire_id = 0;
  while (transverse_position > -plane.getMaxTransverseCoord()) {
    sand_geometry::tracker::WireInfo w;

    sand_geometry::tracker::CellID cell_unique_id = encodeCellId(plane.uId(), sand_geometry::tracker::CellID(wire_id));
    w.setId(sand_geometry::tracker::WireID(cell_unique_id()));
    w.setType(sand_geometry::tracker::WireInfo::Type::kSignal);
    TVector2 local_2d_position = rotatedToLocal(TVector2(0, transverse_position), plane);
    std::vector<TVector2> intersection_points = 
            getGlobalLinePlaneIntersections(local_2d_position, plane);
    for (const auto& point:intersection_points) {
      TVector3 intersection(point.X(), point.Y(), plane.getPosition().Z());
      w.setPoint(intersection);
    }

    if (w.getPoints().size() == 2) {
      w.setCenter((w.getFirstPoint() + w.getSecondPoint()) * 0.5);
      w.setLength((w.getSecondPoint() - w.getFirstPoint()).Mag());
      if (fabs(w.getFirstPoint().Y() - plane.getPosition().Y() + vertices[0].Y()) < 1E-3) {
        w.setReadoutEnd(sand_geometry::tracker::WireInfo::ReadoutEnd::kFirst);
      }
      if (fabs(w.getFirstPoint().X() - plane.getPosition().X() + vertices[0].X()) < 1E-3) {
        w.setReadoutEnd(sand_geometry::tracker::WireInfo::ReadoutEnd::kFirst);
      }
      if (fabs(w.getSecondPoint().Y() - plane.getPosition().Y() + vertices[0].Y()) < 1E-3) {
        w.setReadoutEnd(sand_geometry::tracker::WireInfo::ReadoutEnd::kSecond);
      }
      if (fabs(w.getSecondPoint().X() - plane.getPosition().X() + vertices[0].X()) < 1E-3) {
        w.setReadoutEnd(sand_geometry::tracker::WireInfo::ReadoutEnd::kSecond);
      }
    }

    if (w.getLength() > tracker_module_configuration::drift::id_to_length[std::to_string(plane.lId()())]) {
      plane.addCell(transverse_position, 
                    sand_geometry::tracker::Cell(cell_unique_id, w, 
                    tracker_module_configuration::drift::id_to_offset[std::to_string(plane.lId()())],
                    plane.getDimension().Z(),
                    tracker_module_configuration::drift::id_to_velocity[std::to_string(plane.lId()())]));
      wire_id++;
    }
    transverse_position -= tracker_module_configuration::drift::id_to_spacing[std::to_string(plane.lId()())];

  }

}


void SANDGeoManager::setPlaneInfo(const TGeoHMatrix& matrix)
{
  TGeoNode* node = gGeoManager->GetCurrentNode();
  TString node_name = node->GetName();
  TGeoMatrix* node_matrix = node->GetMatrix();
  TGeoHMatrix node_hmatrix = matrix * (*node_matrix);
  if (isDriftPlane(node_name)) {
    setDriftPlaneInfo(node, node_hmatrix);
  } else if (isSttPlane(node_name)) {
    setSttPlaneInfo(node, node_hmatrix);
  } else {
    for (int i = 0; i < node->GetNdaughters(); i++) {
      gGeoManager->CdDown(i);
      setPlaneInfo(node_hmatrix);
      gGeoManager->CdUp();
    }
  }
}
double SANDGeoManager::getMinDistanceBetweenSegments(TVector3 a, TVector3 b,
                                                     TVector3 c, TVector3 d)
{
  TVector3 s(b - a);
  TVector3 r(d - c);

  TVector3 diff = a - c;
  double A = s.Dot(s);    // s . s
  double B = s.Dot(r);    // s . r
  double C = r.Dot(r);    // r . r
  double D = s.Dot(diff);    // s . (a - c)
  double E = r.Dot(diff);    // r . (a - c)

  double det = A * C - B * B;
  if ( (fabs(det) > 1E-9)) {
    double t = (B * E - C * D) / det;
    double t_prime = (A * E - B * D) / det;

    t = std::max(0.0, std::min(1.0, t));
    t_prime = std::max(0.0, std::min(1.0, t_prime));
    
    TVector3 point1 = a + t * (b - a);
    TVector3 point2 = c + t_prime * (d - c);
    
    if (t == 0 || t == 1) {
      TVector3 AP = point1 - c;
      t_prime = AP.Dot(r) / r.Mag2();
      t_prime = std::max(0.0, std::min(1.0, t_prime));
    }
    
    if (t_prime == 0 || t_prime == 1) {
      TVector3 AP = point2 - a;
      t = AP.Dot(s) / s.Mag2();
      t = std::max(0.0, std::min(1.0, t));
    }

    point1 = a + t * (b - a);
    point2 = c + t_prime * (d - c);

    return (point1 - point2).Mag();
  } else {
    TVector3 w(a - c);
    double t = w.Dot(r) / r.Mag2();
    TVector3 p_closest = c + t * r;
    return (a - p_closest).Mag();
  }


}

// To Do: I don't like geometry being a string.. maybe use an enum?
// Notice: Currently a single dz and dy are considered. If planes will have 
//        different thickness or different wire smaplings, this won't work
void SANDGeoManager::fillAdjacentCells(std::string geometry)
{
  double dz; 
  double dy;
  auto first_cell  = planes_.at(0).getIdToCellMap().begin();
  auto cell_size = first_cell->second.getSize();
  dy = cell_size.h;
  dz = cell_size.w;
  if (geometry == "STT") {
    dz = dz * sqrt(3) / 2.;
  }
  
  double max_distance = sqrt(dy*dy + dz*dz) + 0.1;
  std::cout << "max_distance " << dy << " " << dz << " " << max_distance << std::endl;

  for(auto plane_it = planes_.begin(); plane_it != planes_.end(); plane_it++) {
    // std::cout << "Checking plane " << plane_it->uId()() << std::endl;
    int c = 0;
    for(auto next_plane_it = plane_it; c < 3 && next_plane_it != planes_.end(); next_plane_it++) {
      // std::cout << "with plane " << next_plane_it->uId()() << std::endl;
      c++;
      auto&      plane_it_cells =      plane_it->getIdToCellMap();
      auto& next_plane_it_cells = next_plane_it->getIdToCellMap();

        
      for(auto& plane_cell:plane_it_cells) {
        for(auto& next_plane_cell:next_plane_it_cells) {
          if(plane_cell.first == next_plane_cell.first) {
            continue;
          }

          double distance = getMinDistanceBetweenSegments(plane_cell.second.getWire().getFirstPoint(),
                                                          plane_cell.second.getWire().getSecondPoint(),
                                                          next_plane_cell.second.getWire().getFirstPoint(),
                                                          next_plane_cell.second.getWire().getSecondPoint());
          // std::cout << c << " " << distance << " " << plane_cell.first() << " " << next_plane_cell.first() << std::endl;
          if (distance < max_distance) {
            plane_cell.second.addAdjacentCell(&(next_plane_cell.second));
            next_plane_cell.second.addAdjacentCell(&(plane_cell.second));
          }
        }
      }
    }
    // break;
  }
}

void SANDGeoManager::rearrangePlanes()
{
  id_to_plane_.clear();
  std::sort(planes_.begin(), planes_.end(), 
            [](const sand_geometry::tracker::Plane& p1, const sand_geometry::tracker::Plane& p2)
              {return p1.getPosition().Z() < p2.getPosition().Z();});

  for (auto it = planes_.begin();
            it != planes_.end(); ++it) {
              id_to_plane_[it->uId()] = it;
            }
            
  for (auto it = planes_.begin();
            it != planes_.end(); ++it) {
              it->updateCells();
            }
}

void SANDGeoManager::setTrackerInfo()
{
  geo_->CdTop();
  TGeoHMatrix matrix = *gGeoIdentity;
  setPlaneInfo(matrix);
  if (geo_->FindVolumeFast("STTtracker_PV")) {
    std::cout << "using SAND tracker : STT\n";
  } else {
    std::cout << "using SAND tracker : DRIFT CHAMBER\n";
  }
  rearrangePlanes();
  // fillAdjacentCells(geometry);
  printModulesInfo(0);
}

void SANDGeoManager::printModulesInfo(int verbose)
{
  std::cout << "There are " << planes_.size() << " planes in the geometry:" << std::endl;
  for (const auto& p:planes_) {
    std::cout << "  - Plane " << p.uId()() << std::endl;
    std::cout << "    Wire rotation: " << p.getRotation() << std::endl;
    std::cout << "    Center position: " << p.getPosition().X() << " " 
                                         << p.getPosition().Y() << " " 
                                         << p.getPosition().Z() << std::endl;
    std::cout << "    Dimensions: "    << p.getDimension().X() << " " 
                                       << p.getDimension().Y() << " " 
                                       << p.getDimension().Z() << std::endl;
    std::cout << "    List of cells (" << p.nCells() << "):"  << std::endl;
    if (verbose >= 1) {
      for(const auto& c:p.getIdToCellMap()) {
        std::cout << "      " << c.first() << std::endl;
        std::cout << "        Center: " << c.second.getWire().getCenter().X() << " "
                                        << c.second.getWire().getCenter().Y() << " "
                                        << c.second.getWire().getCenter().Z() << std::endl;
        std::cout << "        Length: " << c.second.getWire().getLength()     << std::endl;
        std::cout << "        Point1: " << c.second.getWire().getFirstPoint().X() << " "
                                        << c.second.getWire().getFirstPoint().Y() << " "
                                        << c.second.getWire().getFirstPoint().Z() << std::endl;
        std::cout << "        Point2: " << c.second.getWire().getSecondPoint().X() << " "
                                        << c.second.getWire().getSecondPoint().Y() << " "
                                        << c.second.getWire().getSecondPoint().Z() << std::endl;
        std::cout << "        Adjacent ids: ";
        for (const auto& adj:c.second.getAdjacentCell()) std::cout << adj->getId()() << " ";
        std::cout << std::endl;

      }
    }
  }
}

void SANDGeoManager::drawModulesInfo()
{
  gStyle->SetOptStat(0);
  TCanvas cc("", "", 1000, 1000);
  cc.cd();

  auto plane = planes_[3];
  
  TH2D h("","", 10, plane.getPosition().X() - plane.getDimension().X() / 2, plane.getPosition().X() + plane.getDimension().X() / 2, 
                10, plane.getPosition().Y() - plane.getDimension().Y() / 2, plane.getPosition().Y() + plane.getDimension().Y() / 2);

  h.Draw();

  for (const auto& c:plane.getIdToCellMap()) {
    TLine* l = new TLine(c.second.getWire().getFirstPoint().X(), c.second.getWire().getFirstPoint().Y(),
                         c.second.getWire().getSecondPoint().X(), c.second.getWire().getSecondPoint().Y());
    l->Draw("same");
  }
  
  cc.SaveAs("plane.png");
}

//  #######################################################
//  ##                     INIT                          ##
//  #######################################################

void SANDGeoManager::init(TGeoManager* const geo)
{
  std::cout << "> Setting geometry info\n";
  geo_ = geo;
  planes_.clear();
  id_to_plane_.clear();
  set_ecal_info();
  setTrackerInfo();
}

void SANDGeoManager::setGeoCurrentPoint(double x, double y, double z) const
{
  double p[3] = {x, y, z};
  geo_->SetCurrentPoint(p);
}

void SANDGeoManager::setGeoCurrentDirection(double x, double y, double z) const
{
  geo_->SetCurrentDirection(x, y, z);
}

void SANDGeoManager::initVolume(volume& v) const
{
  auto p = geo_->GetCurrentPoint();
  v.geo_volume = geo_->FindNode(p[0], p[1], p[2])->GetVolume();
  v.volume_path = geo_->GetPath();
  if (v.volume_path.Contains("Active")) {
    v.IsActive = true;
  } else {
    v.IsActive = false;
  }
}

int SANDGeoManager::get_ecal_cell_id(double x, double y, double z, bool include_passive) const
{
  if (geo_ == 0) {
    std::cout << "ERROR: TGeoManager pointer not initialized" << std::endl;
    return -999;
  }

  /////
  TGeoNode* node = geo_->FindNode(x, y, z);

  if (node == 0) return -999;

  TString volume_name = node->GetName();
  TString volume_path = geo_->GetPath();

  // std::cout << "vol_path: " << volume_path << "\n";
  // std::cout << "vol name: " << volume_name << "\n";

  if (check_and_process_ecal_path(volume_path) == false) return -999;
  //////
  // temporary solution: the detector_id is currently encoded directly in
  // module_id
  int detector_id = 0;
  int module_id;
  int layer_id;
  int cell_local_id;

  // barrel modules
  if (is_ecal_barrel(volume_name, include_passive)) {
    get_ecal_barrel_module_and_layer(volume_name, volume_path, detector_id,
                                     module_id, layer_id);
    get_ecal_barrel_cell_local_id(x, y, z, node, cell_local_id);
  }

  /* So, if all goes well in the check below, node should be a slab in one
  sections of an endcap module. In get_ecal_endcap_cell_local_id one could use
  geo_ to get the outer module node. Inside that: a function that computes all
  the stuff */

  else if (is_ecal_endcap(volume_name, include_passive)) {
    get_ecal_endcap_module_and_layer(volume_name, volume_path, detector_id,
                                     module_id, layer_id);
    get_ecal_endcap_cell_local_id(x, y, z, module_id, cell_local_id);
  } else {
    // std::cout << ">get_ecal_cell_id exiting with error:\n"
    //           << volume_name << "\n";
    return -999;
  }

  int cell_unique_id =
      encode_ecal_cell_id(detector_id, module_id, layer_id, cell_local_id);

  // std::cout << ">>uniqID: " << cell_unique_id <<", detID: " << detector_id
  //           << ", modID: " << module_id << ", planeID: " << layer_id
  //           << ", cell_local_id: " << cell_local_id << "\n";

  return cell_unique_id;
}

double SANDGeoManager::getHitCellDistance(TVector2 rotated_yz_hit_position, 
                                        std::map<sand_geometry::tracker::CellID, sand_geometry::tracker::Cell>::const_iterator cell_it, 
                                        const sand_geometry::tracker::Plane& plane) const
{
  TVector2 global_wire_xy_position(cell_it->second.getWire().getCenter().X(), cell_it->second.getWire().getCenter().Y());
  TVector2 rotated_wire_2d_position = globalToRotated(global_wire_xy_position, plane);
  
  TVector2 rotated_yz_wire_position(rotated_wire_2d_position.Y(), 
                                    cell_it->second.getWire().getCenter().Z() - plane.getPosition().Z());
  return (rotated_yz_hit_position - rotated_yz_wire_position).Mod();
}

// To Do: check why the STT geometry sometimes gives the wrong result when checking
// for the closest cell. This piece of code was a fix but without understanding the root of the problem
sand_geometry::tracker::CellID SANDGeoManager::getClosestCellToHit(TVector3 hit_center, const sand_geometry::tracker::Plane& plane, bool checkCloseCells = false) const
{
  TVector2 global_hit_xy_position(hit_center.X(), hit_center.Y());
  TVector2 rotated_hit_xy_position = globalToRotated(global_hit_xy_position, plane);
  double transverse_coord = rotated_hit_xy_position.Y();
  
  TVector2 rotated_yz_hit_position(transverse_coord, hit_center.Z() - plane.getPosition().Z());

  std::map<sand_geometry::tracker::CellID, sand_geometry::tracker::Cell>::const_iterator cell_it = plane.getLowerBoundCell(transverse_coord);
  if (cell_it == plane.getIdToCellMapEnd()) {
    cell_it = plane.getIdToCellMap().begin();
  }
  std::map<sand_geometry::tracker::CellID, sand_geometry::tracker::Cell>::const_iterator next_cell_it = std::next(cell_it);
  if (next_cell_it == plane.getIdToCellMapEnd()) {
    next_cell_it = plane.getIdToCellMap().begin();
  }
 
  double distance1 = 1E9;
  double distance2 = 1E9;
  while (true) {
    distance1 = getHitCellDistance(rotated_yz_hit_position, cell_it, plane);
    distance2 = getHitCellDistance(rotated_yz_hit_position, next_cell_it, plane);

    if (checkCloseCells) {
      auto cell_size = cell_it->second.getSize();
      double h = cell_size.h;
      double w = cell_size.w;
      w /= 2;
      h /= 2;
      if (distance1 > w && distance2 > w) {
        if (cell_it != plane.getIdToCellMap().cbegin()) {
          cell_it--;
        } 
        if (std::next(next_cell_it) != plane.getIdToCellMap().end()) {
          next_cell_it++;
        }
        if (cell_it == plane.getIdToCellMap().cbegin() && 
            std::next(next_cell_it) == plane.getIdToCellMap().cend()) {
          break;
        }
      } else {
        break;
      }
    } else {
      break;
    }
  }

  return (distance1 < distance2) ? cell_it->first : next_cell_it->first;
}

sand_geometry::tracker::CellID SANDGeoManager::getSttTubeId(double x, double y, double z) const
{
  if (geo_ == 0) {
    std::cout << "ERROR: TGeoManager pointer not initialized" << std::endl;
    return -999;
  }

  TGeoNode* node = geo_->FindNode(x, y, z);

  TString node_path = gGeoManager->GetPath();
  sand_geometry::tracker::PlaneID stt_plane_unique_id = getSttPlaneId(node_path);

  if (stt_plane_unique_id == -999) {
    return sand_geometry::tracker::CellID(-999);
  }

  auto& plane = planes_.at(getPlaneIndex(stt_plane_unique_id)());

  TVector3 hit_center(x, y, z);
  return getClosestCellToHit(hit_center, plane, true);
}

long SANDGeoManager::printSttTubeId(double x, double y, double z) const
{
  if (geo_ == 0) {
    std::cout << "ERROR: TGeoManager pointer not initialized" << std::endl;
    return -999;
  }

  TGeoNode* node = geo_->FindNode(x, y, z);
  TString volume_name = node->GetName();
  // std::cout << volume_name << "\n";
  return -1;
}

std::vector<sand_geometry::tracker::CellID> SANDGeoManager::getSegmentIds(const TG4HitSegment& hseg)
    const
{

  auto middle = (hseg.Start + hseg.Stop) * 0.5;

  TGeoNode* node = geo_->FindNode(middle.X(), middle.Y(), middle.Z());
  TString node_path = gGeoManager->GetPath();
  sand_geometry::tracker::PlaneID drift_plane_unique_id = getDriftPlaneId(node_path);

  // To Do: use the map?
  auto& plane = planes_.at(getPlaneIndex(drift_plane_unique_id)());

  sand_geometry::tracker::CellID cell_id_start = getClosestCellToHit(hseg.Start.Vect(), plane);
  sand_geometry::tracker::CellID cell_id_stop  = getClosestCellToHit(hseg.Stop.Vect(),  plane);

  return {cell_id_start, cell_id_stop};
}
