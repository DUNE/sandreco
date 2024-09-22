#include "SANDGeoManager.h"

#include <iostream>

#include <TGeoTrd2.h>
#include <TGeoTube.h>
#include <TObjString.h>

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

std::pair<int, int> SANDGeoManager::decode_ecal_barrel_cell_local_id(
    int id) const
{
  int cell = id / 100;
  int layer = id % 100;
  return std::make_pair(layer, cell);
}

std::pair<int, int> SANDGeoManager::decode_ecal_endcap_cell_local_id(
    int id) const
{
  int cell = id / 100;
  int layer = id % 100;
  return std::make_pair(layer, cell);
}

std::vector<double> SANDGeoManager::get_levels_z(
    double half_module_height) const
{
  // z edge of the cells
  std::vector<double> zlevel;
  zlevel.push_back(-half_module_height);

  for (int i = 0; i < sand_geometry::ecal::number_of_layers; i++) {
    zlevel.push_back(zlevel.back() + sand_geometry::ecal::layer_thickness[i]);
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
//- I guess this should be updated?
std::map<int, TVector3>
    SANDGeoManager::get_ecal_endcap_cell_center_local_position(
        const std::vector<double>& zlevels, double rmin, double rmax) const
{
  // z position of the center of the cells
  std::map<int, TVector3> ecal_endcap_cell_center_local_positions;
  for (auto i = 0u; i < zlevels.size() - 1u; i++) {
    auto z_this_layer = 0.5 * (zlevels.at(i) + zlevels.at(i + 1));

    // cell width at the z position of the center of the cell
    //- I'm not sure how the cell width is defined --> need to look at the
    // geometry
    //- This bit may need to be changed depending on the shape of the new cells
    // x seems the direction along which the endcap is segmented (| | | |o| | |
    // |)
    double x_cell_width =
        2 * rmax / sand_geometry::ecal::number_of_cells_per_endcap_layer;

    //- the x position of the cells now varies over the endcap
    // position of the center of the cells
    for (int j = 0; j < sand_geometry::ecal::number_of_cells_per_endcap_layer;
         j++) {
      auto x = x_cell_width * (j + 0.5) - rmax;
      auto y = 0.;
      auto z = z_this_layer;
      auto id = encode_ecal_endcap_cell_local_id(i, j);
      ecal_endcap_cell_center_local_positions[id] = TVector3(x, y, z);
    }
  }
  return ecal_endcap_cell_center_local_positions;
}

// this should still work
int SANDGeoManager::encode_ecal_cell_id(int detector_id, int module_id,
                                        int layer_id, int cell_local_id)
{
  return cell_local_id + 100 * layer_id + 1000 * module_id +
         detector_id * 100000;
}

void SANDGeoManager::decode_ecal_cell_id(int cell_global_id, int& detector_id,
                                         int& module_id, int& layer_id,
                                         int& cell_local_id)
{
  detector_id = cell_global_id / 100000;
  cell_global_id -= detector_id * 100000;

  module_id = cell_global_id / 1000;
  cell_global_id -= module_id * 1000;

  layer_id = cell_global_id / 100;
  cell_global_id -= layer_id * 100;

  cell_local_id = cell_global_id;
}

bool SANDGeoManager::is_ecal_barrel(const TString& volume_name) const
{
  // something like: volECALActiveSlab_21_PV_0
  return volume_name.Contains("volECAL") == true &&
         volume_name.Contains("Active") == true &&
         volume_name.Contains("end") == false;
}

bool SANDGeoManager::is_ecal_endcap(const TString& volume_name) const
{
  // something like: endvolECALActiveSlab_0_PV_0
  return volume_name.Contains("endvolECAL") == true &&
         volume_name.Contains("Active") == true;
}

bool SANDGeoManager::is_endcap_mod(const TString& volume_name) const
{
  return volume_name.Contains(endcap_mod_regex_);
}

bool SANDGeoManager::check_and_process_ecal_path(TString& volume_path) const
{
  // this bit seems strange: what if other paths end up having more than 8
  // tokens?

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
  volume_path = ((TObjString*)obj->At(6))->GetString() + "/" +
                ((TObjString*)obj->At(7))->GetString() + "/" +
                ((TObjString*)obj->At(8))->GetString();
  delete obj;

  return true;
}
// Note: I think plane_id should be renamed to layer_id for consistency with
// function calls also plane_id is used in the STT functions
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
  int slab_id = ((TObjString*)obja1->At(1))->GetString().Atoi();

  delete obja1;
  delete obja2;

  // layer_id==0 -> smallest slab -> internal (slab_id==0 ?)
  // layer_id==208 -> biggest slab -> external (slab_id==208 ?)
  // why layer_id = slab_id / 40? Because const int number_of_layers = 5 =
  // 208//40
  layer_id = slab_id / 40;

  if (layer_id > 4) layer_id = 4;
}

// Note: I think plane_id should be renamed to layer_id for consistency with
// function calls
/* void SANDGeoManager::get_ecal_endcap_module_and_layer(
    const TString& volume_name, const TString& volume_path, int& detector_id,
    int& module_id, int& layer_id) const
{
  TObjArray* obja1 =
      volume_name.Tokenize("_");  // ENDCAP => endvolECALActiveSlab_0_PV_0
  TObjArray* obja2 = volume_path.Tokenize("_");  // ENDCAP => ECAL_end_lv_PV_0

  module_id = ((TObjString*)obja2->At(4))->GetString().Atoi();
  int slab_id = ((TObjString*)obja1->At(1))->GetString().Atoi();

  // mod == 40 -> left  -> detID = 1
  // mod == 30 -> right -> detID = 3
  // (see issue: https://baltig.infn.it/dune/sand-reco/-/issues/18)
  if (module_id == 0) {
    detector_id = 1;
    module_id = 40;
  } else if (module_id == 1) {
    detector_id = 3;
    module_id = 30;
  }
  delete obja1;
  delete obja2;

  // layer_id==0 -> internal (slab_id==0 ?)
  // layer_id==208 -> external (slab_id==208 ?)
  layer_id = slab_id / 40;

  if (layer_id > 4) layer_id = 4;
} */

// NEW VERSION!!!
void SANDGeoManager::get_ecal_endcap_module_and_layer(
    const TString& volume_name, const TString& volume_path, int& detector_id,
    int& module_id, int& layer_id) const
{
  TObjArray* obja1 =
      volume_name.Tokenize("_");  // ENDCAP => endvolECALActiveSlab_0_PV_0
  TObjArray* obja2 = volume_path.Tokenize(
      "_");  // ENDCAP =>
             // ECAL_end_lv_PV_0/ECAL_ec_mod_4_lv_PV_1/ECAL_ec_mod_vert_0_lv_PV

  int slab_id = ((TObjString*)obja1->At(1))->GetString().Atoi();

  // I'm fairly sure the module_id should be encoded here to get the endcap_mod
  // from the map in get_ecal_endcap_cell_local_id (which uses the module
  // width)-->extract the side and replica indices
  int side_id = ((TObjString*)obja2->At(4))->GetString().Atoi();
  int mod_id = ((TObjString*)obja2->At(8))->GetString().Atoi();
  int replica_id = ((TObjString*)obja2->At(10))->GetString().Atoi();

  module_id = encode_endcap_mod_id(mod_id, replica_id, side_id);

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
}

/* // cell_local structure should be changed
void SANDGeoManager::get_ecal_endcap_cell_local_id(double x, double y, double z,
                                                   const TGeoNode* const node,
                                                   int& cell_local_id) const
{
  double master[3];
  double local[3];
  master[0] = x;
  master[1] = y;
  master[2] = z;

  geo_->GetCurrentNavigator()->MasterToLocal(master, local);

  TGeoTube* tub = (TGeoTube*)node->GetVolume()->GetShape();

  double rmin = tub->GetRmin();
  double rmax = tub->GetRmax();
  double dz = tub->GetDz();

  cell_local_id =
      int((local[0] / rmax + 1.) *
          sand_geometry::ecal::number_of_cells_per_endcap_layer * 0.5);
} */

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
  endcapmap_.at(endcap_mod_id).mod_hmatrix().MasterToLocal(master, local);

  // double rmin = tub->GetRmin();
  // double rmax = tub->GetRmax();
  // double dz = tub->GetDz();

  // Cell width at z = Plocal[2]
  double cell_width = 2. * endcapmap_.at(endcap_mod_id).width() /
                      sand_geometry::ecal::number_of_cells_per_barrel_layer;

  // cellID = distanza dall'estremo diviso larghezza cella
  cell_local_id =
      (local[0] + endcapmap_.at(endcap_mod_id).width()) / cell_width;
}

int SANDGeoManager::get_hit_path_len(const double& hx, const double& hy,
                                     const double& hz, const int& endcap_mod_id,
                                     double& d1, double& d2) const
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
  // check whether the layer is actually contained inside the module
  if (!volume_path.Contains(ec_mod.path())) return -999;

  // convert to layer local coordinates
  geo_->GetCurrentNavigator()->MasterToLocal(master, local);

  // manage each section separately
  if (volume_path.Contains("vert")) {
    auto depth = (ec_mod.mod_dz() - 2 * ec_mod.Al_dz())-local[2];
    d1 =
        0.5 * ec_mod.l_vert() + ec_mod.get_curv_arc_len(depth) + ec_mod.l_hor();
    d2 = ec_mod.get_cell_tot_len(depth) - d1;
  } else if (volume_path.Contains("hor") && volume_path.Contains("lv_PV_0")) {
    auto depth = (ec_mod.mod_dz() - 2 * ec_mod.Al_dz())-local[2];
    d1 = 0.5 * ec_mod.l_hor() - local[1];
    d2 = ec_mod.get_cell_tot_len(depth) - d1;
  } else if (volume_path.Contains("hor") && volume_path.Contains("lv_PV_1")) {
    auto depth = (ec_mod.mod_dz() - 2 * ec_mod.Al_dz()) - local[2];
    d2 = 0.5 * ec_mod.l_hor() - local[1];
    d1 = ec_mod.get_cell_tot_len(depth) - d2;
  } else if (volume_path.Contains("curv") && volume_path.Contains("lv_PV_0")) {
    auto depth = std::sqrt(std::pow(local[1], 2) + std::pow(local[2], 2)) -
                 ec_mod.rmin();
    d1 = std::sqrt(std::pow(local[1], 2) + std::pow(local[2], 2)) *
             std::atan(std::abs(local[2] / local[1])) +
         ec_mod.l_hor();
    d2 = ec_mod.get_cell_tot_len(depth) - d1;
  } else if (volume_path.Contains("curv") && volume_path.Contains("lv_PV_1")) {
    auto depth = std::sqrt(std::pow(local[1], 2) + std::pow(local[2], 2)) -
                 ec_mod.rmin();
    d2 = std::sqrt(std::pow(local[1], 2) + std::pow(local[2], 2)) *
             std::atan(std::abs(local[2] / local[1])) +
         ec_mod.l_hor();
    d1 = ec_mod.get_cell_tot_len(depth) - d2;
  }

  return 1;
}

void SANDGeoManager::set_ecal_info()
{
  // https://root.cern.ch/root/htmldoc/guides/users-guide/Geometry.html#shapes
  // GetDx1() half length in x at -Dz
  // GetDx2() half length in x at +Dz
  // Dx1 < Dx2 => -Dz corresponds to minor width => internal side
  //-- TGeoTrd2 methods, which defines a trapezoid with up and down faces
  // parallel to zy and oblique sides along xz
  TGeoTrd2* mod =
      (TGeoTrd2*)geo_->FindVolumeFast(sand_geometry::ecal::barrel_module_name)
          ->GetShape();
  double ecal_barrel_xmin = mod->GetDx1();
  double ecal_barrel_xmax = mod->GetDx2();
  double ecal_barrel_dz = mod->GetDz();
  double ecal_barrel_dy = mod->GetDy1();
  //-- the geometry probably does not have an endcap_module_name volume anymore
  //-- the TGeoTube must be split into separate parts, corresponding to the
  // single module sections
  //--
  TGeoTube* ec =
      (TGeoTube*)geo_->FindVolumeFast(sand_geometry::ecal::endcap_module_name)
          ->GetShape();
  double ecal_endcap_rmax = ec->GetRmax();  // Maximum radius = 2000
  double ecal_endcap_rmin = ec->GetRmin();

  // get z of the levels between the layers
  /* This bit seems to be still viable: ECAL_lv is the valid barrel module
   * volume name now (again: volume not node as with the regexs)*/
  auto z_levels = get_levels_z(ecal_barrel_dz);

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
          TString::Format(sand_geometry::ecal::path_barrel_template, module_id)
              .Data());

      geo_->LocalToMaster(local, master);

      // here we create new cellInfo
      int detector_id = 2;
      int cell_unique_id =
          encode_ecal_cell_id(detector_id, module_id, layer_id, cell_local_id);
      cellmap_[cell_unique_id] = SANDECALCellInfo(
          cell_unique_id, master[0], master[1], master[2], ecal_barrel_dy,
          SANDECALCellInfo::Orient::kHorizontal);
    }
  }

  // eval endcap cell center global position
  auto ecal_endcap_cell_center_local_positions =
      get_ecal_endcap_cell_center_local_position(z_levels, ecal_endcap_rmin,
                                                 ecal_endcap_rmax);

  for (auto module_id : sand_geometry::ecal::endcap_module_ids) {
    for (auto cell_position : ecal_endcap_cell_center_local_positions) {
      local[0] = cell_position.second.X();
      local[1] = cell_position.second.Y();
      local[2] = cell_position.second.Z();

      double cell_length = ecal_endcap_rmax *
                           TMath::Sin(TMath::ACos(local[0] / ecal_endcap_rmax));

      auto cell_and_layer_id =
          decode_ecal_endcap_cell_local_id(cell_position.first);
      auto layer_id = cell_and_layer_id.first;
      auto cell_local_id = cell_and_layer_id.second;
      int detector_id = 0;

      if (module_id == 30) {
        detector_id = 3;
        geo_->cd(sand_geometry::ecal::path_endcapR_template);
      } else if (module_id == 40) {
        detector_id = 1;
        geo_->cd(sand_geometry::ecal::path_endcapL_template);
      }
      // as if in the geometry the ECAL endcaps were not divided into modules
      // so one has to compute the local position w.r.t. the endcap and
      // transform
      geo_->LocalToMaster(local, master);

      // here we create new cellInfo
      int cell_unique_id =
          encode_ecal_cell_id(detector_id, module_id, layer_id, cell_local_id);
      cellmap_[cell_unique_id] =
          SANDECALCellInfo(cell_unique_id, master[0], master[1], master[2],
                           cell_length, SANDECALCellInfo::Orient::kVertical);
    }
  }
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
    endcapmap_[mod_id] = SANDENDCAPModInfo(mod_id, node, node_hmatrix);
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
  TGeoHMatrix matrix = *gGeoIdentity;
  std::cout << "> Checking endcap info\n";
  set_ecal_endcap_info(matrix);
}

//  #######################################################
//  ##                     STT ENCODING                  ##
//  #######################################################

int SANDGeoManager::encode_stt_tube_id(int stt_plane_global_id,
                                       int stt_tube_local_id)
{
  return stt_tube_local_id * 100000 + stt_plane_global_id;
}

void SANDGeoManager::decode_stt_tube_id(int stt_tube_global_id,
                                        int& stt_plane_global_id,
                                        int& stt_tube_local_id)
{
  stt_tube_local_id = stt_tube_global_id / 100000;
  stt_plane_global_id = stt_tube_global_id % 100000;  // global id
}

int SANDGeoManager::encode_stt_plane_id(int stt_module_id,
                                        int stt_plane_local_id,
                                        int stt_plane_type)
{
  return stt_module_id * 100 + stt_plane_local_id * 10 + stt_plane_type;
}

void SANDGeoManager::decode_stt_plane_id(int stt_plane_global_id,
                                         int& stt_module_id,
                                         int& stt_plane_local_id,
                                         int& stt_plane_type)
{
  stt_module_id = stt_plane_global_id / 100;
  stt_plane_local_id = (stt_plane_global_id - stt_module_id * 100) / 10;
  stt_plane_type = stt_plane_global_id % 10;
}

bool SANDGeoManager::is_stt_tube(const TString& volume_name) const
{
  return volume_name.Contains(stt_tube_regex_);
}

bool SANDGeoManager::is_stt_plane(const TString& volume_name) const
{
  return volume_name.Contains(stt_plane_regex_);
}

int SANDGeoManager::get_stt_plane_id(const TString& volume_path) const
{
  auto plane_matches = stt_plane_regex_.MatchS(volume_path);
  auto module_matches = stt_module_regex_.MatchS(volume_path);

  // if (plane_matches->GetEntries() == 0) {
  //   delete plane_matches;
  //   delete module_matches;
  //   return 0;
  // }

  if (plane_matches->GetEntries() < 5) {
    // Sometimes the volume path returned by the TGeoManager does not match
    // with the expected one for a tube...to be investigated!!!
    // std::cout << "Error: volume path for STT digit not expected!! returning
    // default value (0) for stt plane id" << std::endl;
    delete plane_matches;
    delete module_matches;
    return 0;
  }

  int module_id =
      (reinterpret_cast<TObjString*>(plane_matches->At(2)))->GetString().Atoi();
  int plane_replica_id =
      (reinterpret_cast<TObjString*>(plane_matches->At(4)))->GetString().Atoi();
  int plane_type = (reinterpret_cast<TObjString*>(plane_matches->At(3)))
                           ->GetString()
                           .EqualTo("XX")
                       ? 2
                       : 1;
  int module_replica_id = (reinterpret_cast<TObjString*>(module_matches->At(3)))
                              ->GetString()
                              .Atoi();

  delete plane_matches;
  delete module_matches;

  return encode_stt_plane_id(module_id * 10 + module_replica_id,
                             2 * plane_replica_id + plane_type, plane_type);
}

//- plane_id is obtained from get_stt_plane_id(node_path)
void SANDGeoManager::set_stt_tube_info(const TGeoNode* const node,
                                       const TGeoHMatrix& matrix,
                                       int stt_plane_id)
{
  int stt_plane_type;
  int stt_module_id;
  int stt_plane_local_id;
  decode_stt_plane_id(stt_plane_id, stt_module_id, stt_plane_local_id,
                      stt_plane_type);

  std::map<double, int> this_plane_stt_tube_tranverse_position_map;

  if (stt_plane_type != 1 && stt_plane_type != 2)
    std::cout << "Error: stt plane type expected 0 or 1 -> " << stt_plane_type
              << std::endl;

  for (int i = 0; i < node->GetNdaughters(); i++) {
    auto tube_node = node->GetDaughter(i);
    auto tube_matches = stt_tube_regex_.MatchS(tube_node->GetName());

    int tube_id = (reinterpret_cast<TObjString*>(tube_matches->At(4)))
                      ->GetString()
                      .Atoi();
    delete tube_matches;

    int tube_unique_id = encode_stt_tube_id(stt_plane_id, tube_id);

    TGeoMatrix* tube_matrix = tube_node->GetMatrix();
    TGeoHMatrix tube_hmatrix = matrix * (*tube_matrix);

    TGeoTube* tube_shape = (TGeoTube*)tube_node->GetVolume()->GetShape();
    double tube_length = 2 * tube_shape->GetDz();
    TString tube_volume_name = tube_node->GetName();

    if (!is_stt_tube(tube_volume_name))
      std::cout << "Error: expected ST but not -> " << tube_volume_name.Data()
                << std::endl;

    TVector3 tube_position;
    tube_position.SetX(tube_hmatrix.GetTranslation()[0]);
    tube_position.SetY(tube_hmatrix.GetTranslation()[1]);
    tube_position.SetZ(tube_hmatrix.GetTranslation()[2]);

    double transverse_coord =
        stt_plane_type == 1 ? tube_position.X() : tube_position.Y();

    this_plane_stt_tube_tranverse_position_map[transverse_coord] =
        tube_unique_id;

    // here we fill STT tube info
    sttmap_[tube_unique_id] = SANDSTTTubeInfo(
        tube_unique_id, tube_position.X(), tube_position.Y(), tube_position.Z(),
        tube_length,
        stt_plane_type == 1 ? SANDSTTTubeInfo::Orient::kVertical
                            : SANDSTTTubeInfo::Orient::kHorizontal,
        stt_plane_type == 1 ? SANDSTTTubeInfo::ReadoutEnd::kPlus
                            : SANDSTTTubeInfo::ReadoutEnd::kPlus);
  }

  stt_tube_tranverse_position_map_[stt_plane_id] =
      this_plane_stt_tube_tranverse_position_map;
}

void SANDGeoManager::set_stt_info(const TGeoHMatrix& matrix)
{
  TGeoNode* node = gGeoManager->GetCurrentNode();
  TString node_path = gGeoManager->GetPath();
  TString node_name = node->GetName();
  TGeoMatrix* node_matrix = node->GetMatrix();
  TGeoHMatrix node_hmatrix = matrix * (*node_matrix);

  // std::cout << "node_name: " << (std::string)node_name << "\n";

  //- this checks the node name with regex (names should be consistent with the
  // new gdmls)
  if (is_stt_plane(node_name)) {
    int plane_id = get_stt_plane_id(node_path);
    set_stt_tube_info(node, node_hmatrix, plane_id);
  } else {
    for (int i = 0; i < node->GetNdaughters(); i++) {
      gGeoManager->CdDown(i);
      set_stt_info(node_hmatrix);
      gGeoManager->CdUp();
    }
  }
}

void SANDGeoManager::set_stt_info()
{
  geo_->CdTop();
  TGeoHMatrix matrix = *gGeoIdentity;
  std::cout << "> Setting stt info\n";
  set_stt_info(matrix);
}

//  #######################################################
//  ##                     INIT                          ##
//  #######################################################

void SANDGeoManager::init(TGeoManager* const geo)
{
  std::cout << "> Setting geometry info\n";
  geo_ = geo;

  cellmap_.clear();
  sttmap_.clear();
  endcapmap_.clear();
  stt_tube_tranverse_position_map_.clear();

  set_ecal_endcap_info();
  std::cout << "endcap module map size: " << endcapmap_.size() << "\n";
  // set_stt_info();
  // set_ecal_info();
}

int SANDGeoManager::get_ecal_cell_id(double x, double y, double z) const
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

  if (check_and_process_ecal_path(volume_path) == false) return -999;
  //////

  // temporary solution: the detector_id is currently encoded directly in
  // module_id
  int detector_id = 0;
  int module_id;
  int layer_id;
  int cell_local_id;

  // barrel modules
  if (is_ecal_barrel(volume_name)) {

    get_ecal_barrel_module_and_layer(volume_name, volume_path, detector_id,
                                     module_id, layer_id);
    get_ecal_barrel_cell_local_id(x, y, z, node, cell_local_id);
  }

  /* So, if all goes well in the check below, node should be a slab in one
  sections of an endcap module. In get_ecal_endcap_cell_local_id one could use
  geo_ to get the outer module node. Inside that: a function that computes all
  the stuff */

  // end cap modules --> NEED TO UPDATE THESE TWO FUNCTIONS!! (Do the new endcap
  // IDs conflict with the barrel IDs?)
  else if (is_ecal_endcap(volume_name)) {

    get_ecal_endcap_module_and_layer(volume_name, volume_path, detector_id,
                                     module_id, layer_id);
    get_ecal_endcap_cell_local_id(x, y, z, module_id, cell_local_id);
  } else {
    return -999;
  }

  int cell_unique_id =
      encode_ecal_cell_id(detector_id, module_id, layer_id, cell_local_id);

  return cell_unique_id;
}

int SANDGeoManager::get_stt_tube_id(double x, double y, double z) const
{
  if (geo_ == 0) {
    std::cout << "ERROR: TGeoManager pointer not initialized" << std::endl;
    return -999;
  }

  TGeoNode* node = geo_->FindNode(x, y, z);
  TString volume_name = node->GetName();

  int plane_id = get_stt_plane_id(geo_->GetPath());
  if (plane_id == 0) return -999;

  int tube_id = -999;
  int module_id;
  int plane_local_id;
  int plane_type;
  decode_stt_plane_id(plane_id, module_id, plane_local_id, plane_type);

  double transverse_coord = 0.;

  if (plane_type == 1)
    transverse_coord = x;
  else
    transverse_coord = y;

  std::map<double, int>::const_iterator it =
      stt_tube_tranverse_position_map_.at(plane_id).lower_bound(
          transverse_coord);

  if (it == stt_tube_tranverse_position_map_.at(plane_id).begin()) {
    tube_id = stt_tube_tranverse_position_map_.at(plane_id).begin()->second;
  } else if (it == stt_tube_tranverse_position_map_.at(plane_id).end()) {
    tube_id = stt_tube_tranverse_position_map_.at(plane_id).rbegin()->second;
  } else {
    SANDSTTTubeInfo tube1 = sttmap_.at(it->second);
    SANDSTTTubeInfo tube2 = sttmap_.at(std::prev(it)->second);

    TVector2 v1;
    TVector2 v2;

    v1.SetX(tube1.z());
    v2.SetX(tube2.z());

    if (plane_type == 1) {
      v1.SetY(tube1.x());
      v2.SetY(tube2.x());
    } else {
      v1.SetY(tube1.y());
      v2.SetY(tube2.y());
    }

    TVector2 v(z, transverse_coord);

    if ((v - v1).Mod() > (v - v2).Mod()) {
      if ((v - v2).Mod() > 5)
        std::cout << "Error: distance grater than ST radius" << std::endl;

      tube_id = std::prev(it)->second;
    } else {
      if ((v - v1).Mod() > 5)
        std::cout << "Error: distance grater than ST radius" << std::endl;

      tube_id = it->second;
    }
  }

  return tube_id;
}