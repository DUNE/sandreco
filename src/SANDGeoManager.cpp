#include "SANDGeoManager.h"

#include <iostream>

#include <TGeoTrd2.h>
#include <TGeoTube.h>
#include <TObjString.h>

bool check_ecal_barrel_geometry_consistency_with_fluka(double xmin, double xmax, double dz) {
    bool condition = (abs(xmin - sand_geometry::ecal::fluka::barrel_module_xmin) > 0.2) ||
        (abs(xmax - sand_geometry::ecal::fluka::barrel_module_xmax) > 0.2) ||
        (abs(dz - sand_geometry::ecal::fluka::barrel_module_thickness) > 0.2);

    if (condition) {
      std::cout << "ERROR ON ECAL GEOMETRY: xmin= " << xmin
                << " instead of what is expected in Fluka"
                << sand_geometry::ecal::fluka::barrel_module_xmin << std::endl;
      std::cout << "ERROR ON ECAL GEOMETRY: xmax= " << xmax
                << " instead of what is expected in Fluka"
                << sand_geometry::ecal::fluka::barrel_module_xmax << std::endl;
      std::cout << "ERROR ON ECAL GEOMETRY: dz= " << dz
                << " instead of what is expected in Fluka"
                << sand_geometry::ecal::fluka::barrel_module_thickness << std::endl;
      // exit(1);
    }

    return condition;
}

bool check_ecal_endcap_geometry_consistency_with_fluka(double ec_r, double ec_dz){

    bool condition = abs(ec_r - sand_geometry::ecal::fluka::endcap_rmax) > 0.2 ||
        (abs(ec_dz - sand_geometry::ecal::fluka::endcap_thickness));

    if (condition) {
      std::cout << "ERROR ON ECAL ENDCAP GEOMETRY: R= " << ec_r
                << " instead of what is expected in Fluka"
                << sand_geometry::ecal::fluka::endcap_rmax << std::endl;
      std::cout << "ERROR ON ECAL ENDCAP GEOMETRY: Thickness= "
                << ec_dz
                << " instead of what is expected in Fluka"
                << sand_geometry::ecal::fluka::endcap_thickness << std::endl;
      //  exit(1);
    }

    return condition;
}

int SANDGeoManager::encode_ecal_barrel_cell_local_id(int layer, int cell)
{
    return cell * 100 + layer;
}

int SANDGeoManager::encode_ecal_endcap_cell_local_id(int layer, int cell)
{
    return cell * 100 + layer;
}

std::pair<int, int> SANDGeoManager::decode_ecal_barrel_cell_local_id(int id)
{
    int cell = id / 100;
    int layer = id % 100;
    return std::make_pair(layer, cell);
}

std::pair<int, int> SANDGeoManager::decode_ecal_endcap_cell_local_id(int id)
{
    int cell = id / 100;
    int layer = id % 100;
    return std::make_pair(layer, cell);
}

std::vector<double> SANDGeoManager::get_levels_z(double half_module_height)
{
  // z edge of the cells
  std::vector<double> zlevel;
  zlevel.push_back(-half_module_height);

  for (int i = 0; i < sand_geometry::ecal::number_of_layers; i++) {
    zlevel.push_back(zlevel.back() + sand_geometry::ecal::layer_thickness[i]);
  }
  return zlevel;
}

std::map<int, TVector3> SANDGeoManager::get_ecal_barrel_cell_center_local_position(const std::vector<double>& zlevels, double m, double q){
  // z position of the center of the cells
  std::map<int, TVector3> ecal_barrel_cell_center_local_positions;
  for (auto i = 0u; i < zlevels.size() - 1u; i++) {
    auto z_this_layer = 0.5 * (zlevels.at(i) + zlevels.at(i + 1));

    // total module width at the z position of the center of the cell
    double x_module_width_at_z = 2 * (m * z_this_layer + q);

    // cell width at the z position of the center of the cell
    double x_cell_width = x_module_width_at_z / sand_geometry::ecal::number_of_cells_per_barrel_layer;

    // x position of the center of the cells
    for (int j = 0; j < sand_geometry::ecal::number_of_cells_per_barrel_layer; j++) {
        auto x = x_cell_width * (j + 0.5) - x_module_width_at_z * 0.5;
        auto y = 0.;
        auto z = z_this_layer;
        auto id = encode_ecal_barrel_cell_local_id(i,j);
        ecal_barrel_cell_center_local_positions[id] = TVector3(x,y,z);
    }
  }
  return ecal_barrel_cell_center_local_positions;
}

std::map<int, TVector3> SANDGeoManager::get_ecal_endcap_cell_center_local_position(const std::vector<double>& zlevels, double rmin, double rmax){
  // z position of the center of the cells
  std::map<int, TVector3> ecal_endcap_cell_center_local_positions;
  for (auto i = 0u; i < zlevels.size() - 1u; i++) {
    auto z_this_layer = 0.5 * (zlevels.at(i) + zlevels.at(i + 1));

    // cell width at the z position of the center of the cell
    double x_cell_width = rmax / sand_geometry::ecal::number_of_cells_per_endcap_layer;

    // x position of the center of the cells
    for (int j = 0; j < sand_geometry::ecal::number_of_cells_per_endcap_layer; j++) {
        auto x = x_cell_width * (j + 0.5) - rmax;
        auto y = 0.;
        auto z = z_this_layer;
        auto id = encode_ecal_endcap_cell_local_id(i,j);
        ecal_endcap_cell_center_local_positions[id] = TVector3(x,y,z);
    }
  }
  return ecal_endcap_cell_center_local_positions;
}

int SANDGeoManager::encode_ecal_cell_id(int detector_id, int module_id, int layer_id, int cell_local_id)
{
  return cell_local_id + 100 * layer_id + 1000 * module_id + detector_id * 100000;
}

void SANDGeoManager::decode_ecal_cell_id(int cell_global_id, int& detector_id, int& module_id, int& layer_id, int& cell_local_id)
{
  detector_id = cell_global_id / 100000;
  cell_global_id -= detector_id * 100000;

  module_id = cell_global_id / 1000;
  cell_global_id -= module_id * 1000;

  layer_id = cell_global_id / 100;
  cell_global_id -= layer_id * 100;

  cell_local_id = cell_global_id;
}

bool SANDGeoManager::is_ecal_barrel(const TString& volume_name)
{
  // something like: volECALActiveSlab_21_PV_0
  return volume_name.Contains("volECAL") == true && volume_name.Contains("Active") == true &&
         volume_name.Contains("end") == false;

}

bool SANDGeoManager::is_ecal_endcap(const TString& volume_name)
{
  // something like: endvolECALActiveSlab_0_PV_0
  return volume_name.Contains("endvolECAL") == true && volume_name.Contains("Active") == true;
}

bool SANDGeoManager::check_and_process_ecal_path(TString& volume_path)
{
  // ENDCAP ==> something like:
  // "/volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_lv_PV_18/volECALActiveSlab_21_PV_0"
  // BARREL ==> something like:
  // "/volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_end_lv_PV_0/endvolECALActiveSlab_0_PV_0"
  TObjArray* obj = volume_path.Tokenize("/");

  int size = obj->GetEntries();
  if (size < 8) {
    return false;
  };

  // BARREL => ECAL_lv_PV_18
  // ENDCAP => ECAL_end_lv_PV_0
  volume_path = ((TObjString*)obj->At(6))->GetString();
  delete obj;

  return true;
}

void get_ecal_barrel_module_and_layer(const TString& volume_name, const TString& volume_path, int& detector_id, int& module_id, int& plane_id)
{
  TObjArray* obja1 = volume_name.Tokenize("_");    // BARERL => volECALActiveSlab_21_PV_0
  TObjArray* obja2 = volume_path.Tokenize("_");  // BARREL => ECAL_lv_PV_18

  // top module => modID == 0
  // increasing modID counterclockwise as seen from positive x
  //(i.e. z(modID==1) < z(modID==0) & z(modID==0) < z(modID==23))
  detector_id = 2;
  module_id = ((TObjString*)obja2->At(3))->GetString().Atoi();
  int slab_id = ((TObjString*)obja1->At(1))->GetString().Atoi();

  delete obja1;
  delete obja2;

  // planeID==0 -> smallest slab -> internal
  // planeID==208 -> biggest slab -> external
  plane_id = slab_id / 40;

  if (plane_id > 4) plane_id = 4;
}

void get_ecal_endcap_module_and_layer(const TString& volume_name, const TString& volume_path, int& detector_id, int& module_id, int& plane_id)
{
  TObjArray* obja1 = volume_name.Tokenize("_");  // ENDCAP => endvolECALActiveSlab_0_PV_0
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

  // planeID==0 -> internal
  // planeID==208 -> external
  plane_id = slab_id / 40;

  if (plane_id > 4) plane_id = 4;
}

void SANDGeoManager::get_ecal_barrel_cell_local_id(double x, double y, double z, const TGeoNode* const node, int& cell_local_id, double& cell_length)
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
  cell_length = trd->GetDy1();

  // http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/html/Detector/Geometry/geomSolids.html
  // if z = -dz -> dx = 2*dx1
  // if z =  dz -> dx = 2*dx2
  // semilarghezza della slab di scintillatore alla quota Plocal[2]
  double dx = 0.5 * local[2] / dz * (dx2 - dx1) + 0.5 * (dx2 + dx1);

  // Cell width at z = Plocal[2]
  double cell_width = 2. * dx / sand_geometry::ecal::number_of_cells_per_barrel_layer;

  // cellID = distanza dall'estremo diviso larghezza cella
  cell_local_id = (local[0] + dx) / cell_width;
}

void SANDGeoManager::get_ecal_endcap_cell_local_id(double x, double y, double z, const TGeoNode* const node, int& cell_local_id, double& cell_length)
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

  cell_length = rmax * TMath::Sin(TMath::ACos(local[0] / rmax));
  cell_local_id = int((local[0] / rmax + 1.) * sand_geometry::ecal::number_of_cells_per_endcap_layer * 0.5);
}

void SANDGeoManager::set_ecal_info()
{
  // https://root.cern.ch/root/htmldoc/guides/users-guide/Geometry.html#shapes
  // GetDx1() half length in x at -Dz
  // GetDx2() half length in x at +Dz
  // Dx1 < Dx2 => -Dz corresponds to minor width => internal side

  double ecal_barrel_xmin;
  double ecal_barrel_xmax;
  double ecal_barrel_dz;
  double ecal_barrel_dy;

  double ecal_endcap_rmin;
  double ecal_endcap_rmax;
  
  if (geo_type_ == SANDGeoType::kFromEdepSim) {
    TGeoTrd2* mod = (TGeoTrd2*)geo_->FindVolumeFast(sand_geometry::ecal::barrel_module_name)->GetShape();
    ecal_barrel_xmin = mod->GetDx1();
    ecal_barrel_xmax = mod->GetDx2();
    ecal_barrel_dz = mod->GetDz();
    ecal_barrel_dy = mod->GetDy1();

    TGeoTube* ec = (TGeoTube*)geo_->FindVolumeFast(sand_geometry::ecal::endcap_module_name)->GetShape();
    ecal_endcap_rmax = ec->GetRmax();  // Maximum radius = 2000
    ecal_endcap_rmin = ec->GetRmin();
  }
  else if (geo_type_ == SANDGeoType::kFromFluka) {
    ecal_barrel_xmin = sand_geometry::ecal::fluka::barrel_module_xmin;;
    ecal_barrel_xmax = sand_geometry::ecal::fluka::barrel_module_xmax;
    ecal_barrel_dz = sand_geometry::ecal::fluka::barrel_module_thickness;
    ecal_barrel_dy = sand_geometry::ecal::fluka::barrel_module_length;

    ecal_endcap_rmax = sand_geometry::ecal::fluka::endcap_rmax;  // Maximum radius = 2000
    ecal_endcap_rmin = 0.;
  }
  else
  {
    std::cout << "ERROR: Unknown Geometry Type" << std::endl;
  }

  // get z of the levels between the layers
  auto z_levels = get_levels_z(ecal_barrel_dz);
  
  // get slop of the edge of the barrel module 
  auto ecal_barrel_edge_slope = 0.5 * (ecal_barrel_xmax - ecal_barrel_xmin) / ecal_barrel_dz;
  auto ecal_barrel_edge_position = 0.5 * (ecal_barrel_xmax + ecal_barrel_xmin);

  // eval barrel cell center global position
  auto ecal_barrel_cell_center_local_positions = get_ecal_barrel_cell_center_local_position(z_levels, ecal_barrel_edge_slope, ecal_barrel_edge_position);
  
  double local[3];
  double master[3];

  for(int module_id = 0; module_id < sand_geometry::ecal::number_of_barrel_modules; module_id++)
  {
    for(auto cell_position: ecal_barrel_cell_center_local_positions)
    {
      local[0] = cell_position.second.X();
      local[1] = cell_position.second.Y();
      local[2] = cell_position.second.Z();

      auto cell_and_layer_id = decode_ecal_barrel_cell_local_id(cell_position.first);
      auto layer_id = cell_and_layer_id.first;
      auto cell_local_id = cell_and_layer_id.second;
      
      geo_->cd(TString::Format(sand_geometry::ecal::path_barrel_template, module_id).Data());
      geo_->LocalToMaster(local, master);

      // here we create new cellInfo
      int detector_id = 2;
      int cell_unique_id = encode_ecal_cell_id(detector_id, module_id, layer_id, cell_local_id);
      cellmap_[cell_unique_id] = SANDECALCellInfo(cell_unique_id, master[0], master[1], master[2], ecal_barrel_dy, SANDECALCellInfo::Orient::kHorizontal);
    }
  }
  
  // eval barrel cell center global position
  auto ecal_endcap_cell_center_local_positions = get_ecal_endcap_cell_center_local_position(z_levels, ecal_endcap_rmin, ecal_endcap_rmax);
  
  for(auto module_id: sand_geometry::ecal::endcap_module_ids)
  {
    for(auto cell_position: ecal_barrel_cell_center_local_positions)
    {
      local[0] = cell_position.second.X();
      local[1] = cell_position.second.Y();
      local[2] = cell_position.second.Z();

      auto cell_and_layer_id = decode_ecal_endcap_cell_local_id(cell_position.first);
      auto layer_id = cell_and_layer_id.first;
      auto cell_local_id = cell_and_layer_id.second;
      int detector_id = 0;
      
      if(module_id == 30)
      {
        detector_id = 3;
        geo_->cd(sand_geometry::ecal::path_endcapR_template);
      }
      else if(module_id == 40)
      {
        detector_id = 1;
        geo_->cd(sand_geometry::ecal::path_endcapL_template);
      }

      geo_->LocalToMaster(local, master);

      // here we create new cellInfo  
      int cell_unique_id = encode_ecal_cell_id(detector_id, module_id, layer_id, cell_local_id);
      cellmap_[cell_unique_id] = SANDECALCellInfo(cell_unique_id, master[0], master[1], master[2],ecal_barrel_dy, SANDECALCellInfo::Orient::kVertical); 

    }
  }
}

int SANDGeoManager::encode_stt_tube_id(int stt_plane_global_id, int stt_tube_local_id)
{
  return stt_tube_local_id * 100000 + stt_plane_global_id;
}

void SANDGeoManager::decode_stt_tube_id(int stt_tube_global_id, int& stt_plane_global_id, int& stt_tube_local_id)
{
  stt_tube_local_id = stt_tube_global_id / 100000;
  stt_plane_global_id = stt_tube_global_id % 100000;  // global id
}

int SANDGeoManager::encode_stt_plane_id(int stt_module_id, int stt_plane_local_id, int stt_plane_type)
{
  return stt_module_id * 100 + stt_plane_local_id * 10 + stt_plane_type;
}

void SANDGeoManager::decode_stt_plane_id(int stt_plane_global_id, int& stt_module_id, int& stt_plane_local_id, int& stt_plane_type)
{
  stt_module_id = stt_plane_global_id / 100;
  stt_plane_local_id = (stt_plane_global_id - stt_module_id * 100) / 10;
  stt_plane_type = stt_plane_global_id % 10;
}

bool SANDGeoManager::is_stt_tube(const TString& volume_name)
{
  return volume_name.Contains(stt_single_tube_regex_);
}

bool SANDGeoManager::is_stt_plane(const TString& volume_name)
{
  return volume_name.Contains(stt_plane_regex_);
}

int SANDGeoManager::get_stt_plane_id(const TString& volume_path)
{
  auto plane_matches = stt_plane_regex_.MatchS(volume_path);
  auto module_matches = stt_module_regex_.MatchS(volume_path);

  if (plane_matches->GetEntries() == 0) {
    delete plane_matches;
    delete module_matches;
    return 0;
  }

  int module_id = (reinterpret_cast<TObjString*>(plane_matches->At(1)))->GetString().Atoi();
  int plane_replica_id = (reinterpret_cast<TObjString*>(plane_matches->At(5)))->GetString().Atoi();
  int plane_type =
      (reinterpret_cast<TObjString*>(plane_matches->At(4)))->GetString().EqualTo("hh")
          ? 2
          : 1;
  int module_replica_id =
      (reinterpret_cast<TObjString*>(module_matches->At(3)))->GetString().Atoi();

  delete plane_matches;
  delete module_matches;

  return encode_stt_plane_id(module_id * 10 + module_replica_id, 2 * plane_replica_id + plane_type, plane_type);
}

void SANDGeoManager::set_stt_tube_info(const TGeoNode* const node, const TGeoHMatrix& matrix, int stt_plane_id)
{
  int stt_plane_type;
  int stt_module_id;
  int stt_plane_local_id;
  decode_stt_plane_id(stt_plane_id, stt_module_id, stt_plane_local_id, stt_plane_type);

  std::map<double, int> this_plane_stt_tube_tranverse_position_map;

  if (stt_plane_type != 1 && stt_plane_type != 2)
    std::cout << "Error: stt plane type expected 0 or 1 -> " << stt_plane_type << std::endl;

  for (int i = 0; i < node->GetNdaughters(); i++) {
    auto two_tubes_node = node->GetDaughter(i);
    auto two_tubes_matches = stt_two_tubes_regex_.MatchS(two_tubes_node->GetName());

    int two_tubes_id =
        (reinterpret_cast<TObjString*>(two_tubes_matches->At(5)))->GetString().Atoi();
    delete two_tubes_matches;

    TGeoMatrix* two_tubes_matrix = two_tubes_node->GetMatrix();
    TGeoHMatrix two_tubes_hmatrix = matrix * (*two_tubes_matrix);

    for (int j = 0; j < two_tubes_node->GetNdaughters(); j++) {
      TGeoNode* tube_node = two_tubes_node->GetDaughter(j);
      TGeoTube* tube_shape = (TGeoTube*)tube_node->GetVolume()->GetShape();
      double tube_length = 2 * tube_shape->GetDz();
      TString tube_volume_name = tube_node->GetName();

      if (!is_stt_tube(tube_volume_name))
        std::cout << "Error: expected ST but not -> " << tube_volume_name.Data()
                  << std::endl;

      TGeoMatrix* tube_matrix = two_tubes_node->GetDaughter(j)->GetMatrix();
      TGeoHMatrix tube_hmatrix = two_tubes_hmatrix * (*tube_matrix);

      auto tube_matches = stt_single_tube_regex_.MatchS(tube_volume_name);

      int tube_id =
          (reinterpret_cast<TObjString*>(tube_matches->At(tube_matches->GetEntries() - 2)))
              ->GetString()
              .Atoi();
      delete tube_matches;

      int tube_local_id = two_tubes_id * 2 + tube_id;
      int tube_unique_id = encode_stt_tube_id(stt_plane_id, tube_local_id);

      int z_coord_index = stt_plane_type % 2;
      int y_coord_index = 1 - z_coord_index;
      TVector3 tube_position;
      tube_position.SetX(tube_hmatrix.GetTranslation()[2]);
      tube_position.SetY(tube_hmatrix.GetTranslation()[y_coord_index]);
      tube_position.SetZ(tube_hmatrix.GetTranslation()[z_coord_index]);

      this_plane_stt_tube_tranverse_position_map[tube_position.Y()] = tube_unique_id;

      // here we fill STT tube info
      sttmap_[tube_unique_id]= SANDSTTTubeInfo(tube_unique_id, 
                                tube_position.X(),
                                tube_position.Y(),
                                tube_position.Z(),
                                tube_length,
                                stt_plane_type == 1 ? SANDSTTTubeInfo::Orient::kVertical : SANDSTTTubeInfo::Orient::kHorizontal,
                                stt_plane_type == 1 ? SANDSTTTubeInfo::ReadoutEnd::kPlus : SANDSTTTubeInfo::ReadoutEnd::kPlus);
    }
  }

  stt_tube_tranverse_position_map_[stt_plane_id] = this_plane_stt_tube_tranverse_position_map;
}

void SANDGeoManager::set_stt_info(const TGeoHMatrix& matrix)
{
  TGeoNode* node = gGeoManager->GetCurrentNode();
  TString node_path = gGeoManager->GetPath();
  TString node_name = node->GetName();
  TGeoMatrix* node_matrix = node->GetMatrix();
  TGeoHMatrix node_hmatrix = matrix * (*node_matrix);

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
    set_stt_info(matrix);
}

void SANDGeoManager::init(TGeoManager* const geo, SANDGeoType geo_type)
{
    geo_ = geo;
    geo_type_ = geo_type;

    set_ecal_info();
    set_stt_info();
}

int SANDGeoManager::get_ecal_cell_id(double x, double y, double z)
{
  /////
  TGeoNode* node = geo_->FindNode(x, y, z);

  if (node == 0) return false;

  TString volume_name = node->GetName();
  TString volume_path = geo_->GetPath();

  if (check_and_process_ecal_path(volume_path) == false)
    return false;
  //////

  int detector_id;
  int module_id;
  int layer_id;
  int cell_local_id;
  double cell_length;

  // barrel modules
  if (is_ecal_barrel(volume_name)) {

    get_ecal_barrel_module_and_layer(volume_name, volume_path, detector_id, module_id, layer_id);
    get_ecal_barrel_cell_local_id(x, y, z, node, cell_local_id, cell_length);

    return true;
  }
  // end cap modules
  else if (is_ecal_endcap(volume_name)) {


    get_ecal_endcap_module_and_layer(volume_name, volume_path, detector_id, module_id, layer_id);
    get_ecal_endcap_cell_local_id(x, y, z, node, cell_local_id, cell_length);

    return true;
  } else {
    return false;
  }

  int cell_unique_id = encode_ecal_cell_id(detector_id, module_id, layer_id, cell_local_id);

  return cell_unique_id;
}

int SANDGeoManager::get_stt_tube_id(double x, double y, double z)
{
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

  std::map<double, int>::iterator it = stt_tube_tranverse_position_map_.at(plane_id).lower_bound(transverse_coord);

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

    if(plane_type == 1)
    {
      v1.SetY(tube1.x());
      v2.SetY(tube2.x());
    }
    else
    {
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