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
  double ecal_endcap_dz;
  
  if (geo_type_ == SANDGeoType::kFromEdepSim) {
    TGeoTrd2* mod = (TGeoTrd2*)geo_->FindVolumeFast(sand_geometry::ecal::barrel_module_name)->GetShape();
    ecal_barrel_xmin = mod->GetDx1();
    ecal_barrel_xmax = mod->GetDx2();
    ecal_barrel_dz = mod->GetDz();
    ecal_barrel_dy = mod->GetDy1();

    TGeoTube* ec = (TGeoTube*)geo_->FindVolumeFast(sand_geometry::ecal::endcap_module_name)->GetShape();
    ecal_endcap_rmax = ec->GetRmax();  // Maximum radius = 2000
    ecal_endcap_rmin = ec->GetRmin();
    ecal_endcap_dz = ec->GetDz();   // half of thickness = 115
  }
  else if (geo_type_ == SANDGeoType::kFromFluka) {
    ecal_barrel_xmin = sand_geometry::ecal::fluka::barrel_module_xmin;;
    ecal_barrel_xmax = sand_geometry::ecal::fluka::barrel_module_xmax;
    ecal_barrel_dz = sand_geometry::ecal::fluka::barrel_module_thickness;
    ecal_barrel_dy = sand_geometry::ecal::fluka::barrel_module_length;

    ecal_endcap_rmax = sand_geometry::ecal::fluka::endcap_rmax;  // Maximum radius = 2000
    ecal_endcap_rmin = 0.;
    ecal_endcap_dz = sand_geometry::ecal::fluka::endcap_thickness;   // half of thickness = 115
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

  for(int mod_id = 0; mod_id < sand_geometry::ecal::number_of_barrel_modules; mod_id++)
  {
    for(auto cell_position: ecal_barrel_cell_center_local_positions)
    {
      local[0] = cell_position.second.X();
      local[1] = cell_position.second.Y();
      local[2] = cell_position.second.Z();

      auto cell_and_layer_id = decode_ecal_barrel_cell_local_id(cell_position.first);
      auto layer_id = cell_and_layer_id.first;
      auto cell_id = cell_and_layer_id.second;
      
      geo_->cd(TString::Format(sand_geometry::ecal::path_barrel_template, mod_id).Data());
      geo_->LocalToMaster(local, master);

      // here we create new cellInfo         
    }
  }
  
  // eval barrel cell center global position
  auto ecal_endcap_cell_center_local_positions = get_ecal_endcap_cell_center_local_position(z_levels, ecal_endcap_rmin, ecal_endcap_rmax);
  
  for(auto endcap_mod_id: sand_geometry::ecal::endcap_module_ids)
  {
    for(auto cell_position: ecal_barrel_cell_center_local_positions)
    {
      local[0] = cell_position.second.X();
      local[1] = cell_position.second.Y();
      local[2] = cell_position.second.Z();

      auto cell_and_layer_id = decode_ecal_endcap_cell_local_id(cell_position.first);
      auto layer_id = cell_and_layer_id.first;
      auto cell_id = cell_and_layer_id.second;

      if(endcap_mod_id == 30)
        geo_->cd(sand_geometry::ecal::path_endcapR_template);
      else if(endcap_mod_id == 40)
        geo_->cd(sand_geometry::ecal::path_endcapL_template);

      geo_->LocalToMaster(local, master);

      // here we create new cellInfo    
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

  if (stt_plane_type != 0 && stt_plane_type != 1)
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

    if(plane_type == 1)
    {
    }
    else
    {
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

  return encode_stt_tube_id(plane_id, tube_id);
}