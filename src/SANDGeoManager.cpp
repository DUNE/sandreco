#include "SANDGeoManager.h"

#include <iostream>

#include <TGeoTrd2.h>
#include <TGeoTube.h>

#include "utils.h"

bool check_ecal_barrel_geometry_consistency_with_fluka(double xmin, double xmax, double dz) {
    bool condition = (abs(xmin - sand_reco::ecal::fluka::xmin_f) > 0.2) ||
        (abs(xmax - sand_reco::ecal::fluka::xmax_f) > 0.2) ||
        (abs(dz - sand_reco::ecal::fluka::dz_f) > 0.2);

    if (condition) {
      std::cout << "ERROR ON ECAL GEOMETRY: xmin= " << xmin
                << " instead of what is expected in Fluka"
                << sand_reco::ecal::fluka::xmin_f << std::endl;
      std::cout << "ERROR ON ECAL GEOMETRY: xmax= " << xmax
                << " instead of what is expected in Fluka"
                << sand_reco::ecal::fluka::xmax_f << std::endl;
      std::cout << "ERROR ON ECAL GEOMETRY: dz= " << dz
                << " instead of what is expected in Fluka"
                << sand_reco::ecal::fluka::dz_f << std::endl;
      // exit(1);
    }

    return condition;
}

bool check_ecal_endcap_geometry_consistency_with_fluka(){

    bool condition = abs(sand_reco::ecal::endcap::ec_r - sand_reco::ecal::fluka::ec_rf) > 0.2 ||
        (abs(sand_reco::ecal::endcap::ec_dz - sand_reco::ecal::fluka::ec_dzf));

    if (condition) {
      std::cout << "ERROR ON ECAL ENDCAP GEOMETRY: R= " << sand_reco::ecal::endcap::ec_r
                << " instead of what is expected in Fluka"
                << sand_reco::ecal::fluka::ec_rf << std::endl;
      std::cout << "ERROR ON ECAL ENDCAP GEOMETRY: Thickness= "
                << sand_reco::ecal::endcap::ec_dz
                << " instead of what is expected in Fluka"
                << sand_reco::ecal::fluka::ec_dzf << std::endl;
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

  for (int i = 0; i < sand_geometry::ecal::nLay; i++) {
    zlevel.push_back(zlevel.back() + sand_geometry::ecal::dzlay[i]);
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
    double x_cell_width = x_module_width_at_z / sand_geometry::ecal::nCel;

    // x position of the center of the cells
    for (int j = 0; j < sand_geometry::ecal::nCel; j++) {
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
    double x_cell_width = rmax / sand_geometry::ecal::nCel_ec;

    // x position of the center of the cells
    for (int j = 0; j < sand_geometry::ecal::nCel_ec; j++) {
        auto x = x_cell_width * (j + 0.5) - rmax;
        auto y = 0.;
        auto z = z_this_layer;
        auto id = encode_ecal_endcap_cell_local_id(i,j);
        ecal_endcap_cell_center_local_positions[id] = TVector3(x,y,z);
    }
  }
  return ecal_endcap_cell_center_local_positions;
}

void SANDGeoManager::set_ecal_cell_center_position(TGeoManager* const geo, SANDGeoType geo_type)
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
  
  if (geo_type == SANDGeoType::kFromEdepSim) {
    TGeoTrd2* mod = (TGeoTrd2*)geo->FindVolumeFast(sand_geometry::ecal::barrel_module_name)->GetShape();
    ecal_barrel_xmin = mod->GetDx1();
    ecal_barrel_xmax = mod->GetDx2();
    ecal_barrel_dz = mod->GetDz();
    ecal_barrel_dy = mod->GetDy1();

    TGeoTube* ec = (TGeoTube*)geo->FindVolumeFast(sand_geometry::ecal::endcap_module_name)->GetShape();
    ecal_endcap_rmax = ec->GetRmax();  // Maximum radius = 2000
    ecal_endcap_rmin = ec->GetRmin();
    ecal_endcap_dz = ec->GetDz();   // half of thickness = 115
  }
  else if (geo_type == SANDGeoType::kFromFluka) {
    ecal_barrel_xmin = sand_reco::ecal::fluka::xmin_f;;
    ecal_barrel_xmax = sand_reco::ecal::fluka::xmax_f;
    ecal_barrel_dz = sand_reco::ecal::fluka::dz_f;
    ecal_barrel_dy = sand_reco::ecal::barrel::lCalBarrel;

    ecal_endcap_rmax = sand_reco::ecal::fluka::ec_rf;  // Maximum radius = 2000
    ecal_endcap_rmin = 0.;
    ecal_endcap_dz = sand_reco::ecal::fluka::ec_dzf;   // half of thickness = 115
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

  for(int mod_id = 0; mod_id < sand_geometry::ecal::nMod; mod_id++)
  {
    for(auto cell_position: ecal_barrel_cell_center_local_positions)
    {
      local[0] = cell_position.second.X();
      local[1] = cell_position.second.Y();
      local[2] = cell_position.second.Z();

      auto cell_and_layer_id = decode_ecal_barrel_cell_local_id(cell_position.first);
      auto layer_id = cell_and_layer_id.first;
      auto cell_id = cell_and_layer_id.second;
      
      geo->cd(TString::Format(sand_geometry::ecal::path_barrel_template, mod_id).Data());
      geo->LocalToMaster(local, master);

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
        geo->cd(sand_geometry::ecal::path_endcapR_template);
      else if(endcap_mod_id == 40)
        geo->cd(sand_geometry::ecal::path_endcapL_template);

      geo->LocalToMaster(local, master);

      // here we create new cellInfo      
    }
  }
}

void SANDGeoManager::init(TGeoManager* const geo, SANDGeoType geo_type)
{
    set_ecal_cell_center_position(geo, geo_type);
}