#include "SANDTrackerPlane.h"

#include <iostream>

namespace sand_geometry
{

namespace tracker
{
  
std::map<CellID, Cell>::iterator Plane::getCell(CellID index)
{
  return id_to_cell_map_.find(index);
}
std::map<CellID, Cell>::const_iterator Plane::getCell(CellID index) const
{
  return id_to_cell_map_.find(index);
}

std::map<CellID, Cell>::iterator Plane::getCell(double coord)
{
  return id_to_cell_map_.find(coord_to_id_map_[coord]);
}
std::map<CellID, Cell>::const_iterator Plane::getCell(double coord) const
{
  return id_to_cell_map_.find(coord_to_id_map_.at(coord));
}

std::map<CellID, Cell>::iterator Plane::getLowerBoundCell(double coord)
{
  return id_to_cell_map_.find(coord_to_id_map_.lower_bound(coord)->second);
}
const std::map<CellID, Cell>::const_iterator Plane::getLowerBoundCell(double coord) const
{
  return id_to_cell_map_.find(coord_to_id_map_.lower_bound(coord)->second);
}



void Plane::computePlaneVertices()
{
  vertices_.push_back(TVector2( dimension_.X() / 2,  dimension_.Y() / 2));
  vertices_.push_back(TVector2(-dimension_.X() / 2,  dimension_.Y() / 2));
  vertices_.push_back(TVector2(-dimension_.X() / 2, -dimension_.Y() / 2));
  vertices_.push_back(TVector2( dimension_.X() / 2, -dimension_.Y() / 2));
}

void Plane::computeMaxTransversePosition()
{
  for (auto v:vertices_) {
    double half_x = v.X();
    double half_y = v.Y();

    double max_y = -half_x * sin(rotation_) + half_y * cos(rotation_);
    if (max_y > max_transverse_position_) {
      max_transverse_position_ = max_y;
    }
  }
  // std::cout << max_transverse_position_ << std::endl;
}
} // namespace tracker
} // namespace sand_geometry
