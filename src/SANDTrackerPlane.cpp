#include "SANDTrackerPlane.h"

#include <iostream>

std::map<SANDTrackerCellID, SANDTrackerCell>::iterator SANDTrackerPlane::getCell(SANDTrackerCellID index)
{
  return _id_to_cell_map.find(index);
}
std::map<SANDTrackerCellID, SANDTrackerCell>::const_iterator SANDTrackerPlane::getCell(SANDTrackerCellID index) const
{
  return _id_to_cell_map.find(index);
}

std::map<SANDTrackerCellID, SANDTrackerCell>::iterator SANDTrackerPlane::getCell(double coord)
{
  return _id_to_cell_map.find(_coord_to_id_map[coord]);
}
std::map<SANDTrackerCellID, SANDTrackerCell>::const_iterator SANDTrackerPlane::getCell(double coord) const
{
  return _id_to_cell_map.find(_coord_to_id_map.at(coord));
}

std::map<SANDTrackerCellID, SANDTrackerCell>::iterator SANDTrackerPlane::getLowerBoundCell(double coord)
{
  return _id_to_cell_map.find(_coord_to_id_map.lower_bound(coord)->second);
}
const std::map<SANDTrackerCellID, SANDTrackerCell>::const_iterator SANDTrackerPlane::getLowerBoundCell(double coord) const
{
  return _id_to_cell_map.find(_coord_to_id_map.lower_bound(coord)->second);
}



void SANDTrackerPlane::computePlaneVertices()
{
  _vertices.push_back(TVector2( _dimension.X() / 2,  _dimension.Y() / 2));
  _vertices.push_back(TVector2(-_dimension.X() / 2,  _dimension.Y() / 2));
  _vertices.push_back(TVector2(-_dimension.X() / 2, -_dimension.Y() / 2));
  _vertices.push_back(TVector2( _dimension.X() / 2, -_dimension.Y() / 2));
}

void SANDTrackerPlane::computeMaxTransversePosition()
{
  for (auto v:_vertices) {
    double half_x = v.X();
    double half_y = v.Y();

    double max_y = -half_x * sin(_rotation) + half_y * cos(_rotation);
    if (max_y > _max_transverse_position) {
      _max_transverse_position = max_y;
    }
  }
  // std::cout << _max_transverse_position << std::endl;
}
