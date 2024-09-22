#include "SANDTrackerPlane.h"

#include <iostream>

SANDTrackerCell &SANDTrackerPlane::getCell(long index)
{
  if (_id_to_cell_map.find(index) != _id_to_cell_map.end())
    return _id_to_cell_map[index];
  else
    throw "cell not found...i should impelent a class for exception and "
          "handling of that";
}
std::map<long, SANDTrackerCell>::const_iterator SANDTrackerPlane::getCell(long index) const
{
  return _id_to_cell_map.find(index);
}

SANDTrackerCell &SANDTrackerPlane::getCell(double coord)
{
  if (_coord_to_id_map.find(coord) != _coord_to_id_map.end())
    return _id_to_cell_map[_coord_to_id_map[coord]];
  else
    throw "cell not found...i should impelent a class for exception and "
          "handling of that";
}
const SANDTrackerCell &SANDTrackerPlane::getCell(double coord) const
{
  if (_coord_to_id_map.find(coord) != _coord_to_id_map.end())
    return _id_to_cell_map.at(_coord_to_id_map.at(coord));
  else
    throw "cell not found...i should impelent a class for exception and "
          "handling of that";
}

std::map<long, SANDTrackerCell>::iterator SANDTrackerPlane::getLowerBoundCell(double coord)
{
  return _id_to_cell_map.find(_coord_to_id_map.lower_bound(coord)->second);
}
const std::map<long, SANDTrackerCell>::const_iterator SANDTrackerPlane::getLowerBoundCell(double coord) const
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
