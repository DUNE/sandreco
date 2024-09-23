#pragma once
#include "SANDTrackerCell.h"
#include <vector>
#include <map>
class SANDTrackerPlane
{
 private:
  long _unique_id;
  long _local_id;
  double _rotation; //rad, mano destra su z
  TVector3 _position;
  TVector3 _dimension;
  std::map<double, long> _coord_to_id_map;
  std::map<long, SANDTrackerCell> _id_to_cell_map;

  std::vector<TVector2> _vertices;

  double _max_transverse_position = 0;

 public:
  SANDTrackerPlane() {};
  SANDTrackerPlane(long uid, long lid)
  {
    _unique_id = uid;
    _local_id  = lid;
  }
  long uid() const
  {
    return _unique_id;
  }
  long lid() const
  {
    return _local_id;
  }
  void addCell(const double transverse_coordinate, SANDWireInfo w)
  {
    if(_coord_to_id_map.find(transverse_coordinate) == _coord_to_id_map.end()) {
      _coord_to_id_map.insert({transverse_coordinate, w.id()});
      _id_to_cell_map.insert({w.id(), SANDTrackerCell(w)});
    }
  }
  void addCell(const double transverse_coordinate, SANDWireInfo w, double width, double eight)
  {
    if(_coord_to_id_map.find(transverse_coordinate) == _coord_to_id_map.end()) {
      _coord_to_id_map.insert({transverse_coordinate, w.id()});
      _id_to_cell_map.insert({w.id(), SANDTrackerCell(w, width, eight)});
    }
  }
  const std::map<long, SANDTrackerCell>& getIdToCellMap() const {return _id_to_cell_map;};
  int nCells() const
  {
    return _coord_to_id_map.size();
  }

  void computePlaneVertices();
  void computeMaxTransversePosition();
  SANDTrackerCell& getCell(long);
  std::map<long, SANDTrackerCell>::const_iterator  getCell(long) const;
  SANDTrackerCell& getCell(double);
  const SANDTrackerCell& getCell(double) const;
  std::map<long, SANDTrackerCell>::iterator getLowerBoundCell(double);
  const std::map<long, SANDTrackerCell>::const_iterator getLowerBoundCell(double) const;
  TVector3 getPosition()  const {return _position;} ;
  TVector3 getDimension() const {return _dimension;} ;
  double getRotation() const {return _rotation;} ;
  double getMaxTransverseCoord() const {return _max_transverse_position;} ;
  std::vector<TVector2> getPlaneVertices() const {return _vertices;} ;
  void setPosition(TVector3 p)  { _position  = p;};
  void setDimension(TVector3 d) { _dimension = d;};
  void setRotation(double r) {_rotation = r;};
};
