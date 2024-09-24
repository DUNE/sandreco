#pragma once
#include "SANDTrackerCell.h"
#include <vector>
#include <map>
class SANDTrackerModule;

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

  SANDTrackerModule* _module = nullptr;

 public:
  SANDTrackerPlane() {};
  SANDTrackerPlane(long u_id, long l_id, SANDTrackerModule* module_ptr)
  {
    _unique_id = u_id;
    _local_id  = l_id;
    _module = module_ptr;
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
      _id_to_cell_map.insert({w.id(), SANDTrackerCell(w, this)});
    }
  }
  void addCell(const double transverse_coordinate, SANDWireInfo w, double width, double eight, double velocity)
  {
    if(_coord_to_id_map.find(transverse_coordinate) == _coord_to_id_map.end()) {
      _coord_to_id_map.insert({transverse_coordinate, w.id()});
      _id_to_cell_map.insert({w.id(), SANDTrackerCell(w, width, eight, velocity, this)});
    }
  }
  const std::map<long, SANDTrackerCell>& getIdToCellMap() const {return _id_to_cell_map;};
  std::map<long, SANDTrackerCell>::const_iterator getIdToCellMapEnd() const {return _id_to_cell_map.end();};
  int nCells() const
  {
    return _coord_to_id_map.size();
  }

  void computePlaneVertices();
  void computeMaxTransversePosition();
  std::map<long, SANDTrackerCell>::iterator getCell(long);
  std::map<long, SANDTrackerCell>::const_iterator  getCell(long) const;
  std::map<long, SANDTrackerCell>::iterator getCell(double);
  std::map<long, SANDTrackerCell>::const_iterator getCell(double) const;
  std::map<long, SANDTrackerCell>::iterator getLowerBoundCell(double);
  const std::map<long, SANDTrackerCell>::const_iterator getLowerBoundCell(double) const;
  TVector3 getPosition()  const {return _position;} ;
  TVector3 getDimension() const {return _dimension;} ;
  double getRotation() const {return _rotation;} ;
  double getMaxTransverseCoord() const {return _max_transverse_position;} ;
  std::vector<TVector2> getPlaneVertices() const {return _vertices;} ;
  SANDTrackerModule* getModule() const {return _module;} ;
  void setPosition(TVector3 p)  { _position  = p;};
  void setDimension(TVector3 d) { _dimension = d;};
  void setRotation(double r) {_rotation = r;};
};
