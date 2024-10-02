#pragma once
#include "SANDTrackerCell.h"
#include <vector>
#include <map>
class SANDTrackerModule;

class SANDTrackerPlaneID : public SingleElStruct<unsigned long>
{
 public:
  SANDTrackerPlaneID(unsigned long id) : SingleElStruct<unsigned long>(id){};
  SANDTrackerPlaneID() : SingleElStruct<unsigned long>(){};
};

class SANDTrackerPlaneIndex : public SingleElStruct<unsigned long>
{
 public:
  SANDTrackerPlaneIndex(unsigned long id) : SingleElStruct<unsigned long>(id){};
  SANDTrackerPlaneIndex() : SingleElStruct<unsigned long>(){};
};

class SANDTrackerPlaneLocalID : public SingleElStruct<unsigned long>
{
 public:
  SANDTrackerPlaneLocalID(unsigned long id) : SingleElStruct<unsigned long>(id){};
  SANDTrackerPlaneLocalID() : SingleElStruct<unsigned long>(){};
};

// To Do: Use the local id

class SANDTrackerPlane
{
 private:
  SANDTrackerPlaneID _unique_id;
  SANDTrackerPlaneID _local_id;
  double _rotation; //rad, mano destra su z
  TVector3 _position;
  TVector3 _dimension;
  std::map<double, SANDTrackerCellID> _coord_to_id_map;
  std::map<SANDTrackerCellID, SANDTrackerCell> _id_to_cell_map;

  std::vector<TVector2> _vertices;

  double _max_transverse_position = 0;

  SANDTrackerModule* _module = nullptr;

 public:
  SANDTrackerPlane() {};
  SANDTrackerPlane(SANDTrackerPlaneID u_id, SANDTrackerPlaneID l_id)
  {
    _unique_id = u_id;
    _local_id  = l_id;
  }
  SANDTrackerPlane(SANDTrackerPlaneID u_id, SANDTrackerPlaneID l_id, SANDTrackerModule* module_ptr)
  {
    _unique_id = u_id;
    _local_id  = l_id;
    _module = module_ptr;
  }
  SANDTrackerPlaneID uid() const
  {
    return _unique_id;
  }
  SANDTrackerPlaneID lid() const
  {
    return _local_id;
  }
  void addCell(const double transverse_coordinate, SANDTrackerCell c)
  {
    if(_coord_to_id_map.find(transverse_coordinate) == _coord_to_id_map.end()) {
      c.setPlane(this);
      _coord_to_id_map.insert({transverse_coordinate, c.id()});
      _id_to_cell_map.insert({c.id(), c});
    }
  }
        std::map<SANDTrackerCellID, SANDTrackerCell>& getIdToCellMap()       {return _id_to_cell_map;};
  const std::map<SANDTrackerCellID, SANDTrackerCell>& getIdToCellMap() const {return _id_to_cell_map;};
        std::map<double, SANDTrackerCellID>& getCoordToIDMap()       {return _coord_to_id_map;};
  const std::map<double, SANDTrackerCellID>& getCoordToIDMap() const {return _coord_to_id_map;};
  std::map<SANDTrackerCellID, SANDTrackerCell>::const_iterator getIdToCellMapEnd() const {return _id_to_cell_map.end();};
  int nCells() const
  {
    return _coord_to_id_map.size();
  }

  void computePlaneVertices();
  void computeMaxTransversePosition();
  std::map<SANDTrackerCellID, SANDTrackerCell>::iterator getCell(SANDTrackerCellID);
  std::map<SANDTrackerCellID, SANDTrackerCell>::const_iterator  getCell(SANDTrackerCellID) const;
  std::map<SANDTrackerCellID, SANDTrackerCell>::iterator getCell(double);
  std::map<SANDTrackerCellID, SANDTrackerCell>::const_iterator getCell(double) const;
  std::map<SANDTrackerCellID, SANDTrackerCell>::iterator getLowerBoundCell(double);
  const std::map<SANDTrackerCellID, SANDTrackerCell>::const_iterator getLowerBoundCell(double) const;
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

using plane_iterator = std::vector<SANDTrackerPlane>::const_iterator;
