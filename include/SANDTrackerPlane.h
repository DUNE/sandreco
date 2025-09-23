#pragma once
#include "SANDTrackerCell.h"
#include <vector>
#include <map>
namespace sand_geometry
{

namespace tracker
{
  
class Module;

class PlaneID : public SingleElStruct<long>
{
 public:
  PlaneID(long id) : SingleElStruct<long>(id){};
  PlaneID() : SingleElStruct<long>(){};
};

class PlaneIndex : public SingleElStruct<long>
{
 public:
  PlaneIndex(long id) : SingleElStruct<long>(id){};
  PlaneIndex() : SingleElStruct<long>(){};
};

class PlaneLocalID : public SingleElStruct<long>
{
 public:
  PlaneLocalID(long id) : SingleElStruct<long>(id){};
  PlaneLocalID() : SingleElStruct<long>(){};
};

// To Do: Use the local id

class Plane
{
 private:
  PlaneID unique_id_;
  PlaneID local_id_;
  double rotation_; //rad, mano destra su z
  TVector3 position_;
  TVector3 dimension_;
  std::map<double, CellID> coord_to_id_map_;
  std::map<CellID, Cell> id_to_cell_map_;

  std::vector<TVector2> vertices_;

  double max_transverse_position_ = 0;

  Module* module_ = nullptr;

 public:
  Plane() {};
  Plane(PlaneID u_id, PlaneID l_id)
  {
    unique_id_ = u_id;
    local_id_  = l_id;
  }
  Plane(PlaneID u_id, PlaneID l_id, Module* module_ptr)
  {
    unique_id_ = u_id;
    local_id_  = l_id;
    module_ = module_ptr;
  }
  PlaneID uId() const
  {
    return unique_id_;
  }
  PlaneID lId() const
  {
    return local_id_;
  }
  void addCell(const double transverse_coordinate, Cell c)
  {
    if(coord_to_id_map_.find(transverse_coordinate) == coord_to_id_map_.end()) {
      coord_to_id_map_.insert({transverse_coordinate, c.getId()});
      id_to_cell_map_.insert({c.getId(), c});
      id_to_cell_map_[c.getId()].setPlane(this);
    }
  }
  void updateCells() {
    for (auto& c:id_to_cell_map_) {
      c.second.setPlane(this);
    }
  }
        std::map<CellID, Cell>& getIdToCellMap()       {return id_to_cell_map_;};
  const std::map<CellID, Cell>& getIdToCellMap() const {return id_to_cell_map_;};
        std::map<double, CellID>& getCoordToIDMap()       {return coord_to_id_map_;};
  const std::map<double, CellID>& getCoordToIDMap() const {return coord_to_id_map_;};
  std::map<CellID, Cell>::const_iterator getIdToCellMapEnd() const {return id_to_cell_map_.end();};
  int nCells() const
  {
    return coord_to_id_map_.size();
  }

  void computePlaneVertices();
  void computeMaxTransversePosition();
  std::map<CellID, Cell>::iterator getCell(CellID);
  std::map<CellID, Cell>::const_iterator  getCell(CellID) const;
  std::map<CellID, Cell>::iterator getCell(double);
  std::map<CellID, Cell>::const_iterator getCell(double) const;
  std::map<CellID, Cell>::iterator getLowerBoundCell(double);
  const std::map<CellID, Cell>::const_iterator getLowerBoundCell(double) const;
  TVector3 getPosition()  const {return position_;} ;
  TVector3 getDimension() const {return dimension_;} ;
  double getRotation() const {return rotation_;} ;
  double getMaxTransverseCoord() const {return max_transverse_position_;} ;
  std::vector<TVector2> getPlaneVertices() const {return vertices_;} ;
  Module* getModule() const {return module_;} ;
  void setPosition(TVector3 p)  { position_  = p;};
  void setDimension(TVector3 d) { dimension_ = d;};
  void setRotation(double r) {rotation_ = r;};
};

using plane_iterator = std::vector<Plane>::const_iterator;

} // namespace tracker
} // namespace sand_geometry
