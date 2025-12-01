#pragma once

#include "SANDWireInfo.h"

namespace sand_geometry
{

namespace tracker
{

class Plane;

class CellID : public SingleElStruct<unsigned long>
{
 public:
  using SingleElStruct<unsigned long>::SingleElStruct;
};

class Cell
{
  CellID id_;
  sand_geometry::tracker::WireInfo wire_;
  double width_;
  double height_;
  double driftVelocity_;  // drift velocity in um/ns

  Plane* plane_;
  std::vector<Cell*> adjacent_cells_;

 public:
   struct CellSize {
    double h;
    double w;
  };

  Cell() {};
  Cell(const CellID cell_id,
                  const sand_geometry::tracker::WireInfo &l, 
                  double w,
                  double h,
                  double v,
                  Plane* plane)
      : id_(cell_id),
        wire_(l),
        width_(w),
        height_(h),
        driftVelocity_(v),
        plane_(plane)
  {
  }
  Cell(const CellID cell_id,
                  const sand_geometry::tracker::WireInfo &l, 
                  double w,
                  double h,
                  double v)
      : id_(cell_id),
        wire_(l),
        width_(w),
        height_(h),
        driftVelocity_(v)
  {
  }

  Cell(const CellID cell_id, const sand_geometry::tracker::WireInfo &l, Plane* plane): 
        id_(cell_id), wire_(l), plane_(plane)
  {
  }

  Cell(const Cell &cell);

  void id(CellID wID)
  {
    id_ = wID;
  }
  void setPlane(Plane* p) 
  {
    plane_ = p;
  }
  Plane* getPlane() const 
  {
    return plane_;
  }
  CellID getId() const
  {
    return id_;
  }
  CellSize getSize() const
  {
    return {height_, width_};
  }
  sand_geometry::tracker::WireInfo getWire()
  {
    return wire_;
  }
  const sand_geometry::tracker::WireInfo& getWire() const
  {
    return wire_;
  }
  double getDriftVelocity() const
  {
    return driftVelocity_;
  }

  void addAdjacentCell(Cell* adj_cell);
  const std::vector<Cell*> getAdjacentCell() const {return adjacent_cells_;};

  bool isAdjacent(const CellID& adj_id) const {
    if (std::find_if(adjacent_cells_.begin(), 
                    adjacent_cells_.end(), 
                    [&adj_id](Cell* cell)
                    { return (cell->getId() == adj_id) ? true : false; })
      != adjacent_cells_.end()) {
        return true;
    } else {
      return false;
    }
  }
};

} // namespace tracker
} // namespace sand_geometry