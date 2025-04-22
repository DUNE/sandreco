#include "SANDTrackerCell.h"
#include "TVector3.h"

namespace sand_geometry
{

namespace tracker
{
 
Cell::Cell(const Cell& cell)
    : wire_(cell.wire_)
{
  driftVelocity_ = cell.driftVelocity_;
  id_ = cell.id_;
  width_ = cell.width_;
  height_ = cell.height_;
}

void Cell::addAdjacentCell(Cell* adj_cell)
{
  if (std::find(adjacent_cells_.begin(), adjacent_cells_.end(), adj_cell) 
      == adjacent_cells_.end()) {
    adjacent_cells_.push_back(adj_cell);
  }
}
} // namespace tracker
} // namespace sand_geometry
