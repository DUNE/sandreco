#include "SANDTrackerCell.h"
#include "TVector3.h"

SANDTrackerCell::SANDTrackerCell(const SANDTrackerCell& cell)
    : _wire(cell._wire)
{
  _driftVelocity = cell._driftVelocity;
  _isFired = cell._isFired;
  _timeResponse = cell._timeResponse;
  _id = cell._id;
  _width = cell._width;
  _height = cell._height;
}

void SANDTrackerCell::addAdjacentCell(SANDTrackerCell* adj_cell)
{
  if (std::find(_adjacent_cells.begin(), _adjacent_cells.end(), adj_cell) 
      == _adjacent_cells.end()) {
    _adjacent_cells.push_back(adj_cell);
  }
}