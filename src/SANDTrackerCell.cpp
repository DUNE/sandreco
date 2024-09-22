#include "SANDTrackerCell.h"
#include "TVector3.h"

SANDTrackerCell::SANDTrackerCell(const SANDTrackerCell& cell)
    : _wire(cell._wire)
{
  _driftVelocity = cell._driftVelocity;
  _isFired = cell._isFired;
  _timeResponse = cell._timeResponse;
  _wireID = cell._wireID;
  _width = cell._width;
  _height = cell._height;
}
