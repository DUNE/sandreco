#include "SANDTrackerPlane.h"
#include <iostream>

SANDTrackerCell &SANDTrackerPlane::getCell(int index)
{
  if (index < (int)_vCells.size())
    return _vCells[index];
  else
    throw "cell not found...i should impelent a class for exception and "
          "handling of that";
}

void SANDTrackerPlane::propagateTrack(const CLine3D track)
{
  int nc = _vCells.size();
  for (int iCell = 0; iCell < nc; iCell++) {
    SANDTrackerCell &c = _vCells[iCell];

    double h, w;
    c.size(h, w);
    double maxDist = 0.5 * sqrt(h * h + w * w);
    double dist = CLine3D::distance(c.wire(), track);
    if (dist < maxDist) {
      c.isFired(true);
      c.timeResponse(dist / c.driftVelocity());
      _vFiredCellsIndex.push_back(iCell);
    }
  }
}

void SANDTrackerPlane::clear()
{
  for (std::vector<long>::iterator cit = _vFiredCellsIndex.begin();
       cit != _vFiredCellsIndex.end(); cit++) {
    _vCells.at(*cit).isFired(false);
  }

  _vFiredCellsIndex.clear();
}
