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

void SANDTrackerPlane::clear()
{
  for (std::vector<long>::iterator cit = _vFiredCellsIndex.begin();
       cit != _vFiredCellsIndex.end(); cit++) {
    _vCells.at(*cit).isFired(false);
  }

  _vFiredCellsIndex.clear();
}
