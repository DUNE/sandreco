#include "SANDTrackerModule.h"
#include <iostream>

SANDTrackerPlane &SANDTrackerModule::getPlane(long index)
{
  if (_vPlanes.find(index) != _vPlanes.end()) {
    return _vPlanes[index];
  } else {
    throw "plane not found...i should impelent a class for exception and "
          "handling of that";
  }
}

bool SANDTrackerModule::addPlane(const SANDTrackerPlane plane)
{
  long ui = plane.uid();
  auto it = _vPlanes.insert({ui, plane});
  return it.second;
}
