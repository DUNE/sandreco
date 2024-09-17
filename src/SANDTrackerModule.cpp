#include "SANDTrackerModule.h"
#include <iostream>

SANDTrackerPlane &SANDTrackerModule::getPlane(long index)
{
  if (_vPlanes_map.find(index) != _vPlanes_map.end()) {
    return _vPlanes_map[index];
  } else {
    throw "plane not found...i should impelent a class for exception and "
          "handling of that";
  }
  // if (index < (int)_vPlanes.size())
    // return _vPlanes[index];
  // else
    // throw "plane not found...i should impelent a class for exception and "
          // "handling of that";
}

bool SANDTrackerModule::addPlane(const SANDTrackerPlane plane)
{
  long i = plane.id();
  auto it = _vPlanes_map.insert({i, plane});
  return it.second;

  // for (SANDTrackerPlane &p : _vPlanes) {
  //   if (i == p.id()) {
  //     throw "plane already in the list...i should impelent a class for "
  //           "exception and handling of that";
  //   }
  // }
  // _vPlanes.push_back(plane);
}

void SANDTrackerModule::clear()
{
  for (std::vector<SANDTrackerPlane>::iterator pit = _vPlanes.begin();
       pit != _vPlanes.end(); pit++)
    pit->clear();
}
