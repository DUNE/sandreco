#include "SANDTrackerModule.h"
#include <iostream>

SANDTrackerPlane &SANDTrackerModule::getPlane(int index)
{
  if (index < (int)_vPlanes.size())
    return _vPlanes[index];
  else
    throw "plane not found...i should impelent a class for exception and "
          "handling of that";
}

void SANDTrackerModule::addPlane(const SANDTrackerPlane plane)
{
  long i = plane.id();
  for (SANDTrackerPlane &p : _vPlanes) {
    if (i == p.id()) {
      throw "plane already in the list...i should impelent a class for "
            "exception and handling of that";
    }
  }
  _vPlanes.push_back(plane);
}

void SANDTrackerModule::propagateTrack(const CLine3D track)
{

  for (std::vector<SANDTrackerPlane>::iterator pit = _vPlanes.begin();
       pit != _vPlanes.end(); pit++) {
    pit->propagateTrack(track);
  }
}

void SANDTrackerModule::clear()
{
  for (std::vector<SANDTrackerPlane>::iterator pit = _vPlanes.begin();
       pit != _vPlanes.end(); pit++)
    pit->clear();
}
