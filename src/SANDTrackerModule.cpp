#include "SANDTrackerModule.h"
#include <iostream>

std::map<long, SANDTrackerPlane>::iterator SANDTrackerModule::getPlane(long index)
{
  return _vPlanes.find(index);
}

std::map<long, SANDTrackerPlane>::const_iterator SANDTrackerModule::getPlane(long index) const
{  
  return _vPlanes.find(index);
}

bool SANDTrackerModule::addPlane(long plane_unique_id, long plane_local_id)
{
  auto it = _vPlanes.insert({plane_unique_id, 
                             SANDTrackerPlane(plane_unique_id, plane_local_id, this)});
  return it.second;
}
