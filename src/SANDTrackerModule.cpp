#include "SANDTrackerModule.h"
#include <iostream>

std::map<SANDTrackerPlaneID, SANDTrackerPlane>::iterator SANDTrackerModule::getPlane(SANDTrackerPlaneID index)
{
  return _vPlanes.find(index);
}

std::map<SANDTrackerPlaneID, SANDTrackerPlane>::const_iterator SANDTrackerModule::getPlane(SANDTrackerPlaneID index) const
{  
  return _vPlanes.find(index);
}

bool SANDTrackerModule::addPlane(SANDTrackerPlaneID plane_unique_id, SANDTrackerPlaneID plane_local_id)
{
  auto it = _vPlanes.insert({plane_unique_id, 
                             SANDTrackerPlane(plane_unique_id, plane_local_id, this)});
  return it.second;
}
