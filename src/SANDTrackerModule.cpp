#include "SANDTrackerModule.h"
#include <iostream>

namespace sand_geometry
{

namespace tracker
{
  
std::map<PlaneID, Plane>::iterator Module::getPlane(PlaneID index)
{
  return planes_.find(index);
}

std::map<PlaneID, Plane>::const_iterator Module::getPlane(PlaneID index) const
{  
  return planes_.find(index);
}

bool Module::addPlane(PlaneID plane_unique_id, PlaneID plane_local_id)
{
  auto it = planes_.insert({plane_unique_id, 
                             Plane(plane_unique_id, plane_local_id, this)});
  return it.second;
}
} // namespace tracker
} // namespace sand_geometry