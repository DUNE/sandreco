#pragma once
#include "SANDTrackerPlane.h"
#include <string>

namespace sand_geometry
{

namespace tracker
{
 
class ModuleID : public SingleElStruct<unsigned long>
{
 public:
  using SingleElStruct<unsigned long>::SingleElStruct;
};

class Module
{
 private:
  std::string target_;
  ModuleID id_;
  std::map<PlaneID, Plane> planes_;

 public:
  Module()
  {
  }
  Module(ModuleID i)
  {
    id_ = i;
  }
  Module(ModuleID i, std::string trg)
  {
    id_ = i;
    target_ = trg;
  }
  
  void Id(const ModuleID idModule)
  {
    id_ = idModule;
  }
  ModuleID Id() const
  {
    return id_;
  }

  void target(const std::string tagetName)
  {
    target_ = tagetName;
  }
  std::string getTarget() const
  {
    return target_;
  }

  bool addPlane(PlaneID plane_unique_id, PlaneID plane_local_id);
  std::map<PlaneID, Plane>::iterator getPlane(PlaneID index);
  std::map<PlaneID, Plane>::const_iterator getPlane(PlaneID index) const;

  int nPlanes() const
  {
    return planes_.size();
  }
  std::map<PlaneID, Plane>& planes()
  {
    return planes_;
  };
  const std::map<PlaneID, Plane>& planes() const
  {
    return planes_;
  };
};
} // namespace tracker
} // namespace sand_geometry
