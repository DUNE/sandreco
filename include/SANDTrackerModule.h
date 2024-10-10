#pragma once
#include "SANDTrackerPlane.h"
#include <string>

class SANDTrackerModuleID : public SingleElStruct<unsigned long>
{
 public:
  SANDTrackerModuleID(unsigned long id) : SingleElStruct<unsigned long>(id){};
  SANDTrackerModuleID() : SingleElStruct<unsigned long>(){};
};

class SANDTrackerModule
{
 private:
  std::string _target;
  SANDTrackerModuleID _id;
  std::map<SANDTrackerPlaneID, SANDTrackerPlane> _vPlanes;

 public:
  SANDTrackerModule()
  {
  }
  SANDTrackerModule(SANDTrackerModuleID i)
  {
    _id = i;
  }
  SANDTrackerModule(SANDTrackerModuleID i, std::string trg)
  {
    _id = i;
    _target = trg;
  }
  
  void id(const SANDTrackerModuleID idModule)
  {
    _id = idModule;
  }
  SANDTrackerModuleID id() const
  {
    return _id;
  }

  void target(const std::string tagetName)
  {
    _target = tagetName;
  }
  std::string target() const
  {
    return _target;
  }

  bool addPlane(SANDTrackerPlaneID plane_unique_id, SANDTrackerPlaneID plane_local_id);
  std::map<SANDTrackerPlaneID, SANDTrackerPlane>::iterator getPlane(SANDTrackerPlaneID index);
  std::map<SANDTrackerPlaneID, SANDTrackerPlane>::const_iterator getPlane(SANDTrackerPlaneID index) const;

  int nPlanes() const
  {
    return _vPlanes.size();
  }
  std::map<SANDTrackerPlaneID, SANDTrackerPlane>& planes()
  {
    return _vPlanes;
  };
  const std::map<SANDTrackerPlaneID, SANDTrackerPlane>& planes() const
  {
    return _vPlanes;
  };
};
