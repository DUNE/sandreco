#pragma once
#include "SANDTrackerPlane.h"
#include <string>
class SANDTrackerModule
{
 private:
  std::string _target;
  long _id;
  std::vector<SANDTrackerPlane> _vPlanes;
  std::map<long, SANDTrackerPlane> _vPlanes_map;

 public:
  SANDTrackerModule()
  {
  }
  SANDTrackerModule(long i)
  {
    _id = i;
  }
  SANDTrackerModule(long i, std::string trg)
  {
    _id = i;
    _target = trg;
  }
  // maybe not needed
  SANDTrackerModule(long i, std::string trg,
                    std::vector<SANDTrackerPlane> vPlanes)
  {
    _id = i;
    _target = trg;
    _vPlanes = vPlanes;
  }

  void id(const long idModule)
  {
    _id = idModule;
  }
  long id() const
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

  void addPlane(const SANDTrackerPlane plane);
  SANDTrackerPlane& getPlane(int index);

  int nPlanes() const
  {
    return _vPlanes.size();
  }
  std::vector<SANDTrackerPlane>& planes()
  {
    return _vPlanes;
  };

  void propagateTrack(const CLine3D track);
  void clear();
};
