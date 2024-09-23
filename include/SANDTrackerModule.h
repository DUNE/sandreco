#pragma once
#include "SANDTrackerPlane.h"
#include <string>
class SANDTrackerModule
{
 private:
  std::string _target;
  long _id;
  std::map<long, SANDTrackerPlane> _vPlanes;

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

  bool addPlane(const SANDTrackerPlane plane);
  SANDTrackerPlane& getPlane(long index);
  const SANDTrackerPlane& getPlane(long index) const;

  int nPlanes()
  {
    return _vPlanes.size();
  }
  std::map<long, SANDTrackerPlane>& planes()
  {
    return _vPlanes;
  };
  const std::map<long, SANDTrackerPlane>& planes() const
  {
    return _vPlanes;
  };
};
