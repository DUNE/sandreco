#pragma once

#include "SANDWireInfo.h"

class SANDTrackerPlane;

class SANDTrackerCellID : public SingleElStruct<unsigned long>
{
 public:
  SANDTrackerCellID(unsigned long id) : SingleElStruct<unsigned long>(id){};
  SANDTrackerCellID() : SingleElStruct<unsigned long>(){};
};

class SANDTrackerCell
{
  SANDTrackerCellID _id;
  SANDWireInfo _wire;
  double _width;
  double _height;

  double _timeResponse;   // time response in nanoseconds
  double _driftVelocity;  // drift velocity in um/ns
  bool _isFired;

  SANDTrackerPlane* _plane;

 public:
  SANDTrackerCell() {};
  SANDTrackerCell(const SANDTrackerCellID cID, const SANDWireInfo &l, const double w, const double h, 
                  const double time, const double vd, const bool fired, SANDTrackerPlane* plane)
      : _id(cID),
        _wire(l),
        _width(w),
        _height(h),
        _timeResponse(time),
        _driftVelocity(vd),
        _isFired(fired),
        _plane(plane)

  {
  }
  
  SANDTrackerCell(const SANDTrackerCellID cID,
                  const SANDWireInfo &l, 
                  const double w,
                  const double h,
                  const double v,
                  SANDTrackerPlane* plane)
      : _id(cID),
        _wire(l),
        _width(w),
        _height(h),
        _driftVelocity(v),
        _plane(plane)

  {
  }
  SANDTrackerCell(const SANDTrackerCellID cID,
                  const SANDWireInfo &l, 
                  const double w,
                  const double h,
                  const double v)
      : _id(cID),
        _wire(l),
        _width(w),
        _height(h),
        _driftVelocity(v)

  {
  }

  SANDTrackerCell(const SANDTrackerCellID cID, const SANDWireInfo &l, SANDTrackerPlane* plane): 
        _id(cID), _wire(l), _plane(plane)
  {
  }

  SANDTrackerCell(const SANDTrackerCell &cell);

  void timeResponse(double time)
  {
    _timeResponse = time;
  }
  void isFired(bool fired)
  {
    _isFired = fired;
  }
  void id(SANDTrackerCellID wID)
  {
    _id = wID;
  }

  double timeResponse() const
  {
    return _timeResponse;
  }
  bool isFired() const
  {
    return _isFired;
  }
  void setPlane(SANDTrackerPlane* p) 
  {
    _plane = p;
  }
  SANDTrackerPlane* getPlane() const 
  {
    return _plane;
  }
  SANDTrackerCellID id() const
  {
    return _id;
  }
  void size(double &h, double &w) const
  {
    h = _height;
    w = _width;
  }
  SANDWireInfo wire()
  {
    return _wire;
  }
  const SANDWireInfo& wire() const
  {
    return _wire;
  }
  double driftVelocity() const
  {
    return _driftVelocity;
  }

  double evaluateDriftRadius() const
  {
    if (_isFired)
      return _timeResponse * _driftVelocity;
    else
      return -1;
  }
};
