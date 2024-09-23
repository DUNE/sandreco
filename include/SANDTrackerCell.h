#pragma once

#include "SANDWireInfo.h"

class SANDTrackerCell
{
  SANDWireInfo _wire;
  long _wireID;
  double _width;
  double _height;

  double _timeResponse;   // time response in nanoseconds
  double _driftVelocity;  // drift velocity in um/ns
  bool _isFired;

 public:
  SANDTrackerCell() {};
  SANDTrackerCell(const SANDWireInfo &l, const long wID, const double w, const double h, 
                  const double time, const double vd, const bool fired)
      : _wire(l),
        _wireID(l.id()),
        _width(w),
        _height(h),
        _timeResponse(time),
        _driftVelocity(vd),
        _isFired(fired)

  {
  }
  
  SANDTrackerCell(const SANDWireInfo &l, 
                  const double w,
                  const double h,
                  const double v)
      : _wire(l),
        _wireID(l.id()),
        _width(w),
        _height(h),
        _driftVelocity(v)

  {
  }

  SANDTrackerCell(const SANDWireInfo &l): _wire(l),
        _wireID(l.id())
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
  void id(long wID)
  {
    _wireID = wID;
  }

  double timeResponse() const
  {
    return _timeResponse;
  }
  bool isFired() const
  {
    return _isFired;
  }
  long id() const
  {
    return _wireID;
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
