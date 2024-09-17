#pragma once

#include "CLine3D.h"

class SANDTrackerCell
{
  CLine3D _wire;
  double _timeResponse;   // time response in nanoseconds
  double _driftVelocity;  // drift velocity in um/ns
  bool _isFired;
  long _wireID;
  double _width;
  double _height;

 public:
  SANDTrackerCell(const CLine3D &l, const double time, const double vd,
                  const bool fired, const long wID, const double w,
                  const double h)
      : _wire(l),
        _timeResponse(time),
        _driftVelocity(vd),
        _isFired(fired),
        _wireID(wID),
        _width(w),
        _height(h)

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
  CLine3D wire() const
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
