#pragma once

#include "TVector3.h"

class CLine3D {
private:
  TVector3 _point;
  TVector3 _direction;
  TVector3 _u;
  TVector3 _v;

public:
  CLine3D() {};
  CLine3D(const TVector3 &p, const TVector3 &d);

  CLine3D(const CLine3D &line);

  static double distance(const CLine3D &line1, const CLine3D &line2);
  TVector3 getPoint() const { return _point; }
  TVector3 getDirection() const { return _direction; }
  TVector3 getU() const { return _u; }
  TVector3 getV() const { return _v; }
};