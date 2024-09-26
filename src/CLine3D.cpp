#include "CLine3D.h"
#include "TVector3.h"

CLine3D::CLine3D(const TVector3& p, const TVector3& d)
{
  _point = p;
  _direction = d * (1. / d.Mag());
  _u = TVector3(-_direction.Y(), _direction.X(), 0);
  _v = _direction.Cross(_u);

  _u = _u * (1. / _u.Mag());
  _v = _v * (1. / _v.Mag());
}

CLine3D::CLine3D(const CLine3D& line)
{
  _point = line._point;
  _direction = line._direction;
  _u = line._u;
  _v = line._v;
}

double CLine3D::distance(const CLine3D& line1, const CLine3D& line2)
{
  TVector3 p1 = line1._point;
  TVector3 d1 = line1._direction;
  TVector3 p2 = line2._point;
  TVector3 d2 = line2._direction;
  TVector3 n = d1.Cross(d2);

  double denom = n.Mag();

  if (denom == 0) {
    // Lines are parallel
    return (p2 - p1).Cross(d1).Mag() / d1.Mag();
  } else {
    // Lines are not parallel
    return std::fabs((p2 - p1).Dot(n)) / denom;
  }
}
