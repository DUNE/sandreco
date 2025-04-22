#include "Line3D.h"
#include "TVector3.h"

Line3D::Line3D(const TVector3& p, const TVector3& d)
{
  point_ = p;
  direction_ = d * (1. / d.Mag());
  u_ = TVector3(-direction_.Y(), direction_.X(), 0);
  v_ = direction_.Cross(u_);

  u_ = u_ * (1. / u_.Mag());
  v_ = v_ * (1. / v_.Mag());
}

Line3D::Line3D(const Line3D& line)
{
  point_ = line.point_;
  direction_ = line.direction_;
  u_ = line.u_;
  v_ = line.v_;
}

// To Do: this is basically the same as the one in the sandgeomanager.
//        Check if this class/function is really needed and merge if possible
double Line3D::distance(const Line3D& line1, const Line3D& line2)
{
  TVector3 p1 = line1.point_;
  TVector3 d1 = line1.direction_;
  TVector3 p2 = line2.point_;
  TVector3 d2 = line2.direction_;
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
