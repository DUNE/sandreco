#pragma once

#include "TVector3.h"

class Line3D
{
 private:
  TVector3 point_;
  TVector3 direction_;
  TVector3 u_;
  TVector3 v_;

 public:
  Line3D() {};
  Line3D(const TVector3& p, const TVector3& d);

  Line3D(const Line3D &line);

  Line3D operator =(const Line3D& line) {
    this->point_ = line.point_;
    this->direction_ = line.direction_;
    this->u_ = line.u_;
    this->v_ = line.v_;
    return *this;
}

  static double distance(const Line3D &line1, const Line3D &line2);
  TVector3 getPoint() const
  {
    return point_;
  }
  TVector3 getDirection() const
  {
    return direction_;
  }
  TVector3 getU() const
  {
    return u_;
  }
  TVector3 getV() const
  {
    return v_;
  }
};