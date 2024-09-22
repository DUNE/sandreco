/*Implementation of the SANDECALCellInfo class,
for storing geometric info of the SAND ECAL cells*/

#include "SANDECALCellInfo.h"

// Default constructor
SANDECALCellInfo::SANDECALCellInfo()
{
}

// Parametric constructor
SANDECALCellInfo::SANDECALCellInfo(int arg_id, double arg_x, double arg_y,
                                   double arg_z, double arg_length,
                                   Orient arg_orientation)
    : id_(arg_id),
      x_(arg_x),
      y_(arg_y),
      z_(arg_z),
      length_(arg_length),
      orientation_(arg_orientation)
{
}

// Setter methods for the attributes
void SANDECALCellInfo::id(int arg_id)
{
  id_ = arg_id;
}
void SANDECALCellInfo::x(double arg_x)
{
  x_ = arg_x;
}
void SANDECALCellInfo::y(double arg_y)
{
  y_ = arg_y;
}
void SANDECALCellInfo::z(double arg_z)
{
  z_ = arg_z;
}
void SANDECALCellInfo::length(double arg_length)
{
  length_ = arg_length;
}
void SANDECALCellInfo::orientation(Orient arg_orientation)
{
  orientation_ = arg_orientation;
}

// Getter methods for the attributes
int SANDECALCellInfo::id()
{
  return id_;
}
double SANDECALCellInfo::x()
{
  return x_;
}
double SANDECALCellInfo::y()
{
  return y_;
}
double SANDECALCellInfo::z()
{
  return z_;
}
double SANDECALCellInfo::length()
{
  return length_;
}
SANDECALCellInfo::Orient SANDECALCellInfo::orientation()
{
  return orientation_;
}