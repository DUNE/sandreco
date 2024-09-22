/*Implementation of the SANDWireInfo class,
for storing geometric info of the SAND Tracker*/

#include "SANDWireInfo.h"

// Default constructor
SANDWireInfo::SANDWireInfo()
{
}

// Parametric constructor
SANDWireInfo::SANDWireInfo(long arg_id, double arg_x, double arg_y, double arg_z,
                           double arg_length, Orient arg_orientation,
                           ReadoutEnd arg_readout_end)
    : id_(arg_id),
      x_(arg_x),
      y_(arg_y),
      z_(arg_z),
      length_(arg_length),
      orientation_(arg_orientation),
      readout_end_(arg_readout_end)
{
}

// Parametric constructor
SANDWireInfo::SANDWireInfo(long arg_id, double arg_x, double arg_y, double arg_z,
                           double arg_length, Orient arg_orientation,
                           ReadoutEnd arg_readout_end, double arg_ax,
                           double arg_ay, double arg_az)
    : id_(arg_id),
      x_(arg_x),
      y_(arg_y),
      z_(arg_z),
      length_(arg_length),
      orientation_(arg_orientation),
      readout_end_(arg_readout_end),
      ax_(arg_ax),
      ay_(arg_ay),
      az_(arg_az)
{
}

// Setter methods for the attributes
void SANDWireInfo::id(long arg_id)
{
  id_ = arg_id;
}
void SANDWireInfo::x(double arg_x)
{
  x_ = arg_x;
}
void SANDWireInfo::y(double arg_y)
{
  y_ = arg_y;
}
void SANDWireInfo::z(double arg_z)
{
  z_ = arg_z;
}
void SANDWireInfo::length(double arg_length)
{
  length_ = arg_length;
}
void SANDWireInfo::ax(double arg_ax)
{
  ax_ = arg_ax;
}
void SANDWireInfo::ay(double arg_ay)
{
  ay_ = arg_ay;
}
void SANDWireInfo::az(double arg_az)
{
  az_ = arg_az;
}
void SANDWireInfo::orientation(Orient arg_orientation)
{
  orientation_ = arg_orientation;
}
void SANDWireInfo::readout_end(ReadoutEnd arg_readout_end)
{
  readout_end_ = arg_readout_end;
}

// Getter methods for the attributes
long SANDWireInfo::id() const
{
  return id_;
}
double SANDWireInfo::x() const
{
  return x_;
}
double SANDWireInfo::y() const
{
  return y_;
}
double SANDWireInfo::z() const
{
  return z_;
}
double SANDWireInfo::length() const
{
  return length_;
}
SANDWireInfo::Orient SANDWireInfo::orientation() const
{
  return orientation_;
}
double SANDWireInfo::ax() const
{
  return ax_;
}
double SANDWireInfo::ay() const
{
  return ay_;
}
double SANDWireInfo::az() const
{
  return az_;
}
SANDWireInfo::ReadoutEnd SANDWireInfo::readout_end() const
{
  return readout_end_;
}