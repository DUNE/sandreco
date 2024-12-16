/*Implementation of the SANDWireInfo class,
for storing geometric info of the SAND Tracker*/

#include "SANDWireInfo.h"

// Default constructor
SANDWireInfo::SANDWireInfo()
{
}

// Parametric constructor
SANDWireInfo::SANDWireInfo(SANDWireID arg_id, double arg_x, double arg_y, double arg_z,
                           double arg_length, ReadoutEnd arg_readout_end)
    : id_(arg_id),
      center_(TVector3(arg_x, arg_y, arg_z)),
      length_(arg_length),
      readout_end_(arg_readout_end)
{
}

// Parametric constructor
SANDWireInfo::SANDWireInfo(SANDWireID arg_id, TVector3 arg_center,
                           double arg_length, ReadoutEnd arg_readout_end)
    : id_(arg_id),
      center_(arg_center),
      length_(arg_length),
      readout_end_(arg_readout_end)
{
}

// Setter methods for the attributes
void SANDWireInfo::id(SANDWireID arg_id)
{
  id_ = arg_id;
}
void SANDWireInfo::length(double arg_length)
{
  length_ = arg_length;
}
void SANDWireInfo::readout_end(ReadoutEnd arg_readout_end)
{
  readout_end_ = arg_readout_end;
}

// Getter methods for the attributes
SANDWireID SANDWireInfo::id() const
{
  return id_;
}
double SANDWireInfo::length() const
{
  return length_;
}
SANDWireInfo::ReadoutEnd SANDWireInfo::readout_end() const
{
  return readout_end_;
}