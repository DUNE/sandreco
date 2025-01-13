/*Implementation of the SANDWireInfo class,
for storing geometric info of the SAND Tracker*/

#include "SANDWireInfo.h"

namespace sand_geometry
{

namespace tracker
{
// Default constructor
WireInfo::WireInfo()
{
}

// Parametric constructor
WireInfo::WireInfo(WireID arg_id, double arg_x, double arg_y, double arg_z,
                           double arg_length, ReadoutEnd arg_readout_end)
    : id_(arg_id),
      center_(TVector3(arg_x, arg_y, arg_z)),
      length_(arg_length),
      readout_end_(arg_readout_end)
{
}

// Parametric constructor
WireInfo::WireInfo(WireID arg_id, TVector3 arg_center,
                           double arg_length, ReadoutEnd arg_readout_end)
    : id_(arg_id),
      center_(arg_center),
      length_(arg_length),
      readout_end_(arg_readout_end)
{
}

// Setter methods for the attributes
void WireInfo::setId(WireID arg_id)
{
  id_ = arg_id;
}
void WireInfo::setLength(double arg_length)
{
  length_ = arg_length;
}
void WireInfo::setReadoutEnd(ReadoutEnd arg_readout_end)
{
  readout_end_ = arg_readout_end;
}

// Getter methods for the attributes
WireID WireInfo::getId() const
{
  return id_;
}
double WireInfo::getLength() const
{
  return length_;
}
WireInfo::ReadoutEnd WireInfo::getReadoutEnd() const
{
  return readout_end_;
}
} // namespace sand_geometry
} // namespace tracker
