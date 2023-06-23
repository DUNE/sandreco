/*Implementation of the SANDWireInfo class,
for storing geometric info of the SAND STTs*/

#include "SANDSTTTubeInfo.h"

// Default constructor
SANDWireInfo::SANDWireInfo() {}

// Parametric constructor
SANDWireInfo::SANDWireInfo(int arg_id, double arg_x, double arg_y,
                                 double arg_z, double arg_length,
                                 Orient arg_orientation,
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

// Setter methods for the attributes
void SANDWireInfo::id(int arg_id) { id_ = arg_id; }
void SANDWireInfo::x(double arg_x) { x_ = arg_x; }
void SANDWireInfo::y(double arg_y) { y_ = arg_y; }
void SANDWireInfo::z(double arg_z) { z_ = arg_z; }
void SANDWireInfo::length(double arg_length) { length_ = arg_length; }
void SANDWireInfo::orientation(Orient arg_orientation)
{
  orientation_ = arg_orientation;
}
void SANDWireInfo::readout_end(ReadoutEnd arg_readout_end)
{
  readout_end_ = arg_readout_end;
}

// Getter methods for the attributes
int SANDWireInfo::id() { return id_; }
double SANDWireInfo::x() { return x_; }
double SANDWireInfo::y() { return y_; }
double SANDWireInfo::z() { return z_; }
double SANDWireInfo::length() { return length_; }
SANDWireInfo::Orient SANDWireInfo::orientation() { return orientation_; }
SANDWireInfo::ReadoutEnd SANDWireInfo::readout_end()
{
  return readout_end_;
}