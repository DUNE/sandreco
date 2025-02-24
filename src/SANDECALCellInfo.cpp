/*Implementation of the SANDECALCellInfo class,
for storing geometric info of the SAND ECAL cells*/

#include "SANDECALCellInfo.h"

namespace sand_geometry
{

namespace ecal
{

// Default constructor
ECALCellInfo::ECALCellInfo() {}

// Parametric constructor
ECALCellInfo::ECALCellInfo(int arg_id, double arg_x, double arg_y,
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
void ECALCellInfo::setId(int arg_id) { id_ = arg_id; }
void ECALCellInfo::setX(double arg_x) { x_ = arg_x; }
void ECALCellInfo::setY(double arg_y) { y_ = arg_y; }
void ECALCellInfo::setZ(double arg_z) { z_ = arg_z; }
void ECALCellInfo::setLength(double arg_length) { length_ = arg_length; }
void ECALCellInfo::setOrientation(Orient arg_orientation)
{
  orientation_ = arg_orientation;
}

// Getter methods for the attributes
int ECALCellInfo::getId() const { return id_; }
double ECALCellInfo::getX() const { return x_; }
double ECALCellInfo::getY() const { return y_; }
double ECALCellInfo::getZ() const { return z_; }
double ECALCellInfo::getLength() const { return length_; }
ECALCellInfo::Orient ECALCellInfo::getOrientation()
{
  return orientation_;
}
} // namespace ecal
} // namesoace sand_geometry