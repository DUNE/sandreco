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
                                   ModuleType arg_module_type)
    : id_(arg_id),
      x_(arg_x),
      y_(arg_y),
      z_(arg_z),
      length_(arg_length),
      module_type_(arg_module_type)
{
}

// Setter methods for the attributes
void ECALCellInfo::setId(int arg_id) { id_ = arg_id; }
void ECALCellInfo::setX(double arg_x) { x_ = arg_x; }
void ECALCellInfo::setY(double arg_y) { y_ = arg_y; }
void ECALCellInfo::setZ(double arg_z) { z_ = arg_z; }
void ECALCellInfo::setLength(double arg_length) { length_ = arg_length; }
void ECALCellInfo::setModuleType(ModuleType arg_module_type)
{
  module_type_ = arg_module_type;
}

// Getter methods for the attributes
int ECALCellInfo::getId() const { return id_; }
double ECALCellInfo::getX() const { return x_; }
double ECALCellInfo::getY() const { return y_; }
double ECALCellInfo::getZ() const { return z_; }
double ECALCellInfo::getLength() const { return length_; }
ECALCellInfo::ModuleType ECALCellInfo::getModuleType()
{
  return module_type_;
}
} // namespace ecal
} // namesoace sand_geometry