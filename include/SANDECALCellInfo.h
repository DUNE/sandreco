#include <TObject.h>

#ifndef ECALCELLINFO_H
#define ECALCELLINFO_H

namespace sand_geometry
{

namespace ecal
{

// class for storing geometric info of the SAND ECAL cells
class ECALCellInfo : public TObject
{
 public:
  enum class ModuleType { kBarrel, kEndcap };

 private:
  int id_;              // id of the cell
  double x_;            // x position of the center of the cell
  double y_;            // y position of the center of the cell
  double z_;            // z position of the center of the cell
  double length_;       // length of the cell
  ModuleType module_type_;  // ModuleType of the cell
  
 public:
  ECALCellInfo();  // Default constructor
  ECALCellInfo(int id, double x, double y, double z, double length,
                   ModuleType module_type);  // parametric constructor

  // Setter methods for the attributes
  void setId(int arg_id);
  void setX(double arg_x);
  void setY(double arg_y);
  void setZ(double arg_z);
  void setLength(double arg_length);
  void setModuleType(ModuleType arg_module_type);
  // Getter methods for the attributes
  int getId() const;
  double getX() const;
  double getY() const;
  double getZ() const;
  double getLength() const;
  ModuleType getModuleType() const;

  ClassDef(ECALCellInfo, 1);
};
} // namespace ecal
} // namesoace sand_geometry

#ifdef __MAKECINT__
#pragma link C++ class sand_geometry::ecal::ECALCellInfo + ;
#endif

#endif
