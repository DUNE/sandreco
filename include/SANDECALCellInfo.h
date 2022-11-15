#include <TObject.h>

#ifndef SANDECALCELLINFO_H
#define SANDECALCELLINFO_H

// class for storing geometric info of the SAND ECAL cells
class SANDECALCellInfo : public TObject
{
 public:
  enum class Orient { kHorizontal, kVertical };

 private:
  int id_;              // id of the cell
  double x_;            // x position of the center of the cell
  double y_;            // y position of the center of the cell
  double z_;            // z position of the center of the cell
  double length_;       // length of the cell
  Orient orientation_;  // orientation of the cell
 public:
  SANDECALCellInfo();  // Default constructor
  SANDECALCellInfo(int id, double x, double y, double z, double length,
                   Orient orientation);  // parametric constructor

  // Setter methods for the attributes
  void id(int arg_id);
  void x(double arg_x);
  void y(double arg_y);
  void z(double arg_z);
  void length(double arg_length);
  void orientation(Orient arg_orientation);
  // Getter methods for the attributes
  int id();
  double x();
  double y();
  double z();
  double length();
  Orient orientation();

  ClassDef(SANDECALCellInfo, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class SANDECALCellInfo + ;
#endif

#endif
