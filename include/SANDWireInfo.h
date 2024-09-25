#include <TObject.h>

#ifndef SANDWireInfo_H
#define SANDWireInfo_H

// class for storing the STT tubes geometrical info
class SANDWireInfo : public TObject
{
 public:
  enum class Orient { kHorizontal, kVertical };
  enum class ReadoutEnd { kPlus, kMinus };

 private:
  int id_;                  // id of tube
  double x_;                // x position of the center of the tube
  double y_;                // y position of the center of the tube
  double z_;                // z position of the center of the tube
  double length_;           // length of the tube
  Orient orientation_;      // orientation of the tube
  ReadoutEnd readout_end_;  // end where signal are read
  double ax_;
  double ay_;
  double az_;

 public:
  SANDWireInfo();  // Default constructor
  SANDWireInfo(int id, double x, double y, double z, double length,
                  Orient orientation, ReadoutEnd readout_end);  // parametric constructor
  SANDWireInfo(int id, double x, double y, double z, double length,
                  Orient orientation, ReadoutEnd readout_end,
                  double arg_ax, double arg_ay, double arg_az);  // parametric constructor

  // Setter methods for the attributes
  void id(int arg_id);
  void x(double arg_x);
  void y(double arg_y);
  void z(double arg_z);
  void length(double arg_length);
  void orientation(Orient arg_orientation);
  void readout_end(ReadoutEnd arg_reaodut_end);
  void ax(double arg_ax);
  void ay(double arg_ay);
  void az(double arg_az);
  // Getter methods for the attributes
  int id();
  double x();
  double y();
  double z();
  double length();
  Orient orientation();
  ReadoutEnd readout_end();
  double ax();
  double ay();
  double az();

  ClassDef(SANDWireInfo, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class SANDWireInfo + ;
#endif

#endif
