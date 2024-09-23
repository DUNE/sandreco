#include <TObject.h>
#include "TVector3.h"

#ifndef SANDWireInfo_H
#define SANDWireInfo_H

// class for storing the STT tubes geometrical info
class SANDWireInfo : public TObject
{
 public:
  enum class Orient {
    kHorizontal,
    kVertical
  };
  enum class ReadoutEnd {
    kPlus,
    kMinus
  };
  enum class Type {
    kSignal,
    kField
  };

 private:
  long id_;                  // id of tube
  double x_;                // x position of the center of the tube
  double y_;                // y position of the center of the tube
  double z_;                // z position of the center of the tube
  TVector3 center_;
  double length_;           // length of the tube
  Orient orientation_;      // orientation of the tube
  ReadoutEnd readout_end_;  // end where signal are read
  Type type_;
  double ax_;
  double ay_;
  double az_;
  std::vector<TVector3> points;

 public:
  SANDWireInfo();  // Default constructor
  SANDWireInfo(long id, double x, double y, double z, double length,
               Orient orientation,
               ReadoutEnd readout_end);  // parametric constructor
  SANDWireInfo(long id, double x, double y, double z, double length,
               Orient orientation, ReadoutEnd readout_end, double arg_ax,
               double arg_ay, double arg_az);  // parametric constructor

  // Setter methods for the attributes
  void id(long arg_id);
  void x(double arg_x);
  void y(double arg_y);
  void z(double arg_z);
  void center(TVector3 c) {center_ = c;};
  void length(double arg_length);
  void orientation(Orient arg_orientation);
  void readout_end(ReadoutEnd arg_reaodut_end);
  void type(Type t) {type_ = t;};
  void ax(double arg_ax);
  void ay(double arg_ay);
  void az(double arg_az);
  void setPoint(TVector3 p) {points.push_back(p);};
  // Getter methods for the attributes
  long id() const;
  double x() const;
  double y() const;
  double z() const;
  TVector3 center() const {return center_;};
  double length() const;
  Orient orientation() const;
  ReadoutEnd readout_end() const;
  Type type() const {return type_;};
  double ax() const;
  double ay() const;
  double az() const;
  std::vector<TVector3> getPoints() {return points;};
  const std::vector<TVector3> getPoints() const {return points;};
  const TVector3 getDirection() const {return (points[0] - points[1]) * (1. / (points[0] - points[1]).Mag());};

  ClassDef(SANDWireInfo, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class SANDWireInfo + ;
#endif

#endif
