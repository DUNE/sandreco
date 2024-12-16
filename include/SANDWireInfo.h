#include <TObject.h>
#include "TVector3.h"

#ifndef SANDWireInfo_H
#define SANDWireInfo_H

template <typename T>
struct SingleElStruct {
 protected:
  T el;

 public:
  SingleElStruct() = default;
  SingleElStruct(T val) { el = val; };
  SingleElStruct(const SingleElStruct &other) { el = other.el; };
  inline T operator()() const { return el; };
  inline bool operator<(const SingleElStruct &other) const
  {
    return el < other.el;
  };
  inline bool operator>(const SingleElStruct &other) const
  {
    return el > other.el;
  };
  inline bool operator==(const SingleElStruct &other) const
  {
    return el == other.el;
  };
  inline bool operator!=(const SingleElStruct &other) const
  {
    return el != other.el;
  };
  inline int operator-(const SingleElStruct &other) const
  {
    return el - other.el;
  }
  inline T operator++(int)
  {
    return el++;
  }
};

class SANDWireID : public SingleElStruct<unsigned long>
{
 public:
  SANDWireID(unsigned long id) : SingleElStruct<unsigned long>(id){};
  SANDWireID() : SingleElStruct<unsigned long>(){};
};


// class for storing the STT tubes geometrical info
class SANDWireInfo : public TObject
{
 public:
  enum class ReadoutEnd {
    kFirst,
    kSecond
  };
  enum class Type {
    kSignal,
    kField
  };

 private:
  SANDWireID id_;                  // id of tube
  TVector3 center_;
  double length_;           // length of the tube
  ReadoutEnd readout_end_;  // end where signal are read
  Type type_;
  std::vector<TVector3> points;

 public:
  SANDWireInfo();  // Default constructor
  SANDWireInfo(SANDWireID id, double x, double y, double z, double length,
               ReadoutEnd readout_end);  // parametric constructor
  SANDWireInfo(SANDWireID id, TVector3 center, double length,
               ReadoutEnd readout_end);  // parametric constructor

  // Setter methods for the attributes
  void id(SANDWireID arg_id);
  void x(double arg_x);
  void y(double arg_y);
  void z(double arg_z);
  void center(TVector3 c) {center_ = c;};
  void length(double arg_length);
  void readout_end(ReadoutEnd arg_reaodut_end);
  void type(Type t) {type_ = t;};
  void setPoint(TVector3 p) {points.push_back(p);};
  // Getter methods for the attributes
  SANDWireID id() const;
  TVector3 center() const {return center_;};
  double length() const;
  ReadoutEnd readout_end() const;
  Type type() const {return type_;};
  std::vector<TVector3> getPoints() {return points;};
  const std::vector<TVector3> getPoints() const {return points;};
  const TVector3 getDirection() const {return (getOppositePointToReadout() - getReadoutPoint());};
  const TVector3 getNormalizedDirection() const {return (getOppositePointToReadout() - getReadoutPoint()) * (1. / (points[1] - points[0]).Mag());};
  const TVector3 getFirstPoint()  const {return points[0];};
  const TVector3 getSecondPoint() const {return points[1];};
  const TVector3& getReadoutPoint() const { return (readout_end_ == ReadoutEnd::kFirst) ? points[0] : points[1];};
  const TVector3& getOppositePointToReadout() const { return (readout_end_ == ReadoutEnd::kFirst) ? points[1] : points[0];};
  ClassDef(SANDWireInfo, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class SANDWireInfo + ;
#endif

#endif
