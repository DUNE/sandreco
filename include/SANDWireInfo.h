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

namespace sand_geometry
{

namespace tracker
{

class WireID : public SingleElStruct<unsigned long>
{
 public:
  WireID(unsigned long id) : SingleElStruct<unsigned long>(id){};
  WireID() : SingleElStruct<unsigned long>(){};
};


// class for storing the STT tubes geometrical info
class WireInfo : public TObject
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
  WireID id_;                  // id of tube
  TVector3 center_;
  double length_;           // length of the tube
  ReadoutEnd readout_end_;  // end where signal are read
  Type type_;
  std::vector<TVector3> points_;

 public:
  WireInfo();  // Default constructor
  WireInfo(WireID id, double x, double y, double z, double length,
               ReadoutEnd readout_end);  // parametric constructor
  WireInfo(WireID id, TVector3 center, double length,
               ReadoutEnd readout_end);  // parametric constructor

  // Setter methods for the attributes
  void setId(WireID arg_id);
  void setX(double arg_x);
  void setY(double arg_y);
  void setZ(double arg_z);
  void setCenter(TVector3 c) {center_ = c;};
  void setLength(double arg_length);
  void setReadoutEnd(ReadoutEnd arg_reaodut_end);
  void setType(Type t) {type_ = t;};
  void setPoint(TVector3 p) {points_.push_back(p);};
  // getter methods for the attributes
  WireID getId() const;
  TVector3 getCenter() const {return center_;};
  double getLength() const;
  ReadoutEnd getReadoutEnd() const;
  Type getType() const {return type_;};
  std::vector<TVector3> getPoints() {return points_;};
  const std::vector<TVector3> getPoints() const {return points_;};
  const TVector3 getDirection() const {return (getOppositePointToReadout() - getReadoutPoint());};
  const TVector3 getNormalizedDirection() const {return (getOppositePointToReadout() - getReadoutPoint()) * (1. / (points_[1] - points_[0]).Mag());};
  const TVector3 getFirstPoint()  const {return points_[0];};
  const TVector3 getSecondPoint() const {return points_[1];};
  const TVector3& getReadoutPoint() const { return (readout_end_ == ReadoutEnd::kFirst) ? points_[0] : points_[1];};
  const TVector3& getOppositePointToReadout() const { return (readout_end_ == ReadoutEnd::kFirst) ? points_[1] : points_[0];};
  ClassDef(WireInfo, 1);
};

} // namespace sand_geometry
} // namespace tracker
#ifdef __MAKECINT__
#pragma link C++ class sand_geometry::tracker::WireInfo + ;
#endif

#endif
