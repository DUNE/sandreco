#ifndef SANDTrackerCLUSTER_H
#define SANDTrackerCLUSTER_H

#include <vector>

#include "SANDTrackerDigitCollection.h"
#include "SANDGeoManager.h"
namespace sand_reco
{
namespace tracker
{
class ClusterID : public SingleElStruct<unsigned long>
{
 public:
  ClusterID(unsigned long id) : SingleElStruct<unsigned long>(id){};
  ClusterID() : SingleElStruct<unsigned long>(){};
};

class Cluster
{
 private:
  ClusterID id_;
  sand_geometry::tracker::plane_iterator plane_;
  std::vector<DigitID> digits_;
  std::vector<DigitID> digits_extended_;

  const SANDGeoManager* sand_geo_;

  static ClusterID counter_;


 public:
  Cluster() = default;
  Cluster(const SANDGeoManager* sand_geo, std::vector<DigitID> &digits);
  Cluster(const SANDGeoManager* sand_geo, std::vector<DigitID> &digits, sand_geometry::tracker::plane_iterator plane);
  enum class RecoAlgo { ELikelihood, EMinuit };
  inline ClusterID getId() const { return id_; };
  inline sand_geometry::tracker::plane_iterator getPlane() const {return plane_;};
  inline sand_geometry::tracker::PlaneID getPlaneId() const 
  { 
    return plane_->uId(); 
  };
  inline double getRotation() const
  {
    return plane_->getRotation();
  };
  inline double getZ() const 
  { 
    return plane_->getPosition().Z(); 
  };
  inline const SANDGeoManager* getSandGeoManager() const
  {
    return sand_geo_;
  }
  inline const std::vector<DigitID> &getDigits() const { return digits_; };
  void getExtendedCluster(int offset);
  inline const std::vector<DigitID> &getExtendedDigits() const { return digits_extended_; };
  static void resetCounter() { counter_ = 0; };

  friend class ClustersInPlane;
};
} // namespace tracker
} // namespace sand_reco

#ifdef __MAKECINT__
#pragma link C++ class ClusterTools::Point + ;
#pragma link C++ class ClusterTools::Tube + ;
#pragma link C++ class ClusterTools::TubeCollection + ;
#pragma link C++ class ClusterTools::Line + ;
#pragma link C++ class std::vector < ClusterTools::Line> + ;
#pragma link C++ class ClusterTools::ClusterParameters + ;
#pragma link C++ class ClusterTools::Cluster + ;
#pragma link C++ class ClusterTools::RecoParams + ;
#pragma link C++ class std::vector < ClusterTools::RecoParams> + ;
#pragma link C++ class ClusterTools::InputParams + ;
#pragma link C++ class Plane + ;
#pragma link C++ class SingleElStruct<unsigned int> + ;
#pragma link C++ class DigitID + ;
#pragma link C++ class PlaneIndex + ;
#pragma link C++ class std::vector < DigitID> + ;
#pragma link C++ class TubeID + ;
#pragma link C++ class Tube + ;
#pragma link C++ class std::map<TubeID,Tube> + ;
#pragma link C++ class PlaneID + ;
#pragma link C++ class Digit + ;
#pragma link C++ class Cluster + ;
#pragma link C++ class std::vector < Cluster> + ;
#endif

#endif