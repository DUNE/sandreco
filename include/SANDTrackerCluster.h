#ifndef SANDTrackerCLUSTER_H
#define SANDTrackerCLUSTER_H

#include <vector>

#include "SANDTrackerDigitCollection.h"
#include "SANDGeoManager.h"

class SANDTrackerClusterID : public SingleElStruct<unsigned long>
{
 public:
  SANDTrackerClusterID(unsigned long id) : SingleElStruct<unsigned long>(id){};
  SANDTrackerClusterID() : SingleElStruct<unsigned long>(){};
};

class SANDTrackerCluster
{
 private:
  SANDTrackerClusterID fId;
  plane_iterator fPlane;
  std::vector<SANDTrackerDigitID> fDigits;
  std::vector<SANDTrackerDigitID> fDigits_extended;

  const SANDGeoManager* _sand_geo;

  static SANDTrackerClusterID fCounter;

  SANDTrackerCluster(const SANDGeoManager* sand_geo, std::vector<SANDTrackerDigitID> &digits, const SANDTrackerPlaneID &id);

 public:
  SANDTrackerCluster() = default;
  enum class RecoAlgo { ELikelihood, EMinuit };
  inline SANDTrackerClusterID GetId() const { return fId; };
  inline SANDTrackerPlaneID GetPlaneId() const 
  { 
    return fPlane->uid(); 
  };
  inline double GetRotation() const
  {
    return fPlane->getRotation();
  };
  inline double GetZ() const 
  { 
    return fPlane->getPosition().Z(); 
  };
  inline const std::vector<SANDTrackerDigitID> &GetDigits() const { return fDigits; };
  void GetExtendedCluster(int offset);
  inline const std::vector<SANDTrackerDigitID> &GetExtendedDigits() const { return fDigits_extended; };
  static void ResetCounter() { fCounter = 0; };

  friend class SANDTrackerClustersInPlane;
};

#ifdef __MAKECINT__
#pragma link C++ class SANDTrackerClusterTools::Point + ;
#pragma link C++ class SANDTrackerClusterTools::Tube + ;
#pragma link C++ class SANDTrackerClusterTools::TubeCollection + ;
#pragma link C++ class SANDTrackerClusterTools::Line + ;
#pragma link C++ class std::vector < SANDTrackerClusterTools::Line> + ;
#pragma link C++ class SANDTrackerClusterTools::ClusterParameters + ;
#pragma link C++ class SANDTrackerClusterTools::Cluster + ;
#pragma link C++ class SANDTrackerClusterTools::RecoParams + ;
#pragma link C++ class std::vector < SANDTrackerClusterTools::RecoParams> + ;
#pragma link C++ class SANDTrackerClusterTools::InputParams + ;
#pragma link C++ class SANDTrackerPlane + ;
#pragma link C++ class SingleElStruct<unsigned int> + ;
#pragma link C++ class SANDTrackerDigitID + ;
#pragma link C++ class SANDTrackerPlaneIndex + ;
#pragma link C++ class std::vector < SANDTrackerDigitID> + ;
#pragma link C++ class SANDTrackerTubeID + ;
#pragma link C++ class SANDTrackerTube + ;
#pragma link C++ class std::map<SANDTrackerTubeID,SANDTrackerTube> + ;
#pragma link C++ class SANDTrackerPlaneID + ;
#pragma link C++ class SANDTrackerDigit + ;
#pragma link C++ class SANDTrackerCluster + ;
#pragma link C++ class std::vector < SANDTrackerCluster> + ;
#endif

#endif