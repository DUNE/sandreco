#ifndef SANDTrackerCLUSTERCOLLECTION_H
#define SANDTrackerCLUSTERCOLLECTION_H

#include <vector>

#include "SANDTrackerClustersInPlane.h"

class SANDTrackerClusterCollection
{
  // The idea is to make this class static and const
 private:
  std::vector<SANDTrackerClustersInPlane> fPlanes;
  const SANDGeoManager* _sand_geo;

 public:
  SANDTrackerClusterCollection(const SANDGeoManager* sand_geo, const std::vector<SANDTrackerDigit> &digits);
  SANDTrackerClusterCollection(const SANDGeoManager* sand_geo);
  ~SANDTrackerClusterCollection(){};
  inline const SANDTrackerClustersInPlane &GetClustersInPlane(
      const SANDTrackerPlaneID &id
      ) const
  {
    return fPlanes.at(_sand_geo->GetPlaneIndex(id)());
  };
  inline const std::vector<SANDTrackerClustersInPlane> &GetPlanes() const
  {
    return fPlanes;
  };
  int GetNClusters() const;
  const SANDTrackerCluster &GetFirstDownstreamCluster();
};
#endif