#ifndef SANDTrackerCLUSTERSINPLANE_H
#define SANDTrackerCLUSTERSINPLANE_H

#include <vector>

#include "SANDTrackerCluster.h"

class SANDTrackerClustersInPlane
{
 private:
  plane_iterator fPlane;
  std::vector<SANDTrackerCluster> fClusters;
  const SANDGeoManager* _sand_geo;

  void Clusterize(const std::vector<SANDTrackerDigitID> &digits);
  inline TVector2 GetDigitCoord(const SANDTrackerDigit *dg) const;

 public:
  SANDTrackerClustersInPlane() {};
  SANDTrackerClustersInPlane(const SANDGeoManager* sand_geo, const SANDTrackerPlaneID &id)
  {
    _sand_geo = sand_geo;
    fPlane = _sand_geo->get_plane_info(id);
  };
  SANDTrackerClustersInPlane(const SANDGeoManager* sand_geo,
                     const SANDTrackerPlaneID &id,
                     const std::vector<SANDTrackerDigitID> &digits)
  {
    _sand_geo = sand_geo;
    fPlane = _sand_geo->get_plane_info(id);
    Clusterize(digits);
  };
  ~SANDTrackerClustersInPlane(){};
  inline const SANDTrackerPlaneID GetId() const 
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
  inline std::vector<SANDTrackerCluster> &GetClusters()
  {
    return fClusters;
  };
  inline const std::vector<SANDTrackerCluster> &GetClusters() const
  {
    return fClusters;
  };
  void AddCluster(const SANDTrackerCluster &clu) { fClusters.push_back(clu); };
  void RemoveCluster(SANDTrackerClusterID cid)
  {
    auto it =
        std::find_if(fClusters.begin(), fClusters.end(),
                     [cid](const SANDTrackerCluster &c) { return (cid == c.GetId()); });
    assert(it != fClusters.end());
    fClusters.erase(it);
  };
  const SANDTrackerCluster &GetNearestCluster(double x, double y) const;
  // int GetBestMatch(const SANDTrackerCluster &cl) const;
};
#endif