#ifndef SANDTrackerCLUSTERSINPLANE_H
#define SANDTrackerCLUSTERSINPLANE_H

#include <vector>

#include "SANDTrackerCluster.h"

class SANDTrackerClustersInPlane
{
 private:
  const SANDTrackerPlane *fPlane;
  std::vector<SANDTrackerCluster> fClusters;
  const SANDGeoManager* _sand_geo;

  void Clusterize(const std::vector<SANDTrackerDigitID> &digits);
  // inline double GetDigitCoord(const SANDTrackerDigit *dg) const;

 public:
  SANDTrackerClustersInPlane() {};
  SANDTrackerClustersInPlane(const SANDGeoManager* sand_geo, const SANDTrackerPlaneID &id)
  {
    _sand_geo = sand_geo;
    fPlane = &(_sand_geo->get_plane_info(id)->second);
  };
  SANDTrackerClustersInPlane(const SANDGeoManager* sand_geo,
                     const SANDTrackerPlaneID &id,
                     const std::vector<SANDTrackerDigitID> &digits)
  {
    _sand_geo = sand_geo;
    fPlane = &(_sand_geo->get_plane_info(id)->second);
    Clusterize(digits);
  };
  ~SANDTrackerClustersInPlane(){};
  // const SANDTrackerPlaneIndex GetIndex() const { return fPlane->GetIndex(); };
  inline const SANDTrackerPlaneID GetId() const 
  { 
    // return fPlane->GetId(); 
  };
  // inline SANDTrackerPlane::EOrientation GetOrientation() const
  // {
    // return fPlane->GetOrientation();
  // };
  inline double GetZ() const 
  { 
    // return fPlane->GetZ(); 
  };
  inline const std::vector<SANDTrackerCluster> &GetClusters() const
  {
    return fClusters;
  };
  void AddCluster(const SANDTrackerCluster &clu) { fClusters.push_back(clu); };
  void RemoveCluster(int cid)
  {
    auto it =
        std::find_if(fClusters.begin(), fClusters.end(),
                     [cid](const SANDTrackerCluster &c) { return (cid == c.GetId()); });
    assert(it != fClusters.end());
    fClusters.erase(it);
  };
  const SANDTrackerCluster &GetNearestCluster(double x) const;
  int GetBestMatch(const SANDTrackerCluster &cl) const;
};
#endif