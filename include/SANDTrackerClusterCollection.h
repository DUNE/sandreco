#ifndef SANDTrackerCLUSTERCOLLECTION_H
#define SANDTrackerCLUSTERCOLLECTION_H

#include <vector>

#include "SANDTrackerClustersContainer.h"

class SANDTrackerClusterCollection
{

  // The idea is to make this class static and const
  public:
  enum class ClusteringMethod {
    kProximityInPlane,
    kCellAdjacency
  };
  SANDTrackerClusterCollection(const SANDGeoManager* sand_geo, const std::vector<SANDTrackerDigit> &digits, ClusteringMethod clu_method);
  ~SANDTrackerClusterCollection(){};

  void ClusterProximityInPlane(const std::vector<SANDTrackerDigit>& digits);
  void ClusterCellAdjacency(const std::vector<SANDTrackerDigit>& digits);
  inline const ClustersContainer* GetClustersInContainerByIndex(const int& index) const
  {
    return containers.at(index);
  };
  inline const ClustersContainer* GetClustersInContainer(const SANDTrackerPlaneID &id) const
  {
    return containers.at(_sand_geo->GetPlaneIndex(id)());
  };
  inline const std::vector<ClustersContainer*> &GetContainer() const
  {
    return containers;
  };
  int GetNClusters() const;
  const SANDTrackerCluster &GetFirstDownstreamCluster();
  
  private:
    std::vector<ClustersContainer*> containers;
    const SANDGeoManager* _sand_geo;
};
#endif