#ifndef SANDTrackerCLUSTERCOLLECTION_H
#define SANDTrackerCLUSTERCOLLECTION_H

#include <vector>

#include "SANDTrackerClustersContainer.h"
namespace sand_reco
{
namespace tracker
{
class ClusterCollection
{

  // The idea is to make this class static and const
  public:
  enum class ClusteringMethod {
    kProximityInPlane,
    kCellAdjacency
  };
  ClusterCollection(const SANDGeoManager* sand_geo, const std::vector<Digit> &digits, ClusteringMethod clu_method);
  ~ClusterCollection(){};

  void ClusterProximityInPlane(const std::vector<Digit>& digits);
  void ClusterCellAdjacency(const std::vector<Digit>& digits);
  inline const ClustersContainer* getClustersInContainerByIndex(const int& index) const
  {
    return containers_.at(index);
  };
  inline const ClustersContainer* getClustersInContainer(const sand_geometry::tracker::PlaneID &id) const
  {
    return containers_.at(*(sand_geo_->getPlaneIndex(id)()));
  };
  inline const std::vector<ClustersContainer*> &getContainers() const
  {
    return containers_;
  };
  int getNClusters() const;
  const Cluster &getFirstDownstreamCluster();
  
  private:
    std::vector<ClustersContainer*> containers_;
    const SANDGeoManager* sand_geo_;
};
} // namespace tracker
} // namespace sand_reco
#endif