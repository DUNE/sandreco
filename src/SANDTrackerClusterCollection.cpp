#include "SANDTrackerClusterCollection.h"
#include "SANDTrackerUtils.h"

#include "utils.h"

void SANDTrackerClusterCollection::ClusterProximityInPlane(const std::vector<SANDTrackerDigit>& digits) {
  std::map<SANDTrackerPlaneIndex, std::vector<SANDTrackerDigitID>> fMapDigits;
  for (auto& dg : digits) {
    SANDTrackerPlaneID plane_global_id;
    SANDTrackerCellID cell_local_id;
    _sand_geo->decode_cell_id(dg.did, plane_global_id, cell_local_id);
    fMapDigits[_sand_geo->GetPlaneIndex(plane_global_id)].push_back(SANDTrackerDigitID(dg.did));
  }

  for (auto& p:fMapDigits) {
    std::sort(p.second.begin(), p.second.end(), [](SANDTrackerDigitID a, SANDTrackerDigitID b)
                                  { return a() > b(); });
    containers.push_back(new SANDTrackerClustersInPlane(_sand_geo, SANDTrackerClustersContainerID(p.first()), p.second));
  }
}

void SANDTrackerClusterCollection::ClusterCellAdjacency(const std::vector<SANDTrackerDigit>& digits) {
  std::vector<SANDTrackerDigitID> digitIds;
  for (auto& dg : digits) {
    digitIds.push_back(SANDTrackerDigitID(dg.did));
  }
  containers.push_back(new SANDTrackerClustersByProximity(_sand_geo, SANDTrackerClustersContainerID(0), digitIds));
}

SANDTrackerClusterCollection::SANDTrackerClusterCollection(const SANDGeoManager* sand_geo, const std::vector<SANDTrackerDigit>& digits, ClusteringMethod clu_method)
{
  _sand_geo = sand_geo;

  if (clu_method == ClusteringMethod::kProximityInPlane) {
    ClusterProximityInPlane(digits);
  }
  if (clu_method == ClusteringMethod::kCellAdjacency) {
    ClusterCellAdjacency(digits);
  }
}

// get number of available dg_tubes
int SANDTrackerClusterCollection::GetNClusters() const
{
  auto n = 0;
  std::for_each(containers.begin(), containers.end(),
                 [&n](const ClustersContainer* p) 
                     { n += p->GetClusters().size(); });
  return n;
}

// get downstream digit
const SANDTrackerCluster &SANDTrackerClusterCollection::GetFirstDownstreamCluster()
{
  auto it = --containers.end();
  while ((*it)->GetClusters().size() == 0) --it;
  return (*it)->GetClusters().front();
}