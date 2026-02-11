#include "SANDTrackerClusterCollection.h"
#include "SANDTrackerUtils.h"

#include "utils.h"

namespace sand_reco
{
namespace tracker
{
void ClusterCollection::ClusterProximityInPlane(const std::vector<sand_reco::tracker::Digit>& digits) {
  std::map<sand_geometry::tracker::PlaneIndex, std::vector<sand_reco::tracker::DigitID>> fMapDigits;
  for (auto& dg : digits) {
    sand_geometry::tracker::PlaneID plane_global_id;
    sand_geometry::tracker::CellID cell_local_id;
    sand_geo_->decodeCellId(dg.did, plane_global_id, cell_local_id);
    fMapDigits[sand_geo_->getPlaneIndex(plane_global_id)].push_back(sand_reco::tracker::DigitID(dg.did));
  }

  for (auto& p:fMapDigits) {
    std::sort(p.second.begin(), p.second.end(), [](sand_reco::tracker::DigitID a, sand_reco::tracker::DigitID b)
                                  { return a() > b(); });
    containers_.push_back(new sand_reco::tracker::ClustersInPlane(sand_geo_, sand_reco::tracker::ClustersContainerID(p.first()), p.second));
  }
}

void ClusterCollection::ClusterCellAdjacency(const std::vector<sand_reco::tracker::Digit>& digits) {
  std::map<sand_geometry::tracker::ModuleID, std::vector<sand_reco::tracker::DigitID>> fMapDigits;
  for (auto& dg : digits) {
    sand_geometry::tracker::ModuleID unique_module_id;
    sand_geometry::tracker::ModuleID supermodule_id;
    sand_geometry::tracker::ModuleID module_id;
    sand_geometry::tracker::ModuleID module_replica_id;
    sand_geometry::tracker::PlaneID plane_global_id;
    sand_geometry::tracker::PlaneID plane_replica_id;
    sand_geometry::tracker::PlaneID plane_type;
    sand_geometry::tracker::CellID cell_local_id;
    sand_geo_->decodeCellId(dg.did, plane_global_id, cell_local_id);

    sand_geo_->decodePlaneId(plane_global_id, unique_module_id, plane_replica_id, plane_type);
    sand_geo_->decodeModuleId(unique_module_id, supermodule_id, module_id, module_replica_id);
    fMapDigits[unique_module_id].push_back(sand_reco::tracker::DigitID(dg.did));
    // std::cout << *unique_module_id() << " " << sand_geo_->getPlanes().at(*(sand_geo_->getPlaneIndex(plane_global_id)())).getPosition().Z() << std::endl;
  }

  for (auto& p:fMapDigits) {
    containers_.push_back(new sand_reco::tracker::ClustersByProximity(sand_geo_, sand_reco::tracker::ClustersContainerID(p.first()), p.second));
  }
}

ClusterCollection::ClusterCollection(const SANDGeoManager* sand_geo, const std::vector<sand_reco::tracker::Digit>& digits, ClusteringMethod clu_method)
{
  sand_geo_ = sand_geo;

  if (clu_method == ClusteringMethod::kProximityInPlane) {
    ClusterProximityInPlane(digits);
  }
  if (clu_method == ClusteringMethod::kCellAdjacency) {
    ClusterCellAdjacency(digits);
  }
}

// get number of available dg_tubes
int ClusterCollection::getNClusters() const
{
  auto n = 0;
  std::for_each(containers_.begin(), containers_.end(),
                 [&n](const ClustersContainer* p) 
                     { n += p->getClusters().size(); });
  return n;
}

// get downstream digit
const sand_reco::tracker::Cluster &ClusterCollection::getFirstDownstreamCluster()
{
  auto it = --containers_.end();
  while ((*it)->getClusters().size() == 0) --it;
  return (*it)->getClusters().front();
}
} // namespace sand_reco
} // namespace tracker