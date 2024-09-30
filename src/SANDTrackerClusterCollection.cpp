#include "SANDTrackerClusterCollection.h"
#include "SANDTrackerUtils.h"

#include "utils.h"

SANDTrackerClusterCollection::SANDTrackerClusterCollection(const SANDGeoManager* sand_geo, const std::vector<SANDTrackerDigit>& digits)
{
  _sand_geo = sand_geo;
  std::map<SANDTrackerPlaneID, std::vector<SANDTrackerDigitID>> fMapDigits;
  for (auto& dg : digits) {
    SANDTrackerPlaneID plane_global_id;
    SANDTrackerCellID cell_local_id;
    _sand_geo->decode_cell_id(dg.did, plane_global_id, cell_local_id);
    fMapDigits[plane_global_id].push_back(SANDTrackerDigitID(dg.did));
  }

  std::map<SANDTrackerPlaneID, SANDTrackerClustersInPlane> fPlanes;
  for (auto& p:fMapDigits) {
    std::sort(p.second.begin(), p.second.end(), [](SANDTrackerDigitID a, SANDTrackerDigitID b)
                                  { return a() > b(); });
    // fPlanes[p.first] = SANDTrackerClustersInPlane(_sand_geo, p.first, p.second);
    SANDTrackerModuleID unique_module_id;
    SANDTrackerPlaneID plane_local_id, plane_type;
    _sand_geo->decode_plane_id(p.first, unique_module_id, 
                               plane_local_id, plane_type);
    fModules[unique_module_id][p.first] = SANDTrackerClustersInPlane(_sand_geo, p.first, p.second);
  }
}

SANDTrackerClusterCollection::SANDTrackerClusterCollection(const SANDGeoManager* sand_geo)
{
  std::map<SANDTrackerPlaneID, SANDTrackerClustersInPlane> fPlanes;
  for (const auto& m:_sand_geo->get_modules_map()) {
    for (const auto& p:m.second.planes()) {
      // fPlanes[p.first] = SANDTrackerClustersInPlane(_sand_geo, p.first);
      SANDTrackerModuleID unique_module_id;
      SANDTrackerPlaneID plane_local_id, plane_type;
      _sand_geo->decode_plane_id(p.first, unique_module_id, 
                                plane_local_id, plane_type);
      fModules[unique_module_id][p.first] = SANDTrackerClustersInPlane(_sand_geo, p.first);
    }
  }
}

// get number of available dg_tubes
int SANDTrackerClusterCollection::GetNClusters() const
{
  auto n = 0;
  for (const auto& fPlanes:fModules) {
    std::for_each(fPlanes.second.begin(), fPlanes.second.end(),
                 [&n](const std::pair<SANDTrackerPlaneID, SANDTrackerClustersInPlane>& p) 
                     { n += p.second.GetClusters().size(); });
  }
  return n;
}

// get downstream digit
const SANDTrackerCluster &SANDTrackerClusterCollection::GetFirstDownstreamCluster() const
{
  // auto it = --fPlanes.end();
  // while (it->second.GetClusters().size() == 0) --it;
  // return it->second.GetClusters().front();
}

// // check if the module has available digits
// bool SANDTrackerClusterCollection::IsPlaneOK(const SANDTrackerPlaneIndex &index) const
// {
//   if (index() >= fPlanes.size())
//     return false;
//   else
//     return (fPlanes.at(index()).GetClusters().size() != 0);
// }