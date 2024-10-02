#include "SANDTrackerClusterCollection.h"
#include "SANDTrackerUtils.h"

#include "utils.h"

SANDTrackerClusterCollection::SANDTrackerClusterCollection(const SANDGeoManager* sand_geo, const std::vector<SANDTrackerDigit>& digits)
{
  _sand_geo = sand_geo;
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
    fPlanes.push_back(SANDTrackerClustersInPlane(_sand_geo, _sand_geo->get_planes().at(p.first()).uid(), p.second));
  }
}

SANDTrackerClusterCollection::SANDTrackerClusterCollection(const SANDGeoManager* sand_geo)
{
  for (const auto& p:_sand_geo->get_planes()) {
    fPlanes.push_back(SANDTrackerClustersInPlane(_sand_geo, p.uid()));
  }
}

// get number of available dg_tubes
int SANDTrackerClusterCollection::GetNClusters() const
{
  auto n = 0;
  std::for_each(fPlanes.begin(), fPlanes.end(),
                 [&n](const SANDTrackerClustersInPlane& p) 
                     { n += p.GetClusters().size(); });
  return n;
}

// get downstream digit
const SANDTrackerCluster &SANDTrackerClusterCollection::GetFirstDownstreamCluster()
{
  auto it = --fPlanes.end();
  while (it->GetClusters().size() == 0) --it;
  return it->GetClusters().front();
}