#include "SANDTrackerCluster.h"
// #include "STTClusterTrackReco.h"
#include "SANDTrackerUtils.h"

#include "Math/Functor.h"
#include <TMath.h>

#include <numeric>

SANDTrackerClusterID SANDTrackerCluster::fCounter(0);

SANDTrackerCluster::SANDTrackerCluster(const SANDGeoManager* sand_geo, std::vector<SANDTrackerDigitID> &digits, const SANDTrackerPlaneID &id)
    : fId(fCounter++), fDigits(digits)
{
  _sand_geo = sand_geo;
  fPlane = _sand_geo->get_plane_info(id);
}

void SANDTrackerCluster::GetExtendedCluster(int offset)
{
  std::vector<ulong> ids;
  for (auto i = 0u; i < fDigits.size(); i++) {
    auto d = SANDTrackerDigitCollection::GetDigit(fDigits.at(i));
    ids.push_back(d.did);
  }

  std::sort(ids.begin(), ids.end());
 
  auto cells_in_plane = fPlane->getIdToCellMap();
  ulong const id_max = std::min(ids.back() + offset,  cells_in_plane.rbegin()->first());
  ulong const id_min = std::max(ids.front() - offset, cells_in_plane.begin()->first());

  for (ulong this_id = id_min; this_id <= id_max; this_id++) {
    if (std::find(ids.begin(), ids.end(), this_id) == ids.end()) {
      auto wire_info = _sand_geo->get_cell_info(this_id)->second.wire();

      SANDTrackerDigit extended_d;
      extended_d.did = this_id;
      extended_d.x = wire_info.center().X();
      extended_d.y = wire_info.center().Y();
      extended_d.z = wire_info.center().Z(); 
      extended_d.tdc = -1; 

      fDigits_extended.push_back(SANDTrackerDigitID(extended_d.did));
    }
  }
}