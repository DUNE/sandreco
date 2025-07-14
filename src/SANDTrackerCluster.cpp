#include "SANDTrackerCluster.h"
// #include "STTClusterTrackReco.h"
#include "SANDTrackerUtils.h"

#include "Math/Functor.h"
#include <TMath.h>

#include <numeric>

namespace sand_reco
{
namespace tracker
{
ClusterID Cluster::counter_(0);

Cluster::Cluster(const SANDGeoManager* sand_geo, 
          std::vector<DigitID> &digits, sand_geometry::tracker::plane_iterator plane)
    : id_(counter_++), digits_(digits)
{
  sand_geo_ = sand_geo;
  plane_ = plane;
}

Cluster::Cluster(const SANDGeoManager* sand_geo, std::vector<DigitID> &digits)
    : id_(counter_++), digits_(digits)
{
  sand_geo_ = sand_geo;
  plane_ = sand_geo_->getPlaneInfo(sand_geometry::tracker::CellID(digits[0]()));

  for (const auto& digitID:digits) {
    auto plane = sand_geo_->getPlaneInfo(sand_geometry::tracker::CellID(digitID()));
    if (plane->getPosition().Z() < plane_->getPosition().Z()) {
      plane_ = plane;
    }
  }
}

void Cluster::getExtendedCluster(int offset)
{

  // To Do: add case for triplet clusters, not only plane ones
  std::vector<long> ids;
  for (auto i = 0u; i < digits_.size(); i++) {
    auto d = DigitCollection::getDigit(digits_.at(i));
    ids.push_back(d.did);
  }

  std::sort(ids.begin(), ids.end());
 
  auto cells_in_plane = plane_->getIdToCellMap();
  long const id_max = std::min(ids.back() + offset,  cells_in_plane.rbegin()->first());
  long const id_min = std::max(ids.front() - offset, cells_in_plane.begin()->first());

  for (long this_id = id_min; this_id <= id_max; this_id++) {
    if (std::find(ids.begin(), ids.end(), this_id) == ids.end()) {
      auto wire_info = sand_geo_->getCellInfo(this_id)->second.getWire();

      Digit extended_d;
      extended_d.did = this_id;
      extended_d.x = wire_info.getCenter().X();
      extended_d.y = wire_info.getCenter().Y();
      extended_d.z = wire_info.getCenter().Z(); 
      extended_d.tdc = -1; 

      digits_extended_.push_back(DigitID(extended_d.did));
    }
  }
}
} // namespace sand_reco
} // namespace tracker