#include "SANDTrackerClustersContainer.h"

#include "SANDTrackerUtils.h"
namespace sand_reco
{
namespace tracker
{
// get digit coordinate according to the plane
inline TVector2 ClustersContainer::getDigitCoord(const sand_reco::tracker::Digit *dg) const
{
  auto plane = sand_geo_->getPlaneInfo(sand_geometry::tracker::CellID(dg->did));
  return getSandGeoManager()->globalToRotated(TVector2(dg->x, dg->y), *plane);
}

bool sand_reco::tracker::ClustersByProximity::isPermutation(const std::vector<sand_reco::tracker::DigitID>& clu)
{
  for (auto& cluster:getClusters()) {
    // Notice: this is O(n^2). is there a better way? i.e sort
    if (std::is_permutation(clu.begin(), clu.end(), cluster.getDigits().begin())) {
      return true;
    }
  }
  return false;
}

void sand_reco::tracker::ClustersByProximity::findCluster(std::vector<sand_reco::tracker::DigitID>& current_cluster, 
                                                 std::map<sand_geometry::tracker::CellID, sand_reco::tracker::DigitID>::iterator it, 
                                                 std::map<sand_geometry::tracker::CellID, sand_reco::tracker::DigitID>& fMap, 
                                                 int cluster_size) {

    for (auto next_it = fMap.begin(); next_it != fMap.end(); next_it++) {
      
      const auto& next_cell = getSandGeoManager()->getCellInfo(next_it->first);
      if (std::find(current_cluster.begin(), current_cluster.end(), next_it->first) 
            != current_cluster.end()) {
        continue;
      }

      int adjacent_count = 0;
      for (const auto& digit : current_cluster) {
        const auto& cluster_cell = getSandGeoManager()->getCellInfo(sand_geometry::tracker::CellID(digit()));

        if (cluster_cell->second.isAdjacent(next_cell->first)) {
          adjacent_count++;
        }
      }
      if (adjacent_count == 0) {
        continue;
      }

      current_cluster.push_back(next_it->second);

      if ((int)current_cluster.size() == cluster_size) {
        if (!isPermutation(current_cluster)) {
          addCluster(sand_reco::tracker::Cluster(getSandGeoManager(), current_cluster));
        }
        current_cluster.pop_back();
      } else {
        findCluster(current_cluster, next_it, fMap, cluster_size);
      }
    }
    current_cluster.pop_back();
}

void sand_reco::tracker::ClustersByProximity::clusterize(const std::vector<sand_reco::tracker::DigitID>& digits)
{
  if (digits.size() > 0) {
    std::map<sand_geometry::tracker::CellID, sand_reco::tracker::DigitID> fMap;
    std::for_each(digits.begin(), digits.end(), [&fMap](const sand_reco::tracker::DigitID &d) {
      fMap[sand_geometry::tracker::CellID(d())] = d;
    });

    // To Do: should be a a config paramenter
    int cluster_size = 2;
    for (auto it = fMap.begin(); it != fMap.end(); it++) {
      std::vector<sand_reco::tracker::DigitID> current_cluster = {it->second};
      findCluster(current_cluster, it, fMap, cluster_size);
    }
  }
}

const sand_reco::tracker::Cluster &sand_reco::tracker::ClustersByProximity::getNearestCluster(double x, double y) const
{
  // To Do: yes
  std::cout << "ERR: Calling sand_reco::tracker::ClustersByProximity::getNearestCluster(double x, double y) "
            << "but it is not implemented yet and you are getting the first cluster of the list." << std::endl;
  throw std::logic_error("Not implemented");
  return getClusters().at(0);
}

void sand_reco::tracker::ClustersInPlane::clusterize(const std::vector<sand_reco::tracker::DigitID>& digits)
{
  if (digits.size() > 0) {
    std::map<sand_geometry::tracker::CellID, sand_reco::tracker::DigitID> fMap;
    std::for_each(digits.begin(), digits.end(), [&fMap](const sand_reco::tracker::DigitID &d) {
      fMap[sand_geometry::tracker::CellID(d())] = d;
    });

    std::vector<sand_reco::tracker::DigitID> clu;
    clu.push_back(fMap.begin()->second);
    auto fThisTube = std::next(fMap.begin());

    while (fThisTube != fMap.end()) {
      if (SANDTrackerUtils::areAdjacent(fThisTube->first,
                                sand_geometry::tracker::CellID(clu.back()()))) {
        clu.push_back(fThisTube->second);
      } else {
        addCluster(sand_reco::tracker::Cluster(getSandGeoManager(), clu, plane_));
        clu.clear();
        clu.push_back(fThisTube->second);
      }
      fThisTube++;
    }
    addCluster(sand_reco::tracker::Cluster(getSandGeoManager(), clu, plane_));
  }
}

const sand_reco::tracker::Cluster &sand_reco::tracker::ClustersInPlane::getNearestCluster(double x, double y) const
{
  std::vector<double> dist;
  TVector2 pos(x, y);

  for (auto const &cl : getClusters()) {
    std::vector<double> dx;
    std::for_each(
        cl.getDigits().cbegin(), cl.getDigits().cend(),
        [this, &dx, pos](const sand_reco::tracker::DigitID &id) {
          dx.push_back(
            (pos - this->getDigitCoord(&sand_reco::tracker::DigitCollection::getDigit(id))).Mod());
        });

    dist.push_back(*std::min_element(dx.begin(), dx.end()));
  }

  return getClusters().at(
      std::distance(dist.begin(), std::min_element(dist.begin(), dist.end())));
}

} // namespace tracker
} // namespace sand_reco