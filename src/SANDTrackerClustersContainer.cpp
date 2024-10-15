#include "SANDTrackerClustersContainer.h"

#include "SANDTrackerUtils.h"

// get digit coordinate according to the plane
inline TVector2 ClustersContainer::GetDigitCoord(const SANDTrackerDigit *dg) const
{
  auto plane = _sand_geo->get_plane_info(SANDTrackerCellID(dg->did));
  return getSandGeoManager()->GlobalToRotated(TVector2(dg->x, dg->y), *plane);
}

bool SANDTrackerClustersByProximity::IsPermutation(const std::vector<SANDTrackerDigitID>& clu)
{
  for (auto& cluster:GetClusters()) {
    // Notice: this is O(n^2). is there a better way? i.e sort
    if (std::is_permutation(clu.begin(), clu.end(), cluster.GetDigits().begin())) {
      return true;
    }
  }
  return false;
}

void SANDTrackerClustersByProximity::findCluster(std::vector<SANDTrackerDigitID>& current_cluster, 
                                                 std::map<SANDTrackerCellID, SANDTrackerDigitID>::iterator it, 
                                                 std::map<SANDTrackerCellID, SANDTrackerDigitID>& fMap, 
                                                 int cluster_size) {

    for (auto next_it = fMap.begin(); next_it != fMap.end(); next_it++) {
      
      const auto& next_cell = getSandGeoManager()->get_cell_info(next_it->first);
      if (std::find(current_cluster.begin(), current_cluster.end(), next_it->first) 
            != current_cluster.end()) {
        continue;
      }

      int adjacent_count = 0;
      for (const auto& digit : current_cluster) {
        const auto& cluster_cell = getSandGeoManager()->get_cell_info(SANDTrackerCellID(digit()));

        if (cluster_cell->second.isAdjacent(next_cell->first)) {
          adjacent_count++;
        }
      }
      if (adjacent_count == 0) {
        continue;
      }

      current_cluster.push_back(next_it->second);

      if ((int)current_cluster.size() == cluster_size) {
        if (!IsPermutation(current_cluster)) {
          AddCluster(SANDTrackerCluster(getSandGeoManager(), current_cluster));
        }
        current_cluster.pop_back();
      } else {
        findCluster(current_cluster, next_it, fMap, cluster_size);
      }
    }
    current_cluster.pop_back();
}

void SANDTrackerClustersByProximity::Clusterize(const std::vector<SANDTrackerDigitID>& digits)
{
  if (digits.size() > 0) {
    std::map<SANDTrackerCellID, SANDTrackerDigitID> fMap;
    std::for_each(digits.begin(), digits.end(), [&fMap](const SANDTrackerDigitID &d) {
      fMap[SANDTrackerCellID(d())] = d;
    });

    int cluster_size = 3;
    for (auto it = fMap.begin(); it != fMap.end(); it++) {
      std::vector<SANDTrackerDigitID> current_cluster = {it->second};
      findCluster(current_cluster, it, fMap, cluster_size);
    }
  }
}

const SANDTrackerCluster &SANDTrackerClustersByProximity::GetNearestCluster(double x, double y) const
{
  // To Do: yes
  return SANDTrackerCluster();
}

void SANDTrackerClustersInPlane::Clusterize(const std::vector<SANDTrackerDigitID>& digits)
{
  if (digits.size() > 0) {
    std::map<SANDTrackerCellID, SANDTrackerDigitID> fMap;
    std::for_each(digits.begin(), digits.end(), [&fMap](const SANDTrackerDigitID &d) {
      fMap[SANDTrackerCellID(d())] = d;
    });

    std::vector<SANDTrackerDigitID> clu;
    clu.push_back(fMap.begin()->second);
    auto fThisTube = std::next(fMap.begin());

    while (fThisTube != fMap.end()) {
      if (SANDTrackerUtils::AreAdjacent(fThisTube->first,
                                SANDTrackerCellID(clu.back()()))) {
        clu.push_back(fThisTube->second);
      } else {
        AddCluster(SANDTrackerCluster(getSandGeoManager(), clu, fPlane));
        clu.clear();
        clu.push_back(fThisTube->second);
      }
      fThisTube++;
    }
    AddCluster(SANDTrackerCluster(getSandGeoManager(), clu, fPlane));
  }
}

const SANDTrackerCluster &SANDTrackerClustersInPlane::GetNearestCluster(double x, double y) const
{
  std::vector<double> dist;
  TVector2 pos(x, y);

  for (auto const &cl : GetClusters()) {
    std::vector<double> dx;
    std::for_each(
        cl.GetDigits().cbegin(), cl.GetDigits().cend(),
        [this, &dx, pos](const SANDTrackerDigitID &id) {
          dx.push_back(
            (pos - this->GetDigitCoord(&SANDTrackerDigitCollection::GetDigit(id))).Mod());
        });

    dist.push_back(*std::min_element(dx.begin(), dx.end()));
  }

  return GetClusters().at(
      std::distance(dist.begin(), std::min_element(dist.begin(), dist.end())));
}

