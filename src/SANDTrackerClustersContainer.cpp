#include "SANDTrackerClustersContainer.h"

#include "SANDTrackerUtils.h"

// get digit coordinate according to the plane
inline TVector2 ClustersContainer::GetDigitCoord(const SANDTrackerDigit *dg) const
{
  auto plane = _sand_geo->get_plane_info(SANDTrackerCellID(dg->did));
  return getSandGeoManager()->GlobalToRotated(TVector2(dg->x, dg->y), *plane);
}

bool SANDTrackerClustersByProximity::CheckTriplet(const std::vector<SANDTrackerDigitID>& clu)
{
  for (auto& cluster:GetClusters()) {
    if(std::find(cluster.GetDigits().begin(), cluster.GetDigits().end(), clu[0]) != cluster.GetDigits().end() &&
       std::find(cluster.GetDigits().begin(), cluster.GetDigits().end(), clu[1]) != cluster.GetDigits().end() &&
       std::find(cluster.GetDigits().begin(), cluster.GetDigits().end(), clu[2]) != cluster.GetDigits().end())
      {
        return true;
      }
  }
  return false;
}

void SANDTrackerClustersByProximity::Clusterize(const std::vector<SANDTrackerDigitID>& digits)
{
  if (digits.size() > 0) {
    std::map<SANDTrackerCellID, SANDTrackerDigitID> fMap;
    std::for_each(digits.begin(), digits.end(), [&fMap](const SANDTrackerDigitID &d) {
      fMap[SANDTrackerCellID(d())] = d;
    });

    std::vector<SANDTrackerDigitID> clu;

    for (auto first_it = fMap.begin(); first_it != fMap.end(); first_it++) {
      const auto& first_cell = getSandGeoManager()->get_cell_info(first_it->first);
      for (auto second_it = fMap.begin(); second_it != fMap.end(); second_it++) {
        if (first_it == second_it) continue;
        const auto& second_cell = getSandGeoManager()->get_cell_info(second_it->first);
        for (auto third_it = fMap.begin(); third_it != fMap.end(); third_it++) {
          if (second_it == third_it || first_it == third_it) continue;
          const auto& third_cell = getSandGeoManager()->get_cell_info(third_it->first);
          
          if (first_cell->second.isAdjacent(second_cell->first) && 
              second_cell->second.isAdjacent(first_cell->first) && 
              second_cell->second.isAdjacent(third_cell->first) && 
              third_cell->second.isAdjacent(second_cell->first))
          {
            clu = {first_it->second, second_it->second, third_it->second};
            if(!CheckTriplet(clu)) {
              AddCluster(SANDTrackerCluster(getSandGeoManager(), clu));
            }
          }
        }
      }  
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

