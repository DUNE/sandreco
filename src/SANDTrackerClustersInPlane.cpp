#include "SANDTrackerClustersInPlane.h"

#include "SANDTrackerUtils.h"

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
        AddCluster(SANDTrackerCluster(_sand_geo, clu, GetId()));
        clu.clear();
        clu.push_back(fThisTube->second);
      }
      fThisTube++;
    }
    AddCluster(SANDTrackerCluster(_sand_geo, clu, GetId()));
  }
}

// get nearest digit in module
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

// int SANDTrackerClustersInPlane::GetBestMatch(const SANDTrackerCluster &cl) const
// {
//   auto z = 0.5 * (cl.GetZ() + GetZ());
//   auto index = -1;
//   for (auto i = 0u; i < fClusters.size(); i++) {
//     auto &this_cl = fClusters.at(i);
//     auto ok1 = cl.GetDigits().size() == 2 ||
//                (cl.GetDigits().size() > 2 &&
//                 cl.GetRecoParameters().front().quality < 10);
//     auto ok2 = this_cl.GetDigits().size() == 2 ||
//                (this_cl.GetDigits().size() > 2 &&
//                 this_cl.GetRecoParameters().front().quality < 10);
//     if (ok1 && ok2) {
//       auto min_dx = 1E200;
//       for (const auto &par1 : cl.GetRecoParameters())
//         for (const auto &par2 : this_cl.GetRecoParameters()) {
//           auto x1 = par1.trk.m * z + par1.trk.q;
//           auto x2 = par2.trk.m * z + par2.trk.q;
//           if (fabs(x1 - x2) < min_dx) {
//             auto ds = fabs(par1.trk.m - par2.trk.m);
//             auto dx = fabs(x1 - x2);
//             if (dx < 5 && ds < 0.4) {
//               min_dx = dx;
//               index = i;
//             }
//           }
//         }
//     }
//   }
//   return index;
// }

// get digit coordinate according to the plane
inline TVector2 SANDTrackerClustersInPlane::GetDigitCoord(const SANDTrackerDigit *dg) const
{
  return _sand_geo->GlobalToRotated(TVector2(dg->x, dg->y), *fPlane);
}