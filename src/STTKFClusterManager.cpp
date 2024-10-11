#include "STTKFClusterManager.h"
#include "SANDTrackerUtils.h"

// std::map<int, STTKFClusterUse> STTKFClusterManager::fClusterUse;
// const STTClusterCollection* STTKFClusterManager::fClusterCollection = 0;

// void STTKFClusterManager::Init(const STTClusterCollection& clusterCollection) {
//   fClusterCollection = &clusterCollection;
//   fClusterUse.clear();

//   // Init cluster usage
//   for(auto& plane: fClusterCollection->GetPlanes())
//   {
//     for(auto& cluster: plane.GetClusters())
//     {
//       fClusterUse[cluster.GetId()] = STTKFClusterUse(cluster);
//     }
//   }
// }

// std::vector<int> STTKFClusterManager::GetClusterIDsForPlane(const STTPlaneIndex &index) {
//   std::vector<int> clusterIDs;

//   for(auto cluster: fClusterCollection->GetClustersInPlane(index).GetClusters())
//   {
//     auto id = cluster.GetId();
//     if(STTKFClusterManager::isClusterUsable(cluster))
//     // if(STTKFClusterManager::isClusterUsable(fClusterUse[id].GetAvailability() == STTKFClusterUse::STTKFClusterUsage::kAvailable)
//     //   if(cluster.GetDigits().size() > 1)
//         clusterIDs.push_back(id);
//   }
//   STTTRACKRECO_LOG("INFO", TString::Format("Found %lu clustera in plane index %u", clusterIDs.size(), index()).Data());
//   return clusterIDs;
// }

// void STTKFClusterManager::SetClusterUsed(int clusterID, int trackID) {
//   fClusterUse[clusterID].SetIsUsed(trackID);
// }

// bool STTKFClusterManager::GetDownPlaneIDAndFirstAvailableCluster(STTPlaneID& planeID, int& clusterID) {
//   auto& planes = fClusterCollection->GetPlanes();
//   for(auto it = planes.rbegin(); it != planes.rend(); it++)
//     if(GetNUsableClusters(it->GetClusters()) != 0)    
//       for(auto& cluster: it->GetClusters())
//         if(STTKFClusterManager::isClusterUsable(cluster))
//         {
//           planeID = it->GetId();
//           clusterID = cluster.GetId();
//           auto& trk = cluster.GetRecoParameters().front().trk;
//           STTTRACKRECO_LOG("INFO", TString::Format("Found downstream cluster %8d in plane %8u with m: %f and q: %f", clusterID, planeID(), trk.m, trk.q).Data());
//           return true;
//         }
//   STTTRACKRECO_LOG("INFO", "No downstream cluster found");
//   return false;
// }

// bool STTKFClusterManager::isClusterUsable(const STTCluster& cluster) {
//   bool condition_1 = fClusterUse[cluster.GetId()].GetAvailability() == STTKFClusterUse::STTKFClusterUsage::kAvailable;
//   bool condition_2 = cluster.GetDigits().size() > 1u;
//   return condition_1 && condition_2;
// }

// int STTKFClusterManager::GetNUsableClusters(const std::vector<STTCluster>& clusters) {
//   int count = 0;
//   for(auto& cluster: clusters)
//     if(isClusterUsable(cluster))
//       count++;
//   STTTRACKRECO_LOG("INFO", TString::Format("%d usable clusters of %lu", count, clusters.size()).Data());
//   return count;
// }