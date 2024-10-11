#ifndef STTKFCLUSTERMANAGER_H
#define STTKFCLUSTERMANAGER_H

#include "SANDTrackerCluster.h"
#include "SANDTrackerClusterCollection.h"

class STTKFClusterUse {
  public:
    enum class STTKFClusterUsage {
      kUsed,
      kAvailable,
    };
    
  private:
    int fTrackID;
    STTKFClusterUsage fUse; 
    SANDTrackerCluster const* fCluster;

  public:
    STTKFClusterUse(const SANDTrackerCluster& cluster): 
      fTrackID(-1),
      fUse(STTKFClusterUsage::kAvailable), 
      fCluster(&cluster) {};
    STTKFClusterUse(): 
      fTrackID(-1),
      fUse(STTKFClusterUsage::kAvailable), 
      fCluster(nullptr) {};
    
    // void SetIsUsed(int trackID) {
    //   fTrackID = trackID;
    //   fUse = STTKFClusterUsage::kUsed;
    // };
    
    // void SetIsAvailable() {
    //   fTrackID = -1;
    //   fUse = STTKFClusterUsage::kAvailable;
    // };

    // const SANDTrackerCluster& GetCluster() const { return *fCluster; };
    // int GetTrackID() { return fTrackID; };
    // const STTKFClusterUsage& GetAvailability() {return fUse; };
};

class STTKFClusterManager {
  private:
    // static std::map<int, STTKFClusterUse> fClusterUse;
    // static const SANDTrackerClusterCollection* fClusterCollection;

  public:
    STTKFClusterManager() {};
    // static void Init(const SANDTrackerClusterCollection& clusterCollection);
    // static const SANDTrackerCluster& GetCluster(int clusterID) {return fClusterUse.at(clusterID).GetCluster(); };
    // static std::vector<int> GetClusterIDsForPlane(const STTPlaneIndex &id);
    // static void SetClusterUsed(int clusterID, int trackID);
    // static bool GetDownPlaneIDAndFirstAvailableCluster(STTPlaneID& planeID, int& clusterID);
    // static bool isClusterUsable(const SANDTrackerCluster& cluster);
    // static int GetNUsableClusters(const std::vector<SANDTrackerCluster>& clusters);
};

#endif