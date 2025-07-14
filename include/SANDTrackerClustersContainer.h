#ifndef SANDTrackerCLUSTERSINPLANE_H
#define SANDTrackerCLUSTERSINPLANE_H

#include <vector>

#include "SANDTrackerCluster.h"
namespace sand_reco
{
namespace tracker
{
class ClustersContainerID : public SingleElStruct<long>
{
 public:
  ClustersContainerID(long id) : SingleElStruct<long>(id){};
  ClustersContainerID() : SingleElStruct<long>(){};
};

class ClustersContainer 
{
  private:
    std::vector<Cluster> clusters_;
    const SANDGeoManager* sand_geo_;
    ClustersContainerID id_;
  
  public:
    virtual void clusterize(const std::vector<DigitID> &digits) = 0;
    virtual ~ClustersContainer(){};
    ClustersContainer() {};
    ClustersContainer(const SANDGeoManager* sand_geo, ClustersContainerID id) : sand_geo_(sand_geo), id_(id) {};

    virtual const Cluster &getNearestCluster(double x, double y) const = 0;
    void addCluster(const Cluster &clu) { clusters_.push_back(clu); };


    const SANDGeoManager* getSandGeoManager() const 
    {
      return sand_geo_;
    }
   
    inline const ClustersContainerID getId() {return id_;};
    
    inline std::vector<Cluster> &getClusters()
    {
      return clusters_;
    };
    inline const std::vector<Cluster> &getClusters() const
    {
      return clusters_;
    };
    inline TVector2 getDigitCoord(const Digit *dg) const;

  protected:
    void removeCluster(ClusterID cid)
    {
      auto it =
          std::find_if(clusters_.begin(), clusters_.end(),
                      [cid](const Cluster &c) { return (cid == c.getId()); });
      assert(it != clusters_.end());
      clusters_.erase(it);
    };
};

class ClustersByProximity : public ClustersContainer
{
 private:
  void clusterize(const std::vector<DigitID> &digits) override;

 public:
  ClustersByProximity() {};
  ClustersByProximity(const SANDGeoManager* sand_geo, const ClustersContainerID &id) : ClustersContainer(sand_geo, id) {};
  ClustersByProximity(const SANDGeoManager* sand_geo, const ClustersContainerID &id, const std::vector<DigitID> &digits) 
    : ClustersContainer(sand_geo, id)
  {
    clusterize(digits);
  };
  ~ClustersByProximity(){};
  
  bool isPermutation(const std::vector<DigitID>& clu);
  const Cluster &getNearestCluster(double x, double y) const override;
  void findCluster(std::vector<DigitID>& current_cluster, 
                                                 std::map<sand_geometry::tracker::CellID, DigitID>::iterator it, 
                                                 std::map<sand_geometry::tracker::CellID, DigitID>& fMap, 
                                                 int cluster_size);
};

class ClustersInPlane : public ClustersContainer
{
 private:
    sand_geometry::tracker::plane_iterator plane_;
    void clusterize(const std::vector<DigitID> &digits) override;

 public:
  ClustersInPlane() {};
  ClustersInPlane(const SANDGeoManager* sand_geo, const ClustersContainerID &id) : ClustersContainer(sand_geo, id), plane_(getSandGeoManager()->getPlaneInfo(sand_geometry::tracker::PlaneID(id()))) {};
  ClustersInPlane(const SANDGeoManager* sand_geo, const ClustersContainerID &id, const std::vector<DigitID> &digits) 
    : ClustersContainer(sand_geo, id), plane_(getSandGeoManager()->getPlaneInfo(sand_geometry::tracker::PlaneID(id())))
  {
    clusterize(digits);
  };
  ~ClustersInPlane(){};

  inline double getRotation() const
  {
    return plane_->getRotation();
  };
  inline double getZ() const 
  { 
    return plane_->getPosition().Z();
  };
  
  const Cluster &getNearestCluster(double x, double y) const override;  
};
} // namespace tracker
} // namespace sand_reco
#endif