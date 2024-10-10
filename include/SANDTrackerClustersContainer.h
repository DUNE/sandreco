#ifndef SANDTrackerCLUSTERSINPLANE_H
#define SANDTrackerCLUSTERSINPLANE_H

#include <vector>

#include "SANDTrackerCluster.h"

class SANDTrackerClustersContainerID : public SingleElStruct<unsigned long>
{
 public:
  SANDTrackerClustersContainerID(unsigned long id) : SingleElStruct<unsigned long>(id){};
  SANDTrackerClustersContainerID() : SingleElStruct<unsigned long>(){};
};

class ClustersContainer 
{
  private:
    std::vector<SANDTrackerCluster> fClusters;
    const SANDGeoManager* _sand_geo;
    SANDTrackerClustersContainerID _id;
  
  public:
    virtual void Clusterize(const std::vector<SANDTrackerDigitID> &digits) = 0;
    virtual ~ClustersContainer(){};
    ClustersContainer() {};
    ClustersContainer(const SANDGeoManager* sand_geo, SANDTrackerClustersContainerID id) : _sand_geo(sand_geo), _id(id) {};

    virtual const SANDTrackerCluster &GetNearestCluster(double x, double y) const = 0;
    void AddCluster(const SANDTrackerCluster &clu) { fClusters.push_back(clu); };


    const SANDGeoManager* getSandGeoManager() const 
    {
      return _sand_geo;
    }
   
    inline const SANDTrackerClustersContainerID GetId() {return _id;};
    
    inline std::vector<SANDTrackerCluster> &GetClusters()
    {
      return fClusters;
    };
    inline const std::vector<SANDTrackerCluster> &GetClusters() const
    {
      return fClusters;
    };
    inline TVector2 GetDigitCoord(const SANDTrackerDigit *dg) const;

  protected:
    void RemoveCluster(SANDTrackerClusterID cid)
    {
      auto it =
          std::find_if(fClusters.begin(), fClusters.end(),
                      [cid](const SANDTrackerCluster &c) { return (cid == c.GetId()); });
      assert(it != fClusters.end());
      fClusters.erase(it);
    };
};

class SANDTrackerClustersByProximity : public ClustersContainer
{
 private:
  void Clusterize(const std::vector<SANDTrackerDigitID> &digits) override;

 public:
  SANDTrackerClustersByProximity() {};
  SANDTrackerClustersByProximity(const SANDGeoManager* sand_geo, const SANDTrackerClustersContainerID &id) : ClustersContainer(sand_geo, id) {};
  SANDTrackerClustersByProximity(const SANDGeoManager* sand_geo, const SANDTrackerClustersContainerID &id, const std::vector<SANDTrackerDigitID> &digits) 
    : ClustersContainer(sand_geo, id)
  {
    Clusterize(digits);
  };
  ~SANDTrackerClustersByProximity(){};
  
  bool IsPermutation(const std::vector<SANDTrackerDigitID>& clu);
  const SANDTrackerCluster &GetNearestCluster(double x, double y) const override;
  void findCluster(std::vector<SANDTrackerDigitID>& current_cluster, 
                                                 std::map<SANDTrackerCellID, SANDTrackerDigitID>::iterator it, 
                                                 std::map<SANDTrackerCellID, SANDTrackerDigitID>& fMap, 
                                                 int cluster_size);
};

class SANDTrackerClustersInPlane : public ClustersContainer
{
 private:
    plane_iterator fPlane;
    void Clusterize(const std::vector<SANDTrackerDigitID> &digits) override;

 public:
  SANDTrackerClustersInPlane() {};
  SANDTrackerClustersInPlane(const SANDGeoManager* sand_geo, const SANDTrackerClustersContainerID &id) : ClustersContainer(sand_geo, id), fPlane(getSandGeoManager()->get_plane_info(SANDTrackerPlaneID(id()))) {};
  SANDTrackerClustersInPlane(const SANDGeoManager* sand_geo, const SANDTrackerClustersContainerID &id, const std::vector<SANDTrackerDigitID> &digits) 
    : ClustersContainer(sand_geo, id), fPlane(getSandGeoManager()->get_plane_info(SANDTrackerPlaneID(id())))
  {
    Clusterize(digits);
  };
  ~SANDTrackerClustersInPlane(){};

  inline double GetRotation() const
  {
    return fPlane->getRotation();
  };
  inline double GetZ() const 
  { 
    return fPlane->getPosition().Z();
  };
  
  const SANDTrackerCluster &GetNearestCluster(double x, double y) const override;  
};
#endif