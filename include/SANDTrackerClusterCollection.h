#ifndef SANDTrackerCLUSTERCOLLECTION_H
#define SANDTrackerCLUSTERCOLLECTION_H

#include <vector>

#include "SANDTrackerClustersInPlane.h"

class SANDTrackerClusterCollection
{
  // The idea is to make this class static and const
 private:
  using PlanesInModule = std::map<SANDTrackerPlaneID, SANDTrackerClustersInPlane>;
  std::map<SANDTrackerModuleID, PlanesInModule> fModules;
  const SANDGeoManager* _sand_geo;

 public:
  SANDTrackerClusterCollection(const SANDGeoManager* sand_geo, const std::vector<SANDTrackerDigit> &digits);
  SANDTrackerClusterCollection(const SANDGeoManager* sand_geo);
  ~SANDTrackerClusterCollection(){};
  inline const SANDTrackerClustersInPlane &GetClustersInPlane(
      const SANDTrackerPlaneID &id
      ) const
  {
    SANDTrackerModuleID unique_module_id;
    SANDTrackerPlaneID plane_local_id, plane_type;
    _sand_geo->decode_plane_id(id, unique_module_id, 
                               plane_local_id, plane_type);
    
    return fModules.at(unique_module_id).at(id);
  };
  inline const std::map<SANDTrackerModuleID, PlanesInModule> &GetModules() const
  {
    return fModules;
  };
  int GetNClusters() const;
  const SANDTrackerCluster &GetFirstDownstreamCluster();
};
#endif