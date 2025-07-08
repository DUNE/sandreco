#include <TVector3.h>
#include <memory>
#include "SANDGeoManager.h"

struct AABB{
  AABB(){};
  AABB(const std::vector<sand_geometry::tracker::CellID>& cells, SANDGeoManager* geo);

  TVector3 min_;
  TVector3 max_;
};

struct Node{
  std::vector<int> indices_;
  std::unique_ptr<Node> left_;
  std::unique_ptr<Node> right_;
  AABB aabb;
};

class BVH{
  public:
  BVH(){};
  Node root_;
  void createTree(const std::vector<sand_geometry::tracker::CellID>& cells, SANDGeoManager* geo);

};