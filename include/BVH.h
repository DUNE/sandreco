#include <TVector3.h>
#include <memory>
#include <map>
#include "SANDGeoManager.h"

struct AABB{
  AABB(){};
  AABB(const sand_geometry::tracker::CellID& cellID, SANDGeoManager* geo);
  void expand(const AABB& second_aabb);
  TVector3 min_;
  TVector3 max_;
};

struct Node{
  sand_geometry::tracker::CellID index_;
  std::unique_ptr<Node> left_;
  std::unique_ptr<Node> right_;
  AABB aabb_;
};

class BVH{
  public:
  BVH(){};
  std::unique_ptr<Node> root_ = std::make_unique<Node>();
  BVH(std::vector<sand_geometry::tracker::CellID>& cells, SANDGeoManager* geo){fillCellAABBMap(cells, geo); createTree(root_, cells.begin(), cells.end(), geo);};
  void createTree(std::unique_ptr<Node>& node, std::vector<sand_geometry::tracker::CellID>::iterator begin, std::vector<sand_geometry::tracker::CellID>::iterator end, SANDGeoManager* geo);


  void fillCellAABBMap(std::vector<sand_geometry::tracker::CellID>& cells, SANDGeoManager* geo);
  std::map<sand_geometry::tracker::CellID, AABB> cellAABBs_;
};