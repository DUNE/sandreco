#include <TVector3.h>
#include <memory>
#include <map>
#include "SANDGeoManager.h"

struct AABB{
  AABB(){};
  AABB(const sand_geometry::tracker::CellID& cellID, SANDGeoManager* geo);
  void expand(const AABB& second_aabb);
  bool isOverlapping(const AABB& second_aabb, double epsilon);
  TVector3 min_;
  TVector3 max_;
};

struct Node{
  sand_geometry::tracker::CellID index_ = -1;
  std::unique_ptr<Node> left_;
  std::unique_ptr<Node> right_;
  AABB aabb_;
};

class BVH{
  public:
    BVH(){};
    BVH(std::vector<sand_geometry::tracker::CellID>& cells, SANDGeoManager* geo) {fillCellAABBMap(cells, geo); createTree(root_, cells.begin(), cells.end(), geo);};
    const std::map<sand_geometry::tracker::CellID,std::vector<sand_geometry::tracker::CellID>>& getAdjacentCells(SANDGeoManager* geo);
    
    private:
    void createTree(std::unique_ptr<Node>& node, std::vector<sand_geometry::tracker::CellID>::iterator begin, std::vector<sand_geometry::tracker::CellID>::iterator end, SANDGeoManager* geo);
    void fillCellAABBMap(std::vector<sand_geometry::tracker::CellID>& cells, SANDGeoManager* geo);
    void searchAdjacentCells(std::unique_ptr<Node>& node, std::unique_ptr<Node>& other_node, SANDGeoManager* geo);
    void getAdjacentCells(std::unique_ptr<Node>& node, SANDGeoManager* geo);
    // std::vector<std::vector<sand_geometry::tracker::CellID>> clusteredCells();
    
    std::map<sand_geometry::tracker::CellID,std::vector<sand_geometry::tracker::CellID>> cellID_to_adjacent_cells_;
    std::map<sand_geometry::tracker::CellID, AABB> cellAABBs_;
    std::unique_ptr<Node> root_ = std::make_unique<Node>();

    SANDGeoManager* geo_ ; 
};