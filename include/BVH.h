#pragma once
#include <TVector3.h>
#include <memory>
#include <map>
#include "SANDGeoManager.h"

class SANDGeoManager;

struct AABB{
  AABB(){};
  AABB(const sand_geometry::tracker::cell_map_iterator cellID, SANDGeoManager* geo);
  void expand(const AABB& second_aabb);
  bool isOverlapping(const AABB& second_aabb, double epsilon = 0) const;
  TVector3 min_;
  TVector3 max_;
};

struct Node{
  sand_geometry::tracker::CellID index_ = std::nullopt;
  sand_geometry::tracker::cell_map_iterator cell_iterator_;
  std::unique_ptr<Node> left_;
  std::unique_ptr<Node> right_;
  AABB aabb_;
};

class BVH{
  public:
    BVH(){};
    BVH(std::vector<sand_geometry::tracker::cell_map_iterator>& cells, SANDGeoManager* geo, double max_distance, double overlap_tolerance) {fillCellAABBMap(cells, geo); createTree(root_, cells.begin(), cells.end(), geo); searchAdjacentCells(root_, root_, geo, max_distance, overlap_tolerance);};
    
    private:
    void createTree(std::unique_ptr<Node>& node, std::vector<sand_geometry::tracker::cell_map_iterator>::iterator begin, std::vector<sand_geometry::tracker::cell_map_iterator>::iterator end, SANDGeoManager* geo);
    void fillCellAABBMap(const std::vector<sand_geometry::tracker::cell_map_iterator>& cells, SANDGeoManager* geo);
    void searchAdjacentCells(std::unique_ptr<Node>& node, std::unique_ptr<Node>& other_node, SANDGeoManager* geo, double max_distance, double overlap_tolerance);
    std::map<sand_geometry::tracker::CellID, AABB> cellAABBs_;
    std::unique_ptr<Node> root_ = std::make_unique<Node>();

    SANDGeoManager* geo_ ; 
};