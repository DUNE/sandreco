#include "BVH.h"
#include <TVector2.h>



AABB::AABB(const std::map<sand_geometry::tracker::CellID, sand_geometry::tracker::Cell>::iterator cell, SANDGeoManager* geo){
    min_ = TVector3(1E9, 1E9, 1E9);
    max_ = TVector3(-1E9, -1E9, -1E9);
    std::vector<TVector3> vertices;
    auto p = cell->second.getWire().getFirstPoint();
    
    auto& plane = *geo->getPlaneInfo(cell->first);
    auto p_rotated = geo->globalToRotated(TVector2(p.X(), p.Y()), plane);  
    
    auto h_2 = cell->second.getSize().h /2.;
    auto w_2 = cell->second.getSize().w /2.;

    TVector2 v1(p_rotated.X(), p_rotated.Y() + h_2);
    TVector2 v2(p_rotated.X(), p_rotated.Y() - h_2);
    
    auto p1 = geo->rotatedToGlobal(TVector2(v1.X(), v1.Y()), plane); 

    auto p2 = geo->rotatedToGlobal(TVector2(v2.X(), v2.Y()), plane);

    vertices.push_back(TVector3(p1.X(), p1.Y(), p.Z() + w_2));
    vertices.push_back(TVector3(p2.X(), p2.Y(), p.Z() + w_2));   

    p = cell->second.getWire().getSecondPoint();
    p_rotated = geo->globalToRotated(TVector2(p.X(), p.Y()), plane); 
    TVector2 v3(p_rotated.X(), p_rotated.Y() + h_2);
    TVector2 v4(p_rotated.X(), p_rotated.Y() - h_2);

    auto p3 = geo->rotatedToGlobal(TVector2(v3.X(), v3.Y()), plane); 
    auto p4 = geo->rotatedToGlobal(TVector2(v4.X(), v4.Y()), plane);

    vertices.push_back(TVector3(p3.X(), p3.Y(), p.Z() - w_2));
    vertices.push_back(TVector3(p4.X(), p4.Y(), p.Z() - w_2));

    for(const auto &v : vertices){
        min_.SetX(std::min(min_.X(), v.X()));
        max_.SetX(std::max(max_.X(), v.X()));
        min_.SetY(std::min(min_.Y(), v.Y()));
        max_.SetY(std::max(max_.Y(), v.Y()));
        min_.SetZ(std::min(min_.Z(), v.Z()));
        max_.SetZ(std::max(max_.Z(), v.Z()));
    }
}

void AABB::expand(const AABB& second_aabb){
    min_.SetX(std::min(min_.X(), second_aabb.min_.X()));
    max_.SetX(std::max(max_.X(), second_aabb.max_.X()));
    min_.SetY(std::min(min_.Y(), second_aabb.min_.Y()));
    max_.SetY(std::max(max_.Y(), second_aabb.max_.Y()));
    min_.SetZ(std::min(min_.Z(), second_aabb.min_.Z()));
    max_.SetZ(std::max(max_.Z(), second_aabb.max_.Z()));
}

bool AABB::isOverlapping(const AABB& second_aabb, double epsilon = 0){
    if(max_.X() + epsilon >= second_aabb.min_.X() && min_.X() - epsilon <= second_aabb.max_.X() && 
       max_.Y() + epsilon >= second_aabb.min_.Y() && min_.Y() - epsilon <= second_aabb.max_.Y() && 
       max_.Z() + epsilon >= second_aabb.min_.Z() && min_.Z() - epsilon <= second_aabb.max_.Z()){
        return true;
    }
    return false;
}

void BVH::fillCellAABBMap(std::vector<std::map<sand_geometry::tracker::CellID, sand_geometry::tracker::Cell>::iterator> cells, SANDGeoManager* geo){
    for(const auto &cell : cells){
        cellAABBs_[cell->first] = AABB(cell, geo);
    }
}

void BVH::createTree(std::unique_ptr<Node>& node, std::vector<sand_geometry::tracker::cell_map_iterator>::iterator begin, std::vector<sand_geometry::tracker::cell_map_iterator>::iterator end, SANDGeoManager* geo){

    node->aabb_ = cellAABBs_[(*begin)->first];
    for (auto it = begin + 1; it != end; ++it) {
        node->aabb_.expand(cellAABBs_[(*it)->first]);
    }

    node->index_ = sand_geometry::tracker::CellID(-1);
    
    if(std::distance(begin, end) == 1){
        node->index_ = (*begin)->first;
        node->cell_iterator_ = (*begin);
        return;
    }
    
    auto sorting_by_z = [geo](const sand_geometry::tracker::cell_map_iterator c1, const sand_geometry::tracker::cell_map_iterator c2){
        const auto& center_1 = c1->second.getWire().getCenter();
        const auto& center_2 = c2->second.getWire().getCenter();
        
        return center_1.Z() < center_2.Z();
    };
    auto sorting_by_y = [geo](const sand_geometry::tracker::cell_map_iterator c1, const sand_geometry::tracker::cell_map_iterator c2){
        const auto& center_1 = c1->second.getWire().getCenter();
        const auto& center_2 = c2->second.getWire().getCenter();
        
        return center_1.Y() < center_2.Y();
    };
    auto sorting_by_x = [geo](const sand_geometry::tracker::cell_map_iterator c1, const sand_geometry::tracker::cell_map_iterator c2){
        const auto& center_1 = c1->second.getWire().getCenter();
        const auto& center_2 = c2->second.getWire().getCenter();
        
        return center_1.X() < center_2.X();
    };
    
    auto sorting_by_id = [geo](const sand_geometry::tracker::cell_map_iterator c1, const sand_geometry::tracker::cell_map_iterator c2){
        const auto& id1 = c1->second.getId();
        const auto& id2 = c2->second.getId();
        
        return id1 < id2;
    };
    
    double deltaZ = node->aabb_.max_.Z() - node->aabb_.min_.Z();
    int middle_point= std::distance(begin, end) / 2;
    const double w = (*begin)->second.getSize().w;
    if (w != deltaZ) {
        std::sort(begin, end, sorting_by_z);
        // std::nth_element(begin, begin + middle_point, end, sorting_by_z);
    } else {
        std::sort(begin, end, sorting_by_id);
        // std::nth_element(begin, begin + middle_point, end, sorting_by_id);
    }
    node->left_ = std::make_unique<Node>();
    createTree(node->left_, begin, begin + middle_point, geo);

    node->right_ = std::make_unique<Node>();
    createTree(node->right_, begin + middle_point, end, geo);
 
    return;
}


void BVH::searchAdjacentCells(std::unique_ptr<Node>& node, std::unique_ptr<Node>& other_node, SANDGeoManager* geo){
    if(!node || !other_node) return;
    
    if(!node->aabb_.isOverlapping(other_node->aabb_, 1)) return;

    if(other_node->index_ != -1 && node->index_ != -1) {
        if (node->index_ == 220006 || other_node->index_ == 220006) {
            std::cout << "comparing " << node->index_() << " with " << other_node->index_() << std::endl;
        }
        if(other_node->index_ == node->index_) {
            return;
        }

        if (node->index_ < other_node->index_) {
            return;
        }
        const auto& wire = geo->getCellInfo(node->index_)->second.getWire();
        const auto& other_wire = geo->getCellInfo(other_node->index_)->second.getWire();
        
        double distance = geo->getMinDistanceBetweenSegments(wire.getFirstPoint(),
        wire.getSecondPoint(),
        other_wire.getFirstPoint(),
        other_wire.getSecondPoint());
        if(distance <10){
            if (node->index_ == 220006 || other_node->index_ == 220006) {
                std::cout << "adding " << node->index_() << " with " << other_node->index_() << std::endl;
            }
            node->cell_iterator_->second.addAdjacentCell(&(other_node->cell_iterator_->second));
            other_node->cell_iterator_->second.addAdjacentCell(&(node->cell_iterator_->second));
        }
        return;
    } 
    
    if (other_node->index_ == -1 && node->index_ == -1){
        searchAdjacentCells(node->left_,  other_node->left_, geo);
        searchAdjacentCells(node->left_,  other_node->right_, geo);
        searchAdjacentCells(node->right_, other_node->right_, geo);
        searchAdjacentCells(node->right_, other_node->left_, geo);
    } else if (node->index_ == -1) {
        searchAdjacentCells(node->left_,  other_node, geo);
        searchAdjacentCells(node->right_, other_node, geo);
    } else if (other_node->index_ == -1) {
        searchAdjacentCells(node, other_node->left_, geo);
        searchAdjacentCells(node, other_node->right_, geo);
    }
}
