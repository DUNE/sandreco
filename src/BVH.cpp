#include "BVH.h"
#include <TVector2.h>



AABB::AABB(const sand_geometry::tracker::CellID& cellID, SANDGeoManager* geo){
    min_ = TVector3(1E9, 1E9, 1E9);
    max_ = TVector3(-1E9, -1E9, -1E9);
      std::vector<TVector3> vertices;
        auto cell = geo->getCellInfo(cellID)->second;
        auto p = cell.getWire().getFirstPoint();
     
        auto& plane = *geo->getPlaneInfo(cellID);
        auto p_rotated = geo->globalToRotated(TVector2(p.X(), p.Y()), plane);  
     
        auto h_2 = cell.getSize().h /2.;
        auto w_2 = cell.getSize().w /2.;

        TVector2 v1(p_rotated.X(), p_rotated.Y() + h_2);
        TVector2 v2(p_rotated.X(), p_rotated.Y() - h_2);
        
        auto p1 = geo->rotatedToGlobal(TVector2(v1.X(), v1.Y()), plane); 

        auto p2 = geo->rotatedToGlobal(TVector2(v2.X(), v2.Y()), plane);

        vertices.push_back(TVector3(p1.X(), p1.Y(), p.Z() + w_2));
        vertices.push_back(TVector3(p2.X(), p2.Y(), p.Z() + w_2));   

        p = cell.getWire().getSecondPoint();
        p_rotated = geo->globalToRotated(TVector2(p.X(), p.Y()), plane); 
        TVector2 v3(p_rotated.X(), p_rotated.Y() + h_2);
        TVector2 v4(p_rotated.X(), p_rotated.Y() - h_2);

        auto p3 = geo->rotatedToGlobal(TVector2(v3.X(), v3.Y()), plane); 
        auto p4 = geo->rotatedToGlobal(TVector2(v4.X(), v4.Y()), plane);

        vertices.push_back(TVector3(p3.X(), p3.Y(), p.Z() - w_2));
        vertices.push_back(TVector3(p4.X(), p4.Y(), p.Z() - w_2));

        for(const auto &v : vertices){

        if(min_.X() > v.X()){
            min_.SetX(v.X());
        }
        if(max_.X() < v.X()){
            max_.SetX(v.X());
        }

        if(min_.Y() > v.Y()){
            min_.SetY(v.Y());
        }
        if(max_.Y() < v.Y()){
            max_.SetY(v.Y());
        }
       
        if(min_.Z() > v.Z()){
            min_.SetZ(v.Z());
        }
        if(max_.Z() < v.Z()){
            max_.SetZ(v.Z());
        }

        }

}

void AABB::expand(const AABB& second_aabb){
    if(min_.X() > second_aabb.min_.X()){
        min_.SetX(second_aabb.min_.X());
    }
    if(max_.X() < second_aabb.max_.X()){
        max_.SetX(second_aabb.max_.X());
    }

    if(min_.Y() > second_aabb.min_.Y()){
        min_.SetY(second_aabb.min_.Y());
    }
    if(max_.Y() < second_aabb.max_.Y()){
        max_.SetY(second_aabb.max_.Y());
    }
   
    if(min_.Z() > second_aabb.min_.Z()){
        min_.SetZ(second_aabb.min_.Z());
    }
    if(max_.Z() < second_aabb.max_.Z()){
        max_.SetZ(second_aabb.max_.Z());
    }
 }

bool AABB::isOverlapping(const AABB& second_aabb, double epsilon = 0){
    if(max_.X() + epsilon >= second_aabb.min_.X() && min_.X() - epsilon <= second_aabb.max_.X() && 
       max_.Y() + epsilon >= second_aabb.min_.Y() && min_.Y() - epsilon <= second_aabb.max_.Y() && 
       max_.Z() + epsilon >= second_aabb.min_.Z() && min_.Z() - epsilon <= second_aabb.max_.Z()){
        return true;
    }

    return false;
}

void BVH::fillCellAABBMap(std::vector<sand_geometry::tracker::CellID>& cells, SANDGeoManager* geo){
    for(const auto &cell : cells){
        cellAABBs_[cell] = AABB(cell, geo);
    }
}

int c = 0;

void BVH::createTree(std::unique_ptr<Node>& node, std::vector<sand_geometry::tracker::CellID>::iterator begin, std::vector<sand_geometry::tracker::CellID>::iterator end, SANDGeoManager* geo){

    node->aabb_ = cellAABBs_[*begin];
    for (auto it = begin + 1; it != end; ++it) {
        node->aabb_.expand(cellAABBs_[*it]);
    }

    // node->indices_ = cells;

    
    if(std::distance(begin, end) == 1){
        node->index_ = *begin;
        // std::cout << "ADDED CELL;, TOT: " << c << " " << (*begin)() << std::endl;
        c++;
        return;
    }
    
    int axis = 0;
    double deltaX = node->aabb_.max_.X() - node->aabb_.min_.X();
    double deltaY = node->aabb_.max_.Y() - node->aabb_.min_.Y();
    double deltaZ = node->aabb_.max_.Z() - node->aabb_.min_.Z();


    if(deltaY > deltaX){
        axis = 1;
    }
    if(deltaZ > deltaY){
        axis = 2;
    } 
    

    auto sorting_function = [axis, geo](const sand_geometry::tracker::CellID& c1, const sand_geometry::tracker::CellID& c2){
        auto center_1 = geo->getCellInfo(c1)->second.getWire().getCenter();
        auto center_2 = geo->getCellInfo(c2)->second.getWire().getCenter();
        
        if(axis==0){
            return center_1.X() < center_2.X();
        }
        if(axis==1){
            return center_1.Y() < center_2.Y();
        }
        if(axis==2){
            return center_1.Z() < center_2.Z();
        }
        return false;
    };

    // for (const auto& cell:cells) {
    //     std::cout << cell() << " ";
    // }
    // std::cout << std::endl;
    // for (const auto& cell:cells) {
        //     std::cout << cell() << " ";
        // }
        // std::cout << std::endl;
        
        
        int middle_point= std::distance(begin, end) / 2;
        // std::nth_element(begin, begin + middle_point, end, sorting_function);
        std::sort(begin, end, sorting_function);
    
    
    
    // std::vector<sand_geometry::tracker::CellID> left_cell_id(cells.begin(), cells.begin() + middle_point);
    node->left_ = std::make_unique<Node>();
    createTree(node->left_, begin, begin + middle_point, geo);
    // std::cout << __LINE__ << " " << size << std::endl;
    
    // std::vector<sand_geometry::tracker::CellID> right_cell_id(begin + middle_point, end); 
    // std::cout << __LINE__ << " " << size << std::endl;
    node->right_ = std::make_unique<Node>();
    createTree(node->right_, begin + middle_point, end, geo);
 
    return;
}

const std::map<sand_geometry::tracker::CellID,std::vector<sand_geometry::tracker::CellID>>& BVH::getAdjacentCells(SANDGeoManager* geo){
    searchAdjacentCells(root_, root_, geo);

    return cellID_to_adjacent_cells_;
}

void BVH::getAdjacentCells(std::unique_ptr<Node>& node, SANDGeoManager* geo){
    // if(node->index_ != -1) {
    //     std::vector<sand_geometry::tracker::CellID> adjacent_cells;
    //     searchAdjacentCells(node,root_, adjacent_cells, geo);
    //     // cellID_to_adjacent_cells_[node->index_] = adjacent_cells; 
    // } else {
    //     getAdjacentCells(node->left_, geo); 
    //     getAdjacentCells(node->right_, geo); 
    // }

    // return;
}

void BVH::searchAdjacentCells(std::unique_ptr<Node>& node, std::unique_ptr<Node>& other_node, SANDGeoManager* geo){
    if(!node || !other_node) return;
    // if(node->index_ == other_node->index_) {return;}
    
    if (node->index_ == 420000 || other_node->index_ == 420000 ) {
        std::cout << "Checking " << node->index_() << " and " << other_node->index_() << std::endl;
    }

    if(!node->aabb_.isOverlapping(other_node->aabb_, 1)) return;

    

    if(other_node->index_ != -1 && node->index_ != -1) {
        if(other_node->index_ == node->index_) {
            return;
        }
        if (node->index_ < other_node->index_) {
            return;
        }
        auto wire = geo->getCellInfo(node->index_)->second.getWire();
        auto other_wire = geo->getCellInfo(other_node->index_)->second.getWire();

        double distance = geo->getMinDistanceBetweenSegments(wire.getFirstPoint(),
                                                            wire.getSecondPoint(),
                                                            other_wire.getFirstPoint(),
                                                            other_wire.getSecondPoint());
        if(distance <10){
            cellID_to_adjacent_cells_[node->index_].push_back(other_node->index_);
            cellID_to_adjacent_cells_[other_node->index_].push_back(node->index_);
            // std::cout << "DDED " << other_node->index_() << " TO " << node->index_() << " and the other way around" << std::endl;
        }
        return;
    } 
    
    if (other_node->index_ == -1 && node->index_ == -1){
        searchAdjacentCells(node->left_,  other_node->left_, geo);
        searchAdjacentCells(node->left_,  other_node->right_, geo);
        searchAdjacentCells(node->right_, other_node->right_, geo);
        // if (other_node->index_ != node->index_) {
            searchAdjacentCells(node->right_, other_node->left_, geo);
        // }
    } else if (node->index_ == -1) {
        searchAdjacentCells(node->left_,  other_node, geo);
        searchAdjacentCells(node->right_, other_node, geo);
    } else if (other_node->index_ == -1) {
        searchAdjacentCells(other_node->left_,  node, geo);
        searchAdjacentCells(other_node->right_, node, geo);
    }
}
