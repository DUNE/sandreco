#include "BVH.h"
#include <TVector2.h>



AABB::AABB(const sand_geometry::tracker::CellID& cellID, SANDGeoManager* geo){
    min_ = TVector3(1E9, 1E9, 1E9);
    max_ = TVector3(-1E9, -1E9, -1E9);
    // std::cout << __LINE__ << std::endl;
    
        std::vector<TVector3> vertices;
        auto cell = geo->getCellInfo(cellID)->second;
        auto p = cell.getWire().getFirstPoint();
        // p.Print();
        auto& plane = *geo->getPlaneInfo(cellID);
        auto p_rotated = geo->globalToRotated(TVector2(p.X(), p.Y()), plane);  
        // p_rotated.Print();
        auto h_2 = cell.getSize().h /2.;
        auto w_2 = cell.getSize().w /2.;
// std::cout << __LINE__ << std::endl;
        TVector2 v1(p_rotated.X(), p_rotated.Y() + h_2);
        TVector2 v2(p_rotated.X(), p_rotated.Y() - h_2);
        
        auto p1 = geo->rotatedToGlobal(TVector2(v1.X(), v1.Y()), plane); 
        // p1.Print();
        // break;
        auto p2 = geo->rotatedToGlobal(TVector2(v2.X(), v2.Y()), plane);
// std::cout << __LINE__ << std::endl;
        vertices.push_back(TVector3(p1.X(), p1.Y(), p.Z() + w_2));
        vertices.push_back(TVector3(p2.X(), p2.Y(), p.Z() + w_2));   

        p = cell.getWire().getSecondPoint();
        p_rotated = geo->globalToRotated(TVector2(p.X(), p.Y()), plane); 
        TVector2 v3(p_rotated.X(), p_rotated.Y() + h_2);
        TVector2 v4(p_rotated.X(), p_rotated.Y() - h_2);
// std::cout << __LINE__ << std::endl;
        auto p3 = geo->rotatedToGlobal(TVector2(v3.X(), v3.Y()), plane); 
        auto p4 = geo->rotatedToGlobal(TVector2(v4.X(), v4.Y()), plane);

        vertices.push_back(TVector3(p3.X(), p3.Y(), p.Z() - w_2));
        vertices.push_back(TVector3(p4.X(), p4.Y(), p.Z() - w_2));
// std::cout << __LINE__ << std::endl;
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

    
// std::cout << __LINE__ << std::endl;
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
    // node->aabb_.min_.Print();
    // node->aabb_.max_.Print();

    // node->indices_ = cells;

    // std::cout << __LINE__ << std::endl;
    if(std::distance(begin, end) == 1){
        node->index_ = *begin;
        std::cout << "ADDED CELL;, TOT: " << c << " " << (*begin)() << " ";
        std::cout << std::endl;
        c++;
        return;
    }
    // std::cout << __LINE__ << std::endl;
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
    // std::cout << __LINE__ << std::endl;

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
        std::nth_element(begin, begin + middle_point, end, sorting_function);
    
    // std::cout << __LINE__ << std::endl;
    
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

