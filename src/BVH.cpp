#include "BVH.h"
#include <TVector2.h>



AABB::AABB(const std::vector<sand_geometry::tracker::CellID>& cells, SANDGeoManager* geo){
    min_ = TVector3(1E9, 1E9, 1E9);
    max_ = TVector3(-1E9, -1E9, -1E9);
    // std::cout << __LINE__ << std::endl;
    for(const auto& cellID : cells){
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

    }
// std::cout << __LINE__ << std::endl;
}



void BVH::createTree(const std::vector<sand_geometry::tracker::CellID>& cells, SANDGeoManager* geo){
    root_.aabb = AABB(cells, geo);
    root_.aabb.min_.Print();
    root_.aabb.max_.Print();

}

