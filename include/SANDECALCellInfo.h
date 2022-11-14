#include <TObject.h>

#ifndef SANDECALCELLINFO_H
#define SANDECALCELLINFO_H

//class for storing 
class SANDECALCellInfo : public TObject{
    public:
        enum class Orient{kHorizontal, kVertical};
        int id_;                                       // id of the cell
        double x_;                                     // x position of the center of the cell
        double y_;                                     // y position of the center of the cell
        double z_;                                     // z position of the center of the cell
        Orient orientation_;                           // orientation of the cell
        double length_;                                // length of the cell

    ClassDef(SANDECALCellInfo, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class SANDECALCellInfo + ;
#endif

#endif 
