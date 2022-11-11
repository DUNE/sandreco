#include <TObject.h>
#include <TVector3.h>

#include <map>

#ifndef CELLINFO_H
#define CELLINFO_H

//class for storing 
class CellInfo : public TObject{
    public:
        enum class Orient{khorizontal, kvertical};
        int fid_;
        double fx_;
        double fy_;
        double fz_;
        Orient forient_;
        double flen_;

    ClassDef(CellInfo, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class CellInfo + ;
#endif

#endif 
