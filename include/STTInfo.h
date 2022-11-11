#include <TObject.h>
#include <TVector3.h>

#include <map>

#ifndef STTINFO_H
#define STTINFO_H

//class for storing the STT tubes geometrical info
class STTInfo : public TObject{
    public:
        enum class Orient{khorizontal, kvertical};
        enum class ReadoutEnd{kplus, kminus};
        int fid_;
        double fx_;
        double fy_;
        double fz_;
        Orient forient_;
        ReadoutEnd freadoutend_;
        double flen_;

    ClassDef(STTInfo, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class STTInfo + ;
#endif


#endif 
