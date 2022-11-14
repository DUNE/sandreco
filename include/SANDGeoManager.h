#include "SANDECALCellInfo.h"
#include "SANDSTTTubeInfo.h"

#include <map>

#ifndef SANDGEOMANAGER_H
#define SANDGEOMANAGER_H

class SANDGeoManager : public TObject {
    private:
        std::map<int, SANDECALCellInfo> cellmap_;
        std::map<int, SANDSTTTubeInfo> sttmap_;

    ClassDef(SANDGeoManager, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class SANDGeoManager + ;
#endif

#endif