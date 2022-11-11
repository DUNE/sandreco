#include "CellInfo.h"
#include "STTInfo.h"

#include <TObject.h>
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TPRegexp.h>
#include <TString.h>
#include <TVector2.h>

#include "struct.h"

#include <TG4Event.h>

#ifndef SANDGEOMANAGER_H
#define SANDGEOMANAGER_H

class SANDGeoManager : public TObject {
    private:
        std::map<int, CellInfo> cellmap_;
        std::map<int, STTInfo> sttmap_;

    ClassDef(SANDGeoManager, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class SANDGeoManager + ;
#endif

#endif