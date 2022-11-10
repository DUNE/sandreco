#ifndef SANDGEOMANAGER_H
#define	SANDGEOMANAGER_H

#include "TObject.h"
#include "TClass.h"
class SANDGeoManager : public TObject {
    ClassDef(SANDGeoManager, 1)
};


#ifdef __MAKECINT__
#pragma link C++ class SANDGeoManager + ;
#endif

#endif
