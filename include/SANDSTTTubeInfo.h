#include <TObject.h>

#ifndef SANDSTTTUBEINFO_H
#define SANDSTTTUBEINFO_H

//class for storing the STT tubes geometrical info
class SANDSTTTubeInfo : public TObject{
    public:
        enum class Orient{kHorizontal, kVertical}; 
        enum class ReadoutEnd{kPlus, kMinus};
        int id_;                                      // id of tube
        double x_;                                    // x position of the center of the tube
        double y_;                                    // y position of the center of the tube
        double z_;                                    // z position of the center of the tube
        Orient orientation_;                          // orientation of the tube
        ReadoutEnd readout_end_;                      // end where signal are read
        double length_;                               // length of the tube

    ClassDef(SANDSTTTubeInfo, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class SANDSTTTubeInfo + ;
#endif


#endif 
