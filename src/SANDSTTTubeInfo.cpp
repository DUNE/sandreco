/*Implementation of the SANDSTTTubeInfo class, 
for storing geometric info of the SAND STTs*/

#include "SANDSTTTubeInfo.h"

//Default constructor
SANDSTTTubeInfo::SANDSTTTubeInfo(){}

//Parametric constructor
SANDSTTTubeInfo::SANDSTTTubeInfo(int id, double x, double y, 
        double z, double length, Orient orientation, ReadoutEnd readout_end)
        : id_(id), x_(x), y_(y), z_(z),length_(length), 
        orientation_(orientation), readout_end_(readout_end){}

//Setter methods for the attributes
void SANDSTTTubeInfo::id(int arg_id){id_=arg_id;}
void SANDSTTTubeInfo::x(double arg_x){x_=arg_x;}
void SANDSTTTubeInfo::y(double arg_y){y_=arg_y;}
void SANDSTTTubeInfo::z(double arg_z){z_=arg_z;}
void SANDSTTTubeInfo::length(double arg_length){length_=arg_length;}
void SANDSTTTubeInfo::orientation(Orient arg_orientation){
        orientation_=arg_orientation;
        }
void SANDSTTTubeInfo::readout_end(ReadoutEnd arg_readout_end){
        readout_end_=arg_readout_end;
        }

//Getter methods for the attributes
int SANDSTTTubeInfo::id(){return id_;}
double SANDSTTTubeInfo::x(){return x_;}
double SANDSTTTubeInfo::y(){return y_;}
double SANDSTTTubeInfo::z(){return z_;}
double SANDSTTTubeInfo::length(){return length_;}
SANDSTTTubeInfo::Orient SANDSTTTubeInfo::orientation(){
    return orientation_;
    }
SANDSTTTubeInfo::ReadoutEnd SANDSTTTubeInfo::readout_end(){
    return readout_end_;
    }