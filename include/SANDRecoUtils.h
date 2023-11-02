#ifndef SANDRECOUTILS_H
#define SANDRECOUTILS_H

#include <cmath>
#include "TVector3.h"

class Helix
{
    public:
        Helix(double R, double dip, double Phi0, double h, TVector3 x0) : R(R), dip(dip), Phi0(Phi0), h(h), x0(x0) {}
        
        // R : radius of the helix
        // pitch : dip angle 
        // Phi0 : is the azimuth angle of the starting point with respect to the helix axis
        // h : +/- sense of rotation
        // x0 : track starting point
        
        double x_h(double s) const{
            return x0.X() + s * sin(dip);
        }

        double y_h(double s) const{
            return x0.Y() + R * (sin(Phi0 + h*s*cos(dip)/R) - sin(Phi0));
        }

        double z_h(double s) const{
            return x0.Z() + R * (cos(Phi0 + h*s*cos(dip)/R) - cos(Phi0));
        }

        double GetImpactParameter(double s, double t, WireLine WL) const{
            // takes to input parameters (s,t) that define a point on the 
            // helix and a point on the wire line WL and return its distance
            double dx = x_h(s) - WL.x_l(t);
            double dy = y_h(s) - WL.y_l(t);
            double dz = z_h(s) - WL.z_l(t);
            return sqrt(dx*dx + dy*dy + dz*dz);
        }

    private:
        double R;
        double dip;
        double Phi0;
        int h;
        TVector3 x0;
    
    ClassDef(Helix, 1);
}

class WireLine : public SANDWireInfo
{
    public:
        WireLine(double dx, double dy) : dx(dx), dy(dy) {}

        // (dx,dy,dz) = vector parallel to wire

        double x_l(double t){
            if(orientation_ == kHorizontal){
                return x_ + t * dx;
            }else{
                return x_;
            }
        }

        double y_l(double t){
            if (orientation_ == kHorizontal)
            {
                return y_;
            }else{
                return y_ + t * dy; 
            }
            
        }

        double z_l(){
            return z_;
        }

    private:
        double dx;
        double dy;
    
    ClassDef(WireLine, 1);
}

#ifdef __MAKECINT__
#pragma link C++ class Helix + ;
#pragma link C++ class WireLine + ;
#endif

// #ifdef __MAKECINT__
// #pragma link C++ class WireLine + ;
// #endif

#endif