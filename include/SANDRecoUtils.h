#ifndef SANDRECOUTILS_H
#define SANDRECOUTILS_H

#include <cmath>
// #include "utils.h"
#include "struct.h"
#include "SANDSTTTubeInfo.h"

#include "TVector3.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

class Helix
{
    public:
        Helix(double arg_R, double arg_dip, double arg_Phi0, int arg_h, TVector3 arg_x0) : 
                    R_(arg_R), dip_(arg_dip), Phi0_(arg_Phi0), h_(arg_h), x0_(arg_x0) {}
        
        // R : radius of the helix
        // pitch : dip angle 
        // Phi0 : is the azimuth angle of the starting point with respect to the helix axis
        // h : +/- sense of rotation
        // x0 : track starting point
        
        double x_h(double s) const{
            return x0_.X() + s * sin(dip_);
        }

        double y_h(double s) const{
            return x0_.Y() + R_ * (sin(Phi0_ + h_*s*cos(dip_)/R_) - sin(Phi0_));
        }

        double z_h(double s) const{
            return x0_.Z() + R_ * (cos(Phi0_ + h_*s*cos(dip_)/R_) - cos(Phi0_));
        }

        TVector3 GetPointAt(double s) const {
            return {x_h(s), y_h(s), z_h(s)};
        }

        void SetHelixParam(double* p){
            R_    = p[0];
            dip_  = p[1];
            Phi0_ = p[2];
            h_    = p[3];
            x0_   = {p[4],p[5],p[6]};
        }

        double R() const {return R_;};
        double dip() const {return dip_;};
        double Phi0() const {return Phi0_;};
        int h() const {return h_;};
        TVector3 x0() const {return x0_;};

    private:
        double R_;
        double dip_;
        double Phi0_;
        int h_;
        TVector3 x0_;
    
    ClassDef(Helix, 1);
};

class Line : public SANDWireInfo
{
    public:
        Line(double arg_dx, double arg_dy, double arg_ax, double arg_ay, double arg_z0) : 
                    dx_(arg_dx), dy_(arg_dy), ax_(arg_ax), ay_(arg_ay), z0_(arg_z0) {}

        // (dx,dy,dz) = vector parallel to wire
        // x = dx * t + ax
        // y = dy * t + ay
        // z = z0
        // a,b : line intercept

        void SetWireLineParam(double* p){
            dx_ = p[0];
            dy_ = p[1];
            ax_ = p[2];
            ay_ = p[3];
            z0_ = p[4];
        }
        double x_l(double t) const {
            return dx_ * t + ax_;
        }
        double y_l(double t) const {
            return dy_ * t + ay_;
        }
        double z_l() const {
            return z0_;
        }

        double dx() const {return dx_;};
        double dy() const {return dy_;};
        double ax() const {return ax_;};
        double ay() const {return ay_;};
        double z0() const {return z0_;};

        TVector3 GetPointAt(double t) const {
            return {x_l(t), y_l(t), z_l()};
        }

    private:
        double dx_;
        double dy_;
        double ax_;
        double ay_;
        double z0_;
    
    ClassDef(Line, 1);
};

namespace RecoUtils{ //RecoUtils

double GetImpactParameter(const Helix& helix, const Line& line, double s, double t);

double FunctorImpactParameter(const double* p);

double GetMinImpactParameter(const Helix& helix, const Line& line);

double GetExpectedRadiusFromDigit(const dg_tube& digit);

Line GetLineFromDigit(const dg_tube& digit);

std::vector<dg_tube> GetEventDigits();

// negative log likelihood

double NLL(const Helix& h, const std::vector<dg_tube>& digits);

double FunctorNLL(const double* p);

const double* GetHelixParameters();

} // RecoUtils

#ifdef __MAKECINT__
#pragma link C++ class Helix + ;
#pragma link C++ class Line + ;
#endif

// #ifdef __MAKECINT__
// #pragma link C++ class Line + ;
// #endif

#endif