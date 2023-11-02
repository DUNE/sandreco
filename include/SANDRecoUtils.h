#ifndef SANDRECOUTILS_H
#define SANDRECOUTILS_H

#include <cmath>
// #include "utils.h"
#include "TVector3.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

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

        TVector3 GetPointAt(double s) const {
            return {x_h(s), y_h(s), z_h(s)};
        }

        void SetHelixParam(double* p){
            R = p[0];
            dip = p[1];
            Phi0 = p[2];
            h = p[3];
            x0 = {p[4],p[5],p[6]};
        }

    private:
        double R;
        double dip;
        double Phi0;
        int h;
        TVector3 x0;
    
    ClassDef(Helix, 1);
};

class Line //: public SANDWireInfo
{
    public:
        Line(double dx, double dy, double ax, double ay) : dx(dx), dy(dy), ax(ax), ay(ay) {}

        // (dx,dy,dz) = vector parallel to wire
        // x = dx * t + ax
        // y = dy * t + ay
        // a,b : line intercept

        void SetWireLineParam(double* p){
            dx = p[0];
            dy = p[1];
            ax = p[2];
            ay = p[3];
        }
        double x_l(double t) const {
            return dx * t + ax;
        }
        double y_l(double t) const {
            return dy * t + ay;
        }
        double z_l() const {
            return 0;
        }
        // double x_l(double t){
        //     if(orientation_ == kHorizontal){
        //         return x_ + t * dx;
        //     }else{
        //         return x_;
        //     }
        // }

        // double y_l(double t){
        //     if (orientation_ == kHorizontal)
        //     {
        //         return y_;
        //     }else{
        //         return y_ + t * dy; 
        //     }
            
        // }

        // double z_l(){
        //     return z_;
        // }

        TVector3 GetPointAt(double t) const {
            return {x_l(t), y_l(t), z_l()};
        }

    private:
        double dx;
        double dy;
        double ax;
        double ay;
    
    ClassDef(Line, 1);
};

// namespace RecoUtils{ //RecoUtils

double GetImpactParameter(const Helix& helix, const Line& line, double s, double t) {

    TVector3 helix_point = helix.GetPointAt(s);
    TVector3 line_point = line.GetPointAt(t);

    return (helix_point - line_point).Mag();
};

double FunctionImpactParameter(const double* p) {
    // first 7 entries of p to define the helix: R, dip, Phi0, h, x0
    // last 4 entries of p to define the line: dx, dy, ax, ay
    TVector3 x0{p[4], p[5], p[6]};

    Helix helix(p[0], p[1], p[2], p[3], x0);
    Line line(p[7],p[8],p[9],p[10]);

    double s = p[11];
    double t = p[12];

    return GetImpactParameter(helix, line, s, t);
}

double GetMinImpactParameter(const Helix& helix, const Line& line){

    ROOT::Math::Functor functor(&FunctionImpactParameter, 13);

    // gDebug=1;

    ROOT::Math::Minimizer * minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");

    minimizer->SetFunction(functor);

    // set boundaries for Phi0
    
    // helix params
    minimizer->SetFixedVariable(0, "R", 1);
    minimizer->SetFixedVariable(1, "dip", 0);
    minimizer->SetFixedVariable(2, "Phi0", 0);
    minimizer->SetFixedVariable(3, "h", -1);
    minimizer->SetFixedVariable(4, "x0_x", 1);
    minimizer->SetFixedVariable(5, "x0_y", 1);
    minimizer->SetFixedVariable(6, "x0_z", 1);
    // line params    
    minimizer->SetFixedVariable(7, "dx", 1);
    minimizer->SetFixedVariable(8, "dy", 1);
    minimizer->SetFixedVariable(9, "ax", 1);
    minimizer->SetFixedVariable(10, "ay", 1);
    // s, t
    minimizer->SetVariable(11, "s", 1, 0.001);
    minimizer->SetVariable(12, "t", 1, 0.001);

    if( minimizer->Minimize()) minimizer->PrintResults();
    return 1.;
}


int SANDRecoUtils(){
    
    TVector3 x0 = {0,0,0};

    // R, dip, Phi0, h, x0
    Helix h(1., 0., 0., 1, x0);
    
    // mx, my, ax, ay -> line that crosses (0,0,0)
    Line l(1., 1., 0., 0.);

    std::cout<< GetMinImpactParameter(h,l) <<"\n";

    return 0;
}

#ifdef __MAKECINT__
#pragma link C++ class Helix + ;
#pragma link C++ class Line + ;
#pragma link C++ class ImpactParameterFunctor + ;
#endif

// #ifdef __MAKECINT__
// #pragma link C++ class Line + ;
// #endif

#endif