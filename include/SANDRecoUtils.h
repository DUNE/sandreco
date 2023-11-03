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

class Line //: public SANDWireInfo
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
    Line line(p[7],p[8],p[9],p[10],p[11]);

    double s = p[12];
    double t = p[13];

    return GetImpactParameter(helix, line, s, t);
}

double GetMinImpactParameter(const Helix& helix, const Line& line){

    ROOT::Math::Functor functor(&FunctionImpactParameter, 13);

    // gDebug=1;

    ROOT::Math::Minimizer * minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");

    minimizer->SetFunction(functor);

    // set boundaries for Phi0
    
    // helix params
    minimizer->SetFixedVariable(0, "R",     helix.R());
    minimizer->SetFixedVariable(1, "dip",   helix.dip());
    minimizer->SetFixedVariable(2, "Phi0",  helix.Phi0());
    minimizer->SetFixedVariable(3, "h",     helix.h());
    minimizer->SetFixedVariable(4, "x0_x",  helix.x0().X());
    minimizer->SetFixedVariable(5, "x0_y",  helix.x0().Y());
    minimizer->SetFixedVariable(6, "x0_z",  helix.x0().Z());
    // line params    
    minimizer->SetFixedVariable(7, "dx",    line.dx());
    minimizer->SetFixedVariable(8, "dy",    line.dy());
    minimizer->SetFixedVariable(9, "ax",    line.ax());
    minimizer->SetFixedVariable(10, "ay",   line.ay());
    minimizer->SetFixedVariable(11, "z0",   line.z0());
    // s, t
    minimizer->SetVariable(12, "s", 1, 0.0001);
    minimizer->SetVariable(13, "t", 1, 0.0001);

    bool HasMinimized = minimizer->Minimize();

    if(HasMinimized) minimizer->PrintResults();
    return HasMinimized;
}


int SANDRecoUtils(){

    double R = 2.3;
    double Phi0 = 1.4;
    int hel = 1;
    double dip = 0.;

    double s_star = -1* Phi0*R/hel/cos(dip); // so that Phi=Phi0
    
    TVector3 x0 = {33., 14., 6.}; //

    // R, dip, Phi0, h, x0
    Helix h(R, dip, Phi0, hel, x0); // this helix crosses (0,0,0) for s=0
    
    auto helix_at_s_star = h.GetPointAt(s_star);
    // mx, my, ax, ay -> line that crosses (0,0,0)
    Line l(2., 3., helix_at_s_star.X(), helix_at_s_star.Y(), helix_at_s_star.Z()); // this line crosses (0,0,0) for t=0

    auto line_at_t_star = l.GetPointAt(0.);

    std::cout<<"s_star : "<<s_star<<"\n";
    std::cout<<"helix at s=s_star hit the point at : "<<helix_at_s_star.X()<<" "<<helix_at_s_star.Y()<<" "<<helix_at_s_star.Z()<<"\n";
    std::cout<<"t_star : 0 \n";
    std::cout<<"line at t=0 hit the point at : "<<line_at_t_star.X()<<" "<<line_at_t_star.Y()<<" "<<line_at_t_star.Z()<<"\n";
    std::cout<<"GetImpactParameter(h, l , s_star, t_star) : "<<GetImpactParameter(h, l , s_star, 0)<<"\n";
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