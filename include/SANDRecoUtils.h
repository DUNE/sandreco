#ifndef SANDRECOUTILS_H
#define SANDRECOUTILS_H

#include <cmath>
#include <complex>
#include <TRandom3.h>
#include <TStyle.h>
#include <type_traits>

// #include "utils.h"
#include "struct.h"
#include "utils.h"
#include "SANDWireInfo.h"
#include "SANDGeoManager.h"

#include "TVector3.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TDatabasePDG.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFitResult.h"

extern TGeoManager* geo;

extern SANDGeoManager geo_manager;

class Helix
{
    public:
        // R : radius of the helix
        // dip : dip angle 
        // Phi0 : is the azimuth angle of the starting point with respect to the helix axis
        // h : +/- sense of rotation
        // x0 : track starting pofint

        // optionsl : low_lim, up_lim 
        // lower and upper limit for linear parameter s to be used in case we want restrict the helix to a portion

        Helix(double arg_R, double arg_dip, double arg_Phi0, int arg_h, TVector3 arg_x0, double low_lim = -999., double up_lim = -999.) : 
                    R_(arg_R), dip_(arg_dip), Phi0_(arg_Phi0), h_(arg_h), x0_(arg_x0), low_lim_(low_lim), up_lim_(up_lim) {}

        Helix(const TG4Trajectory& trj) {
            auto momentum = trj.GetInitialMomentum();
            auto pt       = sqrt(momentum.Z()*momentum.Z() + momentum.Y()*momentum.Y());
            auto pl       = momentum.X();
            int pdg       = trj.GetPDGCode();
            TDatabasePDG    database; // maybe class static member could be better
            int charge    = database.GetParticle(pdg)->Charge();

            x0_           = trj.Points[0].GetPosition().Vect();
            dip_          = TMath::ATan2(pl,pt);
            R_            = (pt/(0.3*0.6)); // [m] = [GeV]/[T] or [mm] = [MeV]/[T] 
            h_            = (charge < 0) ? 1 : -1;
            Phi0_         = TMath::ATan2(momentum.Y(),momentum.Z()) + h_*TMath::Pi()*0.5; // for mu pi/2 should be added

            up_lim_       = 0.; // track starts from s = low_lim to s = 0 (for s<0 tracks proceed upwards according to this parmetrization)
            low_lim_      = -1e3;
        }

        Helix() {
            R_    = 0.;
            dip_  = 0.;
            Phi0_ = 0.;
            h_    = 1;
            x0_   = {0.,0.,0.};
        }

        // copy constructor
        Helix(const Helix& source){
            R_    = source.R();
            dip_  = source.dip();
            Phi0_ = source.Phi0();
            h_    = source.h();
            x0_   = source.x0();
        }

        double x_h(double s) const{
            return x0_.X() - s * sin(dip_);
        }

        double y_h(double s) const{
            return x0_.Y() + R_ * (sin(Phi0_ + h_*s*cos(dip_)/R_) - sin(Phi0_));
        }

        double z_h(double s) const{
            return x0_.Z() + R_ * (cos(Phi0_ + h_*s*cos(dip_)/R_) - cos(Phi0_));
        }

        // first derivatives
        double dx_over_ds() const{
            return -sin(dip_);
        }

        double dy_over_ds(double s) const{
            return cos(Phi0_ + h_*s*cos(dip_)/R_)*h_*cos(dip_);
        }

        double dz_over_ds(double s) const{
            return -sin(Phi0_ + h_*s*cos(dip_)/R_)*h_*cos(dip_);
        }

        TVector3 GetPointAt(double s) const {
            return {x_h(s), y_h(s), z_h(s)};
        }

        TVector3 GetTangentVector(double s) const {
            return {dx_over_ds(), dy_over_ds(s), dz_over_ds(s)};
        }

        double GetPhiFromZ(double z) const {
            return TMath::ACos((z - x0_.Z())/R_ + cos(Phi0_)) - Phi0_;
        }

        double GetSFromPhi(double Phi) const {
            return Phi * R_ / h_ / cos(dip_);
        }

        std::vector<TVector3> GetHelixPoints(double s_min = -1e3, double s_max = 1e3, double step = 10.) const {
            
            if(s_min>=s_max){
                std::cout << "s_min cannot be larger than s_max";
                throw;
            }

            std::vector<TVector3> points;

            while (s_min<s_max)
            {
                points.push_back(GetPointAt(s_min));
                s_min += step;
            }
            
            return points;
        } 

        void SetHelixParam(double* p){
            R_    = p[0];
            dip_  = p[1];
            Phi0_ = p[2];
            h_    = p[3];
            x0_   = {p[4],p[5],p[6]};
        }

        void SetR(double arg){
            R_ = arg;
        }

        void Setdip(double arg){
            dip_ = arg;
        }

        void Seth(double arg){
            h_ = arg;
        }

        void SetPhi0(double arg){
            R_ = arg;
        }

        void Setx0(TVector3 arg){
            x0_ = arg;
        }

        void SetLowLim(double arg_low){
            low_lim_ = arg_low;
        }
        
        void SetUpLim(double arg_up){
            up_lim_ = arg_up;
        }

        void SetHelixRangeFromDigit(const dg_wire& digit){
        // gives the s_min and s_max values that define the portion of the
        // helix in the digit z range

        // mm -> module thickenss is 10 mm, set 8 mm x 2 to be conservative
        double z_min   = digit.z - 8;
        double z_max   = digit.z + 8;

        double Phi_min = GetPhiFromZ(z_max); // Phi is decreasing function of z (goes as acos(z))
        double Phi_max = GetPhiFromZ(z_min);

        SetLowLim(GetSFromPhi(Phi_min));
        SetUpLim(GetSFromPhi(Phi_max));
        }

        void PrintHelixPars() const {
            std::cout<<std::setprecision(8)<<"R   -> "<<R_<<"\n";
            std::cout<<std::setprecision(8)<<"dip -> "<<dip_<<"\n";
            std::cout<<std::setprecision(8)<<"phi -> "<<Phi0_<<"\n";
            std::cout<<std::setprecision(8)<<"h   -> "<<h_<<"\n";
            std::cout<<std::setprecision(8)<<"x0, y0, z0 -> "<< x0_.X() << " " 
                                                             << x0_.Y() << " " 
                                                             << x0_.Z() <<"\n";
            std::cout<<std::setprecision(8)<<"zc, yc     -> "<< x0_.Z() - R_*cos(Phi0_) <<" "
                                                             << x0_.Y() - R_*sin(Phi0_) <<"\n";
            std::cout << "\n";                                                             
        }

        double   R()      const {return R_;};
        double   dip()    const {return dip_;};
        double   Phi0()   const {return Phi0_;};
        int      h()      const {return h_;};
        TVector3 x0()     const {return x0_;};
        double   LowLim() const {return low_lim_;};
        double   UpLim()  const {return up_lim_;};
        TVector2 Center() const {return {x0_.Z()-R_*cos(Phi0_), x0_.Y()-R_*sin(Phi0_)};};

        virtual ~Helix() {}

    private:
        double   R_;
        double   dip_;
        double   Phi0_;
        int      h_;
        TVector3 x0_;
        double   low_lim_;
        double   up_lim_;
    
    ClassDef(Helix, 1);
};

class Circle
{
    public:
        Circle(double arg_center_x, double arg_center_y, double arg_R) :
            center_x_(arg_center_x), center_y_(arg_center_y), R_(arg_R) {}

        Circle() :
            center_x_(0), center_y_(0), R_(1) {}
        
        double x_l(double angle) const {
            // angle should be in [0, 2pi)
            return center_x_ + R_ * cos(angle);
        }
        double y_l(double angle) const {
            // angle should be in [0, 2pi)
            return center_y_ + R_ * sin(angle);
        }

        double dx_derivative(double angle) const {
            return -1. * R_ * sin(angle);
        }
        
        double dy_derivative(double angle) const {
            return R_ * cos(angle);
        }

        TVector2 GetPointAt(double angle) const {
            return {x_l(angle), y_l(angle)};
        }

        TVector2 GetDerivativeAt(double angle) const {
            return {dx_derivative(angle), dy_derivative(angle)};
        }

        TVector2 GetDerivativeAt(double arg_x, double arg_y) const {
            double angle = GetAngleFromPoint(arg_x, arg_y);
            return {dx_derivative(angle), dy_derivative(angle)};
        }

        double GetAngleFromPoint(double arg_x, double arg_y) const {
            // return the angle [0, 2pi) corresponding to the point (x,y)
            if(arg_y -  center_y_ > 0){
                return atan2(arg_y - center_y_, arg_x - center_x_);
            }else{
                return atan2(arg_y - center_y_, arg_x - center_x_) + 2 * TMath::Pi();
            }
        }

        TF1* GetUpperSemiCircle() const {
            // y as function of the x coordinate
            TF1* up_circle = new TF1("UpperCircle", "[0] + sqrt([1]*[1] - (x-[2])*(x-[2]))");
            up_circle->SetParameters(center_y_, R_, center_x_);
            return up_circle;
        }

        TF1* GetLowerSemiCircle() const {
            // y as function of the x coordinate
            TF1* up_circle = new TF1("LowerCircle", "[0] - sqrt([1]*[1] - (x-[2])*(x-[2]))");
            up_circle->SetParameters(center_y_, R_, center_x_);
            return up_circle;
        }

        double Distance2Point(TVector2 point) const {
            return fabs((center() - point).Mod() - R());
        }

        TVector2 center() const {
            TVector2 c = {center_x_, center_y_};
            return c;
        }

        void PrintCircleInfo() const {
            std::cout << "circle center (x,y) : ("
                      << center_x_ << ", "
                      << center_y_ << "), R : "
                      << R_ << "\n";
        }

        double center_x() const {return center_x_;};
        double center_y() const {return center_y_;};
        double R() const {return R_;};

        virtual ~Circle() {}

    private:
        double center_x_;
        double center_y_;
        double R_;
    
    ClassDef(Circle, 1);
};

class Spiral2D
{
    /*
        (arg_center_x, arg_center_y) : spiral center
        arg_R0 : Spiral starting radius
        theta0 : angle offset
        k : rate of radius decrease
    */
    public:
        Spiral2D(double arg_center_x, double arg_center_y, double arg_R0, double arg_k) : 
        center_x_(arg_center_x), center_y_(arg_center_y), R0_(arg_R0), k_(arg_k) {}

        Spiral2D(){
            center_x_ = 0.;
            center_y_ = 0.;
            R0_ = 1.;
            k_ = 0.1;
        }

        double x_l(double angle) const {
            // angle should be in [0, 2pi)
            return center_x_ + (R0_ - k_ * angle) * cos(angle);
        }

        double y_l(double angle) const {
            // angle should be in [0, inf)
            return center_y_ + (R0_ - k_ * angle) * sin(angle);
        }

        TVector2 center() const {
            TVector2 c = {center_x_, center_y_};
            return c;
        }

        double center_x() const {return center_x_;};
        double center_y() const {return center_y_;};
        double R0() const {return R0_;};
        double kequazione() const {return k_;};

    private:
        double center_x_;
        double center_y_;
        double R0_;
        double k_;
};

class Line2D
{
    public:
    /*
        parametric equation : 
            x = dx * t + ax
            y = dy * t + ay
        cartesian equation : 
            y = m * x + q
    */ 
        Line2D(double arg_dx, double arg_dy, double arg_ax, double arg_ay){
            dx_ = arg_dx;
            dy_ = arg_dy;
            ax_ = arg_ax;
            ay_ = arg_ay;
            m_ = dy_ / dx_;
            q_ = arg_ay - m_ * arg_ax;
            direction_ = {dx_, dy_};
            p0_ = {ax_, ay_};
        }

        Line2D(double arg_m, double arg_q, double arg_ax = 0.){
            m_ = arg_m;
            q_ = arg_q;
            dx_ = 1.;
            dy_ = arg_m;
            ax_ = arg_ax;
            ay_ = arg_q + arg_m * arg_ax;
            direction_ = {dx_, dy_};
            p0_ = {ax_, ay_};
        }

        Line2D(TVector2 arg_direction, TVector2 arg_p0){
            direction_ = arg_direction;
            p0_ = arg_p0;
            dx_ = arg_direction.X();
            dy_ = arg_direction.Y();
            ax_ = arg_p0.X();
            ay_ = arg_p0.Y();
            m_ = dy_ / dx_;
            q_ = arg_p0.Y() - m_ * arg_p0.X();
        }

        Line2D(){
            dx_ = 1.;
            dy_ = 1.;
            ax_ = 0.;
            ay_ = 0.;
            m_ = 1.;
            q_ = 0.;
            direction_ = {1.,1.};
            p0_ = {0.,0.};
        }

        double Distance2Point(TVector2 p1) const {
            TVector2 diff = p1 - p0_;

            double projection_length = diff * direction_.Unit();
            TVector2 projection = projection_length * direction_.Unit();

            TVector2 diff_projection = diff - projection;
            return diff_projection.Mod();
        }

        double GetXFromY(double arg_y) const {
            return (arg_y - q_)/m_;
        }

        double dx() const {return dx_;};
        double dy() const {return dy_;};
        double ax() const {return ax_;};
        double ay() const {return ay_;};
        double m() const {return m_;};
        double q() const {return q_;};
        TVector2 direction() const {return direction_;};
        TVector2 p0() const {return p0_;};

        virtual ~Line2D() {}

    private:
        double dx_;
        double dy_;
        double ax_;
        double ay_;
        double m_;
        double q_;
        TVector2 direction_;
        TVector2 p0_;
    
    ClassDef(Line2D, 1);
};

class Line : public SANDWireInfo
{
    public:
        // (dx,dy,dz) = vector parallel to wire
        // x = dx * t + ax
        // y = dy * t + ay
        // z = dz * t + az
        // ax,ay,az : line intercepts

        // optional : low_lim_, up_lim_
        // usefull as we are interest in a portion of the line defined by the wire length
        
        Line(double arg_dx, double arg_dy, double arg_dz, double arg_ax, double arg_ay, double arg_az, double low_lim = -1e6, double up_lim = 1e6) : 
                    dx_(arg_dx), dy_(arg_dy), dz_(arg_dz), ax_(arg_ax), ay_(arg_ay), az_(arg_az), low_lim_(low_lim), up_lim_(up_lim) {}

        Line(const TG4HitSegment& hit){
            // define a line from hit that lies on it
            TVector3 hit_direction = (hit.GetStop() - hit.GetStart()).Vect();
            TVector3 middle_point = ((hit.GetStop() + hit.GetStart())*0.5).Vect();
            // double hit_length = (hit.GetStop() - hit.GetStart()).Vect().Mag();

            dx_ = hit_direction.X();
            dy_ = hit_direction.Y();
            dz_ = hit_direction.Z();

            // X0 of the line is the hit middle point
            ax_ = middle_point.X();
            ay_ = middle_point.Y();
            az_ = middle_point.Z();

            up_lim_ = 0.5;
            low_lim_ = - 0.5; 
        }

        Line(){
            dx_ = 0.;
            dy_ = 0.;
            dz_ = 0.;

            ax_ = 0.;
            ay_ = 0.;
            az_ = 0.;
        }

        void SetWireLineParam(double* p){
            dx_ = p[0];
            dy_ = p[1];
            dz_ = p[2];
            ax_ = p[3];
            ay_ = p[4];
            az_ = p[5];
        }

        void SetLowLim(double arg_low){
            low_lim_ = arg_low;
        }
        
        void SetUpLim(double arg_up){
            up_lim_ = arg_up;
        }

        void SetLineLength(double l){
            SetLowLim(-l/2./GetDirectionVector().Mag());
            SetUpLim(l/2./GetDirectionVector().Mag());
        }
 
        void SetLineRangeFromDigit(const dg_wire& digit){
            SetLowLim(-digit.wire_length/2./GetDirectionVector().Mag());
            SetUpLim(digit.wire_length/2./ GetDirectionVector().Mag());
        }

        double x_l(double t) const {
            return dx_ * t + ax_;
        }
        double y_l(double t) const {
            return dy_ * t + ay_;
        }
        double z_l(double t) const {
            return dz_ * t + az_;
        }

        // first derivatives
        double dx_over_dt() const {
            return dx_;
        }
        double dy_over_dt() const {
            return dy_;
        }
        double dz_over_dt() const {
            return dz_;
        }
        // derivatives

        TVector3 GetDirectionVector() const{
            // return the wire direction vector
            return {dx_over_dt(), dy_over_dt(), dz_over_dt()};
        }

        TVector3 GetPointAt(double t) const {
            return {x_l(t), y_l(t), z_l(t)};
        }

        TVector3 GetLinePointX0() const{
            return GetPointAt(0);
        }

        TVector3 GetLineUpperLimit() const{
            return GetPointAt(up_lim_);
        }

        TVector3 GetLineLowerLimit() const{
            return GetPointAt(low_lim_);
        }

        double GetLineLength() const{
            return (GetPointAt(up_lim_) - GetPointAt(low_lim_)).Mag();
        }

        bool IsInLineLimits(double t) const{
            /*
            find if t is a point on the line
            between line [low_lim, up_lim]
            x----------------x----------------x
            low_lim          X0           up_lim
            */
           auto point = GetPointAt(t);
           auto line_center = GetLinePointX0();
           auto pointDist2Center = (point - line_center).Mag();
           return (pointDist2Center < GetLineLength()/2.);
        }

        void PrintLineInfos(){
            std::cout << "line X0 : (" << ax() <<
                                   ", "<< ay() <<
                                   ", "<< az() <<")\n";
            std::cout << "line direction : (" << dx() <<
                                   ", "<< dy() <<
                                   ", "<< dz() <<")\n";
            std::cout << "line limits t : (" << LowLim() 
                                      << ", "<< UpLim() <<")\n";                                 
        }

        double dx() const {return dx_;};
        double dy() const {return dy_;};
        double dz() const {return dz_;};
        double ax() const {return ax_;};
        double ay() const {return ay_;};
        double az() const {return az_;};
        double LowLim() const {return low_lim_;};
        double UpLim() const {return up_lim_;};

    virtual ~Line() {}

    private:
        double dx_;
        double dy_;
        double dz_;
        double ax_;
        double ay_;
        double az_;
        double low_lim_;
        double up_lim_;
    
    ClassDef(Line, 1);
};

namespace RecoUtils{ //RecoUtils

extern              std::vector<dg_wire>* event_digits;

void                InitWireInfos(TGeoManager* g);

double              GetDistHelix2Line(const Helix& helix, double s, const Line& line, double& t);

double              GetDistHelix2LineDerivative(const Helix& helix, double s, const Line& line); // may be can be remover

double              GetLinePointDistance(const Line& l, TVector3 point);

double              GetLineLineDistance(const Line& l1, const Line& l2, double& closest2line2, double& closest2line1);

double              GetSegmentSegmentDistance(const Line& l1, const Line& l2, double& closest2line2, double& closest2line1);

double              GetImpactParameter(const Helix& helix, const Line& line, double s, double t);

double              FunctorImpactParameter(const double* p);

double              GetMinImpactParameter(const Helix& helix, const Line& line);

double              GetMinImpactParameter(const Helix& helix, const Line& line, double& s_min, double& t_min, bool& HasMinimized);

double              NewtonRaphson2D(TF1* f, TF1* fprime, double& x_guess, double tol, unsigned int max_iterations);

Line                GetLineFromWire(const dg_wire& digit);

Line2D              WiresLinearFit(const std::vector<dg_wire*>& wires);

TF1*                WiresSinFit(const std::vector<dg_wire>& wires);

Circle              WiresCircleFit(const std::vector<dg_wire*>& wires);

TF1*                InvertSin(TF1* fSin);

double              NLL(Helix& h,const std::vector<dg_wire>& digits); // negative log likelihood function

double              FunctorNLL(const double* p);

double              GetDipAngleFromCircleLine(const Circle& circle, 
                                             const Line2D& line, 
                                             double Phi0,
                                             int helicity,
                                             TVector3& Momentum);

Helix               GetHelixFromCircleLine(const Circle& circle, 
                                           const Line2D& line, 
                                           const Helix& true_helix,
                                           TVector3& momentum);

const double*       GetHelixParameters(const Helix& helix_initial_guess, int& TMinuitStatus);

dg_wire             Copy(const dg_wire& wire);

std::vector<double> SmearVariable(double mean, double sigma, int nof_points);

std::vector<Line2D> GetTangentsTo2Circles(const Circle& c1, const Circle& c2);

Line2D              GetTangetTo3Circles(const Circle& c1, const Circle& c2, const Circle& c3); // bug to be fixed

Line2D              GetBestTangent2NCircles(const std::vector<Circle>& circles);

} // RecoUtils

struct Parameter
{
    std::string name;
    int id;
    bool fixed_in_fit;
    double initial_guess;
    double value;
    double error;
    // Parameter() 
    //     : name(""), fixed_in_fit(false), initial_guess(0.0), value(0.0), error(0.0) {}
    // Parameter(std::string n, int i, bool f, double ig, double v, double e) 
    //     : name(n), id(i), fixed_in_fit(f), initial_guess(ig), value(v), error(e) {}
};

struct MinuitFitInfos
{
    // std::string                 Auxiliary_name;
    const char*                 Auxiliary_name;
    int                         TMinuitFinalStatus; // 0 if converged, 4 if falied
    int                         NIterations; // number of iterations to reach the minimum
    double                      MinValue; // value of the func to minimize at its minimum
    std::vector<Parameter>      fitted_parameters;
};

struct RecoObject
{
    // edepsim info of the reconstructed track_____________________
    std::vector<TLorentzVector> trj_points;
    
    // digitization info___________________________________________
    /*
      fired_wires : 
        all the fired wires related to the reconstructed obj.
        These should be provided by some pattern (track) reco algo.
    */
    std::vector<dg_wire>        fired_wires;
    /*
        impact_par_estimated : distance track (helix, 
        circle or sin) to the wire center.  
    */
    std::vector<double>         impact_par_estimated;

    /* local fit of 3 circles (segmented track)____________________   
        track_segments : vector of segments tangent to the
        measured impact parameters
    */
    
    // std::vector<Line2D>         track_segments_ZY;
    // std::vector<Line2D>         track_segments_XZ;
    
    Helix                       true_helix;
    Helix                       reco_helix;

    double                      pt_true;
    double                      pt_reco;

    TVector3                    p_true;
    TVector3                    p_reco;

    // fitting info TMinuit _______________________________________
    MinuitFitInfos              fit_infos_xz;
    MinuitFitInfos              fit_infos_zy;
};

namespace Color {
    enum Code {
        FG_RED      = 31,
        FG_GREEN    = 32,
        FG_BLUE     = 34,
        FG_DEFAULT  = 39,
        BG_RED      = 41,
        BG_GREEN    = 42,
        BG_BLUE     = 44,
        BG_DEFAULT  = 49
    };
    class Modifier {
        Code code;
    public:
        Modifier(Code pCode) : code(pCode) {}
        friend std::ostream&
        operator<<(std::ostream& os, const Modifier& mod) {
            return os << "\033[" << mod.code << "m";
        }
    };
}

#endif