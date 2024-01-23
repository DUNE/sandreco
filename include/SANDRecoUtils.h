#ifndef SANDRECOUTILS_H
#define SANDRECOUTILS_H

#include <cmath>
#include <TRandom3.h>
// #include "utils.h"
#include "struct.h"
#include "utils.h"
#include "SANDSTTTubeInfo.h"
#include "SANDGeoManager.h"

#include "TVector3.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TDatabasePDG.h"

extern TGeoManager* geo;

extern SANDGeoManager geo_manager;

class Helix
{
    public:
        // R : radius of the helix
        // dip : dip angle 
        // Phi0 : is the azimuth angle of the starting point with respect to the helix axis
        // h : +/- sense of rotation
        // x0 : track starting point

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
            h_            = (charge < 0) ? 1 : -1;  // to check muon helicity consistency
            Phi0_         = TMath::ATan2(momentum.Y(),momentum.Z()) + h_*TMath::Pi()*0.5; // for mu pi/2 should be added
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
        double dx_over_ds(double s) const{
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
            return {dx_over_ds(s), dy_over_ds(s), dz_over_ds(s)};
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
            std::cout<<std::setprecision(8)<<"x0, y0, z0 -> "<<x0_.X()<<" "<<x0_.Y()<<" "<<x0_.Z()<<"\n";
        }

        double   R()      const {return R_;};
        double   dip()    const {return dip_;};
        double   Phi0()   const {return Phi0_;};
        int      h()      const {return h_;};
        TVector3 x0()     const {return x0_;};
        double   LowLim() const {return low_lim_;};
        double   UpLim()  const {return up_lim_;};

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

class Line : public SANDWireInfo
{
    public:
        Line(double arg_dx, double arg_dy, double arg_ax, double arg_ay, double arg_z0, double low_lim = -999., double up_lim = -999.) : 
                    dx_(arg_dx), dy_(arg_dy), ax_(arg_ax), ay_(arg_ay), z0_(arg_z0), low_lim_(low_lim), up_lim_(up_lim){}

        // (dx,dy,dz) = vector parallel to wire
        // x = dx * t + ax
        // y = dy * t + ay
        // z = z0
        // a,b : line intercept

        // optional : low_lim_, up_lim_
        // usefull as we are interest in a portion of the line defined by the wire length

        void SetWireLineParam(double* p){
            dx_ = p[0];
            dy_ = p[1];
            ax_ = p[2];
            ay_ = p[3];
            z0_ = p[4];
        }

        void SetLowLim(double arg_low){
            low_lim_ = arg_low;
        }
        
        void SetUpLim(double arg_up){
            up_lim_ = arg_up;
        }

        void SetLineRangeFromDigit(const dg_wire& digit){
            auto wire = geo_manager.get_wire_info(digit.did);
            SetLowLim(-wire.length()/2. * 1.01); // 1.01 to be not too stringent on the limits
            SetUpLim(wire.length()/2. * 1.01);
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

        // first derivatives
        double dx_over_dt() const {
            return dx_;
        }
        double dy_over_dt() const {
            return dy_;
        }
        double dz_over_dt() const {
            return 0.;
        }
        // derivatives

        TVector3 GetDirectionVector() const{
            return {dx_over_dt(), dy_over_dt(), dz_over_dt()};
        }

        TVector3 GetLinePointX0() const{
            return GetPointAt(0);
        }

        double dx() const {return dx_;};
        double dy() const {return dy_;};
        double ax() const {return ax_;};
        double ay() const {return ay_;};
        double z0() const {return z0_;};
        double LowLim() const {return low_lim_;};
        double UpLim() const {return up_lim_;};

        TVector3 GetPointAt(double t) const {
            return {x_l(t), y_l(t), z_l()};
        }

    private:
        double dx_;
        double dy_;
        double ax_;
        double ay_;
        double z0_;
        double low_lim_;
        double up_lim_;
    
    ClassDef(Line, 1);
};

namespace RecoUtils{ //RecoUtils

extern        std::vector<dg_wire>* event_digits;

void          InitWireInfos(TGeoManager* g);

double        GetDistHelix2Line(const Helix& helix, double s, const Line& line);

double        GetDistHelix2LineDerivative(const Helix& helix, double s, const Line& line); // may be can be remover

double        GetImpactParameter(const Helix& helix, const Line& line, double s, double t);

double        FunctorImpactParameter(const double* p);

double        GetMinImpactParameter(const Helix& helix, const Line& line);

double        GetMinImpactParameter(const Helix& helix, const Line& line, double& s_min, double& t_min, bool& HasMinimized);

double        GetExpectedRadiusFromDigit(const dg_wire& digit);

Line          GetLineFromDigit(const dg_wire& digit);

double        NLL(Helix& h,const std::vector<dg_wire>& digits); // negative log likelihood function

double        FunctorNLL(const double* p);

const double* InitHelixPars(const std::vector<dg_wire>& digits);

const double* GetHelixParameters(const Helix& helix_initial_guess, int& TMinuitStatus);

} // RecoUtils

struct MinuitFitInfos
{
    int                         TMinuitFinalStatus; // 0 if converged, 4 if falied
    int                         NIterations; // number of iterations to reach the minimum
    double                      MinValue; // value of the func to minimize at its minimum
    std::vector<double>         parameters_errors; // error associated to the par estimate
};

struct RecoObject
{
    int                         traj_edep_index;
    std::vector<dg_wire>        fired_wires;
    Helix                       true_helix;
    Helix                       helix_first_guess;
    Helix                       reco_helix;
    std::vector<double>         impact_par_from_TDC;
    std::vector<double>         impact_par_estimated;
    std::vector<TLorentzVector> trj_points; // from edepsim
    double                      pt_true;
    double                      pt_reco;
    MinuitFitInfos              fit_infos;
};

struct EventReco
{   
    int                         event_index;
    RecoObject                  reco_object;
    std::vector<dg_wire>        event_fired_wires;
    int                         nof_digits = event_fired_wires.size();
};


#endif