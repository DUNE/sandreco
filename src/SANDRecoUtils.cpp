#include "SANDRecoUtils.h"

// using namespace RecoUtils;

// std::vector<dg_wire>* event_digits;

SANDGeoManager geo_manager;

void RecoUtils::InitWireInfos(TGeoManager* g){
    geo_manager.init(g);
}

double RecoUtils::GetDistHelix2Line(const Helix& helix, double s, const Line& line, double& t){
    // return the distance between a point on the helix at s
    // and a give line. The projection of the helix point on the line occurs at t

    auto point_h = helix.GetPointAt(s);
    auto point_l = line.GetLinePointX0();
    auto hl      = point_h - point_l;
    auto v       = line.GetDirectionVector();
    t = hl.Dot(v)/(v.Mag()*v.Mag());
    
    return (hl.Cross(v)).Mag()/v.Mag();
}

double RecoUtils::GetDistHelix2LineDerivative(const Helix& helix, double s, const Line& line){
    // return the derivateive of distance between a point on the helix
    // and a give line calculated for the helix at s
        
    auto point_h    = helix.GetPointAt(s);
    auto point_l    = line.GetLinePointX0();
    auto hl         = point_h - point_l;
    auto v          = line.GetDirectionVector();
    auto dh_over_ds = helix.GetTangentVector(s);
    auto hl_cross_v = hl.Cross(v);
    
    return hl_cross_v.Dot(dh_over_ds.Cross(v));
}

double RecoUtils::GetImpactParameter(const Helix& helix, const Line& line, double s, double t){

    TVector3 helix_point = helix.GetPointAt(s);
    TVector3 line_point  = line.GetPointAt(t);

    return (helix_point - line_point).Mag();
}

double RecoUtils::GetMinImpactParameter(const Helix& helix, const Line& line, double& s_min, double& t_min, bool& HasMinimized){
    /* given an helix and a wire line, find the impact parameter
       s_min and t_min are the parameter of the helix and line 
       that five the point on the helix and on the line corresponding
       to the impact parameter */
   
    // sample the helix in the drift plane
    const int nof_helix_sampling_points = 1000;

    auto s_step = (helix.UpLim() - helix.LowLim()) / nof_helix_sampling_points;

    // initialize the impact parameter
    double impact_parameter = 1e6;

    for (auto i = 0u; i < nof_helix_sampling_points; i++)
    {
        auto s_i = helix.LowLim() + s_step * i;

        double t_i = -999;

        auto dist2wire = RecoUtils::GetDistHelix2Line(helix, s_i, line, t_i);
        // if(dist2wire < impact_parameter)
        //     std::cout << "t_i : "                     << fabs(t_i) 
        //               << ", "                         << ((line.GetPointAt(t_i)-line.GetLinePointX0()).Mag())
        //               << ", line.GetLineLength()/2. " << line.GetLineLength()/2.
        //               << ", "                         << (fabs(t_i)<=line.GetLineLength()/2.)
        //               << ", "                         << ((dist2wire < impact_parameter)&&((line.GetPointAt(t_i)-line.GetLinePointX0()).Mag()))
        //               <<"\n";
    
        if((dist2wire < impact_parameter) && 
           ((line.GetPointAt(t_i)-line.GetLinePointX0()).Mag()<=line.GetLineLength()/2.) // is within SAND x range
           )
        {
            impact_parameter = dist2wire;
            s_min = s_i;
            t_min = t_i;
        }
    }
    if(s_min!=-999) HasMinimized=1;
    return impact_parameter;
}

double RecoUtils::GetMinImpactParameter(const Helix& helix, const Line& line){
    double s_min = -999;
    double t_min = -999;
    bool HasMinimized = 0;
    return RecoUtils::GetMinImpactParameter(helix, line, s_min, t_min, HasMinimized);
}

double RecoUtils::GetExpectedRadiusFromDigit(const dg_wire& digit){
    // get TDC and convert it to a radius
    TRandom3 rand;
    auto wire_info = geo_manager.get_wire_info(digit.did);
    
    double guess_wire_pmt_dist     = wire_info.length()/2.;
    double signal_propagation_time = guess_wire_pmt_dist/sand_reco::stt::v_signal_inwire;
    double t0                      = 1; // assume 1 ns
    return (digit.tdc - signal_propagation_time - t0)*sand_reco::stt::v_drift;
}

Line RecoUtils::GetLineFromDigit(const dg_wire& digit){
    // get Line from fired wire 
    // x = dx * t + ax
    // y = dy * t + ay
    // z = z0
    if(digit.hor==true)
    {
        // horizontal == line along x axis
        // x = 
        // y = y0
        // z = z0
        Line l(1., 0., digit.x, digit.y, digit.z);
        l.SetLineLength(digit.wire_length);
        return l;
    }else{// vertical
        Line l(0., 1., digit.x, digit.y, digit.z);
        l.SetLineLength(digit.wire_length);
        return l;
    }
}

double RecoUtils::NLL(Helix& h,const std::vector<dg_wire>& digits)
{   
    static int iter    = 0;
    double nll         = 0.;
    const double sigma = 0.2; // 200 mu_m = 0.2 mm
    int digit_index    = 0;

    for(auto& digit : digits) 
    {
        // each digit define a line completly
        Line l = RecoUtils::GetLineFromDigit(digit);

        /* find s_lower and s_upper that gives the portion of the helix 
        in the plane containing the wire */
        h.SetHelixRangeFromDigit(digit);
        
        double r_estimated = RecoUtils::GetMinImpactParameter(h,l);
        
        // double r_measured  = RecoUtils::GetExpectedRadiusFromDigit(digit);
        double r_true = digit.drift_time * sand_reco::stt::v_drift;
        
        // nll += (r_estimated - r_measured) * (r_estimated - r_measured) / (sigma * sigma);
        nll += (r_estimated - r_true) * (r_estimated - r_true) / (sigma * sigma);
        
        digit_index++;
    }
    // if(iter==2) throw "";
    iter ++;
    return sqrt(nll)/digits.size();
}

double RecoUtils::FunctorNLL(const double* p)
{
    double R    = p[0];
    double dip  = p[1];
    double Phi0 = p[2];
    int hel     = p[3];
    TVector3 x0 = {p[4],p[5],p[6]};

    Helix h(R, dip, Phi0, hel, x0);

    auto nll = RecoUtils::NLL(h,  *event_digits);
    
    return nll;
}

const double* RecoUtils::InitHelixPars(const std::vector<dg_wire>& digits)
{
    double* p;
    return p;
}

const double* RecoUtils::GetHelixParameters(const Helix& helix_initial_guess, int& TMinuitStatus)
{
    ROOT::Math::Functor functor(&RecoUtils::FunctorNLL, 7);

    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");

    minimizer->SetPrintLevel(4);

    // minimizer->SetMaxIterations(100);
    // minimizer->SetMaxFunctionCalls(100);

    minimizer->SetFunction(functor);

    // helix params
    // SetLimitedVariable       (ivar, name,   val,  step,  low,   up)
    minimizer->SetLimitedVariable(0,    "R",    helix_initial_guess.R(), 100,  100,     1e5);
    minimizer->SetLimitedVariable(1,    "dip",  helix_initial_guess.dip(), 0.01,  -1.6,  1.6);
    // minimizer->SetLimitedVariable(2,    "Phi0", helix_initial_guess.Phi0(), 0.01,  -3.14, 3.14);
    // minimizer->SetLimitedVariable(3,    "h",    helix_initial_guess.h(), 0.01,  -1.1,  1.1);
    minimizer->SetLimitedVariable(4,    "x0_x", helix_initial_guess.x0().X(), 1,  -1800, 1800);
    // minimizer->SetLimitedVariable(5,    "x0_y", helix_initial_guess.x0().Y(), 1,  -4500, -300);
    // minimizer->SetLimitedVariable(6,    "x0_z", helix_initial_guess.x0().Z(), 1,  21500, 26000);

    // minimizer->SetFixedVariable(0,    "R",    helix_initial_guess.R());
    // minimizer->SetFixedVariable(1,    "dip",    helix_initial_guess.dip());
    minimizer->SetFixedVariable(2,    "Phi0",    helix_initial_guess.Phi0());
    minimizer->SetFixedVariable(3,    "h",    helix_initial_guess.h());
    // minimizer->SetFixedVariable(4,    "x0_x",    helix_initial_guess.x0().X());
    minimizer->SetFixedVariable(5,    "x0_y",    helix_initial_guess.x0().Y());
    minimizer->SetFixedVariable(6,    "x0_z",    helix_initial_guess.x0().Z());

    minimizer->Minimize();
    minimizer->PrintResults();

    TMinuitStatus      = minimizer->Status();  
    const double* pars = minimizer->X();

    return pars;
}