#include "SANDRecoUtils.h"

// using namespace RecoUtils;

// std::vector<dg_tube>* event_digits;

SANDGeoManager geo_manager;

void RecoUtils::InitWireInfos(TGeoManager* g){
    geo_manager.init(g);
}

double RecoUtils::GetDistHelix2Line(const Helix& helix, double s, const Line& line){
    // return the distance between a point on the helix at s
    // and a give line

    auto point_h = helix.GetPointAt(s);
    auto point_l = line.GetLinePointX0();
    auto hl      = point_h - point_l;
    auto v       = line.GetDirectionVector();
    
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

// double RecoUtils::FunctorImpactParameter(const double* p) {
//     // first 7 entries of p to define the helix: R, dip, Phi0, h, x0
//     // last 4 entries of p to define the line: dx, dy, ax, ay
//     TVector3 x0{p[4], p[5], p[6]};

//     Helix helix(p[0], p[1], p[2], p[3], x0);
//     Line line(p[7],p[8],p[9],p[10],p[11]);

//     double s = p[12];
//     double t = p[13];

//     return RecoUtils::GetImpactParameter(helix, line, s, t);
// }

// double RecoUtils::GetMinImpactParameter(const Helix& helix, const Line& line, double& s_min, double& t_min, bool& HasMinimized){

//     ROOT::Math::Functor functor(&RecoUtils::FunctorImpactParameter, 14);

//     ROOT::Math::Minimizer * minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");

//     minimizer->SetFunction(functor);

//     // std::cout<<__FILE__<<" "<<__LINE__<<" RecoUtils::GetMinImpactParameter\n";
//     // helix.PrintHelixPars();
    
//     // helix params
//     minimizer->SetFixedVariable(0, "R",     helix.R());
//     minimizer->SetFixedVariable(1, "dip",   helix.dip());
//     minimizer->SetFixedVariable(2, "Phi0",  helix.Phi0());
//     minimizer->SetFixedVariable(3, "h",     helix.h());
//     minimizer->SetFixedVariable(4, "x0_x",  helix.x0().X());
//     minimizer->SetFixedVariable(5, "x0_y",  helix.x0().Y());
//     minimizer->SetFixedVariable(6, "x0_z",  helix.x0().Z());
//     // line params    
//     minimizer->SetFixedVariable(7, "dx",    line.dx());
//     minimizer->SetFixedVariable(8, "dy",    line.dy());
//     minimizer->SetFixedVariable(9, "ax",    line.ax());
//     minimizer->SetFixedVariable(10, "ay",   line.ay());
//     minimizer->SetFixedVariable(11, "z0",   line.z0());
//     // s, t
//     double s_low = helix.LowLim();
//     double s_up = helix.UpLim();
//     double t_low = line.LowLim();
//     double t_up = line.UpLim();

//     minimizer->SetVariable(12, "s", 1, 0.0001);
//     minimizer->SetVariable(13, "t", 1, 0.0001);
//     // minimizer->SetVariable(12, "s", (s_low+s_up)*0.5, 0.0001);
//     // minimizer->SetVariable(13, "t", (t_low+t_up)*0.5, 0.0001);

//     // set limits
//     // std::cout<<"setting limits on s : ["<<s_low<<","<<s_up<<"]\n";
//     // std::cout<<"setting limits on t : ["<<t_low<<","<<t_up<<"]\n";
//     // minimizer->SetVariableLimits(12, s_low-3, s_low+3);
//     // minimizer->SetVariableLimits(13, t_low, t_up);

//     HasMinimized = minimizer->Minimize();

//     if(HasMinimized) {
//         // minimizer->PrintResults();
//         auto pars = minimizer->X(); 
//         s_min = pars[12]; 
//         t_min = pars[13];
//         return GetImpactParameter(helix,line,s_min,t_min);
//         // std::cout<<"found s : "<<s_min<<" found t : "<<t_min<<" that gives the distance : "<<GetImpactParameter(helix,line,s_min,t_min)<<"\n";
//         // std::cout<<"point on the Helix  : "<<helix.GetPointAt(s_min).X()<<" "<<helix.GetPointAt(s_min).Y()<<" "<<helix.GetPointAt(s_min).Z()<<"\n";
//         // std::cout<<"point on the Line   : "<<line.GetPointAt(t_min).X()<<" "<<line.GetPointAt(t_min).Y()<<" "<<line.GetPointAt(t_min).Z()<<"\n";
//     }else{
//         std::cout<<"wasnt' able to minimize\n";
//         throw "";
//     };
// }

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

        auto dist2wire = RecoUtils::GetDistHelix2Line(helix, s_i, line);

        // std::cout << "dist2wire : " << dist2wire << ", s_i : " << s_i << "\n";
    
        if(dist2wire < impact_parameter)
        {
            impact_parameter = dist2wire;
            s_min = s_i;
            // what about t_min ?
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

double RecoUtils::GetExpectedRadiusFromDigit(const dg_tube& digit){
    // get TDC and convert it to a radius
    TRandom3 rand;
    auto wire_info = geo_manager.get_wire_info(digit.did);

    // std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
    // std::cout<<"wire length : "<<wire_info.length()<<"\n";
    // std::cout<<"signal velocity  : "<<sand_reco::stt::v_signal_inwire<<"\n";
    // std::cout<<"drift velocity  : "<<sand_reco::stt::v_drift<<"\n";
    // std::cout<<"tdc : "<<digit.tdc<<"\n";
    // std::cout<<"t0 : "<<digit.t0<<"\n";

    double guess_wire_pmt_dist = wire_info.length()/2.;
    double signal_propagation_time = guess_wire_pmt_dist/sand_reco::stt::v_signal_inwire;
    double t0 = 1; // assume 1 ns
    return (digit.tdc - signal_propagation_time - t0)*sand_reco::stt::v_drift;
}

Line RecoUtils::GetLineFromDigit(const dg_tube& digit){
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
        return l;
    }else{// vertical
        Line l(0., 1., digit.x, digit.y, digit.z);
        return l;
    }
}

double RecoUtils::NLL(Helix& h,const std::vector<dg_tube>& digits)
{   
    static int iter    = 0;
    double nll         = 0.;
    const double sigma = 0.2; // 200 mu_m = 0.2 mm
    int digit_index    = 0;

    // std::cout <<"\n";
    // std::cout << "iteration number : "<< iter << "\n";
    // std::cout << iter << "," 
    //           << h.R() << ","
    //           << h.dip() << ","
    //           << h.Phi0() << ","
    //           << h.x0().X() << ","
    //           << h.x0().Y() << ","
    //           << h.x0().Z() << "\n";
    // h.PrintHelixPars();

    for(auto& digit : digits) 
    {
        // each digit define a line completly
        Line l = RecoUtils::GetLineFromDigit(digit);

        /* find s_lower and s_upper that gives the portion of the helix 
        in the plane containing the wire */
        h.SetHelixRangeFromDigit(digit);
        
        double r_estimated = RecoUtils::GetMinImpactParameter(h,l);
        
        double r_measured  = RecoUtils::GetExpectedRadiusFromDigit(digit);
        
        nll += sqrt((r_estimated - r_measured) * (r_estimated - r_measured) / (sigma * sigma));
        
        // std::cout << "digit : << "         << digit_index <<
        //             //  ", s_low_lim : "      << h.LowLim() << 
        //             //  ", s_up_lim : "       << h.UpLim() << 
        //              ", r_estimated : "    << r_estimated << 
        //              ", r_measured : "     << r_measured << 
        //              ", cumulative nll : " << nll << 
        //              "\n";
        
        digit_index++;
    }
    // if(iter==2) throw "";
    iter ++;
    return nll;
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

const double* RecoUtils::InitHelixPars(const std::vector<dg_tube>& digits)
{
    double* p;
    return p;
}

// const double* RecoUtils::GetHelixParameters(const double* p, const std::vector<dg_tube>& digits)
const double* RecoUtils::GetHelixParameters(const double* p)
{
    // std::cout<<__FILE__<<" "<<__LINE__<<" RecoUtils::GetHelixParameters\n";

    ROOT::Math::Functor functor(&RecoUtils::FunctorNLL, 7);

    ROOT::Math::Minimizer * minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");

    minimizer->SetPrintLevel(4);

    minimizer->SetFunction(functor);

    // helix params
    // >SetLimitedVariable       (ivar, name,   val,  step,  low,   up)
    minimizer->SetLimitedVariable(0,    "R",    p[0], 100,  100,     1e5);
    minimizer->SetLimitedVariable(1,    "dip",  p[1], 0.01,  -1.6,  1.6);
    // minimizer->SetLimitedVariable(2,    "Phi0", p[2], 0.01,  -3.14, 3.14);
    // minimizer->SetLimitedVariable(3,    "h",    p[3], 0.01,  -1.1,  1.1);
    // minimizer->SetLimitedVariable(4,    "x0_x", p[4], 1,  -1800, 1800);
    // minimizer->SetLimitedVariable(5,    "x0_y", p[5], 1,  -4500, -300);
    // minimizer->SetLimitedVariable(6,    "x0_z", p[6], 1,  21500, 26000);

    // minimizer->SetVariable(0,    "R",    p[0], 1);
    // minimizer->SetFixedVariable(0,    "R",    p[0]);
    // minimizer->SetFixedVariable(1,    "dip",    p[1]);
    minimizer->SetFixedVariable(2,    "Phi0",    p[2]);
    minimizer->SetFixedVariable(3,    "h",    p[3]);
    minimizer->SetFixedVariable(4,    "x0_x",    p[4]);
    minimizer->SetFixedVariable(5,    "x0_y",    p[5]);
    minimizer->SetFixedVariable(6,    "x0_z",    p[6]);

    minimizer->Minimize();

    minimizer->PrintResults();

    const double* pars = minimizer->X();

    return pars;
}