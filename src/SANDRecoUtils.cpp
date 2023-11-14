#include "SANDRecoUtils.h"

std::vector<dg_tube>* event_digits = nullptr;

SANDGeoManager* geo_manager;

double RecoUtils::GetImpactParameter(const Helix& helix, const Line& line, double s, double t) {

    TVector3 helix_point = helix.GetPointAt(s);
    TVector3 line_point = line.GetPointAt(t);

    return (helix_point - line_point).Mag();
}

double RecoUtils::FunctorImpactParameter(const double* p) {
    // first 7 entries of p to define the helix: R, dip, Phi0, h, x0
    // last 4 entries of p to define the line: dx, dy, ax, ay
    TVector3 x0{p[4], p[5], p[6]};

    Helix helix(p[0], p[1], p[2], p[3], x0);
    Line line(p[7],p[8],p[9],p[10],p[11]);

    double s = p[12];
    double t = p[13];

    return RecoUtils::GetImpactParameter(helix, line, s, t);
}

double RecoUtils::GetMinImpactParameter(const Helix& helix, const Line& line){

    ROOT::Math::Functor functor(&RecoUtils::FunctorImpactParameter, 13);

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

    if(HasMinimized) {
        minimizer->PrintResults();
        auto s_min = minimizer->X()[12];
        auto t_min = minimizer->X()[13];
        std::cout<<"found s : "<<s_min<<"\n";
        std::cout<<"found t : "<<t_min<<"\n";
        std::cout<<"that gives the distance : "<<GetImpactParameter(helix,line,s_min,t_min)<<"\n";
    };

    return HasMinimized;
}

double RecoUtils::GetExpectedRadiusFromDigit(const dg_tube& digit){
    // get TDC and convert it to a radius
    TRandom3 rand;
    auto wire_info = geo_manager->get_wire_info(digit.did);
    double guess_wire_pmt_dist = wire_info.length()/2.;
    double signal_propagation_time = guess_wire_pmt_dist/sand_reco::stt::v_signal_inwire;
    double smearing = rand.Gaus(0, sand_reco::stt::tm_stt_smearing);
    
    return (digit.tdc - digit.t0 - signal_propagation_time - smearing)/sand_reco::stt::v_drift;
}

Line RecoUtils::GetLineFromDigit(const dg_tube& digit){
    // get Line from fired wire 
    // x = dx * t + ax
    // y = dy * t + ay
    // z = z0
    if(digit.hor==true){
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

double RecoUtils::NLL(const Helix& h,const std::vector<dg_tube>& digits)
{   
    double nll = 0.;

    const double sigma = 200.*1e-6;

    for(auto& digit : digits)
    {
        Line l = RecoUtils::GetLineFromDigit(digit); // each digit define a line completly
        double r_estimated = RecoUtils::GetMinImpactParameter(h,l);
        double r_expected  = RecoUtils::GetExpectedRadiusFromDigit(digit);
        nll += (r_estimated - r_expected) * (r_estimated - r_expected) / (sigma*sigma); 
    }

    return nll;
}

double RecoUtils::FunctorNLL(const double* p)
{
    double R = p[0];
    double dip = p[1];
    double Phi0 = p[2];
    int hel = p[3];
    TVector3 x0 = {p[4],p[5],p[6]};

    Helix h(R, dip, Phi0, hel, x0);

    return RecoUtils::NLL(h, *event_digits);
}

const double* RecoUtils::InitHelixPars(const std::vector<dg_tube>& digits)
{
    double* p;
    return p;
}

const double* RecoUtils::GetHelixParameters(const double* p)
{
    // pass first guess to the Helix parmeters
    ROOT::Math::Functor functor(&RecoUtils::FunctorNLL, 5);

    ROOT::Math::Minimizer * minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");

    minimizer->SetFunction(functor);

    // helix params
    minimizer->SetVariable(0, "R", p[0], 0.01);
    minimizer->SetVariable(1, "dip", p[1], 0.01);
    minimizer->SetVariable(2, "Phi0", p[2], 0.01);
    minimizer->SetVariable(3, "h", p[3], 1);
    minimizer->SetVariable(4, "x0_x", p[4], 0.01);
    minimizer->SetVariable(5, "x0_y", p[5], 0.01);
    minimizer->SetVariable(6, "x0_z", p[6], 0.01);

    if(minimizer->Minimize()) std::cout<<"minimizer successfull\n";

    return minimizer->X();
}