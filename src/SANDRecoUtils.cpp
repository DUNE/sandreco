#include "SANDRecoUtils.h"
#include <numeric>

SANDGeoManager geo_manager;

// constants___________________________________________________________
const double SAND_CENTER_X = 0.;

const double SAND_CENTER_Y = -2384.73;

const double SAND_CENTER_Z = 23910.;

const double SAND_TRACKER_X_LENGTH = 3220.0; // does not include frames
//______________________________________________________________________

void RecoUtils::InitWireInfos(TGeoManager* g){
    geo_manager.init(g);
}

std::vector<double> RecoUtils::SmearVariable(double mean, double sigma, int nof_points){

    TRandom3 r(0);

    std::vector<double> smeared_points;

    for (int i = 0; i < nof_points; ++i) 
        smeared_points.push_back(r.Gaus(mean, sigma));
    
    return smeared_points;
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

dg_wire  RecoUtils::Copy(const dg_wire& wire){
    dg_wire copy;
    copy.det = wire.did;
    copy.did = wire.did;
    copy.x = wire.x;
    copy.y = wire.y;
    copy.z = wire.z;
    copy.t0 = wire.t0;
    copy.de = wire.de;
    copy.adc = wire.adc;
    copy.tdc = wire.tdc;
    copy.t0 = wire.t0;
    copy.hor = wire.hor;
    copy.hindex = wire.hindex;
    copy.drift_time = wire.drift_time;
    copy.signal_time = wire.signal_time;
    copy.t_hit = wire.t_hit;
    copy.wire_length = wire.wire_length;
    return copy; 
}

double RecoUtils::GetLinePointDistance(const Line& l, TVector3 point){
    auto x = l.GetPointAt(0);
    auto v = l.GetDirectionVector();
    auto xp = x - point;
    return (xp.Cross(v)).Mag()/v.Mag();
}

double RecoUtils::GetLineLineDistance(const Line& l1, const Line& l2, 
                                      double& closest2line2, double& closest2line1){
    /*
        Return the (shortest) distance between 2 lines l1 and l2 and 
        fill the values closest2line2 and closest2line1:
        - closest2line2 : is the double that gives the point on the line 1
                          closest to the line 2
        - closest2line1 : is the double that gives the point on the line 2
                          closest to the line 1
        
        Reference : https://math.stackexchange.com/questions/2213165/find-shortest-distance-between-lines-in-3d

        NOTE : the direction vector of the lines l1 and l2 are not supposed to be 
               necessarely unitary. So closest2line1 and closest2line2 are NOT
               necessarely the distance from the lines X0 point (usually taken
               as the line center.
    */

   // define direction of the line that connects the 2 lines
   auto e1 = l1.GetDirectionVector();
   auto e2 = l2.GetDirectionVector(); 
   auto n = e1.Cross(e2);
//    n *= 1./n.Mag();

   // define vector connecting the 2 lines
   auto l21 = l2.GetLinePointX0() - l1.GetLinePointX0();

   // point on l1 closest to l2
   closest2line2 = (e2.Cross(n)).Dot(l21)/(n.Dot(n));
   
   // point on l2 closest to l1
   closest2line1 = (e1.Cross(n)).Dot(l21)/(n.Dot(n));

   return (l1.GetPointAt(closest2line2) - l2.GetPointAt(closest2line1)).Mag();
}

double RecoUtils::GetSegmentSegmentDistance(const Line& l1, const Line& l2, 
                                      double& closest2line2, double& closest2line1){
    /*
        Return the distance between to lines defined in a given range 
        = distance between two segments
    */                                        
    auto d = GetLineLineDistance(l1, l2, closest2line2, closest2line1);
    double seg_closest2line1, seg_closest2line2;
    
    // point on l1 closest to line l2 is outside l1 range
    if(!l1.IsInLineLimits(closest2line2)){
        // closest point is either up or low lim
        auto d_up  = (l1.GetLineUpperLimit() - l2.GetPointAt(closest2line1)).Mag();
        auto d_low = (l1.GetLineLowerLimit() - l2.GetPointAt(closest2line1)).Mag();
        seg_closest2line2 = (d_up < d_low) ? l1.UpLim() : l1.LowLim();
    }else{
        seg_closest2line2 = closest2line2;
    }

    // point on l2 closest to line l1 is outside l2 range
    if(!l2.IsInLineLimits(closest2line1)){
        // closest point is either up or low lim
        auto d_up  = (l2.GetLineUpperLimit() - l1.GetPointAt(closest2line2)).Mag();
        auto d_low = (l2.GetLineLowerLimit() - l1.GetPointAt(closest2line2)).Mag();
        seg_closest2line1 = (d_up < d_low) ? l2.UpLim() : l2.LowLim();
    }else{
        seg_closest2line1 = closest2line1;
    }
    
    d = (l1.GetPointAt(seg_closest2line2) - l2.GetPointAt(seg_closest2line1)).Mag();
    
    closest2line1 = seg_closest2line1;
    closest2line2 = seg_closest2line2;
    
    return d;
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

double RecoUtils::NewtonRaphson2D(TF1* f, TF1* fprime, double& x_guess, 
                                  double tol, unsigned int max_iterations){
    /*
        given a function f and its derivative fprime
        find the point x_sol such that fprime(x_sol) 
        approx 0 (within a given tollerance tol and 
        with a max number of iterations max_iterations).
        Return f(x_sol)
    */
    // std::cout << "x_guess " << x_guess << "f->Eval(x_guess) "<< f->Eval(x_guess) << "\n";
    double min_sofar = 1e6;
    unsigned int iter = 0;
    for (auto i = 0u; i < max_iterations; i++)
    {
        double f_at_x = f->Eval(x_guess);
        double fprime_at_x = fprime->Eval(x_guess);
        if(f_at_x < min_sofar) min_sofar = f_at_x;
        // std::cout << "i : " << i 
        //   << ", x : " << x_guess
        //   << ", distance : " << f_at_x
        //   << ", fprime_at_x : " << fprime_at_x 
        //   << "\n";
        if(fabs(fprime_at_x) < tol) break;
        x_guess -= f_at_x / fprime_at_x;
        if (fabs(f_at_x) < tol) break;
        iter++;
    }
    if(iter==max_iterations - 1u) return min_sofar;
    return f->Eval(x_guess);
}

double RecoUtils::GetMinImpactParameter(const Helix& helix, const Line& line){
    double s_min = -999;
    double t_min = -999;
    bool HasMinimized = 0;
    return RecoUtils::GetMinImpactParameter(helix, line, s_min, t_min, HasMinimized);
}

Line RecoUtils::GetLineFromWire(const dg_wire& wire){
    // get Line from fired wire 
    // x = dx * t + ax
    // y = dy * t + ay
    // z = z0
    if(wire.hor==true)
    {
        // horizontal == line along x axis
        // x = t + x0
        // y = y0
        // z = z0
        Line l(1., 0., 0., wire.x, wire.y, wire.z);
        l.SetLineLength(wire.wire_length);
        return l;
    }else{
        // vertical == line along y axis
        // x = x0
        // y = t + y0
        // z = z0
        Line l(0., 1., 0., wire.x, wire.y, wire.z);
        l.SetLineLength(wire.wire_length);
        return l;
    }
}

TF1* RecoUtils::InvertSin(TF1* fSin){
    /*
        - expected input z(x) = A *sin(B*x + C) + D
        - output x(z) = [arcsin((z - D)/A) - C] / B
        with output function define in x half period
    */
    double amplitude = fSin->GetParameter(0); // A
    double frequency = fSin->GetParameter(1); // B
    double phase = fSin->GetParameter(2); // C
    double offset = fSin->GetParameter(3); // D

    double z_min = offset - amplitude;// for sin = -1
    double z_max = offset + amplitude;

    // std::cout << "function [x_min, x_max] = [" << fSin->GetXmin() << ", " << fSin->GetXmax() << "]\n"; 
    // std::cout << "function [z_min, z_max] = [" << fSin->Eval(fSin->GetXmin()) << ", " << fSin->Eval(fSin->GetXmax()) << "]\n"; 
    // std::cout << "inverse function [z_min, z_max] = [" << z_min << ", " << z_max << "]\n"; 
    
    // TF1* fAsin = new TF1("fAsin", "(( asin((x - [0])/[1]) - [2] ) / [3] )", z_min, z_max);
    TF1* fAsin = new TF1("fAsin", "((asin((x-[0])/[1]) + 2*TMath::Pi()*[2])-[3])/[4]");
    
    fAsin->SetParameter(0, offset);
    fAsin->SetParameter(1, amplitude);
    fAsin->SetParameter(2, -1);
    fAsin->SetParameter(3, phase);
    fAsin->SetParameter(4, frequency);

    TCanvas *c1 = new TCanvas("c2", "ASine Function", 800, 600);
    fAsin->Draw();
    c1->SaveAs("inverse_fitted_sin.png");
    
    return fAsin;
}

TF1* RecoUtils::WiresSinFit(const std::vector<dg_wire>& wires){
    /*
        Perform a 2D fit with sin function of the XZ coordinates
        of the verticlal wire
        
        PARAMETERS OF THE FITS: [0]*sin([1]*x + [2]) + [3]
        - [0] Amplitude
        - [1] Frequency
        - [2] Phase
        - [3] Offset

        DATA : (X,Z) coordinates of vertical fired wires
    */
    std::vector<double> z,x;
    for(auto& wire : wires){
        x.push_back(wire.x);
        z.push_back(wire.z);
    }
    std::cout << "number of zy hits: " << z.size() << "\n";

    // Create a TGraph object and fill it with the data 
    TGraph *graph = new TGraph(z.size(), &x[0], &z[0]);

    double x_min = *std::min_element(x.begin(), x.end());
    double x_max = *std::max_element(x.begin(), x.end());
    
    double z_min = *std::min_element(z.begin(), z.end());
    double z_max = *std::max_element(z.begin(), z.end());

    double initial_amplitude_guess = (z_max - z_min)/2.;
    double initial_frequency_guess = 2. * M_PI/(x_max - x_min);
    double initial_phase_guess = 0.;
    double initial_offset_guess = std::accumulate(z.begin(), z.end(), 0.0) / z.size();

    // Define a sine function with parameters p0, p1, p2, and p3
    TF1* sineFit = new TF1("sineFit", "[0]*sin([1]*x + [2]) + [3]", x_min, x_max);

    // std::cout << "sinFit [x_min, x_max] : [" << x_min << ", " << x_max << "]\n";

    // Set initial parameters for the fit
    sineFit->SetParameters(initial_amplitude_guess, 
                           initial_frequency_guess, 
                           initial_phase_guess, 
                           initial_offset_guess);

    // Perform the fit
    graph->Fit(sineFit, "Q"); // "Q" option suppresses the fit statistics output

    // Retrieve fit parameters
    double amplitude = sineFit->GetParameter(0);
    double frequency = sineFit->GetParameter(1);
    double phase = sineFit->GetParameter(2);
    double offset = sineFit->GetParameter(3);

    // Visualize the fit
    TCanvas *c1 = new TCanvas("c1", "Sine Function Fit", 800, 600);
    graph->Draw("AP");
    sineFit->Draw("same");

    // Print out the fit parameters
    // std::cout << "Fit parameters:" << std::endl;
    // std::cout << "Amplitude: " << amplitude << std::endl;
    // std::cout << "Frequency: " << frequency << std::endl;
    // std::cout << "Phase: " << phase << std::endl;
    // std::cout << "Offset: " << offset << std::endl;
 
    // Optionally save the plot and parameters to a file
    c1->SaveAs("sine_function_fit.png");
    return sineFit;
}

int fitCircle(int n, const std::vector<double>& x, const std::vector<double>& y,
              double& xc, double& yc, double& r, double& errr, double& chi2)
{
  xc = -999;
  yc = -999;
  r = -999;
  errr = -999;
  chi2 = -999;

  if (x.size() != y.size()) return 1;

  double sumx = 0, sumy = 0;                            // linear    terms
  double sumx2 = 0, sumy2 = 0, sumxy = 0;               // quadratic terms
  double sumxy2 = 0, sumx2y = 0, sumx3 = 0, sumy3 = 0;  // cubic     terms

  for (int i = 0; i < n; i++) {
    double xp = x.at(i);
    double yp = y.at(i);
    sumx += xp;
    sumy += yp;
    sumx2 += xp * xp;
    sumy2 += yp * yp;
    sumxy += xp * yp;
    sumxy2 += xp * yp * yp;
    sumx2y += xp * xp * yp;
    sumx3 += xp * xp * xp;
    sumy3 += yp * yp * yp;
  }

  double a = n * sumx2 - sumx * sumx;
  double b = n * sumxy - sumx * sumy;
  double c = n * sumy2 - sumy * sumy;
  double d = 0.5 * (n * sumxy2 - sumx * sumy2 + n * sumx3 - sumx * sumx2);
  double e = 0.5 * (n * sumx2y - sumy * sumx2 + n * sumy3 - sumy * sumy2);

  if (a * c - b * b == 0.) return 2;

  xc = (d * c - b * e) / (a * c - b * b);
  yc = (a * e - b * d) / (a * c - b * b);

  double rMean = 0;
  double rrms = 0;

  for (int i = 0; i < n; i++) {
    double xp = x.at(i);
    double yp = y.at(i);
    double r2 = (xp - xc) * (xp - xc) + (yp - yc) * (yp - yc);

    rMean += sqrt(r2);
    rrms += r2;
  }

  rMean /= n;
  rrms /= n;
  r = rMean;

  errr = sqrt(rrms - rMean * rMean);

  chi2 = 0.0;

  for (int i = 0; i < n; i++) {
    chi2 += TMath::Abs((y.at(i) - yc) * (y.at(i) - yc) +
                       (x.at(i) - xc) * (x.at(i) - xc) - r * r);
  }

  chi2 /= n;

  return 0;
}

Circle RecoUtils::WiresCircleFit(const std::vector<dg_wire*>& wires){
    /*
        Perform a 2D circular fit of the ZY coordinates
        of the horizontal wires 
    */
    std::vector<double> z,y;
    for(auto& wire : wires){
        z.push_back(wire->z);
        y.push_back(wire->y);
    }

    double centerZ, centerY, radius, err, chi2;

    fitCircle(z.size(), z, y, centerZ, centerY, radius, err, chi2);

    Circle c(centerZ, centerY, radius);

    return c;
}

Line2D RecoUtils::WiresLinearFit(const std::vector<dg_wire*>& wires){
    /*
        Perform a 2D linear fit of the XZ coordinates
        of the vertical wire 
    */
   std::vector<double> z,x;
   for(auto& wire : wires){
    if(wire->hor==false){ // vertical
        x.push_back(wire->x);
        z.push_back(wire->z);
    }
   }
    // Create a TGraph object and fill it with the data
    TGraph* graph = new TGraph(z.size(), &x[0], &z[0]);

    // Define a linear function with parameters p0 (slope) and p1 (intercept)
    TF1* linearFit = new TF1("linearFit", "[0]*x + [1]");

    // Perform the fit using chi-square minimization method
    graph->Fit(linearFit, "Q"); // "Q" option suppresses the fit statistics output

    // Retrieve fit parameters
    double slope = linearFit->GetParameter(0);
    double intercept = linearFit->GetParameter(1);

    // Visualize the fit
    // TCanvas *c1 = new TCanvas("c1", "Linear Fit", 800, 600);
    // graph->Draw("AP");
    // linearFit->Draw("same");

    // Optionally save the plot and parameters to a file
    // c1->SaveAs("linear_fit.png");

    // Print out the fit parameters
    Line2D fitted_line(slope, intercept);
    
    // std::cout << " Slope (m): " << fitted_line.m() << ", Intercept (q): " << fitted_line.q() << std::endl;
    
    return fitted_line;
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
        Line l = RecoUtils::GetLineFromWire(digit);

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

Helix RecoUtils::GetHelixFromCircleLine(const Circle& circle, 
                                        const Line2D& line, 
                                        const Helix& true_helix,
                                        TVector3& momentum){
    /*
        input:
            circle : fitted trajectory in the ZY bending plane
            line : fitted trajectory in the XZ plane
            point x0 : trajectory stariting point (x0 <-> Phi0)
            momentum : passed as reference to be filled
        construct a helix from these inputs
        NOTE: a guess for the initial vertex is needed because related to the dip angle
    */
    auto vertex = true_helix.x0();

    auto helicity = true_helix.h();

    double pt = circle.R() * 0.29979 * 0.6;
    
    double pz = pt * (helicity * (vertex.Y() - circle.center_y()) / circle.R());
    
    double py = pt * (-helicity * (vertex.Z() - circle.center_x()) / circle.R());
    
    double pz_over_px = line.m();
    
    double px = pz / pz_over_px;
    
    double Phi0 = TMath::ATan2(py, pz) + helicity * TMath::Pi()*0.5;;

    momentum = {px, py, pz};

    double dip_angle = TMath::ATan2(px, pt);
    
    return Helix(circle.R(), dip_angle, Phi0, helicity, vertex);
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

std::vector<double> CardanoSolutions(double a, double b, double c, double d) {
    double f, g, h;
    std::complex<double> x1, x2, x3;

    // Calculate intermediate values
    f = ((3 * c / a) - ((b * b) / (a * a))) / 3;
    g = ((2 * (b * b * b) / (a * a * a)) - (9 * b * c / (a * a)) + (27 * d / a)) / 27;
    h = ((g * g) / 4) + ((f * f * f) / 27);
    std::cout << "h : " << h << "\n";
    if (f == 0 && g == 0 && h == 0) {
        // All roots are real and equal
        double root = pow(d / a, 1.0 / 3.0);
        std::cout << "x = " << root << std::endl;
    } else if (h <= 0) {
        // All 3 roots are real
        double i = sqrt(((g * g) / 4) - h);
        double j = pow(i, 1.0 / 3.0);
        double k = acos(-g / (2 * i));
        double l = -j;
        double m = cos(k / 3);
        double n = sqrt(3) * sin(k / 3);
        double p = -b / (3 * a);

        x1 = 2 * j * m - p;
        x2 = l * (m + n) + p;
        x3 = l * (m - n) + p;

        std::cout << "x1 = " << x1 << std::endl;
        std::cout << "x2 = " << x2 << std::endl;
        std::cout << "x3 = " << x3 << std::endl;
    } else if (h > 0) {
        // One real root and two complex roots
        double r = -g / 2 + sqrt(h);
        double s = pow(r, 1.0 / 3.0);
        double t = -g / 2 - sqrt(h);
        double u = pow(t, 1.0 / 3.0);

        x1 = s + u - b / (3 * a);
        x2 = std::complex<double>((s + u) / 2, -(s - u) * sqrt(3) / 2) - b / (3 * a);
        x3 = std::complex<double>((s + u) / 2, (s - u) * sqrt(3) / 2) - b / (3 * a);

        std::cout << "x1 = " << x1 << std::endl;
        std::cout << "x2 = " << x2 << std::endl;
        std::cout << "x3 = " << x3 << std::endl;
    }
    return {x1.real(), x2.real(), x3.real()};
}

std::vector<Line2D> RecoUtils::GetTangentsTo2Circles(const Circle& c1, const Circle& c2){
    /*
        Get the 4 tangents tangent to 2 separate circles;
    */
    std::vector<Line2D> tangents;
    double x1 = c1.center_x(), x2 = c2.center_x();
    double y1 = c1.center_y(), y2 = c2.center_y();
    double r1 = c1.R(), r2 = c2.R();
    double dx = x2 - x1;
    double dy = y2 - y1;
    double d = sqrt(dx*dx + dy*dy);
    for(auto i : {-1., 1.}){ // outer tangents, inner tangents
        double dr = r2 + i*r1;
        double gamma = i*atan2(dy, dx);
        double beta = asin(dr/d);
        for(auto k : {-1.,1.}){
            double alpha = gamma + k*beta;
            // tangent point to circle 1
            double x3 = x1 + k*r1*sin(alpha);
            double y3 = y1 - i*k*r1*cos(alpha);
            // tangent point to circle 2
            double x4 = x2 - i*k*r2*sin(alpha);
            double y4 = y2 + k*r2*cos(alpha);
            // direction, p0
            Line2D tangent({x4-x3, y4-y3}, {x3, y3});
            tangents.push_back(tangent);
        }
    }
    return tangents;
}

Line2D RecoUtils::GetTangetTo3Circles(const Circle& c1, const Circle& c2, const Circle& c3){
    // !! bug to be fixed
    // get 4 tangets to 2 circles
    Line2D tangent_3_circles;
    auto min_radius_difference = 999.;
    auto tangents = RecoUtils::GetTangentsTo2Circles(c1, c2);
    for(auto& t : tangents){
        double expected_radius = t.Distance2Point({c3.center_x(), c3.center_y()});
        double measured_radius = c3.R();
        double radius_difference = fabs(expected_radius - measured_radius);
        if(radius_difference < min_radius_difference){
            tangent_3_circles = t;
            min_radius_difference = radius_difference;
        }
    }
    return tangent_3_circles;
}

Line2D RecoUtils::GetBestTangent2NCircles(const std::vector<Circle>& circles){
  if(circles.size()<4){
    std::cout << "too few circles to fit in function " << __FILE__ << ", " << __LINE__ << "\n";
    std::cout << "try use RecoUtils::GetTangetTo3Circles \n";
    throw ""; 
  }
  auto guess_tangents = RecoUtils::GetTangentsTo2Circles(circles[0], circles[1]);
  double best_score = 999.;
  double sigma = 0.2;
  Line2D best_tangent;
  for(auto& tangent : guess_tangents){
    double tangent_score = 0.;
    for(auto i = 2u; i<circles.size(); i++){
      double expected_radius = tangent.Distance2Point({circles[i].center_x(), circles[i].center_y()});
      double measured_radius = circles[i].R();
      tangent_score += (expected_radius - measured_radius)*(expected_radius - measured_radius) / sigma*sigma;
    }
    tangent_score = sqrt(tangent_score)/(circles.size()-2);
    if(tangent_score < best_score){
      best_tangent = tangent;
      best_score = tangent_score;
    }
  }
  return best_tangent;
}
