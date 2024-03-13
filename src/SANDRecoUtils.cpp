#include "SANDRecoUtils.h"

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
        // x = t + x0
        // y = y0
        // z = z0
        Line l(1., 0., 0., digit.x, digit.y, digit.z);
        l.SetLineLength(digit.wire_length);
        return l;
    }else{
        // vertical == line along y axis
        // x = x0
        // y = t + y0
        // z = z0
        Line l(0., 1., 0., digit.x, digit.y, digit.z);
        l.SetLineLength(digit.wire_length);
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

    std::cout << "sinFit [x_min, x_max] : [" << x_min << ", " << x_max << "]\n";

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
    std::cout << "Fit parameters:" << std::endl;
    std::cout << "Amplitude: " << amplitude << std::endl;
    std::cout << "Frequency: " << frequency << std::endl;
    std::cout << "Phase: " << phase << std::endl;
    std::cout << "Offset: " << offset << std::endl;
 
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

TF1* RecoUtils::WiresCircleFit(const std::vector<dg_wire>& wires){
    /*
        Perform a 2D circular fit of the ZY coordinates
        of the horizontal wires 
    */
    std::vector<double> z,y;
    for(auto& wire : wires){
        z.push_back(wire.z);
        y.push_back(wire.y);
    }
    // Create a TGraph object and fill it with the data
    TGraph* graph = new TGraph(z.size(), &z[0], &y[0]);

    double z_min = *std::min_element(z.begin(), z.end());
    double z_max = *std::max_element(z.begin(), z.end());
    double y_min = *std::min_element(y.begin(), y.end());
    double y_max = *std::max_element(y.begin(), y.end());

    // Define a circle function
    TF1* circleFit = new TF1("circleFit", "[0] + sqrt([1]*[1] - (x-[2])*(x-[2]))", z_min-20, z_max+20); // Adjust range as needed

    double centerZ, centerY, radius, err, chi2;

    fitCircle(z.size(), z, y, centerZ, centerY, radius, err, chi2);
    // Imposta i parametri della funzione circolare
    circleFit->SetParameter(0, centerY);
    circleFit->SetParameter(1, radius);
    circleFit->SetParameter(2, centerZ);

    // Print the fitted parameters
    std::cout << "Fitted Circle Parameters:" << std::endl;
    std::cout << "Center (z, y): (" << centerZ << ", " << centerY << ")" << std::endl;
    std::cout << "Radius: " << radius << std::endl;
    std::cout << "err: " << err << std::endl;
    std::cout << "chi2: " << chi2 << std::endl;

    // Create a canvas to draw the graph and the fitted circle
    TCanvas *canvas = new TCanvas("canvas", "Fitted Circle", 800, 600);
    // gStyle->SetOptFit(1111); // Show fit results on the graph
    graph->Draw("AP");
    circleFit->Draw("same");

    // Save the canvas to a file if needed
    canvas->SaveAs("fitted_circle.png");

    return circleFit;
}

void  RecoUtils::WiresLinearFit(const std::vector<dg_wire>& wires){
    /*
        Perform a 2D linear fit of the XZ coordinates
        of the vertical wire 
    */
   std::vector<double> z,x;
   for(auto& wire : wires){
    if(wire.hor==false){ // vertical
        x.push_back(wire.x);
        z.push_back(wire.z);
    }
   }
   std::cout << "number of hits : " << z.size() << "\n";
    
    // Create a TGraph object and fill it with the data
    TGraph *graph = new TGraph(z.size(), &x[0], &z[0]);

    // Define a linear function with parameters p0 (slope) and p1 (intercept)
    TF1 *linearFit = new TF1("linearFit", "[0]*x + [1]", wires.front().z - 20, wires.front().z + 20); // Adjust xmin and xmax

    // Perform the fit using chi-square minimization method
    graph->Fit(linearFit, "Q"); // "Q" option suppresses the fit statistics output

    // Retrieve fit parameters
    double slope = linearFit->GetParameter(0);
    double intercept = linearFit->GetParameter(1);

    // Visualize the fit
    TCanvas *c1 = new TCanvas("c1", "Linear Fit", 800, 600);
    graph->Draw("AP");
    linearFit->Draw("same");

    // Print out the fit parameters
    std::cout << "Slope: " << slope << std::endl;
    std::cout << "Intercept: " << intercept << std::endl;

    // Optionally save the plot and parameters to a file
    c1->SaveAs("linear_fit.png");
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

void RecoUtils::InitHelixPars(const std::vector<dg_wire>& fired_wires,
                              Helix& helix_initial_guess)
{
    /*
        !!! ASSUMPTION !!! : 
            - WIRES CONFIGURATION XXY
            - WIRES ARE ORDERED
        constrain helix parameter first guess using fired wires
        - x0 form fit of vertical wires coordinates with a sin function 
    */
    std::vector<dg_wire> hor_wires, ver_wires;

    for (auto& wire : fired_wires)
    {
        if(wire.hor==true){ // horizontal
            hor_wires.push_back(wire);
        }else{// (wire.hor==false) // vertical
            ver_wires.push_back(wire);
        }
    }
    
    TF1* SinFit = RecoUtils::WiresSinFit(ver_wires);
    TF1* CircleFit = RecoUtils::WiresCircleFit(hor_wires);
    /*
        first gues of (x0,z0) is the value of the fitted function
        around the first fired wire
    */
    TVector3 x0 = {fired_wires.front().x, 
                   helix_initial_guess.x0().Y(), 
                   helix_initial_guess.x0().Z()};

    helix_initial_guess.Setx0(x0);
    helix_initial_guess.SetR(CircleFit->GetParameter(1));
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