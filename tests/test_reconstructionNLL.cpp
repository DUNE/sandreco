#include <iostream>
#include <fstream>
#include <TAxis.h>
#include <TRandom3.h>

#include "SANDRecoUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TG4Event.h"
#include "TH1D.h"

/*
    CODE SKELETON:
    - DEFAULT SIMULATION SETTINGS
    - CONSTANS
    - GLOBAL VARIABLES
    - INPUT/OUTPUT FUNCTIONS
    - EVENT SELECTION
    - DIGITIZATION
    - EDEP DIGITIZATION
    - FITTING
*/

// DEFAULT SIMULATION SETTINGS_________________________________________________

bool INCLUDE_SIGNAL_PROPAGATION = false;

bool INCLUDE_HIT_TIME = false;

bool INCLUDE_TDC_SMEARING = false;

bool _DEBUG_ = false;

bool USE_NON_SMEARED_TRACK = false;

int FITTING_STRATEGY = 0;
/*
    0 : fit TDCs on the XZ plane
        with an sin function and 
        on the ZY with a circle. 
        The fits are done separa
        tely.
    
    1 : fit TDCs simultaneously 
        on XZ and ZY plane using
        an Helix
*/ 

// CONSTANS____________________________________________________________________

const double SAND_CENTER_X = 0.;

const double SAND_CENTER_Y = -2384.73;

const double SAND_CENTER_Z = 23910.;

// !! HARDCODED (I know you are judging me, stop it)
// does not include frames
const double SAND_TRACKER_X_LENGTH = 3220.0;

// mm (approx), this includes 9 C3H6 mod + 1 C mod + clearances
const double SAND_SUPEMOD_Z_THICK = 364.6;

const double SAND_TRACKER_Z_START = 22952.14;

const std::vector<double> SUPERMOD_LENGTHS = {3129.97277395, // A1
                                              3550.97659024, // B1
                                              3755.16996258, // C1
                                              3755.16996258, // C2
                                              3550.97659024, // B2
                                              3129.97277395, // A2
                                              2466.05424949, // D
                                              1262.30218417};// F

const double DRIFT_CELL_DIM = 10.; // 1 cm

const double TDC_SMEARING = 1.; // ns, default disabled

const double NEUTRINO_INTERACTION_TIME = 1.; // ns

const double LIGHT_VELOCITY = 299.792458; // mm/ns

// GLOBAL VARIABLES____________________________________________________________

std::vector<dg_wire>* RecoUtils::event_digits = nullptr;

std::vector<dg_wire>* vertical_fired_digits = nullptr;

std::vector<dg_wire>* horizontal_fired_digits = nullptr;

dg_wire* first_fired_wire = nullptr;

std::vector<Line2D>* track_segments_ZY = nullptr;

std::vector<Line2D>* track_segments_XZ = nullptr;

TGeoManager* geo = nullptr;

TG4Event* evEdep = nullptr;

// COLORS ______________________________________________________________________

Color::Modifier def(Color::FG_DEFAULT);
Color::Modifier red(Color::FG_RED); // warnings
Color::Modifier green(Color::FG_GREEN); // opeartions
Color::Modifier blue(Color::FG_BLUE); // results

void LOG(TString i, const char* out){
    if(i.Contains("I")){ // info
        std::cout << green << "[INFO] " << out << def << std::endl;
        std::cout << "\n";
    }else if(i.Contains("i")){
        std::cout << green << " |___ " << out << def << std::endl;
        std::cout << "\n";
    }else if(i.Contains("R")){ // result
        std::cout << blue << "[RESULT] " << out << def << std::endl;
        std::cout << "\n";
    }else if(i.Contains("W")){ // waring
        std::cout << red << "[WARNING] " << out << def << std::endl;
        std::cout << "\n";
    }else{
        std::cout << out;
    }
}

// _____________________________________________________________________________

// INPUT/OUTPUT FUNCTIONS_______________________________________________________

void help_input(){
    std::cout << "\n";
    std::cout << red << "./buil/bin/test_reconstructionNLL "
              << "-edep <EDep file> "
            //   << "-digit <digitization file> "
              << "-wireinfo <WireInfo file> "
              << "-o <fOuptut.root> "
              << "[signal_propagation] [hit_time] [debug] [track_no_smear]\n";
    std::cout << "\n";
    std::cout << "<WireInfo file>      : see tests/wireinfo.txt \n";
    std::cout << "--signal_propagation : include signal_propagation in digitization \n";
    std::cout << "--hit_time           : include hit time in digitization \n";
    std::cout << "--track_no_smear     : reconstruct non smeared track (NO E_LOSS NO MCS)\n";
    std::cout << "--debug              : use higher verbosity for TMinuit " << def << std::endl;
}

void PrintEventInfos(const Helix& true_helix,
                     const Helix& helix_first_guess,
                     const std::vector<dg_wire>& fired_wires){
        LOG("", "muon true helix : \n");
        true_helix.PrintHelixPars();
        // const double true_pars[3] = {true_helix.Center().X(), true_helix.Center().Y(), true_helix.R(),};
        // double chi_squared = FunctorNLL_Circle(true_pars);
        // LOG("", TString::Format("Expected Chi2 from true helix : %f \n",chi_squared).Data());
}

void ReadWireInfos(const std::string fInput, std::vector<dg_wire>& wire_infos){
    std::ifstream stream(fInput);
    std::string line;
    std::string wire_value;
    int w_id, w_hor;
    double w_x, w_y, w_z, w_l;
    std::getline(stream, line); // needed to skip header
    while(std::getline(stream, line))
    {
        std::stringstream lineStream(line);
        std::vector<double> wire_values;
        while (getline(lineStream, wire_value, ','))
        {
            double value = std::stod(wire_value);
            wire_values.push_back(value);
        }
        dg_wire wire;
        wire.did = wire_values[0];
        wire.x = wire_values[1];
        wire.y = wire_values[2];
        wire.z = wire_values[3];
        wire.wire_length = wire_values[4];
        /*
            if on txt file you read hor(ientation) == 0,
            means horizontal. Here I invert the values
            to give to hor the mmeaning of horizontal
        */
        wire.hor =  (wire_values[5]==0) ? 1 : 0;
        wire_infos.push_back(wire);
    }
    stream.close();
}

// EVENT SELECTION_____________________________________________________________

bool PassSelectionNofHits(const std::vector<dg_wire>& fired_wires){
    /*
      Select events that have:
      - at least 3 hits on the banding plane (to make circular fit)
      - at least 2 hits on the XZ plane (to make linear fit)
    */
   int nof_ZY_hits = 0;
   int nof_XZ_hits = 0;
    if(fired_wires.size()<5){
        std::cout << "nof fired wires: " << fired_wires.size() << "\n";
        LOG("W", "skipping event");
        std::cout << "\n";
        return false;
    }else{
        for(auto& f : fired_wires){
            if(f.hor == true){
                nof_ZY_hits++;
            }else{
                nof_XZ_hits++;
            }
        }
        if((nof_ZY_hits<5)|(nof_XZ_hits<3)){
            std::cout << "ZY_hits : " << nof_ZY_hits << ", XZ_hits : " << nof_XZ_hits << "\n";
            LOG("W", "skipping event");
            return false;
        }else{
            return true;
        }
    };

}

bool IsInFiducialVolume(const TG4PrimaryVertex& vertex){
    /*
        select events in the sand fiducial volume at 10 cm far from the frames:
        - 5 cm along x from the frame edge;
    */
    auto vertex_position = vertex.GetPosition();
    bool pass_x = (fabs(vertex_position.X() - SAND_CENTER_X) <= SAND_TRACKER_X_LENGTH/2. - 50);
    // get supermod number
    int smod = (vertex_position.Z() - SAND_TRACKER_Z_START)/SAND_SUPEMOD_Z_THICK;
    // event should be 5 cm far from supermod frame
    bool pass_zy = (fabs(vertex_position.Y() - SAND_CENTER_Y) <= SUPERMOD_LENGTHS[smod]/2. - 50);
    return pass_x * pass_zy;
}

// DIGITIZATION________________________________________________________________

void CreateDigitsFromHelix(Helix& h, 
                           std::vector<dg_wire>& wire_infos, 
                           std::vector<dg_wire>& fired_wires){
    /* 
        - run over wires
        - set helix range
        - call RecoUtils::GetMinImpactParameter(const Helix& helix, const Line& line)
        - if it is < 5*sqrt(2) mm ( = half diagonal of a cell size of 10 mm)
          and the helix point is inside SAND tracker x range
            - fill dg_wire info with
            - add to event_wire_digits
    */
    fired_wires.clear();
    // usefull in case t_hit included in TDC
    double hit_timer = NEUTRINO_INTERACTION_TIME;
    double last_s_min = 0.;
    //
   // scan all the wires TDC and select only the onces with the closest impact par
   for(auto& w : wire_infos){

    Line l = RecoUtils::GetLineFromDigit(w);

    h.SetHelixRangeFromDigit(w);

    if(h.LowLim() > 0) continue;

    /* define point on the helix and on the line corresponding to the impact parameter
       t_min is how much do I move from the wire center to get to the line point closest to the helix.
       It has a sign.
    */
    double s_min, t_min;

    bool HasMinimized = false;

    double impact_par = RecoUtils::GetMinImpactParameter(h, l, s_min, t_min, HasMinimized);

    if((impact_par <= DRIFT_CELL_DIM*0.5*sqrt(2)) & (fabs(t_min)*l.GetDirectionVector().Mag() <= w.wire_length)){ // conditions for fired_wire

        w.drift_time = impact_par / sand_reco::stt::v_drift;

        w.tdc = w.drift_time;

        if(INCLUDE_TDC_SMEARING)
            w.tdc = RecoUtils::SmearVariable(w.drift_time, TDC_SMEARING, 1)[0];
        
        if(INCLUDE_SIGNAL_PROPAGATION){
            // ATTENCTION WE ARE ASSUMING t_min CALCULATE WRT TO WIRE CENTER
            // CHECK GetMinImpactParameter
            w.signal_time = (t_min + w.wire_length/2.) / sand_reco::stt::v_signal_inwire;
            w.tdc += w.signal_time;
        }

        if(INCLUDE_HIT_TIME){
            if(last_s_min!=0) // for first fired hit hit_timer=NEUTRINO_INTERACTION_TIME
                hit_timer += (h.GetPointAt(s_min) - h.GetPointAt(last_s_min)).Mag() / LIGHT_VELOCITY;
            w.t_hit = hit_timer;
            w.tdc += w.t_hit;
            last_s_min = s_min;
        }
        fired_wires.push_back(w);
    }
   }
   // oreder the wire hit by z and give the a hit time
//    SortWires(fired_wires);
}

std::vector<TG4HitSegment> FilterHits(const std::vector<TG4HitSegment>& hits, const int PDG){
    /*
    filter hits whose PrimaryId is equal to PDG
    */
    
    std::vector<TG4HitSegment> filtered_hits;

    for (auto& hit : hits)
    {
        bool contribution = false;
        for(auto& id : hit.Contrib){
            if((evEdep->Trajectories[id].GetPDGCode()==PDG)&&(evEdep->Trajectories[id].GetParentId()==-1)){
                contribution = true;
                break;
                }
        }
        if(contribution) filtered_hits.push_back(hit);
    }
    return filtered_hits;
}

double Hit2WireDistance(const dg_wire& wire,
                        const TG4HitSegment& hit,
                        double& closest2Hit,
                        double& closest2Wire){
                        // std::ofstream& fout){
    /*
        given a hit and a wire, find the 3D distance
        and fill:
        - closest2Hit: the point on wire closest to the hit
        - closest2Wire: the point on hit line closest to the wire
    */
    Line wire_line = RecoUtils::GetLineFromDigit(wire);
    Line hit_line = Line(hit);
    
    double d = RecoUtils::GetSegmentSegmentDistance(wire_line, hit_line, closest2Hit, closest2Wire);
    auto w_point = wire_line.GetPointAt(closest2Hit);
    auto h_point = hit_line.GetPointAt(closest2Wire);
    // if(d<=DRIFT_CELL_DIM * 0.5 * sqrt(2)){
    //     fout << evEdep->EventId <<","<< hit.GetStart().X() << "," << hit.GetStart().Y() << "," << hit.GetStart().Z()<< ",";
    //                             fout << hit.GetStop().X()  << "," << hit.GetStop().Y()  << "," << hit.GetStop().Z()<< ",";
    //                             fout << w_point.X()  << "," << w_point.Y()  << "," << w_point.Z()<< ",";
    //                             fout << h_point.X()  << "," << h_point.Y()  << "," << h_point.Z() << std::endl;
    // }
    return d;
}

void UpdateFiredWires(std::vector<dg_wire>& fired_wires, dg_wire new_fired){
    
    bool found = false;
    //(f.did == new_fired.did) && (f.tdc > new_fired.tdc)
    for(auto& f : fired_wires){
        if((f.did == new_fired.did)){ // wire already fired
            found = true;
            // if tdc is shorter than the one saved, update otherwise discard wire
            if(f.tdc > new_fired.tdc){
                f.tdc = new_fired.tdc;
                f.drift_time = new_fired.drift_time;
                f.signal_time = new_fired.signal_time;
                f.t_hit = new_fired.t_hit;
            }
        }
    }
    if(!found) fired_wires.push_back(new_fired);
}

// EDEP DIGITIZATION___________________________________________________________

void SortWiresByTime(std::vector<dg_wire>& wires){
    // sort fired wires by hit time (MC truth): usefull to subract assumed t_hit from measured TDC
    std::sort(wires.begin(), wires.end(), [](const dg_wire &w1, const dg_wire &w2) {
    return w1.t_hit < w2.t_hit;
    });
}

void CreateDigitsFromEDep(const std::vector<TG4HitSegment>& hits,
                          const std::vector<dg_wire>& wire_infos, 
                          std::vector<dg_wire>& fired_wires){
    /*
    Perform digitization of edepsim hits
    */
    fired_wires.clear();
    
    // filter only muon hit segments
    LOG("i", "Select muon hits from MC truth");
    auto muon_hits = FilterHits(hits, 13);

    for(auto& hit : muon_hits){
        
        for(auto& wire : wire_infos){

            auto hit_middle = (hit.GetStart()+hit.GetStop())*0.5;
            
            if(fabs(wire.z - hit_middle.Z()) <= DRIFT_CELL_DIM * 0.5){//first scan along z to find the wire plane

                double closest2Hit, closest2Wire;
                double hit2wire_dist = Hit2WireDistance(wire, hit, closest2Hit, closest2Wire);

                bool pass = (wire.hor==true) ? 
                    (hit2wire_dist <= DRIFT_CELL_DIM * 0.5 * sqrt(2.)) : (hit2wire_dist <= DRIFT_CELL_DIM * sqrt(2.)); // due to an error on the cell dimension of vertical planes
                
                if(pass){// find closest wire on the plane
                    auto fired = RecoUtils::Copy(wire);
                    fired.drift_time = hit2wire_dist / sand_reco::stt::v_drift;
                    fired.signal_time = fabs(closest2Hit - fired.wire_length*0.5) / sand_reco::stt::v_signal_inwire;
                    fired.t_hit = ((hit.GetStop() + hit.GetStart())*0.5 + (hit.GetStop() - hit.GetStart())*closest2Wire).T();
                    fired.tdc = fired.drift_time + fired.signal_time +fired.t_hit;
                    UpdateFiredWires(fired_wires, fired);
                }
            }
        }
    }
}

// FITTING_____________________________________________________________________

Parameter CreateParam(std::string name, 
                      int id,
                      bool fixed_in_fit,
                      double initial_guess,
                      double value,
                      double error){
    Parameter p;
    p.name = name;
    p.id = id;
    p.fixed_in_fit = fixed_in_fit;
    p.initial_guess = initial_guess;
    p.value = value;
    p.error = error;
    return p;
}

Helix GetHelixFirstGuess(const Helix& input_helix){
    Helix first_guess(RecoUtils::SmearVariable(input_helix.R(), 100, 1)[0], // smear of 10 cm around true R
                      RecoUtils::SmearVariable(input_helix.dip(), 0.1, 1)[0], 
                      input_helix.Phi0(), 
                       1, 
                      {RecoUtils::SmearVariable(input_helix.x0().X(), 10, 1)[0],
                       input_helix.x0().Y(),
                       input_helix.x0().Z()});
    return first_guess;
}

double WireZdistance(const dg_wire& wire1, const dg_wire& wire2) {
    return std::fabs(wire1.z - wire2.z);
}

dg_wire FindClosestZWire(const dg_wire& input_wire, 
                                       const std::vector<dg_wire>& wires){
    /*
        return wire with opposit orientation 
        closest to the input_wire along
        the Z axis.
    */
    dg_wire closest;
    double minDistance = std::numeric_limits<double>::max();
    for(const auto& w : wires){
        // compute input_wire - w distance
        double dist = WireZdistance(input_wire, w);
        if(dist < minDistance){
            // closest wire on the right
            minDistance = dist;
            closest = w;
        }
    }
    return closest;
}

void TDC2ImpactPar(dg_wire& wire){
    /*
        Covert the measure TDC into a impact
        parameter (distance track - wire) that
        has produced the measured TDC.
    */
   wire.t_hit_measured = wire.t_hit;

   if(INCLUDE_HIT_TIME){
     wire.t_hit_measured = wire.t_hit; // use MC truth for now
   }

   if(INCLUDE_SIGNAL_PROPAGATION){
    /*
        the coordinate on the (horizontal / vertical) wire from 
        which the signal propagates is approximatively given by  
        the closest (vertical / horizontal) fired wire in space.
    */
        if(wire.hor==true){ // horizontal
            // find closest vertical wire along z
            auto ClosestFiredWire = FindClosestZWire(wire, *vertical_fired_digits);
            
            wire.signal_time_measured = fabs(ClosestFiredWire.x - wire.wire_length*0.5) / sand_reco::stt::v_signal_inwire;
        }else{ // vertical
            // find closest horizontal wires along z
            auto ClosestFiredWire = FindClosestZWire(wire, *horizontal_fired_digits);

            wire.signal_time_measured = fabs(ClosestFiredWire.y + wire.wire_length*0.5) / sand_reco::stt::v_signal_inwire;
        }
   }
   wire.drift_time_measured = wire.tdc - wire.signal_time_measured - wire.t_hit_measured;
}

void TDC2DriftDistance(std::vector<dg_wire>& fired_wires){
    /*
        Convert the observed TDC into a drift distance
        for each fired_wires. Store the infos in the wire
        attributes: t_hit_measured, signal_time_measured,
        drift_time_measured
    */
   for (auto& wire : fired_wires)
   {
        TDC2ImpactPar(wire);
        // std::cout << "wire id " << wire.did << " is horizontal : " << wire.hor <<"\n";
        // LOG("R", TString::Format("drift time: true %f, measured %f", wire.drift_time, wire.drift_time_measured).Data());
        // std::cout << "\n";
   }
}

std::vector<Circle> Wire2DriftCircle(const std::vector<dg_wire>& wires){
    std::vector<Circle> drift_cirlces;
    for(const auto& wire : wires){
        if(wire.hor == true){ // horizontal
            Circle c(wire.z, wire.y, wire.drift_time_measured * sand_reco::stt::v_drift);
            drift_cirlces.push_back(c);
        }else{ // vertical
            Circle c(wire.x, wire.z, wire.drift_time_measured * sand_reco::stt::v_drift);
            drift_cirlces.push_back(c);
        }
    }
    return drift_cirlces;
}

// FITTING STRATEGY 0__________________________________________________________

// SEGMENT FITTING_____________________________________________________________

std::vector<Line2D> GetTrackSegments(const std::vector<Circle>& drift_circles){
    std::vector<Line2D> track_segments;
    unsigned int nof_circles = drift_circles.size();
    for(auto i=0u; i < nof_circles - 2; i++){
        const Circle& first = drift_circles[i];
        const Circle& second = drift_circles[i+1];
        const Circle& third = drift_circles[i+2];
        auto tangent12 = RecoUtils::GetTangetTo3Circles(first, second, third);
        track_segments.push_back(tangent12);
        // if(i==0){
        //     first.PrintCircleInfo();
        //     second.PrintCircleInfo();
        //     third.PrintCircleInfo();
        // }
    }
    return track_segments;
}

// FITTING ON ZY PLANE_________________________________________________________

// CIRLCE FITTING______________________________________________________________

double GetImpactParamiterCircle(Circle c, const dg_wire& wire){
    TVector2 wire_center = {wire.z, wire.y};
    return  c.Distance2Point(wire_center);
}

double FunctorNLL_Circle(const double* p){
    
    double zc = p[0];
    double yc = p[1];
    double R = p[2];

    double nll         = 0.;
    const double sigma = 0.2; // 200 mu_m = 0.2 mm

    Circle c(zc, yc, R);

    // circle has to be "perfectly" tangent to the first 2 circles
    auto first_point = (*track_segments_ZY)[0].p0();
    auto second_point = (*track_segments_ZY)[1].p0();
    auto dist2first = c.Distance2Point({first_point.X(), first_point.Y()});
    auto dist2second = c.Distance2Point({second_point.X(), second_point.Y()});
    double i = 1;
    for(auto& wire : *horizontal_fired_digits){

        // r_estimated = impact parameter estimated as distance sin function - wire
        double r_estimated = GetImpactParamiterCircle(c, wire);

        // r_observed = impact parameter from the observed tdc
        double r_observed = wire.drift_time_measured * sand_reco::stt::v_drift;

        nll += (r_estimated - r_observed) * (r_estimated - r_observed) / (sigma * sigma);

        i++;
    }
    double nll_epsilon = (dist2first*dist2second)/(sigma*sigma);
    return sqrt(nll + nll_epsilon)/(horizontal_fired_digits->size());
}

Circle GetRecoZYTrack(Circle& first_guess,
                      MinuitFitInfos& fit_infos){
    /*
        Fit observed TDCs on ZY plane with a Circle.
        Parameter of the circle:
        - (zc, yc) center of the fitted circle
        - R radius of the circle
    */
    double zc = first_guess.center_x();
    double yc = first_guess.center_y();
    double R = first_guess.R();
 
    ROOT::Math::Functor functor(&FunctorNLL_Circle, 3);
    
    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");

    minimizer->SetFunction(functor);

    // zc
    minimizer->SetVariable(0, "zc", zc, 5);
    // yc
    minimizer->SetVariable(1, "yc", yc, 5);
    // R
    minimizer->SetVariable(2, "R", R, 200); // 200 mm -> 35 MeV difference

    // minimization settings
    if(_DEBUG_) minimizer->SetPrintLevel(4);

    // start minimization
    minimizer->Minimize();

    // retrieve result of the minimization and errors
    const double* pars = minimizer->X();
    const double* parsErrors = minimizer->Errors();
    
    // fill output object fit_infos
    Parameter center_z = CreateParam("zc", 0, false, zc, pars[0], parsErrors[0]);
    Parameter center_y = CreateParam("yc", 1, false, yc, pars[1], parsErrors[1]);
    Parameter radius   = CreateParam("R", 2, false, R, pars[2], parsErrors[2]);
    fit_infos.Auxiliary_name = "Circular_fit_ZY";
    fit_infos.fitted_parameters = {center_z, center_y, radius};
    fit_infos.TMinuitFinalStatus = minimizer->Status();
    fit_infos.NIterations = minimizer->NIterations();
    fit_infos.MinValue = minimizer->MinValue();

    // track from final fit
    Circle c_reco(center_z.value, center_y.value, radius.value);

    // print results of the minimization
    minimizer->PrintResults();

    return c_reco;
}

// FITTING ON XZ PLANE_________________________________________________________

// linear fit
double FunctorNLL_Line(const double* p){
    double m = p[0];
    double q = p[1];

    Line2D test_line(m, q);
    
    double nll = 0.;
    double sigma = 0.2;

    // circle has to be "perfectly" tangent to the first 2 circles
    // auto first_point = (*track_segments_XZ)[0].p0();
    // auto second_point = (*track_segments_XZ)[1].p0();
    // auto dist2first = test_line.Distance2Point({first_point.X(), first_point.Y()});
    // auto dist2second = test_line.Distance2Point({second_point.X(), second_point.Y()});
    // std::cout << "-----------------------------------------------------------------\n";
    double dist2first = 0.;
    double dist2second = 0.;
    double i = 1;
    for(auto& wire : *vertical_fired_digits){

        // r_estimated = impact parameter estimated as 2D distance line - wire
        double r_estimated = test_line.Distance2Point({wire.x, wire.z});

        // r_observed = impact parameter from the observed tdc
        double r_observed = wire.drift_time_measured * sand_reco::stt::v_drift;

        nll += (r_estimated - r_observed) * (r_estimated - r_observed) / (sigma * sigma);

        // std::cout << "wire id : " << wire.did
        //           << ", impact par : true " << wire.drift_time * sand_reco::stt::v_drift
        //           << ", r_observed : " << r_observed
        //           << ", r_estimated : " << r_estimated << "\n";

        i++;
    }
    double nll_epsilon = (dist2first*dist2second)/(sigma*sigma);
    return sqrt(nll+nll_epsilon)/(horizontal_fired_digits->size());
    // return sqrt(nll)/(horizontal_fired_digits->size());
}

Line2D GetLineFitXZPlane(Line2D first_guess,
                       MinuitFitInfos& fit_infos){
    /*
        Fit observed TDCs on XZ plane with a Line.
        Parameter of the Line:
        y = m*x + q
        - m : slope
        - q : intercept
    */
    ROOT::Math::Functor functor(&FunctorNLL_Line, 2);

    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");

    minimizer->SetFunction(functor);

    minimizer->SetVariable(0, "m", first_guess.m(), 0.001);
    
    minimizer->SetVariable(1, "q", first_guess.q(), 1.);

    // minimization settings
    if(_DEBUG_) minimizer->SetPrintLevel(4);

    minimizer->Minimize();

    // retrieve result of the minimization and errors
    const double* pars = minimizer->X();
    const double* parsErrors = minimizer->Errors();

    // fill output object fit_infos
    Parameter m = CreateParam("m", 0, false, first_guess.m(), pars[0], parsErrors[0]);
    Parameter q = CreateParam("q", 1, false, first_guess.q(), pars[1], parsErrors[1]);
    fit_infos.Auxiliary_name = "Linear_fit_XZ";
    fit_infos.fitted_parameters = {m, q};
    fit_infos.TMinuitFinalStatus = minimizer->Status();
    fit_infos.NIterations = minimizer->NIterations();
    fit_infos.MinValue = minimizer->MinValue();
    
    // print results of the minimization
    minimizer->PrintResults();

    Line2D l(m.value, q.value);

    return l;
}

// sin fit

double GetImpactParameterSin(TF1* TestSin, const dg_wire& wire){
    /*
        Use numerical Newton method to find distance sin function to a point 
        in (wire.x, wire.z) in 2D space. TestSin = A sin(Bx + C) + D, functi
        on to be minimized is the euclidian distance: 
        - D(x) = (x_wire - x)**2 + (z_wire - TestSin(x))**2
        - D_prime(x) = 2(x_wire - x) + 2(z_wire - TestSin(x))*TestSin_prime(x)
        - TestSin_prime(x) = A B cos(Bx + C)
    */
    double x_wire = wire.x;
    double z_wire = wire.z;
    double x_min = TestSin->GetXmin();
    double x_max = TestSin->GetXmax();

    TF1* TestSin_prime = new TF1("TestSin_prime", "[0]*[1]*cos([1]*x + [2])", x_min, x_max);

    TestSin_prime->SetParameters(TestSin->GetParameter(0),
                                 TestSin->GetParameter(1),
                                 TestSin->GetParameter(2));
    /*
        distance_sin_wire : distance curve TestSin to a point (x_wire, z_wire)
        - variable x (space coordinate)
        - parameter (x_wire, z_wire) coordinate f the wire to which the distace
          is calculated
    */
    TF1* distance_sin_wire = new TF1("distance_sin_wire", [&](double*x, double *par){
        double sin_at_x = TestSin->Eval(x[0]);
        double x_w = par[0];
        double z_w = par[1];
        return sqrt((x[0] - x_w)*(x[0] - x_w) + (sin_at_x - z_w)*(sin_at_x - z_w));
    }, x_min, x_max, 2);              
    distance_sin_wire->SetParameters(x_wire, z_wire);
    
    // derivative of distance
    TF1* distance_sin_wire_prime = new TF1("distance_sin_wire_prime", [&](double*x, double *par){
        double distance = distance_sin_wire->Eval(x[0]);
        double sin_at_x = TestSin->Eval(x[0]);
        double derivative_at_x = TestSin_prime->Eval(x[0]);
        double x_w = par[0];
        double z_w = par[1];
        return 1./(2.*distance)*(2*(x[0] - x_w) + 2*(sin_at_x - z_w)*derivative_at_x);
    }, x_min, x_max, 2);
    distance_sin_wire_prime->SetParameters(x_wire, z_wire);

    double x_guess = x_wire;

    double impact_par = RecoUtils::NewtonRaphson2D(distance_sin_wire, distance_sin_wire_prime, x_guess, 0.2, 1000);

    return impact_par;
}

double NLL_Sin(TF1* TestSin, std::vector<dg_wire>& wires){
    /*
        run over vertical fired wires
        - get impact par from TDC
        - compare it with wire - track 
          min distance
    */
    static int calls   = 0;
    double nll         = 0.;
    const double sigma = 0.2; // 200 mu_m = 0.2 mm
    
    for(auto& wire : wires){

        // r_estimated = impact parameter estimated as distance sin function - wire
        double r_estimated = GetImpactParameterSin(TestSin, wire);

        // r_observed = impact parameter from the observed tdc
        double r_observed = wire.drift_time_measured * sand_reco::stt::v_drift;

        nll += (r_estimated - r_observed) * (r_estimated - r_observed) / (sigma * sigma);

        // std::cout << "impact par : true " << wire.drift_time * sand_reco::stt::v_drift
        //                << ", r_observed : " << r_observed
        //                << ", r_estimated : " << r_estimated << "\n";

    }
    // throw "";
    calls ++;
    return sqrt(nll)/wires.size();
}

double FunctorNLL_Sin(const double* p){
    
    double amplitude = p[0];
    double frequency = p[1];
    double phase = p[2];
    double offset = p[3];
    // fixed
    double x_min = p[4];
    double x_max = p[5];
    
    TF1* TestSin = new TF1("TestSin", "[0]*sin([1]*x + [2]) + [3]", x_min, x_max);
    
    // Set initial parameters for the fit
    TestSin->SetParameters(amplitude, 
                           frequency, 
                           phase, 
                           offset);
    
    return NLL_Sin(TestSin, *vertical_fired_digits);
}

TF1* GetRecoXZTrack(TF1* first_guess,
                    MinuitFitInfos& fit_infos){
    /*
        Fit TDCs on XZ plane using a sin function.
        Parameter of the sin function:
        A * sin(B * x + C) + D
        - A: amplitude 
        - B: frequency
        - C: phase
        - D: offset
    */
    double amplitude = first_guess->GetParameter(0); // A
    double frequency = first_guess->GetParameter(1); // B
    double phase = first_guess->GetParameter(2); // C
    double offset = first_guess->GetParameter(3); // D
    double x_min = first_guess->GetXmin();
    double x_max = first_guess->GetXmax();
    
    ROOT::Math::Functor functor(&FunctorNLL_Sin, 6);                  

    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");

    minimizer->SetFunction(functor);

    // A
    minimizer->SetLimitedVariable(0, "A", amplitude, amplitude*0.01, amplitude*0.8,amplitude*1.2);
    // B
    minimizer->SetLimitedVariable(1, "B", frequency, frequency*0.01, frequency*0.8, frequency*1.2);
    // C
    // minimizer->SetLimitedVariable(2, "C", phase, phase*0.01, phase*0.8, phase*1.2);
    // D
    // minimizer->SetLimitedVariable(3, "D", offset, offset*0.01, offset*0.8, offset*1.2);

    // extreme [x_min, x_max]
    minimizer->SetFixedVariable(2, "C", phase);
    minimizer->SetFixedVariable(3, "D", offset);
    minimizer->SetFixedVariable(4, "x_min", x_min);
    minimizer->SetFixedVariable(5, "x_max", x_max);

    // minimization settings
    if(_DEBUG_) minimizer->SetPrintLevel(4);

    // start minimization
    minimizer->Minimize();

    // retrieve result of the minimization and errors
    const double* pars = minimizer->X();
    const double* parsErrors = minimizer->Errors();
    
    // fill output object fit_infos
    Parameter A = CreateParam("amplitude", 0, false, amplitude, pars[0], parsErrors[0]);                                        
    Parameter B = CreateParam("frequency", 1, false, frequency, pars[1], parsErrors[1]);                                        
    Parameter C = CreateParam("phase", 2, true, phase, pars[2], parsErrors[2]);                                        
    Parameter D = CreateParam("offset", 3, true, offset, pars[3], parsErrors[3]);
    fit_infos.Auxiliary_name = "Sin_Fit_XZ";                          
    fit_infos.fitted_parameters = {A,B,C,D};                          
    fit_infos.TMinuitFinalStatus = minimizer->Status();
    fit_infos.NIterations = minimizer->NIterations();
    fit_infos.MinValue = minimizer->MinValue();

    // final fit
    TF1* SinFinalFit = new TF1("SinFinalFit", "[0]*sin([1]*x + [2]) + [3]", x_min, x_max);
    SinFinalFit->SetParameters(pars[0],pars[1],pars[2],pars[3]);

    // print results of the minimization
    minimizer->PrintResults();

    return SinFinalFit;
}

// HELIX FITTING (FITTING STRATEGY 1)__________________________________________

double NLL(Helix& h, std::vector<dg_wire>& wires){
    
    static int calls   = 0;
    double nll         = 0.;
    const double sigma = 0.2; // 200 mu_m = 0.2 mm
    for (auto& wire : wires)
    {
        Line l = RecoUtils::GetLineFromDigit(wire);
        
        h.SetHelixRangeFromDigit(wire);
        
        // r_estimated = impact parameter estimated as distance helix - wire
        double r_estimated = RecoUtils::GetMinImpactParameter(h,l);

        // r_observed = impact parameter from the observed tdc
        double r_observed = wire.drift_time_measured * sand_reco::stt::v_drift;
        
        nll += (r_estimated - r_observed) * (r_estimated - r_observed) / (sigma * sigma);
    }
    calls ++;
    return sqrt(nll)/wires.size();
}

double FunctorNLL(const double* p){
    
    double R    = p[0];
    double dip  = p[1];
    double Phi0 = p[2];
    int hel     = p[3];
    TVector3 x0 = {p[4],p[5],p[6]};

    Helix h(R, dip, Phi0, hel, x0);

    return NLL(h, *RecoUtils::event_digits);
}

Helix GetRecoHelix(Helix& helix_initial_guess,
                   std::vector<dg_wire>& fired_wires, 
                   MinuitFitInfos& fit_infos,
                   int FitStrategy = 1){
    /*
        Fit measured TDCs from fired wires assuming
        helicoidal track.
        Helix is described by the parameters :
        - R: helix radius
        - dip: dip angle
        - x0: TVector3 vertex of the track
        - Phi0: angle between x0 and the heli center
    */

    ROOT::Math::Functor functor(&FunctorNLL, 7);

    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");

    minimizer->SetFunction(functor);
    
    minimizer->SetLimitedVariable(0,    "R",    helix_initial_guess.R(), 
                                        10, helix_initial_guess.R()*0.8, helix_initial_guess.R()*1.2); // R fit range
    
    minimizer->SetLimitedVariable(1,    "dip",  helix_initial_guess.dip(), 
                                        0.01, helix_initial_guess.dip()*0.8, helix_initial_guess.dip()*1.2); // dip fit range
    
    minimizer->SetLimitedVariable(4,    "x0_x", helix_initial_guess.x0().X(), 
                                        1, helix_initial_guess.x0().X() - 20, helix_initial_guess.x0().X() + 20); // x0 fit range
    
    minimizer->SetFixedVariable(2,    "Phi0",    helix_initial_guess.Phi0());
    minimizer->SetFixedVariable(3,    "h",       helix_initial_guess.h());
    minimizer->SetFixedVariable(5,    "x0_y",    helix_initial_guess.x0().Y());
    minimizer->SetFixedVariable(6,    "x0_z",    helix_initial_guess.x0().Z());

    // minimization settings
    if(_DEBUG_) minimizer->SetPrintLevel(4);

    if(FitStrategy==1){
    }else{
        // relax default tolerance
        minimizer->SetTolerance(minimizer->Tolerance()*10);
        // change strategy
        minimizer->SetStrategy(2);
        std::cout<< "------------------------------\n";
        std::cout<< "setting fitting strategy to 2 \n";
    }

    // start minimization
    minimizer->Minimize();
    
    // retrieve result of the minimization and errors
    const double* pars = minimizer->X();
    const double* parsErrors = minimizer->Errors();

    // print results of the minimization
    minimizer->PrintResults();
    
    if(minimizer->Status()){ // failed fit
        if(minimizer->Strategy()==2){
            std::cout << "failed both with stratefy 1 and 2 \n";
        }else{
            // try to fit changing fit strategy
            GetRecoHelix(helix_initial_guess, fired_wires, fit_infos, 2);
        }
    }
    
    // fill output object fit_infos
    Parameter R_    = CreateParam("R", 0, false, helix_initial_guess.R(), pars[0], parsErrors[0]);
    Parameter dip_  = CreateParam("dip", 1, false, helix_initial_guess.dip(), pars[1], parsErrors[1]);
    Parameter Phi0_ = CreateParam("Phi0", 2, true, helix_initial_guess.Phi0(), pars[2], parsErrors[2]);
    Parameter h_    = CreateParam("h", 3, true, helix_initial_guess.h(), pars[3], parsErrors[3]);
    Parameter x0_   = CreateParam("x0", 4, false, helix_initial_guess.x0().X(), pars[4], parsErrors[4]);
    Parameter y0_   = CreateParam("y0", 5, true, helix_initial_guess.x0().Y(), pars[5], parsErrors[5]);
    Parameter z0_   = CreateParam("z0", 6, true, helix_initial_guess.x0().Z(), pars[6], parsErrors[6]);
    fit_infos.Auxiliary_name = "Helix_fit";
    fit_infos.fitted_parameters = {R_,dip_,Phi0_,h_,x0_,y0_,z0_};
    fit_infos.TMinuitFinalStatus = minimizer->Status();
    fit_infos.NIterations = minimizer->NIterations();
    fit_infos.MinValue = minimizer->MinValue();

    // create reconstructed helix from fitted parameters
    Helix reco_helix(pars[0], pars[1], pars[2], pars[3], {pars[4], pars[5], pars[6]});

    return reco_helix;
}

Helix Reconstruct(Circle FittedCircle,
                  Line2D FittedLine,
                  const std::vector<dg_wire>& hor_wires,
                  const std::vector<dg_wire>& ver_wires
                  ){
    /*
        Define the reconstructed helix from
        the separate reconstruction of the 
        track in the XZ plane and ZY plane
    */
    double zc = FittedCircle.center_x();
    double yc = FittedCircle.center_y();
    double R = FittedCircle.R();
    double m = FittedLine.m();
    double q = FittedLine.q();
    /*
        vertex of the track on the ZY plane is
        given by the point (z0,y0) on the fitted
        circle where z0 1 cm far from the first wire
        and y0=y0(z0).
    */
    bool forward_track = (hor_wires[0].z < hor_wires[1].z) ? true : false;
    auto first_point_XZ = (*track_segments_XZ)[0];
    TVector2 track_direction_XZ = first_point_XZ.direction() * (1./first_point_XZ.direction().Mod());
    
    // move 1 cm backwards (or forwards) in the direction of the line from the line starting point
    TVector2 vertex_XZ = (forward_track) ? first_point_XZ.p0() - track_direction_XZ * 10. : first_point_XZ.p0() + track_direction_XZ * 10.;
    // std::cout << "first point on XZ plane : (" << first_point_XZ.p0().X() << ", " << first_point_XZ.p0().Y() << ")\n";
    // std::cout << "x0, z0  = (" << vertex_XZ.X() << "," << vertex_XZ.Y() << ")\n";

    // reconstruct vertex yz starting from the tangent to the first 3 circles
    auto first_point_ZY = (*track_segments_ZY)[0];
    auto track_direction_ZY = first_point_ZY.direction() * (1./first_point_ZY.direction().Mod());
    TVector2 vertex_ZY = (forward_track) ? first_point_ZY.p0() - track_direction_ZY * 10. : first_point_ZY.p0() + track_direction_ZY * 10.;
    // std::cout << "z0, y0  = (" << vertex_ZY.X() << "," << vertex_ZY.Y() << ")\n";

    double Phi0 = FittedCircle.GetAngleFromPoint(vertex_ZY.X(), vertex_ZY.Y());
    TVector2 pt = R*0.3*0.6 * track_direction_ZY;
    /*
        reconstructed pt is the value of the tangent0 
        to the circle in the first point of the curve
    */
    double pz = pt.X();
    double py = pt.Y();
    double px = pz / FittedLine.m();
    // std::cout << "pz_reco : " << pz << "\n";
    // std::cout << "pz_reco : " << py << "\n";
    // std::cout << "px_reco : " << px << "\n";
    
    double dip = atan2(px, pt.Mod());
    
    Helix h(R, dip, Phi0, 1, {vertex_XZ.X(), vertex_ZY.Y(), vertex_XZ.Y()});
    
    return h;
}

Helix Reconstruct(MinuitFitInfos& fit_infos,
                  Helix& true_helix,
                  Helix& helix_first_guess,
                  std::vector<dg_wire>& wire_infos,
                  std::vector<dg_wire>& fired_wires){

    Helix reco_helix;

    reco_helix = GetRecoHelix(helix_first_guess, fired_wires, fit_infos);

    std::cout<<"\n";
    std::cout<< "true helix - reco helix \n";
    std::cout<<"R_true - R_reco     : "<< true_helix.R() - reco_helix.R() << "\n";
    std::cout<<"dip_true - dip_reco : "<< true_helix.dip() - reco_helix.dip() << "\n";
    std::cout<<"x0_true - x0_reco   : "<< true_helix.x0().X() - reco_helix.x0().X() << "\n";

    return reco_helix;
}

// MAIN________________________________________________________________________

int main(int argc, char* argv[]){

    if (argc < 3 || argc > 11) {
        help_input();
    return -1;
    }
    /* 
       provide edepsim input to run the reconstruction on non-smeared muon tracks
       (for testing used file "/storage/gpfs_data/neutrino/users/gi/sand-reco/tests/test_drift_edep.root")
    */
    const char* fEDepInput;
    
    /*
       provide digit file input to run the reconstruction on edep smeared muon tracks
    */
    const char* fDigitInput;
    
    const char* fOutput;
    
    /* 
       look up table previously created in the Digitization as csv file that contains
       all info about wires in the geometry
       for testing used file  "/storage/gpfs_data/neutrino/users/gi/sand-reco/tests/wireinfo.txt"
    */
    const char* fWireInfo;

    int index = 1;
    
    LOG("","\n");
    LOG("I","Reading inputs");

    while (index < argc)
    {
        TString opt = argv[index];
        if(opt.CompareTo("-edep")==0){
            try
            {
                fEDepInput = argv[++index];
                std::cout << "edep file       " << fEDepInput << std::endl;
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() << '\n';
                return 1;
            }
        }else if(opt.CompareTo("-wireinfo")==0){
            try
            {
                fWireInfo = argv[++index];
                std::cout << "wire_info file  " << fWireInfo << std::endl;
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() << '\n';
                return 1;
            }
        }else if(opt.CompareTo("-o")==0){
            try
            {
                fOutput = argv[++index];
                std::cout << "Output file     " << fOutput << std::endl;
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() << '\n';
                return 1;
            }
        }else if(opt.CompareTo("--signal_propagation")==0){
            INCLUDE_SIGNAL_PROPAGATION = true;
        }else if(opt.CompareTo("--hit_time")==0){
            INCLUDE_HIT_TIME = true;
        }else if(opt.CompareTo("--track_no_smear")==0){
            USE_NON_SMEARED_TRACK = true;
        }else if(opt.CompareTo("--debug")==0){
            _DEBUG_ = true;
        }
        else{
            auto ui = argv[++index];
            std::cout<<"unknown input "<< ui << "\n";
            return 1; 
        }
        index++;
    }

    std::cout << "Signal propagation..." << (INCLUDE_SIGNAL_PROPAGATION ? "enabled" : "disabled") << std::endl;
    std::cout << "Hit time............." << (INCLUDE_HIT_TIME ? "enabled" : "disabled") << std::endl;
    std::cout << "Debug mode..........." << (_DEBUG_ ? "enabled" : "disabled") << std::endl;
    std::cout << "\n";
    
    TFile fEDep(fEDepInput, "READ");
    
    // TFile fDigit(fDigitInput, "READ");
    
    TFile fout(fOutput, "RECREATE");
    
    TTree* tEdep = (TTree*)fEDep.Get("EDepSimEvents");
    
    // TTree* tDigit = (TTree*)fDigit.Get("tDigit");
    
    TTree tout("tReco", "tReco");
    
    // EventReco event_reco;
    RecoObject reco_object;
    
    tEdep->SetBranchAddress("Event", &evEdep);
    
    // tout.Branch("event_reco", "EventReco", &event_reco);
    tout.Branch("reco_object", "reco_object", &reco_object);
    
    std::vector<dg_wire> wire_infos;
    
    std::vector<dg_wire> fired_wires;

    LOG("I","Loading wires lookup table");
    ReadWireInfos(fWireInfo, wire_infos);

    for(auto i= 992u; i < 993; i++)
    // for(auto i= 9u; i < 10; i++)
    // for(auto i=0u; i < tEdep->GetEntries(); i++)
    {
        tEdep->GetEntry(i);
        
        auto muon_trj = evEdep->Trajectories[0];
        auto vertex = evEdep->Primaries[0];
        auto trj = evEdep->Trajectories;

        if(muon_trj.GetPDGCode()!=13) continue;

        LOG("I","Test muon");
        std::cout << "event number : " << i 
                  << ", muon_momentum : (" 
                  << muon_trj.GetInitialMomentum().X() << ", " 
                  << muon_trj.GetInitialMomentum().Y() << ", " 
                  << muon_trj.GetInitialMomentum().Z() << ", " 
                  << muon_trj.GetInitialMomentum().T() << ")\n";
        LOG("","\n");
        
        LOG("I","Check Event Vertex in Fiducial Volume (10 cm from frames)");
        bool event_in_fiducial_volume = IsInFiducialVolume(vertex);
        
        if(!event_in_fiducial_volume){
            LOG("W","skipping event, not in fiducial volume");
            continue;
        }
        
        // create the non-smeared track (no E loss nor MCS) from initial muon momentum and direction
        Helix true_helix(muon_trj);

        if(USE_NON_SMEARED_TRACK){
            // define fired wires from the muon helix
            LOG("I", "Digitization of non seamred track");
            CreateDigitsFromHelix(true_helix, wire_infos, fired_wires);
        }else{
            // create digits from edepsim selected hits
            LOG("I", "Digitization of edepsim TG4HitSegments");
            CreateDigitsFromEDep(evEdep->SegmentDetectors["DriftVolume"], wire_infos, fired_wires);
            // throw "";
        }

        LOG("i", "Sort wires by true hit time");
        SortWiresByTime(fired_wires);

        RecoUtils::event_digits = &fired_wires;

        first_fired_wire = &fired_wires.front();

        LOG("I", "Check event with enough hits");
        bool good_event = PassSelectionNofHits(fired_wires);

        if(!good_event) continue;

        Helix helix_first_guess = GetHelixFirstGuess(true_helix);

        MinuitFitInfos fit_infos, fit_XZ, fit_ZY;

        LOG("I", "Split wires in vertical and horizontal");
        std::vector<dg_wire> hor_wires, ver_wires;
        for (auto& wire : fired_wires){
            if(wire.hor==true){ // horizontal
                hor_wires.push_back(wire);
            }else{// (wire.hor==false) // vertical
                ver_wires.push_back(wire);
            }}

        // sort wires (they should be already sorted)
        LOG("I","Sorting wire in time (MC truth)");
        SortWiresByTime(hor_wires);
        SortWiresByTime(ver_wires);

        horizontal_fired_digits = &hor_wires;
        vertical_fired_digits = &ver_wires;
        
        // convert measured TDC into drift distance (i.e dift circles)
        LOG("I","Converting measured TDC into drift time ");
        TDC2DriftDistance(fired_wires);
        LOG("ii","Converting horizontal wires");
        TDC2DriftDistance(*horizontal_fired_digits);
        LOG("ii","Converting vertical wires");
        TDC2DriftDistance(*vertical_fired_digits);

        LOG("I", "Creating drift circles from drift time");
        std::vector<Circle> vertical_drift_circles = Wire2DriftCircle(*vertical_fired_digits);
        std::vector<Circle> horizontal_drift_circles = Wire2DriftCircle(*horizontal_fired_digits);

        LOG("I","Find local tangets to drift circles in XZ plane");
        auto segmentsXZ = GetTrackSegments(vertical_drift_circles);
        track_segments_XZ = &segmentsXZ;
        LOG("R", TString::Format("Tangent to first TDC XZ plane m = %f, q = %f", segmentsXZ[0].m(), segmentsXZ[0].q()).Data());

        auto segmentsZY = GetTrackSegments(horizontal_drift_circles);
        track_segments_ZY = &segmentsZY;
        LOG("R", TString::Format("Tangent to first TDC XZ plane m = %f, q = %f", segmentsZY[0].m(), segmentsZY[0].q()).Data());

        Helix reco_helix;

        if(FITTING_STRATEGY==0){
            /*
                Fit the coordinate of the fired wires to
                provide an initial guess of the NLL method
            */
            PrintEventInfos(true_helix, helix_first_guess, fired_wires);
            
            LOG("I", "Fitting drift circles with NLL method on XZ PLANE");
            LOG("i", "Track First Guess : fitting wire coordinates on XZ PLANE");
            Line2D LineFit = RecoUtils::WiresLinearFit(*vertical_fired_digits);
            LOG("R", TString::Format("first guess : Slope (m) = %f, Intercept (q) =  %f", LineFit.m(), LineFit.q()).Data());
            LOG("i", "Fitting drift circles starting from first guess");
            std::cout << blue << "\n";
            Line2D RecoXZTrack = GetLineFitXZPlane(LineFit, fit_XZ);
            std::cout << def << "\n";
            
            LOG("I", "Fitting drift circles with NLL method on ZY PLANE");
            LOG("i", "Track First Guess : fitting wire coordinates on ZY PLANE");
            Circle CircleFit = RecoUtils::WiresCircleFit(*horizontal_fired_digits);
            LOG("R", TString::Format("first guess : Center = (%f, %f), R =  %f", 
                                            CircleFit.center_x(), CircleFit.center_y(), CircleFit.R()).Data());
            LOG("i", "Fitting drift circles starting from first guess");
            std::cout << blue << "\n";
            Circle RecoZYTrack = GetRecoZYTrack(CircleFit, fit_ZY);
            std::cout << def << "\n";

            // get reconstructed helix from separate fit
            reco_object.fit_infos.push_back(fit_XZ);
            reco_object.fit_infos.push_back(fit_ZY);
            reco_object.track_segments_ZY = *track_segments_ZY;
            reco_object.track_segments_XZ = *track_segments_XZ;

            reco_helix = Reconstruct(RecoZYTrack, 
                                     RecoXZTrack,
                                     *horizontal_fired_digits,
                                     *vertical_fired_digits);

            reco_helix.PrintHelixPars();                                     
        }else if(FITTING_STRATEGY==1){
            // RecoUtils::InitHelixPars(fired_wires, helix_first_guess);
            // PrintEventInfos(true_helix, helix_first_guess, fired_wires);
            reco_helix = Reconstruct(fit_infos,         // infos on the fit
                                     true_helix,        // truth
                                     helix_first_guess, // first guess
                                     wire_infos,        // all wire infos
                                     fired_wires);      // wires fired by the test helix
            reco_object.fit_infos.push_back(fit_infos);
        }

        reco_object.fired_wires = *RecoUtils::event_digits;
        reco_object.true_helix = true_helix;
        reco_object.pt_true = true_helix.R()*0.3*0.6;
        reco_object.reco_helix = reco_helix;
        reco_object.pt_reco = reco_helix.R()*0.3*0.6;
        
        reco_object.edep_file_input = fEDepInput;
        reco_object.digit_file_input = fDigitInput;
        reco_object.event_index = i;
        if(USE_NON_SMEARED_TRACK) reco_object.use_track_no_smearing = true;

        tout.Fill();
    }
    fout.cd();

    tout.Write();

    fout.Close();
}