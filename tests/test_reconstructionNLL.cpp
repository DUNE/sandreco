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

std::vector<dg_wire>* RecoUtils::event_digits = nullptr;

TGeoManager* geo = nullptr;

TG4Event* evEdep = nullptr;

dg_wire* first_fired_wire = nullptr;

// SIMULATION SETTINGS ----------------------

bool INCLUDE_SIGNAL_PROPAGATION = false;

bool INCLUDE_HIT_TIME = false;

bool INCLUDE_TDC_SMEARING = false;

bool _DEBUG_ = false;

bool USE_NON_SMEARED_TRACK = false;

// CONSTANS

const double SAND_CENTER_X = 0.;

const double SAND_CENTER_Y = -2384.73;

const double SAND_CENTER_Z = 23910.;

const double DRIFT_CELL_DIM = 10.; // 1 cm

const double TDC_SMEARING = 1.; // ns

const double NEUTRINO_INTERACTION_TIME = 1.; // ns

const double LIGHT_VELOCITY = 299.792458; // mm/ns

// -------------------------------------------

std::vector<double> SmearVariable(double mean, double sigma, int nof_points){

    TRandom3 r(0);

    std::vector<double> smeared_points;

    for (int i = 0; i < nof_points; ++i) 
        smeared_points.push_back(r.Gaus(mean, sigma));
    
    return smeared_points;
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
        wire.hor =  (wire_values[5]==0) ? 1 : 0;
        wire_infos.push_back(wire);
    }
    stream.close();
}

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
            w.tdc = SmearVariable(w.drift_time, TDC_SMEARING, 1)[0];
        
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
    // for (auto& hit : hits)
    // {
    //     if((evEdep->Trajectories[hit.PrimaryId].GetPDGCode()==PDG)
    //     &&(evEdep->Trajectories[hit.PrimaryId].GetParentId()==-1))
    //         filtered_hits.push_back(hit);
    // }
    // return filtered_hits;
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
    // look for already fired wires and update tdc
    for(auto& f : fired_wires){
        if((f.did == new_fired.did) && (f.tdc > new_fired.tdc)){
            found = true;
            f.tdc = new_fired.tdc;
            f.drift_time = new_fired.drift_time;
            f.signal_time = new_fired.signal_time;
            f.t_hit = new_fired.t_hit;
        }
    }
    if(!found) fired_wires.push_back(new_fired);
}

void Hits2File(const std::vector<TG4HitSegment>& hits){
    
    std::ofstream outputFile("tests/this_event_hits.txt");

    outputFile << "event_index,start_x,start_y,start_z,stop_x,stop_y,stop_z" << std::endl;

    for (const auto& hit : hits) {
        outputFile << evEdep->EventId <<","<< hit.GetStart().X() << "," << hit.GetStart().Y() << "," << hit.GetStart().Z() << ",";
        
        outputFile << hit.GetStop().X() << "," << hit.GetStop().Y() << "," << hit.GetStop().Z() << std::endl;
    }
    
    // Chiusura del file
    outputFile.close();

}

void CreateDigitsFromEDep(const std::vector<TG4HitSegment>& hits,
                          const std::vector<dg_wire>& wire_infos, 
                          std::vector<dg_wire>& fired_wires){
    /*
    Perform digitization of edepsim hits
    */
    fired_wires.clear();
    
    // filter only muon hit segments
    auto muon_hits = FilterHits(hits, 13);

    // 4 checks
    // Hits2File(muon_hits);
    // std::ofstream outputFile("tests/this_event_hits.txt"); // checks
    // outputFile << "event_index,start_x,start_y,start_z,stop_x,stop_y,stop_z,wpoint_x,wpoint_y,wpoint_z,hpoint_x,hpoint_y,hpoint_z" << std::endl; // checks
    //

    for(auto& hit : muon_hits){
        
        for(auto& wire : wire_infos){

            auto hit_middle = (hit.GetStart()+hit.GetStop())*0.5;
            
            if(fabs(wire.z - hit_middle.Z()) <= DRIFT_CELL_DIM * 0.5){//first scan along z to find the wire plane

                double closest2Hit, closest2Wire;
                double hit2wire_dist = Hit2WireDistance(wire, hit, closest2Hit, closest2Wire, outputFile);

                bool pass = (wire.hor==1) ? 
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
    outputFile.close();
}

Helix GetHelixFirstGuess(const Helix& input_helix){
    Helix first_guess(SmearVariable(input_helix.R(), 100, 1)[0], // smear of 10 cm around true R
                      SmearVariable(input_helix.dip(), 0.1, 1)[0], 
                      input_helix.Phi0(), 
                       1, 
                      {SmearVariable(input_helix.x0().X(), 10, 1)[0],
                       input_helix.x0().Y(),
                       input_helix.x0().Z()});
    return first_guess;
}

void FillTrueInfos(const Helix& true_helix, RecoObject& reco_obj){
    reco_obj.true_helix = true_helix;
    reco_obj.pt_true = true_helix.R()*0.3*0.6;
}

void FillRecoInfos(const Helix& test_helix, const Helix& reco_helix, RecoObject& reco_object){
    for(auto& w : *RecoUtils::event_digits){
            reco_object.impact_par_from_TDC.push_back(w.tdc * sand_reco::stt::v_drift);
            auto l = RecoUtils::GetLineFromDigit(w);
            reco_object.impact_par_estimated.push_back(RecoUtils::GetMinImpactParameter(test_helix, l));
    }
    reco_object.reco_helix = reco_helix;
    reco_object.pt_reco = reco_helix.R()*0.3*0.6;
}

void FillEventInfos(int ev_index, const RecoObject& reco_object, EventReco& event_reco){
        event_reco.event_index = ev_index;
        event_reco.event_fired_wires = *RecoUtils::event_digits;
        event_reco.nof_digits = RecoUtils::event_digits->size();
        event_reco.reco_object = reco_object;
}

void FillTreeOut(int ev_index, 
                 const Helix& test_helix,
                 const Helix& helix_first_guess,
                 const Helix& reco_helix,
                 RecoObject& reco_object,
                 EventReco& event_reco){
    
    reco_object.fired_wires = *RecoUtils::event_digits;
    reco_object.helix_first_guess = helix_first_guess;
    
    FillTrueInfos(test_helix, reco_object);
    FillRecoInfos(test_helix, reco_helix, reco_object);
    FillEventInfos(ev_index, reco_object, event_reco);
}

double TDC2ImpactPar(const dg_wire& wire){
    
    double assumed_signal_propagation = 0.;
    static double assumed_hit_time = 0.;
    static dg_wire previous_wire = *first_fired_wire;

    if(&wire == first_fired_wire){ // if a new helix is being tested, reset t_hit timer
        assumed_hit_time = NEUTRINO_INTERACTION_TIME;
        previous_wire = *first_fired_wire;
    }

    if(INCLUDE_SIGNAL_PROPAGATION) 
        assumed_signal_propagation = (wire.wire_length/2.) / sand_reco::stt::v_signal_inwire;

    if(INCLUDE_HIT_TIME){
        double step = sqrt((wire.x - previous_wire.x)*(wire.x - previous_wire.x) + 
                           (wire.y - previous_wire.y)*(wire.y - previous_wire.y) +
                           (wire.z - previous_wire.z)*(wire.z - previous_wire.z));
        assumed_hit_time += step / LIGHT_VELOCITY;
    }

    previous_wire = wire;
    return (wire.tdc - assumed_signal_propagation - assumed_hit_time) * sand_reco::stt::v_drift;
}

double NLL(Helix& h,
           std::vector<dg_wire>& digits){
    
    static int calls   = 0;
    double nll         = 0.;
    const double sigma = 0.2; // 200 mu_m = 0.2 mm
    for (auto& digit : digits)
    {
        Line l = RecoUtils::GetLineFromDigit(digit);
        
        h.SetHelixRangeFromDigit(digit);
        
        // r_estimated = impact parameter estimated as distance helix - wire
        double r_estimated = RecoUtils::GetMinImpactParameter(h,l);

        // r_observed = impact parameter from the observed tdc
        double r_observed = TDC2ImpactPar(digit);
        
        nll += (r_estimated - r_observed) * (r_estimated - r_observed) / (sigma * sigma);
    }
    calls ++;
    return sqrt(nll)/digits.size();
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
                   MinuitFitInfos& fit_infos,
                   int FitStrategy = 1){
    
    ROOT::Math::Functor functor(&FunctorNLL, 7);

    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");

    minimizer->SetFunction(functor);
    
    minimizer->SetLimitedVariable(0,    "R",    helix_initial_guess.R(), 100,  100,     1e5);
    minimizer->SetLimitedVariable(1,    "dip",  helix_initial_guess.dip(), 0.01,  -1.6,  1.6);
    minimizer->SetLimitedVariable(4,    "x0_x", helix_initial_guess.x0().X(), 1,  -1800, 1800);
    
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
    
    // fil output object fit_infos
    fit_infos.TMinuitFinalStatus = minimizer->Status();
    fit_infos.NIterations = minimizer->NIterations();
    fit_infos.MinValue = minimizer->MinValue();

    // print results of the minimization
    minimizer->PrintResults();
    
    if(fit_infos.TMinuitFinalStatus==4){ // failed fit
        if(minimizer->Strategy()==2){
            std::cout << "failed both with stratefy 1 and 2 \n";
        }else{
            // try to fit changing fit strategy
            GetRecoHelix(helix_initial_guess, fit_infos, 2);
        }
    }
    
    for (auto i = 0u; i < minimizer->NDim(); i++)
        fit_infos.parameters_errors.push_back(parsErrors[i]);

    // create reconstructed helix from fitted parameters
    Helix reco_helix(pars[0], pars[1], pars[2], pars[3], {pars[4], pars[5], pars[6]});

    return reco_helix;
}

Helix Reconstruct(MinuitFitInfos& fit_infos,
                  Helix& test_helix,
                  Helix& helix_first_guess,
                  std::vector<dg_wire>& wire_infos,
                  std::vector<dg_wire>& fired_wires){
    
    RecoUtils::event_digits = &fired_wires;

    auto reco_helix = GetRecoHelix(helix_first_guess, fit_infos);

    std::cout<<"\n";
    std::cout<< "true helix - reco helix \n";
    std::cout<<"R_true - R_reco     : "<< test_helix.R() - reco_helix.R() << "\n";
    std::cout<<"dip_true - dip_reco : "<< test_helix.dip() - reco_helix.dip() << "\n";
    std::cout<<"x0_true - x0_reco   : "<< test_helix.x0().X() - reco_helix.x0().X() << "\n";

    return reco_helix;

}

bool PassEventSelection(const std::vector<dg_wire>& fired_wires){
    /*
      Select events that have:
      - at least 3 hits on the banding plane (to make circular fit)
      - at least 2 hits on the XZ plane (to make linear fit)
    */
   int nof_ZY_hits = 0;
   int nof_XZ_hits = 0;
    if(fired_wires.size()<5){
        std::cout << "skipping event with nof fired wires: " << fired_wires.size() << "\n";
        std::cout << "\n";
        return false;
    }else{
        for(auto& f : fired_wires){
            if(f.hor == 1){
                nof_ZY_hits++;
            }else{
                nof_XZ_hits++;
            }
        }
        if((nof_ZY_hits<3)|(nof_XZ_hits<2)){
            std::cout << "skipping event, ZY_hits : " << nof_ZY_hits 
                      << ", XZ_hits : " << nof_XZ_hits << "\n";
            std::cout << "\n";
            return false;
        }else{
            return true;
        }
    };

}

void PrintEventInfos(int i,
                     const Helix& test_helix,
                     const Helix& helix_first_guess,
                     const std::vector<dg_wire>& fired_wires){
        std::cout<<"\n";
        std::cout<< "event number : " << i 
                 <<", nof of fired wires : "<< fired_wires.size() 
                 <<"\n";
        std::cout<<"\n";
        std::cout<<"test helix : \n";
        test_helix.PrintHelixPars();
        std::cout<<"\n";
        std::cout<<"first guess : \n";
        helix_first_guess.PrintHelixPars();
        std::cout<<"\n";
        
}

void help_input(){
    std::cout << "\n";
    std::cout << "./buil/bin/test_reconstructionNLL "
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
    std::cout << "--debug              : use higher verbosity for TMinuit \n";
}

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
    
    std::cout << "\n";
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
    
    EventReco event_reco;
    
    tEdep->SetBranchAddress("Event", &evEdep);
    
    tout.Branch("event_reco", "EventReco", &event_reco);
    
    std::vector<dg_wire> wire_infos;
    
    std::vector<dg_wire> fired_wires;

    ReadWireInfos(fWireInfo, wire_infos);

    // for(auto i= 0u; i < 10; i++)
    for(auto i=0u; i < tEdep->GetEntries(); i++)
    {
        // if(i!=78) continue;
        tEdep->GetEntry(i);
        
        // tDigit->GetEntry(i);

        auto muon_trj = evEdep->Trajectories[0];
        auto vertex = evEdep->Primaries[0];
        auto trj = evEdep->Trajectories;

        if(muon_trj.GetPDGCode()!=13) continue;

        // create the non-smeared track (no E loss nor MCS) from initial muon momentum and direction
        std::cout << "event : " << i 
                  << ", muon_momentum : (" 
                  << muon_trj.GetInitialMomentum().X() << ", " 
                  << muon_trj.GetInitialMomentum().Y() << ", " 
                  << muon_trj.GetInitialMomentum().Z() << ", " 
                  << muon_trj.GetInitialMomentum().T() << ")\n";
        
        Helix test_helix(muon_trj);

        Helix helix_first_guess = GetHelixFirstGuess(test_helix);

        if(USE_NON_SMEARED_TRACK){
            // define fired wires from the muon helix
            CreateDigitsFromHelix(test_helix, wire_infos, fired_wires);
        }else{
            // create digits from edepsim selected hits
            CreateDigitsFromEDep(evEdep->SegmentDetectors["DriftVolume"], wire_infos, fired_wires);
        }

        // sort fired wires by hit time (MC truth): usefull to subract assumed t_hit from measured TDC
        std::sort(fired_wires.begin(), fired_wires.end(), [](const dg_wire &a, const dg_wire &b) {
        return a.t_hit < b.t_hit;
        });

        bool good_event = PassEventSelection(fired_wires);

        if(!good_event){ // at least 4 hits required
            continue;
        }else{
            PrintEventInfos(i, test_helix, helix_first_guess, fired_wires);
        }
        
        MinuitFitInfos fit_infos;

        // usefull to subract t_hit from measured TDC
        first_fired_wire = &fired_wires.front();

        auto reco_helix = Reconstruct(fit_infos,         // infos on the fit
                                      test_helix,        // truth
                                      helix_first_guess, // first guess
                                      wire_infos,        // all wire infos
                                      fired_wires);      // wires fired by the test helix

        RecoObject reco_object;

        reco_object.fit_infos = fit_infos;
        event_reco.edep_file_input = fEDepInput;
        event_reco.digit_file_input = fDigitInput;
        if(USE_NON_SMEARED_TRACK) event_reco.use_track_no_smearing = true;

        FillTreeOut(i, test_helix,
                       helix_first_guess,
                       reco_helix,
                       reco_object, 
                       event_reco);

        tout.Fill();

    }
    fout.cd();

    tout.Write();

    fout.Close();
}