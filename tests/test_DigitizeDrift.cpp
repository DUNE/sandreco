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

// DEFAULT SIMULATION SETTINGS_________________________________________________

bool INCLUDE_SIGNAL_PROPAGATION = false;

bool INCLUDE_HIT_TIME = false;

bool INCLUDE_TDC_SMEARING = false;

bool _DEBUG_ = false;

bool USE_NON_SMEARED_TRACK = false;

// CONSTANS____________________________________________________________________

const double SAND_CENTER_X = 0.;

const double SAND_CENTER_Y = -2384.73;

const double SAND_CENTER_Z = 23910.;

// !! HARDCODED (I know you are judging me, stop it)
// does not include frames
const double SAND_TRACKER_X_LENGTH = 3220.0;

// distance in mm from SUPERMOD frames
const double FIDUCIAL_CUT = 50.;

// height in mm
const double SUPERMOD_Y_HEIGHT[5] = {3755.16996258, // A1, A2
                                     3550.97659024, // B1, B2
                                     3129.97277395, // C1, C2
                                     1262.30218417,  // X1    
                                     2466.05424949 // X0   
};

const double MYLAR_2_MYLAR_DIST = 10.02; // mm

const double SENSE_2_SENSE_DIST = 20.28; // mm

const double TDC_SMEARING = 1.; // ns, default disabled

const double NEUTRINO_INTERACTION_TIME = 1.; // ns

const double LIGHT_VELOCITY = 299.792458; // mm/ns

// GLOBAL VARIABLES____________________________________________________________

std::vector<dg_wire>* RecoUtils::event_digits = nullptr;

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

void help_input(){
    std::cout << "\n";
    std::cout << red << "./build/bin/test_DigitizeDrift"
              << "-edep <EDep file> "
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

// std::map<int, std::vector<const TG4HitSegment*>> GroupHitsByTrajectory(const std::vector<TG4HitSegment>& hits)
// {
//     // Define a map to group pointers to hit segments by primary trajectory ID
//     std::map<int, std::vector<const TG4HitSegment*>> grouped_hits;
    
//     // Iterate through each hit segment in the input vector
//     for (const auto& hit : hits)
//     {
//         // Get the primary trajectory ID for the current segment
//         int hit_trj_id = hit.GetPrimaryId();
        
//         // Add a pointer to the current segment to the vector associated with the primary trajectory ID in the map
//         grouped_hits[hit_trj_id].push_back(&hit);
//     }
    
//     // Return the map containing pointers to hit segments grouped by primary trajectory ID
//     return grouped_hits;
// }

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
    // fired_wires.clear();

    // for(auto& hit : hits)
    for(auto i = 0u; i < hits.size(); i++){

        auto hit_middle = (hits[i].GetStop() + hits[i].GetStart())*0.5;
        auto hit_delta =  (hits[i].GetStop() - hits[i].GetStart());
        
        Line hit_line = Line(hits[i]);
        
        for(auto& wire : wire_infos){

            // first scan along z to find the wire plane
            bool is_hit_in_wire_plane = fabs(wire.z - hit_middle.Z()) < MYLAR_2_MYLAR_DIST * 0.5;
            
            if(is_hit_in_wire_plane){ // pass if hit z coodinate is found in the wire plane

                Line wire_line = RecoUtils::GetLineFromWire(wire);

                double closest_2_hit, closest_2_wire;

                double hit_2_wire_dist = RecoUtils::GetSegmentSegmentDistance(wire_line, hit_line, closest_2_hit, closest_2_wire);

                TVector3 w_point = wire_line.GetPointAt(closest_2_hit);
                
                TVector3 h_point = hit_line.GetPointAt(closest_2_wire);

                bool is_hit_in_cell = (wire.hor == true) ? 
                    fabs(w_point.Y() - h_point.Y()) < SENSE_2_SENSE_DIST * 0.5 : fabs(w_point.X() - h_point.X()) <= SENSE_2_SENSE_DIST * 0.5;

                if(is_hit_in_cell){ // pass if hit y (or x) coodinate is found in the hor (or vertical) wire cell

                    auto fired = RecoUtils::Copy(wire);
                    
                    fired.drift_time = hit_2_wire_dist / sand_reco::stt::v_drift;

                    fired.signal_time = (w_point - wire_line.GetLineUpperLimit()).Mag() / sand_reco::stt::v_signal_inwire;
                    
                    fired.t_hit = (hit_middle + hit_delta * closest_2_wire).T();

                    fired.tdc = fired.drift_time + fired.signal_time + fired.t_hit;

                    fired.hindex.push_back(i);

                    UpdateFiredWires(fired_wires, fired);
                }
            }
        }
    }
}

int main(int argc, char* argv[]){
    // std::cout << "DIGITIZE DRIFT \n";

    if (argc < 3 || argc > 11) {
        help_input();
    return -1;
    }

    const char* fEDepInput;

    const char* fDigitOutput;

    const char* fWireInfo;

    int index = 1;

    std::string trackerType = "DriftVolume";

    LOG("","\n");
    LOG("I", "Parsing inputs");

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
                fDigitOutput = argv[++index];
                std::cout << "Output file     " << fDigitOutput << std::endl;
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

    LOG("I", "Reading EDepSim Input");
    TFile fEDep(fEDepInput, "READ");

    TFile fout(fDigitOutput, "RECREATE");

    TTree tout("tDigit", "tDigit");

    TTree* tEdep = (TTree*)fEDep.Get("EDepSimEvents");

    geo = (TGeoManager*)fEDep.Get("EDepSimGeometry");

    // Checks for STT or DRIFT Chaber geometry
    if(geo->FindVolumeFast("STTtracker_PV")){
        std::cout<<"\n--- STT based simulation ---\n";
        trackerType="Straw";
    }
    else if (geo->FindVolumeFast("SANDtracker_PV")){
        std::cout<<"\n--- Drift based simulation ---\n";
        trackerType="DriftVolume";
    }
    else{
        std::cout<<"Error in retriving volume information from Geo Manager, exiting...\n";
        exit(-1);
    }

    std::vector<dg_wire> wire_infos;

    std::vector<dg_wire> fired_wires;

    tout.Branch("edep_file", "edep_file", &fEDepInput);
    
    tout.Branch("fired_wires", "fired_wires", &fired_wires);
    
    LOG("I","Loading wires lookup table");
    ReadWireInfos(fWireInfo, wire_infos);

    LOG("I", "Reading branch Event from EDepFile");
    tEdep->SetBranchAddress("Event", &evEdep);

    for(auto i=0u; i < tEdep->GetEntries(); i++)
    //for(auto i=0u; i < 10; i++)
    {
        LOG("ii", TString::Format("Processing Event %d", i).Data());
        tEdep->GetEntry(i);

        fired_wires.clear();

        CreateDigitsFromEDep(evEdep->SegmentDetectors[trackerType], wire_infos, fired_wires);
        
        LOG("i", "Sort wires by true hit time");
        SortWiresByTime(fired_wires);
        
        tout.Fill();
    }
    fout.cd();

    tout.Write();

    fout.Close();

    return 0;
}