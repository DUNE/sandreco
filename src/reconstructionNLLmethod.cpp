#include <iostream>

#include "SANDRecoUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TG4Event.h"

std::vector<dg_tube>* RecoUtils::event_digits = nullptr;

TG4Event* evEdep = nullptr;

TGeoManager* geo = nullptr;

void FillTrueInfos(RecoObject& reco_obj, const TG4Trajectory& trj){
        // fill the reco_obj true MC infos given by edepsim
        Helix                      true_helix(trj); // true helix
        auto initial_momentum    = trj.GetInitialMomentum();
        reco_obj.traj_edep_index = trj.GetTrackId();
        reco_obj.true_helix      = true_helix;
        reco_obj.pt_true         = sqrt(initial_momentum.Pz()*initial_momentum.Pz()+ 
                                        initial_momentum.Py()*initial_momentum.Py());
        
        for(auto& point : trj.Points) reco_obj.trj_points.push_back(point.GetPosition());
        std::cout << "true helix : \n";
        true_helix.PrintHelixPars();
}

void FillRecoInfos(RecoObject& reco_obj, const Helix& helix_initial_guess){
        auto fitted_pars        = RecoUtils::GetHelixParameters(helix_initial_guess);
        double R_reco                  = fitted_pars[0];
        double dip_reco                = fitted_pars[1];
        double Phi0_reco               = fitted_pars[2];
        int    h_reco                  = 1;
        TVector3 x0_reco               = {fitted_pars[4],
                                          fitted_pars[5],
                                          fitted_pars[6],};
        reco_obj.reco_helix     = Helix(R_reco, dip_reco, Phi0_reco, h_reco, x0_reco);
        reco_obj.pt_reco        = R_reco*0.3*0.6;
        
}

int main(int argc, char* argv[]){
    
    if(argc<4)
    {
        std::cout<<"ReconstructNLLmethod <EDEP FILE> <DIGIT FILE> <FILE OUTPUT NAME>\n";
        throw "";
    }

    const char* InputEDep = argv[1];

    const char* InputDigit = argv[2];

    const char* fOutput = argv[3];

    // read digitization file
    TFile fEDep(InputEDep, "READ");
    TFile fDigit(InputDigit, "READ");
    TFile fout(fOutput, "RECREATE");

    // read geometry
    geo = (TGeoManager*)fEDep.Get("EDepSimGeometry");
    
    // initialize wiremap
    if(geo)
    {
        RecoUtils::InitWireInfos(geo);

    }else{
        std::cout<<"TGeoManager null \n";
        throw "";
    }

    // MC info tree
    TTree* tDigit = (TTree*)fDigit.Get("tDigit");
    TTree* tEdep  = (TTree*)fEDep.Get("EDepSimEvents");
    TTree tout("reco", "t");
    
    const int nev = tDigit->GetEntries();

    // digits
    tDigit->SetBranchAddress("dg_tube", &RecoUtils::event_digits);
    tEdep->SetBranchAddress("Event", &evEdep);

    // output file branches
    EventReco               event_reco;
    std::vector<RecoObject> reco_infos;
    std::vector<double>     impact_par_from_TDC;
    std::vector<double>     impact_par_estimated;
    std::vector<double>     trj_points_x;
    std::vector<double>     trj_points_y;
    std::vector<double>     trj_points_z;
    std::vector<double>     fired_digits_x;
    std::vector<double>     fired_digits_y;
    std::vector<double>     fired_digits_z;
    std::vector<double>     fired_digits_hor;
    double                  R_true, dip_true, Phi0_true, pt_true;
    double                  R_reco, dip_reco, Phi0_reco, pt_reco;
    int                     h_true, h_reco;
    TVector3                x0_true;
    TVector3                x0_reco;

    tout.Branch("event_reco",               "EventReco",                 &event_reco);
    // tout.Branch("ImpactParFromTDC",         "std::vector<double>",       &impact_par_from_TDC);
    // tout.Branch("ImpactParEstimated",       "std::vector<double>",       &impact_par_estimated);
    // tout.Branch("edep_trj_points_x",        "std::vector<double>",       &trj_points_x);
    // tout.Branch("edep_trj_points_y",        "std::vector<double>",       &trj_points_y);
    // tout.Branch("edep_trj_points_z",        "std::vector<double>",       &trj_points_z);
    // tout.Branch("fired_digits_x",           "std::vector<double>",       &fired_digits_x);
    // tout.Branch("fired_digits_y",           "std::vector<double>",       &fired_digits_y);
    // tout.Branch("fired_digits_z",           "std::vector<double>",       &fired_digits_z);
    // tout.Branch("fired_digits_hor",         "std::vector<double>",       &fired_digits_hor);
    // tout.Branch("R_true",                   &R_true,                     "R_true/D");
    // tout.Branch("R_reco",                   &R_reco,                     "R_reco/D");
    // tout.Branch("dip_true",                 &dip_true,                   "dip_true/D");
    // tout.Branch("dip_reco",                 &dip_reco,                   "dip_reco/D");
    // tout.Branch("Phi0_true",                &Phi0_true,                  "Phi0_true/D");
    // tout.Branch("Phi0_reco",                &Phi0_reco,                  "Phi0_reco/D");
    // tout.Branch("pt_true",                  &pt_true,                    "pt_true/D");
    // tout.Branch("pt_reco",                  &pt_reco,                    "pt_reco/D");
    // tout.Branch("h_true",                   &h_true,                     "h_true/I");
    // tout.Branch("h_reco",                   &h_reco,                     "h_reco/I");
    // tout.Branch("x0_true",                  "TVector3",                  &x0_true);
    // tout.Branch("x0_reco",                  "TVector3",                  &x0_reco);
    
    // for (auto i = 0; i < nev; i++)
    for (auto i = 0; i < 3; i++)
    {
        tDigit->GetEntry(i);

        tEdep->GetEntry(i);

        auto muon_trj = evEdep->Trajectories[0];

        if(muon_trj.GetPDGCode()!=13) continue;

        RecoObject reco_object; // here should be a run over trajectories in future, now we have 1 event, 1 mu
        
        Helix true_helix(muon_trj); 

        FillTrueInfos(reco_object, muon_trj);

        FillRecoInfos(reco_object, true_helix); // feeded helix is the initial guess for the minimizer
        
        std::cout<<"\n";

        std::cout<<"TRUE values - FITTED VALUES \n";

        std::cout << "R_true - R_reco             : " << reco_object.true_helix.R() - reco_object.reco_helix.R() << "[mm]\n";
        std::cout << "dip_true - dip_reco         : " << reco_object.true_helix.dip() - reco_object.reco_helix.dip() << "\n";
        std::cout << "Phi0_true - Phi0_reco       : " << reco_object.true_helix.Phi0() - reco_object.reco_helix.Phi0() << "\n";
        std::cout << "x0_true - x0_reco           : " << (reco_object.true_helix.x0() - reco_object.reco_helix.x0()).Mag() << "[mm]\n";
        std::cout << "(pt_true - pt_reco)/pt_true : " << (reco_object.pt_true - reco_object.pt_reco)/reco_object.pt_true << "\n";
        
        std::cout<<"-----------------------------------------\n";

        // run over digits
        for(auto j=0u; j<RecoUtils::event_digits->size(); j++)
        {
            auto digit = RecoUtils::event_digits->at(j);

            /* find s_min and s_max that gives the portion of the helix 
            in the plane containing the wire */
            true_helix.SetHelixRangeFromDigit(digit);

            // std::cout<<"helix in range s : ["<<true_helix.LowLim()<<","<<true_helix.UpLim()<<"]\n";

            Line digit_line           = RecoUtils::GetLineFromDigit(digit);
            auto TDC_impact_par       = RecoUtils::GetExpectedRadiusFromDigit(digit);
            auto estimated_impact_par = RecoUtils::GetMinImpactParameter(true_helix, digit_line);
            
            std::cout<<"impact par from TDC : "<<TDC_impact_par<<", estimated_impact_par (my minimizer) : "<<estimated_impact_par<<" mm \n";

            impact_par_from_TDC.push_back(TDC_impact_par);
            impact_par_estimated.push_back(estimated_impact_par);
            fired_digits_x.push_back(digit.x);
            fired_digits_y.push_back(digit.y);
            fired_digits_z.push_back(digit.z);
            fired_digits_hor.push_back(digit.hor);

            reco_object.impact_par_from_TDC.push_back(TDC_impact_par);
            reco_object.impact_par_estimated.push_back(estimated_impact_par);
            reco_object.fired_wires.push_back(digit);
        }// run over digits
        
        event_reco.event_index = i;
        event_reco.event_fired_wires = *RecoUtils::event_digits;
        event_reco.reco_infos.push_back(reco_object);
    }// run over event

    tout.Fill();
    // // write output
    fout.cd();
    tout.Write();
    fout.Close();

    return 0;
}