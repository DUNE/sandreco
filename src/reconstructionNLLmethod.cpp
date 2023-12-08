#include <iostream>

#include "SANDRecoUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TLegend.h"

// using namespace RecoUtils;

std::vector<dg_tube>* RecoUtils::event_digits = nullptr;

TG4Event* evEdep = nullptr;

TGeoManager* geo = nullptr;

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

    // read geometry
    geo = (TGeoManager*)fEDep.Get("EDepSimGeometry");
    
    // initialize wiremap
    if(geo)
    {
        // std::cout<<__FILE__<<" "<<__LINE__<<"\n";
        RecoUtils::InitWireInfos(geo);

    }else{
        std::cout<<"TGeoManager null \n";
        throw "";
    }

    // MC info tree
    TTree* t = (TTree*)fDigit.Get("tDigit");
    TTree* tEdep = (TTree*)fEDep.Get("EDepSimEvents");
    
    const int nev = t->GetEntries();

    std::vector<dg_tube>* mydigit = nullptr;

    // digits
    // t->SetBranchAddress("dg_tube", &RecoUtils::event_digits);
    t->SetBranchAddress("dg_tube", &mydigit);
    tEdep->SetBranchAddress("Event", &evEdep);

    // output file branches
    std::vector<double>     impact_par_from_TDC;
    std::vector<double>     impact_par_estimated;
    std::vector<double>     trj_points_x;
    std::vector<double>     trj_points_y;
    std::vector<double>     trj_points_z;
    TLorentzVector          initial_momentum;
    TLorentzVector          initial_position;
    std::vector<double>     fired_digits_x;
    std::vector<double>     fired_digits_y;
    std::vector<double>     fired_digits_z;
    std::vector<double>     fired_digits_hor;
    double                  R_true, dip_true, Phi0_true;
    double                  R_reco, dip_reco, Phi0_reco;
    int                     h_true, h_reco;
    TVector3                x0_true;
    TVector3                x0_reco;

    TFile fout(fOutput, "RECREATE");
    TTree tout("reco", "t");

    tout.Branch("ImpactParFromTDC",    "std::vector<double>",            &impact_par_from_TDC);
    tout.Branch("ImpactParEstimated",    "std::vector<double>",          &impact_par_estimated);
    tout.Branch("edep_trj_points_x",        "std::vector<double>",       &trj_points_x);
    tout.Branch("edep_trj_points_y",        "std::vector<double>",       &trj_points_y);
    tout.Branch("edep_trj_points_z",        "std::vector<double>",       &trj_points_z);
    tout.Branch("edep_initial_momentum",    "TLorentzVector",            &initial_momentum);
    tout.Branch("edep_initial_position",    "TLorentzVector",            &initial_position);
    tout.Branch("fired_digits_x",           "std::vector<double>",       &fired_digits_x);
    tout.Branch("fired_digits_y",           "std::vector<double>",       &fired_digits_y);
    tout.Branch("fired_digits_z",           "std::vector<double>",       &fired_digits_z);
    tout.Branch("fired_digits_hor",         "std::vector<double>",       &fired_digits_hor);
    tout.Branch("R_true",                   &R_true,                     "R_true/D");
    tout.Branch("R_reco",                   &R_reco,                     "R_reco/D");
    tout.Branch("dip_true",                 &dip_true,                   "dip_true/D");
    tout.Branch("dip_reco",                 &dip_reco,                   "dip_reco/D");
    tout.Branch("Phi0_true",                &Phi0_true,                  "Phi0_true/D");
    tout.Branch("Phi0_reco",                &Phi0_reco,                  "Phi0_reco/D");
    tout.Branch("h_true",                   &h_true,                     "h_true/I");
    tout.Branch("h_reco",                   &h_reco,                     "h_reco/I");
    tout.Branch("x0_true",                  "TVector3",                  &x0_true);
    tout.Branch("x0_reco",                  "TVector3",                  &x0_reco);
    
    // for (auto i = 0; i < nev; i++)
    for (auto i = 0; i < 1; i++)
    {
        // if(i!=6) continue;
        t->GetEntry(i);
        tEdep->GetEntry(i);
        RecoUtils::event_digits = new std::vector<dg_tube>();
        for (auto& digit : *mydigit)
        {
            RecoUtils::event_digits->push_back(digit);
        }
        
        std::cout<<"evento : "<<i<<" numero di digits : "<<RecoUtils::event_digits->size()<<"\n";
        
        auto muon_trj = evEdep->Trajectories[0];

        if(muon_trj.GetPDGCode()!=13) continue;

        initial_momentum = muon_trj.GetInitialMomentum();
        initial_position = muon_trj.Points[0].GetPosition();
        
        Helix true_helix(muon_trj); // true helix
        
        R_true      = true_helix.R();
        dip_true    = true_helix.dip();
        Phi0_true   = true_helix.Phi0();
        h_true      = 1;            
        x0_true     = {true_helix.x0().X(),
                       true_helix.x0().Y(),
                       true_helix.x0().Z()};

        for(auto& point : muon_trj.Points){
            // point is TG4TRajectoryPoint
            TLorentzVector point_position = point.GetPosition();
            trj_points_x.push_back(point_position.X());
            trj_points_y.push_back(point_position.Y());
            trj_points_z.push_back(point_position.Z());
        }

        std::cout<<"-----------------------------------------\n";
        
        double p[7] = {R_true, 
                       dip_true, 
                       Phi0_true, 
                       h_true, 
                       x0_true.X(),
                       x0_true.Y(), 
                       x0_true.Z()};

        std::cout << "true helix : \n";

        true_helix.PrintHelixPars();

        // check
        // auto digit2 = RecoUtils::event_digits->at(2);
        // Line digit2_line = RecoUtils::GetLineFromDigit(digit2);
        // true_helix.SetHelixRangeFromDigit(digit2);
        // std::cout << "digit 2, RecoUtils::GetMinImpactParameter : " << RecoUtils::GetMinImpactParameter(true_helix, digit2_line);

        // check

        const double* pars = p;  
        
        auto fitted_pars = RecoUtils::GetHelixParameters(pars);

        R_reco = fitted_pars[0];
        dip_reco = fitted_pars[1];
        Phi0_reco = fitted_pars[2];
        h_reco = 1;
        x0_reco = {fitted_pars[4],
                   fitted_pars[5],
                   fitted_pars[6],
                   };
        
        std::cout<<"\n";

        std::cout<<"TRUE values - FITTED VALUES \n";

        std::cout << "R_true - R_reco       : " << R_true - R_reco << "[mm]\n";
        std::cout << "dip_true - dip_reco   : " << dip_true - dip_reco << "\n";
        std::cout << "Phi0_true - Phi0_reco : " << Phi0_true - Phi0_reco << "\n";
        std::cout << "x0_true - x0_reco     : " << (x0_true - x0_reco).Mag() << "[mm]\n";
        
        std::cout<<"-----------------------------------------\n";

        // run over digits
        for(auto j=0u; j<RecoUtils::event_digits->size(); j++)
        {
            auto digit = RecoUtils::event_digits->at(j);

            /* find s_min and s_max that gives the portion of the helix 
            in the plane containing the wire */
            true_helix.SetHelixRangeFromDigit(digit);

            // std::cout<<"helix in range s : ["<<true_helix.LowLim()<<","<<true_helix.UpLim()<<"]\n";

            Line digit_line = RecoUtils::GetLineFromDigit(digit);

            auto TDC_impact_par = RecoUtils::GetExpectedRadiusFromDigit(digit);
            
            auto estimated_impact_par = RecoUtils::GetMinImpactParameter(true_helix, digit_line);
            
            std::cout<<"impact par from TDC : "<<TDC_impact_par<<", estimated_impact_par (my minimizer) : "<<estimated_impact_par<<" mm \n";
            // std::cout<<"\n";
            
            impact_par_from_TDC.push_back(TDC_impact_par);
            impact_par_estimated.push_back(estimated_impact_par);
            fired_digits_x.push_back(digit.x);
            fired_digits_y.push_back(digit.y);
            fired_digits_z.push_back(digit.z);
            fired_digits_hor.push_back(digit.hor);
        }
        
        delete RecoUtils::event_digits;
    }
    tout.Fill();
    // // write output
    fout.cd();
    tout.Write();
    fout.Close();

    return 0;
}