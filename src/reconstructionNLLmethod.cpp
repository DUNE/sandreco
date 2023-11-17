#include <iostream>

#include "SANDRecoUtils.h"

#include "TFile.h"
#include "TTree.h"

std::vector<dg_tube>* event_digits = nullptr;

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
    // TEST________________________________________
    // double R = 2.3;
    // double Phi0 = TMath::Pi()/4.;
    // int hel = 1;
    // double dip = 5.;

    // double s_star = 0; // so that Phi=Phi0
    
    // TVector3 x0 = {0., 0., 0.}; //

    // // R, dip, Phi0, h, x0
    // Helix h(R, dip, Phi0, hel, x0); // this helix crosses (0,0,0) for s=0
    
    // auto helix_at_s_star = h.GetPointAt(s_star);
    // // mx, my, ax, ay, z0 -> line that crosses (0,0,0)
    // Line l(2., 0., 0., R*sin(TMath::Pi()/4.), R*cos(TMath::Pi()/4.)); // this line crosses (0,0,0) for t=0

    // auto line_at_t_star = l.GetPointAt(0.);

    // std::cout<<"s_star : "<<s_star<<"\n";
    // std::cout<<"helix at s=s_star hit the point at : "<<helix_at_s_star.X()<<" "<<helix_at_s_star.Y()<<" "<<helix_at_s_star.Z()<<"\n";
    // std::cout<<"t_star : 0 \n";
    // std::cout<<"line at t=0 hit the point at : "<<line_at_t_star.X()<<" "<<line_at_t_star.Y()<<" "<<line_at_t_star.Z()<<"\n";
    // std::cout<<"GetImpactParameter(h, l , s_star, t_star) : "<<RecoUtils::GetImpactParameter(h, l , s_star, 0)<<"\n";
    // std::cout<< RecoUtils::GetMinImpactParameter(h,l) <<"\n";
    // TEST________________________________________
    
    // read digitization file
    TFile fEDep(InputEDep, "READ");
    TFile fDigit(InputDigit, "READ");

    // read geometry
    geo = (TGeoManager*)fEDep.Get("EDepSimGeometry");
    
    // initialize wiremap
    if(geo)
    {
        std::cout<<__FILE__<<" "<<__LINE__<<"\n";
        RecoUtils::InitWireInfos(geo);

    }else{
        std::cout<<"TGeoManager null \n";
        throw "";
    }

    // MC info tree
    TTree* t = (TTree*)fDigit.Get("tDigit");
    TTree* tEdep = (TTree*)fEDep.Get("EDepSimEvents");
    
    const int nev = t->GetEntries();

    // digits
    t->SetBranchAddress("dg_tube", &event_digits);
    tEdep->SetBranchAddress("Event", &evEdep);

    // check radius of each digit
    std::vector<double> radii;
    std::vector<double> trj_points_x;
    std::vector<double> trj_points_y;
    std::vector<double> trj_points_z;
    std::vector<double> fired_digits_x;
    std::vector<double> fired_digits_y;
    std::vector<double> fired_digits_z;
    std::vector<double> helix_pars_from_edep;

    TFile fout(fOutput, "RECREATE");
    TTree tout("tTDCRadius", "t");
    tout.Branch("ExpectedRadiusFromTDC", "std::vector<double>", &radii);
    tout.Branch("edep_trj_points_x", "std::vector<double>", &trj_points_x);
    tout.Branch("edep_trj_points_y", "std::vector<double>", &trj_points_y);
    tout.Branch("edep_trj_points_z", "std::vector<double>", &trj_points_z);
    tout.Branch("fired_digits_x", "std::vector<double>", &fired_digits_x);
    tout.Branch("fired_digits_y", "std::vector<double>", &fired_digits_y);
    tout.Branch("fired_digits_z", "std::vector<double>", &fired_digits_z);
    tout.Branch("helix_pars_from_edep", "std::vector<double>", &helix_pars_from_edep);

    // for (auto i = 0; i < nev; i++)
    for (auto i = 0; i < 1; i++)
    {
        // if(i!=6) continue;
        t->GetEntry(i);
        tEdep->GetEntry(i);
        std::cout<<"evento : "<<i<<" numero di digits : "<<event_digits->size()<<"\n";
        
        auto muon_trj = evEdep->Trajectories[0];

        if(muon_trj.GetPDGCode()!=13) continue;

        for(auto& point : muon_trj.Points){
            // point is TG4TRajectoryPoint
            TLorentzVector point_position = point.GetPosition();
            trj_points_x.push_back(point_position.X());
            trj_points_y.push_back(point_position.Y());
            trj_points_z.push_back(point_position.Z());
        }

        Helix true_helix(muon_trj); // true helix

        helix_pars_from_edep.push_back(true_helix.R());
        helix_pars_from_edep.push_back(true_helix.Phi0());
        helix_pars_from_edep.push_back(true_helix.dip());
        helix_pars_from_edep.push_back(true_helix.h());
        helix_pars_from_edep.push_back(true_helix.x0().X());
        helix_pars_from_edep.push_back(true_helix.x0().Y());
        helix_pars_from_edep.push_back(true_helix.x0().Z());

        // run over digits
        for(auto j=0u; j<event_digits->size(); j++)
        {
            double s_min, t_min;
            bool has_minimized = 0;
            auto digit = event_digits->at(j);

            Line digit_line = RecoUtils::GetLineFromDigit(digit);
            true_helix.SetHelixRangeFromDigit(digit);
            digit_line.SetLineRangeFromDigit(digit);

            TVector2 helix_center = {true_helix.x0().Z() - true_helix.R()*cos(true_helix.Phi0()),
                                     true_helix.x0().Y() - true_helix.R()*sin(true_helix.Phi0())}; 

            TVector2 diff = {digit.z - helix_center.X(), digit.y - helix_center.Y()};

            std::cout<<"R - diff ; "<<true_helix.R() - diff.Mod()<<"\n";                                           

            auto reco_radius = RecoUtils::GetExpectedRadiusFromDigit(digit);
            auto true_helix_radius = RecoUtils::GetMinImpactParameter(true_helix,  digit_line, s_min, t_min, has_minimized);

            std::cout<<"angle x0-z "<< TMath::ATan2(true_helix.GetPointAt(s_min).Z() - helix_center.X(),true_helix.GetPointAt(s_min).Y() - helix_center.Y()) <<"\n";                                
            std::cout<<"angle digit-z "<< TMath::ATan2(digit.z - helix_center.X(), digit.y - helix_center.Y()) <<"\n";  
            std::cout<<"digit x y z hor did : "<<digit.x<<", "<<digit.y<<", "<<digit.z<<", "<<digit.hor<<", "<<digit.did<<"\n";
            std::cout<<"expected radius : "<<reco_radius<<", true helix radius : "<<true_helix_radius<<" mm \n";
            std::cout<<"\n";
            
            radii.push_back(reco_radius);
            fired_digits_x.push_back(digit.x);
            fired_digits_y.push_back(digit.y);
            fired_digits_z.push_back(digit.z);
        }
    }
    tout.Fill();
    // write output
    fout.cd();
    tout.Write();
    fout.Close();

    return 0;
}