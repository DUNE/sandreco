#include <iostream>

#include "SANDRecoUtils.h"

#include "TFile.h"
#include "TTree.h"

std::vector<dg_tube>* event_digits = nullptr;

TGeoManager* geo = nullptr;

// SANDGeoManager* geo_manager = nullptr;
SANDGeoManager* geo_manager;

int main(int argc, char* argv[]){
    
    if(argc<4)
    {
        std::cout<<"ReconstructNLLmethod <EDEP FILE> <DIGIT FILE> <FILE OUTPUT NAME>\n";
        throw "";
    }

    const char* fInput = argv[1];

    const char* fGeometry = argv[2];

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
    TFile f(fInput, "READ");

    // read geometry
    geo = (TGeoManager*)f.Get("EDepSimGeometry");
    
    // initialize wiremap
    if(geo)
    {
        std::cout<<__FILE__<<" "<<__LINE__<<"\n";
        geo_manager->init(geo);
    }else{
        std::cout<<"TGeoManager null \n";
        throw "";
    }

    // MC info tree
    TTree* t = (TTree*)f.Get("tDigit");
    
    const int nev = t->GetEntries();

    // digits
    t->SetBranchAddress("dg_tube", &event_digits);

    for (auto i = 0; i < nev; i++)
    {
        if(i!=6) continue;
        t->GetEntry(i);
        std::cout<<"evento : "<<i<<" numero di digits : "<<event_digits->size()<<"\n";
        // const double* helix_pars_guess = RecoUtils::InitHelixPars(*event_digits);
        // auto helix_pars = RecoUtils::GetHelixParameters(helix_pars_guess);

        // run over digits
        for(auto j=0u; j<event_digits->size(); j++)
        {
            std::cout<<"fired wire : "<<event_digits->at(j).did<<"\n";
            auto expected_radius = RecoUtils::GetExpectedRadiusFromDigit(event_digits->at(j));
            // std::cout<<"digit "<<j<<" expected radius : "<<expected_radius<<"\n";
        }
    }
    

    return 0;
}