#include <iostream>

#include "SANDRecoUtils.h"

#include "TFile.h"
#include "TTree.h"

std::vector<dg_tube>* event_digits = nullptr;

int main(){

    double R = 2.3;
    double Phi0 = TMath::Pi()/4.;
    int hel = 1;
    double dip = 5.;

    double s_star = 0; // so that Phi=Phi0
    
    TVector3 x0 = {0., 0., 0.}; //

    // R, dip, Phi0, h, x0
    Helix h(R, dip, Phi0, hel, x0); // this helix crosses (0,0,0) for s=0
    
    auto helix_at_s_star = h.GetPointAt(s_star);
    // mx, my, ax, ay, z0 -> line that crosses (0,0,0)
    Line l(2., 0., 0., R*sin(TMath::Pi()/4.), R*cos(TMath::Pi()/4.)); // this line crosses (0,0,0) for t=0

    auto line_at_t_star = l.GetPointAt(0.);

    std::cout<<"s_star : "<<s_star<<"\n";
    std::cout<<"helix at s=s_star hit the point at : "<<helix_at_s_star.X()<<" "<<helix_at_s_star.Y()<<" "<<helix_at_s_star.Z()<<"\n";
    std::cout<<"t_star : 0 \n";
    std::cout<<"line at t=0 hit the point at : "<<line_at_t_star.X()<<" "<<line_at_t_star.Y()<<" "<<line_at_t_star.Z()<<"\n";
    std::cout<<"GetImpactParameter(h, l , s_star, t_star) : "<<RecoUtils::GetImpactParameter(h, l , s_star, 0)<<"\n";
    std::cout<< RecoUtils::GetMinImpactParameter(h,l) <<"\n";

    // read digitization file
    TFile f("/storage/gpfs_data/neutrino/users/gi/sand-reco/test_drift_digitize.root", "READ");

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
        const double* helix_pars_guess = RecoUtils::InitHelixPars(*event_digits);
        auto helix_pars = RecoUtils::GetHelixParameters(helix_pars_guess);
    }
    

    return 0;
}