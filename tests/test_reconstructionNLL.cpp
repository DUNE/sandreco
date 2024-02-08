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

const bool INCLUDE_SIGNAL_PROPAGATION = true;

const bool INCLUDE_HIT_TIME = true;

const bool INCLUDE_TDC_SMEARING = false;

const bool _DEBUG_ = false;

// CONSTANS

const double SAND_CENTER_X = 0;

const double SAND_CENTER_Y = -2384.73;

const double SAND_CENTER_Z = 23910;

const double TDC_SMEARING = 10; // ns

const double NEUTRINO_INTERACTION_TIME = 1; // ns

const double LIGHT_VELOCITY = 299.792458; // mm/ns

// -------------------------------------------

std::vector<double> SmearVariable(double mean, double sigma, int nof_points){

    TRandom3 r(0);

    std::vector<double> smeared_points;

    for (int i = 0; i < nof_points; ++i) 
        smeared_points.push_back(r.Gaus(mean, sigma));
    
    return smeared_points;
}

void ReadWireInfos(const std::string& fInput, std::vector<dg_wire>& wire_infos){
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

// double angleWithCenter(const dg_wire &a) {
//     return std::atan2(a.y - SAND_CENTER_Y, a.z - SAND_CENTER_Z);
// }

// void SortWires(std::vector<dg_wire>& wires){
//     // first wire is the most upstream (the lowest z)
//     std::sort(wire_infos.begin(), wire_infos.end(), [](const dg_wire &a, const dg_wire &b) {
//             return a.z < b.z;});
    
//     dg_wire firstWire = wire_infos.front();

//     // following wires are sorted by clockwise criterion
//     std::sort(wire_infos.begin() + 1, wire_infos.end(), [&firstWire](const dg_wire &a, const dg_wire &b) {
//         return angleWithCenter(a) < angleWithCenter(b);
//     });
// }

void CreateDigitsFromHelix(Helix& h, 
                           std::vector<dg_wire>& wire_infos, 
                           std::vector<dg_wire>& fired_wires){
    /* 
        - run over wires
        - set helix range
        - callRecoUtils::GetMinImpactParameter(const Helix& helix, const Line& line)
        - if it is < 5*sqrt(2) mm (is the half diagonal of a cell size of 10 mm)
          and the helix point is inside SAND x range
            - fill dg_wire info with
            - add to event_wire_digits
    */
    // usefull in case t_hit included in TDC
    double hit_timer = NEUTRINO_INTERACTION_TIME;
    double last_s_min = 0.;
    //
   // scan all the wires TDC and select only the onces with the closest impact par
   for(auto& w : wire_infos){

    Line l = RecoUtils::GetLineFromDigit(w);

    h.SetHelixRangeFromDigit(w);

    /* define point on the helix and on the line corresponding to the impact parameter
       t_min is how much do I move from the wire center to get to the line point closest to the helix.
       It has a sign.
    */
    double s_min, t_min;

    bool HasMinimized = false;

    double impact_par = RecoUtils::GetMinImpactParameter(h, l, s_min, t_min, HasMinimized);

    if((impact_par <= 5*sqrt(2)) & (fabs(t_min) <= w.wire_length)){

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
            std::cout << "hit_timer : " << hit_timer << " ns \n";
            last_s_min = s_min;
        }
        fired_wires.push_back(w);
    }
   }
   // oreder the wire hit by z and give the a hit time
//    SortWires(fired_wires);
}

Helix GetHelixFirstGuess(const Helix& input_helix){
    Helix first_guess(SmearVariable(input_helix.R(), 100, 1)[0], 
                      SmearVariable(input_helix.dip(), 0.1, 1)[0], 
                      input_helix.Phi0(), 
                       1, 
                      {SmearVariable(input_helix.x0().X(), 10, 1)[0],
                       input_helix.x0().Y(),
                       input_helix.x0().Z()});
    return first_guess;
}

void FillTrueInfos(RecoObject& reco_obj, const Helix& true_helix){
    reco_obj.true_helix = true_helix;
    reco_obj.pt_true = true_helix.R()*0.3*0.6;
}

void FillRecoInfos(RecoObject& reco_obj, const Helix& reco_helix){
    reco_obj.reco_helix = reco_helix;
    reco_obj.pt_reco = reco_helix.R()*0.3*0.6;
}

void Plot(TCanvas *canvas,
          std::vector<dg_wire>& wire_infos, 
          std::vector<dg_wire>& fired_wires,
          Helix& test_helix,
          Helix& reco_helix){
    
    TGraph *wire_infoZY = new TGraph();
    TGraph *wire_infoXZ = new TGraph();

    TGraph *fired_wiresZY = new TGraph();
    TGraph *fired_wiresXZ = new TGraph();

    TGraph *test_helixZY = new TGraph();
    TGraph *test_helixXZ = new TGraph();

    TGraph *reco_helixZY = new TGraph();
    TGraph *reco_helixXZ = new TGraph();

    
    for (auto i = 0u; i < wire_infos.size(); i++)
    {
        if(wire_infos[i].hor==1){
            wire_infoZY->SetPoint(i, wire_infos[i].z, wire_infos[i].y);
        }else{
            wire_infoXZ->SetPoint(i, wire_infos[i].x, wire_infos[i].z);
        }
    }

    for (auto i = 0u; i < fired_wires.size(); i++)
    {
        if(fired_wires[i].hor==1){
            fired_wiresZY->SetPoint(i, fired_wires[i].z, fired_wires[i].y);
        }else{
            fired_wiresXZ->SetPoint(i, fired_wires[i].x, fired_wires[i].z);
        }
    }

    auto test_helix_points  = test_helix.GetHelixPoints(-3*1e4,3*1e4);

    for (auto i = 0u; i < test_helix_points.size(); i++)
    {
        test_helixZY->SetPoint(i, test_helix_points[i].Z(), test_helix_points[i].Y());
        test_helixXZ->SetPoint(i, test_helix_points[i].X(), test_helix_points[i].Z());
    }

    auto reco_helix_points  = reco_helix.GetHelixPoints(-3*1e4,3*1e4);

    for (auto i = 0u; i < test_helix_points.size(); i++)
    {
        reco_helixZY->SetPoint(i, reco_helix_points[i].Z(), reco_helix_points[i].Y());
        reco_helixXZ->SetPoint(i, reco_helix_points[i].X(), reco_helix_points[i].Z());
    }
    
    canvas->Divide(2,1);

    canvas->cd(1);
    
    wire_infoZY->SetMarkerStyle(20);
    wire_infoZY->SetMarkerColor(kBlue);
    wire_infoZY->SetMarkerSize(0.3);
    
    fired_wiresZY->SetMarkerStyle(20);
    fired_wiresZY->SetMarkerColor(kRed);
    fired_wiresZY->SetMarkerSize(0.3);

    test_helixZY->SetMarkerStyle(20);
    test_helixZY->SetMarkerColor(kGreen);
    test_helixZY->SetMarkerSize(0.3);

    reco_helixZY->SetMarkerStyle(20);
    reco_helixZY->SetMarkerColor(kOrange);
    reco_helixZY->SetMarkerSize(0.3);

    wire_infoZY->GetXaxis()->SetLimits(21500,26000);
    wire_infoZY->GetYaxis()->SetRangeUser(-4500, -300);

    wire_infoZY->Draw("AP");
    fired_wiresZY->Draw("P");
    test_helixZY->Draw("PL");
    reco_helixZY->Draw("PL");

    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(wire_infoZY, "all wires", "p");
    legend->AddEntry(fired_wiresZY, "fired wires", "p");
    legend->AddEntry(test_helixZY, "test helix", "p");
    legend->AddEntry(reco_helixZY, "reco helix", "p");
    legend->SetTextSize(0.04);
    legend->Draw();

    canvas->cd(2);

    wire_infoXZ->SetMarkerStyle(20);
    wire_infoXZ->SetMarkerColor(kBlue);
    wire_infoXZ->SetMarkerSize(0.3);

    fired_wiresXZ->SetMarkerStyle(20);
    fired_wiresXZ->SetMarkerColor(kRed);
    fired_wiresXZ->SetMarkerSize(0.3);

    test_helixXZ->SetMarkerStyle(20);
    test_helixXZ->SetMarkerColor(kGreen);
    test_helixXZ->SetMarkerSize(0.3);

    reco_helixXZ->SetMarkerStyle(20);
    reco_helixXZ->SetMarkerColor(kOrange);
    reco_helixXZ->SetMarkerSize(0.3);

    wire_infoXZ->GetXaxis()->SetLimits(-1700, 1700);
    wire_infoXZ->GetYaxis()->SetRangeUser(21500,26000);

    wire_infoXZ->Draw("AP");
    fired_wiresXZ->Draw("P");
    test_helixXZ->Draw("PL");
    reco_helixXZ->Draw("PL");

    canvas->Print("tests/test_helix_vs_reco_helix.pdf");

}

double TDC2ImpactPar(const dg_wire& wire){
    
    double assumed_signal_propagation = 0.;
    static double assumed_hit_time = 0.;
    static dg_wire previous_wire = *first_fired_wire;

    if(&wire == first_fired_wire){ // if a new helix is being tested, reset t_hit timer
        assumed_hit_time = NEUTRINO_INTERACTION_TIME;
        previous_wire = *first_fired_wire;}

    if(INCLUDE_SIGNAL_PROPAGATION) 
        assumed_signal_propagation = (wire.wire_length/2.) / sand_reco::stt::v_signal_inwire;

    if(INCLUDE_HIT_TIME)
        assumed_hit_time += fabs(wire.z - previous_wire.z) / LIGHT_VELOCITY;
    
    std::cout
            //  << "previous_wire.t_hit: "<<previous_wire.t_hit
             << " assumed_hit_time: "<<assumed_hit_time
             << " true t_hit: "<<wire.t_hit<<"\n";
    previous_wire = wire;
    return (wire.tdc - assumed_signal_propagation - assumed_hit_time) * sand_reco::stt::v_drift;
}

double NLL(Helix& h,
           std::vector<dg_wire>& digits){
    
    static int calls   = 0;
    double nll         = 0.;
    const double sigma = 0.2; // 200 mu_m = 0.2 mm
    std::cout << "call number : "<<calls<<"\n";
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

void PrintEventInfos(int i,
                     const Helix& test_helix,
                     const Helix& helix_first_guess,
                     const std::vector<dg_wire>& fired_wires){
        std::cout<<"\n";
        std::cout<< "event number : " << i <<", nof of fired wires : "<<fired_wires.size()<<"\n";
        std::cout<<"\n";
        std::cout<<"test helix : \n";
        test_helix.PrintHelixPars();
        std::cout<<"\n";
        std::cout<<"first guess : \n";
        helix_first_guess.PrintHelixPars();
        std::cout<<"\n";
        
}

int main(int argc, char* argv[]){

    // read edepsim input to create realistic muon tracks
    TFile fEDep("/storage/gpfs_data/neutrino/users/gi/sand-reco/tests/test_drift_edep.root", "READ");
    
    // file output to store reconstructed Helix
    TFile fout("tests/test_101_reconstruct_NLLmethod.root", "RECREATE");

    // text file with wire infos extrapolated from digitization
    std::string fWireInfo = "/storage/gpfs_data/neutrino/users/gi/sand-reco/tests/wireinfo.txt";
    
    TTree* tEdep  = (TTree*)fEDep.Get("EDepSimEvents");

    TTree tout("tReco", "tReco");

    EventReco event_reco;
    
    tEdep->SetBranchAddress("Event", &evEdep);

    tout.Branch("event_reco", "EventReco", &event_reco);

    // TCanvas *canvas = new TCanvas("c","c" , 1000, 800);

    // canvas->Print("tests/test_helix_vs_reco_helix.pdf[");
    
    std::vector<dg_wire> wire_infos;
    
    ReadWireInfos(fWireInfo, wire_infos);

    // for(auto i=0u; i < tEdep->GetEntries(); i++)
    for(auto i=0u; i < 2; i++)
    {
        // if(i!=78) continue;
        
        tEdep->GetEntry(i);

        auto muon_trj = evEdep->Trajectories[0];

        if(muon_trj.GetPDGCode()!=13) continue;

        // create the ideal helix (no E loss nor MCS) from initial muon momentum and direction
        Helix test_helix(muon_trj);
        
        Helix helix_first_guess = GetHelixFirstGuess(test_helix);

        std::vector<dg_wire> fired_wires;

        // define fired wires from the muon helix
        CreateDigitsFromHelix(test_helix, wire_infos, fired_wires);

        // // sort fired wires by z coordinate
        // std::sort(fired_wires.begin(), fired_wires.end(), [](const dg_wire &a, const dg_wire &b) {
        // return a.z < b.z;
        // });

        first_fired_wire = &fired_wires.front();

        if(fired_wires.size()<3){ // at least 4 hits required
            continue;
        }else{
            PrintEventInfos(i, test_helix, helix_first_guess, fired_wires);
        }

        int TMinuitStatus = -1;

        MinuitFitInfos fit_infos;

        auto reco_helix = Reconstruct(fit_infos,         // infos on the fit
                                      test_helix,        // truth
                                      helix_first_guess, // first gues
                                      wire_infos,        // all wire infos
                                      fired_wires);      // wires fired by the test helix

        RecoObject reco_object;

        reco_object.fired_wires = fired_wires;

        reco_object.helix_first_guess = helix_first_guess;

        reco_object.fit_infos = fit_infos;

        for(auto& w : fired_wires){
            reco_object.impact_par_from_TDC.push_back(w.tdc * sand_reco::stt::v_drift);
            auto l = RecoUtils::GetLineFromDigit(w);
            reco_object.impact_par_estimated.push_back(RecoUtils::GetMinImpactParameter(test_helix, l));
        }

        FillTrueInfos(reco_object, test_helix);

        FillRecoInfos(reco_object, reco_helix);

        event_reco.event_index = i;

        event_reco.event_fired_wires = fired_wires;

        event_reco.nof_digits = fired_wires.size();

        event_reco.reco_object = reco_object;

        tout.Fill();                          

        // TCanvas *canvas_i = new TCanvas(Form("event%d",i), Form("event%d",i), 800, 600);

        // canvas_i->SetTitle(Form("event%d",i));

        // Plot(canvas_i, wire_infos, fired_wires, test_helix, reco_helix);

        // delete canvas_i;

    }
    // canvas->Print("tests/test_helix_vs_reco_helix.pdf]");
    fout.cd();

    tout.Write();

    fout.Close();
}