#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TDirectory.h>

#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>


int makePlot(const std::string& input_root = "kf_trees.root",
                   const std::string& output_root = "kf_plots.root");


struct PointInfos{
  int step{};
  double z{};
  double x_pred{}, y_pred{}; // predicted
  double x_filt{}, y_filt{}; // filtered 
  double x_smooth{}, y_smooth{}; // smoothed 
  double x_meas{}, y_meas{}; // measured (orientation??)
};

 // ---------------- plot TRACCE (reco/smooth/filt/meas) ------------------------
static void makeTrackGraph(TTree* steps, TDirectory* outDir) {

  Int_t run=0, event=0, track_id=0, step_idx=0;
  steps->SetBranchAddress("run",&run);
  steps->SetBranchAddress("event",&event);
  steps->SetBranchAddress("track_id",&track_id);
  steps->SetBranchAddress("step_idx",&step_idx);

  Double_t z=0;
  Double_t x_meas=0, y_meas=0;
  Double_t x_pred=0, y_pred=0;
  Double_t x_filt=0, y_filt=0;
  Double_t x_smooth=0, y_smooth=0;

  steps->SetBranchAddress("z", &z);

  steps->SetBranchAddress("x_meas", &x_meas);
  steps->SetBranchAddress("y_meas", &y_meas);

  steps->SetBranchAddress("x_pred", &x_pred);
  steps->SetBranchAddress("y_pred", &y_pred);
  steps->SetBranchAddress("x_filt", &x_filt);
  steps->SetBranchAddress("y_filt", &y_filt);
  steps->SetBranchAddress("x_smooth",    &x_smooth);
  steps->SetBranchAddress("y_smooth",    &y_smooth);

  std::map<int, std::vector<PointInfos>> byTrack;
  const auto n = steps->GetEntries();
  for (auto i=0;i<n;++i){
    steps->GetEntry(i);
    PointInfos p{};
    p.step = step_idx; p.z = z;
    p.x_pred = x_pred; 
    p.y_pred = y_pred;
    p.x_filt = x_filt;
    p.y_filt = y_filt;
    p.x_smooth = x_smooth;
    p.y_smooth = y_smooth;
    p.x_meas = x_meas;
    p.y_meas = y_meas;
    byTrack[track_id].push_back(p);
  }

  outDir->cd();
  
  for (auto& kv : byTrack) {
    const int tid = kv.first;
    auto& v = kv.second;
    std::sort(v.begin(), v.end(), [](const PointInfos& a, const PointInfos& b){ return a.step < b.step; });
    const int N = (int)v.size();

    auto* yz_predicted = new TGraph(N);
    auto* yz_filtered  = new TGraph(N);
    auto* yz_smoothed  = new TGraph(N);
    auto* yz_measured  = new TGraph(N);

    auto* xz_predicted = new TGraph(N);
    auto* xz_filtered  = new TGraph(N);
    auto* xz_smoothed  = new TGraph(N);
    auto* xz_measured  = new TGraph(N);

    for (int i=0;i<N;++i){
      const auto& p = v[i];
      yz_predicted->SetPoint(i, p.z, p.y_pred);
      yz_filtered ->SetPoint(i, p.z, p.y_filt);
      yz_smoothed ->SetPoint(i, p.z, p.y_smooth);
      if (std::isfinite(p.y_meas)) yz_measured->SetPoint(yz_measured->GetN(), p.z, p.y_meas);

      xz_predicted->SetPoint(i, p.z, p.x_pred);
      xz_filtered ->SetPoint(i, p.z, p.x_filt);
      xz_smoothed ->SetPoint(i, p.z, p.x_smooth);
      if (std::isfinite(p.x_meas)) xz_measured->SetPoint(xz_measured->GetN(), p.z, p.x_meas);
    }

const std::string tidS = std::to_string(tid);

  std::string name = "track_" + tidS + "_yz_predicted";
  yz_predicted->SetName(name.c_str());
  name = "track_" + tidS + "_yz_filtered";
  yz_filtered->SetName(name.c_str());
  name = "track_" + tidS + "_yz_smoothed";
  yz_smoothed->SetName(name.c_str());
  name = "track_" + tidS + "_yz_measured";
  yz_measured->SetName(name.c_str());

  name = "track_" + tidS + "_xz_predicted";
  xz_predicted->SetName(name.c_str());
  name = "track_" + tidS + "_xz_filtered";
  xz_filtered->SetName(name.c_str());
  name = "track_" + tidS + "_xz_smoothed";
  xz_smoothed->SetName(name.c_str());
  name = "track_" + tidS + "_xz_measured";
  xz_measured->SetName(name.c_str());

  std::string title = "Track " + tidS + " YZ predicted; z [mm]; y [mm]";
  yz_predicted->SetTitle(title.c_str());
  title = "Track " + tidS + " YZ filtered;  z [mm]; y [mm]";
  yz_filtered->SetTitle(title.c_str());
  title = "Track " + tidS + " YZ smoothed;  z [mm]; y [mm]";
  yz_smoothed->SetTitle(title.c_str());
  title = "Track " + tidS + " YZ measured;  z [mm]; y [mm]";
  yz_measured->SetTitle(title.c_str());

  title = "Track " + tidS + " XZ predicted; z [mm]; x [mm]";
  xz_predicted->SetTitle(title.c_str());
  title = "Track " + tidS + " XZ filtered;  z [mm]; x [mm]";
  xz_filtered->SetTitle(title.c_str());
  title = "Track " + tidS + " XZ smoothed;  z [mm]; x [mm]";
  xz_smoothed->SetTitle(title.c_str());
  title = "Track " + tidS + " XZ measured;  z [mm]; x [mm]";
  xz_measured->SetTitle(title.c_str());

  std::string dname = "track_" + tidS;
  TDirectory* dtr = outDir->mkdir(dname.c_str());
  dtr->cd();

  yz_predicted->Write();
  yz_filtered->Write();
  yz_smoothed->Write();
  yz_measured->Write();
  xz_predicted->Write();
  xz_filtered->Write();
  xz_smoothed->Write();
  xz_measured->Write();

  outDir->cd();

    delete yz_predicted;
    delete yz_filtered;
    delete yz_smoothed;
    delete yz_measured;
    delete xz_predicted;
    delete xz_filtered;
    delete xz_smoothed;
    delete xz_measured;
  }
}

int makePlot(const std::string& input_root,
                   const std::string& output_root)
{
  TFile fin(input_root.c_str(), "READ");
  if (fin.IsZombie()) { std::cerr<<"Cannot open "<<input_root<<"\n"; return 1; }

  TTree* steps  = (TTree*)fin.Get("steps");
  TTree* tracks = (TTree*)fin.Get("tracks");
  if (!steps || !tracks) { std::cerr<<"Missing 'steps' or 'tracks' trees\n"; return 2; }

    Int_t run=0, event=0, track_id=0, step_idx=0;
    Double_t z=0.0;
    steps->SetBranchAddress("run",&run);
    steps->SetBranchAddress("event",&event);
    steps->SetBranchAddress("track_id",&track_id);
    steps->SetBranchAddress("step_idx",&step_idx);
      steps->SetBranchAddress("z", &z);

    //---true---
    Double_t x_true=0.0;
    Double_t y_true=0.0;
    Double_t invR_true=0.0;
    Double_t tanL_true=0.0;
    Double_t phi_true=0.0;
    steps->SetBranchAddress("x_true",   &x_true);
    steps->SetBranchAddress("y_true",   &y_true);
    steps->SetBranchAddress("invR_true",   &invR_true);
    steps->SetBranchAddress("tanL_true",   &tanL_true);
    steps->SetBranchAddress("phi_true",    &phi_true);

    //---smooth---
    Double_t x_smooth=0.0;
    Double_t y_smooth=0.0;
    Double_t invR_smooth=0.0;
    Double_t tanL_smooth=0.0;
    Double_t phi_smooth=0.0;
    steps->SetBranchAddress("x_smooth", &x_smooth);
    steps->SetBranchAddress("y_smooth", &y_smooth);
    steps->SetBranchAddress("invR_smooth", &invR_smooth);
    steps->SetBranchAddress("tanL_smooth", &tanL_smooth);
    steps->SetBranchAddress("phi_smooth",  &phi_smooth);

    //sigma smoothed
    Double_t sigma_x_smooth=0.0;
    Double_t sigma_y_smooth=0.0;
    Double_t sigma_invR_smooth=0.0;
    Double_t sigma_tanL_smooth=0.0;
    Double_t sigma_phi_smooth=0.0;
    if (steps->GetBranch("sigma_x_smooth"))    steps->SetBranchAddress("sigma_x_smooth",&sigma_x_smooth);
    if (steps->GetBranch("sigma_y_smooth"))    steps->SetBranchAddress("sigma_y_smooth",&sigma_y_smooth);
    if (steps->GetBranch("sigma_invR_smooth")) steps->SetBranchAddress("sigma_invR_smooth",&sigma_invR_smooth);
    if (steps->GetBranch("sigma_tanL_smooth")) steps->SetBranchAddress("sigma_tanL_smooth",&sigma_tanL_smooth);
    if (steps->GetBranch("sigma_phi_smooth"))  steps->SetBranchAddress("sigma_phi_smooth",&sigma_phi_smooth);

    Double_t x_meas=0.0, y_meas=0.0;
    Double_t x_pred=0.0, y_pred=0.0; 
    Double_t x_filt=0.0, y_filt=0.0;
    steps->SetBranchAddress("x_meas", &x_meas);
    steps->SetBranchAddress("y_meas", &y_meas);
    steps->SetBranchAddress("x_pred", &x_pred);
    steps->SetBranchAddress("y_pred", &y_pred);
    steps->SetBranchAddress("x_filt", &x_filt);
    steps->SetBranchAddress("y_filt", &y_filt);

    Double_t p_true=0.0; //MeV
    Double_t p_smooth=0.0; //MeV
    Double_t innov_pos=0.0, innov_ang=0.0, chi2=0.0;
    steps->SetBranchAddress("p_true", &p_true);
    steps->SetBranchAddress("p_smooth", &p_smooth);
    if (steps->GetBranch("innov_pos")) steps->SetBranchAddress("innov_pos",&innov_pos);
    if (steps->GetBranch("innov_ang")) steps->SetBranchAddress("innov_ang",&innov_ang);
    if (steps->GetBranch("chi2"))      steps->SetBranchAddress("chi2",&chi2);

  //---smoooth - true ---
  TH1D x_res_step("(x_smooth - x_true)","(x_{smooth}-x_{true}) [mm];mm;Entries",200,-2,2);
  TH1D y_res_step("(y_smooth - y_true)","(y_{smooth}-y_{true}) [mm];mm;Entries",200,-2,2);
  TH1D invR_res_step("(invR_smooth - invR_true)","invR_{smooth}-invR_{true};;Entries",200,-1e-2,1e-2);
  TH1D tanL_res_step ("(tanL_smooth - tanL_true)","tan#lambda_{smooth}-tan#lambda_{true};;Entries",200,-0.1,0.1);
  TH1D phi_res_step  ("(phi_smooth - phi_true)","#phi_{smooth}-#phi_{true} [rad];rad;Entries",200,-0.1,0.1);
  TH1D smooth_x_res("(x_smooth - x_true)","(x_{smooth}-x_{true}) [mm];mm;Entries",200,-2,2);
  TH1D smooth_y_res("(y_smooth - y_true)","(y_{smooth}-y_{true}) [mm];mm;Entries",200,-2,2);

  TH1D x_pull_step  ("x_pull_step","(x_s-x_true)/#sigma_{x};;Entries",200,-5,5);
  TH1D y_pull_step  ("y_pull_step","(y_s-y_true)/#sigma_{y};;Entries",200,-5,5);
  TH1D invR_pull_step("invR_pull_step","(invR_s-invR_t)/#sigma_{invR};;Entries",200,-5,5);
  TH1D tanL_pull_step("tanL_pull_step","(tanL_s-tanL_t)/#sigma_{tanL};;Entries",200,-5,5);
  TH1D phi_pull_step ("phi_pull_step","(phi_s-phi_t)/#sigma_{phi};;Entries",200,-5,5);

  //--- misura - true ---
  TH1D measx_res_step("(x_meas - x_true)","(x_{meas}-x_{true}) [mm];mm;Entries",200,-2,2);
  TH1D measy_res_step("(y_meas - y_true)","(y_{meas}-y_{true}) [mm];mm;Entries",200,-2,2);

  // --- momentum---
  TH1D p_res_step("p_res_step","p_{smooth} - p_{true} [MeV];MeV;Entries",200,-500,500);
  TH1D p_rel_res_step("(p_smooth - p_true)/p_true","(p_{s}-p_{t})/p_{t};;Entries",200,-1,1);
  
  TH2D p_true_vs_p_smooth_step("p_true_vs_p_smooth_step","True p vs smoothed p; p_{true}[MeV]; p_{reco}[MeV]",200,0,5000,200,0,5000);
  TH2D dp_vs_ptrue_step("dp_vs_ptrue_step","#Delta p vs p_{true};p_{true} [MeV];#Delta p [MeV]",200,0,5000,200,-500,500);
  TH2D dprel_vs_ptrue_step("dprel_vs_ptrue_step","#Delta p/p vs p_{true};p_{true} [MeV];(p_{s}-p_{t})/p_{t}",200,0,5000,200,-1,1);


  TH1D h_gpos("innovation_pos","Innovation (pos);;Entries",100,-3,3);
  TH1D h_gang("innovation_ang","Innovation (ang);;Entries",100,-3,3);
  TH1D chi2_h("chi2","chi2;;Entries",1000,0,1e5);
  

  // ---------------- plot per STEP ------------------------
  const Long64_t ns = steps->GetEntries(); //numero step
  for (Long64_t i=0;i<ns;++i){
    steps->GetEntry(i);

    // smoothed vs. true
    x_res_step.Fill(x_smooth - x_true);
    y_res_step.Fill(y_smooth - y_true);
    invR_res_step.Fill(invR_smooth - invR_true);
    tanL_res_step.Fill (tanL_smooth - tanL_true);
    phi_res_step.Fill  (phi_smooth  - phi_true);

    if (std::isfinite(sigma_x_smooth) && sigma_x_smooth>0) x_pull_step.Fill((x_smooth - x_true)/sigma_x_smooth);
    if (std::isfinite(sigma_y_smooth) && sigma_y_smooth>0) y_pull_step.Fill((y_smooth - y_true)/sigma_y_smooth);
    if (std::isfinite(sigma_invR_smooth) && sigma_invR_smooth>0) invR_pull_step.Fill((invR_smooth - invR_true)/sigma_invR_smooth);
    if (std::isfinite(sigma_tanL_smooth) && sigma_tanL_smooth>0) tanL_pull_step.Fill((tanL_smooth - tanL_true)/sigma_tanL_smooth);
    if (std::isfinite(sigma_phi_smooth)  && sigma_phi_smooth>0)  phi_pull_step.Fill((phi_smooth - phi_true)/sigma_phi_smooth);

    // measured vs. truth (per ora solo smearing di differenza)
    measx_res_step.Fill(x_meas - x_true);
    measy_res_step.Fill(y_meas - y_true);


    // momentum per step
    if (std::isfinite(p_true) && std::isfinite(p_smooth)) {
      p_res_step.Fill(p_smooth - p_true);
      p_true_vs_p_smooth_step.Fill(p_true, p_smooth);
    }
    
    if (std::isfinite(p_true) && p_true>0 && std::isfinite(p_smooth)) {
      const double dp = p_smooth - p_true;
      p_rel_res_step.Fill(dp / p_true);
      dp_vs_ptrue_step.Fill(p_true, dp);
      dprel_vs_ptrue_step.Fill(p_true, dp / p_true);
    }

    // innovation e chi2
    if (std::isfinite(innov_pos)) h_gpos.Fill(innov_pos);
    if (std::isfinite(innov_ang)) h_gang.Fill(innov_ang);
    chi2_h.Fill(chi2);

  }


// ------ SEED CHECKS------

std::map<std::tuple<int,int,int>, std::pair<int, Long64_t>> lastEntry;
for (Long64_t ie=0; ie<ns; ++ie) {
  steps->GetEntry(ie);
  auto key = std::make_tuple(run, event, track_id);
  auto it = lastEntry.find(key);
  if (it==lastEntry.end() || step_idx > it->second.first) {
    lastEntry[key] = { step_idx, ie };
  }
}

TH1D seed_dx   ("seed_dx","Seed: x_{s}-x_{t} [mm];mm;Entries",200,-2,2);
TH1D seed_dy   ("seed_dy","Seed: y_{s}-y_{t} [mm];mm;Entries",200,-2,2);
TH1D seed_dinvR("seed_dinvR","Seed: invR_{s}-invR_{t};;Entries",200,-1e-2,1e-2);
TH1D seed_dtanL("seed_dtanL","Seed: tan#lambda_{s}-tan#lambda_{t};;Entries",200,-0.1,0.1);
TH1D seed_dphi ("seed_dphi","Seed: #phi_{s}-#phi_{t} [rad];rad;Entries",200,-0.1,0.1);

TH1D seed_pullx("seed_pullx","Seed: (x_{s}-x_{t})/#sigma_{x};pull_{x};Entries",100,-5,5);
TH1D seed_pully("seed_pully","Seed: (y_{s}-y_{t})/#sigma_{y};pull_{y};Entries",100,-5,5);
TH1D seed_invR_pull("seed_invR_pull","Seed: (invR_{s}-invR_{t})/#sigma_{invR};;Entries",100,-5,5);
TH1D seed_tanL_pull("seed_tanL_pull","Seed: (tan#lambda_{s}-tan#lambda_{t})/#sigma_{tanL};;Entries",100,-5,5);
TH1D seed_phi_pull ("seed_phi_pull","Seed: (#phi_{s}-#phi_{t})/#sigma_{phi};;Entries",100,-5,5);

TH1D seed_dp("seed_dp","p_{s}^{last} - p_{t}^{first} [MeV];MeV;Entries",200,-1000,1000);
TH1D seed_dprel("seed_dprel","(p_{s}^{last}-p_{t}^{first})/p_{t}^{first};;Entries",200,-1,1);
TH1D seed_ratio_ps_over_pfirst("seed_ratio_ps_over_pfirst","p_{s}^{last}/p_{t}^{first};;Entries",200,0,2);

Int_t tr_run=0, tr_event=0, tr_tid=0, n_steps=0;
Double_t p_true_first=0, p_true_last=0, p_smooth_last=0;
tracks->SetBranchAddress("run",&tr_run);
tracks->SetBranchAddress("event",&tr_event);
tracks->SetBranchAddress("track_id",&tr_tid);
tracks->SetBranchAddress("n_steps",&n_steps);
tracks->SetBranchAddress("p_true_first",&p_true_first);

tracks->SetBranchAddress("p_true_last",&p_true_last);
tracks->SetBranchAddress("p_smooth_last",&p_smooth_last);


std::map<std::tuple<int,int,int>, double> pfirst_by_key;
const Long64_t nt = tracks->GetEntries();
for (Long64_t it=0; it<nt; ++it) {
  tracks->GetEntry(it);
  pfirst_by_key[ std::make_tuple(tr_run,tr_event,tr_tid) ] = p_true_first;
}

for (const auto& kv : lastEntry) {
  const auto key = kv.first;
  const Long64_t ie = kv.second.second;
  steps->GetEntry(ie);

  const double dx = x_smooth - x_true;
  const double dy = y_smooth - y_true;
  seed_dx.Fill(dx);
  seed_dy.Fill(dy);
  seed_dinvR.Fill(invR_smooth - invR_true);
  seed_dtanL.Fill (tanL_smooth - tanL_true);
  seed_dphi.Fill  (phi_smooth  - phi_true);

  if (std::isfinite(sigma_x_smooth) && sigma_x_smooth>0) seed_pullx.Fill(dx / sigma_x_smooth);
  if (std::isfinite(sigma_y_smooth) && sigma_y_smooth>0) seed_pully.Fill(dy / sigma_y_smooth);
  if (std::isfinite(sigma_invR_smooth) && sigma_invR_smooth>0) seed_invR_pull.Fill((invR_smooth - invR_true)/sigma_invR_smooth);
  if (std::isfinite(sigma_tanL_smooth) && sigma_tanL_smooth>0) seed_tanL_pull.Fill((tanL_smooth - tanL_true)/sigma_tanL_smooth);
  if (std::isfinite(sigma_phi_smooth)  && sigma_phi_smooth>0)  seed_phi_pull.Fill((phi_smooth - phi_true)/sigma_phi_smooth);


  const auto itp = pfirst_by_key.find(key);
  if (itp != pfirst_by_key.end() && std::isfinite(itp->second) && itp->second>0 && std::isfinite(p_smooth)) {
    const double p_first = itp->second;
    const double dp_first = p_smooth - p_first;
    seed_dp.Fill(dp_first);
    seed_dprel.Fill(dp_first / p_first);
    seed_ratio_ps_over_pfirst.Fill(p_smooth / p_first);
  }
}



  TFile fout(output_root.c_str(),"RECREATE");
  if (fout.IsZombie()) { std::cerr<<"Cannot create "<<output_root<<"\n"; return 3; }

    x_res_step.Write();
    y_res_step.Write();
    invR_res_step.Write();
    tanL_res_step.Write();
    phi_res_step.Write();

    measx_res_step.Write();
    measy_res_step.Write();

  p_res_step.Write();
  p_true_vs_p_smooth_step.Write();

  p_rel_res_step.Write();
  dp_vs_ptrue_step.Write();
  dprel_vs_ptrue_step.Write();

  h_gpos.Write();
  h_gang.Write();
  chi2_h.Write();

  seed_dx.Write();
  seed_dy.Write();
  seed_dinvR.Write();
  seed_dtanL.Write();
  seed_dphi.Write();
  seed_pullx.Write();
  seed_pully.Write();
  seed_invR_pull.Write();
  seed_tanL_pull.Write();
  seed_phi_pull.Write();

  seed_dp.Write();
  seed_dprel.Write();
  seed_ratio_ps_over_pfirst.Write();

  TDirectory* dGraphs = fout.mkdir("graphs");
  makeTrackGraph(steps, dGraphs);

  // ---- Istogrammi per track  ----
  TH2D p_init_vs_reco_smooth("p_init_vs_reco_smooth","Initial true p vs last smoothed reco p; p_{true}^{init} [MeV]; p_{reco}^{smooth,last} [MeV]",200,0,5000,200,0,5000);//???
  TH2D p_true_vs_reco_smooth("p_true_vs_reco_smooth","True p at last step vs last smoothed reco p; p_{true}^{last} [MeV]; p_{reco}^{smooth,last} [MeV]",200,0,5000,200,0,5000);
  TH2D p_diff_vs_points("p_diff_vs_points","DeltaP vs nPoints; nPoints; (p_{reco}^{last} - p_{true}^{last}) [MeV]",200,0,200,200,-500,500);

  for (Long64_t i=0;i<nt;++i){
    tracks->GetEntry(i);
    if (std::isfinite(p_true_first) && std::isfinite(p_smooth_last))
      p_init_vs_reco_smooth.Fill(p_true_first, p_smooth_last);
      
    if (std::isfinite(p_true_last) && std::isfinite(p_smooth_last)) {
      //p_first_vs_reco_smooth.Fill(p_true_last, p_smooth_last);
      p_diff_vs_points.Fill(n_steps, p_smooth_last - p_true_last);
    }
  }

  p_init_vs_reco_smooth.Write();
  //p_first_vs_reco_smooth.Write();
  p_diff_vs_points.Write();

  fout.Close();
  return 0;
}


int main(int argc, char** argv) {
  const char* in  = (argc > 1) ? argv[1] : "kf_trees.root";
  const char* out = (argc > 2) ? argv[2] : "kf_plots.root";
  return makePlot(in, out);
}