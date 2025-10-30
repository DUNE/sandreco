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

  Int_t track_id=0;
  step_idx=0 /*orientation=0*/;
  Double_t z=0;
  Double_t x_meas=0, y_meas=0;
  Double_t x_pred=0, y_pred=0;
  Double_t x_filt=0, y_filt=0;
  Double_t x_smooth=0, y_smooth=0;

  steps->SetBranchAddress("track_id", &track_id);
  steps->SetBranchAddress("step_idx", &step_idx);
  //steps->SetBranchAddress("orientation", &orientation);
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

  Double_t x_true=0, y_true=0;
  Double_t x_meas=0, y_meas=0;
  Double_t x_smooth=0, y_smooth=0;
  Double_t p_true=0, p_smooth=0;
  Double_t var_x=NAN, var_y=NAN;
  Double_t innov_pos=0, innov_ang=0, chi2=0;

  steps->SetBranchAddress("x_true", &x_true);
  steps->SetBranchAddress("y_true", &y_true);
  steps->SetBranchAddress("x_meas", &x_meas);
  steps->SetBranchAddress("y_meas", &y_meas);
  steps->SetBranchAddress("x_smooth", &x_smooth);
  steps->SetBranchAddress("y_smooth", &y_smooth);
  steps->SetBranchAddress("p_true", &p_true);
  steps->SetBranchAddress("p_smooth",&p_smooth);
  if (steps->GetBranch("var_x"))     steps->SetBranchAddress("var_x",&var_x);
  if (steps->GetBranch("var_y"))     steps->SetBranchAddress("var_y",&var_y);
  if (steps->GetBranch("innov_pos")) steps->SetBranchAddress("innov_pos",&innov_pos);
  if (steps->GetBranch("innov_ang")) steps->SetBranchAddress("innov_ang",&innov_ang);
  if (steps->GetBranch("chi2"))      steps->SetBranchAddress("chi2",&chi2);

  TH1D meas_x_res("(x_meas - x_true)","(x_{meas}-x_{true}) [mm];mm;Entries",200,-2,2);
  TH1D meas_y_res("(y_meas - y_true)","(y_{meas}-y_{true}) [mm];mm;Entries",200,-2,2);
  TH1D smooth_x_res("(x_smooth - x_true)","(x_{smooth}-x_{true}) [mm];mm;Entries",200,-2,2);
  TH1D smooth_y_res("(y_smooth - y_true)","(y_{smooth}-y_{true}) [mm];mm;Entries",200,-2,2);

  TH1D p_res_step("p_res_step","p_{smooth} - p_{true} [MeV];MeV;Entries",200,-500,500);
  TH2D p_true_vs_p_smooth_step("p_true_vs_p_smooth_step","True p vs smoothed p; p_{true}[MeV]; p_{reco}[MeV]",200,0,5000,200,0,5000);

  TH1D h_gpos("innovation_pos","Innovation (pos);;Entries",100,-3,3);
  TH1D h_gang("innovation_ang","Innovation (ang);;Entries",100,-3,3);
  TH1D x_pull("x_pull_smooth","(x_s-x_true)/#sigma_x;;Entries",200,-5,5);
  TH1D y_pull("y_pull_smooth","(y_s-y_true)/#sigma_y;;Entries",200,-5,5);
  TH1D chi2_h("chi2","chi2;;Entries",1000,0,1e5);
  TH1D h_pull_x_first("x_pull_first", "(x_s - x_true)/#sigma_{x} at first step; pull_{x}; entries", 100, -5, 5);
  TH1D h_pull_y_first("y_pull_first", "(y_s - y_true)/#sigma_{y} at first step; pull_{y}; entries", 100, -5, 5);
  
//--------PULL SEED------------
std::unordered_map<int, Long64_t> firstEntry;  
std::unordered_map<int, int>      firstStep; 

const Long64_t nSteps = steps->GetEntries();
for (Long64_t ie = 0; ie < nSteps; ++ie) {
  steps->GetEntry(ie);
  auto it = firstStep.find(trk);
  if (it == firstStep.end() || step < it->second) {
    firstStep[trk]  = step;
    firstEntry[trk] = ie;
  }
}

for (const auto& kv : firstEntry) {
  steps->GetEntry(kv.second);

  const double sigma_pos = 200E-6;
  
  if (std::isfinite(sigma_pos) && sigma_pos > 0) {
    const double pull_x = (x_smoothed - x_true) / sigma_pos;
    if (std::isfinite(pull_x)) h_pull_x_first.Fill(pull_x);
  }
  if (std::isfinite(sigma_pos) && sigma_pos > 0) {
    const double pull_y = (y_smoothed - y_true) / sigma_pos;
    if (std::isfinite(pull_y)) h_pull_y_first.Fill(pull_y);
  }
}

h_pull_x_first.Write();
h_pull_y_first.Write();



  // ---------------- plot per STEP ------------------------
  const Long64_t ns = steps->GetEntries(); //numero step
  for (Long64_t i=0;i<ns;++i){
    steps->GetEntry(i);

    // measured vs. truth (per ora solo smearing di differenza)
    meas_x_res.Fill(x_meas - x_true);
    meas_y_res.Fill(y_meas - y_true);

    // smoothed vs. true
    smooth_x_res.Fill(x_smooth - x_true);
    smooth_y_res.Fill(y_smooth - y_true);

    // momentum per step
    if (std::isfinite(p_true) && std::isfinite(p_smooth)) {
      p_res_step.Fill(p_smooth - p_true);
      p_true_vs_p_smooth_step.Fill(p_true, p_smooth);
    }

    // innovation e chi2
    if (std::isfinite(innov_pos)) h_gpos.Fill(innov_pos);
    if (std::isfinite(innov_ang)) h_gang.Fill(innov_ang);
    chi2_h.Fill(chi2);

    // pulls: var_x/var_y in m^2, pos in mm → dividere per 1e3??
    if (std::isfinite(var_x) && var_x>0) x_pull.Fill((x_smooth - x_true)/(std::sqrt(var_x)*1e3));
    if (std::isfinite(var_y) && var_y>0) y_pull.Fill((y_smooth - y_true)/(std::sqrt(var_y)*1e3));
  }


  TFile fout(output_root.c_str(),"RECREATE");
  if (fout.IsZombie()) { std::cerr<<"Cannot create "<<output_root<<"\n"; return 3; }

  meas_x_res.Write(); 
  meas_y_res.Write();
  smooth_x_res.Write(); 
  smooth_y_res.Write();
  p_res_step.Write();  
  p_true_vs_p_smooth_step.Write();
  h_gpos.Write(); 
  h_gang.Write();
  x_pull.Write(); 
  y_pull.Write();
  chi2_h.Write();

  TDirectory* dGraphs = fout.mkdir("graphs");
  makeTrackGraph(steps, dGraphs);

  // ---- Istogrammi per track  ----
  Int_t tr_track_id=0, tr_run=0, tr_event=0, n_steps=0;
  Double_t p_init_true=0;
  Double_t p_true_last=0;
  Double_t p_smooth_last=0;

  tracks->SetBranchAddress("track_id",&tr_track_id);
  tracks->SetBranchAddress("run",&tr_run);
  tracks->SetBranchAddress("event",&tr_event);
  tracks->SetBranchAddress("n_steps",&n_steps);
  tracks->SetBranchAddress("p_init_true",&p_init_true);
  tracks->SetBranchAddress("p_true_last",&p_true_last);
  tracks->SetBranchAddress("p_smooth_last",&p_smooth_last);

  TH2D p_init_vs_reco_smooth("p_init_vs_reco_smooth","Initial true p vs last smoothed reco p; p_{true}^{init} [MeV]; p_{reco}^{smooth,last} [MeV]",200,0,5000,200,0,5000);//???
  TH2D p_true_vs_reco_smooth("p_true_vs_reco_smooth","True p at last step vs last smoothed reco p; p_{true}^{last} [MeV]; p_{reco}^{smooth,last} [MeV]",200,0,5000,200,0,5000);
  TH2D p_diff_vs_points("p_diff_vs_points","DeltaP vs nPoints; nPoints; (p_{reco}^{last} - p_{true}^{last}) [MeV]",200,0,200,200,-500,500);

  const Long64_t nt = tracks->GetEntries();
  for (Long64_t i=0;i<nt;++i){
    tracks->GetEntry(i);
    if (std::isfinite(p_init_true) && std::isfinite(p_smooth_last))
      p_init_vs_reco_smooth.Fill(p_init_true, p_smooth_last);
      
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