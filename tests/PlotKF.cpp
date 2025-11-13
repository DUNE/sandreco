#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TTree.h>
#include <string>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

int makePlot(const std::string& input_root = "kf_trees.root",
             const std::string& output_root = "kf_plots.root");

struct PointInfos {
  int step{};
  double z{};
  double x_pred{}, y_pred{};      // predicted
  double x_filt{}, y_filt{};      // filtered
  double x_smooth{}, y_smooth{};  // smoothed
  double x_meas{}, y_meas{};      // measured (orientation??)
};

// ---------------- plot TRACCE (reco/smooth/filt/meas) ------------------------
static void makeTrackGraph(TTree* steps, TDirectory* outDir)
{

  Int_t run = 0, event = 0, track_id = 0, step_idx = 0;
  steps->SetBranchAddress("run", &run);
  steps->SetBranchAddress("event", &event);
  steps->SetBranchAddress("track_id", &track_id);
  steps->SetBranchAddress("step_idx", &step_idx);

  Double_t z = 0;
  Double_t x_meas = 0, y_meas = 0;
  Double_t x_pred = 0, y_pred = 0;
  Double_t x_filt = 0, y_filt = 0;
  Double_t x_smooth = 0, y_smooth = 0;

  steps->SetBranchAddress("z", &z);

  steps->SetBranchAddress("x_meas", &x_meas);
  steps->SetBranchAddress("y_meas", &y_meas);

  steps->SetBranchAddress("x_pred", &x_pred);
  steps->SetBranchAddress("y_pred", &y_pred);
  steps->SetBranchAddress("x_filt", &x_filt);
  steps->SetBranchAddress("y_filt", &y_filt);
  steps->SetBranchAddress("x_smooth", &x_smooth);
  steps->SetBranchAddress("y_smooth", &y_smooth);

  std::map<int, std::vector<PointInfos>> byTrack;
  const auto n = steps->GetEntries();
  for (auto i = 0; i < n; ++i) {
    steps->GetEntry(i);
    PointInfos p{};
    p.step = step_idx;
    p.z = z;
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
    std::sort(v.begin(), v.end(), [](const PointInfos& a, const PointInfos& b) {
      return a.step < b.step;
    });
    const int N = (int)v.size();

    auto* yz_predicted = new TGraph(N);
    auto* yz_filtered = new TGraph(N);
    auto* yz_smoothed = new TGraph(N);
    auto* yz_measured = new TGraph(N);

    auto* xz_predicted = new TGraph(N);
    auto* xz_filtered = new TGraph(N);
    auto* xz_smoothed = new TGraph(N);
    auto* xz_measured = new TGraph(N);

    for (int i = 0; i < N; ++i) {
      const auto& p = v[i];
      yz_predicted->SetPoint(i, p.z, p.y_pred);
      yz_filtered->SetPoint(i, p.z, p.y_filt);
      yz_smoothed->SetPoint(i, p.z, p.y_smooth);
      if (std::isfinite(p.y_meas))
        yz_measured->SetPoint(yz_measured->GetN(), p.z, p.y_meas);

      xz_predicted->SetPoint(i, p.z, p.x_pred);
      xz_filtered->SetPoint(i, p.z, p.x_filt);
      xz_smoothed->SetPoint(i, p.z, p.x_smooth);
      if (std::isfinite(p.x_meas))
        xz_measured->SetPoint(xz_measured->GetN(), p.z, p.x_meas);
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

int makePlot(const std::string& input_root, const std::string& output_root)
{
  TFile fin(input_root.c_str(), "READ");
  if (fin.IsZombie()) {
    std::cerr << "Cannot open " << input_root << "\n";
    return 1;
  }

  TTree* steps = (TTree*)fin.Get("steps");
  TTree* tracks = (TTree*)fin.Get("tracks");
  if (!steps || !tracks) {
    std::cerr << "Missing 'steps' or 'tracks' trees\n";
    return 2;
  }

  Int_t run = 0, event = 0, track_id = 0, step_idx = 0;
  Double_t z = 0.0;
  steps->SetBranchAddress("run", &run);
  steps->SetBranchAddress("event", &event);
  steps->SetBranchAddress("track_id", &track_id);
  steps->SetBranchAddress("step_idx", &step_idx);
  steps->SetBranchAddress("z", &z);

  //---true---
  Double_t x_true = 0.0;
  Double_t y_true = 0.0;
  Double_t invR_true = 0.0;
  Double_t tanL_true = 0.0;
  Double_t phi_true = 0.0;
  steps->SetBranchAddress("x_true", &x_true);
  steps->SetBranchAddress("y_true", &y_true);
  steps->SetBranchAddress("invR_true", &invR_true);
  steps->SetBranchAddress("tanL_true", &tanL_true);
  steps->SetBranchAddress("phi_true", &phi_true);

  //---smooth---
  Double_t x_smooth = 0.0;
  Double_t y_smooth = 0.0;
  Double_t invR_smooth = 0.0;
  Double_t tanL_smooth = 0.0;
  Double_t phi_smooth = 0.0;
  steps->SetBranchAddress("x_smooth", &x_smooth);
  steps->SetBranchAddress("y_smooth", &y_smooth);
  steps->SetBranchAddress("invR_smooth", &invR_smooth);
  steps->SetBranchAddress("tanL_smooth", &tanL_smooth);
  steps->SetBranchAddress("phi_smooth", &phi_smooth);
  //---sigma smoothed----
  Double_t sigma_x_smooth = 0.0;
  Double_t sigma_y_smooth = 0.0;
  Double_t sigma_invR_smooth = 0.0;
  Double_t sigma_tanL_smooth = 0.0;
  Double_t sigma_phi_smooth = 0.0;
  if (steps->GetBranch("sigma_x_smooth"))
    steps->SetBranchAddress("sigma_x_smooth", &sigma_x_smooth);
  if (steps->GetBranch("sigma_y_smooth"))
    steps->SetBranchAddress("sigma_y_smooth", &sigma_y_smooth);
  if (steps->GetBranch("sigma_invR_smooth"))
    steps->SetBranchAddress("sigma_invR_smooth", &sigma_invR_smooth);
  if (steps->GetBranch("sigma_tanL_smooth"))
    steps->SetBranchAddress("sigma_tanL_smooth", &sigma_tanL_smooth);
  if (steps->GetBranch("sigma_phi_smooth"))
    steps->SetBranchAddress("sigma_phi_smooth", &sigma_phi_smooth);
  //---measurements---
  Double_t x_meas = 0.0, y_meas = 0.0;
  Double_t x_pred = 0.0, y_pred = 0.0;
  Double_t x_filt = 0.0, y_filt = 0.0;
  steps->SetBranchAddress("x_meas", &x_meas);
  steps->SetBranchAddress("y_meas", &y_meas);
  steps->SetBranchAddress("x_pred", &x_pred);
  steps->SetBranchAddress("y_pred", &y_pred);
  steps->SetBranchAddress("x_filt", &x_filt);
  steps->SetBranchAddress("y_filt", &y_filt);

  Double_t p_true = 0.0;    // MeV
  Double_t p_smooth = 0.0;  // MeV
  Double_t innov_pos = 0.0, innov_ang = 0.0, chi2 = 0.0;
  steps->SetBranchAddress("p_true", &p_true);
  steps->SetBranchAddress("p_smooth", &p_smooth);
  if (steps->GetBranch("innov_pos"))
    steps->SetBranchAddress("innov_pos", &innov_pos);
  if (steps->GetBranch("innov_ang"))
    steps->SetBranchAddress("innov_ang", &innov_ang);
  if (steps->GetBranch("chi2")) steps->SetBranchAddress("chi2", &chi2);

  //---smoooth vs. true ---
  TH1D x_res_step(
      "(x_smooth - x_true)",
      "Residuals of reconstructed x;(x_{smooth}-x_{true}) [mm];Entries", 500,
      -1, 1);
  TH1D y_res_step(
      "(y_smooth - y_true)",
      "Residuals of reconstruced y;(y_{smooth}-y_{true}) [mm];Entries", 500, -1,
      1);
  TH1D invR_res_step(
      "(invR_smooth - invR_true)",
      "Residuals of reconstructed (1/R)_{smooth}-(1/R)_{true};Entries", 500,
      -0.1, 0.1);
  TH1D tanL_res_step("(tanL_smooth - tanL_true)",
                     "Residuals of reconstructed "
                     "tan#lambda;tan#lambda_{smooth}-tan#lambda_{true};Entries",
                     500, -20, 20);
  TH1D phi_res_step("(phi_smooth - phi_true)",
                    "Residuals of reconstructed #Phi;#phi_{smooth}-#phi_{true} "
                    "[mrad];Entries",
                    500, -20, 20);
  //--- pull test on smooth ---
  TH1D x_smooth_pull("x_smooth_pull",
                     "Pull test on reco x; #it{g(x)_{smooth}}; Entries", 500,
                     -10, 10);
  TH1D y_smooth_pull("y_smooth_pull",
                     "Pull test on reco y; #it{g(y)_{smooth}}; Entries", 500,
                     -10, 10);
  TH1D invR_smooth_pull("invR_smooth_pull",
                        "Pull test on reco 1/R; #it{g(1/R)_{smooth}};Entries",
                        500, -10, 10);
  TH1D tanL_smooth_pull(
      "tanL_smooth_pull",
      "Pull test on reco tan#lambda; #it{g(tan#lambda)_{smooth}}; Entries", 500,
      -5, 5);
  TH1D phi_smooth_pull("phi_smooth_pull",
                       "Pull test on reco #Phi; #it{g(#Phi)_{smooth}}; Entries",
                       500, -5, 5);
  //--- misura vs. true ---
  TH1D measx_res_step(
      "(x_meas - x_true)",
      "Residuals of measured x; (x_{meas}-x_{true}) [mm]; Entries", 500, -2, 2);
  TH1D measy_res_step(
      "(y_meas - y_true)",
      "Residuals of measured y; (y_{meas}-y_{true}) [mm]; Entries", 500, -2, 2);
  //--- pull tests on measurements ---
  TH1D pos_meas_pull("innovation_pos",
                     "Pull test on measured position; #it{g(pos)_{k}}; Entries",
                     100, -3, 3);
  TH1D ang_meas_pull(
      "innovation_ang",
      "Pull test on measured direction; #it{g(#theta)_{k}}; Entries", 100, -3,
      3);
  // DA AGGIUNGERE: residui sugli altri parametri di stato dalla misura e pull
  // misura su altri parametri stato

  TH1D chi2_h("chi2", "chi2;;Entries", 1000, 0, 1e5);

  // --- momentum---
  TH1D p_res_step("p_res_step",
                  "Residuals of reconstructed momentum; "
                  "p_{reco}-p_{true}/p_{true}; Entries",
                  500, -0.5, 0.5);
  TH2D h2_p_reco_vs_p_true("h2_p_reco_vs_p_true",
                           "p_{reco} vs p_{true};p_{true} [MeV];p_{reco} [MeV]",
                           200, 0, 3000, 200, 0, 3000);
  TH2D dp_vs_ptrue_step("dp_vs_ptrue_step",
                        "#Delta p vs p_{true};p_{true} [MeV];#Delta p [MeV]",
                        200, 0, 3000, 200, -1500, 1500);

  //--- sigma distribution ---
  TH1D h_sigma_x_smooth("sigma_x_smooth",
                        "#sigma_{x}^{smooth};#sigma_{x}^{smooth} [mm];Entries",
                        500, 0, 0.5);
  TH1D h_sigma_y_smooth("sigma_y_smooth",
                        "#sigma_{y}^{smooth};#sigma_{y}^{smooth} [mm];Entries",
                        500, 0, 0.5);
  TH1D h_sigma_invR_smooth(
      "sigma_invR_smooth",
      "#sigma_{1/R}^{smooth};#sigma_{1/R}^{smooth} [1/m];Entries", 500, 0,
      20);  // 1/[0.5, 0.5]
  TH1D h_sigma_tanL_smooth(
      "sigma_tanL_smooth",
      "#sigma_{tan#lambda}^{smooth};#sigma_{tan#lambda}^{smooth};Entries", 500,
      0, 0.01);
  TH1D h_sigma_phi_smooth(
      "sigma_phi_smooth",
      "#sigma_{#phi}^{smooth};#sigma_{#phi}^{smooth} [rad];Entries", 500, 0,
      0.05);

  // ---------------- plot per STEP ------------------------
  const Long64_t ns = steps->GetEntries();  // numero step
  for (Long64_t i = 0; i < ns; ++i) {
    steps->GetEntry(i);

    // smoothed vs. true
    x_res_step.Fill(x_smooth - x_true);
    y_res_step.Fill(y_smooth - y_true);
    invR_res_step.Fill(invR_smooth - invR_true);
    tanL_res_step.Fill((tanL_smooth - tanL_true) * 1000);
    phi_res_step.Fill((phi_smooth - phi_true) * 1000);
    // pull on smooth
    if (std::isfinite(sigma_x_smooth) && sigma_x_smooth > 0)
      x_smooth_pull.Fill((x_smooth - x_true) / (sigma_x_smooth));
    if (std::isfinite(sigma_y_smooth) && sigma_y_smooth > 0)
      y_smooth_pull.Fill((y_smooth - y_true) / (sigma_y_smooth));
    if (std::isfinite(sigma_invR_smooth) && sigma_invR_smooth > 0)
      invR_smooth_pull.Fill((invR_smooth - invR_true) / (sigma_invR_smooth));
    if (std::isfinite(sigma_tanL_smooth) && sigma_tanL_smooth > 0)
      tanL_smooth_pull.Fill((tanL_smooth - tanL_true) / (sigma_tanL_smooth));
    if (std::isfinite(sigma_phi_smooth) && sigma_phi_smooth > 0)
      phi_smooth_pull.Fill((phi_smooth - phi_true) / (sigma_phi_smooth));
    // measured vs. truth (per ora solo smearing di differenza)
    measx_res_step.Fill(x_meas - x_true);
    measy_res_step.Fill(y_meas - y_true);
    // pull on measurements
    if (std::isfinite(innov_pos)) pos_meas_pull.Fill(innov_pos);
    if (std::isfinite(innov_ang)) ang_meas_pull.Fill(innov_ang);

    chi2_h.Fill(chi2);

    // momentum per step
    if (std::isfinite(p_true) && std::isfinite(p_smooth)) {
      p_res_step.Fill((p_smooth - p_true) / p_true);
      h2_p_reco_vs_p_true.Fill(p_true, p_smooth);
    }

    if (std::isfinite(p_true) && p_true > 0 && std::isfinite(p_smooth)) {
      const double dp = p_smooth - p_true;
      dp_vs_ptrue_step.Fill(p_true, dp);
    }

    //--- sigma ---
    if (std::isfinite(sigma_x_smooth) && sigma_x_smooth > 0)
      h_sigma_x_smooth.Fill(sigma_x_smooth);
    if (std::isfinite(sigma_y_smooth) && sigma_y_smooth > 0)
      h_sigma_y_smooth.Fill(sigma_y_smooth);
    if (std::isfinite(sigma_invR_smooth) && sigma_invR_smooth > 0)
      h_sigma_invR_smooth.Fill((sigma_invR_smooth) * 1000);
    if (std::isfinite(sigma_tanL_smooth) && sigma_tanL_smooth > 0)
      h_sigma_tanL_smooth.Fill((sigma_tanL_smooth));
    if (std::isfinite(sigma_phi_smooth) && sigma_phi_smooth > 0)
      h_sigma_phi_smooth.Fill((sigma_phi_smooth));
  }

  // ---SEED----
  std::map<std::tuple<int, int, int>, std::pair<int, Long64_t>> lastEntry;
  for (Long64_t ie = 0; ie < ns; ++ie) {
    steps->GetEntry(ie);
    auto key = std::make_tuple(run, event, track_id);
    auto it = lastEntry.find(key);
    if (it == lastEntry.end() || step_idx > it->second.first) {
      lastEntry[key] = {step_idx, ie};
    }
  }

  TH1D seed_dx("seed_dx",
               "Residuals of x at seed;x_{smooth}-x_{true} [mm];Entries", 200,
               -2, 2);
  TH1D seed_dy("seed_dy",
               "Residuals of y at seed;y_{smooth}-y_{true} [mm];Entries", 200,
               -2, 2);
  TH1D seed_dinvR("seed_dinvR",
                  "Residuals of 1/R at seed;1/R_{smooth}-1/R_{true};Entries",
                  200, -0.2, 0.2);
  TH1D seed_dtanL("seed_dtanL",
                  "Residuals of tan#lambda at "
                  "seed;tan#lambda_{smooth}-tan#lambda_{true};Entries",
                  200, -0.02, 0.02);
  TH1D seed_dphi(
      "seed_dphi",
      "Residuals of #Phi at seed;#phi_{smooth}-#phi_{true} [rad];Entries", 200,
      -0.04, 0.04);

  TH1D x_seed_pull(
      "x_seed_pull",
      "Pull test on reco x at seeding point ; #it{g(x)_{seed}}; Entries", 500,
      -10, 10);
  TH1D y_seed_pull(
      "y_seed_pull",
      "Pull test on reco y at seeding point ; #it{g(y)_{seed}}; Entries", 500,
      -10, 10);
  TH1D invR_seed_pull(
      "invR_seed_pull",
      "Pull test on reco 1/R at seeding point ; #it{g(1/R)_{seed}};Entries",
      500, -10, 10);
  TH1D tanL_seed_pull("tanL_seed_pull",
                      "Pull test on reco tan#lambda at seeding point ; "
                      "#it{g(tan#lambda)_{seed}}; Entries",
                      500, -10, 10);
  TH1D phi_seed_pull(
      "phi_seed_pull",
      "Pull test on reco #Phi at seeding point ; #it{g(#Phi)_{seed}}; Entries",
      500, -10, 10);

  //---TRACK---
  Int_t tr_run = 0, tr_event = 0, tr_tid = 0, n_steps = 0;
  Double_t p_true_first = 0, p_true_last = 0, p_smooth_last = 0;
  tracks->SetBranchAddress("run", &tr_run);
  tracks->SetBranchAddress("event", &tr_event);
  tracks->SetBranchAddress("track_id", &tr_tid);
  tracks->SetBranchAddress("n_steps", &n_steps);
  tracks->SetBranchAddress("p_true_first", &p_true_first);
  tracks->SetBranchAddress("p_true_last", &p_true_last);
  tracks->SetBranchAddress("p_smooth_last", &p_smooth_last);

  std::map<std::tuple<int, int, int>, double> pfirst_by_key;
  const Long64_t nt = tracks->GetEntries();
  for (Long64_t it = 0; it < nt; ++it) {
    tracks->GetEntry(it);
    pfirst_by_key[std::make_tuple(tr_run, tr_event, tr_tid)] = p_true_first;
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
    seed_dtanL.Fill(tanL_smooth - tanL_true);
    seed_dphi.Fill(phi_smooth - phi_true);

    if (std::isfinite(sigma_x_smooth) && sigma_x_smooth > 0)
      x_seed_pull.Fill((x_smooth - x_true) / (sigma_x_smooth));
    if (std::isfinite(sigma_y_smooth) && sigma_y_smooth > 0)
      y_seed_pull.Fill((y_smooth - y_true) / (sigma_y_smooth));
    if (std::isfinite(sigma_invR_smooth) && sigma_invR_smooth > 0)
      invR_seed_pull.Fill((invR_smooth - invR_true) / (sigma_invR_smooth));
    if (std::isfinite(sigma_tanL_smooth) && sigma_tanL_smooth > 0)
      tanL_seed_pull.Fill((tanL_smooth - tanL_true) / (sigma_tanL_smooth));
    if (std::isfinite(sigma_phi_smooth) && sigma_phi_smooth > 0)
      phi_seed_pull.Fill((phi_smooth - phi_true) / (sigma_phi_smooth));

    // const auto itp = pfirst_by_key_pull.find(key);
    // if (itp != pfirst_by_key.end() && std::isfinite(itp->second) &&
    // itp->second>0 && std::isfinite(p_smooth)) {
    //   const double p_first = itp->second;
    //   const double dp_first = p_smooth - p_first;
    //   seed_dp.Fill(dp_first);
    //   seed_dprel.Fill(dp_first / p_first);
    //   seed_ratio_ps_over_pfirst.Fill(p_smooth / p_first);
    // }
  }

  TFile fout(output_root.c_str(), "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "Cannot create " << output_root << "\n";
    return 3;
  }

  x_res_step.Write();
  y_res_step.Write();
  invR_res_step.Write();
  tanL_res_step.Write();
  phi_res_step.Write();

  x_smooth_pull.Write();
  y_smooth_pull.Write();
  invR_smooth_pull.Write();
  tanL_smooth_pull.Write();
  phi_smooth_pull.Write();

  measx_res_step.Write();
  measy_res_step.Write();

  pos_meas_pull.Write();
  ang_meas_pull.Write();
  chi2_h.Write();

  p_res_step.Write();
  h2_p_reco_vs_p_true.Write();
  dp_vs_ptrue_step.Write();

  h_sigma_x_smooth.Write();
  h_sigma_y_smooth.Write();
  h_sigma_invR_smooth.Write();
  h_sigma_tanL_smooth.Write();
  h_sigma_phi_smooth.Write();

  seed_dx.Write();
  seed_dy.Write();
  seed_dinvR.Write();
  seed_dtanL.Write();
  seed_dphi.Write();

  x_seed_pull.Write();
  y_seed_pull.Write();
  invR_seed_pull.Write();
  tanL_seed_pull.Write();
  phi_seed_pull.Write();

  // seed_dp.Write();
  // seed_dprel.Write();
  // seed_ratio_ps_over_pfirst.Write();

  TDirectory* dGraphs = fout.mkdir("graphs");
  // makeTrackGraph(steps, dGraphs);

  // ---- plot per TRACK  ----
  TH2D p_init_vs_reco_smooth(
      "p_init_vs_reco_smooth",
      "Initial true p vs last smoothed reco p; p_{true}^{init} [MeV]; "
      "p_{reco}^{smooth,last} [MeV]",
      200, 0, 3000, 200, 0, 3000);  //???
  TH2D p_res_last("p_res_last",
                  "True p at last step vs last smoothed reco p; "
                  "p_{true}^{last} [MeV]; p_{reco}^{smooth,last} [MeV]",
                  200, 0, 3000, 200, 0, 3000);
  TH2D p_diff_vs_points(
      "p_diff_vs_points",
      "DeltaP vs nPoints; nPoints; (p_{reco}^{last} - p_{true}^{last}) [MeV]",
      200, 0, 200, 200, -500, 500);

  for (Long64_t i = 0; i < nt; ++i) {
    tracks->GetEntry(i);
    if (std::isfinite(p_true_first) && std::isfinite(p_smooth_last))
      p_init_vs_reco_smooth.Fill(p_true_first, p_smooth_last);

    if (std::isfinite(p_true_last) && std::isfinite(p_smooth_last)) {
      p_res_last.Fill(p_true_last, p_smooth_last);
      p_diff_vs_points.Fill(n_steps, p_smooth_last - p_true_last);
    }
  }

  p_init_vs_reco_smooth.Write();
  p_res_last.Write();
  p_diff_vs_points.Write();

  fout.Close();
  return 0;
}

int main(int argc, char** argv)
{
  const char* in = (argc > 1) ? argv[1] : "kf_trees.root";
  const char* out = (argc > 2) ? argv[2] : "kf_plots.root";
  return makePlot(in, out);
}