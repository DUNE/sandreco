#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"

#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
  gStyle->SetOptStat(1110);

  std::string filename = "triplets.root";
  if (argc > 1) {
    filename = argv[1];
  }

  TFile* f = TFile::Open(filename.c_str(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "[PlotTriplets] Error opening file " << filename << std::endl;
    return 1;
  }

  TTree* t = dynamic_cast<TTree*>(f->Get("triplets"));
  if (!t) {
    std::cerr << "[PlotTriplets] Cannot find TTree 'triplets' in "
              << filename << std::endl;
    return 1;
  }

  int    run = 0;
  int    event = 0;
  int    moduleID = 0;

  double x_true = 0.0;
  double y_true = 0.0;
  double z_true = 0.0;
  double x      = 0.0;
  double y      = 0.0;
  double z      = 0.0;

  double score  = 0.0;
  double D      = 0.0;

  double dx     = 0.0;
  double dy     = 0.0;

  int nU = 0;
  int nV = 0;
  int nY = 0;

  t->SetBranchAddress("run",      &run);
  t->SetBranchAddress("event",    &event);
  t->SetBranchAddress("moduleID", &moduleID);

  t->SetBranchAddress("x_true",   &x_true);
  t->SetBranchAddress("y_true",   &y_true);
  t->SetBranchAddress("z_true",   &z_true);

  t->SetBranchAddress("z",        &z);
  t->SetBranchAddress("x",        &x);
  t->SetBranchAddress("y",        &y);

  t->SetBranchAddress("score",    &score);
  t->SetBranchAddress("D",        &D);

  t->SetBranchAddress("dx",       &dx);
  t->SetBranchAddress("dy",       &dy);

  t->SetBranchAddress("nU",       &nU);
  t->SetBranchAddress("nV",       &nV);
  t->SetBranchAddress("nY",       &nY);


  // NB: score = |D|, sempre >= 0; D è con segno; dx,dy residui in mm
  TH1D hScore("hScore", "Triplet score;score = |D|;entries", 100, 0.0, 10.0);
  TH1D hD    ("hD",     "D = U+V-2Ycos#theta;D;entries",      100, -10.0, 10.0);
  TH1D hDx   ("hDx",    "x_{reco} - x_{true};#Delta x [mm];entries", 100, -100.0, 100.0);
  TH1D hDy   ("hDy",    "y_{reco} - y_{true};#Delta y [mm];entries", 100, -100.0, 100.0);

  Long64_t nentries = t->GetEntries();
  std::cout << "[PlotTriplets] Entries in tree: " << nentries << std::endl;

  for (Long64_t i = 0; i < nentries; ++i) {
    t->GetEntry(i);

    hScore.Fill(score);
    hD.Fill(D);

    hDx.Fill(dx);
    hDy.Fill(dy);

  }


  TFile fout("triplets_plots.root", "RECREATE");
  hScore.Write();
  hD.Write();
  hDx.Write();
  hDy.Write();
  fout.Close();

  return 0;
}
