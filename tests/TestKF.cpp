#include <TDatabasePDG.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <unordered_map>
#include <vector>

#include "SANDGeoManager.h"
#include "SANDKalmanFilter.h"
#include "SANDProcessTracklets.h"
#include "SANDTrackerClusterCollection.h"
#include "SANDTrackerDigitCollection.h"
#include "SANDTrackerUtils.h"
#include "SANDTrackletFinder.h"
#include "TryCompleteManager.h"
#include "utils.h"

void processEventWithKF(SANDGeoManager* sand_geo, TG4Event* mc_event,
                        std::vector<dg_wire>* digits, StepsTree* steps,
                        TracksTree* tracks, int run_number, int event_number)
{
  // unknown parameters
  int p[9] = {100, -2000, 2000, 100, -4000, -0, 100, 22500, 26000};

  // fill ID vs index map of DigitCollection
  sand_reco::tracker::DigitCollection::fillMap(digits);
  // get the vector of digits
  auto digit_vec = sand_reco::tracker::DigitCollection::getDigits();
  // check if the vector of digits is empty
  if (digit_vec.empty()) {
    return;
  }

  // get the name of the tracker
  std::string tracker_name =
      sand_reco::tracker::DigitCollection::getDigits().begin()->det;

  // build the clusters starting from digits and return them as a collection
  sand_reco::tracker::ClusterCollection clusters(
      sand_geo, digit_vec,
      sand_reco::tracker::ClusterCollection::ClusteringMethod::kCellAdjacency);

  // init the SAND Tracker Utils
  SANDTrackerUtils::init(sand_geo->getTGeoManager());

  TRandom3 rand(0);
  EDEPTree tree;
  tree.InizializeFromEdep(*mc_event, sand_geo->getTGeoManager());

  std::vector<EDEPTrajectory> primaryTrj;
  tree.Filter(
      std::back_insert_iterator<std::vector<EDEPTrajectory>>(primaryTrj),
      [](const EDEPTrajectory& trj) { return trj.GetParentId() == -1; });

  TDatabasePDG pdg_db;
  std::vector<SParticleInfo> particleInfos;

  std::vector<int> indeces;
  int ii = -1;
  for (auto trj : primaryTrj) {
    ii++;

    if (trj.GetHitMap().find(string_to_component[tracker_name]) ==
        trj.GetHitMap().end()) {
      continue;
    }

    if (trj.GetTrajectoryPoints().find(string_to_component[tracker_name]) ==
        trj.GetTrajectoryPoints().end()) {
      continue;
    }

    auto particle = pdg_db.GetParticle(trj.GetPDGCode());

    if (!particle) {
      continue;
    }

    if (particle->Mass() == 0 || particle->Charge() == 0) {
      continue;
    }

    SParticleInfo pi;
    pi.pdg_code = trj.GetPDGCode();
    pi.id = trj.GetId();
    pi.mass = particle->Mass();
    pi.charge = particle->Charge() / 3;

    double max_z = 0;
    bool to_be_reconstructed = false;

    for (auto& point :
         trj.GetTrajectoryPoints().at(string_to_component[tracker_name])) {
      if (point.GetPosition().Z() > max_z && point.GetMomentum().Z() > 100) {
        max_z = point.GetPosition().Z();
        pi.pos = point.GetPosition().Vect();
        pi.mom = point.GetMomentum();
        to_be_reconstructed = true;
      }
    }

    if (!to_be_reconstructed) continue;

    double sigma_pos = SANDTrackerUtils::getSigmaPositionMeasurement() * 1E3;
    double sigma_mom = 0.05;

    double x_smeared = rand.Gaus(pi.pos.X(), sigma_pos);
    double y_smeared = rand.Gaus(pi.pos.Y(), sigma_pos);
    double px_smeared = pi.mom.X() * rand.Gaus(1, sigma_mom);
    double py_smeared = pi.mom.Y() * rand.Gaus(1, sigma_mom);
    double pz_smeared = pi.mom.Z() * rand.Gaus(1, sigma_mom);

    pi.pos = TVector3(x_smeared, y_smeared, pi.pos.Z());
    pi.mom = TVector3(px_smeared, py_smeared, pz_smeared);
    pi.initial_pos = trj.GetTrajectoryPoints()
                         .at(string_to_component[tracker_name])[0]
                         .GetPosition()
                         .Vect();
    pi.initial_mom = trj.GetTrajectoryPoints()
                         .at(string_to_component[tracker_name])[0]
                         .GetMomentum();
    particleInfos.push_back(pi);
    indeces.push_back(ii);

    //---------------------------------------------------------------------------
    //  TrajectoryPoints option: obtains the measurements as
    //  the smearing of the true tracklet computed from trajectory points at
    //  defined steps.There is no clustering here, the z is ssociated directly
    //  to the closest avalaible MC-trajectory point.
    //---------------------------------------------------------------------------
    std::map<double, std::vector<Tracklet>> z_to_tracklets;
    auto points =
        trj.GetTrajectoryPoints().at(string_to_component[tracker_name]);
    const double step = 1.5;
    auto z_truth = z_to_truth(points, step);

    for (const auto& kv : z_truth) {
      const double z = kv.first;
      const Truth& t = kv.second;
      Tracklet measurements = makeMeasurementTrackletFromTruth(t, rand);
      z_to_tracklets[z].push_back(measurements);
    }

    // It keeps only one measurement per module, it has been used to test hypothesis on pull tests
    for (auto& kv : z_to_tracklets) {
      auto& vec = kv.second;
      if (vec.size() > 1) {
        int idx = rand.Integer(static_cast<int>(vec.size()));
        Tracklet keep = vec[idx];
        vec.clear();
        vec.push_back(keep);
      }
    }

    if (z_to_tracklets.empty()) continue;
    
    for (int ip = 0; ip < (int)indeces.size(); ++ip) {
      tryCompleteManager(z_to_tracklets, particleInfos[ip],
                         /*mg*/ nullptr, /*mgx*/ nullptr, steps, tracks,
                         run_number, event_number);
    }
  }

  int nParticles = particleInfos.size();

  if (nParticles == 0) {
    std::cerr << "no particles to be reconstructed...process aborted"
              << std::endl;
    return;
  }

}

//-----------------------------------------------------------------------------
//                                  MAIN 
//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  gStyle->SetOptStat(0);

  TFile f(argv[1], "READ");
  TGeoManager* geo = 0;
  geo = (TGeoManager*)f.Get("EDepSimGeometry");

  // MC info tree
  TTree* t_h = (TTree*)f.Get("EDepSimEvents");
  TG4Event* ev = new TG4Event;
  t_h->SetBranchAddress("Event", &ev);

  TFile f_d(argv[2], "READ");
  TTree* t = (TTree*)f_d.Get("tDigit");

  std::vector<dg_wire>* digits = 0;
  t->SetBranchAddress("dg_wire", &digits);

  bool plots = false;

  SANDGeoManager sand_geo;
  sand_geo.init(geo);

  std::string geometry;
  if (geo->FindVolumeFast("STTtracker_PV")) {
    geometry = "STT";
  } else if (geo->FindVolumeFast("SANDtracker_PV")) {
    geometry = "DRIFT";
  }
  sand_geo.fillAdjacentCells(geometry);

  StepsTree steps;
  TracksTree tracks;
  create_kf_trees("kf_steps.root", steps, tracks);

  int nev = t->GetEntries();
  for (int i = 0; i < nev; i++) {
    t_h->GetEntry(i);
    t->GetEntry(i);

    int run_number = 0;
    int event_number = i;

    processEventWithKF(&sand_geo, ev, digits, &steps, &tracks, run_number,
                       event_number);
  }

  close_kf_trees(steps, tracks);
}