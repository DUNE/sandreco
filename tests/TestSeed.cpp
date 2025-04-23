#include <TVector3.h>
#include <TTree.h>
#include <TFile.h>
#include <TMarker.h>
#include <TArrow.h>
#include <TDatabasePDG.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <random>

#include "SANDGeoManager.h"
#include "SANDTrackletFinder.h"
#include "SANDTrackerClusterCollection.h"
#include "SANDTrackerDigitCollection.h"
#include "SANDKalmanFilter.h"
#include "utils.h"

#include "EDEPTree.h"

// Generate helix trajectory points
sand_reco::kf::TrackletMap GenerateHelixZY(
  double z0, double y0, double x0,
  double invR, double tanLambda,
  double phi0,
  int nPoints,
  double stepX,
  double sx,
  double sy  // distance between x samples
) {

  double R = 1.0 / invR;
  double omega = invR;  // curvature = 1/R = delta(phi)/ds
  TRandom3 randGen(0);
  sand_reco::kf::TrackletMap three_tracklets;

  for (int i = 0; i < nPoints; ++i) {
      double x = x0 + i * stepX;
      double s = (x - x0) / tanLambda;  // arc length along helix

      double phi = phi0 + omega * s;  // angle swept
      double z = z0 + R * ( std::cos(phi) -  std::cos(phi0));
      double y = y0 + R * ( std::sin(phi) -  std::sin(phi0));

      TVectorD trklet(8);
      trklet[0] = x+randGen.Gaus(0, sx);
      trklet[1] = y+randGen.Gaus(0, sy);
      std::vector<TVectorD> trklet_vec;
      trklet_vec.push_back(trklet);
      three_tracklets[z]= trklet_vec;
  }

  return three_tracklets;
}

sand_reco::kf::TrackletMap Find3fromTrajectory(std::vector<EDEPTrajectoryPoint> trj_points, double sx = 0.0, double sy = 0.0) {

  sand_reco::kf::TrackletMap three_tracklets;
  TRandom3 randGen(0);

  TVectorD trklet(8);
  trklet[0] = trj_points[0].GetPosition().X()+randGen.Gaus(0, sx*1E3);
  trklet[1] = trj_points[0].GetPosition().Y()+randGen.Gaus(0, sy*1E3);
  std::vector<TVectorD> trklet_vec;
  trklet_vec.push_back(trklet);
  three_tracklets[trj_points[0].GetPosition().Z()]= trklet_vec;

  TVectorD trklet2(8);
  trklet2[0] = trj_points[trj_points.size()/2].GetPosition().X()+randGen.Gaus(0, sx*1E3);
  trklet2[1] = trj_points[trj_points.size()/2].GetPosition().Y()+randGen.Gaus(0, sy*1E3);
  std::vector<TVectorD> trklet_vec2;
  trklet_vec2.push_back(trklet2);
  three_tracklets[trj_points[trj_points.size()/2].GetPosition().Z()]= trklet_vec2;

  TVectorD trklet3(8);
  trklet3[0] = trj_points[trj_points.size()-1].GetPosition().X()+randGen.Gaus(0, sx*1E3);
  trklet3[1] = trj_points[trj_points.size()-1].GetPosition().Y()+randGen.Gaus(0, sy*1E3);
  std::vector<TVectorD> trklet_vec3;
  trklet_vec3.push_back(trklet3);
  three_tracklets[trj_points[trj_points.size()-1].GetPosition().Z()]= trklet_vec3;

  return three_tracklets;
  
}

void FlattenMatrix(const TMatrixD* matrix, std::vector<double>* output) {
  if (!matrix || !output) return;

  int nRows = matrix->GetNrows();
  int nCols = matrix->GetNcols();
  output->clear();
  output->reserve(nRows * nCols);

  for (int i = 0; i < nRows; ++i) {
      for (int j = 0; j < nCols; ++j) {
          output->push_back((*matrix)(i, j));
      }
  }
}

void trySeedManager(sand_reco::kf::TrackletMap z_to_tracklets, 
                    SParticleInfo particleInfo, 
                    std::vector<EDEPTrajectoryPoint> trj_points,
                    TMatrixD & StateVectorMC,
                    TMatrixD & StateCovMC,
                    TMatrixD & StateVectorSeed,
                    TMatrixD & StateCovSeed) {

  double sx = SANDTrackerUtils::getSigmaPositionMeasurement();
  double sy = SANDTrackerUtils::getSigmaPositionMeasurement();
  sand_reco::kf::Manager managerSeed;
  //auto closest= managerSeed.FindSeedPoints_MCstart(&z_to_tracklets, particleInfo, 200);
  auto closest = Find3fromTrajectory(trj_points,sx,sy);
  for (auto el:closest) {
    std::cout << "Z: " << el.first*1E-3 << std::endl;
    for (auto el2:el.second) {
      std::cout << "X: " << el2[0]*1E-3 << " Y: " << el2[1]*1E-3 << std::endl;
    }
  }
  
  managerSeed.initFromSeed(&closest,&z_to_tracklets, particleInfo, sx, sy);
  

  sand_reco::kf::Manager managerMC;
  managerMC.initFromMC(&z_to_tracklets, particleInfo);
  StateVectorMC = managerMC.getTrack().getStep(0).getStage(sand_reco::kf::TrackStep::TrackStateStage::kFiltering).getStateVector()();
  StateCovMC = managerMC.getTrack().getStep(0).getStage(sand_reco::kf::TrackStep::TrackStateStage::kFiltering).getStateCovMatrix();
  StateVectorSeed = managerSeed.getTrack().getStep(0).getStage(sand_reco::kf::TrackStep::TrackStateStage::kFiltering).getStateVector()();
  StateCovSeed = managerSeed.getTrack().getStep(0).getStage(sand_reco::kf::TrackStep::TrackStateStage::kFiltering).getStateCovMatrix();

  auto dir=-1.0;
  auto simple_helix = GenerateHelixZY(
  particleInfo.pos.Z()*1E-3, StateVectorMC[1][0], StateVectorMC[0][0],
  dir*StateVectorMC[2][0], dir*StateVectorMC[3][0],
  StateVectorMC[4][0],
  3,  // number of points
  dir*0.6,    // stepX
  0,0);

  for (auto el:simple_helix) {
    std::cout << "Z: " << el.first << std::endl;
    for (auto el2:el.second) {
      std::cout << "X: " << el2[0] << " Y: " << el2[1] << std::endl;
    }
  }
  
  return;
}

void processEventWithSeed(SANDGeoManager* sand_geo, 
                          TG4Event* mc_event, 
                          std::vector<dg_wire>* digits,
                          TMatrixD & StateVectorMC,
                          TMatrixD & StateCovMC,
                          TMatrixD & StateVectorSeed,
                          TMatrixD & StateCovSeed)
{
  
  int p[9] = {100, -2000, 2000, 100, -4000, -1000, 100, 23800, 26000};

  sand_reco::tracker::DigitCollection::fillMap(digits);
  auto digit_map =  sand_reco::tracker::DigitCollection::getDigits();
  if (sand_reco::tracker::DigitCollection::getDigits().empty()) {
    return;
  }
  std::string tracker_name = sand_reco::tracker::DigitCollection::getDigits().begin()->det;
  sand_reco::tracker::ClusterCollection clusters(sand_geo, sand_reco::tracker::DigitCollection::getDigits(), sand_reco::tracker::ClusterCollection::ClusteringMethod::kCellAdjacency);
  
  TrackletFinder traklet_finder;
  traklet_finder.setVolumeParameters(p);
  traklet_finder.setSigmaPosition(0.2);
  traklet_finder.setSigmaAngle(0.2);

  std::map<double, std::vector<TVectorD>> z_to_tracklets;

  SANDTrackerUtils::init(sand_geo->getTGeoManager());

  int gg = 0;
  for (const auto& container:clusters.getContainers()) {
    for (const auto& cluster_in_container:container->getClusters()) {

      traklet_finder.setCells(cluster_in_container);
      auto minima = traklet_finder.findTracklets();
      double z_start = cluster_in_container.getZ();
      for (uint trk = 0; trk < minima.size(); trk++) {
        if (minima[trk][4] < 1E-2) {
          z_to_tracklets[cluster_in_container.getZ()].push_back(minima[trk]);
        }
      }
      traklet_finder.clear();
    }
  }
  
  int sum = 0;
  for (auto el:z_to_tracklets) {
    sum += el.second.size();
  }
  if (sum == 0) {
    return;
  }

  EDEPTree tree;
  tree.InizializeFromEdep(*mc_event, sand_geo->getTGeoManager());
  
  std::vector<EDEPTrajectory> primaryTrj;
  tree.Filter(std::back_insert_iterator<std::vector<EDEPTrajectory>>(primaryTrj), 
    [](const EDEPTrajectory& trj) { return trj.GetParentId() == -1;} );

  TDatabasePDG pdg_db;
  std::vector<SParticleInfo> particleInfos;
  std::vector<std::vector<EDEPTrajectoryPoint>> trj_points;
  for (auto trj:primaryTrj) {
    SParticleInfo pi;
    pi.pdg_code = trj.GetPDGCode();
    pi.id       = trj.GetId();
    auto particle = pdg_db.GetParticle(pi.pdg_code);
    if (!particle) continue;
    pi.mass = particle->Mass();
    pi.charge = particle->Charge() / 3;

    pi.pos = trj.GetTrajectoryPoints().at(string_to_component[tracker_name]).back().GetPosition().Vect();
    pi.mom = trj.GetTrajectoryPoints().at(string_to_component[tracker_name]).back().GetMomentum();
    particleInfos.push_back(pi);

    auto points = trj.GetTrajectoryPoints().at(string_to_component[tracker_name]);
    trj_points.push_back(points);

    std::cout << "Initial Momentum " << trj.GetInitialMomentum().Vect().Mag() << std::endl;
  }

  int nParticles = particleInfos.size();
  
  if (nParticles == 0) {
    std::cerr << "no particles to be reconstructed...process aborted"
              << std::endl;
    return;
  }

  for (int ip = 0; ip < nParticles; ip++) {
    trySeedManager(z_to_tracklets, particleInfos[ip], trj_points[ip],
                   StateVectorMC, StateCovMC,
                   StateVectorSeed, StateCovSeed);
    std::cout << "Number of trajectory points inside the tracker: " << trj_points[ip].size() << std::endl;
  }
}

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

    
  TFile* h_out = new TFile("seed_out.root", "RECREATE");

  // Create a TTree
  TTree *tree = new TTree("MatrixTree", "Tree with TMatrixD branches");

  // Create pointers to TMatrixD objects
  TMatrixD *StateVectorMC_mat   = new TMatrixD(5, 1);  // Example: 5x1 vector
  TMatrixD *StateVectorSeed_mat = new TMatrixD(5, 1);
  TMatrixD *StateCovMC_mat      = new TMatrixD(5, 5);  // Example: 5x5 covariance
  TMatrixD *StateCovSeed_mat    = new TMatrixD(5, 5);

  std::vector<double>* StateVectorMC = new std::vector<double>(5);
  std::vector<double>* StateVectorSeed = new std::vector<double>(5);
  std::vector<double>* StateCovMC = new std::vector<double>(25);
  std::vector<double>* StateCovSeed = new std::vector<double>(25);

  // Create branches
  tree->Branch("StateVectorMC", &StateVectorMC);
  tree->Branch("StateVectorSeed", &StateVectorSeed);
  tree->Branch("StateCovMC", &StateCovMC);
  tree->Branch("StateCovSeed", &StateCovSeed);

  SANDGeoManager sand_geo;
  sand_geo.init(geo);

  for (int i = 0; i < 1; i++) {
    t_h->GetEntry(i);
    t->GetEntry(i);

    processEventWithSeed(&sand_geo, ev, digits,
                         *StateVectorMC_mat, *StateCovMC_mat,
                         *StateVectorSeed_mat, *StateCovSeed_mat);

    FlattenMatrix(StateVectorMC_mat, StateVectorMC);
    FlattenMatrix(StateVectorSeed_mat, StateVectorSeed);
    FlattenMatrix(StateCovMC_mat, StateCovMC);
    FlattenMatrix(StateCovSeed_mat, StateCovSeed);
    // Fill the tree with the current TMatrixD objects
    tree->Fill();

    StateVectorMC->clear();
    StateVectorSeed->clear();
    StateCovMC->clear();
    StateCovSeed->clear();


  }

  tree->Write();
  h_out->Close();

  std::cout << "Tree saved to tree_with_matrices.root" << std::endl;

  // Clean up
  delete StateVectorMC_mat;
  delete StateVectorSeed_mat;
  delete StateCovMC_mat;
  delete StateCovSeed_mat;
  delete StateVectorMC;
  delete StateVectorSeed;
  delete StateCovMC;
  delete StateCovSeed;

}