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

sand_reco::kf::TrackletMap Find3fromTrajectory(std::vector<EDEPTrajectoryPoint> trj_points){

  sand_reco::kf::TrackletMap three_tracklets;

  TVectorD trklet(8);
  trklet[0] = trj_points[0].GetPosition().X();
  trklet[1] = trj_points[0].GetPosition().Y();
  std::vector<TVectorD> trklet_vec;
  trklet_vec.push_back(trklet);
  three_tracklets[trj_points[0].GetPosition().Z()]= trklet_vec;

  TVectorD trklet2(8);
  trklet2[0] = trj_points[trj_points.size()/2].GetPosition().X();
  trklet2[1] = trj_points[trj_points.size()/2].GetPosition().Y();
  std::vector<TVectorD> trklet_vec2;
  trklet_vec2.push_back(trklet2);
  three_tracklets[trj_points[trj_points.size()/2].GetPosition().Z()]= trklet_vec2;

  TVectorD trklet3(8);
  trklet3[0] = trj_points[trj_points.size()-1].GetPosition().X();
  trklet3[1] = trj_points[trj_points.size()-1].GetPosition().Y();
  std::vector<TVectorD> trklet_vec3;
  trklet_vec3.push_back(trklet3);
  three_tracklets[trj_points[trj_points.size()-1].GetPosition().Z()]= trklet_vec3;

  return three_tracklets;
  
}

void trySeedManager(sand_reco::kf::TrackletMap z_to_tracklets, SParticleInfo particleInfo, std::vector<EDEPTrajectoryPoint> trj_points) {
  sand_reco::kf::Manager managerSeed;
  //auto closest= managerSeed.FindSeedPoints_MCstart(&z_to_tracklets, particleInfo, 200);
  auto closest = Find3fromTrajectory(trj_points);
  for (auto el:closest) {
    std::cout << "Z: " << el.first << std::endl;
    for (auto el2:el.second) {
      std::cout << "X: " << el2[0] << " Y: " << el2[1] << std::endl;
    }
  }
  
  managerSeed.initFromSeed(&closest,&z_to_tracklets, particleInfo);

  sand_reco::kf::Manager managerMC;
  managerMC.initFromMC(&z_to_tracklets, particleInfo);
  return;
}

void processEventWithSeed(SANDGeoManager* sand_geo, TG4Event* mc_event, std::vector<dg_wire>* digits)
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
    trySeedManager(z_to_tracklets, particleInfos[ip], trj_points[ip]);
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

    
  TFile* h_out = new TFile("h_out.root", "RECREATE");

  SANDGeoManager sand_geo;
  sand_geo.init(geo);

  for (int i = 0; i < 20; i++) {
    t_h->GetEntry(i);
    t->GetEntry(i);

    processEventWithSeed(&sand_geo, ev, digits);


  }

}