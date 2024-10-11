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
#include "STTKFKalmanFilter.h"
#include "utils.h"

#include "EDEPTree.h"

struct SParticleInfo {
  int charge;
  double mass;
  int pdg_code;
  int id;

  TVector3 pos;
  TVector3 mom;
};

void ProcessEventWithKF(SANDGeoManager* sand_geo, TG4Event* mc_event, std::vector<dg_wire>* digits)
{
  
  int p[9] = {100, -2000, 2000, 100, -3200, -2300, 100, 23800, 26000};

  SANDTrackerDigitCollection::FillMap(digits);
  auto digit_map =  SANDTrackerDigitCollection::GetDigits();
  std::string tracker_name = SANDTrackerDigitCollection::GetDigits().begin()->det;
  // std::cout << "TRACKER: " << tracker_name << std::endl;

  SANDTrackerClusterCollection clusters(sand_geo, SANDTrackerDigitCollection::GetDigits(), SANDTrackerClusterCollection::ClusteringMethod::kCellAdjacency);
  
  TrackletFinder traklet_finder;
  traklet_finder.SetVolumeParameters(p);
  traklet_finder.SetSigmaPosition(0.2);
  traklet_finder.SetSigmaAngle(0.2);

  std::map<double, std::vector<TVectorD>> z_to_tracklets;

  SANDTrackerUtils::Init(sand_geo->GetTGeoManager());
  STTKFGeoManager::Init();

  int gg = 0;
  for (const auto& container:clusters.GetContainers()) {
    for (const auto& cluster_in_container:container->GetClusters()) {
      // std::cout << 100 * gg / container->GetClusters().size() << std::endl;
      // gg++;
      // if (gg == 100) break;
      // if(cluster_in_container.GetZ() < 25650) continue;

      traklet_finder.SetCells(cluster_in_container);
      auto minima = traklet_finder.FindTracklets();
      double z_start = cluster_in_container.GetZ();
      for (uint trk = 0; trk < minima.size(); trk++) {
        if (minima[trk][4] < 1E-2) {
          z_to_tracklets[cluster_in_container.GetZ()].push_back(minima[trk]);
        }
      }
      traklet_finder.Clear();
    }
  }
  
  int sum = 0;
  for (auto el:z_to_tracklets) {
    // std::cout << "At z = " << el.first << " there are " << el.second.size() << " tracklets" << std::endl;
    sum += el.second.size();
  }
  std::cout << "Total tracklets: " << sum << std::endl;
  if (sum == 0) {
    return;
  }

  EDEPTree tree;
  tree.InizializeFromEdep(*mc_event, sand_geo->GetTGeoManager());
  
  std::vector<EDEPTrajectory> primaryTrj;
  tree.Filter(std::back_insert_iterator<std::vector<EDEPTrajectory>>(primaryTrj), 
    [](const EDEPTrajectory& trj) { return trj.GetParentId() == -1;} );

  // std::string print;
  // for (auto trj:primaryTrj) trj.Print(print);

  TDatabasePDG pdg_db;
  std::vector<SParticleInfo> particleInfos;
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

    // std::cout << pi.pdg_code << " " << pi.id << " " << pi.mass << " " << pi.charge << std::endl;
    std::cout << "Initial Momentum " << trj.GetInitialMomentum().Vect().Mag() << std::endl;
    // pi.pos.Print();
    // pi.mom.Print();
  }

  int nParticles = particleInfos.size();
  
  if (nParticles == 0) {
    std::cerr << "no particles to be reconstructed...process aborted"
              << std::endl;
    return;
  }

  for (int ip = 0; ip < nParticles; ip++) {
    int charge = particleInfos[ip].charge;
    double particle_mass = particleInfos[ip].mass;
    int pdg = particleInfos[ip].pdg_code;
  
    STTKFTrack this_track;
    STTKFStateVector current_state;
    STTKFKalmanFilterManager manager;
    STTKFKalmanFilterManager::Orientation current_orientation = manager.GetOrientation();
    
    // To Do: find a smarter algorithm to compute this
    TMatrixD initial_cov_matrix(5, 5);
    initial_cov_matrix[0][0] = pow(200E-6, 2);
    initial_cov_matrix[1][1] = pow(200E-6, 2);
    initial_cov_matrix[2][2] = pow(0.1, 2);
    initial_cov_matrix[3][3] = pow(0.01, 2);
    initial_cov_matrix[4][4] = pow(0.01, 2);

   
    // To Do: implment a seeding algorithm
    STTKFStateVector initial_state_vector = STTKFCheck::get_state_vector(particleInfos[ip].mom * 1E-3, // GeV
                                                                         particleInfos[ip].pos * 1E-3, // m 
                                                                         particleInfos[ip].charge);


    // initial_cov_matrix.Print();


    std::vector<STTKFStateCovarianceMatrix> propagator_matrices;

    STTKFTrackStep trackStep;
    trackStep.SetStage(STTKFTrackStep::STTKFTrackStateStage::kPrediction,
                        STTKFState(initial_state_vector, initial_cov_matrix));
    trackStep.SetStage(STTKFTrackStep::STTKFTrackStateStage::kFiltering,
                        STTKFState(initial_state_vector, initial_cov_matrix));

    propagator_matrices.push_back(initial_cov_matrix);

    this_track.AddStep(trackStep);


    double previous_z = particleInfos[ip].pos.Z();
    auto current_z_it = std::prev(z_to_tracklets.lower_bound(particleInfos[ip].pos.Z()));
    for ( ; current_z_it != z_to_tracklets.begin(); current_z_it--) {

      // Get previous kf step
      auto previous_step =
            this_track.GetStep(this_track.GetSteps().size() - 1);
      auto previousStateVector =
          previous_step
              .GetStage(STTKFTrackStep::STTKFTrackStateStage::kFiltering)
              .GetStateVector();
      auto previousCovMatrix =
          previous_step
              .GetStage(STTKFTrackStep::STTKFTrackStateStage::kFiltering)
              .GetStateCovMatrix();

      
      // Compute energy loss for the new step
      auto dir = -1. * manager.GetDirectiveCosinesFromStateVector(
                               previousStateVector);
      
      // To Do: check if this is still valid and add a real fix if needed
      if (dir.Z() > 0) {
        dir *= -1;
      }

      auto current_mom = SANDTrackerUtils::GetMomentumInMeVFromRadiusInMM(
                              previousStateVector.Radius(),
                              previousStateVector.TanLambda()) / 1000;
      double gamma = sqrt(current_mom * current_mom + particle_mass * particle_mass) /
                     particle_mass;
      double beta_from_gamma = sqrt(1 - pow(1 / gamma, 2));

      // std::cout << previous_z << " "
      //           << previousStateVector.X() << " "
      //           << previousStateVector.Y() << " "
      //           << current_z_it->first << " " << std::endl;
      // dir.Print();

      // To Do: check all units
      auto de_step = STTKFGeoManager::GetDE(
                          current_z_it->first, 
                          1000 * previousStateVector.X(), 1000 * previousStateVector.Y(), previous_z, 
                          dir.X(), dir.Y(), dir.Z(),
                          beta_from_gamma, particle_mass, charge) / 1000;
      // std::cout << "de_step " << de_step << std::endl;

      double dz = (current_z_it->first - previous_z) / 1000;
      // std::cout << dz << std::endl;

      // Propagation of state and cov matrix
      auto predictedStateVector = manager.PropagateState(
          previousStateVector, dz, de_step, particle_mass);
      auto nextPhi = predictedStateVector.Phi();
      TMatrixD covariance_noise(5, 5);
      covariance_noise = manager.GetProcessNoiseMatrix(
            previousStateVector, nextPhi, dz, de_step, previous_z,
            particle_mass);
      
      auto propagatorMatrix = manager.GetPropagatorMatrix(
          previousStateVector, nextPhi, dz, de_step, particle_mass);
      
      auto predictedCovMatrix = manager.PropagateCovMatrix(
          previousCovMatrix, propagatorMatrix, covariance_noise);
        

      // Filling prediction step of track
      STTKFTrackStep currentTrackState;
      currentTrackState.SetStage(
          STTKFTrackStep::STTKFTrackStateStage::kPrediction,
          STTKFState(predictedStateVector, predictedCovMatrix));

      // previousStateVector().Print();
      // predictedStateVector().Print();

      // Prediction and filtering on the next plane
      STTKFMeasurement prediction =
          manager.GetPrediction(current_orientation, predictedStateVector);
      auto projectionMatrix =
          manager.GetProjectionMatrix(current_orientation, predictedStateVector);
      auto measurementNoiseMatrix = manager.GetMeasurementNoiseMatrix();

      auto kalmanGainMatrix = manager.GetKalmanGainMatrix(
          predictedCovMatrix, projectionMatrix, measurementNoiseMatrix);

      TMatrixD projectionMatrixTransposed(TMatrixD::kTransposed,
                                          projectionMatrix);
      auto Sk = measurementNoiseMatrix + projectionMatrix *
                                              predictedCovMatrix *
                                              projectionMatrixTransposed;
      
      
      // prediction.Print();
      double best_chi = 1E9;
      STTKFMeasurement best_measurement(2, 1);
      for (const auto& tracklet:current_z_it->second) {
        STTKFMeasurement measurement(2, 1);
        // To Do: vertical and horizontal are outdated and confusing. Replace with something more meaningful.
        // Notice: vertical planes means horizontal measurements and the opposite
        if (current_orientation == STTKFKalmanFilterManager::Orientation::kVertical) {
          measurement[0][0] = tracklet[0] / 1000.;
          measurement[1][0] = M_PI_2 - tracklet[2];
        } else {
          measurement[0][0] = tracklet[1] / 1000.;
          measurement[1][0] = tracklet[3];
        }
        // measurement[0][0] += rand->Gaus(0, sigmaPos);
        // measurement[1][0] += rand->Gaus(0, sigmaAng);
        // measurement.Print();
        // std::cout << tracklet[0] << " " << tracklet[1] << " " << tracklet[2] << " " << tracklet[3] << std::endl;
        auto chi2 = manager.EvalChi2(measurement, prediction, Sk);
        if (chi2 < best_chi) {
          best_chi = chi2;
          best_measurement = measurement;
        }
      }
      if (best_chi > 1.5) {
        continue;
      }
      propagator_matrices.push_back(propagatorMatrix);

      // std::cout << "chi2 " << best_chi << std::endl;
      // best_measurement.Print();

      auto filteredStateVector =
          manager.FilterState(predictedStateVector, kalmanGainMatrix,
                              best_measurement, prediction);
      // filteredStateVector().Print();
      
      auto filteredCovMatrix = manager.FilterCovMatrix(
            predictedCovMatrix, projectionMatrix, measurementNoiseMatrix);

      currentTrackState.SetStage(
          STTKFTrackStep::STTKFTrackStateStage::kFiltering,
          STTKFState(filteredStateVector, filteredCovMatrix));
      this_track.AddStep(currentTrackState);


      if (current_orientation == STTKFKalmanFilterManager::Orientation::kVertical) {
        current_orientation = STTKFKalmanFilterManager::Orientation::kHorizontal;
      } else {
        current_orientation = STTKFKalmanFilterManager::Orientation::kVertical;
      }

      previous_z = current_z_it->first;
    }

    auto reco_state =
          this_track.GetSteps()[this_track.GetSteps().size() - 1]
              .GetStage(STTKFTrackStep::STTKFTrackStateStage::kFiltering)
              .GetStateVector();
      auto reco_mom = SANDTrackerUtils::GetMomentumInMeVFromRadiusInMM(
                                 reco_state.Radius(), reco_state.TanLambda());
    std::cout << "Initial Reco Momentum " << reco_mom << std::endl;
  }

}

int main(int argc, char* argv[])
{
  gStyle->SetOptStat(0);  

  TFile f(argv[2], "READ");
  TGeoManager* geo = 0;
  geo = (TGeoManager*)f.Get("EDepSimGeometry");


  // MC info tree
  TTree* t_h = (TTree*)f.Get("EDepSimEvents");
  TG4Event* ev = new TG4Event;
  t_h->SetBranchAddress("Event", &ev);
  
  TFile f_d(argv[3], "READ");
  TTree* t = (TTree*)f_d.Get("tDigit");

  std::vector<dg_wire>* digits = 0;
  t->SetBranchAddress("dg_wire", &digits);

    
  TFile* h_out = new TFile("h_out.root", "RECREATE");
  TH1D*  h_res = new TH1D("h_res", "h_res", 1000,-100,100);
  TH1D*  h_minima1000 = new TH1D("minima1000", "minima1000", 1000,0,100000);
  TH1D*  h_minima_100 = new TH1D("minima100", "minima100", 1000,0,100);
  TH1D*  h_minima_0_1 = new TH1D("minima0.1", "minima0.1", 1000,0,0.1);
  TH1D*  h_minima_0_0001 = new TH1D("minima0.0001", "minima0.0001", 1000,0,0.0001);

  SANDGeoManager sand_geo;
  sand_geo.init(geo);

  for (int i = 2; i < 3; i++) {
    t_h->GetEntry(i);
    t->GetEntry(i);

    ProcessEventWithKF(&sand_geo, ev, digits);


    // int p[9] = {100, -2000, 2000, 100, -3200, -2300, 100, 23800, 26000};

    // SANDTrackerDigitCollection::FillMap(digits);
    // SANDTrackerClusterCollection clusters(&sand_geo, SANDTrackerDigitCollection::GetDigits(), SANDTrackerClusterCollection::ClusteringMethod::kCellAdjacency);
    // auto digit_map =  SANDTrackerDigitCollection::GetDigits();
    
    // TrackletFinder traklet_finder;
    // traklet_finder.SetVolumeParameters(p);
    // traklet_finder.SetSigmaPosition(0.2);
    // traklet_finder.SetSigmaAngle(0.2);
    

    // TCanvas* canvas_cluster = new TCanvas("canvas_cluster","canvas_cluster",2000,1000);
    // canvas_cluster->Divide(2,1);
    
    // TH2D* h_cluster_yz = new TH2D("h","h", p[6],p[7], p[8], p[3],p[4], p[5]);
    // TH2D* h_cluster_xz = new TH2D("h","h", p[6],p[7], p[8], p[0],p[1], p[2]);
    // canvas_cluster->cd(1);
    // h_cluster_yz->Draw();
    // canvas_cluster->cd(2);
    // h_cluster_xz->Draw();

    // canvas_cluster->Print("clu.pdf(","pdf");

    // std::map<double, std::vector<TVectorD>> z_to_tracklets;

    // int color = 2;
    // for (const auto& container:clusters.GetContainers()) {
    //   int gg = 0;
    //   for (const auto& cluster_in_container:container->GetClusters()) {
    //     std::cout << (double)gg / container->GetClusters().size() * 100 << std::endl;
    //     gg++;
    //     if (gg == 500) break;
    //     if (color > 9) color = 2;
        

    //     traklet_finder.SetCells(cluster_in_container);
    //     auto minima = traklet_finder.FindTracklets();
    //     // Draw tracklets
    //     if (minima.size() != 0) {
    //       canvas_cluster->cd();
    //       std::sort(minima.begin(), minima.end(),
    //                 [](TVectorD v1, TVectorD v2){ return v1[4] < v2[4];});
    //       double z_start = cluster_in_container.GetZ();
    //       for (uint trk = 0; trk < minima.size(); trk++) {
    //         h_minima1000->Fill(minima[trk][4]);
    //         h_minima_100->Fill(minima[trk][4]);
    //         h_minima_0_1->Fill(minima[trk][4]);
    //         h_minima_0_0001->Fill(minima[trk][4]);
            
    //         if (minima[trk][4] < 1E-2) {
    //           // std::cout << minima[trk][0] << " " << minima[trk][2] << std::endl;
              
    //           z_to_tracklets[cluster_in_container.GetZ()].push_back(minima[trk]);

    //           TVector2 start_tracklet_yz(z_start, minima[trk][1]);
    //           TVector2 start_tracklet_xz(z_start, minima[trk][0]);
    //           double z_end = z_start + 5 * cos(minima[trk][3]);
    //           double y_end = minima[trk][1] + 5 * sin(minima[trk][3]);
    //           double x_end = minima[trk][0] + 5 * cos(minima[trk][2]);
    //           TVector2 end_tracklet_yz(z_end, y_end);
    //           TVector2 end_tracklet_xz(z_end, x_end);
              
    //           TLine* line_yz_tracklet = new TLine(start_tracklet_yz.X(), start_tracklet_yz.Y(), end_tracklet_yz.X(), end_tracklet_yz.Y());
    //           TLine* line_xz_tracklet = new TLine(start_tracklet_xz.X(), start_tracklet_xz.Y(), end_tracklet_xz.X(), end_tracklet_xz.Y());
    //           line_yz_tracklet->SetLineColor(color);
    //           line_yz_tracklet->SetLineWidth(1);
    //           line_xz_tracklet->SetLineColor(color);
    //           line_xz_tracklet->SetLineWidth(1);
              
    //           canvas_cluster->cd(1);
    //           line_yz_tracklet->Draw();
    //           canvas_cluster->cd(2);
    //           line_xz_tracklet->Draw();
    //         }
    //       }
    //     }

    //     auto digitId_to_drift_time = traklet_finder.GetDigitToDriftTimeMap();
        
    //     // Draw cells of all digits
    //     for (auto digit:digit_map) {
    //       auto cell = sand_geo.get_cell_info(SANDTrackerCellID(digit.did));
    //       double h,w;
    //       cell->second.size(w,h);
    //       TBox* box_yz = new TBox(cell->second.wire().center().Z() - h/2., cell->second.wire().center().Y() - w/2., cell->second.wire().center().Z() + h/2., cell->second.wire().center().Y() + w/2.);
    //       TBox* box_xz = new TBox(cell->second.wire().center().Z() - h/2., cell->second.wire().center().X() - w/2., cell->second.wire().center().Z() + h/2., cell->second.wire().center().X() + w/2.);
    //       box_yz->SetFillStyle(0);
    //       box_yz->SetLineColor(1);
    //       box_yz->SetLineWidth(1);
    //       box_xz->SetFillStyle(0);
    //       box_xz->SetLineColor(1);
    //       box_xz->SetLineWidth(1);
    //       canvas_cluster->cd(1);
    //       box_yz->Draw();
    //       canvas_cluster->cd(2);
    //       box_xz->Draw();
    //     }

    //     std::vector<SANDTrackerDigitID> digits_cluster = cluster_in_container.GetDigits();
    //     for (uint d = 0; d < digits_cluster.size(); d++) {
    //       canvas_cluster->cd();

    //       auto digit = SANDTrackerDigitCollection::GetDigit(digits_cluster[d]);
    //       auto cell = sand_geo.get_cell_info(SANDTrackerCellID(digit.did));

    //       // Draw lines connecting cells in cluster
    //       // if (d < digits_cluster.size() - 1) {
    //       //   auto next_digit = SANDTrackerDigitCollection::GetDigit(digits_cluster[d+1]);
    //       //   auto next_cell  = sand_geo.get_cell_info(SANDTrackerCellID(next_digit.did));
    //       //   TLine* line_yz1 = new TLine(cell->second.wire().center().Z(), cell->second.wire().center().Y(), next_cell->second.wire().center().Z(), next_cell->second.wire().center().Y());
    //       //   canvas_cluster->cd();
    //       //   line_yz1->SetLineWidth(1);
    //       //   line_yz1->SetLineColor(color);
    //       //   line_yz1->Draw();
    //       // }
          
    //       // Draw cells of cluster
    //       double h,w;
    //       cell->second.size(w,h);
    //       TBox* box_yz = new TBox(cell->second.wire().center().Z() - h/2., cell->second.wire().center().Y() - w/2., cell->second.wire().center().Z() + h/2., cell->second.wire().center().Y() + w/2.);
    //       TBox* box_xz = new TBox(cell->second.wire().center().Z() - h/2., cell->second.wire().center().X() - w/2., cell->second.wire().center().Z() + h/2., cell->second.wire().center().X() + w/2.);
    //       box_yz->SetFillStyle(0);
    //       box_yz->SetLineColor(color);
    //       box_yz->SetLineWidth(1);
    //       box_xz->SetFillStyle(0);
    //       box_xz->SetLineColor(color);
    //       box_xz->SetLineWidth(1);
    //       canvas_cluster->cd(1);
    //       box_yz->Draw();
    //       canvas_cluster->cd(2);
    //       box_xz->Draw();
          

    //       // Draw reco drift time of digits in cluster
    //       TEllipse* el_yz_comp = new TEllipse(cell->second.wire().center().Z(), cell->second.wire().center().Y(), 
    //                             sand_reco::stt::wire_radius + cell->second.driftVelocity() * digitId_to_drift_time[digits_cluster[d]]);
    //       TEllipse* el_xz_comp = new TEllipse(cell->second.wire().center().Z(), cell->second.wire().center().X(), 
    //                             sand_reco::stt::wire_radius + cell->second.driftVelocity() * digitId_to_drift_time[digits_cluster[d]]);
    //       el_yz_comp->SetFillStyle(0);
    //       el_yz_comp->SetLineColor(color);
    //       el_yz_comp->SetLineWidth(1);
    //       el_xz_comp->SetFillStyle(0);
    //       el_xz_comp->SetLineColor(color);
    //       el_xz_comp->SetLineWidth(1);
    //       canvas_cluster->cd(1);
    //       el_yz_comp->Draw();
    //       canvas_cluster->cd(2);
    //       el_xz_comp->Draw();
          
    //       // Draw true drift time of digits in cluster
    //       TEllipse* el_yz = new TEllipse(cell->second.wire().center().Z(), cell->second.wire().center().Y(), 
    //                             sand_reco::stt::wire_radius + cell->second.driftVelocity() * digit.drift_time);
    //       TEllipse* el_xz = new TEllipse(cell->second.wire().center().Z(), cell->second.wire().center().X(), 
    //                             sand_reco::stt::wire_radius + cell->second.driftVelocity() * digit.drift_time);
    //       el_yz->SetFillStyle(0);
    //       el_yz->SetLineWidth(1);
    //       el_yz->SetLineColor(1);
    //       el_xz->SetFillStyle(0);
    //       el_xz->SetLineWidth(1);
    //       el_xz->SetLineColor(1);
    //       canvas_cluster->cd(1);
    //       el_yz->Draw();
    //       canvas_cluster->cd(2);
    //       el_xz->Draw();

    //       // Draw hit segments for the cluster
    //       for (auto& kk:digit.hindex) {
    //         const TG4HitSegment& hseg = ev->SegmentDetectors[digit.det].at(kk);
    //         TLine* l_yz = new TLine(hseg.Start.Z(), hseg.Start.Y(), hseg.Stop.Z(), hseg.Stop.Y());
    //         TLine* l_xz = new TLine(hseg.Start.Z(), hseg.Start.X(), hseg.Stop.Z(), hseg.Stop.X());
    //         l_yz->SetLineColor(1);
    //         l_xz->SetLineColor(1);
    //         canvas_cluster->cd(1);
    //         l_yz->Draw();
    //         canvas_cluster->cd(2);
    //         l_xz->Draw();
    //       }
          
    //       // h_res->Fill(digitId_to_drift_time[d] - digit.drift_time);
    //     }
    //     color++;
    //     canvas_cluster->Write();
    //     canvas_cluster->Print("clu.pdf","pdf");
    //     canvas_cluster->Clear();

    //     canvas_cluster->Divide(2,1);
    //     canvas_cluster->cd(1);
    //     h_cluster_yz->Draw();
    //     canvas_cluster->cd(2);
    //     h_cluster_xz->Draw();
    //     traklet_finder.Clear();

    //   }
    // }
    // canvas_cluster->Print("clu.pdf)","pdf");
    // h_minima1000->Write();
    // h_minima_100->Write();
    // h_minima_0_1->Write();
    // h_minima_0_0001->Write();

    // int sum = 0;
    // for (auto el:z_to_tracklets) {
    //   std::cout << "At z = " << el.first << " there are " << el.second.size() << " tracklets" << std::endl;
    //   sum += el.second.size();
    // }
    // std::cout << "Total tracklets: " << sum << std::endl;




    // TCanvas* canvas_digitization = new TCanvas("canvas_digitization","canvas_digitization",2000,1000);
    // canvas_digitization->Divide(2,1);
    // TH2D* h_digitization_yz = new TH2D("h","h", p[6],p[7], p[8], p[3],p[4], p[5]);
    // TH2D* h_digitization_xz = new TH2D("h","h", p[6],p[7], p[8], p[0],p[1], p[2]);
    // canvas_digitization->cd(1);
    // h_digitization_yz->Draw();
    // canvas_digitization->cd(2);
    // h_digitization_xz->Draw();
    // for (const auto& digit:digit_map) {
    //   auto cell = sand_geo.get_cell_info(SANDTrackerCellID(digit.did));
          
    //   TEllipse* el_yz = new TEllipse(cell->second.wire().center().Z(), cell->second.wire().center().Y(), sand_reco::stt::wire_radius + cell->second.driftVelocity() * digit.drift_time);
    //   TEllipse* el_xz = new TEllipse(cell->second.wire().center().Z(), cell->second.wire().center().X(), sand_reco::stt::wire_radius + cell->second.driftVelocity() * digit.drift_time);
    //   double h,w;
    //   cell->second.size(w,h);
    //   TBox* box_yz = new TBox(cell->second.wire().center().Z() - h/2., cell->second.wire().center().Y() - w/2., cell->second.wire().center().Z() + h/2., cell->second.wire().center().Y() + w/2.);
    //   TBox* box_xz = new TBox(cell->second.wire().center().Z() - h/2., cell->second.wire().center().X() - w/2., cell->second.wire().center().Z() + h/2., cell->second.wire().center().X() + w/2.);

    //   canvas_digitization->cd(1);
    //   el_yz->SetFillStyle(0);
    //   el_yz->Draw();
      
    //   box_yz->SetFillStyle(0);
    //   box_yz->SetLineWidth(1);
    //   canvas_digitization->cd(1);
    //   box_yz->Draw();

    //   canvas_digitization->cd(2);
    //   el_xz->SetFillStyle(0);
    //   el_xz->Draw();
      
    //   box_xz->SetFillStyle(0);
    //   box_xz->SetLineWidth(1);
    //   box_xz->Draw();

    //   for (auto& hi:digit.hindex) {
    //     const TG4HitSegment& hseg = ev->SegmentDetectors[digit.det].at(hi);
    //     TLine* l_yz = new TLine(hseg.Start.Z(), hseg.Start.Y(), hseg.Stop.Z(), hseg.Stop.Y());
    //     canvas_digitization->cd(1);
    //     l_yz->Draw();
    //     TLine* l_xz = new TLine(hseg.Start.Z(), hseg.Start.X(), hseg.Stop.Z(), hseg.Stop.X());
    //     canvas_digitization->cd(2);
    //     l_xz->Draw();
    //   }

    //   TMarker* mark_yz = new TMarker(cell->second.wire().center().Z(), cell->second.wire().center().Y(), 5);
    //   mark_yz->SetMarkerColor(1);
    //   mark_yz->SetMarkerSize(0.5);
    //   canvas_digitization->cd(1);
    //   mark_yz->Draw();

    //   TMarker* mark_xz = new TMarker(cell->second.wire().center().Z(), cell->second.wire().center().X(), 5);
    //   mark_xz->SetMarkerColor(1);
    //   mark_xz->SetMarkerSize(0.5);
    //   canvas_digitization->cd(2);
    //   mark_xz->Draw();
    // }
    // canvas_digitization->SaveAs("./c2D.png");
    // canvas_digitization->SaveAs("./c2D.C");



  }
  h_res->Write();

}