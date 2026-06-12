bool test_dE_comparison(std::string inFile, int charge, double particle_mass, double pdg_code, double dz, TFile& out_file)
{
  STTKFKalmanFilterManager manager;
  std::string name = "DdE_" + std::to_string(pdg_code);
  TH1D* h = new TH1D(name.c_str(), name.c_str(), 100, -5, 5);
  name = "BBdE_" + std::to_string(pdg_code);
  TH1D* hh = new TH1D(name.c_str(), name.c_str(), 1000, 0, 10);
  name = "TdE_" + std::to_string(pdg_code);
  TH1D* hh2 = new TH1D(name.c_str(), name.c_str(), 1000, 0, 10);
  name = "TdE_vs_BBdE_" + std::to_string(pdg_code);
  TH2D* h2 = new TH2D(name.c_str(), name.c_str(), 1000, 0, 10, 1000, 0, 10);
  name = "comp_BBdE_" + std::to_string(pdg_code);
  TH2D* h3 = new TH2D(name.c_str(), name.c_str(), 1000, 0, 10, 1000, 0, 10);
  name = "edep_BBdE_" + std::to_string(pdg_code);
  TH2D* h4 = new TH2D(name.c_str(), name.c_str(), 1000, 0, 10, 1000, 0, 10);
  name = "res_vs_mom_" + std::to_string(pdg_code);
  TH2D* h5 = new TH2D(name.c_str(), name.c_str(), 1000, 0, 10, 100, -5, 5);
  CheckPoints check_points;
  std::string inFileFull = "/storage/gpfs_data/neutrino/users/vpia/edepsim-prod/numu_default/production/" + inFile;
  check_points.Init(inFileFull, charge, particle_mass, pdg_code);
  out_file.cd();

  for (int ev = 0; ev < 250; ev++) {
    check_points.InitTree(ev);
    std::vector<EDEPTrajectory> primaryTrj;
    check_points.tree.Filter(std::back_insert_iterator<std::vector<EDEPTrajectory>>(primaryTrj), 
    [](const EDEPTrajectory& trj) { return trj.GetParentId() == -1;} );
    
    TDatabasePDG pdg_db;
    // std::vector<SParticleInfo> particleInfos;
    for (auto trj:primaryTrj) {
      // SParticleInfo pi;
      // pi.pdg_code = trj.GetPDGCode();
      // pi.id       = trj.GetId();
      if (trj.GetPDGCode() == pdg_code) {
        check_points.charge = charge;
        check_points.pdg_code = pdg_code;
        check_points.particle_mass = particle_mass;
        check_points.particle_id = trj.GetId();
        break;
      }
    }


    if (!check_points.InitEvent(ev, dz, 200E-6, 0.2)) {
      continue;
    } 
    int nPoints = check_points.GetNPoints();
    std::vector<point> points = check_points.GetPoints();


    for (int i = 1; i < nPoints; i++) {
      STTKFCheck::ParticleState previous_particle_state({points[i - 1].x, points[i - 1].y, points[i - 1].z},
                                                        {points[i - 1].px / 1000, points[i - 1].py / 1000, points[i - 1].pz / 1000,});
      STTKFCheck::ParticleState current_particle_state({points[i].x, points[i].y, points[i].z},
                                                      {points[i].px / 1000, points[i].py / 1000, points[i].pz / 1000,});
      double previous_mom = STTUtils::GetMomentumInMeVFromRadiusInMM(previous_particle_state.get_state_vector(charge).Radius(), 
                                                                    previous_particle_state.get_state_vector(charge).TanLambda()) / 1000;
      double current_mom = STTUtils::GetMomentumInMeVFromRadiusInMM(current_particle_state.get_state_vector(charge).Radius(), 
                                                                    current_particle_state.get_state_vector(charge).TanLambda()) / 1000;
      double current_energy = sqrt(current_mom*current_mom + particle_mass*particle_mass);
      double previous_energy = sqrt(previous_mom*previous_mom + particle_mass*particle_mass);
      
      double gamma = sqrt(previous_mom*previous_mom + particle_mass*particle_mass) / particle_mass;
      double beta = sqrt( 1 - pow(1/gamma, 2));

      auto dir = -1. * manager.GetDirectiveCosinesFromStateVector(previous_particle_state.get_state_vector(charge));
      if (dir.Z() > 0) dir *= -1;
      auto BBdE = STTKFGeoManager::GetDE(1000*points[i].z, 
                                        1000*points[i - 1].x, 
                                        1000*points[i - 1].y, 
                                        1000*points[i - 1].z,
                                        dir.X(), 
                                        dir.Y(), 
                                        dir.Z(),
                                        beta,
                                        particle_mass,
                                        charge);

      double true_dE = (current_energy - previous_energy)*1000;

      double crossedMaterial = STTKFGeoManager::GetCrossedMaterialInGCM2(1000*points[i].z, 
                                                                1000*points[i - 1].x, 
                                                                1000*points[i - 1].y, 
                                                                1000*points[i - 1].z,
                                                                dir.X(), 
                                                                dir.Y(), 
                                                                dir.Z());
      h->Fill((BBdE-true_dE)/true_dE);
      hh->Fill(BBdE);
      hh2->Fill(true_dE);
      h2->Fill(true_dE, BBdE);
      h3->Fill(previous_mom, BBdE/crossedMaterial);
      h4->Fill(previous_mom, true_dE/crossedMaterial);
      h5->Fill(previous_mom, (BBdE-true_dE)/true_dE);
    }
  }
  h->Write();
  hh->Write();
  hh2->Write();
  h2->Write();
  h3->Write();
  h4->Write();
  h5->Write();


}

void tryCompleteManager(sand_reco::kf::utils::TrackletMap z_to_tracklets, SParticleInfo particle, TH1D* h_gpos_distribution, TH1D* h_gang_distribution, TH1D* x_res, TH1D* y_res, TH1D* theta_y_res, TH1D* theta_x_res, TH1D* mom_res, TMultiGraph* mg, TMultiGraph* mgx) {
  sand_reco::kf::Manager manager;
  manager.initFromMC(&z_to_tracklets, particle);
  manager.run();

  auto track = manager.getTrack();
  if (track.getSteps().size() > 0) { // was commented, with > 3
    auto last_step = track.getSteps().back(); //crash if empty due to the .back().
    auto reco_state =
          last_step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kSmoothing).getStateVector();
    auto reco_mom = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(
                                  reco_state.radius(), reco_state.tanLambda());

    auto initial_state = sand_reco::kf::utils::getStateVector(particle.initial_mom * 1E-3, particle.initial_pos * 1E-3, particle.charge);
    auto initial_mom = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(initial_state.radius(), initial_state.tanLambda());
 

    std::cout << "Initial Momentum " << initial_mom << std::endl;
    std::cout << "Initial Smoothed Reco Momentum " << reco_mom << std::endl;
    std::cout << "track.getSteps() " << track.getSteps().size() << std::endl;
    
    TGraph* yz_predicted = new TGraph(track.getSteps().size());
    TGraph* yz_filtered = new TGraph(track.getSteps().size());
    TGraph* yz_smoothed = new TGraph(track.getSteps().size());
    TGraph* yz_measured = new TGraph(track.getSteps().size());
    TGraph* xz_predicted = new TGraph(track.getSteps().size());
    TGraph* xz_filtered = new TGraph(track.getSteps().size());
    TGraph* xz_smoothed = new TGraph(track.getSteps().size());
    TGraph* xz_measured = new TGraph(track.getSteps().size());

    x_res->Fill(initial_state.x() - reco_state.x()); 
    y_res->Fill(initial_state.y() - reco_state.y()); 
    theta_y_res->Fill(initial_state.phi() - reco_state.phi()); 
    theta_x_res->Fill(initial_state.tanLambda() - reco_state.tanLambda()); 
    mom_res->Fill(initial_mom - reco_mom); 


    
    int i = 0;
    for (auto& step : track.getSteps()) {
      auto prediction = step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kPrediction).getStateVector();
      auto filtering = step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kFiltering).getStateVector();
      auto smoothing =  step.getStage(sand_reco::kf::TrackStep::TrackStateStage::kSmoothing).getStateVector();
      
      yz_predicted->SetPoint(i, step.getZ(), prediction.y()*1000 );
      yz_filtered->SetPoint(i, step.getZ() , filtering.y()*1000);
      yz_smoothed->SetPoint(i, step.getZ() , smoothing.y()*1000);
      yz_measured->SetPoint(i, step.getZ() , step.getY());
      xz_predicted->SetPoint(i, step.getZ(), prediction.x()*1000 );
      xz_filtered->SetPoint(i, step.getZ() , filtering.x()*1000);
      xz_smoothed->SetPoint(i, step.getZ() , smoothing.x()*1000);
      xz_measured->SetPoint(i, step.getZ() , step.getX());
      i++;
      
      auto& innovation = step.getInnovation();
      if (innovation.empty()) {
        continue;
      }
  
      h_gpos_distribution->Fill(innovation[0]);
      h_gang_distribution->Fill(innovation[1]);
    
    }
      
    yz_predicted->SetLineColor(3);
    yz_predicted->SetMarkerStyle(3);
    mg->Add(yz_predicted);
    yz_filtered->SetLineColor(4);
    yz_filtered->SetMarkerStyle(4);
    mg->Add(yz_filtered);
    yz_smoothed->SetLineColor(6);
    yz_smoothed->SetMarkerStyle(5);
    mg->Add(yz_smoothed);
    yz_measured->SetLineColor(2);
    yz_measured->SetMarkerStyle(2);
    mg->Add(yz_measured);

    xz_predicted->SetLineColor(3);
    xz_predicted->SetMarkerStyle(3);
    mgx->Add(xz_predicted);
    xz_filtered->SetLineColor(4);
    xz_filtered->SetMarkerStyle(4);
    mgx->Add(xz_filtered);
    xz_smoothed->SetLineColor(6);
    xz_smoothed->SetMarkerStyle(5);
    mgx->Add(xz_smoothed);
    xz_measured->SetLineColor(2);
    xz_measured->SetMarkerStyle(2);
    mgx->Add(xz_measured);
  }

  return;
}

void processEventWithKF(SANDGeoManager* sand_geo, TG4Event* mc_event, std::vector<dg_wire>* digits, TH1D* h_gpos_distribution,TH1D* h_gang_distribution,
                        TH1D* h_x_diff, TH1D* h_y_diff, TH1D* h_theta_x_diff, TH1D* h_theta_y_diff, TH1D* x_res, TH1D* y_res, TH1D* theta_y_res, TH1D* theta_x_res, TH1D* mom_res)
{
  
  int p[9] = {100, -2000, 2000, 100, -4000, -0, 100, 22500, 26000};

  sand_reco::tracker::DigitCollection::fillMap(digits);
  auto digit_map =  sand_reco::tracker::DigitCollection::getDigits();
  if (sand_reco::tracker::DigitCollection::getDigits().empty()) {
    return;
  }
  std::string tracker_name = sand_reco::tracker::DigitCollection::getDigits().begin()->det;
  sand_reco::tracker::ClusterCollection clusters(sand_geo, sand_reco::tracker::DigitCollection::getDigits(), sand_reco::tracker::ClusterCollection::ClusteringMethod::kCellAdjacency);
  
  std::map<double, std::vector<Tracklet>> z_to_tracklets;

  SANDTrackerUtils::init(sand_geo->getTGeoManager());

  TRandom3 rand(0);
  for (const auto& container:clusters.getContainers()) {
    for (const auto& cluster_in_container:container->getClusters()) {

      TVector3 first_point;
      TVector3 last_point;
      double min_z = 10e8;
      double max_z = -10e8;
      for (uint d = 0; d < cluster_in_container.getDigits().size(); d++) {
        auto digit = sand_reco::tracker::DigitCollection::getDigit(cluster_in_container.getDigits()[d]);
        
        if (digit.z > max_z) {
          max_z = digit.z;
          last_point = TVector3(digit.x, digit.y, digit.z);
        }
        if (digit.z < min_z) {
          min_z = digit.z;
          first_point = TVector3(digit.x, digit.y, digit.z);
        }
      }

      auto true_tracklet = getTrueTrackletOfCluster(first_point, last_point, cluster_in_container.getZ());
      
      TVector3 true_pos = true_tracklet[0];
      TVector3 true_dir = true_tracklet[1];     
      double true_theta_yz = atan(true_dir.Y() / true_dir.Z());
      double true_theta_xz = atan(true_dir.X() / true_dir.Z());
      if (true_theta_xz > M_PI_2) true_theta_xz -= M_PI;

      Tracklet measurement_from_true_tracklet;
      measurement_from_true_tracklet.x = true_pos.X()  + rand.Gaus(0, SANDTrackerUtils::getSigmaPositionMeasurement() * 1E3);
      measurement_from_true_tracklet.y = true_pos.Y()  + rand.Gaus(0, SANDTrackerUtils::getSigmaPositionMeasurement() * 1E3);
      measurement_from_true_tracklet.theta_xz = true_theta_xz + rand.Gaus(0, SANDTrackerUtils::getSigmaAngleMeasurement());
      measurement_from_true_tracklet.theta_yz = true_theta_yz + rand.Gaus(0, SANDTrackerUtils::getSigmaAngleMeasurement());
      for (uint d = 0; d < cluster_in_container.getDigits().size(); d++) {
        auto digit = sand_reco::tracker::DigitCollection::getDigit(cluster_in_container.getDigits()[d]);
        measurement_from_true_tracklet.digits.push_back(digit);
      }
      // std::cout << true_pos.X() << std::endl;

      z_to_tracklets[cluster_in_container.getZ()].push_back(measurement_from_true_tracklet);
    }
  }
  
  if (z_to_tracklets.empty()) {
    return;
  }

  EDEPTree tree;
  tree.InizializeFromEdep(*mc_event, sand_geo->getTGeoManager());
  
  std::vector<EDEPTrajectory> primaryTrj;
  tree.Filter(std::back_insert_iterator<std::vector<EDEPTrajectory>>(primaryTrj), 
    [](const EDEPTrajectory& trj) { return trj.GetParentId() == -1;} );

  TDatabasePDG pdg_db;
  std::vector<SParticleInfo> particleInfos;
  std::map<double, std::vector<TVectorD>> z_to_best_tracklet;

  double sigma_pos = 0;
  double sigma_mom = 0;
  std::vector<int> indeces;
  int ii = -1;
  for (auto trj:primaryTrj) {
    ii++;

    if (trj.GetHitMap().find(string_to_component[tracker_name]) == trj.GetHitMap().end()) {
      continue;
    }
    
    if (trj.GetTrajectoryPoints().find(string_to_component[tracker_name]) == trj.GetTrajectoryPoints().end()) {
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
    pi.id       = trj.GetId();
    pi.mass     = particle->Mass();
    pi.charge   = particle->Charge() / 3;

    double max_z = 0;
    bool to_be_reconstructed = false;
    for (auto& point : trj.GetTrajectoryPoints().at(string_to_component[tracker_name])) {
      if (point.GetPosition().Z() > max_z && point.GetMomentum().Z() > 100) {
        max_z = point.GetPosition().Z();
        pi.pos = point.GetPosition().Vect();
        pi.mom = point.GetMomentum();
        to_be_reconstructed = true;
      }
    }

    if (!to_be_reconstructed) continue;
    double x_smeared = rand.Gaus(pi.pos.X(), sigma_pos);
    double y_smeared = rand.Gaus(pi.pos.Y(), sigma_pos);
    double px_smeared = pi.mom.X() * rand.Gaus(1, sigma_mom);
    double py_smeared = pi.mom.Y() * rand.Gaus(1, sigma_mom);
    double pz_smeared = pi.mom.Z() * rand.Gaus(1, sigma_mom);

    pi.pos = TVector3(x_smeared, y_smeared, pi.pos.Z());
    pi.mom = TVector3(px_smeared, py_smeared, pz_smeared);
    pi.initial_pos = trj.GetTrajectoryPoints().at(string_to_component[tracker_name])[0].GetPosition().Vect();
    pi.initial_mom = trj.GetTrajectoryPoints().at(string_to_component[tracker_name])[0].GetMomentum();
    particleInfos.push_back(pi);
    indeces.push_back(ii);
  }

   
  
  int nParticles = particleInfos.size();
  
  if (nParticles == 0) {
    std::cerr << "no particles to be reconstructed...process aborted"
              << std::endl;
    return;
  }

  for (int ip = 0; ip < (int)indeces.size(); ip++) {
    std::string name_mg = "YZ_" + std::to_string(indeces[ip]);
    TMultiGraph* mg = new TMultiGraph(name_mg.c_str(), name_mg.c_str());
    TGraph* yz_true = new TGraph(primaryTrj[indeces[ip]].GetTrajectoryPoints().at(string_to_component[tracker_name]).size());
    std::string name_mgx = "XZ_" + std::to_string(indeces[ip]);
    TMultiGraph* mgx = new TMultiGraph(name_mgx.c_str(), name_mgx.c_str());
    TGraph* xz_true = new TGraph(primaryTrj[indeces[ip]].GetTrajectoryPoints().at(string_to_component[tracker_name]).size());

    for (uint i = 0; i <  primaryTrj[indeces[ip]].GetTrajectoryPoints().at(string_to_component[tracker_name]).size(); i++){
      auto point = primaryTrj[indeces[ip]].GetTrajectoryPoints().at(string_to_component[tracker_name])[i];
       yz_true->SetPoint(i, point.GetPosition().Z() , point.GetPosition().Y());
       xz_true->SetPoint(i, point.GetPosition().Z() , point.GetPosition().X());
    }


    tryCompleteManager(z_to_tracklets, particleInfos[ip], h_gpos_distribution, h_gang_distribution, x_res, y_res, theta_y_res, theta_x_res, mom_res, mg, mgx);

    std::string title = name_mg + "; z [mm]; y [mm]";
    mg->SetTitle(title.c_str());
    yz_true->SetMarkerStyle(4);
    mg->Add(yz_true);
    mg->Write();
    title = name_mgx + "; z [mm]; x [mm]";
    mgx->SetTitle(title.c_str());
    xz_true->SetMarkerStyle(4);
    mgx->Add(xz_true);
    mgx->Write();
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

  SANDGeoManager sand_geo;
  sand_geo.init(geo);
  
  std::string geometry;
  if (geo->FindVolumeFast("STTtracker_PV")) {
    geometry = "STT";
  } else if (geo->FindVolumeFast("SANDtracker_PV")) {
    geometry = "DRIFT";
  } 
  sand_geo.fillAdjacentCells(geometry);

  TFile* h_out = new TFile("h_out.root", "RECREATE");

  int nev = t->GetEntries();
  for (int i = 0; i < nev; i++) {
    t_h->GetEntry(i);
    t->GetEntry(i);

    if (!plots) {
      innovation_test->cd();
      processEventWithKF(&sand_geo, ev, digits);
    }
  }
}