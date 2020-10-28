void extract_info()
{
  gSystem->Load("/data/mt/reco/kloe-simu/lib/libStruct.so");
  
  TFile f("/data/mt/reco/files/reco/numu_internal_10k.0.reco.root");
  TTree* t = (TTree*) f.Get("tDigit");
  TGeoManager* geo = (TGeoManager*) f.Get("EDepSimGeometry");
  
  TFile fmc("/data/mt/reco/files/edep-sim/numu_internal_10k.0.edep-sim.root");
  TTree* tmc = (TTree*) fmc.Get("EDepSimEvents");
  TTree* tgn = (TTree*) fmc.Get("DetSimPassThru/gRooTracker");
  
  std::vector<digit>* digits = new std::vector<digit>();
  t->SetBranchAddress("Stt",&digits);

  TG4Event* ev = new TG4Event();
  tmc->SetBranchAddress("Event",&ev);

  const int max_mcpart = 5000;

  int gn_pdg[max_mcpart];
  double gn_p4[max_mcpart][4];
  
  tgn->SetBranchAddress("StdHepPdg",&gn_pdg);
  tgn->SetBranchAddress("StdHepP4",&gn_p4);
  
  TFile fout("/data/mt/reco/files/reco/ev_info.root","RECREATE");
  TTree tout("tinfo","tinfo");

  double ev_x;
  double ev_y;
  double ev_z;
  double ev_t;
  double ev_px;
  double ev_py;
  double ev_pz;
  double ev_E;
  int ev_nupdg;
  int ev_tgpdg;
  int ev_CC;
  
  const int max_digit = 5000;
  int dg_n;
  double dg_x[max_digit];
  double dg_y[max_digit];
  double dg_z[max_digit];
  double dg_t[max_digit];
  double dg_de[max_digit];
  int dg_did[max_digit];
  int dg_hor[max_digit];
  int dg_tid[max_digit];
  
  const int max_particle = 5000;
  int pp_n;
  double pp_px[max_particle];
  double pp_py[max_particle];
  double pp_pz[max_particle];
  double pp_E[max_particle];
  int pp_pdg[max_particle];
  int pp_tid[max_particle];
  
  tout.Branch("ev_x",&ev_x,"ev_x/D");
  tout.Branch("ev_y",&ev_y,"ev_y/D");
  tout.Branch("ev_z",&ev_z,"ev_z/D");
  tout.Branch("ev_t",&ev_t,"ev_t/D");
  tout.Branch("ev_px",&ev_px,"ev_px/D");
  tout.Branch("ev_py",&ev_py,"ev_py/D");
  tout.Branch("ev_pz",&ev_pz,"ev_pz/D");
  tout.Branch("ev_nupdg",&ev_nupdg,"ev_nupdg/I");
  tout.Branch("ev_tgpdg",&ev_tgpdg,"ev_tgpdg/I");
  tout.Branch("ev_CC",&ev_CC,"ev_CC/I");
  
  tout.Branch("dg_n",&dg_n,"dg_n/I");
  tout.Branch("dg_x",&dg_x,"dg_x[dg_n]/D");
  tout.Branch("dg_y",&dg_y,"dg_y[dg_n]/D");
  tout.Branch("dg_z",&dg_z,"dg_z[dg_n]/D");
  tout.Branch("dg_t",&dg_t,"dg_t[dg_n]/D");
  tout.Branch("dg_de",&dg_de,"dg_de[dg_n]/D");
  tout.Branch("dg_did",&dg_did,"dg_did[dg_n]/I");
  tout.Branch("dg_hor",&dg_hor,"dg_hor[dg_n]/I");
  tout.Branch("dg_tid",&dg_tid,"dg_tid[dg_n]/I");
  
  tout.Branch("pp_n",&pp_n,"pp_n/I");
  tout.Branch("pp_px",&pp_px,"pp_px[pp_n]/D");
  tout.Branch("pp_py",&pp_py,"pp_py[pp_n]/D");
  tout.Branch("pp_pz",&pp_pz,"pp_pz[pp_n]/D");
  tout.Branch("pp_E",&pp_E,"pp_E[pp_n]/D");
  tout.Branch("pp_pdg",&pp_pdg,"pp_pdg[pp_n]/I");
  tout.Branch("pp_tid",&pp_tid,"pp_tid[pp_n]/I");
  
  for(int i = 0; i < t->GetEntries(); i++)
  {
      t->GetEntry(i);
      tmc->GetEntry(i);
      tgn->GetEntry(i);

      ev_x = ev->Primaries.at(0).GetPosition().X();
      ev_y = ev->Primaries.at(0).GetPosition().Y();
      ev_z = ev->Primaries.at(0).GetPosition().Z();
      ev_t = ev->Primaries.at(0).GetPosition().T();

      ev_px = gn_p4[0][0];
      ev_py = gn_p4[0][1];
      ev_pz = gn_p4[0][2];
      ev_E = gn_p4[0][3];

      TString rec(ev->Primaries.at(0).GetReaction());
      std::map<std::string,std::string> fields;

      TObjArray *obj = rec.Tokenize(";");
      for(int j = 0; j < obj->GetEntries(); j++)
      {
        TObjString* str = (TObjString*) obj->At(j);
        TObjArray *obj2 = str->GetString().Tokenize(":");

        std::string key = ((TObjString*) obj2->At(0))->GetString().Data();
        std::string value = ((TObjString*) obj2->At(1))->GetString().Data();

        fields[key] = value;

        delete obj2;
      }

      delete obj;

      ev_nupdg = std::stoi(fields["nu"]);
      ev_tgpdg = std::stoi(fields["tgt"]);
      ev_CC = fields["proc"].find("CC") != std::string::npos ? 1 : 0;
      
      dg_n = digits->size();
      for(int j = 0; j < dg_n; j++)
      {
          dg_x[j] = digits->at(j).x;
          dg_y[j] = digits->at(j).y;
          dg_z[j] = digits->at(j).z;
          dg_t[j] = digits->at(j).t;
          dg_de[j] = digits->at(j).de;
          dg_did[j] = digits->at(j).did;
          dg_hor[j] = digits->at(j).hor;

          int hit_id = digits->at(j).hindex.at(0);
          dg_tid[j] = ev->SegmentDetectors["Straw"][hit_id].GetPrimaryId();
      }

      pp_n = ev->Primaries.at(0).Particles.size();
      for(int j = 0; j < pp_n; j++)
      {
          pp_px[j] = ev->Primaries.at(0).Particles.at(j).GetMomentum().X();
          pp_py[j] = ev->Primaries.at(0).Particles.at(j).GetMomentum().Y();
          pp_pz[j] = ev->Primaries.at(0).Particles.at(j).GetMomentum().Z();
          pp_E[j] = ev->Primaries.at(0).Particles.at(j).GetMomentum().T();
          pp_pdg[j] = ev->Primaries.at(0).Particles.at(j).GetPDGCode();
          pp_pdg[j] = ev->Primaries.at(0).Particles.at(j).GetTrackId();
      }

      tout.Fill();
  }
  
  fout.cd();
  geo->Write("geo");
  tout.Write("tdig");

  f.Close();
  fout.Close();
}