#include <TGeoManager.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>

#include "../include/struct.h"

#include "/opt/exp_software/neutrino/EDEPSIM/include/EDepSim/TG4Event.h"

void extract_info()
{
  gSystem->Load("/data/mt/reco/kloe-simu/lib/libStruct.so");
  
  TFile f("/data/mt/reco/files/reco/numu_internal_10k.1.reco.root");
  TTree* t = (TTree*) f.Get("tDigit");
  TGeoManager* geo = (TGeoManager*) f.Get("EDepSimGeometry");
  
  TFile fmc("/data/mt/reco/files/edep-sim/numu_internal_10k.1.edep-sim.root");
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
  
  TFile fout("/data/mt/reco/files/reco/event_info.1.root","RECREATE");
  TTree teve("teve","teve");
  TTree tgeo("tgeo","tgeo");

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
  
  double pl_z;
  double pl_id;
  double pl_uid;
  double pl_hor;
  
  teve.Branch("ev_x",&ev_x,"ev_x/D");
  teve.Branch("ev_y",&ev_y,"ev_y/D");
  teve.Branch("ev_z",&ev_z,"ev_z/D");
  teve.Branch("ev_t",&ev_t,"ev_t/D");
  teve.Branch("ev_px",&ev_px,"ev_px/D");
  teve.Branch("ev_py",&ev_py,"ev_py/D");
  teve.Branch("ev_pz",&ev_pz,"ev_pz/D");
  teve.Branch("ev_nupdg",&ev_nupdg,"ev_nupdg/I");
  teve.Branch("ev_tgpdg",&ev_tgpdg,"ev_tgpdg/I");
  teve.Branch("ev_CC",&ev_CC,"ev_CC/I");
  
  teve.Branch("dg_n",&dg_n,"dg_n/I");
  teve.Branch("dg_x",&dg_x,"dg_x[dg_n]/D");
  teve.Branch("dg_y",&dg_y,"dg_y[dg_n]/D");
  teve.Branch("dg_z",&dg_z,"dg_z[dg_n]/D");
  teve.Branch("dg_t",&dg_t,"dg_t[dg_n]/D");
  teve.Branch("dg_de",&dg_de,"dg_de[dg_n]/D");
  teve.Branch("dg_did",&dg_did,"dg_did[dg_n]/I");
  teve.Branch("dg_hor",&dg_hor,"dg_hor[dg_n]/I");
  teve.Branch("dg_tid",&dg_tid,"dg_tid[dg_n]/I");
  
  teve.Branch("pp_n",&pp_n,"pp_n/I");
  teve.Branch("pp_px",&pp_px,"pp_px[pp_n]/D");
  teve.Branch("pp_py",&pp_py,"pp_py[pp_n]/D");
  teve.Branch("pp_pz",&pp_pz,"pp_pz[pp_n]/D");
  teve.Branch("pp_E",&pp_E,"pp_E[pp_n]/D");
  teve.Branch("pp_pdg",&pp_pdg,"pp_pdg[pp_n]/I");
  teve.Branch("pp_tid",&pp_tid,"pp_tid[pp_n]/I");
  
  tgeo.Branch("pl_z",&pl_z,"pl_z/D");
  tgeo.Branch("pl_id",&pl_id,"pl_id/I");
  tgeo.Branch("pl_uid",&pl_uid,"pl_uid/I");
  tgeo.Branch("pl_hor",&pl_hor,"pl_hor/I");
  
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
          pp_tid[j] = ev->Primaries.at(0).Particles.at(j).GetTrackId();
      }

      teve.Fill();
  }
  
  TString stt_path = "volWorld_PV/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/MagIntVol_volume_PV_0/volSTTFULL_PV_0/";
  TString full_path;
  TGeoNode* pla;
  TString pla_name;
  
  double master[3];
  double local[] = {0., 0., 0.};
  
  TGeoVolume* stt = geo->FindVolumeFast("volSTTFULL_PV");
  
  for(int i = 0; i < stt->GetNdaughters(); i++)
  {
    TGeoNode* mod = stt->GetNode(i);
    TString mod_name = mod->GetName();
    
    if(mod_name.Contains("volfrontST") == true)
    {
      TObjArray* oba = mod_name.Tokenize("_");
      
      TString* str = (TString*) oba->At(0);
      int mid = str->ReplaceAll("volfrontST","").Atoi();
    
      for(int j = 0; j < mod->GetVolume()->GetNdaughters(); j++)
      {
         pla = mod->GetVolume()->GetNode(j);
         pla_name = pla->GetName();
         
         if(pla_name.Contains("hor") == true)
         {
           full_path = stt_path + mod_name + "/" + pla_name;
           
           geo->cd(full_path.Data());
           
           geo->LocalToMaster(local,master);
           
           pl_z = master[2];
           pl_id = mid + 90;
           pl_uid = pl_id * 10 + 2;
           pl_hor = 1;
           
           tgeo.Fill();
         }
         else if(pla_name.Contains("ver") == true)
         {
           full_path = stt_path + mod_name + "/" + pla_name;
           
           geo->cd(full_path.Data());
           
           geo->LocalToMaster(local,master);
           
           pl_z = master[2];
           pl_id = mid;
           pl_uid = mid * 10 + 1;
           pl_hor = 0;
           
           tgeo.Fill();
         }
      }
      
      delete oba;     
    }
    else if(mod_name.Contains("sttmod") == true)
    {
      TObjArray* oba = mod_name.Tokenize("_");
      
      TString* str = (TString*) oba->At(0);
      int mid = str->ReplaceAll("volfrontST","").Atoi();
    
      for(int j = 0; j < mod->GetVolume()->GetNdaughters(); j++)
      {
         pla = mod->GetVolume()->GetNode(j);
         pla_name = pla->GetName();
         
         if(pla_name.Contains("hor") == true)
         {
           full_path = stt_path + mod_name + "/" + pla_name;
           
           geo->cd(full_path.Data());
           
           geo->LocalToMaster(local,master);
           
           pl_z = master[2];
           pl_id = mid;
           pl_uid = mid * 10 + 2;
           pl_hor = 1;
           
           tgeo.Fill();
         }
         else if(pla_name.Contains("ver") == true)
         {
           full_path = stt_path + mod_name + "/" + pla_name;
           
           geo->cd(full_path.Data());
           
           geo->LocalToMaster(local,master);
           
           pl_z = master[2];
           pl_id = mid;
           pl_uid = pl_id * 10 + 1;
           pl_hor = 0;
           
           tgeo.Fill();
         }
      } 
      delete oba;    
    }
  }
  
  fout.cd();
  tgeo.Write("tgeo");
  teve.Write("teve");

  f.Close();
  fout.Close();
}