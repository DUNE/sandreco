#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TGeoManager.h>

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#include "/opt/exp_software/neutrino/EDEPSIM/include/EDepSim/TG4Event.h"
#include "/opt/exp_software/neutrino/EDEPSIM/include/EDepSim/TG4HitSegment.h"

#include "../include/utils.h"

void extract_MC_info(const char* filename)
{
    TString fname(filename);
    TFile fin(fname.Data());
    TTree* tMC = (TTree*) fin.Get("EDepSimEvents");
    TTree* tgn = (TTree*) fin.Get("DetSimPassThru/gRooTracker");
    TGeoManager* geo = (TGeoManager*) fin.Get("EDepSimGeometry");
    
    kloe_simu::init(geo);

    TG4Event* ev = new TG4Event();

    tMC->SetBranchAddress("Event",&ev);
    
    const int max_mcpart = 5000;

    int gn_pdg[max_mcpart];
    double gn_p4[max_mcpart][4];
  
    tgn->SetBranchAddress("StdHepPdg",&gn_pdg);
    tgn->SetBranchAddress("StdHepP4",&gn_p4);

    TFile fout(fname.ReplaceAll(".root",".MCinfo.root").Data(),"RECREATE");
    TTree thit("thit","thit");
    TTree ttub("ttub","ttub");
    TTree tpla("tpla","tpla");
    TTree tmod("tmod","tmod");
    TTree tpar("tpar","tpar");
    TTree teve("teve","teve");

    int eventID = 0;

    const int max_hit = 10000;
    int h_n;
    int h_id[max_hit];
    int h_tid[max_hit];
    int h_did[max_hit];
    double h_de[max_hit];
    double h_x1[max_hit];
    double h_y1[max_hit];
    double h_z1[max_hit];
    double h_t1[max_hit];
    double h_x2[max_hit];
    double h_y2[max_hit];
    double h_z2[max_hit];
    double h_t2[max_hit];

    int tb_id;
    int tb_hor;
    double tb_x;
    double tb_y;
    double tb_z;

    int pl_id;
    int pl_hor;
    double pl_z;

    int md_id;
    double md_z;
    
    const int max_particle = 5000;
    int pp_n;
    double pp_px[max_particle];
    double pp_py[max_particle];
    double pp_pz[max_particle];
    double pp_E[max_particle];
    int pp_pdg[max_particle];
    int pp_tid[max_particle];
    
    double ev_x;
    double ev_y;
    double ev_z;
    double ev_t;
    double ev_px;
    double ev_py;
    double ev_pz;
    double ev_e;
    int ev_nupdg;
    int ev_tgpdg;
    int ev_CC;

    thit.Branch("eventID",&eventID,"eventID/I");
    thit.Branch("h_n",&h_n,"h_n/I");
    thit.Branch("h_id",&h_id,"h_id[h_n]/I");
    thit.Branch("h_tid",&h_tid,"h_tid[h_n]/I");
    thit.Branch("h_did",&h_did,"h_did[h_n]/I");
    thit.Branch("h_de",&h_de,"h_de[h_n]/D");
    thit.Branch("h_x1",&h_x1,"h_x1[h_n]/D");
    thit.Branch("h_y1",&h_y1,"h_y1[h_n]/D");
    thit.Branch("h_z1",&h_z1,"h_z1[h_n]/D");
    thit.Branch("h_t1",&h_t1,"h_t1[h_n]/D");
    thit.Branch("h_x2",&h_x2,"h_x2[h_n]/D");
    thit.Branch("h_y2",&h_y2,"h_y2[h_n]/D");
    thit.Branch("h_z2",&h_z2,"h_z2[h_n]/D");
    thit.Branch("h_t2",&h_t2,"h_t2[h_n]/D");
    
    ttub.Branch("tb_id",&tb_id,"tb_id/I");
    ttub.Branch("tb_hor",&tb_hor,"tb_hor/I");
    ttub.Branch("tb_x",&tb_x,"tb_x/D");
    ttub.Branch("tb_y",&tb_y,"tb_y/D");
    ttub.Branch("tb_z",&tb_z,"tb_z/D");
    
    tpla.Branch("pl_id",&pl_id,"pl_id/I");
    tpla.Branch("pl_hor",&pl_hor,"pl_hor/I");
    tpla.Branch("pl_z",&pl_z,"pl_z/D");
    
    tmod.Branch("md_id",&md_id,"md_id/I");
    tmod.Branch("md_z",&md_z,"md_z/D");

    tpar.Branch("pp_n",&pp_n,"pp_n/I");
    tpar.Branch("pp_px",&pp_px,"pp_px[pp_n]/D");
    tpar.Branch("pp_py",&pp_py,"pp_py[pp_n]/D");
    tpar.Branch("pp_pz",&pp_pz,"pp_pz[pp_n]/D");
    tpar.Branch("pp_E",&pp_E,"pp_E[pp_n]/D");
    tpar.Branch("pp_pdg",&pp_pdg,"pp_pdg[pp_n]/I");
    tpar.Branch("pp_tid",&pp_tid,"pp_tid[pp_n]/I");
    
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

    int nevent = tMC->GetEntries();
    //int nevent = 10;
    
    std::cout << "Events: " << nevent << " [";
    std::cout << std::setw(3) << int(0) << "%]" << std::flush;

    for(int i = 0; i < nevent; i++)
    {
        tMC->GetEntry(i);
        tgn->GetEntry(i);

        std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i) / nevent * 100)
              << "%]" << std::flush;
        

        ev_x = ev->Primaries.at(0).GetPosition().X();
        ev_y = ev->Primaries.at(0).GetPosition().Y();
        ev_z = ev->Primaries.at(0).GetPosition().Z();
        ev_t = ev->Primaries.at(0).GetPosition().T();

        ev_px = gn_p4[0][0];
        ev_py = gn_p4[0][1];
        ev_pz = gn_p4[0][2];
        ev_e = gn_p4[0][3];

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

        teve.Fill();
        
        h_n = 0;
        eventID = i;

        h_n = ev->SegmentDetectors["Straw"].size();

        for(int j = 0; j < h_n; j++)
        {
            h_id[j] = j;
            h_tid[j] = ev->SegmentDetectors["Straw"].at(j).GetPrimaryId();
            h_de[j] = ev->SegmentDetectors["Straw"].at(j).GetEnergyDeposit();
            h_x1[j] = ev->SegmentDetectors["Straw"].at(j).GetStart().X();
            h_y1[j] = ev->SegmentDetectors["Straw"].at(j).GetStart().Y();
            h_z1[j] = ev->SegmentDetectors["Straw"].at(j).GetStart().Z();
            h_t1[j] = ev->SegmentDetectors["Straw"].at(j).GetStart().T();
            h_x2[j] = ev->SegmentDetectors["Straw"].at(j).GetStop().X();
            h_y2[j] = ev->SegmentDetectors["Straw"].at(j).GetStop().Y();
            h_z2[j] = ev->SegmentDetectors["Straw"].at(j).GetStop().Z();
            h_t2[j] = ev->SegmentDetectors["Straw"].at(j).GetStop().T();
            h_did[j] = kloe_simu::getSTUniqID(geo, 0.5*(h_x1[j]+h_x2[j]), 0.5*(h_y1[j]+h_y2[j]), 0.5*(h_z1[j]+h_z2[j]));
        }
        thit.Fill();

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
        tpar.Fill();
    }

    std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
    std::cout << std::endl;

    for(std::map<int, std::map<int, TVector2> >::iterator it = kloe_simu::stPos.begin(); it != kloe_simu::stPos.end(); it++)
    {
        std::map<int, TVector2>& plane = it->second;
        pl_id = it->first;
        pl_hor = pl_id % 2 == 0 ? 1 : 0;
        pl_z = 0.5*(plane.at(0).Y() + plane.at(1).Y());
        tpla.Fill();

        if(pl_id % 10 == 2)
        {
            md_id = int(pl_id/10);
            md_z = 0.5*(kloe_simu::stPos.at(pl_id).at(1).X() + kloe_simu::stPos.at(pl_id-1).at(0).X());
            tmod.Fill();
        }

        for(std::map<int, TVector2>::iterator ite = plane.begin(); ite != plane.end(); ite++)
        {
            tb_id = 1000 * ite->first + pl_id;
            tb_hor = pl_hor;
            tb_x = tb_hor == 1 ? 0. : ite->second.Y();
            tb_y = tb_hor == 0 ? 0. : ite->second.Y();
            tb_z = ite->second.X();
            ttub.Fill();
        }
    }

    fout.cd();
    thit.Write();
    ttub.Write();
    tpla.Write();
    tmod.Write();
    tpar.Write();
    teve.Write();
    fout.Close();

    fin.Close();
}