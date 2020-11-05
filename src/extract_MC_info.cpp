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
    TGeoManager* geo = (TGeoManager*) fin.Get("EDepSimGeometry");

    TG4Event* ev = new TG4Event();

    tMC->SetBranchAddress("Event",&ev);

    TFile fout(fname.ReplaceAll(".root",".MCinfo.root").Data(),"RECREATE");
    TTree thit("thit","thit");
    TTree ttub("ttub","ttub");
    TTree tpla("tpla","tpla");
    TTree tmod("tmod","tmod");

    int eventID = 0;

    const int max_hit = 10000;
    int h_n;
    int h_id[max_hit];
    int h_tid[max_hit];
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

    thit.Branch("eventID",&eventID,"eventID/I");
    thit.Branch("h_n",&h_n,"h_n/I");
    thit.Branch("h_id",&h_id,"h_id[h_n]/I");
    thit.Branch("h_tid",&h_tid,"h_tid[h_n]/I");
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

    int nevent = tMC->GetEntries();
    //int nevent = 10;
    
    std::cout << "Events: " << nevent << " [";
    std::cout << std::setw(3) << int(0) << "%]" << std::flush;

    for(int i = 0; i < nevent; i++)
    {
        tMC->GetEntry(i);

        std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i) / nevent * 100)
              << "%]" << std::flush;
        
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
        }
        thit.Fill();
    }

    std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
    std::cout << std::endl;
    
    kloe_simu::init(geo);

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
    fout.Close();

    fin.Close();
}