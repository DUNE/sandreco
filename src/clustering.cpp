#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>


#include "../include/Linkdef.h"
#include "../include/struct.h"
#include "../include/utils.h"
#include "TSystem.h"

using namespace std;

std::ofstream infofile;
std::ofstream info_new; 
std::tuple<double, double, double, double> fit_ls(int, double[], double[],
                                                  double[]);

std::pair<std::vector<dg_cell>, std::vector<dg_cell>> ProcessMultiHits(
    std::vector<dg_cell>, std::vector<dg_cell>);
std::pair<std::vector<dg_cell>, std::vector<int>> GetNeighbours(
    std::vector<dg_cell>, int, std::vector<int>, std::vector<dg_cell>);

std::vector<cluster> Clusterize(std::vector<dg_cell>*);
std::vector<cluster> TrackFit(std::vector<cluster>);
std::vector<cluster> Merge(std::vector<cluster>);
std::vector<cluster> Split(std::vector<cluster>);
std::vector<cluster> RecoverIncomplete(std::vector<cluster>,
                                       std::vector<dg_cell>);
cluster Calc_variables(std::vector<dg_cell>);

bool RepetitionCheck(std::vector<int>, int);
bool isNeighbour(int, int);

double TfromTDC(double t1, double t2, double L);
double AttenuationFactor(double d, int planeID);
double EfromADC(double adc1, double adc2, double d1, double d2, int planeID);
double EfromADCsingle(double adc, double f);
double DfromTDC(double, double);

void Clust_info(cluster);

int clustering(std::string const& input)
{
    // TH1F* ClusterEnergy = new TH1F("ClusterEnergy", "Total Energy reconstructed
    // in the ECal;Energy [MeV];Entries/bins", 100, 0., 1500.);
    //gROOT->ProcessLine("gInterpreter->AddIncludePath("-I/opt/exp_software/neutrino/EDEPSIM/include/EDepSim")");
    gSystem->Load("libStruct.so");
    const char* finname = input.c_str();
    TFile f(finname, "READ");
    TTree* t = (TTree*)f.Get("tDigit");
    int nEvents = t->GetEntries();
    std::vector<dg_cell>* cell = new std::vector<dg_cell>;
    std::vector<cluster> f_clust, og_clust;
    TString output = input;
    output.ReplaceAll(".root", ".Clusters.root");
    TFile fout(output, "RECREATE");
    TTree tout("tCluster", "Clustering");
    tout.Branch("cluster", "std::vector<cluster>", &f_clust);
    t->SetBranchAddress("dg_cell", &cell);
    //infofile.open("info.txt");
    info_new.open("muon_6marzo.txt");
    //nEvents = 10; //DC
    for (int i = 0; i < nEvents; i++) {
        //infofile<< "-----------" << endl;
        infofile<< "Entry " << i << endl;
        t->GetEntry(i);
        std::vector<cluster> clust = Clusterize(std::move(cell));
        double CluEn = 0;
        for (auto const& clu_info : clust) {
            Clust_info(clu_info);
            
            if(clu_info.sy >0.2){
                info_new<< "Entry with cluster.sy> 0.2 : " << i << std::endl;  
            }
        }
        f_clust = clust;
        tout.Fill();
        clust.clear();
        f_clust.clear();
        og_clust.clear();
    }
    fout.cd();
    tout.Write();
    fout.Close();
    delete cell;
    return 0;
}

std::vector<cluster> Clusterize(std::vector<dg_cell>* vec_cellraw)
{
    std::vector<dg_cell> complete_cells, broken_cells, multicomplete_cells;
    std::vector<cluster> vec_clust;
    // Create vector of complete (signal on both the photosensors) and incomplete cells
    for (auto const& cell : *vec_cellraw) {
        if (cell.ps1.size() == 0 && cell.ps2.size() == 0) {
            ////infofile<< "In Clusterize: Found empty cell with id: " << cell.id << endl;
            continue;
        } else if (cell.ps1.size() == 0 || cell.ps2.size() == 0) {
            ////infofile<< "In Clusterize: Found borken cell with id: " << cell.id << " and ps size:" << cell.ps1.size() << " "<<  cell.ps2.size() << endl;
            broken_cells.push_back(cell);
        } else if((cell.ps1.size() != 0 && cell.ps2.size() != 0)) {
            // Complete cell
            ////infofile<< "In Clusterize: Found complete cell with id: " << cell.id << " and ps size:" << cell.ps1.size() << " "<<  cell.ps2.size() << endl;
            complete_cells.push_back(cell);
        }
    }
    // Funzione che prende dentro complete_cells e cerca se ci sono celle complete
    // con multiple hits. L'input deve essere vettore celle in cui la size di
    // ps1/ps2 deve essere >= 1, l'output deve essere un vettore di celle in cui
    // la size di ps1/ps2 � sempre 1. In aggiunta a questo ci deve essere un
    // vettore di broken cells in cui si aggiungoino i multiple hits che non hanno
    // un corrispettivo sull'altro ps.
    
    
    //Celle multicomplete possono esserci per pile up, se una sola particella passa nel 
    //detector un solo ps1/2 viene generato
    std::pair<std::vector<dg_cell>, std::vector<dg_cell>> processed_cells =
    ProcessMultiHits(complete_cells, broken_cells);
    multicomplete_cells = processed_cells.first;
    broken_cells = processed_cells.second;
    
    
    std::vector<int> checked_array;
    std::vector<int> vec_cell;
    
    std::vector<int> chck;
    ////infofile<<"Fill vec_clust, Get Neighbours:"<<std::endl;
    for (int i = 0; i < multicomplete_cells.size(); i++) {
        ////infofile<<"Position in multicomplete_cells:"<< i <<std::endl;
        std::vector<dg_cell> v_cell;
        // if (i == 0) {
            //   chck.push_back(i);
            // } else if (RepetitionCheck(chck, i) == true)
            // continue;
            if (RepetitionCheck(chck, i) == true){
                continue;
            } else {
                chck.push_back(i);
            }
            ////infofile<<"inserted id cell in vec_cell (seed): "<< multicomplete_cells.at(i).id << std::endl;
            v_cell.push_back(multicomplete_cells.at(i));
            //DC add chck.push_back(i);?? chck viene riempito in GetNeighb con "i+1"
            std::pair<std::vector<dg_cell>, std::vector<int>> Neighbours =
            GetNeighbours(multicomplete_cells, i, chck, v_cell);
            v_cell = Neighbours.first;
            chck = Neighbours.second;
            struct cluster Clust;
            ////infofile<< "cells given to Calc_variables(v_cell)" <<std::endl;
            for(int k=0; k<v_cell.size();k++){
                //infofile<< v_cell.at(k).id << " with PS size: " << v_cell.at(k).ps1.size() << ", " << v_cell.at(k).ps2.size() <<std::endl;
            }
            Clust = Calc_variables(v_cell);
            vec_clust.push_back(Clust);
            //here 
            //Clust_info(Clust);
            
        }
        // SPLIT
        //infofile<<"split: "<<std::endl;
        vec_clust = Split(vec_clust);
        //for (auto const& clu_info : vec_clust) {
          //  Clust_info(clu_info);
        //}
        //infofile<<"merge: "<<std::endl;
        // MERGE
        vec_clust = Merge(vec_clust);
        //for (auto const& clu_info : vec_clust) {
          //  Clust_info(clu_info);
        //}
        //infofile<<"track_fit: "<<std::endl;
        // Track Fit
        vec_clust = TrackFit(vec_clust);
        //for (auto const& clu_info : vec_clust) {
            // Clust_info(clu_info);
            //}
            //infofile<<"recover incomplete"<<std::endl;
            
            // Recover Incomplete cells
            vec_clust = RecoverIncomplete(vec_clust, broken_cells);
            //for (auto const& clu_info : vec_clust) {
                //  Clust_info(clu_info);
                //}
                return vec_clust;
            }
            
            // per tutte le celle complete cerca se esistono e quali sono i  ps1.tdc - ps2.tdc < delta (distanza temporale massima possibile tra t_a e t_b)
            //in modo da creare complete cells anche per eventi di  multiple hits. 
            
            std::pair<std::vector<dg_cell>, std::vector<dg_cell>> ProcessMultiHits(
                std::vector<dg_cell> og_cell, std::vector<dg_cell> incomplete_cells)
                {
                    
        ////infofile<< "IN PROCESS MULTI HITS"<< endl;
        std::vector<dg_cell> complete_cells;
        for (auto const& cell : og_cell) {
            double delta = cell.l * kloe_simu::vlfb / kloe_simu::m_to_mm; 
            ////infofile<< delta << " -> "<<cell.l<< endl;
            for (int i = 0; i < cell.ps1.size(); i++) {
                int found = 0;
                for (int j = 0; j < cell.ps2.size(); j++) {
                    if (abs(cell.ps1.at(i).tdc - cell.ps2.at(j).tdc) < delta) {
                        dg_cell good_cell;
                        good_cell.id = cell.id;
                        good_cell.z = cell.z;
                        good_cell.x = cell.x;
                        good_cell.y = cell.y;
                        good_cell.l = cell.l;
                        good_cell.mod = cell.mod;
                        good_cell.lay = cell.lay;
                        good_cell.cel = cell.cel;
                        good_cell.ps1.push_back(cell.ps1.at(i));
                        good_cell.ps2.push_back(cell.ps2.at(j));
                        complete_cells.push_back(good_cell);
                        //infofile<< " cell id " << good_cell.id << " [ " << good_cell.ps1.at(0).tdc << " " << good_cell.ps2.at(0).tdc << " ] " << endl;
                        found++;
                        break;
                    }
                }
                if (found == 0) {
                    // buttare in broken cells
                    dg_cell ps1bad_cell;
                    ps1bad_cell.id = cell.id;
                    ps1bad_cell.z = cell.z;
                    ps1bad_cell.x = cell.x;
                    ps1bad_cell.y = cell.y;
                    ps1bad_cell.l = cell.l;
                    ps1bad_cell.mod = cell.mod;
                    ps1bad_cell.lay = cell.lay;
                    ps1bad_cell.cel = cell.cel;
                    ps1bad_cell.ps1.push_back(cell.ps1.at(i));
                    incomplete_cells.push_back(ps1bad_cell);
                    //infofile<<" ps1 - Brok cell id "<<ps1bad_cell.id<<" ("  <<ps1bad_cell.ps1.at(0).tdc <<" ) "<<endl;
                }
            }
            
            
            //perch� farlo due volte?? 
            for (int k = 0; k < cell.ps2.size(); k++) {
                int found = 0;
                for (int l = 0; l < cell.ps1.size(); l++) {
                    if (abs(cell.ps1.at(l).tdc - cell.ps2.at(k).tdc) < delta) {
                        found++;
                    }
                }
                if (found == 0) {
                    // buttare in broken cells
                    dg_cell ps2bad_cell;
                    ps2bad_cell.id = cell.id;
                    ps2bad_cell.z = cell.z;
                    ps2bad_cell.x = cell.x;
                    ps2bad_cell.y = cell.y;
                    ps2bad_cell.l = cell.l;
                    ps2bad_cell.mod = cell.mod;
                    ps2bad_cell.lay = cell.lay;
                    ps2bad_cell.cel = cell.cel;
                    ps2bad_cell.ps2.push_back(cell.ps2.at(k));
                    incomplete_cells.push_back(ps2bad_cell);
                    //infofile<< " ps2 - Brok cell id " << ps2bad_cell.id << " ( " << ps2bad_cell.ps2.at(0).tdc << " ) " << endl;
                }
            }
        }
        return std::make_pair(complete_cells, incomplete_cells);
    } 
    
    std::vector<cluster> RecoverIncomplete(std::vector<cluster> clus,
    std::vector<dg_cell> incomplete_cells)
    {
        
        for (auto const& brok_cells : incomplete_cells) {
            //infofile<< "Broken cells found: " << brok_cells.id << " and ps size: " << brok_cells.ps1.size() << " " << brok_cells.ps2.size() << endl;
            if(brok_cells.ps1.size() > 1){
                for(auto ps1: brok_cells.ps1){
                    //infofile<< "tdc for ps1: " << ps1.tdc << endl;
                }
            } else if(brok_cells.ps2.size() > 1) {
                for(auto ps2: brok_cells.ps2){
                    //infofile<< "tdc for ps2: " << ps2.tdc << endl;
                }
            }
        }
        vector<double> tAvec;
        vector<double> tBvec;
        for(auto cl: clus){
            int splitted = 0;
            double tA = 0, tB = 0, tA2 = 0, tB2 = 0;
            double EA, EB, EAtot = 0, EBtot = 0, EA2tot = 0, EB2tot = 0;
            double tRMS_A, tRMS_B, dist;
            vector<dg_cell>all_cells = cl.cells;
            for (int q = 0; q < all_cells.size(); q++) {
                EA = all_cells.at(q).ps1.at(0).adc;
                EB = all_cells.at(q).ps2.at(0).adc;
                EAtot += EA;
                EA2tot += EA * EA;
                EBtot += EB;
                EB2tot += EB * EB;
                tA += (all_cells.at(q).ps1.at(0).tdc)*
                EA;
                tA2 += std::pow(all_cells.at(q).ps1.at(0).tdc,2) *
                EA;
                tB += (all_cells.at(q).ps2.at(0).tdc) *
                EB;
                tB2 += std::pow(all_cells.at(q).ps2.at(0).tdc,2) *
                EB;
            }
            tA = tA / EAtot;
            tAvec.push_back(tA);
            tA2 = tA2 / EAtot;
            tB = tB / EBtot;
            tBvec.push_back(tB);
            tB2 = tB2 / EBtot;
        }
        
        
        for (auto const& brok_cells : incomplete_cells) {
            int isbarrel = 0;
            if (brok_cells.mod == 30) isbarrel = 1;
            if (brok_cells.mod == 40) isbarrel = 2;
            double cell_phi =
            atan((brok_cells.z - 23910.00) / (brok_cells.y + 2384.73)) * 180 /
            TMath::Pi();
            double cell_theta =
            atan((brok_cells.z - 23910.00) / (brok_cells.x)) * 180 / TMath::Pi();
            int minentry = 0;
            int found = 0;
            //infofile<< "Itering with broken cell " << brok_cells.id << " with tdc: ";
            if(brok_cells.ps1.size() > 0){
                for(auto ps1: brok_cells.ps1){
                    //infofile<< "ps1 ->" << ps1.tdc << endl;
                }
            } else {
                for(auto ps2: brok_cells.ps2){
                    //infofile<< "ps2 ->" << ps2.tdc << endl;
                }
            }
            vector<int> clus_index;
            for (int j = 0; j < clus.size(); j++) {
                double rec_en = 0;
                bool hasNeigh = false;
                int nclusterC = 0;
                //infofile<< "Cluster "<< j <<" tdc ps1 mean " << tAvec.at(j) << " and tdc ps2 mean " << tBvec.at(j) << endl; 
                for(int i = 0; i < clus.at(j).cells.size(); i++){
                    
                    if(isNeighbour(brok_cells.id, clus.at(j).cells.at(i).id)){
                        
                        hasNeigh = true;
                        nclusterC++;
                        
                        if(find(clus_index.begin(), clus_index.end(), j) == clus_index.end()) {
                            clus_index.push_back(j);
                        }
                        
                        
                        //infofile<< "Broken Cell with id: " << brok_cells.id << " is neighbour of cell with id: " <<  clus.at(j).cells.at(i).id << " with PS size: " << clus.at(j).cells.at(i).ps1.size() << " ," <<clus.at(j).cells.at(i).ps2.size() <<" in cluster " << j << endl;
                        
                        for(auto ps1: brok_cells.ps1){
                            for(auto cps1:clus.at(j).cells.at(i).ps1){
                                //infofile<< "ps1 size and time of broken cell: " << brok_cells.ps1.size() << ", " <<  ps1.tdc << " time of complete neighbour cell: " << cps1.tdc;
                                for(auto cps2: clus.at(j).cells.at(i).ps2){
                                    //infofile<< ", for ps2 " << cps2.tdc << endl;
                                }
                            }
                        }
                        for(auto ps2: brok_cells.ps2){
                            for(auto cps1:clus.at(j).cells.at(i).ps1){
                                //infofile<< "ps2 size and time of broken cell: " << brok_cells.ps2.size() << ", " <<  ps2.tdc << " time of complete neighbour cell: " << cps1.tdc;
                                for(auto cps2: clus.at(j).cells.at(i).ps2){
                                    //infofile<< ", for ps2 " << cps2.tdc << endl;
                                }
                            }
                        }
                        
                    }
                }
                
                if(hasNeigh){
                    
                    minentry = j;
                    found = 1;
                    //infofile<< endl;
                    //infofile<< "RECOVER INCOMPLETE -> BROKEN CELL WITH ID: " << brok_cells.id << " associated with cluster: ";
                    //for(auto index: clus_index){
                        //infofile<< index << endl;
                        //}
                        //if(clus_index.size() > 1){
                            //infofile<< "HELP! MORE CLUSTER CANDIDATES FOR THE BROKE CELL!" << endl;
                            //}
                            
                        } else if(found == 0) {
                            double clus_phi =
                            atan((clus.at(j).z - 23910.00) / (clus.at(j).y + 2384.73)) * 180 /
                            TMath::Pi();
                            // if (isbarrel != 0) cout << "Cluster Phi: " << clus_phi << endl;
                            double clus_theta =
                            atan((clus.at(j).z - 23910.00) / (clus.at(j).x)) * 180 / TMath::Pi();
                            // if (isbarrel == 0) cout << "Cluster Theta: " << clus_theta << endl;
                            double minphi = 999, mintheta = 999;
                            int isbarrelc = 0;
                            if (clus.at(j).cells[0].mod == 30) isbarrelc = 1;
                            if (clus.at(j).cells[0].mod == 40) isbarrelc = 2;
                            if (abs(cell_phi - clus_phi) < 3 && (isbarrelc == isbarrel) &&
                            isbarrelc == 0) {
                                found = 1;
                                double dist =
                                sqrt((brok_cells.z - clus.at(j).z) * (brok_cells.z - clus.at(j).z) +
                                (brok_cells.y - clus.at(j).y) * (brok_cells.y - clus.at(j).y));
                                if (abs(cell_phi - clus_phi) < minphi && dist < 2000) {
                                    // cout << "Barrel found 1: Clus Theta/Phi " << clus_theta << " / " <<
                                    // clus_phi << " VS Cell Theta/Phi: " << cell_theta << " " << cell_phi
                                    // << endl; cout << " Cluster: " << clus.at(j).x << " " << clus.at(j).y
                                    // << " " << clus.at(j).z << " VS Cell: " << brok_cells.x << " " <<
                                    // brok_cells.y << " " << brok_cells.z << "  ->  Clus-cell: " <<
                                    // clus.at(j).cells.at(0).id << " and cell id: " << brok_cells.id << "
                                    // -> " << dist << endl;
                                    minphi = abs(cell_phi - clus_phi);
                                    minentry = j;
                                }
                                continue;
                            }
                            if (isbarrel == isbarrelc && isbarrel != 0) {
                                if (abs(cell_theta - clus_theta) < 3) {
                                    // cout << "EndCap found 1: Clus Theta/Phi " << clus_theta << " / " <<
                                    // clus_phi << " VS Cell Theta/Phi: " << cell_theta << " " << cell_phi
                                    // << endl; cout << " Cluster: " << clus.at(j).x << " " << clus.at(j).y
                                    // << " " << clus.at(j).z << " VS Cell: " << brok_cells.x << " " <<
                                    // brok_cells.y << " " << brok_cells.z << "  ->  Clus-cell: " <<
                                    // clus.at(j).cells.at(0).id << "-" << clus.at(j).cells.at(0).mod << "
                                    // and cell id: " << brok_cells.id << endl;
                                    found = 1;
                        mintheta = abs(cell_theta - clus_theta);
                        minentry = j;
                    }
                    continue;
                }
            }
        }
        if (found == 1 && isbarrel == 0) {
            double rec_en = 0;
            double DpmA = clus.at(minentry).x / 10 + 215;
            double DpmB = -clus.at(minentry).x / 10 + 215;
            if (brok_cells.ps1.size() != 0 && brok_cells.ps2.size() != 0) {
                double Ea = brok_cells.ps1.at(0).adc;
                double Eb = brok_cells.ps2.at(0).adc;
                rec_en = kloe_simu::EfromADC(Ea, Eb, DpmA, DpmB, brok_cells.lay);
                clus.at(minentry).e = clus.at(minentry).e + rec_en;
                clus.at(minentry).cells.push_back(brok_cells);
            } else if (brok_cells.ps1.size() != 0) {
                int laycell = brok_cells.lay;
                double f = kloe_simu::AttenuationFactor(DpmA, laycell);
                rec_en = EfromADCsingle(brok_cells.ps1.at(0).adc, f);
                clus.at(minentry).e = clus.at(minentry).e + rec_en;
                clus.at(minentry).cells.push_back(brok_cells);
            } else if (brok_cells.ps2.size() != 0) {
                int laycell = brok_cells.lay;
                double f = kloe_simu::AttenuationFactor(DpmB, laycell);
                rec_en = EfromADCsingle(brok_cells.ps2.at(0).adc, f);
                clus.at(minentry).e = clus.at(minentry).e + rec_en;
                clus.at(minentry).cells.push_back(brok_cells);
            }
        }
        if (found == 1 && isbarrel != 0) {
            double rec_en = 0;
            double ecl = brok_cells.l;
            double DpmA = clus.at(minentry).z / 10 + ecl / 20;
            double DpmB = -clus.at(minentry).z / 10 + ecl / 20;
            if (brok_cells.ps1.size() != 0 && brok_cells.ps2.size() != 0) {
                double Ea = brok_cells.ps1.at(0).adc;
                double Eb = brok_cells.ps2.at(0).adc;
                rec_en = kloe_simu::EfromADC(Ea, Eb, DpmA, DpmB, brok_cells.lay);
                clus.at(minentry).e = clus.at(minentry).e + rec_en;
                clus.at(minentry).cells.push_back(brok_cells);
            } else if (brok_cells.ps1.size() != 0) {
                int laycell = brok_cells.lay;
                double f = kloe_simu::AttenuationFactor(DpmA, laycell);
                rec_en = EfromADCsingle(brok_cells.ps1.at(0).adc, f);
                clus.at(minentry).e = clus.at(minentry).e + rec_en;
                clus.at(minentry).cells.push_back(brok_cells);
            } else if (brok_cells.ps2.size() != 0) {
                int laycell = brok_cells.lay;
                double f = kloe_simu::AttenuationFactor(DpmB, laycell);
                rec_en = EfromADCsingle(brok_cells.ps2.at(0).adc, f);
                clus.at(minentry).e = clus.at(minentry).e + rec_en;
                clus.at(minentry).cells.push_back(brok_cells);
            }
        }
    }
    
    
    return clus;
}

/*void Clust_info(cluster clus)
{
    cout << "=O=O=O=O=O=O=O=O=O=O=O=O=O" << endl;
    cout << "cluster: Energy " << clus.e << " MeV" << endl;
  cout << "Coordinate centroide: " << clus.x << " [X] " << clus.y << " [Y] "
       << clus.z << " [z] e tempo di arrivo medio " << clus.t << " ns" << endl;
  cout << "Varianza: " << clus.varx << " [X] " << clus.vary << " [Y] "
       << clus.varz << " [z]" << endl;
  cout << "Composto dalle seguenti celle: ";
  for (int i = 0; i < clus.cells.size(); i++) {
    cout << "Cell: " << clus.cells.at(i).id << " X: " << clus.cells.at(i).x
         << " Y: " << clus.cells.at(i).y << " Z: " << clus.cells.at(i).z
         << endl;
  }
  cout << endl;
}
*/

void Clust_info(cluster clus){
    infofile<< "Energy " << clus.e << " MeV" << endl; 
    infofile<< "Coordinate centroide: " << clus.x << " [X] " << clus.y << " [Y] "
    << clus.z << " [z] e tempo di arrivo medio " << clus.t << " ns" << endl;
    infofile<< "Varianza: " << clus.varx << " [X] " << clus.vary << " [Y] "
    << clus.varz << " [z]" << endl;
    infofile<< "Composto dalle seguenti celle: ";
    //for (int i = 0; i < clus.cells.size(); i++) {
      //  infofile<< "Cell: " << clus.cells.at(i).id << " X: " << clus.cells.at(i).x
        //<< " Y: " << clus.cells.at(i).y << " Z: " << clus.cells.at(i).z
        //<< endl;
   //}
    infofile<< endl;    
}

std::vector<cluster> Split(std::vector<cluster> original_clu_vec)
{
  std::vector<cluster> clu_vec;
  std::vector<dg_cell> all_cells;
  for (auto const& clus : original_clu_vec) {
    int splitted = 0;
    double tA = 0, tB = 0, tA2 = 0, tB2 = 0;
    double EA, EB, EAtot = 0, EBtot = 0, EA2tot = 0, EB2tot = 0;
    double tRMS_A, tRMS_B, dist;
    all_cells = clus.cells;
    for (int j = 0; j < all_cells.size(); j++) {
      EA = all_cells.at(j).ps1.at(0).adc;
      EB = all_cells.at(j).ps2.at(0).adc;
      EAtot += EA;
      EA2tot += EA * EA;
      EBtot += EB;
      EB2tot += EB * EB;
      double d =
          DfromTDC(all_cells[j].ps1.at(0).tdc, all_cells[j].ps2.at(0).tdc);
      double d1, d2, d3;
      d1 = 0.5 * all_cells[j].l + d;
      d2 = 0.5 * all_cells[j].l - d;
      double cell_E = kloe_simu::EfromADC(all_cells[j].ps1.at(0).adc,
                                          all_cells[j].ps2.at(0).adc, d1, d2,
                                          all_cells[j].lay);
      tA += (all_cells.at(j).ps1.at(0).tdc -
             kloe_simu::vlfb * d1 / kloe_simu::m_to_mm) *
            EA;
      tA2 += std::pow(all_cells.at(j).ps1.at(0).tdc -
                          kloe_simu::vlfb * d1 / kloe_simu::m_to_mm,
                      2) *
             EA;
      tB += (all_cells.at(j).ps2.at(0).tdc -
             kloe_simu::vlfb * d2 / kloe_simu::m_to_mm) *
            EB;
      tB2 += std::pow(all_cells.at(j).ps2.at(0).tdc -
                          kloe_simu::vlfb * d2 / kloe_simu::m_to_mm,
                      2) *
             EB;
    }
    tA = tA / EAtot;
    //infofile<< "tA: "<< tA << std::endl;
    tA2 = tA2 / EAtot;
    tB = tB / EBtot;
    //infofile<< "tB: "<< tB << std::endl;
    tB2 = tB2 / EBtot;
    tRMS_A = sqrt( abs((tA2 - tA * tA) * (EA2tot - EAtot * EAtot) / EA2tot) );
    //infofile<< "tRMS_A: "<< tRMS_A << std::endl;
    tRMS_B = sqrt( abs((tB2 - tB * tB) * (EB2tot - EBtot * EBtot) / EB2tot));
    //infofile<< "tRMS_B: "<< tRMS_B << std::endl;
    dist = std::sqrt(tRMS_A * tRMS_A + tRMS_B * tRMS_B);
    //infofile<< "dist: "<< dist << std::endl; 
    if (dist > 5) {
      std::vector<dg_cell> q1_cells, q2_cells, q3_cells, q4_cells;
      for (auto const& a_cells : all_cells) {
        double t_difA = a_cells.ps1.at(0).tdc - tA;
        double t_difB = a_cells.ps2.at(0).tdc - tB;
        if (t_difA > 0) {
          if (t_difB > 0) {
            q1_cells.push_back(a_cells);
          } else {
            q2_cells.push_back(a_cells);
          }
        } else {
          if (t_difB > 0) {
            q3_cells.push_back(a_cells);
          } else {
            q4_cells.push_back(a_cells);
          }
        }
      }
      if (q1_cells.size() != 0) {
        cluster clus = Calc_variables(q1_cells);
        clu_vec.push_back(clus);
        splitted++;
      }
      if (q2_cells.size() != 0) {
        cluster clus = Calc_variables(q2_cells);
        clu_vec.push_back(clus);
        splitted++;
      }
      if (q3_cells.size() != 0) {
        cluster clus = Calc_variables(q3_cells);
        clu_vec.push_back(clus);
        splitted++;
      }
      if (q4_cells.size() != 0) {
        cluster clus = Calc_variables(q4_cells);
        clu_vec.push_back(clus);
        splitted++;
      }
      q1_cells.clear();
      q2_cells.clear();
      q3_cells.clear();
      q4_cells.clear();
    } else {
      clu_vec.push_back(clus);
    }
    all_cells.clear();
  }
  original_clu_vec.clear();
  return clu_vec;
}

std::vector<cluster> Merge(std::vector<cluster> Og_cluster)
{
  std::vector<cluster> mgd_cluster;
  std::vector<int> checked;
  for (int i = 0; i < Og_cluster.size(); i++) {
    double xi = Og_cluster.at(i).x;
    double yi = Og_cluster.at(i).y;
    double zi = Og_cluster.at(i).z;
    double varxi = Og_cluster.at(i).varx;
    double varyi = Og_cluster.at(i).vary;
    double varzi = Og_cluster.at(i).varz;
    double ti = Og_cluster.at(i).t;
    double ei = Og_cluster.at(i).e;
    bool RepCheck = RepetitionCheck(checked, i);
    if (RepCheck == true) {
      continue;
    }
    checked.push_back(i);
    cluster clust;
    clust.x = Og_cluster.at(i).x;
    clust.y = Og_cluster.at(i).y;
    clust.z = Og_cluster.at(i).z;
    clust.t = Og_cluster.at(i).t;
    clust.e = Og_cluster.at(i).e;
    clust.varx = Og_cluster.at(i).varx;
    clust.vary = Og_cluster.at(i).vary;
    clust.varz = Og_cluster.at(i).varz;
    clust.cells = Og_cluster.at(i).cells;
    for (int j = i; j < Og_cluster.size(); j++) {
      RepCheck = RepetitionCheck(checked, j);
      if (RepCheck == true) {
        continue;
      }
      double xj = Og_cluster.at(j).x;
      double yj = Og_cluster.at(j).y;
      double zj = Og_cluster.at(j).z;
      double tj = Og_cluster.at(j).t;
      double ej = Og_cluster.at(j).e;
      double varxj=Og_cluster.at(j).varx;
      double varyj=Og_cluster.at(j).vary;
      double varzj=Og_cluster.at(j).varz;

      double D = sqrt((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj) +
                      (zi - zj) * (zi - zj));
      double DT = abs(ti - tj);
      double dx = sqrt(std::pow(xi - xj, 2));
      double dy = sqrt(std::pow(yi-yj, 2));
      double dz = sqrt(std::pow(zi -zj, 2));
      // //infofile<< dx << " " << dy << " " << dz << endl;
      // //infofile<< "Variance " << 2*(varxi + varxj) << " " <<  2*(varyi + varyj) << " " << 2*(varzi + varzj) << endl;

      ////infofile<< "D " << D << " DT " << DT << endl; 

      if (D < 40 && DT < 2.5) {
      // if (dx <= 3*(varxi + varxj) && dy <= 3*(varyi + varyj) && dz <= 3*(varzi + varzj) && DT < 2.5) {
        //infofile<< "IT'S MERGIN' TIME!" << endl;
        bool endcap = false;
        if (clust.cells[0].id > 25000) {
          endcap = true;
        }
        if (endcap == true) {
          double Dz_ec = abs(yi - yj);
          D = sqrt((xi - xj) * (xi - xj) + (zi - zj) * (zi - zj));
          if (Dz_ec < 30 && D < 40) {
            std::vector<dg_cell> vec_cells_j = Og_cluster.at(j).cells;
            for (int k = 0; k < vec_cells_j.size(); k++) {
              clust.cells.push_back(vec_cells_j.at(k));
            }
            clust = Calc_variables(clust.cells);
            checked.push_back(j);
          }
        } else if (endcap == false) {
          double Dz_bar = abs(xi - xj);
          D = sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj));
          if (Dz_bar < 30 && D < 40) {
            std::vector<dg_cell> vec_cells_j = Og_cluster.at(j).cells;
            for (int k = 0; k < vec_cells_j.size(); k++) {
              clust.cells.push_back(vec_cells_j.at(k));
            }
            clust = Calc_variables(clust.cells);
            checked.push_back(j);
          }
        }
      }
    }
    mgd_cluster.push_back(clust);
  }
  return mgd_cluster;
}

std::vector<cluster> TrackFit(std::vector<cluster> clu_vec)
{
  const double xl[5] = {4.44, 4.44, 4.44, 4.44, 5.24};
  for (int i = 0; i < clu_vec.size(); i++) {
    double apx[3], eapx[3] = {0, 0, 0}, ctrk[3] = {0, 0, 0},
                   ectrk[3] = {0, 0, 0};
    std::vector<dg_cell> cell_vec_0, cell_vec_1, cell_vec_2, cell_vec_3,
        cell_vec_4;
    for (int j = 0; j < clu_vec.at(i).cells.size(); j++) {
      if (clu_vec.at(i).cells.at(j).lay == 0) {
        cell_vec_0.push_back(clu_vec.at(i).cells.at(j));
      } else if (clu_vec.at(i).cells.at(j).lay == 1) {
        cell_vec_1.push_back(clu_vec.at(i).cells.at(j));
      } else if (clu_vec.at(i).cells.at(j).lay == 2) {
        cell_vec_2.push_back(clu_vec.at(i).cells.at(j));
      } else if (clu_vec.at(i).cells.at(j).lay == 3) {
        cell_vec_3.push_back(clu_vec.at(i).cells.at(j));
      } else if (clu_vec.at(i).cells.at(j).lay == 4) {
        cell_vec_4.push_back(clu_vec.at(i).cells.at(j));
      }
    }
    cluster Lay0, Lay1, Lay2, Lay3, Lay4;
    Lay0 = Calc_variables(cell_vec_0);
    Lay1 = Calc_variables(cell_vec_1);
    Lay2 = Calc_variables(cell_vec_2);
    Lay3 = Calc_variables(cell_vec_3);
    Lay4 = Calc_variables(cell_vec_4);
    double LayE[5] = {Lay0.e, Lay1.e, Lay2.e, Lay3.e, Lay4.e};
    infofile << "Layer 0 : "<< std::endl;
    Clust_info(Lay0); 
    
    infofile << "Layer 1 : "<< std::endl;
    Clust_info(Lay1); 

    infofile << "Layer 2 : "<< std::endl;
    Clust_info(Lay2);

    infofile << "Layer 3 : "<< std::endl;
    Clust_info(Lay3); 

    infofile << "Layer 4 : "<< std::endl;
    Clust_info(Lay4); 

    bool isBarrel = true;
    double yx[5] = {0, 0, 0, 0, 0}, yy[5] = {0, 0, 0, 0, 0},
           yz[5] = {0, 0, 0, 0, 0}, wx[5] = {0, 0, 0, 0, 0},
           wy[5] = {0, 0, 0, 0, 0}, wz[5] = {0, 0, 0, 0, 0};
    double X[5] = {0, 0, 0, 0, 0}, D = 0;
    if (clu_vec.at(i).cells[0].id > 25000) {
      isBarrel = false;
    }
    int lay_cross = 0, first_lay = 0;
    if (Lay0.e > 0) {
      lay_cross++;
      if (lay_cross == 1) {
        first_lay = 1;
      }
      yx[lay_cross - 1] = Lay0.x;
      yy[lay_cross - 1] = Lay0.y;
      yz[lay_cross - 1] = Lay0.z;
      wz[lay_cross - 1] = 0.6; 
      wy[lay_cross - 1] = 0.6;
      wx[lay_cross - 1] = 0.001 * Lay0.e; 
      if (isBarrel == false) {
        wy[lay_cross - 1] = wx[lay_cross - 1];
        wx[lay_cross - 1] = 0.6;
      }
      // D = D + xl[0];
      // if (first_lay == 1) X[0] = 0.5 * xl[0];
    }
    if (Lay1.e > 0) {
      lay_cross++;
      yx[lay_cross - 1] = Lay1.x;
      yy[lay_cross - 1] = Lay1.y;
      yz[lay_cross - 1] = Lay1.z;
      wx[lay_cross - 1] = 0.001 * Lay1.e;
      wy[lay_cross - 1] = 0.6;
      wz[lay_cross - 1] = 0.6;
      if (lay_cross == 1) {
        first_lay = 2;
      }
      if (isBarrel == false) {
        wy[lay_cross - 1] = wx[lay_cross - 1];
        wx[lay_cross - 1] = 0.6;
      }
      // D = D + xl[1];
      // if (first_lay == 2)
      //  X[0] = 0.5 * xl[1];
      // else if (first_lay == 1)
      //  X[1] = xl[0] + 0.5 * xl[1];
    }
    if (Lay2.e > 0) {
      lay_cross++;
      yx[lay_cross - 1] = Lay2.x;
      yy[lay_cross - 1] = Lay2.y;
      yz[lay_cross - 1] = Lay2.z;
      wx[lay_cross - 1] = 0.001 * Lay2.e;
      wy[lay_cross - 1] = 0.6;
      wz[lay_cross - 1] = 0.6;
      if (lay_cross == 1) {
        first_lay = 3;
      }
      if (isBarrel == false) {
        wy[lay_cross - 1] = wx[lay_cross - 1];
        wx[lay_cross - 1] = 0.6;
      }
      // D = D + xl[2];
      // if (first_lay == 3)
      // X[0] = 0.5 * xl[2];
      // else if (first_lay == 2)
      //   X[1] = xl[1] + 0.5 * xl[2];
      // else if (first_lay == 1)
      //  X[2] = 2 * xl[0] + 0.5 * xl[2];
    }
    if (Lay3.e > 0) {
      lay_cross++;
      yx[lay_cross - 1] = Lay3.x;
      yy[lay_cross - 1] = Lay3.y;
      yz[lay_cross - 1] = Lay3.z;
      wx[lay_cross - 1] = 0.001 * Lay3.e;
      wy[lay_cross - 1] = 0.6;
      wz[lay_cross - 1] = 0.6;
      if (lay_cross == 1) {
        first_lay = 4;
      }
      if (isBarrel == false) {
        wy[lay_cross - 1] = wx[lay_cross - 1];
        wx[lay_cross - 1] = 0.6;
      }
      // D = D + xl[3];
      // if (first_lay == 4)
      //  X[0] = 0.5 * xl[3];
      // else if (first_lay == 3)
      //  X[1] = xl[2] + 0.5 * xl[3];
      // else if (first_lay == 2)
      //  X[2] = 2 * xl[0] + 0.5 * xl[3];
      // else if (first_lay == 1)
      //  X[3] = 3 * xl[0] + 0.5 * xl[3];
    }
    if (Lay4.e > 0) {
      lay_cross++;
      yx[lay_cross - 1] = Lay4.x;
      yy[lay_cross - 1] = Lay4.y;
      yz[lay_cross - 1] = Lay4.z;
      wx[lay_cross - 1] = 0.001 * Lay4.e;
      wy[lay_cross - 1] = 0.6;
      wz[lay_cross - 1] = 0.6;
      if (lay_cross == 1) {
        first_lay = 5;
      }
      if (isBarrel == false) {
        wy[lay_cross - 1] = wx[lay_cross - 1];
        wx[lay_cross - 1] = 0.6;
      }
      // D = D + xl[4];
      // if (first_lay == 5)
      //   X[0] = 0.5 * xl[4];
      // else if (first_lay == 4)
      //   X[1] = xl[3] + 0.5 * xl[4];
      // else if (first_lay == 3)
      //   X[2] = 2 * xl[3] + 0.5 * xl[4];
      // else if (first_lay == 2)
      //   X[3] = 3 * xl[3] + 0.5 * xl[4];
      // else if (first_lay == 1)
      //   X[4] = 4 * xl[3] + 0.5 * xl[4];
    }
    if (lay_cross == 0) { 
      continue;
    }

    infofile << "lay cross: " << lay_cross << std::endl; 
    infofile << "first layer: " << first_lay << std::endl; 

    // D = D - 0.5 * xl[lay_cross];
    int Q = 0, L = 0;
    for (int k = first_lay; k <= 5; k++) {
      if (Q == 0 && (LayE[k - 1] >= 0.05 * clu_vec.at(i).e)) { //or 0.5? 
        Q = k; // DC the first layer with E > 5/50 % of the total cluster energy
      }
      L++; // DC, L is the number of layers with some deposited energy. 
    }
    infofile << "Q: " << Q << std::endl; 
    double E1 = 0, E2 = 0;
    if (lay_cross > 1) {
      for (int k_i = 5; k_i >= Q; k_i--) {
        E2 = E1;
        E1 = E1 + LayE[k_i - 1];
      }
      infofile << "E1: " << E1 << ", E2: " << E2 <<std::endl; 

      double Rk = E2 / E1;
      double B = 3;
      if (clu_vec.at(i).e > 16.5) B = 1.5 / log(clu_vec.at(i).e / 10);
      double Zmin = 0;
      double Zapx = 0.5 * B * xl[Q - 1]; 
      double Zmax = B * xl[Q - 1];
      double R1 = 0;
      for (int j = 0; j < 4; j++) {
        R1 = exp(-Zapx) * (1 + Zapx);
        if (R1 > Rk) {
          Zmin = Zapx;
          Zapx = 0.5 * (Zmax + Zmin);
        } else if (R1 < Rk) {
          Zmax = Zapx;
          Zapx = 0.5 * (Zmax + Zmin);
        }
      }
      Zapx = -Zapx / B;
      for (int j = 0; j < Q; j++) {
        Zapx = Zapx + xl[j];
      }
      double XFix[5] = {0, 0, 0, 0, 0};
      if (first_lay == 1) {
        XFix[0] = 2.22;
        XFix[1] = 6.66;
        XFix[2] = 11.1;
        XFix[3] = 15.54;
        //XFix[4] = 18.16;
        XFix[4] = 20.38;
      } else if (first_lay == 2) {
        XFix[0] = 6.66;
        XFix[1] = 11.1;
        XFix[2] = 15.54;
        //XFix[3] = 18.16;
        XFix[3] = 20.38;
      } else if (first_lay == 3) {
        XFix[0] = 11.1;
        XFix[1] = 15.54;
        //XFix[2] = 18.16;
        XFix[2] = 20.38;
      } else if (first_lay == 4) {
        XFix[0] = 15.54;
        XFix[1] = 20.38;
        //XFix[1] = 18.16;
      } else if (first_lay == 5) {
        //XFix[0] = 18.16;
        XFix[0] = 20.38;
      }
      for (int j = 0; j < lay_cross; j++) {
        X[j] = XFix[j] - Zapx; //apex position in module height direction???
        infofile << "X[" << j << "]: " << X[j] << ", yx[" << j << "] : " << yx[j] << ", yy[" << j << "] : " << yy[j] << ", yz[" << j << "] : " << yz[j] << ", Zapx: " << Zapx << ", wx[" << j << "] : " << wx[j] << ", wy[" << j << "] : " << wy[j] << ", wz[" << j << "] : " << wz[j] << std::endl; 
      }

      infofile << "fit_ls x: " << std::endl; 
      std::tuple<double, double, double, double> fit_varx =
          fit_ls(lay_cross, X, yx, wx); //num layer crossed, apex position in module height direction??, positions of layer centroids in x coordinate, constant to understand
      infofile << "fit_ls y: " << std::endl;
      std::tuple<double, double, double, double> fit_vary =
          fit_ls(lay_cross, X, yy, wy);
          infofile << "fit_ls z: " << std::endl;
      std::tuple<double, double, double, double> fit_varz =
          fit_ls(lay_cross, X, yz, wz);
      double trktot = sqrt(get<1>(fit_varx) * get<1>(fit_varx) +
                           get<1>(fit_vary) * get<1>(fit_vary) +
                           get<1>(fit_varz) * get<1>(fit_varz));
      ctrk[0] = get<1>(fit_varx) / trktot;
      ctrk[1] = get<1>(fit_vary) / trktot;
      ctrk[2] = get<1>(fit_varz) / trktot;
      ectrk[0] = get<3>(fit_varx) / trktot;
      ectrk[1] = get<3>(fit_vary) / trktot;
      ectrk[2] = get<3>(fit_varz) / trktot;
      apx[0] = get<0>(fit_varx);
      apx[1] = get<0>(fit_vary);
      apx[2] = get<0>(fit_varz);
      eapx[0] = get<2>(fit_varx);
      eapx[1] = get<2>(fit_vary);
      eapx[2] = get<2>(fit_varz);
      // cout << "Apx: [" << apx[0] << "; " << apx[1] << "; " << apx[2] << "] "
      //<< " e direzione vettore: [" << ctrk[0] << "; " << ctrk[1] << "; " <<
      // ctrk[2] << "] " << endl;
    }
    if (lay_cross == 1) {
      // cout << "apx(x)=" << yx[0] << " apx(y)=" << yy[0] << " apx(z)=" <<
      // yz[0] << endl;
      // cout << "e_apx(x)=" << sqrt(1/wx[0]) << " e_apx(y)=" << sqrt(1 / wy[0])
      //<< " e_apx(z)=" << sqrt(1 / wz[0]) << endl;
      apx[0] = yx[0];
      apx[1] = yy[0];
      apx[2] = yz[0];
    }

    // DISCLAIMER: se vogliamo possiamo utilizzare l'apice del cluster, basta
    // scommentare 
    
    clu_vec.at(i).x = apx[0]; 
    clu_vec.at(i).y = apx[1];
    clu_vec.at(i).z = apx[2];

    clu_vec.at(i).sx = ctrk[0];
    clu_vec.at(i).sy = ctrk[1];
    clu_vec.at(i).sz = ctrk[2];
  }
  return clu_vec;
}

std::tuple<double, double, double, double> fit_ls(int lay, double* X,double* Y, double* W)
{
  double norm = 0, xa = 0, ya = 0, xya = 0, x2a = 0;
  double det, A, B, dA, dB;
  for (int i = 0; i < lay; i++) {
    norm = norm + W[i]; 
    xa = xa + W[i] * X[i];
    ya = ya + W[i] * Y[i];
    xya = xya + W[i] * Y[i] * X[i];
    x2a = x2a + W[i] * X[i] * X[i];
  }
  norm = norm / lay;
  xa = xa / lay;
  ya = ya / lay;
  xya = xya / lay;
  x2a = x2a / lay;
  det = x2a * norm - xa * xa;
  B = (norm * xya - xa * ya) / (norm * x2a - xa * xa);
  A = ya / norm - B * xa / norm;
  dB = 1 / sqrt(lay * det);
  dA = sqrt(x2a / (lay * det));
  infofile << "xa: " << xa << "ya: "<< ya << std::endl; 
  infofile << "A : " << A <<  ", B : " << B << ", dA : " << dA <<  ", dB : " << dB <<  std::endl; 
  return std::make_tuple(A, B, dA, dB);
}

//Calcolo variabili solo per il primo pulse 

cluster Calc_variables(std::vector<dg_cell> cells)
{
  ////infofile<<"Calc_variables" << std::endl;
  double x_weighted = 0, y_weighted = 0, z_weighted = 0, t_weighted = 0,
         x2_weighted = 0, y2_weighted = 0, z2_weighted = 0, Etot = 0, E2tot,
         EvEtot = 0, EA, EAtot = 0, EB, EBtot = 0, TA = 0, TB = 0;
  for (auto const& cell : cells) {
    double d = DfromTDC(cell.ps1.at(0).tdc, cell.ps2.at(0).tdc); // distance module center - photosignal 
    double d1, d2, d3;
    d1 = 0.5 * cell.l + d; // distance pmt1 - module center -> middle of the layer + distance module center-hit 
    d2 = 0.5 * cell.l - d;
    double cell_E = kloe_simu::EfromADC(cell.ps1.at(0).adc, cell.ps2.at(0).adc,
                                        d1, d2, cell.lay);
    double cell_T =
        kloe_simu::TfromTDC(cell.ps1.at(0).tdc, cell.ps2.at(0).tdc, cell.l); // time of the photosignal 
        ////infofile<<"time of the photosignal: " << cell_T << "in the cell with id: " << cell.id << std::endl;
    if (cell.mod > 25) { //endcap  
      d3 = cell.y - d;
      y_weighted = y_weighted + (d3 * cell_E);
      y2_weighted = y2_weighted + (d3 * d3 * cell_E);
      x_weighted = x_weighted + (cell.x * cell_E);
      x2_weighted = x2_weighted + (cell.x * cell.x * cell_E);
    } else { //barrel 
      d3 = cell.x - d; //posizione ricostruita dell'hit lungo il modulo a partire dal centro, cell.x = module.x !  
      x_weighted = x_weighted + (d3 * cell_E); 
      x2_weighted = x2_weighted + (d3 * d3 * cell_E);
      y_weighted = y_weighted + (cell.y * cell_E);
      y2_weighted = y2_weighted + (cell.y * cell.y * cell_E);
    }
    t_weighted = t_weighted + cell_T * cell_E;
    ////infofile<< "t_weighted: " << t_weighted << std::endl; 
    z_weighted = z_weighted + (cell.z * cell_E);
    z2_weighted = z2_weighted + (cell.z * cell.z * cell_E);
    Etot = Etot + cell_E;
    E2tot = E2tot + cell_E * cell_E;
  }
  x_weighted = x_weighted / Etot;
  x2_weighted = x2_weighted / Etot;
  y_weighted = y_weighted / Etot;
  y2_weighted = y2_weighted / Etot;
  z_weighted = z_weighted / Etot;
  z2_weighted = z2_weighted / Etot;
  t_weighted = t_weighted / Etot;
  if (x_weighted > -0.000001 && x_weighted < 0.000001) x_weighted = 0;
  if (y_weighted > -0.000001 && y_weighted < 0.000001) y_weighted = 0;
  if (z_weighted > -0.000001 && z_weighted < 0.000001) z_weighted = 0;
  double dx, dy, dz;
  double neff = Etot * Etot / E2tot;
  double dum = neff / (neff - 1);
  if (cells.size() == 1) {
    dx = 0;
    dy = 0;
    dz = 0;
  } else {
    if (x2_weighted - x_weighted * x_weighted < 0) {
      dx = 0;
    } else {
      dx = sqrt(dum * (x2_weighted - x_weighted * x_weighted));
    }
    if (y2_weighted - y_weighted * y_weighted < 0) {
      dy = 0;
    } else {
      dy = sqrt(dum * (y2_weighted - y_weighted * y_weighted));
    }
    if (z2_weighted - z_weighted * z_weighted < 0) {
      dz = 0;
    } else {
      dz = sqrt(dum * (z2_weighted - z_weighted * z_weighted));
    }
  }
  cluster clust;
  clust.e = Etot;
  clust.x = x_weighted;
  clust.y = y_weighted;
  clust.z = z_weighted;
  clust.t = t_weighted;
  clust.varx = dx;
  clust.vary = dy;
  clust.varz = dz;
  clust.cells = cells;
  return clust;
}

bool RepetitionCheck(std::vector<int> v, int check)
{
  if (std::find(v.begin(), v.end(), check) != v.end()) {
    return true;
  } else {
    return false;
  }
}

bool isNeighbour(int id, int c_id)
{
  // Middle of the module
  if (c_id == id + 101 || c_id == id + 100 || c_id == id + 99 ||
      c_id == id + 1 || c_id == id - 1 || c_id == id - 99 || c_id == id - 100 ||
      c_id == id - 101) {
    return true;
  }
  // Right edge of the module
  else if (c_id == id - 889 || c_id == id - 989 || c_id == id - 1089) {
    return true;
  }
  // Right edge of 0 module
  else if ((c_id == id + 23111 || c_id == id + 23011 || c_id == id + 22911) &&
           id < 25000) {
    return true;
  }
  // Left edge of the module
  else if (c_id == id + 1089 || c_id == id + 989 || c_id == id + 889) {
    return true;
  }
  // Left edge of the 23 module
  else if ((c_id == id - 23011 || c_id == id - 22911 || c_id == id - 23111) &&
           id < 25000) {
    return true;
  }
  // Multiple hit on same cell 
  else if (c_id == id) {
    return true;
  } else
    return false;
}



// std::vector<int> GetNeighCells(int id)
// {
//     std::string strNumber = std::to_string(id);
//     int module;
//     int layer;
//     int cell;
//     // Extract substrings
//     if(id < 1000){
//         module = std::stoi("00");
//         layer = std::stoi(strNumber.substr(0, 1));
//         cell = std::stoi(strNumber.substr(1, 2));
        
//     } else {
//         module = std::stoi(strNumber.substr(0, 2));
//         layer = std::stoi(strNumber.substr(2, 1));
//         cell = std::stoi(strNumber.substr(3, 2));
//     }
//     std::vector<int> NCells;
    
//     // Middle of the module
//     if(cell < 11 && cell > 0){
//         NCells.push_back(id+101);
//         NCells.push_back(id+100);
//         NCells.push_back(id+99);
//         NCells.push_back(id+1);
//         NCells.push_back(id - 1);
//         NCells.push_back(id-99);
//         NCells.push_back(id-100);
//         NCells.push_back(id-101);
//     }
    
    
//     if(module != 0 && module != 23){
//         // Right edge of the module
//         if(cell == 0){
//             NCells.push_back(id-889);
//             NCells.push_back(id - 989);
//             NCells.push_back(id - 1089);
//         }
//         // Left edge of the module
//         else if (cell == 11 ) {
//             NCells.push_back(id+889);
//             NCells.push_back(id + 989);
//             NCells.push_back(id + 1089);
//         }
//     } else if (module == 0 || module == 23){
//         // Right edge of 0 module
//         if (cell == 11) {
//             NCells.push_back(id+23111);
//             NCells.push_back(id+23011);
//             NCells.push_back(id+22911);
//         } 
//         // Left edge of the 23 module
//         else if(cell == 0){
//            NCells.push_back(id-23111);
//             NCells.push_back(id-23011);
//             NCells.push_back(id-22911);
//         }
//     }
//     return NCells
// }








std::pair<std::vector<dg_cell>, std::vector<int>> GetNeighbours(
    std::vector<dg_cell> cells, int start, std::vector<int> checked,
    std::vector<dg_cell> neigh_chain)
{
  for (int i = 0; i < cells.size(); i++) {
    // //infofile<<"Id Cells start: " << cells.at(start).id << ", id Cell i: " << cells.at(i).id <<std::endl;
    // //infofile<< "Position in original multicell vector -> start: " << start << " , i: " << i << std::endl;
    if (RepetitionCheck(checked, i) == true) continue;
    bool check = isNeighbour(cells.at(start).id, cells.at(i).id);
    
    if(check ==true)
    // //infofile<< cells.at(start).id << " isNeighbour of cell " << cells.at(i).id <<std::endl;
    
    if (check == true) {
      neigh_chain.push_back(cells.at(i));
      checked.push_back(i);
      // //infofile<< "Cell with id: " << cells.at(i).id << " was put in the checked vector." << std::endl;
      std::pair<std::vector<dg_cell>, std::vector<int>> find_chain = GetNeighbours(cells, i, checked, neigh_chain);
      neigh_chain = find_chain.first;
      checked = find_chain.second;
    }
  }
  return std::make_pair(neigh_chain, checked);
}

double kloe_simu::AttenuationFactor(double d, int planeID)
{
  /*
       dE/dx attenuation - Ea=p1*exp(-d/atl1)+(1.-p1)*exp(-d/atl2)
         d    distance from photocatode - 2 cells/cell; d1 and d2
        atl1  50. cm
        atl2  430 cm planes 1-2    innermost plane is 1
              380 cm plane 3
              330 cm planes 4-5
         p1   0.35
  */
  double atl2 = 0.0;

  switch (planeID) {
    case 0:
    atl2 = kloe_simu::atl2_01;
    break; //TO CHECK DC 
    case 1:
    atl2 = kloe_simu::atl2_01;
    break;
    
    case 2:
    atl2 = kloe_simu::atl2_2;
    break;
    
    case 3:
    atl2 = kloe_simu::atl2_34;
    break;
    
    case 4:
    atl2 = kloe_simu::atl2_34;
    break;
    
    default:
    // std::cout << "planeID out if range" << std::endl;
      atl2 = -999.0;
      break;
  }
  
  return kloe_simu::p1 * TMath::Exp(-d / kloe_simu::atl1) +
         (1. - kloe_simu::p1) * TMath::Exp(-d / atl2);
}

// reconstruct t of the hit from tdc1 and tdc2
double kloe_simu::TfromTDC(double t1, double t2, double L)
{
  return 0.5 * (t1 + t2 - kloe_simu::vlfb * L / kloe_simu::m_to_mm);
}

// energy deposit of the hit from adc1 and adc2 and
// reconstructed longidutinal coordinate
double kloe_simu::EfromADC(double adc1, double adc2, double d1, double d2,
                           int planeID)
{
  double f1 = AttenuationFactor(d1, planeID);
  double f2 = AttenuationFactor(d2, planeID);

  double const attpassratio = 0.187;
  return 0.5 * (adc1 / f1 + adc2 / f2) /
         (attpassratio * kloe_simu::pe2ADC * kloe_simu::e2p2);
}

double EfromADCsingle(double adc, double f)
{
  double const attpassratio = 0.187;
  return adc / (f * attpassratio * kloe_simu::pe2ADC * kloe_simu::e2p2);
}

double DfromTDC(double ta, double tb)
{
  return 0.5 * (ta - tb) / kloe_simu::vlfb * kloe_simu::m_to_mm;
}

bool endsWith(const std::string& fullString, const std::string& ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}
bool initializeFiles(int argc, char* argv[], std::string& digitFileName) {
    // Process command-line arguments to build the file names
    for (int i = 1; i < argc; i += 2) {
        std::string flag = argv[i];
        std::string fileName = argv[i + 1];

        if (flag == "-d") {
            digitFileName = fileName;
        } else {
            std::cerr << "Error: Unknown flag: " << flag << std::endl;
            return false;
        }
    }

    // Check if the required files are present
    if (digitFileName.empty()) {
        std::cerr << "Error: Missing required flags. Please use '-d' with corresponding file name." << std::endl;
        return false;
    }

    // Validate file names
    if (!endsWith(digitFileName, "edep.digit.root")) {
        std::cerr << "Error: Invalid arguments. Please use '-d' before the edep.digit.root file" << std::endl;
        return false;
    }

    return true; // Files initialized successfully
}



int main(int argc, char* argv[]){
     if (argc < 2) {
        std::cerr << "Usage: Clustering" << " -d <file.edep.digit.root>" << std::endl;
        return 1; // Exit with an error code
    }
    string digitFileName;
    if (!initializeFiles(argc, argv, digitFileName)) {
        return 1; // Exit with an error code
    }
    
    clustering(digitFileName);
  return 0;
}