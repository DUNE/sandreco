#include <vector>
#include <utility>
#include <string>
#include <complex>
#include <tuple>
#include <iostream>
#include <fstream>
#include "TLeaf.h"
#include "TFile.h"
#include "TTree.h"

#include "struct2.h"
#include "utils.h"
#include "Linkdef.h"


using namespace std;

bool RepetitionCheck(std::vector<int>, int);

std::pair<vector<int>, vector<int>>FindNeighbours(int[], std::vector<int>, int, int, std::vector<int>);

std::pair<int, int> Next(std::vector<int>, int[], int);
std::tuple<double, double, double, double> fit_ls(int, double[], double[], double[]);
void TrackFit(std::vector<cluster2>);
void Clust_info(cluster2);
std::vector<cluster2> Merge(std::vector<cluster2>);
std::vector<cluster2> Split(std::vector<cluster2>);
std::vector<cluster2> RecoverIncomplete(std::vector<cluster2>, std::vector<dg_cell>);
cluster2 Calc_variables(std::vector<dg_cell>);

double TfromTDC(double t1, double t2, double L);
double AttenuationFactor(double d, int planeID);
double EfromADC(double adc1, double adc2, double d1, double d2, int planeID);
double DfromADC(double, double);

std::pair<std::vector<cluster2>, std::vector<cluster2>> Preclustering(std::vector<dg_cell> *vec_cellraw) {
    int dg_size=vec_cellraw->size();
    int dg_idvec[dg_size],i_id=0;
    std::vector<dg_cell> incomplete_cells, complete_cells;
    std::vector<cluster2> vec_clust;
    //Create vector of complete and incomplete cells
    for (int i = 0; i < vec_cellraw->size(); i++) {
        if (vec_cellraw->at(i).adc1 == 0 || vec_cellraw->at(i).adc2 == 0 || vec_cellraw->at(i).tdc1 == 0 || vec_cellraw->at(i).tdc1 == 0) {
            //Incomplete cell
            incomplete_cells.push_back(vec_cellraw->at(i));
        }
        else {
            //Complete cell
            complete_cells.push_back(vec_cellraw->at(i));
            dg_idvec[i_id] = vec_cellraw->at(i).id;
            i_id++;
        }   
    }
    int cluster[i_id], n_cluster = 0;
    std::vector<int> checked_array;
    //Star preclustering on incomplete cells
    std::vector<int> vec_cell;
    for (int j = 0; j < complete_cells.size(); j++) {
        if (j == 0) {
            cluster[j] = dg_idvec[0];
        }
        else {
            std::pair<int, int> clu_nclu = Next(checked_array, dg_idvec, complete_cells.size());
            cluster[j] = clu_nclu.first;
            if (cluster[j] == 0) {
                break;
            }
            n_cluster = clu_nclu.second;
        }
        vec_cell.push_back(n_cluster);
        checked_array.push_back(cluster[j]);
        std::pair < vector<int>, vector<int>> check_cell = FindNeighbours(dg_idvec, checked_array, i_id, cluster[j], vec_cell);
        checked_array = check_cell.first;
        vec_cell = check_cell.second;
        std::vector<int> vec_cellid;
        for (std::vector<int>::const_iterator it = vec_cell.begin(); it != vec_cell.end(); it++) {
            vec_cellid.push_back(*it);
        }
        std::vector<dg_cell> cluster_cells;
        for (int its = 0; its < vec_cellid.size(); its++) {
            int cell = vec_cellid[its];
            cluster_cells.push_back(complete_cells.at(cell));
        }
        vec_cellid.clear();
        cluster2 clust;
        clust = Calc_variables(cluster_cells);
        vec_clust.push_back(clust);
        vec_cell.clear();
    }
    //SPLIT
    vec_clust = Split(vec_clust);
    //MERGE
    vec_clust = Merge(vec_clust);
    for (int cl = 0; cl < vec_clust.size(); cl++) {
        //Clust_info(vec_clust.at(cl));
    }
    //cout << "=O=O=O=O=O=O=O=O=O=O=O=O=O" << endl;
    //cout << "Incomplete cells: " << incomplete_cells.size() << endl;
    cout << endl;
    //cout << "--------------------------" << endl;
    //Track Fit
    TrackFit(vec_clust);
    std::vector<cluster2> Og_clus = vec_clust;
    vec_clust = RecoverIncomplete(vec_clust, incomplete_cells);
    for (int cl = 0; cl < vec_clust.size(); cl++) {
        //Clust_info(vec_clust.at(cl));
    }
    return std::make_pair(vec_clust, Og_clus);
}

int Testing(std::string input) {
    const char* finname = input.c_str();
    TFile f(finname, "READ");
    TTree* t = (TTree*)f.Get("tDigit");
    int TotHits = t->GetEntries();
    std::vector<dg_cell>* cell = new std::vector<dg_cell>;
    std::vector<cluster2> clust, og_clust;
    TH1F* ClusterEnergy = new TH1F("ClusterEnergy", "Total Energy reconstructed in the ECal", 100, 0., 3000.);
    TH1F* Og_ClusterEnergy = new TH1F("ClusterEnergy", "Total Energy reconstructed-corrected in the ECal", 100, 0., 3000.);
    TH1F* Xdistri = new TH1F("Xdistribution", "X distribution", 100, -2150., 2150.);
    t->SetBranchAddress("dg_cell", &cell);
    for (int i = 0; i < TotHits; i++) {
        cout << "--------------" << endl;
        cout << "Entry " << i << endl;
        t->GetEntry(i);
        std::pair<std::vector<cluster2>, std::vector<cluster2>>clustered=Preclustering(cell);
        clust = clustered.first;
        og_clust = clustered.second;
        double CluEn = 0, Og_CluEn=0, Xdis=0;
        for (int j = 0; j < clust.size(); j++) {
            CluEn = CluEn + clust.at(j).e;
            Og_CluEn = Og_CluEn + og_clust.at(j).e;
            Xdis = Xdis + clust.at(j).x * clust.at(j).e;
        }
        Xdis = Xdis / CluEn;
        ClusterEnergy->Fill(CluEn);
        Og_ClusterEnergy->Fill(Og_CluEn);
        Xdistri->Fill(Xdis);
        clust.clear();
        og_clust.clear();
    }
    ClusterEnergy->SetDirectory(gROOT);
    Og_ClusterEnergy->SetDirectory(gROOT);
    Xdistri->SetDirectory(gROOT);
    //Og_ClusterEnergy->SetFillColor(kRed);
    //ClusterEnergy->SetFillColor(kBlue);
    ClusterEnergy->SetLineColorAlpha(kBlue, 1);
    Og_ClusterEnergy->SetLineColorAlpha(kRed, 1);
    ClusterEnergy->Draw();
    Og_ClusterEnergy->Draw("same");
    //Xdistri->Draw();
    return 0;
}

std::vector<cluster2> RecoverIncomplete(std::vector<cluster2> clus, std::vector<dg_cell> incomplete_cells) {
    for (int i = 0; i < incomplete_cells.size(); i++) {
        double cell_phi = atan((incomplete_cells.at(i).z - 23910.00) / (incomplete_cells.at(i).y + 2384.73)) * 180 / 3.14159265358979323846;
        for (int j = 0; j < clus.size(); j++) {
            double rec_en = 0;
            double clus_phi = atan((clus.at(j).z - 23910.00) / (clus.at(j).y + 2384.73)) * 180 / 3.14159265358979323846;
            if (abs(cell_phi - clus_phi) < 3 || abs(clus_phi)!=90) {
                cout << "Aggiungere la cella " << incomplete_cells.at(i).id << " al cluster " << j << endl;
                double DpmA = clus.at(j).x/10+215;
                double DpmB = -clus.at(j).x/10+215;
                cout << "DpmA " << DpmA << " DpmB  " << DpmB << endl;
                double AttA = 0.75 * exp(-DpmA / 450) + 0.35 * exp(-DpmA / 50);
                double AttB = 0.75 * exp(-DpmB / 450) + 0.35 * exp(-DpmB / 50);
                double Scala = 0.75 * exp(-215 / 450) + 0.35 * exp(-215 / 50);
                AttA = Scala * 0.9 / (AttA * 0.96);
                AttB = Scala * 0.9 / (AttB * 0.96);
                cout << "AttA " << AttA << " AttB " << AttB << endl;
                cout << "Cluster [Y,Z]= [" << clus.at(j).y << "," << clus.at(j).z << "] vs Cell [Y,Z]= [" << incomplete_cells.at(i).y << ", " << incomplete_cells.at(i).z << "]" << endl;
                cout << "Cluster Energy Pre: " << clus.at(j).e << endl;
                if (incomplete_cells.at(i).adc1 != 0 && incomplete_cells.at(i).adc2) {
                    double Ea = incomplete_cells.at(i).adc1;
                    double Eb = incomplete_cells.at(i).adc2;
                    rec_en = kloe_simu::EfromADC(Ea, Eb, DpmA, DpmB, incomplete_cells.at(i).lay);
                    cout << "EA e EB > 0, rec_En " << rec_en << endl;
                    clus.at(j).e = clus.at(j).e + rec_en;
                    clus.at(j).cells.push_back(incomplete_cells.at(i));
                }
                else if (incomplete_cells.at(i).adc1 != 0 && incomplete_cells.at(i).tdc1) {
                    rec_en = incomplete_cells.at(i).adc1 * AttA/0.145;
                    cout << "EA > 0, rec_En " << rec_en << endl;
                    clus.at(j).e = clus.at(j).e + rec_en;
                    clus.at(j).cells.push_back(incomplete_cells.at(i));
                }
                else if (incomplete_cells.at(i).adc2 != 0 && incomplete_cells.at(i).tdc2) {
                    rec_en = incomplete_cells.at(i).adc2 * AttB / 0.145;
                    cout << "EB > 0, rec_En " << rec_en << endl;
                    clus.at(j).e = clus.at(j).e + rec_en;
                    clus.at(j).cells.push_back(incomplete_cells.at(i));
                }
                cout << "Cluster Energy Post: " << clus.at(j).e << endl;
                continue;
            }
        }
    }
    return clus;
}

void Clust_info(cluster2 clus) {
    cout << "=O=O=O=O=O=O=O=O=O=O=O=O=O" << endl;
    cout << "New cluster: Energy " << clus.e << " MeV" << endl;
    cout << "Coordinate centroide: " << clus.x << " [X] " << clus.y << " [Y] " << clus.z << " [z] e tempo di arrivo medio " << clus.t << " ns" << endl;
    double phi = atan((clus.z - 23910.00)/ (clus.y + 2384.73)) * 180 / 3.14159265358979323846;
    cout << "Angolo azimutale: " << phi << endl;
    cout << "Composto dalle seguenti celle: ";
    for (int i = 0; i < clus.cells.size(); i++) {
        cout << clus.cells.at(i).id << " ";
    }
    cout << endl;
    
}

std::vector<cluster2> Split(std::vector<cluster2> original_clu_vec) {
    std::vector<cluster2> clu_vec;
    std::vector <dg_cell> all_cells;
    cout << "Presplit size: " << original_clu_vec.size() << endl;
    for (int k = 0; k < original_clu_vec.size(); k++) {
        //cout << "Cluster " << k << " con " << original_clu_vec.at(k).cells.size() << " celle ed energia " << original_clu_vec.at(k).e << endl;
        int splitted = 0;
        double tA = 0, tB = 0, tA2 = 0, tB2 = 0;
        double EA, EB, EAtot = 0, EBtot = 0, EA2tot = 0, EB2tot = 0;
        double tRMS_A, tRMS_B, dist;
        all_cells = original_clu_vec.at(k).cells;
        for (int j = 0; j < all_cells.size(); j++) {
            EA = all_cells.at(j).adc1;
            EB = all_cells.at(j).adc2;
            EAtot += EA;
            EA2tot += EA * EA;
            EBtot += EB;
            EB2tot += EB * EB;
            double d = DfromADC(all_cells[j].tdc1, all_cells[j].tdc2);
            double d1, d2, d3;
            d1 = 0.5 * all_cells[j].l + d;
            d2 = 0.5 * all_cells[j].l - d;
            double cell_E = kloe_simu::EfromADC(all_cells[j].adc1, all_cells[j].adc2, d1, d2, all_cells[j].lay);
            //double d1 = DfromADC(all_cells.at(j).tdc1, all_cells.at(j).tdc2);
            //double d2;
            //if (all_cells.at(j).tdc1 <= all_cells.at(j).tdc2) {
            //    d1 = 0.5 * (all_cells.at(j).l - d1);
            //    d2 = 0.5 * (all_cells.at(j).l + d1);
           // }
            //else {
             //   d2 = 0.5 * (all_cells.at(j).l - d1);
             //   d1 = 0.5 * (all_cells.at(j).l + d1);
            //}
            tA += (all_cells.at(j).tdc1 - kloe_simu::vlfb * d1 / kloe_simu::m_to_mm) * EA;
            tA2 += std::pow(all_cells.at(j).tdc1 - kloe_simu::vlfb * d1 / kloe_simu::m_to_mm, 2) * EA;
            tB += (all_cells.at(j).tdc2 - kloe_simu::vlfb * d2 / kloe_simu::m_to_mm) * EB;
            tB2 += std::pow(all_cells.at(j).tdc2 - kloe_simu::vlfb * d2 / kloe_simu::m_to_mm, 2) * EB;
        }
        tA = tA / EAtot;
        tA2 = tA2 / EAtot;
        tB = tB / EAtot;
        tB2 = tB2 / EBtot;
        tRMS_A = (tA2 - tA * tA) * (EA2tot - EAtot * EAtot) / EA2tot;
        tRMS_B = (tB2 - tB * tB) * (EB2tot - EBtot * EBtot) / EB2tot;
        dist = std::sqrt(tRMS_A * tRMS_A + tRMS_B * tRMS_B);
        if (dist > 5) {
            std::vector <dg_cell> q1_cells, q2_cells, q3_cells, q4_cells;
            for (unsigned int i = 0; i < all_cells.size(); i++) {
                double t_difA = all_cells.at(i).tdc1 - tA;
                double t_difB = all_cells.at(i).tdc2 - tB;
                if (t_difA > 0) {
                    if (t_difB > 0) { q1_cells.push_back(all_cells.at(i)); }
                    else { q2_cells.push_back(all_cells.at(i)); }
                }
                else {
                    if (t_difB > 0) { q3_cells.push_back(all_cells.at(i)); }
                    else { q4_cells.push_back(all_cells.at(i)); }
                }
            }
            if (q1_cells.size() != 0) {
                //cout << "Splitting su quadrante 1" << endl;
                cluster2 clus = Calc_variables(q1_cells);
                //cout << "SPLIT: Cluster con " << clus.cells.size() << " celle ed energia " << clus.e << endl;
                clu_vec.push_back(clus);
                splitted++;
            }
            if (q2_cells.size() != 0) {
                //cout << "Splitting su quadrante 2" << endl;
                cluster2 clus = Calc_variables(q2_cells);
                //cout << "SPLIT: Cluster con " << clus.cells.size() << " celle ed energia " << clus.e << endl;
                clu_vec.push_back(clus);
                splitted++;
            }
            if (q3_cells.size() != 0) {
                //cout << "Splitting su quadrante 3" << endl;
                cluster2 clus = Calc_variables(q3_cells);
                //cout << "SPLIT: Cluster con " << clus.cells.size() << " celle ed energia " << clus.e << endl;
                clu_vec.push_back(clus);
                splitted++;
            }
            if (q4_cells.size() != 0) {
                //cout << "Splitting su quadrante 4" << endl;
                cluster2 clus = Calc_variables(q4_cells);
                //cout << "SPLIT: Cluster con " << clus.cells.size() << " celle ed energia " << clus.e << endl;
                clu_vec.push_back(clus);
                splitted++;
            }
            //if (splitted >2) cout << "We got " << splitted << " splits." << endl;
            q1_cells.clear();
            q2_cells.clear();
            q3_cells.clear();
            q4_cells.clear();
        }
        else {
            clu_vec.push_back(original_clu_vec.at(k));
            //cout << "// Cluster con " << original_clu_vec.at(k).cells.size() << " celle ed energia " << original_clu_vec.at(k).e << endl;
        }
        all_cells.clear();
    }
    original_clu_vec.clear();
    cout << "Post-Split size: " << clu_vec.size() << endl;
    return clu_vec;
}

std::vector<cluster2> Merge(std::vector<cluster2> Og_cluster) {
    std::vector<cluster2> mgd_cluster;
    std::vector<int> checked;
    for (int i = 0; i < Og_cluster.size(); i++) {
        double xi = Og_cluster.at(i).x;
        double yi = Og_cluster.at(i).y;
        double zi = Og_cluster.at(i).z;
        double ti = Og_cluster.at(i).t;
        double ei = Og_cluster.at(i).e;
        bool RepCheck = RepetitionCheck(checked, i);
        if (RepCheck == true) {
            continue;
        }
        checked.push_back(i);
        cluster2 clust;
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
            double D = sqrt((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj) + (zi - zj) * (zi - zj));
            double DT = abs(ti - tj);
            if (D < 40 && DT < 2.5) {
                bool endcap = false;
                if (clust.cells[0].id > 25000) {
                    endcap = true;
                }
                if (endcap == true) {
                    //siamo sull'endcap
                    double Dz_ec = abs(yi - yj);
                    D = sqrt((xi - xj) * (xi - xj) + (zi - zj) * (zi - zj));
                    if (Dz_ec < 30 && D < 40) {
                        //Unisce il cluster j ad i
                        //cout << "Unire il cluster j ad i .2 endcap" << endl;
                        std::vector<dg_cell> vec_cells_j = Og_cluster.at(j).cells;
                        //loop su vec_cells_j e push_back su clust_cells
                        for (int k = 0; k < vec_cells_j.size(); k++) {
                            clust.cells.push_back(vec_cells_j.at(k));
                        }
                        clust = Calc_variables(clust.cells);
                        checked.push_back(j);
                    }
                }
                else if (endcap == false) {
                    double Dz_bar = abs(xi - xj);
                    D = sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj));
                    if (Dz_bar < 30 && D < 40) {
                        //Unisce il cluster j ad i
                        //cout << "Unire il cluster j ad i .3 barrel" << endl;
                        std::vector<dg_cell> vec_cells_j = Og_cluster.at(j).cells;
                        //loop su vec_cells_j e push_back su clust_cells
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
    cout << "Post-Merge size: " << mgd_cluster.size() << endl;
    return mgd_cluster;
}

void TrackFit(std::vector<cluster2> clu_vec) {
    const double xl[5] = { 4.44, 4.44, 4.44, 4.44, 5.24 };
    for (int i = 0; i < clu_vec.size(); i++) {
        double apx[3], eapx[3] = { 0,0,0 }, ctrk[3], ectrk[3] = {0,0,0};
        std::vector<dg_cell> cell_vec_0, cell_vec_1, cell_vec_2, cell_vec_3, cell_vec_4;
        for (int j = 0; j < clu_vec.at(i).cells.size(); j++) {
            if (clu_vec.at(i).cells.at(j).lay == 0) {
                //aggiungere la cella al layer 0
                cell_vec_0.push_back(clu_vec.at(i).cells.at(j));
            }
            else if (clu_vec.at(i).cells.at(j).lay == 1) {
                //aggiungere la cella al layer 1
                cell_vec_1.push_back(clu_vec.at(i).cells.at(j));
            }
            else if (clu_vec.at(i).cells.at(j).lay == 2) {
                //aggiungere la cella al layer 2
                cell_vec_2.push_back(clu_vec.at(i).cells.at(j));
            }
            else if (clu_vec.at(i).cells.at(j).lay == 3) {
                //aggiungere la cella al layer 3
                cell_vec_3.push_back(clu_vec.at(i).cells.at(j));
            }
            else if (clu_vec.at(i).cells.at(j).lay == 4) {
                //aggiungere la cella al layer 4
                cell_vec_4.push_back(clu_vec.at(i).cells.at(j));
            }
        }
        cluster2 Lay0, Lay1, Lay2, Lay3, Lay4;
        Lay0 = Calc_variables(cell_vec_0);
        Lay1 = Calc_variables(cell_vec_1);
        Lay2 = Calc_variables(cell_vec_2);
        Lay3 = Calc_variables(cell_vec_3);
        Lay4 = Calc_variables(cell_vec_4);
        double LayE[5] = { Lay0.e,Lay1.e, Lay2.e, Lay3.e, Lay4.e };
        bool isBarrel = true;
        double yx[5] = { 0,0,0,0,0 }, yy[5] = { 0,0,0,0,0 }, yz[5] = { 0,0,0,0,0 }, wx[5] = { 0,0,0,0,0 }, wy[5] = { 0,0,0,0,0 }, wz[5] = { 0,0,0,0,0 };
        double X[5] = { 0,0,0,0,0 }, D;
        if (clu_vec.at(i).cells[0].id > 25000) {
            isBarrel = false;
        }
        int lay_cross = 0, first_lay = 0;
        if (Lay0.e > 0) {
            lay_cross++;
            if (lay_cross == 1) {
                first_lay = 1;
                //cout << "/=/=/=/=/=/=/=/=/=/=" << endl;
                //cout << Lay0.x << " " << Lay0.y << " " << Lay0.z << endl;
                //cout << "/=/=/=/=/=/=/=/=/=/=" << endl;
            }
            yx[lay_cross-1] = Lay0.x;
            wx[lay_cross - 1] = 0.6;
            wy[lay_cross - 1] = 0.6;
            wz[lay_cross - 1] = 0.001 * Lay0.e;
            yy[lay_cross - 1] = Lay0.y;
            yz[lay_cross - 1] = Lay0.z;
            if (isBarrel == false) {
                wy[lay_cross - 1] = wz[lay_cross - 1];
                wz[lay_cross - 1] = 0.6;
            }
            D = D + xl[0];
            if (first_lay == 1) X[0] = 0.5 * xl[0];
        }
        if (Lay1.e > 0) {
            lay_cross++;
            yx[lay_cross - 1] = Lay1.x;
            yy[lay_cross - 1] = Lay1.y;
            yz[lay_cross - 1] = Lay1.z;
            wx[lay_cross - 1] = 0.6;
            wy[lay_cross - 1] = 0.6;
            wz[lay_cross - 1] = 0.001 * Lay1.e;
            if (lay_cross == 1) {
                first_lay = 2;
                //cout << "/=/=/=/=/=/=/=/=/=/=" << endl;
                //cout << Lay1.x << " " << Lay1.y << " " << Lay1.z << endl;
                //cout << "/=/=/=/=/=/=/=/=/=/=" << endl;
            }
            if (isBarrel == false) {
                wy[lay_cross - 1] = wz[lay_cross - 1];
                wz[lay_cross - 1] = 0.6;
            }
            D = D + xl[1];
            if (first_lay == 2) X[0] = 0.5 * xl[1];
            else if (first_lay == 1) X[1] = xl[0] + 0.5 * xl[1];
        }
        if (Lay2.e > 0) {
            lay_cross++;
            yx[lay_cross - 1] = Lay2.x;
            yy[lay_cross - 1] = Lay2.y;
            yz[lay_cross - 1] = Lay2.z;
            wx[lay_cross - 1] = 0.6;
            wy[lay_cross - 1] = 0.6;
            wz[lay_cross - 1] = 0.001 * Lay2.e;
            if (lay_cross == 1) {
                first_lay = 3;
                //cout << "/=/=/=/=/=/=/=/=/=/=" << endl;
                //cout << Lay2.x << " " << Lay2.y << " " << Lay2.z << endl;
                //cout << "/=/=/=/=/=/=/=/=/=/=" << endl;
            }
            if (isBarrel == false) {
                wy[lay_cross - 1] = wz[lay_cross - 1];
                wz[lay_cross - 1] = 0.6;
            }
            D = D + xl[2];
            if (first_lay == 3) X[0] = 0.5 * xl[2];
            else if (first_lay == 2) X[1] = xl[0] + 0.5 * xl[2];
            else if (first_lay == 1) X[2] = 2 * xl[0] + 0.5 * xl[3];
        }
        if (Lay3.e > 0) {
            lay_cross++;
            yx[lay_cross - 1] = Lay3.x;
            yy[lay_cross - 1] = Lay3.y;
            yz[lay_cross - 1] = Lay3.z;
            wx[lay_cross - 1] = 0.6;
            wy[lay_cross - 1] = 0.6;
            wz[lay_cross - 1] = 0.001 * Lay3.e;
            if (lay_cross == 1) {
                first_lay = 4;
                //cout << "/=/=/=/=/=/=/=/=/=/=" << endl;
                //cout << Lay3.x << " " << Lay3.y << " " << Lay3.z << endl;
                //cout << "/=/=/=/=/=/=/=/=/=/=" << endl;
            }
            if (isBarrel == false) {
                wy[lay_cross - 1] = wz[lay_cross - 1];
                wz[lay_cross - 1] = 0.6;
            }
            D = D + xl[3];
            if (first_lay == 4) X[0] = 0.5 * xl[3];
            else if (first_lay == 3) X[1] = xl[0] + 0.5 * xl[3];
            else if (first_lay == 2) X[2] = 2 * xl[0] + 0.5 * xl[3];
            else if (first_lay == 1) X[3] = 3 * xl[0] + 0.5 * xl[3];
        }
        if (Lay4.e > 0) {
            lay_cross++;
            yx[lay_cross - 1] = Lay4.x;
            yy[lay_cross - 1] = Lay4.y;
            yz[lay_cross - 1] = Lay4.z;
            wx[lay_cross - 1] = 0.6;
            wy[lay_cross - 1] = 0.6;
            wz[lay_cross - 1] = 0.001 * Lay4.e;
            if (lay_cross == 1) {
                first_lay = 5;
                //cout << "/=/=/=/=/=/=/=/=/=/=" << endl;
                //cout << Lay4.x << " " << Lay4.y << " " << Lay4.z << endl;
                //cout << "/=/=/=/=/=/=/=/=/=/=" << endl;
            }
            if (isBarrel == false) {
                wy[lay_cross - 1] = wz[lay_cross - 1];
                wz[lay_cross - 1] = 0.6;
            }
            D = D + xl[4];
            if (first_lay == 5) X[0] = 0.5 * xl[4];
            else if (first_lay == 4) X[1] = xl[3] + 0.5 * xl[4];
            else if (first_lay == 3) X[2] = 2 * xl[3] + 0.5 * xl[4];
            else if (first_lay == 2) X[3] = 3 * xl[3] + 0.5 * xl[4];
            else if (first_lay == 1) X[4] = 4 * xl[3] + 0.5 * xl[4];
        }
        if (lay_cross == 0) {
            continue;
        }
        D = D - 0.5 * xl[lay_cross];
        int Q = 0, L=0;
        for (int k = first_lay; k <= 5; k++) {
            //cout << "Layer " << k - 1 << " Energy " << LayE[k - 1] <<" Xpos "<<yx[L]<< " Ypos " << yy[L] << " Zpos " << yz[L] << endl;
            if (Q == 0 && (LayE[k - 1] >= 0.05 * clu_vec.at(i).e)) {
                Q = k;
            }
            L++;
        }
        double E1 = 0, E2 = 0;
        if (lay_cross > 1) {
            for (int k_i = 5; k_i >= Q; k_i--) {
                E2 = E1;
                E1 = E1 + LayE[k_i - 1];
            }
            double Rk = E2 / E1;
            //cout << "Parte dal layer: " << first_lay - 1 << " che vede " << LayE[Q - 1] << " MeV su " << clu_vec.at(i).e << " MeV. Energia residua: " << E2 << " MeV" << endl;
            //cout << "Rk=" << Rk << endl;
            double B = 3;
            if (clu_vec.at(i).e > 16.5) B = 1.5 / log(clu_vec.at(i).e / 10);
            double Zmin = 0;
            double Zapx = 0.5 * B * xl[Q - 1];
            double Zmax = B * xl[Q - 1];
            //cout << "Provvisorio: Zapx=" << Zapx <<" B="<<B<<" and xl="<<xl[Q-1]<<endl;
            double R1;
            for (int j = 0; j < 4; j++) {
                R1 = exp(-Zapx) * (1 + Zapx);
             //   cout << "R1=" << R1 <<" Zapx="<<Zapx<< endl;
                if (R1 > Rk) {
                    Zmin = Zapx;
                    Zapx = 0.5 * (Zmax + Zmin);
                }
                else if (R1 < Rk) {
                    Zmax = Zapx;
                    Zapx = 0.5 * (Zmax + Zmin);
                }
            }
            //cout << "R1=" << R1 << " Zapx=" << Zapx << endl;
            Zapx = -Zapx / B;
            for (int j = 0; j < Q; j++) {
                Zapx = Zapx + xl[j];
             //   cout << Zapx << " " << xl[j] << endl;
            }
            for (int j = 0; j < lay_cross; j++) {
                X[j] = X[j] - Zapx;
            //    cout << "Zapx: " << X[j] << endl;
            }
            std::tuple<double, double, double, double> fit_varx = fit_ls(lay_cross, X, yx, wx);
            //cout <<"Ax="<< get<0>(fit_varx) << " Bx=" << get<1>(fit_varx)<< " dAx=" << get<2>(fit_varx)<<" dBx=" << get<3>(fit_varx)<< endl;
            std::tuple<double, double, double, double> fit_vary = fit_ls(lay_cross, X, yy, wy);
            //cout << "Ay=" << get<0>(fit_vary) << " By=" << get<1>(fit_vary) << " dAy=" << get<2>(fit_vary) << " dBy=" << get<3>(fit_vary) << endl;
            std::tuple<double, double, double, double> fit_varz = fit_ls(lay_cross, X, yz, wz);
            //cout << "Az=" << get<0>(fit_varz) << " Bz=" << get<1>(fit_varz) << " dAz=" << get<2>(fit_varz) << " dBz=" << get<3>(fit_varz) << endl;
            double trktot = sqrt(get<1>(fit_varx) * get<1>(fit_varx) + get<1>(fit_vary) * get<1>(fit_vary) + get<1>(fit_varz) * get<1>(fit_varz));
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
            cout << "Apx: [" << apx[0] << "; " << apx[1] << "; " << apx[2] << "] " << " e direzione vettore: [" << ctrk[0] << "; " << ctrk[1] << "; " << ctrk[2] << "] " << endl;
        }
        if (lay_cross == 1) {
            cout << "apx(x)=" << yx[0] << " apx(y)=" << yy[0] << " apx(z)=" << yz[0] << endl;
            cout << "e_apx(x)=" << sqrt(1/wx[0]) << " e_apx(y)=" << sqrt(1 / wy[0]) << " e_apx(z)=" << sqrt(1 / wz[0]) << endl;
            cout << "Parte dal layer: " << Q - 1 << " che vede " << LayE[Q - 1] << " MeV su " << clu_vec.at(i).e << endl;
            apx[0] = yx[0];
            apx[1] = yy[0];
            apx[2] = yz[0];
        }
        else cout << "Parte dal layer: " << Q - 1 << " che vede " << LayE[Q - 1] << " MeV su " << clu_vec.at(i).e << endl;
    }
}

std::tuple<double, double, double, double> fit_ls(int lay, double X[lay], double Y[lay], double W[lay]) {
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
    B = (norm * xya - xa * ya) / (norm*x2a-xa*xa);
    A = ya / norm - B * xa / norm;
    dB = 1 / sqrt(lay * det);
    dA = sqrt(x2a / (lay * det));
    return std::make_tuple(A, B, dA, dB);
}

cluster2 Calc_variables(std::vector<dg_cell> cells)
{
    double x_weighted = 0, y_weighted = 0, z_weighted = 0, t_weighted = 0, x2_weighted = 0, y2_weighted = 0, z2_weighted = 0, Etot = 0, E2tot, EvEtot = 0, EA, EAtot=0, EB, EBtot=0, TA=0, TB=0;
    for (int j=0; j < cells.size(); j++) {
        double d = DfromADC(cells[j].tdc1, cells[j].tdc2);
        double d1,d2, d3;
        d1 = 0.5 * cells[j].l + d;
        d2 = 0.5 * cells[j].l - d;
        double cell_E = kloe_simu::EfromADC(cells[j].adc1, cells[j].adc2, d1, d2, cells[j].lay);
        double cell_T = kloe_simu::TfromTDC(cells[j].tdc1, cells[j].tdc2, cells[j].l);
        if (cells[j].mod > 25) {
            d3 = cells[j].y - d;
            y_weighted = y_weighted - (d3 * cell_E);
            y2_weighted = y2_weighted - (d3 * d3 * cell_E);
            x_weighted = x_weighted + (cells[j].x * cell_E);
            x2_weighted = x2_weighted + (cells[j].x * cells[j].x * cell_E);
        }
        else {
            d3 = cells[j].x - d;
            x_weighted = x_weighted + (d3 * cell_E);
            x2_weighted = x2_weighted + (d3* d3 * cell_E);
            y_weighted = y_weighted + (cells[j].y * cell_E);
            y2_weighted = y2_weighted + (cells[j].y * cells[j].y * cell_E);
        }
        t_weighted = t_weighted + cell_T * cell_E;
        //x_weighted = x_weighted + (cells[j].x * cell_E);
        //x2_weighted = x2_weighted + (cells[j].x * cells[j].x * cell_E);
        //y_weighted = y_weighted + (cells[j].y * cell_E);
        //y2_weighted = y2_weighted + (cells[j].y * cells[j].y * cell_E);
        z_weighted = z_weighted + (cells[j].z * cell_E);
        z2_weighted = z2_weighted + (cells[j].z * cells[j].z * cell_E);
        Etot = Etot + cell_E;
        E2tot = E2tot + cell_E * cell_E;
    }
    x_weighted = x_weighted / Etot;
    y_weighted = y_weighted / Etot;
    z_weighted = z_weighted / Etot;
    t_weighted = t_weighted / Etot;
    if (x_weighted > -0.000001 && x_weighted < 0.000001) x_weighted = 0;
    if (y_weighted > -0.000001 && y_weighted < 0.000001) y_weighted = 0;
    if (z_weighted > -0.000001 && z_weighted < 0.000001) z_weighted = 0;
    double dx, dy, dz, ta, tb;
    double neff = -1;
    double dum = neff / (neff - 1);
    if (neff == 1 && cells.size() == 1) {
        dx = 0;
        dy = 0;
        dz = 0;
    }
    else {
        if (x2_weighted - x_weighted * x_weighted < 0) {
            dx = 0;
        }
        else {
            dx = sqrt(dum * (x2_weighted - x_weighted * x_weighted));
        }
        if (y2_weighted - y_weighted * y_weighted < 0) {
            dy = 0;
        }
        else {
            dy = sqrt(dum * (y2_weighted - y_weighted * y_weighted));
        }
        if (z2_weighted - z_weighted * z_weighted < 0) {
            dz = 0;
        }
        else {
            dz = sqrt(dum * (z2_weighted - z_weighted * z_weighted));
        }
    }
    cluster2 clust;
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

bool RepetitionCheck(std::vector<int> v, int check) {
    if (std::find(v.begin(), v.end(), check) != v.end()) {
        return true;
    }
    else {
        return false;
    }
}

std::pair<std::vector<int>, std::vector<int>> FindNeighbours(int digits[], std::vector<int> already_checked, int size, int check, std::vector<int> vec_cell) {
    for (int i = 0; i < size; i++) {
        int id_check = digits[i];
        if (id_check == check) {
            continue;
        }
        // Middle of the module
        else if (id_check == check + 101 || id_check == check + 100 || id_check == check + 99 || id_check == check + 1 || id_check == check - 1 || id_check == check - 99 || id_check == check - 100 || id_check == check - 101) {
            bool RepCheck = RepetitionCheck(already_checked, id_check);
            if (RepCheck == true) {
                continue;
            }
            else {
                vec_cell.push_back(i);
                already_checked.push_back(id_check);
                std::pair < vector<int>, vector<int>> check_cell = FindNeighbours(digits, already_checked, size, id_check, vec_cell);
                already_checked = check_cell.first;
                vec_cell = check_cell.second;
            }
        } 
        // Right edge of the module
        else if (id_check == check - 889 || id_check == check - 989 || id_check == check - 1089) {
            bool RepCheck = RepetitionCheck(already_checked, id_check);
            if (RepCheck == true) {
                continue;
            }
            else {
                vec_cell.push_back(i);
                already_checked.push_back(id_check);
                std::pair < vector<int>, vector<int>> check_cell = FindNeighbours(digits, already_checked, size, id_check, vec_cell);
                already_checked = check_cell.first;
                vec_cell = check_cell.second;
            }
        }
        // Right edge of 0 module
        else if (id_check == check + 23111 || id_check == check + 23011 || id_check == check + 22911) {
            bool RepCheck = RepetitionCheck(already_checked, id_check);
            if (RepCheck == true) {
                continue;
            }
            else {
                vec_cell.push_back(i);
                already_checked.push_back(id_check);
                std::pair < vector<int>, vector<int>> check_cell = FindNeighbours(digits, already_checked, size, id_check, vec_cell);
                already_checked = check_cell.first;
                vec_cell = check_cell.second;
            }
        }
        // Left edge of the module
        else if (id_check == check + 1089 || id_check == check + 989 || id_check == check + 889) {
            bool RepCheck = RepetitionCheck(already_checked, id_check);
            if (RepCheck == true) {
                continue;
            }
            else {
                vec_cell.push_back(i);
                already_checked.push_back(id_check);
                std::pair < vector<int>, vector<int>> check_cell = FindNeighbours(digits, already_checked, size, id_check, vec_cell);
                already_checked = check_cell.first;
                vec_cell = check_cell.second;
            }
        }
        // Left edge of the 23 module
        else if (id_check == check - 23011 || id_check == check - 22911 || id_check == check - 23111) {
            bool RepCheck = RepetitionCheck(already_checked, id_check);
            if (RepCheck == true) {
                continue;
            }
            else {
                vec_cell.push_back(i);
                already_checked.push_back(id_check);
                std::pair < vector<int>, vector<int>> check_cell = FindNeighbours(digits, already_checked, size, id_check, vec_cell);
                already_checked = check_cell.first;
                vec_cell = check_cell.second;
            }
        }
    }
    return std::make_pair(already_checked, vec_cell);
}

std::pair<int, int> Next(std::vector<int> already_checked, int digits[], int size) {
    int next = 0;
    int entry = 0;
    for (int i = 0; i < size; i++) {
        bool RepCheck = RepetitionCheck(already_checked, digits[i]);
        if (RepCheck == true) {
            continue;
        }
        else {
            next = digits[i];
            entry = i;
            //cout << i << endl;
            break;
        }
    }
    return std::make_pair(next, entry);
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
    case 1:
        atl2 = kloe_simu::atl2_01;
        break;

    case 2:
        atl2 = kloe_simu::atl2_2;
        break;

    case 3:
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

    return 0.5 * (adc1 / f1 + adc2 / f2) * kloe_simu::adc2MeV;
}

double DfromADC(double ta, double tb) {
    return 0.5 * (ta - tb) / kloe_simu::vlfb * kloe_simu::m_to_mm;
}