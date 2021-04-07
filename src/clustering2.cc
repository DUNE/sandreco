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
void TrackFit(std::vector<cluster2>);
std::vector<cluster2> Merge(std::vector<cluster2>);
cluster2 Calc_variables(std::vector<dg_cell>);
std::vector<cluster2> Split(std::vector<cluster2>);

double TfromTDC(double t1, double t2, double L);
double AttenuationFactor(double d, int planeID);
double EfromADC(double adc1, double adc2, double d1, double d2, int planeID);
double DfromADC(double, double);

std::vector<cluster2> Preclustering(std::vector<dg_cell> *vec_cellraw) {
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
        //cout << "x: " << clust.x << " y: " << clust.y << " z: " << clust.z << " t: " << clust.t << " E: " << clust.e << " [";
        //for (int k = 0; k < clust.cells.size(); k++) {
        //    cout << clust.cells.at(k).id << " ";
        //}
        //cout << "]" << endl;
        vec_clust.push_back(clust);
        vec_cell.clear();
    }

    //TrackFit(vec_clust);
    return vec_clust;
}

//METTERE LA FUNZIONE MERGE E DEBUGGARLA

int Testing(std::string input) {
    const char* finname = input.c_str();
    TFile f(finname, "READ");
    TTree* t = (TTree*)f.Get("tDigit");
    int TotHits = t->GetEntries();
    std::vector<dg_cell>* cell = new std::vector<dg_cell>;
    std::vector<cluster2> clust;
    t->SetBranchAddress("dg_cell", &cell);
    for (int i = 0; i < 5; i++) {
        cout << "Entry " << i << endl;
        t->GetEntry(i);
        clust=Preclustering(cell);
        clust.clear();
    }
    return 0;
}

void TrackFit(std::vector<cluster2> clu_vec) {
    const double xl[5] = { 4.44, 4.44, 4.44, 4.44, 5.24 };
    for (int i = 0; i < clu_vec.size(); i++) {
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
        double W0[3], W1[3], W2[3], W3[3], W4[3];
        double X[5] = { 0,0,0,0,0 }, D;
        if (clu_vec.at(i).cells[0].id > 25000) {
            isBarrel = false;
        }
        int lay_cross = 0, first_lay = 0;
        if (Lay0.e > 0) {
            lay_cross++;
            if (lay_cross == 1) first_lay = 1;
            W0[0] = 0.6;
            W0[1] = 0.6;
            W0[2] = 0.001 * Lay0.e;
            if (isBarrel == false) {
                W0[1] = W0[2];
                W0[2] = 0.6;
            }
            D = D + xl[0];
            if (first_lay == 1) X[0] = 0.5 * xl[0];
        }
        if (Lay1.e > 0) {
            lay_cross++;
            if (lay_cross == 1) first_lay = 2;
            W1[0] = 0.6;
            W1[1] = 0.6;
            W1[2] = 0.001 * Lay2.e;
            if (isBarrel == false) {
                W1[1] = W1[2];
                W1[2] = 0.6;
            }
            D = D + xl[1];
            if (first_lay == 2) X[0] = 0.5 * xl[1];
            else if (first_lay == 1) X[1] = xl[0] + 0.5 * xl[1];
        }
        if (Lay2.e > 0) {
            lay_cross++;
            if (lay_cross == 1) first_lay = 3;
            W2[0] = 0.6;
            W2[1] = 0.6;
            W2[2] = 0.001 * Lay2.e;
            if (isBarrel == false) {
                W2[1] = W2[2];
                W2[2] = 0.6;
            }
            D = D + xl[2];
            if (first_lay == 3) X[0] = 0.5 * xl[2];
            else if (first_lay == 2) X[1] = xl[0] + 0.5 * xl[2];
            else if (first_lay == 1) X[2] = 2 * xl[0] + 0.5 * xl[3];
        }
        if (Lay3.e > 0) {
            lay_cross++;
            if (lay_cross == 1) first_lay = 4;
            W3[0] = 0.6;
            W3[1] = 0.6;
            W3[2] = 0.001 * Lay3.e;
            if (isBarrel == false) {
                W3[1] = W3[2];
                W3[2] = 0.6;
            }
            D = D + xl[3];
            if (first_lay == 4) X[0] = 0.5 * xl[3];
            else if (first_lay == 3) X[1] = xl[0] + 0.5 * xl[3];
            else if (first_lay == 2) X[2] = 2 * xl[0] + 0.5 * xl[3];
            else if (first_lay == 1) X[3] = 3 * xl[0] + 0.5 * xl[3];
        }
        if (Lay4.e > 0) {
            lay_cross++;
            if (lay_cross == 1) first_lay = 5;
            W4[0] = 0.6;
            W4[1] = 0.6;
            W4[2] = 0.001 * Lay4.e;
            if (isBarrel == false) {
                W4[1] = W4[2];
                W4[2] = 0.6;
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
        int Q = 0;
        for (int k = first_lay; k <= 5; k++) {
            cout << "Layer " << k - 1 << " Energy " << LayE[k - 1] << endl;
            if (Q == 0 && (LayE[k - 1] >= 0.05 * clu_vec.at(i).e)) {
                Q = k;
            }
        }
        double E1 = 0, E2 = 0;
        if (lay_cross > 1) {
            for (int k_i = 5; k_i > Q; k_i--) {
                E2 = E1;
                E1 = E1 + LayE[k_i - 1];
            }
            double Rk = E2 / E1;
            cout << "Parte dal layer: " << first_lay - 1 << " che vede " << LayE[Q - 1] << " MeV su " << clu_vec.at(i).e << " MeV. Energia residua: " << E1 << " MeV" << endl;
            double B = 3;
            if (clu_vec.at(i).e > 16.5) B = 1.5 / log(clu_vec.at(i).e / 10);
            double Zmin = 0;
            double Zapx = 0.5 * B * xl[Q - 1];
            double Zmax = B * xl[Q - 1];
            double R1;
            for (int j = 0; j < 4; j++) {
                R1 = exp(-Zapx) * (1 + Zapx);
                if (R1 > Rk) {
                    Zmin = Zapx;
                    Zapx = 0.5 * (Zmax + Zmin);
                }
                else if (R1 < Rk) {
                    Zmax = Zapx;
                    Zapx = 0.5 * (Zmax + Zmin);
                }
            }
            Zapx = -Zapx / B;
            for (int j = 0; j < Q; j++) {
                Zapx = Zapx + xl[j];
                cout << Zapx << " " << xl[j] << endl;
            }
            for (int j = 0; j < lay_cross; j++) {
                X[j] = X[j] - Zapx;
                cout << "Zapx: " << X[j] << endl;
            }
        }
        else cout << "Parte dal layer: " << Q - 1 << " che vede " << LayE[Q - 1] << " MeV su " << clu_vec.at(i).e << endl;
    }
}

cluster2 Calc_variables(std::vector<dg_cell> cells)
{
    double x_weighted = 0, y_weighted = 0, z_weighted = 0, t_weighted = 0, x2_weighted = 0, y2_weighted = 0, z2_weighted = 0, t2_weighted = 0, Etot = 0, E2tot, EvEtot = 0;
    for (int j=0; j < cells.size(); j++) {
        double d1 = DfromADC(cells[j].tdc1, cells[j].tdc2);
        double d2;
        if (cells[j].tdc1 <= cells[j].tdc2) {
            d1 = 0.5 * (4300 - d1);
            d2 = 0.5 * (4300 + d1);
        }
        else {
            d2 = 0.5 * (4300 - d1);
            d1 = 0.5 * (4300 + d1);
        }
        double cell_E = kloe_simu::EfromADC(cells[j].adc1, cells[j].adc2, d1, d2, cells[j].lay);
        double cell_T = kloe_simu::TfromTDC(cells[j].tdc1, cells[j].tdc2, cells[j].l);
        t_weighted = t_weighted + cell_T * cell_E;
        t2_weighted = t2_weighted + cell_T * cell_T * cell_E;
        x_weighted = x_weighted + (cells[j].x * cell_E);
        x2_weighted = x2_weighted + (cells[j].x * cells[j].x * cell_E);
        y_weighted = y_weighted + (cells[j].y * cell_E);
        y2_weighted = y2_weighted + (cells[j].y * cells[j].y * cell_E);
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

    double dx, dy, dz, dt;
    double neff = -1;
    double dum = neff / (neff - 1);
    if (neff == 1 && cells.size() == 1) {
        dx = 0;
        dy = 0;
        dz = 0;
        dt = 0;
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
        if (t2_weighted - t_weighted * t_weighted < 0) {
            dt = 0;
        }
        else {
            dt = sqrt(dum * (t2_weighted - t_weighted * t_weighted));
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
    clust.vart = dt;
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

std::vector<cluster2> Split(std::vector<cluster2> original_clu_vec){
    std::vector<cluster2> clu_vec;
    std::vector <dg_cell> all_cells;
    for(int k=0; k < original_clu_vec.size(); k++){
        double tA = 0, tB = 0, tA2 = 0, tB2 = 0;
        double EA, EB, EAtot = 0, EBtot = 0, EA2tot = 0, EB2tot = 0;
        double tRMS_A, tRMS_B, dist;
        all_cells = original_clu_vec.at(k).cells;
        for(int j=0; j < all_cells.size(); j++){
            EA = all_cells.at(j).adc1;
            EB = all_cells.at(j).adc2;
            EAtot += EA;
            EA2tot += EA * EA;
            EBtot += EB;
            EB2tot += EB * EB;
            double d1 = DfromADC(all_cells.at(j).tdc1, all_cells.at(j).tdc2);
            double d2;
            if (all_cells.at(j).tdc1 <= all_cells.at(j).tdc2) {
                d1 = 0.5 * (all_cells.at(j).l - d1);
                d2 = 0.5 * (all_cells.at(j).l + d1);
            }
            else {
                d2 = 0.5 * (all_cells.at(j).l - d1);
                d1 = 0.5 * (all_cells.at(j).l + d1);
            }
            tA += (all_cells.at(j).tdc1 - kloe_simu::vlfb *d1/ kloe_simu::m_to_mm) * EA;
            tA2 += std::pow( all_cells.at(j).tdc1 - kloe_simu::vlfb *d1/ kloe_simu::m_to_mm, 2) * EA;
            tB += (all_cells.at(j).tdc2 - kloe_simu::vlfb *d2/ kloe_simu::m_to_mm) * EB;
            tB2 += std::pow(all_cells.at(j).tdc2 - kloe_simu::vlfb *d2/ kloe_simu::m_to_mm, 2) * EB;
        }
        tA = tA/EAtot;
        tA2 = tA2/EAtot;
        tB = tB/EAtot;
        tB2 = tB2/EBtot;
        tRMS_A = (tA2 - tA*tA) * (EA2tot - EAtot*EAtot)/EA2tot;
        tRMS_B = (tB2 - tB*tB) * (EB2tot - EBtot*EBtot)/EB2tot;
        dist = std::sqrt(tRMS_A * tRMS_A + tRMS_B * tRMS_B);
        
        
        if(dist > 5){
            std::vector <dg_cell> q1_cells, q2_cells, q3_cells, q4_cells;        
            for(unsigned int i =0; i < all_cells.size(); i++){     
                double t_difA = all_cells.at(i).tdc1 - tA;
                double t_difB = all_cells.at(i).tdc2 - tB;
                if (t_difA > 0){
                  if(t_difB > 0){q1_cells.push_back(all_cells.at(i));}
                  else{q2_cells.push_back(all_cells.at(i));}
                }
                else{
                  if(t_difB > 0){q3_cells.push_back(all_cells.at(i));}
                  else{q4_cells.push_back(all_cells.at(i));}
                }          
            }
            if(q1_cells.size() != 0){
                cluster2 clus = Calc_variables(q1_cells);
                clu_vec.push_back(clus);
            }
            if(q2_cells.size() != 0){
                cluster2 clus = Calc_variables(q2_cells);
                clu_vec.push_back(clus);
            }
            if(q3_cells.size() != 0){
                cluster2 clus = Calc_variables(q3_cells);
                clu_vec.push_back(clus);
            }
            if(q4_cells.size() != 0){
                cluster2 clus = Calc_variables(q4_cells);
                clu_vec.push_back(clus);
            }
            q1_cells.clear();
            q2_cells.clear();
            q3_cells.clear();
            q4_cells.clear();            
        }
        else{  
           clu_vec.push_back(original_clu_vec.at(k));  
        }
        all_cells.clear();        
    }
    original_clu_vec.clear();
    return clu_vec;
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
