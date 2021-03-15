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

void Merge(std::vector<cluster2>, std::vector<cluster2>);

double TfromTDC(double t1, double t2, double L);
double AttenuationFactor(double d, int planeID);
double EfromADC(double adc1, double adc2, double d1, double d2, int planeID);
double DfromADC(double , double);

//Preclustering (Clustering senza merging)


/*
struct cluster2
{
    double x;
    double y;
    double z;
    double t;
    double e;
    double sx;
    double sy;
    double sz;
    double varx;
    double vary;
    double varz;
    double vart;
    std::vector<dg_cells> cells;
};

struct dg_cell
{
  int id;
  double z;
  double y;
  double x;
  double l;
  double adc1;
  double tdc1;
  double adc2;
  double tdc2;
  int mod;
  int lay;
  int cel;
  std::vector<double> pe_time1;
  std::vector<int> hindex1;
  std::vector<double> pe_time2;
  std::vector<int> hindex2;
};

*/

void Preclustering(std::vector<dg_cell> *vec_cell_tree, std::vector<cluster2>& vec_clust)
{
	Double_t x_weighted = 0, y_weighted = 0, z_weighted = 0, t_weighted = 0, x2_weighted = 0, y2_weighted = 0, z2_weighted = 0, t2_weighted = 0, Etot = 0, E2tot, EvEtot = 0, d1, d2;
    std::vector<int> cell_id, cell_lay;
	std::vector<double> cell_x, cell_y, cell_z, cell_l, cell_adc1, cell_tdc1, cell_adc2, cell_tdc2;
    std::vector <dg_cell> evt_cells;
	for(unsigned int i; i < vec_cell_tree->size(); i++){
        evt_cells.push_back(vec_cell_tree->at(i));
        cell_id.push_back(vec_cell_tree->at(i).id);
        cell_x.push_back(vec_cell_tree->at(i).x);
        cell_y.push_back(vec_cell_tree->at(i).y);
        cell_z.push_back(vec_cell_tree->at(i).z);
        cell_l.push_back(vec_cell_tree->at(i).l);
        cell_adc1.push_back(vec_cell_tree->at(i).adc1);
        cell_tdc1.push_back(vec_cell_tree->at(i).tdc1);
        cell_adc2.push_back(vec_cell_tree->at(i).adc2);
        cell_tdc2.push_back(vec_cell_tree->at(i).tdc2);
        cell_lay.push_back(vec_cell_tree->at(i).lay);
    }
    int lentry = cell_id.size();
    int digits[lentry];
    std::vector<int> checked_array;
    int n_cluster = 0;
    int cluster[lentry];
    // Loop over event digits
    for (int j = 0; j < lentry; j++) {
        digits[j] = cell_id.at(j);
    }
    std::vector<int> vec_cell;
    for (int k = 0; k < lentry; k++) {
        cluster2 clust;
        if (k == 0) {
            cluster[k] = digits[0];
        }   else {
            std::pair<int, int> clu_nclu = Next(checked_array, digits, lentry);
            cluster[k] = clu_nclu.first;
            if (cluster[k] == 0) {
                break;
            }
            n_cluster = clu_nclu.second;
        }
        vec_cell.push_back(n_cluster);
        checked_array.push_back(cluster[k]);
        std::pair < vector<int>, vector<int>> check_cell = FindNeighbours(digits, checked_array, lentry, cluster[k], vec_cell);
        checked_array = check_cell.first;
        vec_cell = check_cell.second;
        std::vector<int> vec_cellid;
        for (std::vector<int>::const_iterator it = vec_cell.begin(); it != vec_cell.end(); it++) {   
            vec_cellid.push_back(*it);
            if (cell_tdc1.at(*it) == 0 || cell_tdc2.at(*it) == 0) {
                continue;
            }
            else {
                d1 = DfromADC(cell_tdc1.at(*it), cell_tdc2.at(*it));
                if (cell_tdc1.at(*it) <= cell_tdc2.at(*it)) {
                    d2 = 0.5 * (4300 - d1);
                    d1 = 0.5 * (4300 + d1);
                }
                else {
                    d2 = 0.5 * (4300 - d1);
                    d1 = 0.5 * (4300 + d1);
                }
                double cell_E = kloe_simu::EfromADC(cell_adc1.at(*it), cell_adc2.at(*it), d1, d2, cell_lay.at(*it));
                double cell_T = kloe_simu::TfromTDC(cell_tdc1.at(*it), cell_tdc2.at(*it), cell_l.at(*it));
                t_weighted = t_weighted + cell_T * cell_E;
                t2_weighted = t2_weighted + cell_T * cell_T * cell_E;
                x_weighted = x_weighted + (cell_x.at(*it) * cell_E);
                x2_weighted = x2_weighted + (cell_x.at(*it) * cell_x.at(*it) * cell_E);
                y_weighted = y_weighted + (cell_y.at(*it) * cell_E);
                y2_weighted = y2_weighted + (cell_y.at(*it) * cell_y.at(*it) * cell_E);
                z_weighted = z_weighted + (cell_z.at(*it) * cell_E);
                z2_weighted = z2_weighted + (cell_z.at(*it) * cell_z.at(*it) * cell_E);
                Etot = Etot + abs(cell_E);
                E2tot = E2tot + cell_E * cell_E;
            }
        } 
        double dx, dy, dz, dt;
        double neff = -1;
        double dum = neff / (neff - 1);
        if (Etot != 0 && E2tot != 0) {
            neff = Etot * Etot / E2tot;
            t_weighted = t_weighted / Etot;
            t2_weighted = t2_weighted / Etot;
            x_weighted = x_weighted / Etot;
            x2_weighted = x2_weighted / Etot;
            y_weighted = y_weighted / Etot;
            y2_weighted = y2_weighted / Etot;
            z_weighted = z_weighted / Etot;
            z2_weighted = z2_weighted / Etot;
            if (neff == 1 && vec_cellid.size()==1) {
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
                //cout << "dum: " << dum << " x_w: " << x_weighted << " x2_w: " << x2_weighted << " Diff: " << x2_weighted - x_weighted * x_weighted << endl;
            /*}if (dum * (x2_weighted - x_weighted * x_weighted) >= 0) {
                dx = sqrt(dum * (x2_weighted - x_weighted * x_weighted));
                cout << "neff: " << neff << " x_w: " << x_weighted << " x2_w: " << x2_weighted << " Diff: "<< x2_weighted - x_weighted * x_weighted <<endl;
            }
            if (dum *(y2_weighted - y_weighted * y_weighted) >= 0) {
                dy = sqrt(dum * (y2_weighted - y_weighted * y_weighted));
                cout << "neff: " << neff << " y_w: " << y_weighted << " y2_w: " << y2_weighted << " Diff: " << y2_weighted - y_weighted * y_weighted << endl;
            }
            if (dum * (z2_weighted - z_weighted * z_weighted) >= 0) {
                dz = sqrt(neff / (neff - 1) * (z2_weighted - z_weighted * z_weighted));
                cout << "neff: " << neff << " z_w: " << z_weighted << " z2_w: " << z2_weighted << " Diff: " << z2_weighted - z_weighted * z_weighted << endl;
            }
            if (dum * (t2_weighted - t_weighted * t_weighted) >= 0) {
                dt = sqrt(neff / (neff - 1) * (t2_weighted - t_weighted * t_weighted));
                cout << "neff: " << neff << " t_w: " << t_weighted << " t2_w: " << t2_weighted << " Diff: " << t2_weighted - t_weighted * t_weighted << endl; */
            }
        }
        // sx - sy - sz sarebbe meglio dopo ever fatto il mergin dei cluster che appertengono alla stessa shower
        std::vector<dg_cell> cluster_cells;
        for (int its = 0; its < vec_cellid.size(); its++) {
            int cell = vec_cellid[its];
            cluster_cells.push_back(vec_cell_tree->at(cell));
        }
        vec_cellid.clear();
        clust.e = Etot;
        clust.x = x_weighted;
        clust.y = y_weighted;
        clust.z = z_weighted;
        clust.t = t_weighted;
        clust.varx = dx;
        clust.vary = dy;
        clust.varz = dz;
        clust.vart = dt;
        clust.cells = cluster_cells;
        //cout << "dx: " << clust.varx << " dy: " << clust.vary << " dz: " << clust.varz << " dt: " << clust.vart << " E: " << clust.e<<" [";
        //for (int k = 0; k < clust.cells.size(); k++) {
        //    cout<<clust.cells.at(k).id<< " ";
        //}
        //cout <<"]"<< endl;
        vec_clust.push_back(clust);
        vec_cell.clear();
        x_weighted = 0;
        y_weighted = 0;
        z_weighted = 0;
        t_weighted = 0;
        x2_weighted = 0;
        y2_weighted = 0;
        z2_weighted = 0;
        t2_weighted = 0;
        Etot = 0;
        E2tot = 0;
        dx = 0;
        dy = 0;
        dz = 0;
        dt = 0;
    }
    // Una qualche funzione che prende vec_clust e fa il merging
    std::vector<cluster2> mgd_clust;
    Merge(vec_clust, mgd_clust);

    // Una qualche funzione che prende il vec_clust merged e sputa fuori sx sy sz
//    cout << "////" << endl;
}

void Merge(std::vector<cluster2> Og_cluster, std::vector<cluster2> mgd_clust) {
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
            if (clust.e == 0) {
                continue;
            }
            double xj = Og_cluster.at(j).x;
            double yj = Og_cluster.at(j).y;
            double zj = Og_cluster.at(j).z;
            double tj = Og_cluster.at(j).t;
            double ej = Og_cluster.at(j).e;
            double D = sqrt((xi - xj)*  (xi - xj) + (yi - yj)* (yi - yj) + (zi - zj)* (zi - zj));
            double DT = abs(ti - tj);
            if (D < 40 && DT < 2.5) {
                bool endcap = false;
                if (clust.cells[0].id > 25000) {
                    endcap = true;
                }
                /*if (ti < tj) {
                    cout << xi << " " << yi << " " << zi << " vs " << xj << " " << yj << " " << zj << endl;
                    cout << "Unire il cluster i ad j .1" << endl;
                    std::vector<dg_cell> vec_cells_j = Og_cluster.at(j).cells;
                    //loop su vec_cells_j e push_back su clust_cells
                    for (int k = 0; k < vec_cells_j.size(); k++) {
                        clust.cells.push_back(vec_cells_j.at(k));
                    }
                    //se questo loop funziona, fare loop sulle celle di clust.cells per costruire il nuovo cluster 
                    clust.x = ((xi * ei) + (xj * ej)) / (ei + ej);
                    clust.y = ((yi * ei) + (yj * ej)) / (ei + ej);
                    clust.z = ((zi * ei) + (zj * ej)) / (ei + ej);
                    clust.t = ((ti * ti) + (tj * tj)) / (ei + ej);
                    clust.e = ei + ej;
                    //funzione per calcolare la nuova varianza

                    checked.push_back(j);
                }*/
                //trovare i valori di destro e sinistro (i limiti dell'endcap)
                if (endcap == true) {
                    //siamo sull'endcap
                    double Dz_ec = abs(yi - yj);
                    D = sqrt((xi - xj)* (xi - xj) + (zi - zj)* (zi - zj));
                    if (Dz_ec < 30 && D < 40) {
                        //Unisce il cluster j ad i
                        cout << "Unire il cluster j ad i .2 endcap" << endl;
                        std::vector<dg_cell> vec_cells_j = Og_cluster.at(j).cells;
                        //loop su vec_cells_j e push_back su clust_cells
                        for (int k = 0; k < vec_cells_j.size(); k++) {
                            clust.cells.push_back(vec_cells_j.at(k));
                        }
                        for (int l = 0; l < clust.cells.size(); l++) {
                            cout << clust.cells[l].id << endl;
                        }
                        clust.x = ((xi * ei) + (xj * ej)) / (ei + ej);
                        clust.y = ((yi * ei) + (yj * ej)) / (ei + ej);
                        clust.z = ((zi * ei) + (zj * ej)) / (ei + ej);
                        clust.t = ((ti * ti) + (tj * tj)) / (ei + ej);
                        clust.e = ei + ej;
                        //funzione per calcolare la nuova varianza

                        checked.push_back(j);
                    }
                }
                else if (endcap == false) {
                    double Dz_bar = abs(zi - zj);
                    D = sqrt((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj));
                    if (Dz_bar< 30 && D < 40) {
                            //Unisce il cluster j ad i
                            cout << "Unire il cluster j ad i .3 barrel" << endl;
                            std::vector<dg_cell> vec_cells_j = Og_cluster.at(j).cells;
                            //loop su vec_cells_j e push_back su clust_cells
                            for (int k = 0; k < vec_cells_j.size(); k++) {
                                clust.cells.push_back(vec_cells_j.at(k));
                            }
                            double e_m = 0, e2_m=0, t_m = 0, t2_m = 0, x_m = 0, x2_m = 0, y_m = 0, y2_m = 0, z_m = 0, z2_m = 0, dx_m, dy_m, dz_m, dt_m, neff;
                            for (int l = 0; l < clust.cells.size(); l++) {
                                cout << clust.cells[l].id << endl;
                                if (clust.cells[l].tdc1 == 0 || clust.cells[l].tdc2 == 0) {
                                    continue;
                                }
                                else {
                                    double d1_m = DfromADC(clust.cells[l].tdc1, clust.cells[l].tdc2);
                                    double d2_m;
                                    if (clust.cells[l].tdc1 <= clust.cells[l].tdc2) {
                                        d2_m = 0.5 * (4300 - d1_m);
                                        d1_m = 0.5 * (4300 + d1_m);
                                    }
                                    else {
                                        d2_m = 0.5 * (4300 - d1_m);
                                        d1_m = 0.5 * (4300 + d1_m);
                                    }
                                    double cell_E = kloe_simu::EfromADC(clust.cells[l].adc1, clust.cells[l].adc2, d1_m, d2_m, clust.cells[l].lay);
                                    double cell_T = kloe_simu::TfromTDC(clust.cells[l].tdc1, clust.cells[l].tdc2, clust.cells[l].l);
                                    if (cell_E != 0) {
                                        t_m = t_m + cell_T * cell_E;
                                        t2_m = t2_m + cell_T * cell_T * cell_E;
                                        x_m = x_m + (clust.cells[l].x * cell_E);
                                        x2_m = x2_m + (clust.cells[l].x * clust.cells[l].x * cell_E);
                                        y_m = y_m + (clust.cells[l].y * cell_E);
                                        y2_m = y2_m + (clust.cells[l].y * clust.cells[l].y * cell_E);
                                        z_m = z_m + (clust.cells[l].z * cell_E);
                                        z2_m = z2_m + (clust.cells[l].z * clust.cells[l].z * cell_E);
                                        e_m = e_m + cell_E;
                                        e2_m = e2_m + cell_E * cell_E;
                                    }
                                }
                            }
                            if (e2_m != 0) {
                                neff = e_m * e_m / e2_m;
                                if (neff / (neff - 1) * (x2_m - x_m * x_m) >= 0) {
                                    dx_m = sqrt(neff / (neff - 1) * (x2_m - x_m * x_m));
                                }
                                if (neff / (neff - 1) * (y2_m - y_m * y_m) >= 0) {
                                    dy_m = sqrt(neff / (neff - 1) * (y2_m - y_m * y_m));
                                }
                                if (neff / (neff - 1) * (z2_m - z_m * z_m) >= 0) {
                                    dz_m = sqrt(neff / (neff - 1) * (z2_m - z_m * z_m));
                                }
                                if (neff / (neff - 1) * (t2_m - t_m * t_m) >= 0) {
                                    dt_m = sqrt(neff / (neff - 1) * (t2_m - t_m * t_m));
                                }
                            }
                            clust.x = x_m/e_m;
                            clust.y = y_m/e_m;
                            clust.z = z_m/e_m;
                            clust.t = t_m/e_m;
                            clust.e = e_m;
                            clust.varx = dx_m;
                            clust.vary = dy_m;
                            clust.varz = dz_m;
                            clust.vart = dt_m;
                            checked.push_back(j);
                    }
                }

            }
        }
//        cout << "dx: " << clust.varx << " dy: " << clust.vary << " dz: " << clust.varz << " dt: " << clust.vart << " E: " << clust.e << endl;
        mgd_clust.push_back(clust);
    }
    //return mgd_clust;
}

std::pair<std::vector<int>, std::vector<int>> FindNeighbours(int digits[], std::vector<int> already_checked, int size, int check, std::vector<int> vec_cell) {
    for (int i = 0; i < size; i++) {
        int id_check = digits[i];
        if (id_check == check) {
            continue;
        }
        else if (id_check == check + 101 || id_check == check + 100 || id_check == check + 99 || id_check == check + 1 || id_check == check - 1 || id_check == check - 99 || id_check == check - 100 || id_check == check - 101) {
            bool RepCheck = RepetitionCheck(already_checked, id_check);
            if (RepCheck == true) {
                continue;
            } 
            else {
//                cout << id_check <<"  ("<<i<<") ";
                vec_cell.push_back(i);
                //cout <<i<<" ";
                already_checked.push_back(id_check);
                std::pair < vector<int>, vector<int>> check_cell = FindNeighbours(digits, already_checked, size, id_check, vec_cell);
                already_checked = check_cell.first;
                vec_cell = check_cell.second;
            }
        }
        else if (id_check == check - 889 || id_check == check - 989 || id_check == check - 1089) {
            bool RepCheck = RepetitionCheck(already_checked, id_check);
            if (RepCheck == true) {
                continue;
            } 
            else {
 //               cout << id_check << " ("<<i<<") ";
                vec_cell.push_back(i);
                //cout << i << " ";
                already_checked.push_back(id_check);
                std::pair < vector<int>, vector<int>> check_cell = FindNeighbours(digits, already_checked, size, id_check, vec_cell);
                already_checked = check_cell.first;
                vec_cell = check_cell.second;
            }
        }
        else if ( id_check == check + 23111 ||id_check == check + 23011 || id_check == check + 22911) {
            bool RepCheck = RepetitionCheck(already_checked, id_check);
            if (RepCheck == true) {
                continue;
            }
            else {
//                cout << id_check << " ("<<i<<") ";
                vec_cell.push_back(i);
//                out << i << " ";
                already_checked.push_back(id_check);
                std::pair < vector<int>, vector<int>> check_cell = FindNeighbours(digits, already_checked, size, id_check, vec_cell);
                already_checked = check_cell.first;
                vec_cell = check_cell.second;
            }
        }
        else if (id_check == check + 1089 || id_check == check + 989 || id_check == check + 889) {
            bool RepCheck = RepetitionCheck(already_checked, id_check);
            if (RepCheck == true) {
                continue;
            }
            else {
//                cout << id_check << " (" << i << ") ";
                vec_cell.push_back(i);
                //cout << i << " ";
                already_checked.push_back(id_check);
                std::pair < vector<int>, vector<int>> check_cell = FindNeighbours(digits, already_checked, size, id_check, vec_cell);
                already_checked = check_cell.first;
                vec_cell = check_cell.second;
            }
        }
        else if (id_check == check -23011 || id_check == check -22911 || id_check == check -23111) {
            bool RepCheck = RepetitionCheck(already_checked, id_check);
            if (RepCheck == true) {
                continue;
            }
            else {
//                cout << id_check << " ("<<i<<") ";
                vec_cell.push_back(i);
                //cout << i << " ";
                already_checked.push_back(id_check);
                std::pair < vector<int>, vector<int>> check_cell = FindNeighbours(digits, already_checked, size, id_check, vec_cell);
                already_checked = check_cell.first;
                vec_cell = check_cell.second;
            }
        }
    }
    return std::make_pair(already_checked, vec_cell);
}

bool RepetitionCheck(std::vector<int> v, int check) {
    if (std::find(v.begin(), v.end(), check) != v.end()) {
        return true;
    }
    else {
        return false;
    }
}

std::pair<int,int> Next(std::vector<int> already_checked, int digits[], int size) {
    int next=0;
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

int Testing(std::string input) {
    const char* finname = input.c_str();
    TFile f(finname, "READ");
    TTree* t = (TTree*)f.Get("tDigit");
    int TotHits = t->GetEntries();
    std::vector<dg_cell>* cell=new std::vector<dg_cell>;
    std::vector<cluster2> clust;
    t->SetBranchAddress("dg_cell", &cell);
    for (int i = 0; i <TotHits; i++) {
        clust.clear();
        cout << "Entry " << i << endl;
        t->GetEntry(i);
        Preclustering(cell, clust);
    }
    return 0;
}

/*
int main(){
    Preclustering("../Output_Test2.reco.root");
    return 0;
}
*/
