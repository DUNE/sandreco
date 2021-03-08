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

using namespace std;

bool RepetitionCheck(std::vector<int>, int);

std::pair<vector<int>, vector<int>>FindNeighbours(int[], std::vector<int>, int, int, std::vector<int>);

std::pair<int, int> Next(std::vector<int>, int[], int);

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
    std::vector<int> cells;
};

*/

void Preclustering(std::vector<dg_cell>* vec_cell_tree, std::vector<cluster2>& vec_clust)
{
   /*
    const char* finname = input.c_str();
    TFile f(finname, "READ");


    TTree* t = (TTree*)f.Get("tDigit");
    Int_t TotHits = t->GetEntries();
    // Access Calorimeter leaves
    TLeaf* cell_id = (TLeaf*)t->GetLeaf("dg_cell.id");

    TLeaf* cell_x = (TLeaf*)t->GetLeaf("dg_cell.x");
    TLeaf* cell_y = (TLeaf*)t->GetLeaf("dg_cell.y");
    TLeaf* cell_z = (TLeaf*)t->GetLeaf("dg_cell.z");
    TLeaf* cell_l = (TLeaf*)t->GetLeaf("dg_cell.l");
    TLeaf* cell_lay = (TLeaf*)t->GetLeaf("dg_cell.lay");
    TLeaf* cell_adc1 = (TLeaf*)t->GetLeaf("dg_cell.adc1");
    TLeaf* cell_tdc1 = (TLeaf*)t->GetLeaf("dg_cell.tdc1");
    TLeaf* cell_adc2 = (TLeaf*)t->GetLeaf("dg_cell.adc2");
    TLeaf* cell_tdc2 = (TLeaf*)t->GetLeaf("dg_cell.tdc2");
   
   
    //loop over events
    for (int i = 0; i <5; i++) {
        std::vector<cluster2> vec_clust; //This must be changed in reconstruct.cpp! 
        int pc = 0;
        t->GetEntry(i);

*/
	Double_t x_weighted = 0, y_weighted = 0, z_weighted = 0, t_weighted = 0, x2_weighted = 0, y2_weighted = 0, z2_weighted = 0, t2_weighted = 0, Etot = 0, E2tot, EvEtot = 0, d1, d2;
	std::vector<int> *cell_id, *cell_lay;
	std::vector<double> *cell_x, *cell_y, *cell_z, *cell_l, *cell_adc1, *cell_tdc1, *cell_adc2, *cell_tdc2;
	for(unsigned int j; j < vec_cell_tree->size(); j++){
            cell_id->push_back(vec_cell_tree->at(j).id);
            cell_x->push_back(vec_cell_tree->at(j).x);
            cell_y->push_back(vec_cell_tree->at(j).y);
            cell_z->push_back(vec_cell_tree->at(j).z);
            cell_l->push_back(vec_cell_tree->at(j).l);
            cell_adc1->push_back(vec_cell_tree->at(j).adc1);
            cell_tdc1->push_back(vec_cell_tree->at(j).tdc1);
            cell_adc2->push_back(vec_cell_tree->at(j).adc2);
            cell_tdc2->push_back(vec_cell_tree->at(j).tdc2);
            cell_lay->push_back(vec_cell_tree->at(j).lay);
        }
	int pc = 0;
        int lentry = cell_id->size();
        int digits[lentry];
        std::vector<int> checked_array;
        int n_cluster = 0;
        int cluster[lentry];
        // Loop over event digits
        for (int j = 0; j < lentry; j++) {
            digits[j] = cell_id->at(j);
            //cout << digits[j] << " ";
        }
        //cout << endl;
//        out<<"999"<<endl;
        std::vector<int> vec_cell;
        for (int k = 0; k < lentry; k++) {
            cluster2 clust;
            //std::vector<int> vec_cell;
            if (k == 0) {
                cluster[k] = digits[0];
            }
            else {
                std::pair<int, int> clu_nclu = Next(checked_array, digits, lentry);
                cluster[k] = clu_nclu.first;
                if (cluster[k] == 0) {
                    break;
                }
                n_cluster = clu_nclu.second;
            }
           //cout << n_cluster <<" ";
            vec_cell.push_back(n_cluster);
            //cout << n_cluster << endl;
            checked_array.push_back(cluster[k]);
            std::pair < vector<int>, vector<int>> check_cell = FindNeighbours(digits, checked_array, lentry, cluster[k], vec_cell);
            checked_array = check_cell.first;
            vec_cell = check_cell.second;
            //vec_clust.at(pc).cells.push_back(vec_cell);
            std::vector<int> vec_cellid;
            for (std::vector<int>::const_iterator it = vec_cell.begin(); it != vec_cell.end(); it++) {
                //std::cout << *it << ' ';    
                vec_cellid.push_back(cell_id->at(*it));
                if (cell_tdc1->at(*it) == 0 || cell_tdc2->at(*it) == 0) {
                    continue;
                }   
                else {
                    d1 = DfromADC(cell_tdc1->at(*it), cell_tdc2->at(*it));
                    //cout << d1 << " ";
                    if (cell_tdc1->at(*it) <= cell_tdc2->at(*it)) {
                        d2 = 0.5 * (4300 - d1);
                        d1 = 0.5 * (4300 + d1);
                    }
                    else {
                        d2 = 0.5 * (4300 - d1);
                        d1 = 0.5 * (4300 + d1);
                    }
                    double cell_E = kloe_simu::EfromADC(cell_adc1->at(*it), cell_adc2->at(*it), d1, d2, cell_lay->at(*it));
                    double cell_T = kloe_simu::TfromTDC(cell_tdc1->at(*it), cell_tdc2->at(*it), cell_l->at(*it));
                    t_weighted = t_weighted + cell_T * cell_E;
                    t2_weighted = t2_weighted + cell_T * cell_T * cell_E;
                    x_weighted = x_weighted + (cell_x->at(*it) * cell_E);
                    x2_weighted = x2_weighted + (cell_x->at(*it) * cell_x->at(*it) * cell_E);
                    y_weighted = y_weighted + (cell_y->at(*it) * cell_E);
                    y2_weighted = y2_weighted + (cell_y->at(*it) * cell_y->at(*it) * cell_E);
                    z_weighted = z_weighted + (cell_z->at(*it) * cell_E);
                    z2_weighted = z2_weighted + (cell_z->at(*it) * cell_z->at(*it) * cell_E);
                    Etot = Etot + cell_E;
                    E2tot = E2tot + cell_E * cell_E;
                }
            } 
            double dx = -1;
            double dy = -1;
            double dz = -1;
            double dt = -1;
            double neff = -1;
            if (E2tot != 0) {
                neff= Etot * Etot / E2tot;
            }
            if (Etot != 0) {
                t_weighted = t_weighted / Etot;
                t2_weighted = t2_weighted / Etot;
                x_weighted = x_weighted / Etot;
                x2_weighted = x2_weighted / Etot;
                y_weighted = y_weighted / Etot;
                y2_weighted = y2_weighted / Etot;
                z_weighted = z_weighted / Etot;
                z2_weighted = z2_weighted / Etot;
           
                if (neff / (neff - 1)*(x2_weighted - x_weighted * x_weighted) >= 0) {
                dx = sqrt(neff / (neff - 1) * (x2_weighted - x_weighted * x_weighted));
                }

                if (neff / (neff - 1) *(y2_weighted - y_weighted * y_weighted) >= 0) {
                    dy = sqrt(neff / (neff - 1) * (y2_weighted - y_weighted * y_weighted));
                }

                if (neff / (neff - 1) *(z2_weighted - z_weighted * z_weighted) >= 0) {
                    dz = sqrt(neff / (neff - 1) * (z2_weighted - z_weighted * z_weighted));
                }

                if (neff / (neff - 1) *(t2_weighted - t_weighted * t_weighted) >= 0) {
                    dt = sqrt(neff / (neff - 1) * (t2_weighted - t_weighted * t_weighted));
                }
            }

            // sx - sy - sz sarebbe meglio dopo ever fatto il mergin dei cluster che appertengono alla stessa shower
            clust.cells = vec_cellid;
            clust.e = Etot;
            clust.x = x_weighted;
            clust.y = y_weighted;
            clust.z = z_weighted;
            clust.t = t_weighted;
            clust.varx = dx;
            clust.vary = dy;
            clust.varz = dz;
            clust.vart = dt;

            //vec_clust.push_back(clust);
            cout << "X: " << clust.x << " Y: " << clust.y << " Z: " << clust.z << " T: " << clust.t << " E: " << clust.e<<" [";
            for (std::vector<int>::const_iterator it = clust.cells.begin(); it != clust.cells.end(); it++) {
                std::cout << *it << ' ';
            }
            cout <<"]"<< endl;
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
            //            out << "888"<<endl;
        }
        // Una qualche funzione che prende vec_clust e fa il merging
        // Una qualche funzione che prende il vec_clust merged e sputa fuori sx sy sz
        cout << "////" << endl;
   
//    out << "777" << endl;
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
        else if (id_check == check + 1089 || id_check == check + 989 || id_check == check + 911) {
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
/*
int main(){
    Preclustering("../Output_Test2.reco.root");
    return 0;
}
*/
