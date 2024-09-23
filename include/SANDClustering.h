#ifndef SANDCLUSTERING_H
#define SANDCLUSTERING_H

#include <complex>      // check if needed
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TSystem.h"

//#include "SANDClusteringLinkDef.h"
#include "struct.h"
//#include "utils.h"

std::tuple<double, double, double, double> fit_ls(int, double[], double[],
                                                  double[]);

std::pair<std::vector<dg_cell>, std::vector<dg_cell>> ProcessMultiHits(
    std::vector<dg_cell>, std::vector<dg_cell>);
std::pair<std::vector<dg_cell>, std::vector<int>> GetNeighbours(
    std::vector<dg_cell>, int, std::vector<int>, std::vector<dg_cell>);


std::vector<cluster> Clusterize(std::vector<dg_cell>*);
void Clust_info(cluster);
std::vector<cluster> TrackFit(std::vector<cluster>);
std::vector<cluster> Merge(std::vector<cluster>);
std::vector<cluster> Split(std::vector<cluster>, bool&);
std::vector<cluster> RecoverIncomplete(std::vector<cluster>,
                                       std::vector<dg_cell>);
cluster Calc_variables(std::vector<reco_cell>);
cluster Create_cluster(std::vector<dg_cell>);

bool RepetitionCheck(std::vector<int>, int);
bool isNeighbour(int, int);

//
/// Maybe duplication of identical functions
//
double TfromTDC(double t1, double t2, double L);  
double AttenuationFactor(double d, int planeID);
double EfromADC(double adc1, double adc2, double d1, double d2, int planeID);
double EfromADCsingle(double adc, double f);
double DfromTDC(double, double);

bool endsWith(const std::string& fullString, const std::string& ending);

// void Clust_info(cluster);

#endif
