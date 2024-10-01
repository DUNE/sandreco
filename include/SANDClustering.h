#ifndef SANDCLUSTERING_H
#define SANDCLUSTERING_H

#include <fstream>
#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"

#include "struct.h"

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

double EfromADCsingle(double adc, double f);
double DfromTDC(double, double);

bool endsWith(const std::string& fullString, const std::string& ending);


#endif
