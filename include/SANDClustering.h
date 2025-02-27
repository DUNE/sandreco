#ifndef SANDCLUSTERING_H
#define SANDCLUSTERING_H

#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>

#include "TFile.h"
#include "TTree.h"

#include "struct.h"

#include "SANDGeoManager.h"

std::tuple<double, double, double, double> fitLs(int, double[], double[],
                                                  double[]);

std::pair<std::vector<dg_cell>, std::vector<dg_cell>> processMultiHits(const SANDGeoManager* sand_geo,
    const std::vector<dg_cell>&);
std::pair<std::vector<dg_cell>, std::vector<int>> getNeighbours(
    const std::vector<dg_cell>&, int, std::vector<int>, std::vector<dg_cell>);

std::vector<cluster> clusterize(const SANDGeoManager* sand_geo, const std::vector<dg_cell>&);
std::vector<cluster> merge(const std::vector<cluster>&);
std::vector<cluster> split(const SANDGeoManager* sand_geo, const std::vector<cluster>&, bool&);

void recoverIncomplete(const SANDGeoManager* sand_geo, std::vector<cluster>&, const std::vector<dg_cell>&);
void clustInfo(cluster);
void trackFit(std::vector<cluster>&);

cluster calcVariables(const std::vector<reco_cell>&);
cluster createCluster(const SANDGeoManager* sand_geo, const std::vector<dg_cell>&);

bool repetitionCheck(std::vector<int>, int);
bool isNeighbour(const dg_cell&, const dg_cell&);

bool endsWith(const std::string& fullString, const std::string& ending);


#endif
