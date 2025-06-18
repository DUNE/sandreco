#pragma once

#include <TVector3.h>
#include <TVectorD.h>
#include <TTree.h>
#include <TFile.h>
#include <TMarker.h>
#include <TArrow.h>
#include <TDatabasePDG.h>
#include <TGraph.h>
#include <TMultiGraph.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <random>
#include <limits>

#include "utils.h"

#include "EDEPTree.h"

double ComputeStd(const std::vector<double>& values, double mean);

std::map<double, std::vector<TVector3>> getInterpolatedZ(const std::vector<EDEPTrajectoryPoint>& trj_points,
                                                         const std::map<double, std::vector<TVectorD>>& z_to_tracklets);

std::vector<double> computeZDistance(const std::vector<EDEPTrajectoryPoint>& trj_points,
                                                         const std::map<double, std::vector<TVectorD>>& z_to_tracklets);

std::map<double, std::vector<TVectorD>> findBestTracklet(const std::map<double, std::vector<TVectorD>>& z_to_tracklets,
                                                         const std::map<double, std::vector<TVector3>>& z_to_interpolated_tracklets);
                                                         

