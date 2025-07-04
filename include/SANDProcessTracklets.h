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
#include "SANDTrackerUtils.h"

double ComputeStd(const std::vector<double>& values, double mean);

std::map<double, std::vector<TVector3>> getInterpolatedZ(const std::vector<EDEPTrajectoryPoint>& trj_points,
                                                         const sand_reco::kf::utils::TrackletMap& z_to_tracklets);

std::vector<double> computeZDistance(const std::vector<EDEPTrajectoryPoint>& trj_points,
                                                         const sand_reco::kf::utils::TrackletMap& z_to_tracklets);

sand_reco::kf::utils::TrackletMap findBestTracklet(const sand_reco::kf::utils::TrackletMap& z_to_tracklets,
                                                         const std::map<double, std::vector<TVector3>>& z_to_interpolated_tracklets);
std::vector<TVector3> getTrueTrackletOfCluster(TVector3 start, TVector3 stop, double z);          
double getScore(const TVectorD& tracklet, const std::vector<TVector3>& true_tracklet);
