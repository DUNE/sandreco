#pragma once

#include <TArrow.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TGraph.h>
#include <TMarker.h>
#include <TMultiGraph.h>
#include <TTree.h>
#include <TVector3.h>
#include <TVectorD.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <unordered_map>
#include <vector>

#include "utils.h"

#include "EDEPTree.h"
#include "SANDTrackerUtils.h"

// ------------------------------------------------------------------------------
// New type introduced as the truth at each z of the cluster, note that "truth"
// does not referr to a unique trajectory point from MC info but to the computed
// tracklet (pos, dir, mom) at that cluster obtained from the interpolation
// between start and stop position of the HitSegment.
// The corresponding z-position of the truth is used to map the tracklets along
// the trajectory.
// ------------------------------------------------------------------------------
struct Truth {
  TVector3 pos_;
  TVector3 dir_;
  TVector3 mom_;
};

double ComputeStd(const std::vector<double>& values, double mean);

Truth getTrueTrackletOfCluster(TVector3 start_pos, TVector3 stop_pos,
                               TVector3 start_mom, TVector3 stop_mom, double z);

Truth getTrueTrackletFromTrajectoryPoint(
    const std::vector<EDEPTrajectoryPoint>& points, double z);

std::map<double, Truth> z_to_truth(
    const std::vector<EDEPTrajectoryPoint>& points, double step);

double getScore(const TVectorD& tracklet,
                const std::vector<TVector3>& true_tracklet);
std::map<double, std::vector<TVector3>> getInterpolatedZ(
    const std::vector<EDEPTrajectoryPoint>& trj_points,
    const sand_reco::kf::utils::TrackletMap& z_to_tracklets);

std::vector<double> computeZDistance(
    const std::vector<EDEPTrajectoryPoint>& trj_points,
    const sand_reco::kf::utils::TrackletMap& z_to_tracklets);

sand_reco::kf::utils::TrackletMap findBestTracklet(
    const sand_reco::kf::utils::TrackletMap& z_to_tracklets,
    const std::map<double, std::vector<TVector3>>& z_to_interpolated_tracklets);
