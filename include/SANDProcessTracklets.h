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

// ------------------------------------------------------------------------------------------------------------------------
// New type introduced as the truth at each z of the cluster, note that "truth" does not referr to a unique 
// trajectory point from MC info but to the computed tracklet (pos, dir, mom) at that cluster obtained from the interpolation
// betweena start and stop position of the HitSegment.
// The corresponding z-position of the truth is used to map the tracklets along the trajectory.
// -------------------------------------------------------------------------------------------------------------------------
struct Truth {
    TVector3 pos_; 
    TVector3 dir_;
    TVector3 mom_;
};

double ComputeStd(const std::vector<double>& values, double mean);

//Truth of each clsuter at z, not sure about this.. it should be a segment as TG4HitSegment or just a tracklet with pos and dir
Truth getTrueTrackletOfCluster(
    TVector3 start_pos, TVector3 stop_pos,
    TVector3 start_mom, TVector3 stop_mom,
    double z);

Truth getTrueTrackletFromTrajectoryPoint(const std::vector<EDEPTrajectoryPoint>& points, double z);

double getScore(const TVectorD& tracklet, const std::vector<TVector3>& true_tracklet);
std::map<double, std::vector<TVector3>> getInterpolatedZ(const std::vector<EDEPTrajectoryPoint>& trj_points,
                                                         const sand_reco::kf::utils::TrackletMap& z_to_tracklets);

std::vector<double> computeZDistance(const std::vector<EDEPTrajectoryPoint>& trj_points,
                                                         const sand_reco::kf::utils::TrackletMap& z_to_tracklets);

sand_reco::kf::utils::TrackletMap findBestTracklet(const sand_reco::kf::utils::TrackletMap& z_to_tracklets,
                                                         const std::map<double, std::vector<TVector3>>& z_to_interpolated_tracklets);

// std::vector<TVector3> getTrueTrackletOfCluster(TVector3 start, TVector3 stop, double z);

