#pragma once


#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <unordered_map>
#include <vector>
#include <array>

#include "TVector3.h"
#include "TVectorD.h"

#include "utils.h"
#include "SANDTrackerUtils.h"
#include "SANDTrackerCluster.h"
#include "SANDClustering.h"
#include "SANDWireInfo.h"
#include "SANDGeoManager.h"
#include "SANDKFTrack.h"


class SANDGeoManager;

namespace sand_reco
{
namespace kf
{
namespace trackbuilding
{

static const double stereoAngleInRad = 5.0*M_PI/180.0;

//a Triplet is defined as a sapce point(x,y,z) reconstructed from one cluster of the three drift chamber plane (U,V,Y)
struct Triplet{
    int moduleID;
    TVector3 triplet_pos_;
    double score;
    double D;

    int nU = 0;
    int nV = 0;
    int nY = 0;

    const sand_reco::tracker::Cluster* clusterU =   nullptr;
    const sand_reco::tracker::Cluster* clusterV =   nullptr;
    const sand_reco::tracker::Cluster* clusterY =   nullptr;   
};

// std::map<double , Triplet> z_to_triplet(const std::vector<EDEPTrajectoryPoint>& points,
//                                         const SANDGeoManager* sand_geo,
//                                         const sand_reco::tracker::ClusterCollection& clusters,
//                                         double step);                                           

using TripletMap = std::map<double, Triplet>;

struct TrackSeed{
    sand_reco::kf::State  initial_state_;
    std::array<const Triplet* ,3> triplet_seed_; //to build the seed we use three triplet from three consecutive modules
    double chi2_ = 0.0;
};

// a Track Candidate is defined as a collection of triplets ordered along z with its seed 
struct TrackCandidate{
    TrackSeed seed_;
    std::vector<const Triplet*> triplets_;
    //we should chooose some criteria to choose if a track candidate is valid ot not (chi2, number of missing planes, minimum number of hits?)
    double chi2_ = 0.0;
    int n_missing_pmales_;
};


//From the available clusters in the tracker it returns a map of triplets along for each of them (i.e. a triplet for each module)
TripletMap buildTripletFromClusters(const SANDGeoManager* sand_geo,
                                     const std::vector<sand_reco::tracker::Cluster>& clusters,
                                     double stereo_angle);
//It finds a seed from 3 consecutive triplets (ideal case: one triplet for three adjacent module)
std::vector<TrackSeed> findSeedFromTriplets(const TripletMap& z_to_triplet,
                                            double chi2, 
                                            double n_planes
                                            );

TrackCandidate buildtrackFromSeed(const TripletMap& z_to_triplet,
                                  const TrackSeed& seed,
                                  double chi2
                                );


}//sandreco
}//kf
}//trackbuilding