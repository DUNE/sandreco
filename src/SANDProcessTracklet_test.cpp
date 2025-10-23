#include "SANDProcessTracklets.h"

double ComputeStd(const std::vector<double>& v){
  auto [sum_x, sum_x2] = std::accumulate(
    v.begin(), v.end(),
    std::make_pair(0LL, 0LL),
    [](auto a, auto e) {
      return std::make_pair(a.first + e, a.second + e * e);
    }
  );
  sum_x /= static_cast<double>(v.size());
  sum_x2 /= static_cast<double>(v.size());

  return std::sqrt(sum_x2 - sum_x * sum_x);
}

// ---------------------------------------------------------------------------------------------------------------------------
// Compute a score to select the best tracklet among all based on position and direction
// ---------------------------------------------------------------------------------------------------------------------------
double getScore(const TVectorD& tracklet, const Truth& truth){

    const TVector3& true_pos = truth.pos; //best trajectory position
    const TVector3& true_dir = truth.dir; //best trajectory direction

    const double true_theta_yz = atan(true_dir.Y() / true_dir.Z());
    const double true_theta_xz = atan(true_dir.Z() / true_dir.X());
    if (true_theta_xz < 0) { true_theta_xz = M_PI + true_theta_xz;} 

    const double x_trk = tracklet[0];
    const double y_trk = tracklet[1];
    const double theta_xz = tracklet[2];
    const double theta_yz = tracklet[3];

    //Find the closest (x,y)
    double position_distance = 
        std::sqrt( (x_trk - pos_true.X())*(x_trk - pos_true.X())
            + (y_trk - pos_true.Y())*(y_trk - pos_true.Y()) );
    
    //Find the best direction
    const double px = std::cos(theta_xz);
    const double py = std::sin(theta_yz);
    double pz_sq = 1.0 - px*px - py*py;
    if (pz_sq < 0.0) pz_sq = 0.0;
    const double pz = std::sqrt(pz_sq);
    const TVector3 dir_trk(px, py, pz); //traklet direction
    //angle between traklet and trajectory
    double dot = dir_trk.Dot(dir_true); //rad (spero)
    if (dot >  1.0) dot = 1.0;
    if (dot < -1.0) dot = -1.0;

    const double sigma_pos = SANDTrackerUtils::getSigmaPositionMeasurement() * 1E3; // mm
    const double sigma_ang = SANDTrackerUtils::getSigmaAngleMeasurement();          // rad

    const double score = position_distance / sigma_pos + std::acos(dot) / sigma_ang;

    return score;
    
}

// ---------------------------------------------------------------------------------------------------------------------------
// Takes the start/stop position and momenta of the deposited segment and the z of the cluster to find a truth at that z by 
// interpolating position and momentum
// ---------------------------------------------------------------------------------------------------------------------------
Truth getTrueTrackletOfCluster(TVector3 start_pos, TVector3 stop_pos,
                               TVector3 start_mom, TVector3 stop_mom,
                               double z)
{
  // Interpolation of z coordinate
  auto first_point = start_pos;
  auto second_point = stop_pos;

   if (first_point.Z() > second_point.Z()) {
    std::swap(first_point, second_point);
    std::swap(start_mom, stop_mom);
  }

  double dz = (second_point.Z() - z);
  double alpha = dz / (second_point.Z() - first_point.Z());

  //Interpolation of position and momentum
  TVector3 pos = first_point * (1.0 - alpha) + second_point * alpha;
  TVector3 mom = start_mom * (1.0 - alpha) + stop_mom * alpha
  TVector3 dir = mom_.Unit();

  return Truth{pos, dir, mom;};
}


// ---------------------------------------------------------------------------------------------------------------------------
// Create a map of with key z and value a Tracklet obtained from TG4HitSegment 
// ---------------------------------------------------------------------------------------------------------------------------
std::map<double, Truth> getInterpolatedZFromSegments(TVector3 start_pos, TVector3 stop_pos,
                                                    TVector3 start_mom, TVector3 stop_mom,
                                                     const sand_reco::kf::utils::TrackletMap& z_to_tracklets)
{

  std::map<double, Truth> z_to_truth_from_seg;

  auto first_point = start_pos;
  auto second_point = stop_pos;

   if (first_point.Z() > second_point.Z()) {
    std::swap(first_point, second_point);
    std::swap(start_mom, stop_mom);
  }

  double dz = (second_point.Z() - first_point.Z());
  double alpha = dz / (second_point.Z() - first_point.Z());

  for (auto it = z_to_tracklets.begin(); it != z_to_tracklets.end(); ++it) {
      const double z_trk = z.first;
      
      //Interpolation of position and momentum
      TVector3 pos = first_point * (1.0 - alpha) + second_point * alpha;
      TVector3 mom = start_mom * (1.0 - alpha) + stop_mom * alpha
      TVector3 dir = mom_.Unit();

      z_to_truth_from_segments[z] = Truth{pos, dir, mom};

      }
  return z_to_truth_from_segments;
}



// ---------------------------------------------------------------------------------------------------------------------------
// Create a map of with key z and value a Tracklet obtained from the trajectory points 
// // ---------------------------------------------------------------------------------------------------------------------------
// std::map<double, Truth> getInterpolatedZFromPoints(
//     const std::vector<EDEPTrajectoryPoint>& trj_points,
//     const sand_reco::kf::utils::TrackletMap& z_to_tracklets) //Tracklet mao è ottenuta dai clsuter dei digit
// {

//   std::map<double, Truth> interpolated_z;

//   //Find the closest z coordinates of the MC trajectory to the tracklet
//   for (const auto& z : z_to_tracklets) {
//     double z_trk = z.first;
//     uint i_min = 0;
//     double z_min = 1e30;

//     for(uint i = 0; i < trj_points.size(); i++ ){ 
//       auto point = trj_points.at(i);
//       double z_trj = point.GetPosition().Vect().Z();
//       double dz = std::fabs(z_trk - z_trj);

//       if(dz < z_min){z_min = dz; i_min = i;}
//     }
  
//     //Second point closest 
//     int i_second_min = 10E8;
//     if(i_min == trj_points.size() - 1 ){
//       i_second_min = i_min - 1;
//     } else if (i_min == 0) {
//       i_second_min = i_min + 1;
//     } else {
//       auto prev_point = trj_points.at(i_min - 1);
//       auto next_point = trj_points.at(i_min + 1);
//       double prev_z = prev_point.GetPosition().Vect().Z();
//       double next_z = next_point.GetPosition().Vect().Z();
//       double prev_dz = std::fabs(z_trk - prev_z);
//       double next_dz = std::fabs(z_trk - next_z);

//       if (prev_dz < next_dz) {
//         i_second_min = i_min - 1;
//       } else {
//         i_second_min = i_min +   1;
//       }
//     }

//     //Ordering points along z direction (pA=point_min pB=point_second_min)
//     auto pA = trj_points[i_min];
//     auto pB = trj_points[i_second_min];
    
//     if (pA.GetPosition().Z() > pB.GetPosition().Z()) std::swap(pA, pB);

//     double zA = pA.GetPosition().Z();
//     double zB = pB.GetPosition().Z();
//     double dz = z_trk - zA;
//     double alpha = dz / std::fabs(zB -zA);
    
//     TVector3 posA = pA.GetPosition().Vect();
//     TVector3 posB = pB.GetPosition().Vect();
//     TVector3 momA = pA.GetMomentum();
//     TVector3 momB = pB.GetMomentum();

//     TVector3 interpolated_pos = posA * (1.0 - alpha) + posB * alpha;
//     TVector3 interpolated_mom = momA * (1.0 - alpha) + momB * alpha;
//     TVector3 interpolated_dir = interpolated_mom.Unit();

//     interpolated_z[z_trk] = Truth{interpolated_pos, interpolated_dir, interpolated_mom};  
//   }
//   return interpolated_z;
// }


//----------------------------------------------------------------------------------------------------------------------------
// Find the closest z coordinates of the MC trajectory to the tracklet
// ---------------------------------------------------------------------------------------------------------------------------
std::vector<double> computeZDistance(const std::vector<EDEPTrajectoryPoint>& trj_points,
                                                         const std::map<double, std::vector<TVectorD>>& z_to_tracklets){

  std::vector<double> dz_vec;
  double dz;
  for (const auto& z : z_to_tracklets) {
    double z_trk = z.first;
    uint i_min = 0;
    double z_min = 10E8;

    for(uint i = 0; i < trj_points.size(); i++ ){ 
      auto point = trj_points.at(i);
      double z_trj = point.GetPosition().Vect().Z();
      double current_dz = fabs(z_trk - z_trj);

      if(current_dz < z_min){
        z_min = current_dz;
        i_min = i;
      }
    }

    int i_second_min = 10E8;
    if(i_min == trj_points.size() - 1 ){
      i_second_min = i_min - 1;
    } else if (i_min == 0) {
      i_second_min = i_min + 1;
    } else {
      auto prev_point = trj_points.at(i_min - 1);
      auto next_point = trj_points.at(i_min + 1);
      double prev_z = prev_point.GetPosition().Vect().Z();
      double next_z = next_point.GetPosition().Vect().Z();

      double prev_dz = fabs(z_trk - prev_z);
      double next_dz = fabs(z_trk - next_z);

      if (prev_dz < next_dz) {
        i_second_min = i_min - 1;
      } else {
        i_second_min = i_min +   1;
      }
    }

    //Interpolation of z coordinate
    auto pA = trj_points[i_min];
    auto pB = trj_points[i_second_min];
    
    auto first_point = pA;
    auto second_point = pB;
    if (pA.GetPosition().Z() < pB.GetPosition().Z()) {
      first_point  = pA;
      second_point = pB;
    } else {
      first_point  = pB;
      second_point = pA;
    }

    dz = fabs(first_point.GetPosition().Z() - second_point.GetPosition().Z()); 
  }
 dz_vec.push_back(dz);
return dz_vec;
}

// ---------------------------------------------------------------------------------------------------------------------------
// Select the best tracklet for each cluster based on the best score 
// ---------------------------------------------------------------------------------------------------------------------------
sand_reco::kf::utils::TrackletMap findBestTracklet(
    const sand_reco::kf::utils::TrackletMap& z_to_tracklets,
    const std::map<double, Truth>& z_to_truth)
{
  sand_reco::kf::utils::TrackletMap z_to_best_tracklet;
  if (z_to_tracklets.empty() || z_to_truth.empty())
    return z_to_best_tracklet;

  for (std::map<double, Truth>::const_iterator it_truth = z_to_truth.begin();
       it_truth != z_to_truth.end(); ++it_truth)
  {
    const double z = it_truth->first;
    const Truth& truth = it_truth->second;

    sand_reco::kf::utils::TrackletMap::const_iterator it = z_to_tracklets.find(z);
    if (it == z_to_tracklets.end() || it->second.empty())
      continue;

    const std::vector<TVectorD>& candidates = it->second;

    double best_score = 1e30;
    TVectorD best = candidates.front();

    for (std::size_t i = 0; i < candidates.size(); ++i) {
      const double s = getScore(candidates[i], truth);
      if (s < best_score) {
        best_score = s;
        best = candidates[i];
      }
    }

    z_to_best_tracklet[z].push_back(best);
  }

  return z_to_best_tracklet;
}