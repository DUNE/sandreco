#include "SANDProcessTracklets.h"

double ComputeStd(const std::vector<double>& values, double mean) {
  double squared_difference = 0.0;
  for (double v : values) {
      squared_difference += (v - mean) * (v - mean);
  }
  return std::sqrt(squared_difference / values.size());
}

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

  double dz = z - first_point.Z();
  double alpha = dz / (second_point.Z() - first_point.Z());

  //Interpolation of position and momentum
  TVector3 pos = first_point * (1.0 - alpha) + second_point * alpha;
  TVector3 mom = start_mom * (1.0 - alpha) + stop_mom * alpha;
  TVector3 dir = mom.Unit();

  return Truth{pos, dir, mom};
}



Truth getTrueTrackletFromTrajectoryPoint(const std::vector<EDEPTrajectoryPoint>& points, double z){
  const EDEPTrajectoryPoint* best_point = &points.front();
  double best_z = std::abs(points.front().GetPosition().Z() - z);

  for (const auto& tp : points) {
    double dz = std::abs(tp.GetPosition().Z() - z);
    if (dz < best_z) {
      best_z = dz;
      best_point = &tp;
    }
  }
  TVector3 true_pos = best_point->GetPosition().Vect();
  TVector3 true_dir = best_point->GetMomentum().Unit();
  TVector3 true_mom = best_point->GetMomentum();

  return Truth{true_pos, true_dir, true_mom};
}

std::map<double, Truth> z_to_truth(const std::vector<EDEPTrajectoryPoint>& points, double step){
  std::map<double, Truth> z_to_trajectory_point;
  if (points.empty()) return z_to_trajectory_point;

  double z_min = points[0].GetPosition().Z();
  double z_max = z_min;
  for (std::size_t i = 1; i < points.size(); ++i) {
    const double z = points[i].GetPosition().Z();
    if (z < z_min) z_min = z;
    if (z > z_max) z_max = z;
  }

  double z = z_min;
  while (z <= z_max) {
    z_to_trajectory_point[z] = getTrueTrackletFromTrajectoryPoint(points, z);
    z += step;
  }
  if (z_to_trajectory_point.find(z_max) == z_to_trajectory_point.end()) {
    z_to_trajectory_point[z_max] = getTrueTrackletFromTrajectoryPoint(points, z_max);
  }
  return z_to_trajectory_point;
}



double getScore(const TVectorD& tracklet, const std::vector<TVector3>& true_tracklet){

  //Select the best traklet based on position and direction
  TVector3 best_trj_point = true_tracklet.at(0);
  TVector3 p_trj_dir = true_tracklet.at(1).Unit();
  
  double true_theta_yz = atan(p_trj_dir.Y() / p_trj_dir.Z());
  double true_theta_xz = atan(p_trj_dir.Z() / p_trj_dir.X());
  if (true_theta_xz < 0) {
    true_theta_xz = M_PI + true_theta_xz;
  } 

  double x_trk = tracklet[0];
  double y_trk = tracklet[1];
  double theta_xz = tracklet[2];
  double theta_yz = tracklet[3];

  //Find the closest (x,y)
  double position_distance = sqrt(pow(x_trk - best_trj_point.X(), 2) + pow(y_trk - best_trj_point.Y(), 2));
  double anglular_distance_xz = fabs(true_theta_xz - theta_xz);
  double anglular_distance_yz = fabs(true_theta_yz - theta_yz);

  //Find the best direction
  // double px_trk = cos(theta_xz);
  // double py_trk = sin(theta_yz);
  // double pz_trk = sqrt(1 - px_trk * px_trk - py_trk * py_trk);

  // TVector3 p_trk_dir(px_trk, py_trk, pz_trk);    

  //to be precise this is the cos of the angle between the trajectory and the traklet    
  // p_trk_dir.Print();
  // p_trj_dir.Print();
  // double direction = p_trk_dir.Dot(p_trj_dir); 

  double score = position_distance / (SANDTrackerUtils::getSigmaPositionMeasurement() * 1E3)  // mm
               + anglular_distance_xz / SANDTrackerUtils::getSigmaAngleMeasurement()          // rad
               + anglular_distance_yz / SANDTrackerUtils::getSigmaAngleMeasurement();         // rad

  return score;
}


std::map<double, std::vector<TVector3>> getInterpolatedZ(const std::vector<EDEPTrajectoryPoint>& trj_points,
                                                         const std::map<double, std::vector<TVectorD>>& z_to_tracklets){

  std::map<double, std::vector<TVector3>> interpolated_z;

  //Find the closest z coordinates of the MC trajectory to the tracklet
  for (const auto& z : z_to_tracklets) {
    double z_trk = z.first;
    uint i_min = 0;
    double z_min = 10E8;

    for(uint i = 0; i < trj_points.size(); i++ ){ 
      auto point = trj_points.at(i);
      double z_trj = point.GetPosition().Vect().Z();
      double z_distance = fabs(z_trk - z_trj);

      if(z_distance < z_min){
        z_min = z_distance;
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

      double prev_z_distance = fabs(z_trk - prev_z);
      double next_z_distance = fabs(z_trk - next_z);

      if (prev_z_distance < next_z_distance) {
        i_second_min = i_min - 1;
      } else {
        i_second_min = i_min +   1;
      }
    }

    //Interpolation of z coordinate
    auto point_min        = trj_points[i_min];
    auto point_second_min = trj_points[i_second_min];
    
    auto first_point = point_min;
    auto second_point = point_second_min;
    if (point_min.GetPosition().Z() < point_second_min.GetPosition().Z()) {
      first_point  = point_min;
      second_point = point_second_min;
    } else {
      first_point  = point_second_min;
      second_point = point_min;
    }

    double z_diff = z_trk - first_point.GetPosition().Z();
    double alpha = z_diff / fabs(second_point.GetPosition().Z() - first_point.GetPosition().Z());

     

    TVector3 interpolated_pos = first_point.GetPosition().Vect() * (1 - alpha) + second_point.GetPosition().Vect() * alpha;
    TVector3 interpolated_mom = first_point.GetMomentum() * (1 - alpha) + second_point.GetMomentum() * alpha;
    TVector3 interpolated_dir = interpolated_mom.Unit();

    interpolated_z[z.first] = { interpolated_pos, interpolated_dir };  
  }
  return interpolated_z;
}

std::vector<double> computeZDistance(const std::vector<EDEPTrajectoryPoint>& trj_points,
                                                         const std::map<double, std::vector<TVectorD>>& z_to_tracklets){

  std::vector<double> z_distance_vec;
  double z_distance;
  //Find the closest z coordinates of the MC trajectory to the tracklet
  for (const auto& z : z_to_tracklets) {
    double z_trk = z.first;
    uint i_min = 0;
    double z_min = 10E8;

    for(uint i = 0; i < trj_points.size(); i++ ){ 
      auto point = trj_points.at(i);
      double z_trj = point.GetPosition().Vect().Z();
      double current_z_distance = fabs(z_trk - z_trj);

      if(current_z_distance < z_min){
        z_min = current_z_distance;
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

      double prev_z_distance = fabs(z_trk - prev_z);
      double next_z_distance = fabs(z_trk - next_z);

      if (prev_z_distance < next_z_distance) {
        i_second_min = i_min - 1;
      } else {
        i_second_min = i_min +   1;
      }
    }

    //Interpolation of z coordinate
    auto point_min        = trj_points[i_min];
    auto point_second_min = trj_points[i_second_min];
    
    auto first_point = point_min;
    auto second_point = point_second_min;
    if (point_min.GetPosition().Z() < point_second_min.GetPosition().Z()) {
      first_point  = point_min;
      second_point = point_second_min;
    } else {
      first_point  = point_second_min;
      second_point = point_min;
    }

    z_distance = fabs(first_point.GetPosition().Z() - second_point.GetPosition().Z()); 
  }
 z_distance_vec.push_back(z_distance);
return z_distance_vec;
}



std::map<double, std::vector<TVectorD>> findBestTracklet(const std::map<double, std::vector<TVectorD>>& z_to_tracklets,
                                       const std::map<double, std::vector<TVector3>>& z_to_interpolated_tracklets){

  std::map<double, std::vector<TVectorD>> z_to_best_tracklets;

  for (auto& current_z:z_to_interpolated_tracklets) {
    const auto& tracklets_at_current_z = z_to_tracklets.at(current_z.first);


    //Select the best traklet based on position and direction
    TVector3 best_trj_point = current_z.second.at(0);
    TVectorD best_tracklet(tracklets_at_current_z.at(0).GetNrows());
    TVector3 p_trj_dir = current_z.second.at(1).Unit();

    std::vector<double> position_errors, direction_errors;
    
    double best_score = 1e8;
    for (const auto& tracklet : tracklets_at_current_z) {
      double x_trk = tracklet[0];
      double y_trk = tracklet[1];
      double theta_xz = tracklet[2];
      double theta_yz = tracklet[3];

      //Find the closest (x,y)
      double position_distance = sqrt(pow(x_trk - best_trj_point.X(), 2) + pow(y_trk - best_trj_point.Y(), 2));
      
      //Find the best direction
      double px_trk = cos(theta_xz);
      double py_trk = sin(theta_yz);
      double pz_trk = sqrt(1 - px_trk * px_trk - py_trk * py_trk);

      TVector3 p_trk_dir(px_trk, py_trk, pz_trk);    
      
      //to be precise this is the cos of the angle between the trajectory and the traklet              
      double direction = p_trk_dir.Dot(p_trj_dir); 

      position_errors.push_back(position_distance);
      direction_errors.push_back(direction);

      double score = position_distance / 200E-3 + acos(direction) / 0.02;

      // std::cout << "z: " << current_z.first<< " position distance  "  
      //             << position_distance 
      //             << " angular distance  "  
      //             << direction
      //             << " SCORE " 
      //             << score
      //             << " x distance: " << x_trk - best_trj_point.X()
      //             << " y distance: " << y_trk - best_trj_point.Y()
      //             << std::endl;
      if (score < best_score) {
        best_score = score;
        best_tracklet = tracklet;
        // std::cout << "BEST SCORE " 
        //           << best_score
        //           << std::endl;
      }
    }

    z_to_best_tracklets[current_z.first].push_back(best_tracklet);
  }
  return z_to_best_tracklets;
}

