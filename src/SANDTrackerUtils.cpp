#include "SANDTrackerUtils.h"

TGeoManager* SANDTrackerUtils::fGeo = 0;
const double SANDTrackerUtils::kMagneticFieldInT = 0.6; /*T*/

const double SANDTrackerUtils::k = 0.299792458;
const double SANDTrackerUtils::c = SANDTrackerUtils::k * 1E3;  // mm/ns

void SANDTrackerUtils::Clear()
{
  // sand_reco::stt::stL.clear();
  // sand_reco::stt::stX.clear();
  // sand_reco::stt::stPos.clear();
  // sand_reco::stt::t0.clear();
  // sand_reco::stt::tubePos.clear();
}

// check if tubes are adjacent
// using tube id
bool SANDTrackerUtils::AreAdjacent(const SANDTrackerCellID &tub1, const SANDTrackerCellID &tub2)
{
  // To Do: this doesn't work for staggered stt tubes
  return (abs(long(tub1()) - long(tub2())) <= 1);
}


double SANDTrackerUtils::GetPerpMomentumInGeVFromRadiusInMM(double radius) {
  return GetRadiusInMMToMomentumInGeVConstant() * radius * GetMagneticField(); 
}

double SANDTrackerUtils::GetRadiusInMMFromPerpMomentumInGeV(double perpMom) {
  return perpMom / (GetRadiusInMMToMomentumInGeVConstant() * GetMagneticField()); 
}

TString SANDTrackerUtils::PrintMatrix(const TMatrixD& m) {
  TString str = "\n";
  for(int i = 0; i < m.GetNrows(); i++)
  {
    for(int j = 0; j < m.GetNcols(); j++)
    {
      str += TString::Format("%20.5f  ",m[i][j]);
    }
    str += "\n";
  }
  str += "\n";

  return str;
}

// TString SANDTrackerUtils::PrintStateVector(const SANDTrackerKFStateVector& v) {
//   TString str = "\n";
//   str += TString::Format("%20.5f\n",v.X());
//   str += TString::Format("%20.5f\n",v.Y());
//   str += TString::Format("%20.5f\n",v.SignedInverseRadius());
//   str += TString::Format("%20.5f\n",v.TanLambda());
//   str += TString::Format("%20.5f\n",v.Phi());
//   return str;
// }

// const SANDTrackerCluster* SANDTrackerUtils::GetClusterPointer(int clusterID, const std::vector<SANDTrackerCluster>& clusters)
// {
//   for(auto& cl: clusters)
//     if(cl.GetId() == clusterID)
//       return &cl;
  
//   return 0;
// }

TVector3 SANDTrackerUtils::GetCartesianCoordinateFromCylindrical(double radius, double angle, double x)
{
  auto sandCenter = SANDTrackerUtils::GetSANDInnerVolumeCenterPosition();
  auto xSandCenter = sandCenter[0];
  auto ySandCenter = sandCenter[1];
  auto zSandCenter = sandCenter[2];
  
  auto y = ySandCenter + radius * sin(angle);
  auto z = zSandCenter + radius * cos(angle);

  return {x,y,z};
}













namespace STTKFCheck {

TVector2 get_Bfield_perp(const TVector3& v)
{
    return TVector2(v.Z(), v.Y());
}

double get_Bfield_long(const TVector3& v)
{
    return v.X();
}

double get_tan_of_dip_angle(const TVector3& mom)
{
    auto perp = get_Bfield_perp(mom);
    auto para = get_Bfield_long(mom);
    return para/perp.Mod();
}

int get_rotation_versus(int charge)
{
    return -charge;
}

int get_charge(int versus)
{
    return -versus;
}

double get_direction_angle(const TVector2& dir)
{
    // dir.Print();
    return TMath::ATan2(dir.Y(), dir.X());
}

double get_rotation_angle(const TVector2& dir, int charge)
{
    auto theta = get_direction_angle(dir);
    // std::cout << "THETA " << theta << std::endl;
    //if (theta < 0 ) theta +=  2 * TMath::Pi();
    auto versus = get_rotation_versus(charge);
    //return (theta + versus * 0.5 * TMath::Pi());
    return theta + versus * 0.5 * TMath::Pi();
}

double get_radius(double perp_mom)
{
    return SANDTrackerUtils::GetRadiusInMMFromPerpMomentumInGeV(perp_mom);
}

TVector2 get_circle_center(const TVector2& perp_mom, const TVector2& position, int charge)
{
    auto phi = get_rotation_angle(perp_mom, charge);
    auto radius = get_radius(perp_mom.Mod());
    // std::cout << "weew: " << phi << " " << radius << std::endl;
    return TVector2(position.X() - radius * TMath::Cos(phi), position.Y() - radius * TMath::Sin(phi));
}

std::pair<double, double> get_circle_ys(double z, double radius, const TVector2& center)
{
    auto y = TMath::Sqrt(radius * radius - (z - center.X()) * (z - center.X()));
    return {center.Y() + y, center.Y() - y};
}

double get_rotation_angle(double z, double y, const TVector2& center)
{
    return TMath::ATan2(y - center.Y(), z - center.X());
}

double get_delta_phi(double phi, double previous_phi, int versus)
{

    // da riguardare e sempilificare
    phi = phi < 0. ? 2 * TMath::Pi() + phi : phi;
    previous_phi = previous_phi < 0. ? 2 * TMath::Pi() + previous_phi : previous_phi;
    auto delta_phi = phi - previous_phi;

    auto dphi = delta_phi * versus >= 0. ? delta_phi : 2 * TMath::Pi() - delta_phi;
    // return dphi;
    return fmod(dphi, M_PI);
}

double get_x(double radius, double x_0, double delta_phi, double tan_lambda, int versus_of_rot)
{
    auto charge = get_charge(versus_of_rot);
    return x_0 + radius * charge * tan_lambda * delta_phi;
}

double get_y(const TVector2& center, double radius, double phi)
{   
    // std::cout << "deltaY: " << radius * TMath::Sin(phi) << std::endl;
    return center.Y() + radius * TMath::Sin(phi);
}

double get_z(const TVector2& center, double radius, double phi)
{
    // std::cout << "deltaZ: " << radius * TMath::Cos(phi) << std::endl;
    return center.X() + radius * TMath::Cos(phi);
}

TVector3 get_vector_momentum(double radius, double phi, double tan_lambda, int versus)
{
    auto perp_mon = SANDTrackerUtils::GetPerpMomentumInGeVFromRadiusInMM(radius);
    auto theta = phi - versus * 0.5 * TMath::Pi();
    return TVector3(perp_mon * tan_lambda, perp_mon * TMath::Sin(theta), perp_mon * TMath::Cos(theta));
}

STTKFStateVector get_state_vector(TVector3 mom, TVector3 pos, int charge)
{
    auto perp_mom = get_Bfield_perp(mom);
    auto perp_pos = get_Bfield_perp(pos);
    auto radius = get_radius(perp_mom.Mod());
    // perp_mom.Print();
    auto tan_lambda = get_tan_of_dip_angle(mom);
    auto phi = get_rotation_angle(perp_mom, charge);

    return STTKFStateVector(pos.X(), pos.Y(), charge/radius, tan_lambda, phi);
}

ParticleState::ParticleState(const STTKFStateVector& vector, double z)
{
    _position = TVector3(vector.X(), vector.Y(), z);

    auto charge = vector.Charge();
    auto versus = get_rotation_versus(charge);
    auto radius = vector.Radius();
    
    _momentum = get_vector_momentum(radius, vector.Phi(), vector.TanLambda(), versus);
}

TrajectoryParameters ParticleState::get_trajectory_parameter(int charge) const
{
    auto perp_mom = get_Bfield_perp(_momentum);
    auto perp_pos = get_Bfield_perp(_position);
    auto radius = get_radius(perp_mom.Mod());
    auto versus_of_rot = get_rotation_versus(charge);
    auto tan_lambda = get_tan_of_dip_angle(_momentum);
    auto center_of_rot = get_circle_center(perp_mom, perp_pos, charge);
    auto phi_0 = get_rotation_angle(perp_mom, charge);
    auto x_0 = _position.X();

    return TrajectoryParameters(radius, versus_of_rot, tan_lambda, center_of_rot, phi_0, x_0);
}

STTKFStateVector ParticleState::get_state_vector(int charge) const
{
    auto perp_mom = get_Bfield_perp(_momentum);
    auto perp_pos = get_Bfield_perp(_position);
    auto radius = get_radius(perp_mom.Mod());
    // perp_mom.Print();
    auto tan_lambda = get_tan_of_dip_angle(_momentum);
    auto phi = get_rotation_angle(perp_mom, charge);

    return STTKFStateVector(_position.X(), _position.Y(), charge/radius, tan_lambda, phi);
}

std::pair<double, double> TrajectoryParameters::get_phi_pair(double z) const {
    double dy = TMath::Sqrt(_radius * _radius - (_center_of_rot.X() - z) * (_center_of_rot.X() - z));
    return std::pair<double, double>(atan2(dy, (z - _center_of_rot.X())), atan2(-dy, (z - _center_of_rot.X())));
}

ParticleState TrajectoryParameters::get_particle_state(double delta_phi) const
{
    auto phi = _phi_0 + delta_phi;
    auto x = get_x(_radius, _x_0, delta_phi, _tan_lambda, _versus_of_rot);
    auto y = get_y(_center_of_rot, _radius, phi);
    auto z = get_z(_center_of_rot, _radius, phi);
    auto p = TVector3{x,y,z};
    auto m = get_vector_momentum(_radius, phi, _tan_lambda, _versus_of_rot);

    return ParticleState(p, m);
}

std::vector<ParticleState> TrajectoryParameters::get_particle_states_from_delta_phi(std::vector<double> delta_phis) const
{
    std::vector<ParticleState> particle_states;
    auto total_delta_phi = 0.;
    for(auto delta_phi: delta_phis)
    {
        total_delta_phi += delta_phi;
        particle_states.emplace_back(get_particle_state(total_delta_phi));
    }
    return particle_states;
}

double TrajectoryParameters::get_smallest_delta_phi(double z, double last_phi) const
{
    auto phi_pair = get_phi_pair(z);
    auto delta_phi_1 = get_delta_phi(phi_pair.first , last_phi, _versus_of_rot);
    auto delta_phi_2 = get_delta_phi(phi_pair.second, last_phi, _versus_of_rot);

    // debug ////////////////////////////////////////
    // std::cout << z << " " << last_phi << " " << _versus_of_rot << " "
    //           << phi_pair.first << " " << phi_pair.second << " "
    //           << delta_phi_1 << " " << delta_phi_2 << " "
    //           << (delta_phi_1 < delta_phi_2 ? delta_phi_1 : delta_phi_2) << std::endl; 
    /////////////////////////////////////////////////

    return fabs(delta_phi_1) < fabs(delta_phi_2) ? delta_phi_1 : delta_phi_2;
}

std::vector<double> TrajectoryParameters::get_delta_phis(std::vector<double> zs) const
{
    std::vector<double> delta_phis;

    auto last_phi = _phi_0;
    auto delta_phi = 0.;

    for(auto& z: zs)
    {
        delta_phi = get_smallest_delta_phi(z, last_phi);

        // debug //////////////////////////////////////
        // std::cout << z << " " << last_phi << " " << delta_phi << std::endl;
        //////////////////////////////////////////////

        delta_phis.push_back(delta_phi);
        last_phi += delta_phi;
    }
    return delta_phis;
}

std::vector<ParticleState> TrajectoryParameters::get_particle_states_from_z(std::vector<double> zs) const
{
        // get delta phi
        auto delta_phis = get_delta_phis(zs);

        // debug /////////////////////////////
        // for(auto i = 0u; i < zs.size(); i++) std:: cout << zs.at(i) << " " << delta_phis.at(i) << std::endl;
        /////////////////////////////////////

        return get_particle_states_from_delta_phi(delta_phis);
}

ParticleState operator -(const ParticleState& p1, const ParticleState& p2) {
    ParticleState p;
    p.set_position(p1.get_position() - p2.get_position());
    p.set_momentum(p1.get_momentum() - p2.get_momentum());
    return p;
}
bool operator ==(const ParticleState& p1, const ParticleState& p2) {
    
    return (p1.get_position() == p2.get_position() && p1.get_momentum() == p2.get_momentum());
}

}