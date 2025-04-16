#include "SANDTrackerUtils.h"

TGeoManager* SANDTrackerUtils::geo_ = 0;
const double SANDTrackerUtils::kMagneticFieldInT_ = 0.6; /*T*/
const double SANDTrackerUtils::kEdepSimDensityToGCM3_ = 6.24E18;

const double SANDTrackerUtils::k_ = 0.299792458;
const double SANDTrackerUtils::c_ = SANDTrackerUtils::k_ * 1E3;  // mm/ns

void SANDTrackerUtils::clear()
{
  // sand_reco::stt::stL.clear();
  // sand_reco::stt::stX.clear();
  // sand_reco::stt::stPos.clear();
  // sand_reco::stt::t0.clear();
  // sand_reco::stt::tubePos.clear();
}

// check if tubes are adjacent
// using tube id
// To Do: replace this with the adjacency of cells defined during
// SANDGeoManager construction
bool SANDTrackerUtils::areAdjacent(const sand_geometry::tracker::CellID &tub1, const sand_geometry::tracker::CellID &tub2)
{
  // To Do: this doesn't work for staggered stt tubes
  return (abs(long(tub1()) - long(tub2())) <= 1);
}


double SANDTrackerUtils::getPerpMomentumInGeVFromRadiusInMM(double radius) {
  return getRadiusInMMToMomentumInGeVConstant() * radius * getMagneticField(); 
}

double SANDTrackerUtils::getRadiusInMMFromPerpMomentumInGeV(double perpMom) {
  return perpMom / (getRadiusInMMToMomentumInGeVConstant() * getMagneticField()); 
}

TString SANDTrackerUtils::printMatrix(const TMatrixD& m) {
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

double SANDTrackerUtils::getCrossedMaterialInGCM2(double z, 
                                           double px, double py, double pz,
                                           double sx, double sy, double sz) {
  // sx, sy, sz represent the directional cosine 
  // and should be correctly normalized:
  // sx^2 + sy^2 + sz^2 = 1
  // The function return the ammount of
  // crossed material in g/cm2
  geo_->InitTrack(px, py, pz, sx, sy, sz);
  const double* lastPosition = 0;

  double crossedMaterial = 0.;
  // std::cout << __LINE__ << std::endl;
  // std::cout << sx << " " << sy << " " << sz << std::endl;
  // std::cout << geo_->GetCurrentPoint()[2] << " " << z << std::endl;
  while((lastPosition = geo_->GetCurrentPoint()) && lastPosition[2] > z) {
    auto density = getDensityInGCM3();
    auto pathLength = getPathLengthInCM();
    crossedMaterial += density * pathLength;
    // std::cout << density << " " << pathLength << std::endl;
    geo_->Step();
  }
  return crossedMaterial;
}

double SANDTrackerUtils::getPathLengthInX0(double z, 
                                           double px, double py, double pz,
                                           double sx, double sy, double sz) {

  geo_->InitTrack(px, py, pz, sx, sy, sz);
  const double* lastPosition = 0;

  double pathLengthInX0 = 0.;
  int count = 0;
  while((lastPosition = geo_->GetCurrentPoint()) && lastPosition[2] > z) {
    auto Z = static_cast<int>(geo_->GetCurrentNode()->GetVolume()->GetMaterial()->GetZ());
    auto A = static_cast<int>(geo_->GetCurrentNode()->GetVolume()->GetMaterial()->GetA());
    auto name = static_cast<std::string>(geo_->GetCurrentNode()->GetVolume()->GetMaterial()->GetName());
    count++;
    auto X0 = SANDTrackerUtils::getX0(Z, A);
    auto density = getDensityInGCM3();
    auto pathLength = getPathLengthInCM();
    pathLengthInX0 += pathLength * density / X0;
    geo_->Step();
  }
  return pathLengthInX0;
}

double SANDTrackerUtils::getPathLengthInCM(double z, 
                                           double px, double py, double pz,
                                           double sx, double sy, double sz) {
  // sx, sy, sz represent the directional cosine 
  // and should be correctly normalized:
  // sx^2 + sy^2 + sz^2 = 1
  // The function return the ammount of
  // crossed material in g/cm2
  geo_->InitTrack(px, py, pz, sx, sy, sz);
  const double* lastPosition = 0;

  double pathLengthInCM = 0.;
  while((lastPosition = geo_->GetCurrentPoint()) && lastPosition[2] > z) {
    auto pathLength = getPathLengthInCM();
    pathLengthInCM += pathLength;
    geo_->Step();
  }
  return pathLengthInCM;
}


double SANDTrackerUtils::getDE(double z, 
                              double px, double py, double pz,
                              double sx, double sy, double sz,
                              double beta, double mass, int charge) {

  double K = 0.307075; // MeV mol−1 cm2
  double m_e = 0.5109989461; // MeV  
  double gamma = 1 / sqrt(1 - beta*beta);
  double W_max = (2 * m_e * beta * beta * gamma* gamma) /
                 (1 + 2 * gamma * m_e / mass + pow(m_e / mass, 2));

  geo_->InitTrack(px, py, pz, sx, sy, sz);
  const double* lastPosition = 0;

  double dE = 0.;
  while((lastPosition = geo_->GetCurrentPoint()) && lastPosition[2] > z) {
    geo_->FindNextBoundary();
    auto density = getDensityInGCM3();
    auto pathLength = getPathLengthInCM();

    auto current_material = geo_->GetCurrentNode()->GetVolume()->GetMaterial();
    auto name = static_cast<std::string>(current_material->GetName());
    // std::cout << " name " << name << " lastPosition[2] " << lastPosition[2] << std::endl;
    // std::cout << " name " << geo_->GetCurrentNode()->GetVolume()->GetName() << std::endl;

    int n_elements = current_material->GetNelements();
    double A = 0;
    double Z = 0;
    double w = 0;
    double material_I = 0;
    double num = 0;
    double den = 0;
    for (int i = 0; i < n_elements; i++) {
      current_material->GetElementProp(A, Z, w, i);
      int I = 0;
      if (Z < 13) {
        I = 12 * Z + 7;
      } else {
        I = 9.76 * Z + 58.8 * pow(Z, -1.19) * Z;
      }

      num += w * Z / A * log(I);
      den += w * Z / A;

    }

    material_I = exp(num / den) * 1E-6;
    
    // std::cout << "density " << density << std::endl;
    // std::cout << "pathLength " << pathLength << std::endl;
    // std::cout << "material_I " << material_I << std::endl;

    auto matZ = static_cast<double>(current_material->GetZ());
    auto matA = static_cast<double>(current_material->GetA());

    double plasma_energy = sqrt(density*fabs(matZ/matA))*28.816E-6;
    double delta = log(plasma_energy/material_I)+log(beta*gamma)-0.5;     
    double coeff = K * charge * charge * matZ / matA / (beta * beta);
    double BB_log = log(2 * m_e * beta * beta * gamma * gamma * W_max / (material_I * material_I));
    double dEdx = coeff * (0.5 * BB_log - beta * beta - 0.5*delta);

    dE += dEdx * density * pathLength;
    geo_->Step(true, true);
  }
  return dE;
}












namespace sand_reco
{

namespace kf
{

namespace utils 
{

TVector2 getBFieldPerp(const TVector3& v)
{
    return TVector2(v.Z(), v.Y());
}

double getBFieldLong(const TVector3& v)
{
    return v.X();
}

double getTanOfDipAngle(const TVector3& mom)
{
    auto perp = getBFieldPerp(mom);
    auto para = getBFieldLong(mom);
    return para/perp.Mod();
}

int getRotationVersus(int charge)
{
    return -charge;
}

int getCharge(int versus)
{
    return -versus;
}

double getDirectionAngle(const TVector2& dir)
{
    // dir.Print();
    return TMath::ATan2(dir.Y(), dir.X());
}

double getRotationAngle(const TVector2& dir, int charge)
{
    auto theta = getDirectionAngle(dir);
    // std::cout << "THETA " << theta << std::endl;
    //if (theta < 0 ) theta +=  2 * TMath::Pi();
    auto versus = getRotationVersus(charge);
    //return (theta + versus * 0.5 * TMath::Pi());
    return theta + versus * 0.5 * TMath::Pi();
}
double getRadius(double perp_mom)
{
    return SANDTrackerUtils::getRadiusInMMFromPerpMomentumInGeV(perp_mom);
}

TVector2 getCircleCenter(const TVector2& perp_mom, const TVector2& position, int charge)
{
    auto phi = getRotationAngle(perp_mom, charge);
    auto radius = getRadius(perp_mom.Mod());
    // std::cout << "weew: " << phi << " " << radius << std::endl;
    return TVector2(position.X() - radius * TMath::Cos(phi), position.Y() - radius * TMath::Sin(phi));
}

std::pair<double, double> getCircleYs(double z, double radius, const TVector2& center)
{
    auto y = TMath::Sqrt(radius * radius - (z - center.X()) * (z - center.X()));
    return {center.Y() + y, center.Y() - y};
}

double getRotationAngle(double z, double y, const TVector2& center)
{
    return TMath::ATan2(y - center.Y(), z - center.X());
}

double getDeltaPhi(double phi, double previous_phi, int versus)
{

    // da riguardare e sempilificare
    phi = phi < 0. ? 2 * TMath::Pi() + phi : phi;
    previous_phi = previous_phi < 0. ? 2 * TMath::Pi() + previous_phi : previous_phi;
    auto delta_phi = phi - previous_phi;

    auto dphi = delta_phi * versus >= 0. ? delta_phi : 2 * TMath::Pi() - delta_phi;
    // return dphi;
    return fmod(dphi, M_PI);
}

double getX(double radius, double x_0, double delta_phi, double tan_lambda, int versus_of_rot)
{
    auto charge = getCharge(versus_of_rot);
    return x_0 + radius * charge * tan_lambda * delta_phi;
}

double getY(const TVector2& center, double radius, double phi)
{   
    // std::cout << "deltaY: " << radius * TMath::Sin(phi) << std::endl;
    return center.Y() + radius * TMath::Sin(phi);
}

double getZ(const TVector2& center, double radius, double phi)
{
    // std::cout << "deltaZ: " << radius * TMath::Cos(phi) << std::endl;
    return center.X() + radius * TMath::Cos(phi);
}

TVector3 getVectorMomentum(double radius, double phi, double tan_lambda, int versus)
{
    auto perp_mon = SANDTrackerUtils::getPerpMomentumInGeVFromRadiusInMM(radius);
    auto theta = phi - versus * 0.5 * TMath::Pi();
    return TVector3(perp_mon * tan_lambda, perp_mon * TMath::Sin(theta), perp_mon * TMath::Cos(theta));
}

sand_reco::kf::StateVector getStateVector(TVector3 mom, TVector3 pos, int charge)
{
    auto perp_mom = getBFieldPerp(mom);
    auto perp_pos = getBFieldPerp(pos);
    auto radius = getRadius(perp_mom.Mod());
    // perp_mom.Print();
    auto tan_lambda = getTanOfDipAngle(mom);
    auto phi = getRotationAngle(perp_mom, charge);

    return sand_reco::kf::StateVector(pos.X(), pos.Y(), charge/radius, tan_lambda, phi);
}


Double_t makeC(Double_t z1,Double_t y1, Double_t z2,Double_t y2, Double_t z3,Double_t y3){
  //-----------------------------------------------------------------
  // Initial approzimation of the track curvature
  //-----------------------------------------------------------------
  z3 -=z1;
  z2 -=z1;
  y3 -=y1;
  y2 -=y1;
  //  
  Double_t det = z3*y2-z2*y3;
  if (TMath::Abs(det)<1e-10){
    return 100;
  }
  //
  Double_t u = 0.5* (z2*(z2-z3)+y2*(y2-y3))/det;
  Double_t z0 = z3*0.5-y3*u;
  Double_t y0 = y3*0.5+z3*u;
  Double_t c2 = 1/TMath::Sqrt(z0*z0+y0*y0);
  if (det<0) c2*=-1;
  return c2;
}


Double_t makeSnp(Double_t z1,Double_t y1, Double_t z2,Double_t y2, Double_t z3,Double_t y3){
  //-----------------------------------------------------------------
  // Initial approzimation of the track snp at position z1
  //-----------------------------------------------------------------
  z3 -=z1;
  z2 -=z1;
  y3 -=y1;
  y2 -=y1;
  //  
  Double_t det = z3*y2-z2*y3;
  if (TMath::Abs(det)<1e-10) {
    return 100;
  }
  //
  Double_t u = 0.5* (z2*(z2-z3)+y2*(y2-y3))/det;
  Double_t z0 = z3*0.5-y3*u; 
  Double_t y0 = y3*0.5+z3*u;
  Double_t c2 = 1/TMath::Sqrt(z0*z0+y0*y0);
  if (det>0) c2*=-1;
  z0*=c2;  
  return z0;
}

Double_t makeYC(Double_t z1,Double_t y1, Double_t z2,Double_t y2, Double_t z3,Double_t y3){
  //-----------------------------------------------------------------
  // Initial approzimation of the y coordinate of the center of the track circumference 
  // in the zy plane, with respects to the first point (z1,y1). Used to check consistency 
  // between points (i.e. if they are in the the same semiplane), not in the seeding itself. 
  // If the sign of yC is the same, the points are in the same semiplane.
  //-----------------------------------------------------------------
  z3 -=z1;
  z2 -=z1;
  y3 -=y1;
  y2 -=y1;
  //  
  Double_t det = z3*y2-z2*y3;
  if (TMath::Abs(det)<1e-10) {
    return 100;
  }
  //
  Double_t u = 0.5* (z2*(z2-z3)+y2*(y2-y3))/det;
  Double_t y0 = y3*0.5+z3*u;
  return y0;
}

//_____________________________________________________________________________
Double_t makeTgln(Double_t z1,Double_t y1, Double_t z2,Double_t y2,Double_t x1,Double_t x2,Double_t c){
  //-----------------------------------------------------------------
  // Initial approzimation of the tangent of the track dip angle
  //-----------------------------------------------------------------
  Double_t d  =  TMath::Sqrt((z1-z2)*(z1-z2)+(y1-y2)*(y1-y2));
  if (TMath::Abs(d*c*0.5)>1) return 0;
  Double_t   angle2    = asin(d*c*0.5);

  angle2  = (x1-x2)*c/(angle2*2.);    //dz /(R*dPhi)
  return angle2;
  //return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

sand_reco::kf::State Seed3Points(std::array<double,3> xyz0, std::array<double,3> xyz1, std::array<double,3> xyz2, double sy, double sx, float bz){
  Double_t sy2=sy*sy;
  Double_t sx2=sx*sx;
  TMatrixD param(5,1);
  TMatrixD c(5,5);
  // calculate initial param
  param[0][0]=xyz0[0];              
  param[1][0]=xyz0[1];
  param[2][0]=makeC(xyz0[2],xyz0[1],xyz1[2],xyz1[1],xyz2[2],xyz2[1]);
  auto versus = param[2][0] > 0 ? -1 : 1; 
  param[3][0]=makeTgln(xyz0[2],xyz0[1],xyz1[2],xyz1[1],xyz0[0],xyz1[0],param[2][0]);
  param[4][0]=std::asin(makeSnp(xyz0[2],xyz0[1],xyz1[2],xyz1[1],xyz2[2],xyz2[1])) + versus * 0.5 * TMath::Pi();
  sand_reco::kf::StateVector vec(param);

  //
  Double_t f20=(makeC(xyz0[2],xyz0[1]+sy,xyz1[2],xyz1[1],xyz2[2],xyz2[1])-param[2][0])/sy;
  Double_t f22=(makeC(xyz0[2],xyz0[1],xyz1[2],xyz1[1]+sy,xyz2[2],xyz2[1])-param[2][0])/sy;
  Double_t f23=(makeC(xyz0[2],xyz0[1],xyz1[2],xyz1[1],xyz2[2],xyz2[1]+sy)-param[2][0])/sy;
  //
  Double_t f40=(makeSnp(xyz0[2],xyz0[1]+sy,xyz1[2],xyz1[1],xyz2[2],xyz2[1])-param[4][0])/sy;
  Double_t f42=(makeSnp(xyz0[2],xyz0[1],xyz1[2],xyz1[1]+sy,xyz2[2],xyz2[1])-param[4][0])/sy;
  Double_t f43=(makeSnp(xyz0[2],xyz0[1],xyz1[2],xyz1[1],xyz2[2],xyz2[1]+sy)-param[4][0])/sy;
  //
  //  makeTgln(xyz0[2],xyz0[1],xyz1[2],xyz1[1],xyz0[2],xyz1[2],param[0][4]);
  Double_t f30=(makeTgln(xyz0[2],xyz0[1]+sy,xyz1[2],  xyz1[1],xyz0[0],xyz1[0],param[2][0])-param[3][0])/sy;
  Double_t f31=(makeTgln(xyz0[2],xyz0[1],xyz1[2],  xyz1[1],xyz0[0]+sx,xyz1[0],param[2][0])-param[3][0])/sx;
  Double_t f32=(makeTgln(xyz0[2],xyz0[1],xyz1[2],  xyz1[1]+sy,xyz0[0],xyz1[0],param[2][0])-param[3][0])/sy;
  Double_t f34=(makeTgln(xyz0[2],xyz0[1],xyz1[2],   xyz1[1],xyz0[0],xyz1[0]+sx,param[2][0])-param[3][0])/sx;
  c[0][0] = sx2;
  c[1][0] = 0.;        c[1][1] = sy2;
  c[2][0] = f20*sy2;   c[2][1] = 0.;       c[2][2] = f20*sy2*f20 + f22*sy2*f22 + f23*sy2*f23;
  c[3][0] = f30*sy2;   c[3][1] = f31*sx2;  c[3][2] = f30*sy2*f20 + f32*sy2*f22;             c[3][3] = f30*sy2*f30 + f31*sx2*f31 + f32*sy2*f32 + f34*sx2*f34;
  c[4][0] = f40*sy2;   c[4][1] = 0.;       c[4][2] = f40*sy2*f20 + f42*sy2*f22 + f43*sy2*f23; c[4][3] = f30*sy2*f40 + f32*sy2*f42; c[4][4] = f40*sy2*f40 + f42*sy2*f42 + f43*sy2*f43;
  

  int n = c.GetNrows(); // Assuming square matrix
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
        c[i][j] = c[j][i]; // Mirror lower triangle to upper
    }
  }

  param.Print();

  sand_reco::kf::State StateSeed(param,c);
  return StateSeed;
}

ParticleState::ParticleState(const sand_reco::kf::StateVector& vector, double z)
{
    position_ = TVector3(vector.x(), vector.y(), z);

    auto charge = vector.charge();
    auto versus = getRotationVersus(charge);
    auto radius = vector.radius();
    
    momentum_ = getVectorMomentum(radius, vector.phi(), vector.tanLambda(), versus);
}

TrajectoryParameters ParticleState::getTrajectoryParameter(int charge) const
{
    auto perp_mom = getBFieldPerp(momentum_);
    auto perp_pos = getBFieldPerp(position_);
    auto radius = getRadius(perp_mom.Mod());
    auto versus_of_rot = getRotationVersus(charge);
    auto tan_lambda = getTanOfDipAngle(momentum_);
    auto center_of_rot = getCircleCenter(perp_mom, perp_pos, charge);
    auto phi_0 = getRotationAngle(perp_mom, charge);
    auto x_0 = position_.X();

    return TrajectoryParameters(radius, versus_of_rot, tan_lambda, center_of_rot, phi_0, x_0);
}

sand_reco::kf::StateVector ParticleState::getStateVector(int charge) const
{
    auto perp_mom = getBFieldPerp(momentum_);
    auto perp_pos = getBFieldPerp(position_);
    auto radius = getRadius(perp_mom.Mod());
    // perp_mom.Print();
    auto tan_lambda = getTanOfDipAngle(momentum_);
    auto phi = getRotationAngle(perp_mom, charge);

    return sand_reco::kf::StateVector(position_.X(), position_.Y(), charge/radius, tan_lambda, phi);
}

std::pair<double, double> TrajectoryParameters::getPhiPair(double z) const {
    double dy = TMath::Sqrt(radius_ * radius_ - (center_of_rot_.X() - z) * (center_of_rot_.X() - z));
    return std::pair<double, double>(atan2(dy, (z - center_of_rot_.X())), atan2(-dy, (z - center_of_rot_.X())));
}

ParticleState TrajectoryParameters::getParticleState(double delta_phi) const
{
    auto phi = phi_0_ + delta_phi;
    auto x = getX(radius_, x_0_, delta_phi, tan_lambda_, versus_of_rot_);
    auto y = getY(center_of_rot_, radius_, phi);
    auto z = getZ(center_of_rot_, radius_, phi);
    auto p = TVector3{x,y,z};
    auto m = getVectorMomentum(radius_, phi, tan_lambda_, versus_of_rot_);

    return ParticleState(p, m);
}

std::vector<ParticleState> TrajectoryParameters::getParticleStatesFromDeltaPhi(std::vector<double> delta_phis) const
{
    std::vector<ParticleState> particle_states;
    auto total_delta_phi = 0.;
    for(auto delta_phi: delta_phis)
    {
        total_delta_phi += delta_phi;
        particle_states.emplace_back(getParticleState(total_delta_phi));
    }
    return particle_states;
}

double TrajectoryParameters::getSmallestDeltaPhi(double z, double last_phi) const
{
    auto phi_pair = getPhiPair(z);
    auto delta_phi_1 = getDeltaPhi(phi_pair.first , last_phi, versus_of_rot_);
    auto delta_phi_2 = getDeltaPhi(phi_pair.second, last_phi, versus_of_rot_);

    // debug ////////////////////////////////////////
    // std::cout << z << " " << last_phi << " " << versus_of_rot_ << " "
    //           << phi_pair.first << " " << phi_pair.second << " "
    //           << delta_phi_1 << " " << delta_phi_2 << " "
    //           << (delta_phi_1 < delta_phi_2 ? delta_phi_1 : delta_phi_2) << std::endl; 
    /////////////////////////////////////////////////

    return fabs(delta_phi_1) < fabs(delta_phi_2) ? delta_phi_1 : delta_phi_2;
}

std::vector<double> TrajectoryParameters::getDeltaPhis(std::vector<double> zs) const
{
    std::vector<double> delta_phis;

    auto last_phi = phi_0_;
    auto delta_phi = 0.;

    for(auto& z: zs)
    {
        delta_phi = getSmallestDeltaPhi(z, last_phi);

        // debug //////////////////////////////////////
        // std::cout << z << " " << last_phi << " " << delta_phi << std::endl;
        //////////////////////////////////////////////

        delta_phis.push_back(delta_phi);
        last_phi += delta_phi;
    }
    return delta_phis;
}

std::vector<ParticleState> TrajectoryParameters::getParticleStatesFromZ(std::vector<double> zs) const
{
        // get delta phi
        auto delta_phis = getDeltaPhis(zs);

        // debug /////////////////////////////
        // for(auto i = 0u; i < zs.size(); i++) std:: cout << zs.at(i) << " " << delta_phis.at(i) << std::endl;
        /////////////////////////////////////

        return getParticleStatesFromDeltaPhi(delta_phis);
}

ParticleState operator -(const ParticleState& p1, const ParticleState& p2) {
    ParticleState p;
    p.setPosition(p1.getPosition() - p2.getPosition());
    p.setMomentum(p1.getMomentum() - p2.getMomentum());
    return p;
}
bool operator ==(const ParticleState& p1, const ParticleState& p2) {
    
    return (p1.getPosition() == p2.getPosition() && p1.getMomentum() == p2.getMomentum());
}

} // namespace utils 
} // namespace kf
} // namespace sand_reco