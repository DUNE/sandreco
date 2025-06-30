#include "SANDKalmanFilter.h"
#include "SANDTrackerClusterCollection.h"
#include "SANDTrackerUtils.h"


#include <vector>

#include "TVectorD.h"

namespace sand_reco
{

namespace kf
{

TVector3 Manager::getDirectiveCosinesFromStateVector(
    const sand_reco::kf::StateVector& stateVector)
{
  auto theta = getThetaFromPhi(stateVector);
  if (theta > M_PI_2) theta -= M_PI;
  TVector3 dir(
      stateVector.tanLambda() / sqrt(1 + pow(stateVector.tanLambda(), 2)),
      sin(theta), cos(theta));
  dir *= 1. / dir.Mag();
  return dir;
}

sand_reco::kf::Measurement Manager::getMeasurementFromCluster(
    int clusterID)
{
  sand_reco::kf::Measurement measurement(2, 1);
  // auto cluster = sand_reco::kf::ClusterManager::getCluster(clusterID);
  // auto recoTrkParameter = cluster.getRecoParameters().front().trk;
  // measurement[0][0] = recoTrkParameter.m * cluster.getZ() + recoTrkParameter.q;
  // measurement[1][0] = atan(recoTrkParameter.m);
  return measurement;
}

// TMatrixD Manager::getInitialCovMatrix(
//     const sand_reco::kf::StateVector& stateVector,
//     const STTPlane::EOrientation& orientation)
// {

//   TMatrixD initialCovMatrix(5, 5);

//   auto initialPositionSigma = SANDTrackerUtils::getSigmaPositionMeasurement();
//   auto initialRadiusSigma = 5000.; /* mm */
//   auto initialAngleSigma = SANDTrackerUtils::getSigmaAngleMeasurement();

//   if (orientation == STTPlane::EOrientation::kHorizontal) {
//     initialCovMatrix[0][0] = pow(SANDTrackerUtils::getSANDInnerVolumeLength(), 2) / 12.;
//     initialCovMatrix[1][1] = initialPositionSigma * initialPositionSigma;
//     initialCovMatrix[2][2] =
//         pow(initialRadiusSigma, 2) * pow(stateVector.signedInverseRadius(), 4);
//     initialCovMatrix[3][3] =
//         pow(initialAngleSigma, 2) * pow(1 + pow(stateVector.tanLambda(), 2), 2);
//     initialCovMatrix[4][4] = pow(initialAngleSigma, 2);
//   } else {
//     initialCovMatrix[0][0] = initialPositionSigma * initialPositionSigma;
//     initialCovMatrix[1][1] = pow(SANDTrackerUtils::getSANDInnerVolumeRadius(), 2) / 12.;
//     initialCovMatrix[2][2] =
//         pow(initialRadiusSigma, 2) * pow(stateVector.signedInverseRadius(), 4);
//     initialCovMatrix[3][3] =
//         pow(initialAngleSigma, 2) * pow(1 + pow(stateVector.tanLambda(), 2), 2);
//     initialCovMatrix[4][4] = pow(initialAngleSigma, 2);
//   }

//   return initialCovMatrix;
// }

TMatrixD Manager::getPropagatorMatrix(
    const sand_reco::kf::StateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass)
{
  // PRIMO INDICE = RIGA
  TMatrixD propagatorMatrix(5, 5);
  propagatorMatrix[0][0] = dxDx(stateVector, nextPhi, dZ, dE, particle_mass);         
  propagatorMatrix[0][1] = dxDy(stateVector, nextPhi, dZ, dE, particle_mass);         
  propagatorMatrix[0][2] = dxDInvCR(stateVector, nextPhi, dZ, dE, particle_mass);     
  propagatorMatrix[0][3] = dxDTanl(stateVector, nextPhi, dZ, dE, particle_mass);      
  propagatorMatrix[0][4] = dxDPhi(stateVector, nextPhi, dZ, dE, particle_mass);       
  propagatorMatrix[1][0] = dyDx(stateVector, nextPhi, dZ, dE, particle_mass);         
  propagatorMatrix[1][1] = dyDy(stateVector, nextPhi, dZ, dE, particle_mass);         
  propagatorMatrix[1][2] = dyDInvCR(stateVector, nextPhi, dZ, dE, particle_mass);     
  propagatorMatrix[1][3] = dyDTanl(stateVector, nextPhi, dZ, dE, particle_mass);      
  propagatorMatrix[1][4] = dyDPhi(stateVector, nextPhi, dZ, dE, particle_mass);       
  propagatorMatrix[2][0] = dInvCRDx(stateVector, nextPhi, dZ, dE, particle_mass);     
  propagatorMatrix[2][1] = dInvCRDy(stateVector, nextPhi, dZ, dE, particle_mass);     
  propagatorMatrix[2][2] = dInvCRDInvCR(stateVector, nextPhi, dZ, dE, particle_mass); 
  propagatorMatrix[2][3] = dInvCRDTanl(stateVector, nextPhi, dZ, dE, particle_mass);  
  propagatorMatrix[2][4] = dInvCRDPhi(stateVector, nextPhi, dZ, dE, particle_mass);   
  propagatorMatrix[3][0] = dTanlDx(stateVector, nextPhi, dZ, dE, particle_mass);      
  propagatorMatrix[3][1] = dTanlDy(stateVector, nextPhi, dZ, dE, particle_mass);      
  propagatorMatrix[3][2] = dTanlDInvCR(stateVector, nextPhi, dZ, dE, particle_mass);
  propagatorMatrix[3][3] = dTanlDTanl(stateVector, nextPhi, dZ, dE, particle_mass);
  propagatorMatrix[3][4] = dTanlDPhi(stateVector, nextPhi, dZ, dE, particle_mass);
  propagatorMatrix[4][0] = dPhiDx(stateVector, nextPhi, dZ, dE, particle_mass);
  propagatorMatrix[4][1] = dPhiDy(stateVector, nextPhi, dZ, dE, particle_mass);
  propagatorMatrix[4][2] = dPhiDInvCR(stateVector, nextPhi, dZ, dE, particle_mass);
  propagatorMatrix[4][3] = dPhiDTanl(stateVector, nextPhi, dZ, dE, particle_mass);
  propagatorMatrix[4][4] = dPhiDPhi(stateVector, nextPhi, dZ, dE, particle_mass);
  return propagatorMatrix;
}

TMatrixD Manager::getProcessNoiseMatrix(
    const sand_reco::kf::StateVector& stateVector, double nextPhi, double dZ, double dE,
    double z, double particle_mass)
{
  TVectorD processNoiseTanlDerivative(5);
  TVectorD processNoisePhiDerivative(5);

  processNoiseTanlDerivative[0] = dxDTanl(stateVector, nextPhi, dZ, dE, particle_mass);
  processNoiseTanlDerivative[1] = dyDTanl(stateVector, nextPhi, dZ, dE, particle_mass);
  processNoiseTanlDerivative[2] = dInvCRDTanl(stateVector, nextPhi, dZ, dE, particle_mass);
  processNoiseTanlDerivative[3] = dTanlDTanl(stateVector, nextPhi, dZ, dE, particle_mass);
  processNoiseTanlDerivative[4] = dPhiDTanl(stateVector, nextPhi, dZ, dE, particle_mass);

  processNoisePhiDerivative[0] = dxDPhi(stateVector, nextPhi, dZ, dE, particle_mass);
  processNoisePhiDerivative[1] = dyDPhi(stateVector, nextPhi, dZ, dE, particle_mass);
  processNoisePhiDerivative[2] = dInvCRDPhi(stateVector, nextPhi, dZ, dE, particle_mass);
  processNoisePhiDerivative[3] = dTanlDPhi(stateVector, nextPhi, dZ, dE, particle_mass);
  processNoisePhiDerivative[4] = dPhiDPhi(stateVector, nextPhi, dZ, dE, particle_mass);

  auto momentumInMeV = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(stateVector.radius(), stateVector.tanLambda());

  auto dir = -1. * getDirectiveCosinesFromStateVector(stateVector);
  if (dir.Z() > 0) dir *= -1;
  auto pathLengthInX0 = dE!=0 ? SANDTrackerUtils::getPathLengthInX0(
      (z + dZ)*1000, stateVector.x()*1000, stateVector.y()*1000, z*1000, dir.X(), dir.Y(), dir.Z()):0;

  // MCS angle
  double radius = stateVector.radius();
  double tan    = stateVector.tanLambda();
  double mom    = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(radius, tan) / 1000;
  
  double constant = SANDTrackerUtils::getRadiusInMMToMomentumInGeVConstant() * SANDTrackerUtils::getMagneticField() / sqrt(1 + pow(tan, 2));
  
  double gamma = sqrt(mom*mom + particle_mass*particle_mass) / particle_mass;
  double beta = sqrt( 1 - pow(1/gamma, 2));
  auto sigmaMCSAngle = dE!=0 ? SANDTrackerUtils::getMCSSigmaAngleFromMomentumInMeV(
      momentumInMeV, beta, pathLengthInX0):0;
  auto sigmaMCSAngleSquared = sigmaMCSAngle * sigmaMCSAngle;

  auto factor = pow(1 + pow(stateVector.tanLambda(), 2), 2);

  TMatrixD processNoiseMatrix(5, 5);
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 5; j++)
      processNoiseMatrix[i][j] =
          sigmaMCSAngleSquared *
          (processNoisePhiDerivative[i] * processNoisePhiDerivative[j] +
           factor * processNoiseTanlDerivative[i] *
               processNoiseTanlDerivative[j]);

  return processNoiseMatrix;
}

TMatrixD Manager::getMeasurementNoiseMatrix()
{
  TMatrixD measurementNoiseMatrix(2, 2);
  measurementNoiseMatrix[0][0] = SANDTrackerUtils::getSigmaPositionMeasurement() *
                                 SANDTrackerUtils::getSigmaPositionMeasurement();

  measurementNoiseMatrix[1][1] = SANDTrackerUtils::getSigmaAngleMeasurement()*
                                 SANDTrackerUtils::getSigmaAngleMeasurement();
//   switch (orientation) {
//     case Orientation::kVertical:
//         measurementNoiseMatrix[0][0] = sigma_x * sigma_x;
//         measurementNoiseMatrix[1][1] = sigma_theta_x * sigma_theta_x;
//         break;
    
//     case Orientation::kHorizontal:
//         measurementNoiseMatrix[0][0] = sigma_y * sigma_y;
//         measurementNoiseMatrix[1][1] = sigma_theta_y * sigma_theta_y;
//         break;
// }
  return measurementNoiseMatrix;
}

TMatrixD Manager::getProjectionMatrix(
    Orientation orientation, const sand_reco::kf::StateVector& stateVector)
{
  TMatrixD projectionMatrix(2, 5);

  switch (orientation) {
    case Orientation::kVertical: {
      auto denominator =
          sin(stateVector.phi()) *
          (1 + pow(stateVector.tanLambda() / sin(stateVector.phi()), 2));
      projectionMatrix[0][0] = 1.;
      projectionMatrix[1][3] = -stateVector.charge() / denominator;
      projectionMatrix[1][4] = stateVector.charge() * stateVector.tanLambda() /
                               tan(stateVector.phi()) / denominator;
      break;
    }
    case Orientation::kHorizontal:
      projectionMatrix[0][1] = 1.;
      projectionMatrix[1][4] = 1.;
      break;
  }
  return projectionMatrix;
}

TMatrixD Manager::getKalmanGainMatrix(
    const TMatrixD& covarianceMatrix, const TMatrixD& projectionMatrix,
    const TMatrixD& measurementNoisMatrix)
{
  TMatrixD projectionMatrixTransposed(TMatrixD::kTransposed, projectionMatrix);
  TMatrixD kalmanGainInvertedDenominator(
      TMatrixD::kInverted,
      measurementNoisMatrix +
          projectionMatrix * covarianceMatrix * projectionMatrixTransposed);
  return covarianceMatrix * projectionMatrixTransposed * kalmanGainInvertedDenominator;
}

TMatrixD Manager::getAMatrix(
    const TMatrixD& covarianceMatrixFiltered,
    const TMatrixD& covarianceMatrixNextPredicted,
    const TMatrixD& propagatorMatrix)
{
  TMatrixD covarianceMatrixNextPredictedInverted(TMatrixD::kInverted,
                                                 covarianceMatrixNextPredicted);
  TMatrixD propagatorMatrixTransposed(TMatrixD::kTransposed, propagatorMatrix);
  return covarianceMatrixFiltered * propagatorMatrixTransposed *
         covarianceMatrixNextPredictedInverted;
}

double Manager::deltaRadius(
    const sand_reco::kf::StateVector& stateVector, double nextPhi, double dZ,
    double dE, double particle_mass) const
{
  // To Do: check all units
  double comp_delta_inv_radius_from_de = 0;

  double radius = stateVector.radius();
  double tan    = stateVector.tanLambda();
  double mom    = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(radius, tan) / 1000;
  // std::cout << mom << std::endl;
  
  double constant = SANDTrackerUtils::getRadiusInMMToMomentumInGeVConstant() * SANDTrackerUtils::getMagneticField() / sqrt(1 + pow(tan, 2));
  double step = 10E-8;
  double tmp_inv_r = 1. / radius;
  // Notice: this is a shortcut. Second order derivative?
  for (int k = 0; k < int(fabs(dE) / step); k++) {
    double delta_inv_r = -pow(tmp_inv_r, 3) / pow(constant, 2) * sqrt(pow(constant, 2) / pow(tmp_inv_r, 2) + pow(particle_mass, 2)) * step;
    comp_delta_inv_radius_from_de += delta_inv_r;
    tmp_inv_r += delta_inv_r;
  }
  comp_delta_inv_radius_from_de  += -pow(tmp_inv_r, 3) / pow(constant, 2) * sqrt(pow(constant, 2) / pow(tmp_inv_r, 2) + pow(particle_mass, 2)) * std::fmod(fabs(dE), step);
        
  return comp_delta_inv_radius_from_de;
}

sand_reco::kf::StateVector Manager::propagateState(
    const sand_reco::kf::StateVector& stateVector, double dZ, double dE, double particle_mass)
{
  auto nexttanLambda = stateVector.tanLambda();

  auto cosNextPhi = cos(stateVector.phi()) + dZ / stateVector.radius();
  if (cosNextPhi > 1.)
    cosNextPhi = 1.;
  else if (cosNextPhi < -1.)
    cosNextPhi = -1.;

  auto nextPhi = acos(cosNextPhi);
  if (stateVector.phi() < 0) nextPhi *= -1.;

  auto nextSignedInverseRadius =
      stateVector.signedInverseRadius() +
      (dE!=0 ? stateVector.charge() * deltaRadius(stateVector, nextPhi, dZ, dE, particle_mass) : 0);

  auto cosNextPhi_corr = cos(stateVector.phi()) + dZ * stateVector.charge() * nextSignedInverseRadius;
  if (cosNextPhi_corr > 1.)
    cosNextPhi_corr = 1.;
  else if (cosNextPhi_corr < -1.)
    cosNextPhi_corr = -1.;

  auto nextPhi_corr = acos(cosNextPhi_corr);
  if (stateVector.phi() < 0) nextPhi_corr *= -1.;

  nextPhi = 0.5 * (nextPhi + nextPhi_corr);

  auto nextX = stateVector.x() + stateVector.chargedRadius() *
                                     stateVector.tanLambda() *
                                     (nextPhi - stateVector.phi());

  auto nextY = stateVector.y() +
               stateVector.radius() * (sin(nextPhi) - sin(stateVector.phi()));
  return sand_reco::kf::StateVector(nextX, nextY, nextSignedInverseRadius, nexttanLambda,
                          nextPhi);
}

sand_reco::kf::StateCovarianceMatrix Manager::propagateCovMatrix(
    const TMatrixD& covarianceMatrix, const TMatrixD& propagatorMatrix,
    const TMatrixD& processNoiseMatrix)
{
  TMatrixD propagatorMatrixTransported(TMatrixD::kTransposed, propagatorMatrix);
  
  return propagatorMatrix * covarianceMatrix * propagatorMatrixTransported +
         processNoiseMatrix;
}

sand_reco::kf::Measurement Manager::getPrediction(
    Orientation orientation, const sand_reco::kf::StateVector& stateVector)
{
  sand_reco::kf::Measurement projector(2, 1);
  if (orientation == Orientation::kHorizontal) {
    projector[0][0] = stateVector.y();
    projector[1][0] =
        stateVector.phi() + stateVector.charge() * 0.5 * TMath::Pi();
  } else {
    projector[0][0] = stateVector.x();
    projector[1][0] = -stateVector.charge() *
                      atan(stateVector.tanLambda() / sin(stateVector.phi()));
  }
  return projector;
}

sand_reco::kf::StateVector Manager::filterState(
    const sand_reco::kf::StateVector& stateVector, const TMatrixD& kalmanGainMatrix,
    const sand_reco::kf::Measurement& observed, const sand_reco::kf::Measurement& predicted)
{
  return sand_reco::kf::StateVector(stateVector() +
                          kalmanGainMatrix * (observed - predicted));
}

sand_reco::kf::StateCovarianceMatrix Manager::filterCovMatrix(
    const TMatrixD& covarianceMatrix, const TMatrixD& projectionMatrix,
    const TMatrixD& measurementNoiseMatrix)
{
  TMatrixD covarianceMatrixInverted(TMatrixD::kInverted, covarianceMatrix);
  TMatrixD projectionMatrixTransposed(TMatrixD::kTransposed, projectionMatrix);
  TMatrixD measurementNoiseMatrixInverted(TMatrixD::kInverted,
                                          measurementNoiseMatrix);
  auto nextCovarianceMatrixInverted =
      covarianceMatrixInverted + projectionMatrixTransposed *
                                     measurementNoiseMatrixInverted *
                                     projectionMatrix;
  return TMatrixD(TMatrixD::kInverted, nextCovarianceMatrixInverted);
}

sand_reco::kf::StateVector Manager::smoothState(
    const sand_reco::kf::StateVector& stateVectorFiltered,
    const sand_reco::kf::StateVector& stateVectorPreviousSmoothed,
    const sand_reco::kf::StateVector& stateVectorPreviousPredicted,
    const TMatrixD& theAMatrix)
{
  return sand_reco::kf::StateVector(stateVectorFiltered() +
                          theAMatrix * (stateVectorPreviousSmoothed() -
                                        stateVectorPreviousPredicted()));
}

sand_reco::kf::StateCovarianceMatrix Manager::smoothCovMatrix(
    const TMatrixD& covarianceMatrixFiltered,
    const TMatrixD& covarianceMatrixPreviousSmoothed,
    const TMatrixD& covarianceMatrixPreviousPredicted,
    const TMatrixD& theAMatrix)
{
  TMatrixD theAMatrixTransposed(TMatrixD::kTransposed, theAMatrix);
  return covarianceMatrixFiltered + theAMatrix *
                                        (covarianceMatrixPreviousSmoothed -
                                         covarianceMatrixPreviousPredicted) *
                                        theAMatrixTransposed;
}

void Manager::propagate(double& dE,
                                         double& dZ,
                                         double& beta)
{

  auto currentState = this_track_.getStep(current_step_);
  auto currentStage =
      currentState.getStage(sand_reco::kf::TrackStep::TrackStateStage::kFiltering);

  auto currentStateVector = currentStage.getStateVector();
  auto predictedStateVector = propagateState(currentStateVector, dZ, dE, particleInfo_.mass);
  auto nextPhi = predictedStateVector.phi();
  
  auto processNoiseMatrix =
      getProcessNoiseMatrix(currentStateVector, nextPhi, dZ, dE, current_z_, particleInfo_.mass);
  auto propagatorMatrix =
      getPropagatorMatrix(currentStateVector, nextPhi, dZ, dE, particleInfo_.mass);
  auto predictedCovMatrix = propagateCovMatrix(
      currentStage.getStateCovMatrix(), propagatorMatrix, processNoiseMatrix);

  sand_reco::kf::TrackStep predictedTrackStep;
  predictedTrackStep.setStage(
      sand_reco::kf::TrackStep::TrackStateStage::kPrediction,
      sand_reco::kf::State(predictedStateVector, predictedCovMatrix));
  predictedTrackStep.setPropagatorMatrix(propagatorMatrix);
  
  this_track_.addStep(predictedTrackStep);

  current_step_++;
  current_stage_ = sand_reco::kf::TrackStep::TrackStateStage::kPrediction;
}

double Manager::evalChi2(
    const sand_reco::kf::Measurement& observation, const sand_reco::kf::Measurement& prediction,
    const TMatrixD& measurementNoiseMatrix)
{
  auto residualVector = observation - prediction;
  TMatrixD residualVectorTransposed(TMatrixD::kTransposed, residualVector);
  TMatrixD measurementNoiseMatrixInverted(TMatrixD::kInverted,
    measurementNoiseMatrix);
    auto chi2Matrix = residualVectorTransposed * measurementNoiseMatrixInverted *
    residualVector;
  return chi2Matrix[0][0];
}

sand_reco::kf::Measurement Manager::getMeasurementFromTracklet(const TVectorD& tracklet)
{
  sand_reco::kf::Measurement measurement(2, 1);
  // To Do: vertical and horizontal are outdated and confusing. Replace with something more meaningful.
  // Notice: vertical planes means horizontal measurements and the opposite
  if (current_orientation_ == Orientation::kVertical) {
    // To Do: Check units!
    measurement[0][0] = tracklet[0] / 1000.;
    measurement[1][0] = tracklet[2];
  } else {
    measurement[0][0] = tracklet[1] / 1000.;
    measurement[1][0] = tracklet[3];
  }

  return measurement;
}

int Manager::findBestMatch(double& nextZ, const sand_reco::kf::Measurement& prediction,
    const TMatrixD& Sk)
{
  double best_chi = 1E9;
  auto& next_tracklets = z_to_tracklets_->at(nextZ);
  auto best_tracklet_index = -1;

  // std::cout << "next_tracklets.size(): " << next_tracklets.size() << std::endl;
  
  for (int i = 0; i < (int)next_tracklets.size(); i++) {
    sand_reco::kf::Measurement measurement = getMeasurementFromTracklet(next_tracklets[i]);
    
    auto chi2 = evalChi2(measurement, prediction, Sk);
    // measurement.Print();
    // prediction.Print();
    // std::cout << "chi2: " << chi2 << std::endl;
    if (chi2 < best_chi) {
      best_chi = chi2;
      best_tracklet_index = i;
    }
  }
  if (best_chi < 3) {
    // std::cout << "best_chi: " << best_chi << std::endl;
    return best_tracklet_index;
  } else {
    return -1;
  }
}

void Manager::setNextOrientation()
{
  if (current_orientation_ == Orientation::kVertical) {
    current_orientation_ = Orientation::kHorizontal;
  } else {
    current_orientation_ = Orientation::kVertical;
  }
}

void Manager::EvaluateInnovation(const SANDKFMeasurement& measurement, 
                                                                const SANDKFMeasurement& prediction,
                                                                const TMatrixD&  Sk)
{       
  auto innovation = measurement - prediction;
  std::vector<double> g(innovation.GetNrows());

  for (int i = 0; i < innovation.GetNrows(); i++) {  
    double r = innovation[i][0];
    double C = Sk[i][i];
    g[i] = r/sqrt(C);
  }
  this_track_.setInnovation(current_step_, g);

}


void Manager::filter(const sand_reco::kf::Measurement& measurement,
  const sand_reco::kf::Measurement& prediction)
{

  // measurement.Print();
  // prediction.Print();

  auto currentState = this_track_.getStep(current_step_);
  auto predictedStage =
      currentState.getStage(sand_reco::kf::TrackStep::TrackStateStage::kPrediction);
  auto predictedStateVector = predictedStage.getStateVector();
  auto predictedCovMatrix = predictedStage.getStateCovMatrix();

  auto measurementNoiseMatrix = getMeasurementNoiseMatrix();
  auto projectionMatrix = getProjectionMatrix(current_orientation_, predictedStateVector);
  auto kalmanGainMatrix = getKalmanGainMatrix(
      predictedCovMatrix, projectionMatrix, measurementNoiseMatrix);
  auto filteredStateVector = filterState(predictedStateVector, kalmanGainMatrix,
                                         measurement, prediction);
  auto filteredCovMatrix = filterCovMatrix(predictedCovMatrix, projectionMatrix,
                                           measurementNoiseMatrix);

  this_track_.setStage(current_step_,
                      sand_reco::kf::TrackStep::TrackStateStage::kFiltering,
                      sand_reco::kf::State(filteredStateVector, filteredCovMatrix));
  current_stage_ = sand_reco::kf::TrackStep::TrackStateStage::kFiltering;

  setNextOrientation();

  auto predictionStateCovMatrix = this_track_.getStep(current_step_).getStage(sand_reco::kf::TrackStep::TrackStateStage::kPrediction)
            .getStateCovMatrix();
  TMatrixD projectionMatrixTransposed(TMatrixD::kTransposed,
                                          projectionMatrix);
  TMatrixD Sk = measurementNoiseMatrix + projectionMatrix *
                                              predictionStateCovMatrix *
                                              projectionMatrixTransposed;
  EvaluateInnovation(measurement, prediction, Sk);
}
 

void Manager::smooth()
{

  auto currentState = this_track_.getStep(current_step_);
  auto filteredState =
      currentState.getStage(sand_reco::kf::TrackStep::TrackStateStage::kFiltering);
  auto filteredStateVector = filteredState.getStateVector();
  auto filteredCovMatrix = filteredState.getStateCovMatrix();

  if (current_step_ == int(this_track_.getSteps().size()) - 1) {
    this_track_.setStage(current_step_,
                        sand_reco::kf::TrackStep::TrackStateStage::kSmoothing,
                        sand_reco::kf::State(filteredStateVector, filteredCovMatrix));
  } else {
    // previous state
    auto previousState = this_track_.getStep(current_step_ + 1);

    // previous smoothed
    auto previousSmoothedState = previousState.getStage(
        sand_reco::kf::TrackStep::TrackStateStage::kSmoothing);
    auto previousSmoothedStateVector = previousSmoothedState.getStateVector();
    auto previousSmoothedCovMatrix = previousSmoothedState.getStateCovMatrix();
    
    // previous predicted
    auto previousPredictedState = previousState.getStage(
        sand_reco::kf::TrackStep::TrackStateStage::kPrediction);
    auto previousPredictedStateVector = previousPredictedState.getStateVector();
    auto previousPredictedCovMatrix =
        previousPredictedState.getStateCovMatrix();

    // current predicted
    auto predictedState = currentState.getStage(
        sand_reco::kf::TrackStep::TrackStateStage::kPrediction);
    auto predictedStateVector = predictedState.getStateVector();
    auto predictedCovMatrix = predictedState.getStateCovMatrix();

    auto nextPhi = previousPredictedStateVector.phi();
    auto propagatorMatrix = currentState.getPropagatorMatrix();

    auto theAMatrix = getAMatrix(filteredCovMatrix, previousPredictedCovMatrix,
                                 propagatorMatrix);

    auto smoothedStateVector =
        smoothState(filteredStateVector, previousSmoothedStateVector,
                    previousPredictedStateVector, theAMatrix);
    auto smoothedCovMatrix =
        smoothCovMatrix(filteredCovMatrix, previousSmoothedCovMatrix,
                        previousPredictedCovMatrix, theAMatrix);

    this_track_.setStage(current_step_,
                        sand_reco::kf::TrackStep::TrackStateStage::kSmoothing,
                        sand_reco::kf::State(smoothedStateVector, smoothedCovMatrix));
    // currentState.setStage(sand_reco::kf::TrackStep::TrackStateStage::kSmoothing,
    // sand_reco::kf::State(smoothedStateVector, smoothedCovMatrix));

  }
  current_stage_ = sand_reco::kf::TrackStep::TrackStateStage::kFiltering;
  current_step_--;
}

void Manager::initFromMC(TrackletMap* z_to_tracklets, const SParticleInfo& particleInfo)
{

  TMatrixD initial_cov_matrix(5, 5);
  initial_cov_matrix[0][0] = 5*pow(200E-6, 2);
  initial_cov_matrix[1][1] = 5*pow(200E-6, 2);
  initial_cov_matrix[2][2] = 5*pow(0.1, 2);
  initial_cov_matrix[3][3] = 5*pow(0.1, 2);
  initial_cov_matrix[4][4] = 5*pow(0.1, 2);

  sand_reco::kf::StateVector initial_state_vector = sand_reco::kf::utils::getStateVector(particleInfo.mom * 1E-3,  // GeV
                                                                      particleInfo.pos * 1E-3,  // m
                                                                      particleInfo.charge);

  sand_reco::kf::TrackStep trackStep;
  trackStep.setStage(sand_reco::kf::TrackStep::TrackStateStage::kPrediction,
                      sand_reco::kf::State(initial_state_vector, initial_cov_matrix));
  trackStep.setStage(sand_reco::kf::TrackStep::TrackStateStage::kFiltering,
                      sand_reco::kf::State(initial_state_vector, initial_cov_matrix));


  trackStep.setPropagatorMatrix(initial_cov_matrix);
  

  TMatrixD vectorMC(5,1);
  vectorMC[0][0] = initial_state_vector.x();
  vectorMC[1][0] = initial_state_vector.y();
  vectorMC[2][0] = initial_state_vector.signedInverseRadius();
  vectorMC[3][0] = initial_state_vector.tanLambda();
  vectorMC[4][0] = initial_state_vector.phi();

  // vectorMC.Print();

  particleInfo_        = particleInfo;
  z_to_tracklets_      = z_to_tracklets;
  current_stage_       = sand_reco::kf::TrackStep::TrackStateStage::kFiltering;
  current_step_        = 0u;
  current_z_           = particleInfo.pos.Z(); //Notice: UNITS!!  mm, why?
  current_orientation_ = Orientation::kHorizontal;
  
  this_track_.Clear();
  this_track_.addStep(trackStep);

  this_track_.setZ(current_step_, particleInfo.pos.Z());
  this_track_.setX(current_step_, particleInfo.pos.X());
  this_track_.setY(current_step_, particleInfo.pos.Y());

}


double Manager::findClosestNonEmptyKey(const TrackletMap& myMap, double target) {
  if (myMap.empty()) {
      throw std::runtime_error("Map is empty!");
  }

  // Find the first key that is >= target
  auto it = myMap.lower_bound(target);

  // Ensure we get a key <= target
  if (it == myMap.end() || it->first > target) {
      if (it == myMap.begin()) {
          throw std::runtime_error("No valid key found that is <= target with a non-empty vector!");
      }
      --it; // Move back to ensure key ≤ target
  }

  // Move back until we find a non-empty vector or reach the beginning
  while (it->first <= target) {
      if (!it->second.empty()) {
          return it->first; // Found a valid key
      }
      if (it == myMap.begin()) {
          break; // Stop if we are at the beginning
      }
      --it;
  }

  throw std::runtime_error("No valid key found that is <= target with a non-empty vector!");
}

TrackletMap Manager::FindSeedPoints_MCstart(TrackletMap* z_to_tracklets, const SParticleInfo& particloInfo, int maxSteps){
  
  double closest_key;  // Variable to store the closest key
  TrackletMap tracklet_map; // Output tracklet map
  
  closest_key = findClosestNonEmptyKey(*z_to_tracklets, particloInfo.pos.Z());
  
  // std::cout << "Closest key in z: " << closest_key << std::endl;
  // std::cout << "Number of traclets in closest key: " << z_to_tracklets->at(closest_key).size() << std::endl;


  // Add the last tracklet to the map
  std::vector<TVectorD> last_traclet;
  last_traclet.push_back(z_to_tracklets->at(closest_key)[0]);
  tracklet_map[closest_key] = last_traclet;
  // std::cout << "First tracklet found (x,y,z) : (" << z_to_tracklets->at(closest_key)[0][0] << " , "; 
  // std::cout << z_to_tracklets->at(closest_key)[0][1] << " , "<<closest_key<<" )" << std::endl;


  // Find tracklet maxStep steps below or as close as possible without an empty vector
  auto trl = z_to_tracklets->find(closest_key);

  int realStep = 0;

  for (int step = 0; step < maxSteps; ++step) {
    if (trl != z_to_tracklets->begin()) {
        --trl;
        realStep++;
    } else {
        break; // Stop if we reach the beginning
    }
  }

  auto finalStep = realStep;
  for (int step = 0; step < realStep; ++step) {
    if (!trl->second.empty()) {
        // std::cout << "Found valid far key: " << trl->first << " for now choosing first value (x,y,z): (";
        // std::cout << trl->second[0][0] << " , " << trl->second[0][1] << " , " << trl->first << " )" << std::endl;
        std::vector<TVectorD> first_tracklet;
        first_tracklet.push_back(trl->second[0]);
        tracklet_map[trl->first] = first_tracklet;
        break; // Stop searching after finding a valid key
    } else if (trl->first != closest_key){
        ++trl; 
        finalStep--;
    } else{
      // std::cout << "No valid point found within the given range." << std::endl;
      return tracklet_map; // Return empty map if no valid point is found
    }
  }

  for (int step = 0; step < finalStep; ++step) {
    if (!trl->second.empty() && step>= finalStep/2) {
        // std::cout << "Found valid middle key: " << trl->first << " for now choosing first value (x,y,z): (";
        // std::cout << trl->second[0][0] << " , " << trl->second[0][1] << " , " << trl->first << " )" << std::endl;
        std::vector<TVectorD> first_tracklet;
        first_tracklet.push_back(trl->second[0]);
        tracklet_map[trl->first] = first_tracklet;
        break; // Stop searching after finding a valid key
    } else if (trl->first != closest_key){
        ++trl; // Stop if we reach the beginning
    } else{
      // std::cout << "No valid middle point found." << std::endl;
      return tracklet_map; // Return empty map if no valid point is found
    }
  }

  return tracklet_map;


}

// To Do: implement a seeding algorithm
void Manager::initFromSeed(TrackletMap* three_tracklets, TrackletMap* z_to_tracklets, const SParticleInfo& particleInfo, double sx, double sy)
{
  try{
    if(three_tracklets->size()==3){}
    else{throw three_tracklets->size();}
  }
  catch (size_t size_tr) {
    std::cerr << "Error, number of tracklets for Seeding is 3, but given: " << size_tr << std::endl;
    return;
  }

  auto it = three_tracklets->begin();
  std::vector<std::array<double,3>> xyz;
  for (int i = 0; i < 3; ++i) {
    std::array<double,3> xyzi = {1E-3*it->second[0][0], 1E-3*it->second[0][1], 1E-3*it->first};
    xyz.push_back(xyzi);
    ++it;
  }

  auto state = sand_reco::kf::utils::Seed3Points(xyz[2],xyz[1],xyz[0],sx,sy);

  sand_reco::kf::TrackStep trackStep;
  trackStep.setStage(sand_reco::kf::TrackStep::TrackStateStage::kPrediction,
                      sand_reco::kf::State(state.getStateVector(), state.getStateCovMatrix()));
  trackStep.setStage(sand_reco::kf::TrackStep::TrackStateStage::kFiltering,
                      sand_reco::kf::State(state.getStateVector(), state.getStateCovMatrix()));


  trackStep.setPropagatorMatrix(state.getStateCovMatrix());
  
  this_track_.addStep(trackStep);

  particleInfo_       = particleInfo;
  z_to_tracklets_     = z_to_tracklets;
  current_stage_       = sand_reco::kf::TrackStep::TrackStateStage::kFiltering;
  current_step_        = 0u;
  current_z_           = xyz[2][2]; //Notice: UNITS!!  mm, why?
  current_orientation_ = Orientation::kVertical;
}


void Manager::run()
{

  // criterio per quando fermare la ricerca
  int stepLength = 1;
  if (z_to_tracklets_->lower_bound(current_z_) == z_to_tracklets_->begin()) {
    this_track_.removeLastStep();
    current_step_--;
    return;
  }


  // Notice: if currentZ is not in the map, the second condition is always true
  bool in_range = true;
  while (stepLength < 10 && std::distance(z_to_tracklets_->begin(), z_to_tracklets_->lower_bound(current_z_)) >= stepLength) {
    // 1- propagate to [currentPlaneID - step]
    auto it = (z_to_tracklets_->lower_bound(current_z_));
    for(int i= 0; i < stepLength; i++){
      if(it == z_to_tracklets_->begin()){
        in_range = false;
        break;
      }
      --it;
    }
    
    if (!in_range) {
      break;
    }
    auto nextZ = std::prev(z_to_tracklets_->lower_bound(current_z_), stepLength)->first;   

    auto currentStep = this_track_.getStep(current_step_);
    auto filteredStateVector =
        currentStep.getStage(sand_reco::kf::TrackStep::TrackStateStage::kFiltering)
            .getStateVector();

    auto dir = -1. * getDirectiveCosinesFromStateVector(filteredStateVector);

    // To Do: check if this is still valid and add a real fix if needed
    if (dir.Z() > 0) {
      dir *= -1;
    }

    auto current_mom = SANDTrackerUtils::getMomentumInMeVFromRadiusInMM(
                              filteredStateVector.radius(),
                              filteredStateVector.tanLambda()) / 1000;
    double gamma = sqrt(current_mom * current_mom + particleInfo_.mass * particleInfo_.mass) /
                    particleInfo_.mass;
    double beta = sqrt(1 - pow(1 / gamma, 2));

    // To Do: check all units
    auto dE = SANDTrackerUtils::getDE(
                        nextZ, 
                        1000 * filteredStateVector.x(), 1000 * filteredStateVector.y(), current_z_, 
                        dir.X(), dir.Y(), dir.Z(),
                        beta, particleInfo_.mass, particleInfo_.charge) / 1000;

    double dZ = (nextZ - current_z_) / 1000;

    propagate(dE, dZ, beta);

    // 2- Search best match
    auto predictionStateVector = this_track_.getStep(current_step_)
            .getStage(sand_reco::kf::TrackStep::TrackStateStage::kPrediction).getStateVector();
    auto predictionStateCovMatrix = this_track_.getStep(current_step_).getStage(sand_reco::kf::TrackStep::TrackStateStage::kPrediction)
            .getStateCovMatrix();
    auto prediction = getPrediction(current_orientation_, predictionStateVector);

    auto measurementNoiseMatrix = getMeasurementNoiseMatrix();
    auto projectionMatrix = getProjectionMatrix(current_orientation_, 
                                                predictionStateVector);
    TMatrixD projectionMatrixTransposed(TMatrixD::kTransposed,
                                          projectionMatrix);
    auto Sk = measurementNoiseMatrix + projectionMatrix * predictionStateCovMatrix *
                                              projectionMatrixTransposed;

    // measurementNoiseMatrix.Print();
    // projectionMatrix.Print();
    // predictionStateCovMatrix.Print();
    int tracklet_index = findBestMatch(nextZ, prediction, measurementNoiseMatrix);

    // // 3- If it is found: step = 1
    // //    else step++

    if (tracklet_index != -1) {
      stepLength = 1;
      auto measurement = getMeasurementFromTracklet(z_to_tracklets_->at(nextZ)[tracklet_index]);
      this_track_.setZ(current_step_, nextZ);
      this_track_.setX(current_step_, z_to_tracklets_->at(nextZ)[tracklet_index][0]);
      this_track_.setY(current_step_, z_to_tracklets_->at(nextZ)[tracklet_index][1]);

      filter(measurement, prediction);

      current_z_ = nextZ;
    } else {
      stepLength++;
      this_track_.removeLastStep();
      current_step_--;
      current_stage_ = sand_reco::kf::TrackStep::TrackStateStage::kFiltering;
    }
  }

  while (current_step_ >= 0) smooth();
}

} // namespace kf
} // namespace sand_reco