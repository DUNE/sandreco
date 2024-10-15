#include "STTKFKalmanFilter.h"
#include "SANDTrackerClusterCollection.h"
#include "SANDTrackerUtils.h"

#include <vector>

#include "TVectorD.h"

TVector3 STTKFKalmanFilterManager::GetDirectiveCosinesFromStateVector(
    const STTKFStateVector& stateVector)
{
  auto theta = GetThetaFromPhi(stateVector);
  if (theta > M_PI_2) theta -= M_PI;
  TVector3 dir(
      stateVector.TanLambda() / sqrt(1 + pow(stateVector.TanLambda(), 2)),
      sin(theta), cos(theta));
  dir *= 1. / dir.Mag();
  return dir;
}

STTKFMeasurement STTKFKalmanFilterManager::GetMeasurementFromCluster(
    int clusterID)
{
  STTKFMeasurement measurement(2, 1);
  // auto cluster = STTKFClusterManager::GetCluster(clusterID);
  // auto recoTrkParameter = cluster.GetRecoParameters().front().trk;
  // measurement[0][0] = recoTrkParameter.m * cluster.GetZ() + recoTrkParameter.q;
  // measurement[1][0] = atan(recoTrkParameter.m);
  return measurement;
}

// TMatrixD STTKFKalmanFilterManager::GetInitialCovMatrix(
//     const STTKFStateVector& stateVector,
//     const STTPlane::EOrientation& orientation)
// {

//   TMatrixD initialCovMatrix(5, 5);

//   auto initialPositionSigma = SANDTrackerUtils::GetSigmaPositionMeasurement();
//   auto initialRadiusSigma = 5000.; /* mm */
//   auto initialAngleSigma = SANDTrackerUtils::GetSigmaAngleMeasurement();

//   if (orientation == STTPlane::EOrientation::kHorizontal) {
//     initialCovMatrix[0][0] = pow(SANDTrackerUtils::GetSANDInnerVolumeLength(), 2) / 12.;
//     initialCovMatrix[1][1] = initialPositionSigma * initialPositionSigma;
//     initialCovMatrix[2][2] =
//         pow(initialRadiusSigma, 2) * pow(stateVector.SignedInverseRadius(), 4);
//     initialCovMatrix[3][3] =
//         pow(initialAngleSigma, 2) * pow(1 + pow(stateVector.TanLambda(), 2), 2);
//     initialCovMatrix[4][4] = pow(initialAngleSigma, 2);
//   } else {
//     initialCovMatrix[0][0] = initialPositionSigma * initialPositionSigma;
//     initialCovMatrix[1][1] = pow(SANDTrackerUtils::GetSANDInnerVolumeRadius(), 2) / 12.;
//     initialCovMatrix[2][2] =
//         pow(initialRadiusSigma, 2) * pow(stateVector.SignedInverseRadius(), 4);
//     initialCovMatrix[3][3] =
//         pow(initialAngleSigma, 2) * pow(1 + pow(stateVector.TanLambda(), 2), 2);
//     initialCovMatrix[4][4] = pow(initialAngleSigma, 2);
//   }

//   return initialCovMatrix;
// }

TMatrixD STTKFKalmanFilterManager::GetPropagatorMatrix(
    const STTKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass)
{
  // PRIMO INDICE = RIGA
  TMatrixD propagatorMatrix(5, 5);
  propagatorMatrix[0][0] = DxDx(stateVector, nextPhi, dZ, dE, particle_mass);         
  propagatorMatrix[0][1] = DxDy(stateVector, nextPhi, dZ, dE, particle_mass);         
  propagatorMatrix[0][2] = DxDInvCR(stateVector, nextPhi, dZ, dE, particle_mass);     
  propagatorMatrix[0][3] = DxDTanl(stateVector, nextPhi, dZ, dE, particle_mass);      
  propagatorMatrix[0][4] = DxDPhi(stateVector, nextPhi, dZ, dE, particle_mass);       
  propagatorMatrix[1][0] = DyDx(stateVector, nextPhi, dZ, dE, particle_mass);         
  propagatorMatrix[1][1] = DyDy(stateVector, nextPhi, dZ, dE, particle_mass);         
  propagatorMatrix[1][2] = DyDInvCR(stateVector, nextPhi, dZ, dE, particle_mass);     
  propagatorMatrix[1][3] = DyDTanl(stateVector, nextPhi, dZ, dE, particle_mass);      
  propagatorMatrix[1][4] = DyDPhi(stateVector, nextPhi, dZ, dE, particle_mass);       
  propagatorMatrix[2][0] = DInvCRDx(stateVector, nextPhi, dZ, dE, particle_mass);     
  propagatorMatrix[2][1] = DInvCRDy(stateVector, nextPhi, dZ, dE, particle_mass);     
  propagatorMatrix[2][2] = DInvCRDInvCR(stateVector, nextPhi, dZ, dE, particle_mass); 
  propagatorMatrix[2][3] = DInvCRDTanl(stateVector, nextPhi, dZ, dE, particle_mass);  
  propagatorMatrix[2][4] = DInvCRDPhi(stateVector, nextPhi, dZ, dE, particle_mass);   
  propagatorMatrix[3][0] = DTanlDx(stateVector, nextPhi, dZ, dE, particle_mass);      
  propagatorMatrix[3][1] = DTanlDy(stateVector, nextPhi, dZ, dE, particle_mass);      
  propagatorMatrix[3][2] = DTanlDInvCR(stateVector, nextPhi, dZ, dE, particle_mass);
  propagatorMatrix[3][3] = DTanlDTanl(stateVector, nextPhi, dZ, dE, particle_mass);
  propagatorMatrix[3][4] = DTanlDPhi(stateVector, nextPhi, dZ, dE, particle_mass);
  propagatorMatrix[4][0] = DPhiDx(stateVector, nextPhi, dZ, dE, particle_mass);
  propagatorMatrix[4][1] = DPhiDy(stateVector, nextPhi, dZ, dE, particle_mass);
  propagatorMatrix[4][2] = DPhiDInvCR(stateVector, nextPhi, dZ, dE, particle_mass);
  propagatorMatrix[4][3] = DPhiDTanl(stateVector, nextPhi, dZ, dE, particle_mass);
  propagatorMatrix[4][4] = DPhiDPhi(stateVector, nextPhi, dZ, dE, particle_mass);
  return propagatorMatrix;
}

TMatrixD STTKFKalmanFilterManager::GetProcessNoiseMatrix(
    const STTKFStateVector& stateVector, double nextPhi, double dZ, double dE,
    double z, double particle_mass)
{
  TVectorD processNoiseTanlDerivative(5);
  TVectorD processNoisePhiDerivative(5);

  processNoiseTanlDerivative[0] = DxDTanl(stateVector, nextPhi, dZ, dE, particle_mass);
  processNoiseTanlDerivative[1] = DyDTanl(stateVector, nextPhi, dZ, dE, particle_mass);
  processNoiseTanlDerivative[2] = DInvCRDTanl(stateVector, nextPhi, dZ, dE, particle_mass);
  processNoiseTanlDerivative[3] = DTanlDTanl(stateVector, nextPhi, dZ, dE, particle_mass);
  processNoiseTanlDerivative[4] = DPhiDTanl(stateVector, nextPhi, dZ, dE, particle_mass);

  processNoisePhiDerivative[0] = DxDPhi(stateVector, nextPhi, dZ, dE, particle_mass);
  processNoisePhiDerivative[1] = DyDPhi(stateVector, nextPhi, dZ, dE, particle_mass);
  processNoisePhiDerivative[2] = DInvCRDPhi(stateVector, nextPhi, dZ, dE, particle_mass);
  processNoisePhiDerivative[3] = DTanlDPhi(stateVector, nextPhi, dZ, dE, particle_mass);
  processNoisePhiDerivative[4] = DPhiDPhi(stateVector, nextPhi, dZ, dE, particle_mass);

  auto momentumInMeV = SANDTrackerUtils::GetMomentumInMeVFromRadiusInMM(stateVector.Radius(), stateVector.TanLambda());

  auto dir = -1. * GetDirectiveCosinesFromStateVector(stateVector);
  if (dir.Z() > 0) dir *= -1;
  auto pathLengthInX0 = STTKFGeoManager::GetPathLengthInX0(
      (z + dZ)*1000, stateVector.X()*1000, stateVector.Y()*1000, z*1000, dir.X(), dir.Y(), dir.Z());

  // MCS angle
  double radius = stateVector.Radius();
  double tan    = stateVector.TanLambda();
  double mom    = SANDTrackerUtils::GetMomentumInMeVFromRadiusInMM(radius, tan) / 1000;
  
  double constant = SANDTrackerUtils::GetRadiusInMMToMomentumInGeVConstant() * SANDTrackerUtils::GetMagneticField() / sqrt(1 + pow(tan, 2));
  
  double gamma = sqrt(mom*mom + particle_mass*particle_mass) / particle_mass;
  double beta = sqrt( 1 - pow(1/gamma, 2));
  auto sigmaMCSAngle = SANDTrackerUtils::GetMCSSigmaAngleFromMomentumInMeV(
      momentumInMeV, beta, pathLengthInX0);
  auto sigmaMCSAngleSquared = sigmaMCSAngle * sigmaMCSAngle;

  auto factor = pow(1 + pow(stateVector.TanLambda(), 2), 2);

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

TMatrixD STTKFKalmanFilterManager::GetMeasurementNoiseMatrix()
{
  TMatrixD measurementNoiseMatrix(2, 2);
  measurementNoiseMatrix[0][0] = SANDTrackerUtils::GetSigmaPositionMeasurement() *
                                 SANDTrackerUtils::GetSigmaPositionMeasurement();
  measurementNoiseMatrix[1][1] = SANDTrackerUtils::GetSigmaAngleMeasurement() *
                                 SANDTrackerUtils::GetSigmaAngleMeasurement();
  return measurementNoiseMatrix;
}

TMatrixD STTKFKalmanFilterManager::GetProjectionMatrix(
    Orientation orientation, const STTKFStateVector& stateVector)
{
  TMatrixD projectionMatrix(2, 5);

  switch (orientation) {
    case Orientation::kVertical: {
      auto denominator =
          sin(stateVector.Phi()) *
          (1 + pow(stateVector.TanLambda() / sin(stateVector.Phi()), 2));
      projectionMatrix[0][0] = 1.;
      projectionMatrix[1][3] = -stateVector.Charge() / denominator;
      projectionMatrix[1][4] = stateVector.Charge() * stateVector.TanLambda() /
                               tan(stateVector.Phi()) / denominator;
      break;
    }
    case Orientation::kHorizontal:
      projectionMatrix[0][1] = 1.;
      projectionMatrix[1][4] = 1.;
      break;
  }
  return projectionMatrix;
}

TMatrixD STTKFKalmanFilterManager::GetKalmanGainMatrix(
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

TMatrixD STTKFKalmanFilterManager::GetAMatrix(
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

double STTKFKalmanFilterManager::DeltaRadius(
    const STTKFStateVector& stateVector, double nextPhi, double dZ,
    double dE, double particle_mass) const
{
  // To Do: check all units
  double comp_delta_inv_radius_from_de = 0;

  double radius = stateVector.Radius();
  double tan    = stateVector.TanLambda();
  double mom    = SANDTrackerUtils::GetMomentumInMeVFromRadiusInMM(radius, tan) / 1000;
  // std::cout << mom << std::endl;
  
  double constant = SANDTrackerUtils::GetRadiusInMMToMomentumInGeVConstant() * SANDTrackerUtils::GetMagneticField() / sqrt(1 + pow(tan, 2));
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

STTKFStateVector STTKFKalmanFilterManager::PropagateState(
    const STTKFStateVector& stateVector, double dZ, double dE, double particle_mass)
{
  auto nextTanLambda = stateVector.TanLambda();

  auto cosNextPhi = cos(stateVector.Phi()) + dZ / stateVector.Radius();
  if (cosNextPhi > 1.)
    cosNextPhi = 1.;
  else if (cosNextPhi < -1.)
    cosNextPhi = -1.;

  auto nextPhi = acos(cosNextPhi);
  if (stateVector.Phi() < 0) nextPhi *= -1.;

  auto nextSignedInverseRadius =
      stateVector.SignedInverseRadius() +
      stateVector.Charge() * DeltaRadius(stateVector, nextPhi, dZ, dE, particle_mass);

  auto cosNextPhi_corr = cos(stateVector.Phi()) + dZ * stateVector.Charge() * nextSignedInverseRadius;
  if (cosNextPhi_corr > 1.)
    cosNextPhi_corr = 1.;
  else if (cosNextPhi_corr < -1.)
    cosNextPhi_corr = -1.;

  auto nextPhi_corr = acos(cosNextPhi_corr);
  if (stateVector.Phi() < 0) nextPhi_corr *= -1.;

  nextPhi = 0.5 * (nextPhi + nextPhi_corr);

  auto nextX = stateVector.X() + stateVector.ChargedRadius() *
                                     stateVector.TanLambda() *
                                     (nextPhi - stateVector.Phi());

  auto nextY = stateVector.Y() +
               stateVector.Radius() * (sin(nextPhi) - sin(stateVector.Phi()));
  return STTKFStateVector(nextX, nextY, nextSignedInverseRadius, nextTanLambda,
                          nextPhi);
}

STTKFStateCovarianceMatrix STTKFKalmanFilterManager::PropagateCovMatrix(
    const TMatrixD& covarianceMatrix, const TMatrixD& propagatorMatrix,
    const TMatrixD& processNoiseMatrix)
{
  TMatrixD propagatorMatrixTransported(TMatrixD::kTransposed, propagatorMatrix);
  
  return propagatorMatrix * covarianceMatrix * propagatorMatrixTransported +
         processNoiseMatrix;
}

STTKFMeasurement STTKFKalmanFilterManager::GetPrediction(
    Orientation orientation, const STTKFStateVector& stateVector)
{
  STTKFMeasurement projector(2, 1);
  if (orientation == Orientation::kHorizontal) {
    projector[0][0] = stateVector.Y();
    projector[1][0] =
        stateVector.Phi() + stateVector.Charge() * 0.5 * TMath::Pi();
  } else {
    projector[0][0] = stateVector.X();
    projector[1][0] = -stateVector.Charge() *
                      atan(stateVector.TanLambda() / sin(stateVector.Phi()));
  }
  return projector;
}

STTKFStateVector STTKFKalmanFilterManager::FilterState(
    const STTKFStateVector& stateVector, const TMatrixD& kalmanGainMatrix,
    const STTKFMeasurement& observed, const STTKFMeasurement& predicted)
{
  return STTKFStateVector(stateVector() +
                          kalmanGainMatrix * (observed - predicted));
}

STTKFStateCovarianceMatrix STTKFKalmanFilterManager::FilterCovMatrix(
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

STTKFStateVector STTKFKalmanFilterManager::smoothState(
    const STTKFStateVector& stateVectorFiltered,
    const STTKFStateVector& stateVectorPreviousSmoothed,
    const STTKFStateVector& stateVectorPreviousPredicted,
    const TMatrixD& theAMatrix)
{
  return STTKFStateVector(stateVectorFiltered() +
                          theAMatrix * (stateVectorPreviousSmoothed() -
                                        stateVectorPreviousPredicted()));
}

STTKFStateCovarianceMatrix STTKFKalmanFilterManager::smoothCovMatrix(
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

void STTKFKalmanFilterManager::Propagate(double& dE,
                                         double& dZ,
                                         double& beta)
{

  auto currentState = fThisTrack.GetStep(fCurrentStep);
  auto currentStage =
      currentState.GetStage(STTKFTrackStep::STTKFTrackStateStage::kFiltering);

  auto currentStateVector = currentStage.GetStateVector();
  auto predictedStateVector = PropagateState(currentStateVector, dZ, dE, particleInfo_.mass);
  auto nextPhi = predictedStateVector.Phi();
  
  auto processNoiseMatrix =
      GetProcessNoiseMatrix(currentStateVector, nextPhi, dZ, dE, fCurrentZ, particleInfo_.mass);
  auto propagatorMatrix =
      GetPropagatorMatrix(currentStateVector, nextPhi, dZ, dE, particleInfo_.mass);
  auto predictedCovMatrix = PropagateCovMatrix(
      currentStage.GetStateCovMatrix(), propagatorMatrix, processNoiseMatrix);

  STTKFTrackStep predictedTrackState;
  predictedTrackState.SetStage(
      STTKFTrackStep::STTKFTrackStateStage::kPrediction,
      STTKFState(predictedStateVector, predictedCovMatrix));
  predictedTrackState.SetPropagatorMatrix(propagatorMatrix);
  
  fThisTrack.AddStep(predictedTrackState);

  fCurrentStep++;
  fCurrentStage = STTKFTrackStep::STTKFTrackStateStage::kPrediction;
}

double STTKFKalmanFilterManager::EvalChi2(
    const STTKFMeasurement& observation, const STTKFMeasurement& prediction,
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

STTKFMeasurement STTKFKalmanFilterManager::GetMeasurementFromTracklet(const TVectorD& tracklet)
{
  STTKFMeasurement measurement(2, 1);
  // To Do: vertical and horizontal are outdated and confusing. Replace with something more meaningful.
  // Notice: vertical planes means horizontal measurements and the opposite
  if (fCurrentOrientation == Orientation::kVertical) {
    // To Do: Check units!
    measurement[0][0] = tracklet[0] / 1000.;
    measurement[1][0] = M_PI_2 - tracklet[2];
  } else {
    measurement[0][0] = tracklet[1] / 1000.;
    measurement[1][0] = tracklet[3];
  }

  return measurement;
}

int STTKFKalmanFilterManager::FindBestMatch(double& nextZ, const STTKFMeasurement& prediction,
    const TMatrixD& Sk)
{
  double best_chi = 1E9;
  auto& next_tracklets = z_to_tracklets_->at(nextZ);
  auto best_tracklet_index = -1;

  for (int i = 0; i < next_tracklets.size(); i++) {
    STTKFMeasurement measurement = GetMeasurementFromTracklet(next_tracklets[i]);

    auto chi2 = EvalChi2(measurement, prediction, Sk);
    if (chi2 < best_chi) {
      best_chi = chi2;
      best_tracklet_index = i;
    }
  }
  if (best_chi < 1.5) {
    return best_tracklet_index;
  } else {
    return -1;
  }
}

void STTKFKalmanFilterManager::SetNextOrientation()
{
  if (fCurrentOrientation == Orientation::kVertical) {
    fCurrentOrientation = Orientation::kHorizontal;
  } else {
    fCurrentOrientation = Orientation::kVertical;
  }
}

void STTKFKalmanFilterManager::Filter(const STTKFMeasurement& measurement,
                                      const STTKFMeasurement& prediction)
{

  auto currentState = fThisTrack.GetStep(fCurrentStep);
  auto predictedStage =
      currentState.GetStage(STTKFTrackStep::STTKFTrackStateStage::kPrediction);
  auto predictedStateVector = predictedStage.GetStateVector();
  auto predictedCovMatrix = predictedStage.GetStateCovMatrix();

  auto measurementNoiseMatrix = GetMeasurementNoiseMatrix();
  auto projectionMatrix = GetProjectionMatrix(fCurrentOrientation, predictedStateVector);
  auto kalmanGainMatrix = GetKalmanGainMatrix(
      predictedCovMatrix, projectionMatrix, measurementNoiseMatrix);
  auto filteredStateVector = FilterState(predictedStateVector, kalmanGainMatrix,
                                         measurement, prediction);
  auto filteredCovMatrix = FilterCovMatrix(predictedCovMatrix, projectionMatrix,
                                           measurementNoiseMatrix);

  fThisTrack.SetStage(fCurrentStep,
                      STTKFTrackStep::STTKFTrackStateStage::kFiltering,
                      STTKFState(filteredStateVector, filteredCovMatrix));
  fCurrentStage = STTKFTrackStep::STTKFTrackStateStage::kFiltering;

  SetNextOrientation();
}

void STTKFKalmanFilterManager::Smooth()
{

  auto currentState = fThisTrack.GetStep(fCurrentStep);
  auto filteredState =
      currentState.GetStage(STTKFTrackStep::STTKFTrackStateStage::kFiltering);
  auto filteredStateVector = filteredState.GetStateVector();
  auto filteredCovMatrix = filteredState.GetStateCovMatrix();

  if (fCurrentStep == int(fThisTrack.GetSteps().size()) - 1) {
    fThisTrack.SetStage(fCurrentStep,
                        STTKFTrackStep::STTKFTrackStateStage::kSmoothing,
                        STTKFState(filteredStateVector, filteredCovMatrix));
  } else {
    // previous state
    auto previousState = fThisTrack.GetStep(fCurrentStep + 1);

    // previous smoothed
    auto previousSmoothedState = previousState.GetStage(
        STTKFTrackStep::STTKFTrackStateStage::kSmoothing);
    auto previousSmoothedStateVector = previousSmoothedState.GetStateVector();
    auto previousSmoothedCovMatrix = previousSmoothedState.GetStateCovMatrix();
    
    // previous predicted
    auto previousPredictedState = previousState.GetStage(
        STTKFTrackStep::STTKFTrackStateStage::kPrediction);
    auto previousPredictedStateVector = previousPredictedState.GetStateVector();
    auto previousPredictedCovMatrix =
        previousPredictedState.GetStateCovMatrix();

    // current predicted
    auto predictedState = currentState.GetStage(
        STTKFTrackStep::STTKFTrackStateStage::kPrediction);
    auto predictedStateVector = predictedState.GetStateVector();
    auto predictedCovMatrix = predictedState.GetStateCovMatrix();

    auto nextPhi = previousPredictedStateVector.Phi();
    auto propagatorMatrix = currentState.GetPropagatorMatrix();

    auto theAMatrix = GetAMatrix(filteredCovMatrix, previousPredictedCovMatrix,
                                 propagatorMatrix);

    auto smoothedStateVector =
        smoothState(filteredStateVector, previousSmoothedStateVector,
                    previousPredictedStateVector, theAMatrix);
    auto smoothedCovMatrix =
        smoothCovMatrix(filteredCovMatrix, previousSmoothedCovMatrix,
                        previousPredictedCovMatrix, theAMatrix);

    fThisTrack.SetStage(fCurrentStep,
                        STTKFTrackStep::STTKFTrackStateStage::kSmoothing,
                        STTKFState(smoothedStateVector, smoothedCovMatrix));
    // currentState.SetStage(STTKFTrackStep::STTKFTrackStateStage::kSmoothing,
    // STTKFState(smoothedStateVector, smoothedCovMatrix));

  }
  fCurrentStage = STTKFTrackStep::STTKFTrackStateStage::kFiltering;
  fCurrentStep--;
}

void STTKFKalmanFilterManager::InitFromMC(TrackletMap* z_to_tracklets, const SParticleInfo& particleInfo)
{

  TMatrixD initial_cov_matrix(5, 5);
  initial_cov_matrix[0][0] = pow(200E-6, 2);
  initial_cov_matrix[1][1] = pow(200E-6, 2);
  initial_cov_matrix[2][2] = pow(0.1, 2);
  initial_cov_matrix[3][3] = pow(0.01, 2);
  initial_cov_matrix[4][4] = pow(0.01, 2);

  STTKFStateVector initial_state_vector = STTKFCheck::get_state_vector(particleInfo.mom * 1E-3,  // GeV
                                                                       particleInfo.pos * 1E-3,  // m
                                                                       particleInfo.charge);

  STTKFTrackStep trackStep;
  trackStep.SetStage(STTKFTrackStep::STTKFTrackStateStage::kPrediction,
                      STTKFState(initial_state_vector, initial_cov_matrix));
  trackStep.SetStage(STTKFTrackStep::STTKFTrackStateStage::kFiltering,
                      STTKFState(initial_state_vector, initial_cov_matrix));


  trackStep.SetPropagatorMatrix(initial_cov_matrix);
  
  fThisTrack.AddStep(trackStep);

  particleInfo_       = particleInfo;
  z_to_tracklets_     = z_to_tracklets;
  fCurrentStage       = STTKFTrackStep::STTKFTrackStateStage::kFiltering;
  fCurrentStep        = 0u;
  fCurrentZ           = particleInfo.pos.Z(); //Notice: UNITS!!  mm, why?
  fCurrentOrientation = Orientation::kVertical;
}

// To Do: implment a seeding algorithm
// void STTKFKalmanFilterManager::Init(const STTPlaneID& planeID, int clusterID)
  // {
  //   auto cluster = STTKFClusterManager::GetCluster(clusterID);
  //   auto trkParameter = cluster.GetRecoParameters().at(0).trk;

  //   double x, y, invR, tanL, phi;
  //   auto plane = STTStrawTubeTracker::GetPlane(planeID);
  //   auto defaultCharge = -1;
  //   auto planeOrientation = plane.GetOrientation();

  //   if (planeOrientation == STTPlane::EOrientation::kHorizontal) {
  //     x = SANDTrackerUtils::GetSANDInnerVolumeCenterPosition()[0];
  //     y = trkParameter.m * plane.GetZ() + trkParameter.q;
  //     invR = defaultCharge /
  //            SANDTrackerUtils::GetRadiusInMMFromPerpMomentumInGeV(1. /*GeV*/);
  //     tanL = 0.;
  //     phi = GetPhiFromTheta(atan(trkParameter.m), defaultCharge);
  //   } else {
  //     x = trkParameter.m * plane.GetZ() + trkParameter.q;
  //     ;
  //     y = SANDTrackerUtils::GetSANDInnerVolumeCenterPosition()[1];
  //     invR = defaultCharge /
  //            SANDTrackerUtils::GetRadiusInMMFromPerpMomentumInGeV(1. /*GeV*/);
  //     tanL = trkParameter.m;
  //     phi = 0.5 * TMath::Pi();
  //   }

  //   STTKFTrackStep trackStep;
  //   trackStep.SetPlaneID(planeID);
  //   trackStep.SetClusterIDForThisState(clusterID);
  //   STTKFStateVector stateVector(x, y, invR, tanL, phi);

  //   auto initialCovMatrix = GetInitialCovMatrix(stateVector, planeOrientation);

  //   trackStep.SetStage(STTKFTrackStep::STTKFTrackStateStage::kPrediction,
  //                      STTKFState(stateVector, initialCovMatrix));
  //   trackStep.SetStage(STTKFTrackStep::STTKFTrackStateStage::kFiltering,
  //                      STTKFState(stateVector, initialCovMatrix));
  //   fThisTrack.AddStep(trackStep);

  //   fCurrentStage = STTKFTrackStep::STTKFTrackStateStage::kFiltering;
  //   fCurrentStep = 0u;

  //   STTTRACKRECO_LOG(
  //       "INFO", TString::Format("State Vector     : %s",
  //                               SANDTrackerUtils::PrintStateVector(stateVector).Data())
  //                   .Data());
  //   STTTRACKRECO_LOG(
  //       "INFO", TString::Format("Covariance Matrix: %s",
  //                               SANDTrackerUtils::PrintMatrix(initialCovMatrix).Data())
  //                   .Data());
// }

void STTKFKalmanFilterManager::Run()
{
  // criterio per quando fermare la ricerca
  int stepLength = 1;

  // Notice: if currentZ is not in the map, the second condition is always true
  while (stepLength < 100 && std::distance(z_to_tracklets_->begin(), z_to_tracklets_->find(fCurrentZ)) >= stepLength) {
    // 1- propagate to [currentPlaneID - step]
    auto nextZ = std::prev(z_to_tracklets_->lower_bound(fCurrentZ), stepLength)->first;

    auto currentStep = fThisTrack.GetStep(fCurrentStep);
    auto filteredStateVector =
        currentStep.GetStage(STTKFTrackStep::STTKFTrackStateStage::kFiltering)
            .GetStateVector();

    auto dir = -1. * GetDirectiveCosinesFromStateVector(filteredStateVector);

    // To Do: check if this is still valid and add a real fix if needed
    if (dir.Z() > 0) {
      dir *= -1;
    }

    auto current_mom = SANDTrackerUtils::GetMomentumInMeVFromRadiusInMM(
                              filteredStateVector.Radius(),
                              filteredStateVector.TanLambda()) / 1000;
    double gamma = sqrt(current_mom * current_mom + particleInfo_.mass * particleInfo_.mass) /
                    particleInfo_.mass;
    double beta = sqrt(1 - pow(1 / gamma, 2));

    // To Do: check all units
    auto dE = STTKFGeoManager::GetDE(
                        nextZ, 
                        1000 * filteredStateVector.X(), 1000 * filteredStateVector.Y(), fCurrentZ, 
                        dir.X(), dir.Y(), dir.Z(),
                        beta, particleInfo_.mass, particleInfo_.charge) / 1000;

    double dZ = (nextZ - fCurrentZ) / 1000;

    Propagate(dE, dZ, beta);

    // 2- Search best match
    auto predictionStateVector = fThisTrack.GetStep(fCurrentStep)
            .GetStage(STTKFTrackStep::STTKFTrackStateStage::kPrediction).GetStateVector();
    auto predictionStateCovMatrix = fThisTrack.GetStep(fCurrentStep).GetStage(STTKFTrackStep::STTKFTrackStateStage::kPrediction)
            .GetStateCovMatrix();
    auto prediction = GetPrediction(fCurrentOrientation, predictionStateVector);

    auto measurementNoiseMatrix = GetMeasurementNoiseMatrix();
    auto projectionMatrix = GetProjectionMatrix(fCurrentOrientation, 
                                                predictionStateVector);
    TMatrixD projectionMatrixTransposed(TMatrixD::kTransposed,
                                          projectionMatrix);
    auto Sk = measurementNoiseMatrix + projectionMatrix *
                                              predictionStateCovMatrix *
                                              projectionMatrixTransposed;

    int tracklet_index = FindBestMatch(nextZ, prediction, Sk);

    // // 3- If it is found: step = 1
    // //    else step++
    if (tracklet_index != -1) {
      stepLength = 1;
      auto measurement = GetMeasurementFromTracklet(z_to_tracklets_->at(nextZ)[tracklet_index]);
      Filter(measurement, prediction);
      fCurrentZ = nextZ;
    } else {
      stepLength++;
      fThisTrack.RemoveLastStep();
      fCurrentStep--;
      fCurrentStage = STTKFTrackStep::STTKFTrackStateStage::kFiltering;
    }
  }

  while (fCurrentStep >= 0) Smooth();
}
