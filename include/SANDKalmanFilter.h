#ifndef SANDKALMANFILTER_H
#define SANDKALMANFILTER_H

#include "TGeoManager.h"
#include "TMatrixD.h"

#include "SANDTrackerClusterCollection.h"
#include "SANDKFTrack.h"
#include "SANDTrackerUtils.h"

struct SParticleInfo {
  int charge;
  double mass;
  int pdg_code;
  int id;

  TVector3 pos;
  TVector3 mom;
};

class SANDKFUtils;

using SANDKFStateCovarianceMatrix = TMatrixD;
using SANDKFMeasurement = TMatrixD;
using TrackletMap = std::map<double, std::vector<TVectorD>>;

class SANDKalmanFilterManager {

  private:
    SANDKFTrack fThisTrack;

    SANDKFTrackStep::SANDKFTrackStateStage fCurrentStage; // forward or backward
    int fCurrentStep; // index of the SANDKFTrackStep in SANDKFTrack
    double fCurrentZ; 
    TrackletMap* z_to_tracklets_;
    SParticleInfo particleInfo_;


  public:
    enum class Orientation {
      kVertical,
      kHorizontal
    };
    Orientation GetOrientation() {return fCurrentOrientation;};
    TVector3 GetDirectiveCosinesFromStateVector(const SANDKFStateVector& stateVector);
    double GetPhiFromTheta(double theta, int charge) { return theta - charge * 0.5*TMath::Pi(); };
    double GetThetaFromPhi(double phi, int charge) { return phi + charge * 0.5*TMath::Pi(); };
    double GetThetaFromPhi(const SANDKFStateVector& stateVector) { return GetThetaFromPhi(stateVector.Phi(), stateVector.Charge()); };
    SANDKFMeasurement GetMeasurementFromCluster(int clusterID);

    // TMatrixD GetInitialCovMatrix(const SANDKFStateVector& stateVector, const Orientation& orientation);
    TMatrixD GetPropagatorMatrix(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass);  // Propagate and Smooth
    TMatrixD GetProcessNoiseMatrix(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double z, double particle_mass);          // Propagate
    TMatrixD GetMeasurementNoiseMatrix();               // Filter
    TMatrixD GetProjectionMatrix(Orientation orientation, const SANDKFStateVector& stateVector);                     // Filter
    TMatrixD GetKalmanGainMatrix(const TMatrixD& covarianceMatrix,
                                 const TMatrixD& projectionMatrix,
                                 const TMatrixD& measurementNoiseMatrix);                     // Filter
    TMatrixD GetAMatrix(const TMatrixD& covarianceMatrixFiltered,
                        const TMatrixD& covarianceMatrixNextPredicted,
                        const TMatrixD& propagatorMatrix);                     // Smooth
  
  // private:
  public:
    Orientation fCurrentOrientation = Orientation::kVertical;

    SANDKFMeasurement GetMeasurementFromTracklet(const TVectorD& tracklet);
    double DeltaRadius(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const;
    inline double DEDTanl(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { auto tan = stateVector.TanLambda(); return dE * tan / (1 + tan*tan); };
    inline double DEDPhi(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { 
      if (dE == 0) {
        return 0.;
      } else {
        return -dE / tan(stateVector.Phi());
      }
    };
    inline double DDeltaInvRDx(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { return 0.; };
    inline double DDeltaInvRDy(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { return 0.; };
    inline double DDeltaInvRDInvCR(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { 
      double tan = stateVector.TanLambda();
      double constant_tan = 0.299792458 * 0.6 / sqrt(1 + pow(tan, 2));
      double x = 1./stateVector.Radius();
      double squared_root = sqrt(pow(constant_tan, 2) / pow(x, 2) + pow(particle_mass, 2));
      double new_derivatives_r = stateVector.Charge() * (-dE * (2 * pow(constant_tan, 2) + 3 * pow(particle_mass, 2) * pow(x, 2)) / pow(constant_tan, 2) / squared_root);
      return new_derivatives_r;
    };

    inline double DDeltaInvRDTanl(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { 
      auto tan = stateVector.TanLambda(); 
      if (dE != 0) {
        double old_derivatives = DeltaRadius(stateVector, nextPhi, dZ, dE, particle_mass) * (1./dE * DEDTanl(stateVector, nextPhi, dZ, dE, particle_mass) + tan / (1 + tan*tan));
        double c = -(pow(1./stateVector.Radius(), 3)) / pow(0.299792458 * 0.6, 2);
        double f = tan*tan + 1;
        double g = sqrt(pow(0.299792458 * 0.6, 2) /
                         pow(1./stateVector.Radius(), 2) / (1 + tan*tan) + pow(particle_mass, 2));
        double h = dE;
        double constant_r = 0.299792458 * 0.6 / (1. / stateVector.Radius());
        double f_derivatives = 2 * tan;
        double g_derivatives = -pow(constant_r, 2) * tan / pow(tan*tan + 1, 2) / sqrt(pow(constant_r, 2) / (tan*tan + 1) + pow(particle_mass, 2));
        double h_derivatives = dE * tan / (1 + tan*tan);
        double new_derivatives_tan = c * (f_derivatives*g*h + f*g_derivatives*h + f*g*h_derivatives);
        return new_derivatives_tan;
      } else {
        return DeltaRadius(stateVector, nextPhi, dZ, dE, particle_mass) * ( tan / (1 + tan*tan)); 
      }
    };
    
    inline double DDeltaInvRDPhi(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { 
      if (dE == 0) {
        return 0.;
      } else {
        return DeltaRadius(stateVector, nextPhi, dZ, dE, particle_mass) / dE * DEDPhi(stateVector, nextPhi, dZ, dE, particle_mass); 
      }
    };

    inline double DPhiDCosPhi(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { 
      return stateVector.Charge()/sqrt(1 - pow(cos(stateVector.Phi()) + dZ/stateVector.Radius(),2));
    };

    inline double DxDx(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { return 1; }
    inline double DyDx(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { return 0.; }
    inline double DInvCRDx(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { return stateVector.Charge() * DDeltaInvRDx(stateVector, nextPhi, dZ, dE, particle_mass); }
    inline double DTanlDx(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { return 0.; }
    inline double DPhiDx(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { return 0.; }

    inline double DxDy(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { return 0.; }
    inline double DyDy(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { return 1.; }
    inline double DInvCRDy(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { return stateVector.Charge() * DDeltaInvRDy(stateVector, nextPhi, dZ, dE, particle_mass); }
    inline double DTanlDy(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { return 0.; }
    inline double DPhiDy(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { return 0.; }

    inline double DPhiDInvCR(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const 
    {
      return stateVector.Charge() * dZ * DPhiDCosPhi(stateVector, nextPhi, dZ, dE, particle_mass);
    }
    inline double DxDInvCR(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const
    {
      return -pow(stateVector.ChargedRadius(),2) * stateVector.TanLambda() * (nextPhi - stateVector.Phi()) + stateVector.ChargedRadius() * stateVector.TanLambda() * DPhiDInvCR(stateVector, nextPhi, dZ, dE, particle_mass);
    }
    inline double DyDInvCR(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const
    {
      return -pow(stateVector.Radius(), 2) * stateVector.Charge() * (sin(nextPhi) - sin(stateVector.Phi())) 
             + stateVector.Radius() * cos(nextPhi) * DPhiDInvCR(stateVector, nextPhi, dZ, dE, particle_mass);
    }
    inline double DInvCRDInvCR(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const
    {
      return 1. + stateVector.Charge() * DDeltaInvRDInvCR(stateVector, nextPhi, dZ, dE, particle_mass);
    }
    inline double DTanlDInvCR(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { return 0.; }

    inline double DxDTanl(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const
    {
      return stateVector.ChargedRadius() * (nextPhi - stateVector.Phi());
    }
    inline double DyDTanl(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { return 0.; }
    inline double DInvCRDTanl(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const
    {
      return stateVector.Charge() * DDeltaInvRDTanl(stateVector, nextPhi, dZ, dE, particle_mass);
    }
    inline double DTanlDTanl(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { return 1.; }
    inline double DPhiDTanl(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { return 0.; }

     inline double DPhiDPhi(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const
    {
      int DPhiDPhi_sign = ((nextPhi - stateVector.Phi()) * stateVector.Charge()) > 0 ? 1 : -1;
      return -sin(stateVector.Phi()) * DPhiDCosPhi(stateVector, nextPhi, dZ, dE, particle_mass);
    }
     inline double DxDPhi(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const
    {
      return stateVector.ChargedRadius() * stateVector.TanLambda() * (DPhiDPhi(stateVector, nextPhi, dZ, dE, particle_mass) - 1.);
    }
     inline double DyDPhi(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const
    {
      return stateVector.Radius() * (cos(nextPhi) * DPhiDPhi(stateVector, nextPhi, dZ, dE, particle_mass) - cos(stateVector.Phi()));
    }
    inline double DInvCRDPhi(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { return stateVector.Charge() * DDeltaInvRDPhi(stateVector, nextPhi, dZ, dE, particle_mass); }
    inline double DTanlDPhi(const SANDKFStateVector& stateVector, double nextPhi, double dZ, double dE, double particle_mass) const { return 0.; }
  
  public:
    SANDKFStateVector PropagateState(const SANDKFStateVector& stateVector, double dZ, double dE, double particle_mass);
    SANDKFStateCovarianceMatrix PropagateCovMatrix(const TMatrixD& covarianceMatrix,
                                                  const TMatrixD& propagatorMatrix,
                                                  const TMatrixD& processNoiseMatrix);

    SANDKFStateVector FilterState(const SANDKFStateVector& stateVector, 
                                 const TMatrixD& kalmanGainMatrix,
                                 const SANDKFMeasurement& observed, 
                                 const SANDKFMeasurement& predicted);

    SANDKFStateCovarianceMatrix FilterCovMatrix(const TMatrixD& covarianceMatrix,
                                               const TMatrixD& projectionMatrix,
                                               const TMatrixD& measurementNoiseMatrix);

    SANDKFStateVector smoothState(const SANDKFStateVector& stateVectorFiltered, 
                                 const SANDKFStateVector& stateVectorPreviousSmoothed,
                                 const SANDKFStateVector& stateVectorPreviousPredicted, 
                                 const TMatrixD& theAMatrix);

    SANDKFStateCovarianceMatrix smoothCovMatrix(const TMatrixD& covarianceMatrixFiltered,
                                               const TMatrixD& covarianceMatrixPreviousSmoothed,
                                               const TMatrixD& covarianceMatrixPreviousPredicted,
                                               const TMatrixD& theAMatrix);
    
    SANDKFMeasurement GetPrediction(Orientation orientation, const SANDKFStateVector& stateVector);

    void Propagate(double& dE, double& dZ, double& beta);
    double EvalChi2(const SANDKFMeasurement& measurement, const SANDKFMeasurement& prediction, const TMatrixD& measurementNoiseMatrix);
    int FindBestMatch(double& nextZ, const SANDKFMeasurement& prediction, const TMatrixD& measurementNoiseMatrix);
    void SetNextOrientation();
    void Filter(const SANDKFMeasurement& measurement, const SANDKFMeasurement& prediction);
    void Smooth();
    void InitFromMC(TrackletMap* z_to_tracklets, const SParticleInfo& particloInfo);
    void Run();
    const SANDKFTrack& GetTrack() {return fThisTrack; };
  
  friend class SANDKFUtils;
};

#endif