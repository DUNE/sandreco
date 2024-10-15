#ifndef SANDKFTRACK_H
#define SANDKFTRACK_H

#include "TMatrixD.h"

// #include "SANDStrawTubeTracker.h"
#include "SANDGeoManager.h"

using SANDKFStateCovarianceMatrix = TMatrixD;
using SANDKFMeasurement = TMatrixD;

class SANDKFStateVector {

 private:
  TMatrixD fVector;
  
 public:
  // default constructors ... probably to be delete in the future
  SANDKFStateVector(): fVector(5,1) {};

  // constructors
  SANDKFStateVector(double x, double y, double signed_inv_radius, double tanlambda, double phi): fVector(5,1) {
    fVector(0,0) = x;
    fVector(1,0) = y;
    fVector(2,0) = signed_inv_radius;
    fVector(3,0) = tanlambda;
    fVector(4,0) = phi;
  };

  // constructors
  SANDKFStateVector(TMatrixD vector): fVector(vector) {};
  // SANDKFStateVector(const TMatrixD& vec): fVector(vec) {};
  
  // copy constructor
  // SANDKFStateVector(const SANDKFStateVector& other): fVector(other.fVector) {};
  // SANDKFStateVector(SANDKFStateVector other): fVector(other.fVector) {};
  
  // destructor
  ~SANDKFStateVector() {};

  // copy assignment
  SANDKFStateVector& operator=(const TMatrixD& vec) {
    fVector = vec;
    return *this;
  };
  SANDKFStateVector& operator=(TMatrixD vec) {
    std::swap(fVector, vec);
    return *this;
  };
  SANDKFStateVector operator+ (SANDKFStateVector p2) {
    this->fVector += p2();
    return *this;
  };

  SANDKFStateVector operator- (SANDKFStateVector p2) {
    this->fVector -= p2();
    return *this;
  };

  SANDKFStateVector operator* (SANDKFStateVector p2) {
    this->fVector = ElementMult(this->fVector, p2());
    return *this;
  };

  // operator()
  const TMatrixD& operator()() const { return fVector; };

  // Getters
  inline double X() const { return fVector(0, 0); };
  inline double Y() const { return fVector(1, 0); };
  inline double SignedInverseRadius() const { return fVector(2, 0); };
  inline double TanLambda() const { return fVector(3, 0); };
  inline double Phi() const { return fVector(4, 0); };

  // usefull function
  inline int Charge() const { return std::signbit(SignedInverseRadius()) == false ? +1 : -1; };
  inline double ChargedRadius() const { return 1./SignedInverseRadius(); };
  inline double Radius() const { return ChargedRadius() * Charge(); };

  // Setters... to be removed
  // void X(double val) { fVector(0, 0) = val; };
  // void Y(double val) { fVector(1, 0) = val; };
  // void SignedInverseRadius(double val) { fVector(2, 0) = val; };
  // void TanLambda(double val) { fVector(3, 0) = val; };
  // void Phi(double val) { fVector(4, 0) = val; };
};

class SANDKFState {
  SANDKFStateVector fVector;
  SANDKFStateCovarianceMatrix fCovMatrix;

  public:
  SANDKFState(): fVector(), fCovMatrix(5,5) {};
  SANDKFState(SANDKFStateVector vector, SANDKFStateCovarianceMatrix matrix): fVector(vector), fCovMatrix(matrix) {};
  const SANDKFStateVector& GetStateVector() const {return fVector; };
  const SANDKFStateCovarianceMatrix& GetStateCovMatrix() const {return fCovMatrix; };
};

class SANDKFTrackStep {

  public:
    enum class SANDKFTrackStateStage {
      kPrediction,
      kFiltering,
      kSmoothing,
    };

  private:
    SANDKFState fPrediction;
    SANDKFState fFiltered;
    SANDKFState fSmoothed;

    // the propagation that bring the vector in this state
    TMatrixD fPropagatorMatrix; 
    // TMatrixD fProjectionMatrix; 
    // TMatrixD fProcessNoiseMatrix; 
    // TMatrixD fMeasurementNoiseMatrix; 
    // TMatrixD fKalmanGainMatrix; 
    // TMatrixD fTheAMatrix; 

    // SANDKFTrackStateStage fStage; // meglio  enumerato

    // ID piano di misura;
    SANDTrackerPlaneID fPlaneID;
    int fClusterID;

  public:
    SANDKFTrackStep(): fPropagatorMatrix(5,5) {}; 
    //                   fProjectionMatrix(2,5),
    //                   fProcessNoiseMatrix(5,5),
    //                   fMeasurementNoiseMatrix(2,2),
    //                   fKalmanGainMatrix(5,2),
    //                   fTheAMatrix(5,5) {};
    void SetPlaneID(const SANDTrackerPlaneID& planeID) {fPlaneID = planeID; };
    const SANDTrackerPlaneID& GetPlaneID() const {return fPlaneID; };
    void SetClusterIDForThisState(int clusterID) { fClusterID = clusterID; };
    int GetClusterIDForThisState() const { return fClusterID; }
    void SetStage(SANDKFTrackStateStage stage, SANDKFState state);
    const SANDKFState& GetStage(SANDKFTrackStateStage stage) const;
    void SetPropagatorMatrix(TMatrixD pMatrix) { fPropagatorMatrix = pMatrix; };
    const TMatrixD GetPropagatorMatrix() { return fPropagatorMatrix; };
};

class SANDKFTrack {
  private:
    std::vector<SANDKFTrackStep> fSteps;
  public:
    const std::vector<SANDKFTrackStep>& GetSteps() const {return fSteps; };
    const SANDKFTrackStep& GetStep(int index) const {return fSteps.at(index); };
    void AddStep(SANDKFTrackStep state) { fSteps.push_back(state); };
    void SetStage(int index, SANDKFTrackStep::SANDKFTrackStateStage stage, SANDKFState state) { fSteps.at(index).SetStage(stage, state); };
    void SetClusterIDForState(int index, int clusterID) { fSteps.at(index).SetClusterIDForThisState(clusterID); };
    void RemoveLastStep() { fSteps.erase(fSteps.end()-1); };
};

#endif