#ifndef STTKFTRACK_H
#define STTKFTRACK_H

#include "TMatrixD.h"

// #include "STTStrawTubeTracker.h"
#include "SANDGeoManager.h"

using STTKFStateCovarianceMatrix = TMatrixD;
using STTKFMeasurement = TMatrixD;

class STTKFStateVector {

 private:
  TMatrixD fVector;
  
 public:
  // default constructors ... probably to be delete in the future
  STTKFStateVector(): fVector(5,1) {};

  // constructors
  STTKFStateVector(double x, double y, double signed_inv_radius, double tanlambda, double phi): fVector(5,1) {
    fVector(0,0) = x;
    fVector(1,0) = y;
    fVector(2,0) = signed_inv_radius;
    fVector(3,0) = tanlambda;
    fVector(4,0) = phi;
  };

  // constructors
  STTKFStateVector(TMatrixD vector): fVector(vector) {};
  // STTKFStateVector(const TMatrixD& vec): fVector(vec) {};
  
  // copy constructor
  // STTKFStateVector(const STTKFStateVector& other): fVector(other.fVector) {};
  // STTKFStateVector(STTKFStateVector other): fVector(other.fVector) {};
  
  // destructor
  ~STTKFStateVector() {};

  // copy assignment
  STTKFStateVector& operator=(const TMatrixD& vec) {
    fVector = vec;
    return *this;
  };
  STTKFStateVector& operator=(TMatrixD vec) {
    std::swap(fVector, vec);
    return *this;
  };
  STTKFStateVector operator+ (STTKFStateVector p2) {
    this->fVector += p2();
    return *this;
  };

  STTKFStateVector operator- (STTKFStateVector p2) {
    this->fVector -= p2();
    return *this;
  };

  STTKFStateVector operator* (STTKFStateVector p2) {
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

class STTKFState {
  STTKFStateVector fVector;
  STTKFStateCovarianceMatrix fCovMatrix;

  public:
  STTKFState(): fVector(), fCovMatrix(5,5) {};
  STTKFState(STTKFStateVector vector, STTKFStateCovarianceMatrix matrix): fVector(vector), fCovMatrix(matrix) {};
  const STTKFStateVector& GetStateVector() const {return fVector; };
  const STTKFStateCovarianceMatrix& GetStateCovMatrix() const {return fCovMatrix; };
};

class STTKFTrackStep {

  public:
    enum class STTKFTrackStateStage {
      kPrediction,
      kFiltering,
      kSmoothing,
    };

  private:
    STTKFState fPrediction;
    STTKFState fFiltered;
    STTKFState fSmoothed;

    // the propagation that bring the vector in this state
    TMatrixD fPropagatorMatrix; 
    // TMatrixD fProjectionMatrix; 
    // TMatrixD fProcessNoiseMatrix; 
    // TMatrixD fMeasurementNoiseMatrix; 
    // TMatrixD fKalmanGainMatrix; 
    // TMatrixD fTheAMatrix; 

    // STTKFTrackStateStage fStage; // meglio  enumerato

    // ID piano di misura;
    SANDTrackerPlaneID fPlaneID;
    int fClusterID;

  public:
    STTKFTrackStep(): fPropagatorMatrix(5,5) {}; 
    //                   fProjectionMatrix(2,5),
    //                   fProcessNoiseMatrix(5,5),
    //                   fMeasurementNoiseMatrix(2,2),
    //                   fKalmanGainMatrix(5,2),
    //                   fTheAMatrix(5,5) {};
    void SetPlaneID(const SANDTrackerPlaneID& planeID) {fPlaneID = planeID; };
    const SANDTrackerPlaneID& GetPlaneID() const {return fPlaneID; };
    void SetClusterIDForThisState(int clusterID) { fClusterID = clusterID; };
    int GetClusterIDForThisState() const { return fClusterID; }
    void SetStage(STTKFTrackStateStage stage, STTKFState state);
    const STTKFState& GetStage(STTKFTrackStateStage stage) const;
    void SetPropagatorMatrix(TMatrixD pMatrix) { fPropagatorMatrix = pMatrix; };
    const TMatrixD GetPropagatorMatrix() { return fPropagatorMatrix; };
};

class STTKFTrack {
  private:
    std::vector<STTKFTrackStep> fSteps;
  public:
    const std::vector<STTKFTrackStep>& GetSteps() const {return fSteps; };
    const STTKFTrackStep& GetStep(int index) const {return fSteps.at(index); };
    void AddStep(STTKFTrackStep state) { fSteps.push_back(state); };
    void SetStage(int index, STTKFTrackStep::STTKFTrackStateStage stage, STTKFState state) { fSteps.at(index).SetStage(stage, state); };
    void SetClusterIDForState(int index, int clusterID) { fSteps.at(index).SetClusterIDForThisState(clusterID); };
    void RemoveLastStep() { fSteps.erase(fSteps.end()-1); };
};

#endif