#ifndef SANDKFTRACK_H
#define SANDKFTRACK_H

#include "TMatrixD.h"

// #include "SANDStrawTubeTracker.h"
#include "SANDGeoManager.h"

namespace sand_reco
{

namespace kf
{

using StateCovarianceMatrix = TMatrixD;
using Measurement = TMatrixD;

class StateVector {

 private:
  TMatrixD vector_;
  
 public:
  // default constructors ... probably to be delete in the future
  StateVector(): vector_(5,1) {};

  // constructors
  StateVector(double arg_x, double arg_y, double arg_signed_inv_radius, double arg_tan_lambda, double arg_phi): vector_(5,1) {
    vector_(0,0) = arg_x;
    vector_(1,0) = arg_y;
    vector_(2,0) = arg_signed_inv_radius;
    vector_(3,0) = arg_tan_lambda;
    vector_(4,0) = arg_phi;
  };

  // constructors
  StateVector(TMatrixD vector): vector_(vector) {};
  // StateVector(const TMatrixD& vec): vector_(vec) {};
  
  // copy constructor
  // StateVector(const StateVector& other): vector_(other.vector_) {};
  // StateVector(StateVector other): vector_(other.vector_) {};
  
  // destructor
  ~StateVector() {};

  // copy assignment
  StateVector& operator=(const TMatrixD& vec) {
    vector_ = vec;
    return *this;
  };
  StateVector& operator=(TMatrixD vec) {
    std::swap(vector_, vec);
    return *this;
  };
  StateVector operator+ (StateVector p2) {
    this->vector_ += p2();
    return *this;
  };

  StateVector operator- (StateVector p2) {
    this->vector_ -= p2();
    return *this;
  };

  StateVector operator* (StateVector p2) {
    this->vector_ = ElementMult(this->vector_, p2());
    return *this;
  };

  // operator()
  const TMatrixD& operator()() const { return vector_; };

  // getters
  inline double x() const { return vector_(0, 0); };
  inline double y() const { return vector_(1, 0); };
  inline double signedInverseRadius() const { return vector_(2, 0); };
  inline double tanLambda() const { return vector_(3, 0); };
  inline double phi() const { return vector_(4, 0); };

  // usefull function
  inline int charge() const { return std::signbit(signedInverseRadius()) == false ? +1 : -1; };
  inline double chargedRadius() const { return 1./signedInverseRadius(); };
  inline double radius() const { return chargedRadius() * charge(); };

  // Setters... to be removed
  // void X(double val) { vector_(0, 0) = val; };
  // void Y(double val) { vector_(1, 0) = val; };
  // void signedInverseRadius(double val) { vector_(2, 0) = val; };
  // void TanLambda(double val) { vector_(3, 0) = val; };
  // void Phi(double val) { vector_(4, 0) = val; };
};

class State {
  StateVector vector_;
  StateCovarianceMatrix cov_matrix_;

  public:
  State(): vector_(), cov_matrix_(5,5) {};
  State(StateVector vector, StateCovarianceMatrix matrix): vector_(vector), cov_matrix_(matrix) {};
  const StateVector& getStateVector() const {return vector_; };
  const StateCovarianceMatrix& getStateCovMatrix() const {return cov_matrix_; };
};

class TrackStep {

  public:
    enum class TrackStateStage {
      kPrediction,
      kFiltering,
      kSmoothing,
    };

  private:
    State prediction_;
    State filtered_;
    State smoothed_;
    std::vector<double> innovation_;
    std::vector<dg_wire> digits_;
    double z_;
    double x_;
    double y_;
    Measurement measurement_;
    double chi2_;

    // the propagation that bring the vector in this state
    TMatrixD propagator_matrix_; 
    // TMatrixD fProjectionMatrix; 
    // TMatrixD fProcessNoiseMatrix; 
    // TMatrixD fMeasurementNoiseMatrix; 
    // TMatrixD fKalmanGainMatrix; 
    // TMatrixD fTheAMatrix; 

    // TrackStateStage fStage; // meglio  enumerato

    // ID piano di misura;
    sand_geometry::tracker::PlaneID plane_id_;
    int clusterid_;

  public:
    TrackStep(): propagator_matrix_(5,5) {}; 
    //                   fProjectionMatrix(2,5),
    //                   fProcessNoiseMatrix(5,5),
    //                   fMeasurementNoiseMatrix(2,2),
    //                   fKalmanGainMatrix(5,2),
    //                   fTheAMatrix(5,5) {};
    void setPlaneID(const sand_geometry::tracker::PlaneID& plane_id) {plane_id_ = plane_id; };
    const sand_geometry::tracker::PlaneID& getPlaneID() const {return plane_id_; };
    void setClusterIDForThisState(int cluster_id) { clusterid_ = cluster_id; };
    int getClusterIDForThisState() const { return clusterid_; }
    void setStage(TrackStateStage stage, State state);
    const State& getStage(TrackStateStage stage) const;
    void setPropagatorMatrix(TMatrixD propagator_matrix) { propagator_matrix_ = propagator_matrix; };
    const TMatrixD getPropagatorMatrix() { return propagator_matrix_; };
    void setInnovation(std::vector<double> innovation) { innovation_ = innovation; };
    const std::vector<double>& getInnovation() const { return innovation_ ;};
    void setZ(double z){z_ = z;};
    double getZ() const {return z_;};
    void setX(double x){x_ = x;};
    double getX() const {return x_;};
    void setY(double y){y_ = y;};
    double getY() const {return y_;};
    void setMeasurement(Measurement measurement) {measurement_ = measurement;};
    const Measurement& getMeasurement() const {return measurement_;};
    void setChi2(double chi2) {chi2_ = chi2;};
    double getChi2() {return chi2_;};
    void addDigits(std::vector<dg_wire> digits ){digits_ = digits;};
    std::vector<dg_wire> getDigits() const {return digits_;};
    

};

class Track {
  private:
    std::vector<TrackStep> steps_;
  public:
    const std::vector<TrackStep>& getSteps() const {return steps_; };
    const TrackStep& getStep(int index) const {return steps_.at(index); };
    void addStep(TrackStep state) { steps_.push_back(state); };
    void setStage(int index, TrackStep::TrackStateStage stage, State state) { steps_.at(index).setStage(stage, state); };
    void setInnovation(int index, std::vector<double> innovation) { steps_.at(index).setInnovation(innovation); };
    void setZ(int index, double z){steps_.at(index).setZ(z); };
    void setX(int index, double x){steps_.at(index).setX(x); };
    void setY(int index, double y){steps_.at(index).setY(y); };
    void setMeasurement(int index, Measurement measurement){steps_.at(index).setMeasurement(measurement); };
    void setChi2(int index, double chi2){steps_.at(index).setChi2(chi2); };
    void addDigits(int index, std::vector<dg_wire> digits) {steps_.at(index).addDigits(digits);};
    void setClusterIDForState(int index, int cluster_id) { steps_.at(index).setClusterIDForThisState(cluster_id); };
    void removeLastStep() { steps_.erase(steps_.end()-1); };
    void Clear() {steps_.clear();}
};
} // namespace kf
} // namespace sand_reco
#endif