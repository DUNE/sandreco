#ifndef SANDTrackerUTILS_H
#define SANDTrackerUTILS_H

#include <vector>
#include <math.h>

#include "TTreeReader.h"
#include "TMatrixD.h"
#include <TDecompChol.h>
#include <TRandom3.h>

#include "utils.h"

#include "SANDTrackerDigitCollection.h"
#include "SANDKFTrack.h"

using SANDTrackANDKFStateCovarianceMatrix = TMatrixD;
using SANDTrackANDKFMeasurement = TMatrixD;

#define IS_VERBOSE false

#if (IS_VERBOSE)
#define SANDTrackerTRACKRECO_LOG(level, message)                                                                                                       \
{                                                                                                                                              \
  std::cout << std::setw(40) << std::left << __FUNCTION__ << " in " <<                                     \
               std::setw(100) << std::left << TString::Format("%s:%d",__FILE__,__LINE__).Data() << " --> " \
                              << level << ": " << message << std::endl;                                    \
}
#else
#define SANDTrackerTRACKRECO_LOG(level, message)
#endif

class SANDTrackerUtils
{
 private:
  static TGeoManager* geo_;
  static const double kMagneticFieldInT_;
  static const double k_;
  static const double c_;

  static const double kEdepSimDensityToGCM3_;
    // https://github.com/ClarkMcGrew/edep-sim/blob/master/README.md#reading-the-output
    // Be aware that in the saved TGeoManager object, the masses and densities are 
    // also in CLHEP units, so that 1 kilogram equals 6.24x10^24^ MeV ns^2^ mm^-2^, 
    // and densities are in units of 6.24x^24^ MeV ns^2^ mm^-5^.

  static double getDensityInGCM3() {return geo_->GetCurrentNode()->GetVolume()->GetMaterial()->GetDensity()/kEdepSimDensityToGCM3_; };
  static double getPathLengthInCM() {return geo_->GetStep() * 0.1; };

 public:
  SANDTrackerUtils(){};
  ~SANDTrackerUtils(){};
  static void clear();
  static bool areAdjacent(const sand_geometry::tracker::CellID &tub1, const sand_geometry::tracker::CellID &tub2);
  static inline double getTubeRadius() { return 2.5; };
  static void init(TGeoManager *geo) {geo_ = geo;};
  static TGeoManager* getGeoManager() {return geo_; };
  

  static inline double getX0(int Z, int A) {
    //https://cds.cern.ch/record/1279627/files/PH-EP-Tech-Note-2010-013.pdf
    //https://halldweb.jlab.org/DocDB/0004/000439/001/radlen.pdf
    //  The equation is an approximation while
    //  the second a result of a compact fit 
    //  to the data with an accuracy of better
    //  than 2.5% for all elements except He
    return 716.408 /*g/cm2*/ * A / (Z * (Z+1) * log(287/sqrt(Z)));
  }

  static inline double getDEInMeV(double crossedMaterialInGCM2) {
    return crossedMaterialInGCM2 * 2. /*MeV/(g/cm2)*/;
  }

  static inline double getDEInGeV(double crossedMaterialInGCM2) {
    return getDEInMeV(crossedMaterialInGCM2) * 1E-3;
  }

  static inline double getMCSSigmaAngleFromMomentumInMeV(double p, 
                                       double beta, 
                                       double pathLengthInX0) {
    // The charge of the particle is assumed
    // to be one unit of electron charge
    return 13.6 /*MeV*/ / (p*beta) * sqrt(pathLengthInX0) * (1 + 0.038 * log(pathLengthInX0/(beta*beta)));
  }

  // To Do: check all units
  static inline double getRadiusInMMToMomentumInGeVConstant() {return 0.299792458; /* GeV/(m*T) */ };
  // static inline double getPerpMomentumInGeVFromRadiusInMM(double radius) {return getRadiusInMMToMomentumInGeVConstant() * radius * SANDTrackANDKFGeoManager::getMagneticField(); };
  static double getPerpMomentumInGeVFromRadiusInMM(double radius);
  static double getRadiusInMMFromPerpMomentumInGeV(double perpMom);
  static inline double getMomentumInGeVFromRadiusInMM(double radius, double tanl) {return getPerpMomentumInGeVFromRadiusInMM(radius) * sqrt(1 + tanl*tanl); };
  static inline double getMomentumInMeVFromRadiusInMM(double radius, double tanl) {return 1.E3 * getMomentumInGeVFromRadiusInMM(radius, tanl); };
  static inline double getSigmaPositionMeasurement() {return 200E-6 /*m*/; };
  static inline double getSigmaAngleMeasurement() {return 0.2 /*rad*/; };
  static double getMagneticField() { return kMagneticFieldInT_; };
  static double getk() { return k_; };
  static double getc() { return c_; };

  static TString printMatrix(const TMatrixD& m);

  static TVector3 getCartesianCoordinateFromCylindrical(double radius, double angle, double x);


  static double getCrossedMaterialInGCM2(double z, 
                              double px, double py, double pz,
                              double sx, double sy, double sz);
    static double getPathLengthInX0(double z, 
                            double px, double py, double pz,
                            double sx, double sy, double sz);

    static double getPathLengthInCM(double z, 
                            double px, double py, double pz,
                            double sx, double sy, double sz);

    static double getDE(double z, 
                                double px, double py, double pz,
                                double sx, double sy, double sz,
                                double beta, double mass, int charge);

  friend class SANDTrackerStrawTubeTracker;
};






namespace sand_reco
{

namespace kf
{

namespace utils 
{

TVector2 getBFieldPerp(const TVector3& v);

double getBFieldLong(const TVector3& v);

double getTanOfDipAngle(const TVector3& mom);

int getRotationVersus(int charge);

int getCharge(int versus);

double getDirectionAngle(const TVector2& dir);

double getRotationAngle(const TVector2& dir, int charge);

double getRadius(double perp_mom);

TVector2 getCircleCenter(const TVector2& momentum, const TVector2& position, int charge);

std::pair<double, double> getCircleYs(double z, double radius, const TVector2& center);

double getRotationAngle(double z, double y, const TVector2& center);

double getDeltaPhi(double phi, double previous_phi, int versus);

double getX(double radius, double x_0, double delta_phi, double tan_lambda, int charge);

double getY(const TVector2& center, double radius, double phi);

double getZ(const TVector2& center, double radius, double phi);

TVector3 getVectorMomentum(double radius, double phi, double tan_lambda, int versus);

sand_reco::kf::StateVector getStateVector(TVector3 mom, TVector3 pos, int charge);

sand_reco::kf::State Seed3Points(std::array<double,3> xyz0, std::array<double,3> xyz1, std::array<double,3> xyz2, double sy, double sx);

// Notice: this was used to validate the KF. It is used to generate
// ideal trajectories and to store their information
class TrajectoryParameters;

class ParticleState {
    private:
        TVector3 position_;
        TVector3 momentum_;
    public:
        ParticleState(): position_(), momentum_() {};
        ParticleState(const TVector3& p, const TVector3& m): position_(p), momentum_(m) {};
        ParticleState(const ParticleState& particle_state): position_(particle_state.position_), momentum_(particle_state.momentum_) {}
        ParticleState(const sand_reco::kf::StateVector& vector, double z);

        TrajectoryParameters getTrajectoryParameter(int charge) const;
        sand_reco::kf::StateVector getStateVector(int charge) const;
        const TVector3& getPosition() const { return position_; };
        const TVector3& getMomentum() const { return momentum_; };
        TVector3& getPosition() { return position_; };
        TVector3& getMomentum() { return momentum_; };
        void setPosition(const TVector3 v) { position_ = v; };
        void setMomentum(const TVector3 v) { momentum_ = v; };

        ParticleState operator =(const ParticleState& p) {
            this->setPosition(p.getPosition());
            this->setMomentum(p.getMomentum());
            return *this;
        }
};





class TrajectoryParameters {
    public:
        double radius_;
        double versus_of_rot_; // right-hand rule in right-handed coordinate system, z == beam, y == vertical
        double tan_lambda_;
        TVector2 center_of_rot_;
        double phi_0_;
        double x_0_;
        std::pair<double, double> getPhiPair(double z) const;
        double getSmallestDeltaPhi(double z, double last_phi) const;
        std::vector<double> getDeltaPhis(std::vector<double> zs) const;
        ParticleState getParticleState(double delta_phi) const;
        std::vector<ParticleState> getParticleStatesFromDeltaPhi(std::vector<double> delta_phi) const;
    public:
        TrajectoryParameters() {};
        TrajectoryParameters(double r, double v, double t, const TVector2& c_, double p, double x): 
            radius_(r), versus_of_rot_(v), tan_lambda_(t), center_of_rot_(c_), phi_0_(p), x_0_(x) {};
        std::vector<ParticleState> getParticleStatesFromZ(std::vector<double> zs) const;
};

class CovMatrixPropCheckOutput {
    public:
        double dz;
        TMatrixD initial_state_propagated_;
        TMatrixD initial_covariance_propagated_;
        std::vector<sand_reco::kf::StateVector> propagated_states_;
        TMatrixD mean_of_propagated_states_;
        TMatrixD covariance_of_propagated_states_;
        TMatrixD variance_of_propagated_states_;
        double distance_;
        CovMatrixPropCheckOutput(): 
            initial_state_propagated_(1,5),
            initial_covariance_propagated_(5,5),
            mean_of_propagated_states_(1,5),
            covariance_of_propagated_states_(5,5),
            variance_of_propagated_states_(5,5) {};
};

std::vector<sand_reco::kf::StateVector> GenerateStateVectors(const sand_reco::kf::StateVector& state, const TMatrixD& cov, int n);

TMatrixD getMean(const std::vector<sand_reco::kf::StateVector>& states);

TMatrixD getCov(const std::vector<sand_reco::kf::StateVector>& states, const TMatrixD& mean);

TMatrixD getVar(const std::vector<sand_reco::kf::StateVector>& states, const TMatrixD& cov);

void getMeanAndCov(const std::vector<sand_reco::kf::StateVector>& states, TMatrixD& mean, TMatrixD& cov);

using propagation = std::vector<ParticleState>;
} // namespace utils
} // namespace kf
} // namespace sand_reco

#endif