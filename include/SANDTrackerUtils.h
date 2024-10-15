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
  static TGeoManager* fGeo;
  static const double kMagneticFieldInT;
  static const double k;
  static const double c;

  static const double kEdepSimDensityToGCM3;
    // https://github.com/ClarkMcGrew/edep-sim/blob/master/README.md#reading-the-output
    // Be aware that in the saved TGeoManager object, the masses and densities are 
    // also in CLHEP units, so that 1 kilogram equals 6.24x10^24^ MeV ns^2^ mm^-2^, 
    // and densities are in units of 6.24x^24^ MeV ns^2^ mm^-5^.

  static double GetDensityInGCM3() {return fGeo->GetCurrentNode()->GetVolume()->GetMaterial()->GetDensity()/kEdepSimDensityToGCM3; };
  static double GetPathLengthInCM() {return fGeo->GetStep() * 0.1; };

 public:
  SANDTrackerUtils(){};
  ~SANDTrackerUtils(){};
  static void Clear();
  static bool AreAdjacent(const SANDTrackerCellID &tub1, const SANDTrackerCellID &tub2);
  static const SANDTrackerCellID GetTubeID(const SANDTrackerDigitID &id);
  static const SANDTrackerPlaneID GetPlaneID(const SANDTrackerDigitID &id);
  static const SANDTrackerPlaneLocalID GetPlaneLocalID(const SANDTrackerPlaneID &id);
  static const SANDTrackerModuleID GetModelID(const SANDTrackerPlaneID &id);
  static inline std::vector<double> GetSANDInnerVolumeCenterPosition()
  {
    return std::vector<double>{sand_reco::stt::stt_center[0],
                               sand_reco::stt::stt_center[1],
                               sand_reco::stt::stt_center[2]};
  };
  // static inline double GetSANDInnerVolumeRadius() { return sand_reco::ecal::endcap::ec_r; };
  // static inline double GetSANDInnerVolumeLength() { return 2 * 1690.; };
  // static inline double GetSANDTrackerElectronDriftVelocity()
  // {
    // return sand_reco::stt::v_drift;
  // };
  static inline double GetTubeRadius() { return 2.5; };
  // static inline double GetTubeMaxDriftTime()
  // {
  //   return GetTubeRadius() / GetSANDTrackerElectronDriftVelocity();
  // };
  // static inline SANDTrackerDigitID const EncodeTubeId(const SANDTrackerPlaneID &pid,
  //                                             const SANDTrackerCellID &tid)
  // {
  //   return SANDTrackerDigitID(sand_reco::stt::encodeSTID(pid(), tid()));
  // };
  static void Init(TGeoManager *geo) {fGeo = geo;};
  static TGeoManager* GetGeoManager() {return fGeo; };
  

  static inline double GetX0(int Z, int A) {
    //https://cds.cern.ch/record/1279627/files/PH-EP-Tech-Note-2010-013.pdf
    //https://halldweb.jlab.org/DocDB/0004/000439/001/radlen.pdf
    //  The equation is an approximation while
    //  the second a result of a compact fit 
    //  to the data with an accuracy of better
    //  than 2.5% for all elements except He
    return 716.408 /*g/cm2*/ * A / (Z * (Z+1) * log(287/sqrt(Z)));
  }

  static inline double GetDEInMeV(double crossedMaterialInGCM2) {
    return crossedMaterialInGCM2 * 2. /*MeV/(g/cm2)*/;
  }

  static inline double GetDEInGeV(double crossedMaterialInGCM2) {
    return GetDEInMeV(crossedMaterialInGCM2) * 1E-3;
  }

  static inline double GetMCSSigmaAngleFromMomentumInMeV(double p, 
                                       double beta, 
                                       double pathLengthInX0) {
    // The charge of the particle is assumed
    // to be one unit of electron charge
    // std::cout << p << " "
    //           << beta << " "
    //           << pathLengthInX0 << " "
    //           << sqrt(pathLengthInX0) << " "
    //           << log(pathLengthInX0) << std::endl;
    return 13.6 /*MeV*/ / (p*beta) * sqrt(pathLengthInX0) * (1 + 0.038 * log(pathLengthInX0/(beta*beta)));
  }

  // To Do: check all units
  static inline double GetRadiusInMMToMomentumInGeVConstant() {return 0.299792458; /* GeV/(m*T) */ };
  // static inline double GetPerpMomentumInGeVFromRadiusInMM(double radius) {return GetRadiusInMMToMomentumInGeVConstant() * radius * SANDTrackANDKFGeoManager::GetMagneticField(); };
  static double GetPerpMomentumInGeVFromRadiusInMM(double radius);
  static double GetRadiusInMMFromPerpMomentumInGeV(double perpMom);
  static inline double GetMomentumInGeVFromRadiusInMM(double radius, double tanl) {return GetPerpMomentumInGeVFromRadiusInMM(radius) * sqrt(1 + tanl*tanl); };
  static inline double GetMomentumInMeVFromRadiusInMM(double radius, double tanl) {return 1.E3 * GetMomentumInGeVFromRadiusInMM(radius, tanl); };
  static inline double GetSigmaPositionMeasurement() {return 200E-6 /*m*/; };
  static inline double GetSigmaAngleMeasurement() {return 0.2 /*rad*/; };
  static double GetMagneticField() { return kMagneticFieldInT; };
  static double Getk() { return k; };
  static double Getc() { return c; };

  static TString PrintMatrix(const TMatrixD& m);
  // static TString PrintStateVector(const SANDTrackANDKFStateVector& v); 
  // static const SANDTrackerCluster* GetClusterPointer(int clusterID, const std::vector<SANDTrackerCluster>& clusters);
  static TVector3 GetCartesianCoordinateFromCylindrical(double radius, double angle, double x);


  static double GetCrossedMaterialInGCM2(double z, 
                              double px, double py, double pz,
                              double sx, double sy, double sz);
    static double GetPathLengthInX0(double z, 
                            double px, double py, double pz,
                            double sx, double sy, double sz);

    static double GetPathLengthInCM(double z, 
                            double px, double py, double pz,
                            double sx, double sy, double sz);

    static double GetDE(double z, 
                                double px, double py, double pz,
                                double sx, double sy, double sz,
                                double beta, double mass, int charge);

  friend class SANDTrackerStrawTubeTracker;
};






namespace SANDKFUtils {

TVector2 get_Bfield_perp(const TVector3& v);

double get_Bfield_long(const TVector3& v);

double get_tan_of_dip_angle(const TVector3& mom);

int get_rotation_versus(int charge);

int get_charge(int versus);

double get_direction_angle(const TVector2& dir);

double get_rotation_angle(const TVector2& dir, int charge);

double get_radius(double perp_mom);

TVector2 get_circle_center(const TVector2& momentum, const TVector2& position, int charge);

std::pair<double, double> get_circle_ys(double z, double radius, const TVector2& center);

double get_rotation_angle(double z, double y, const TVector2& center);

double get_delta_phi(double phi, double previous_phi, int versus);

double get_x(double radius, double x_0, double delta_phi, double tan_lambda, int charge);

double get_y(const TVector2& center, double radius, double phi);

double get_z(const TVector2& center, double radius, double phi);

TVector3 get_vector_momentum(double radius, double phi, double tan_lambda, int versus);

SANDKFStateVector get_state_vector(TVector3 mom, TVector3 pos, int charge);

class TrajectoryParameters;

class ParticleState {
    private:
        TVector3 _position;
        TVector3 _momentum;
    public:
        ParticleState(): _position(), _momentum() {};
        ParticleState(const TVector3& p, const TVector3& m): _position(p), _momentum(m) {};
        ParticleState(const SANDKFStateVector& vector, double z);
        TrajectoryParameters get_trajectory_parameter(int charge) const;
        SANDKFStateVector get_state_vector(int charge) const;
        const TVector3& get_position() const { return _position; };
        const TVector3& get_momentum() const { return _momentum; };
        TVector3& get_position() { return _position; };
        TVector3& get_momentum() { return _momentum; };
        void set_position(const TVector3 v) { _position = v; };
        void set_momentum(const TVector3 v) { _momentum = v; };

        ParticleState operator =(const ParticleState& p) {
            this->set_position(p.get_position());
            this->set_momentum(p.get_momentum());
            return *this;
        }
};





class TrajectoryParameters {
    public:
        double _radius;
        double _versus_of_rot; // right-hand rule in right-handed coordinate system, z == beam, y == vertical
        double _tan_lambda;
        TVector2 _center_of_rot;
        double _phi_0;
        double _x_0;
        std::pair<double, double> get_phi_pair(double z) const;
        double get_smallest_delta_phi(double z, double last_phi) const;
        std::vector<double> get_delta_phis(std::vector<double> zs) const;
        ParticleState get_particle_state(double delta_phi) const;
        std::vector<ParticleState> get_particle_states_from_delta_phi(std::vector<double> delta_phi) const;
    public:
        TrajectoryParameters() {};
        TrajectoryParameters(double r, double v, double t, const TVector2& c, double p, double x): 
            _radius(r), _versus_of_rot(v), _tan_lambda(t), _center_of_rot(c), _phi_0(p), _x_0(x) {};
        std::vector<ParticleState> get_particle_states_from_z(std::vector<double> zs) const;
};

class CovMatrixPropCheckOutput {
    public:
        double dz;
        TMatrixD initial_state_propagated;
        TMatrixD initial_covariance_propagated;
        std::vector<SANDKFStateVector> propagated_states;
        TMatrixD mean_of_propagated_states;
        TMatrixD covariance_of_propagated_states;
        TMatrixD variance_of_propagated_states;
        double distance;
        CovMatrixPropCheckOutput(): 
            initial_state_propagated(1,5),
            initial_covariance_propagated(5,5),
            mean_of_propagated_states(1,5),
            covariance_of_propagated_states(5,5),
            variance_of_propagated_states(5,5) {};
};

std::vector<SANDKFStateVector> generate_state_vectors(const SANDKFStateVector& state, const TMatrixD& cov, int n);

TMatrixD get_mean(const std::vector<SANDKFStateVector>& states);

TMatrixD get_cov(const std::vector<SANDKFStateVector>& states, const TMatrixD& mean);

TMatrixD get_var(const std::vector<SANDKFStateVector>& states, const TMatrixD& cov);

void get_mean_and_cov(const std::vector<SANDKFStateVector>& states, TMatrixD& mean, TMatrixD& cov);

using propagation = std::vector<ParticleState>;
}

#endif