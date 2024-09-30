#ifndef SANDTrackerUTILS_H
#define SANDTrackerUTILS_H

#include <vector>
#include <math.h>

#include "TTreeReader.h"
#include "TMatrixD.h"

#include "utils.h"

#include "SANDTrackerDigitCollection.h"
// #include "SANDTrackerStrawTubeTracker.h"
// #include "SANDTrackerKFKalmanFilter.h"

using SANDTrackerKFStateCovarianceMatrix = TMatrixD;
using SANDTrackerKFMeasurement = TMatrixD;

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

  
  static void InitGeo(TGeoManager *geo) { fGeo = geo; };
  // static inline const std::map<int, std::map<int, TVector2>> &GetTubeMap()
  // {
    // return sand_reco::stt::stPos;
  // };
  static std::tuple<SANDTrackerPlaneID, SANDTrackerCellID> GetPlaneAndTubeIds(
      const SANDTrackerDigitID &id);
  static std::tuple<SANDTrackerModuleID, SANDTrackerPlaneLocalID> GetModuleIdAndPlaneLocalId(
      const SANDTrackerPlaneID &id);

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
  static void Init(TGeoManager *geo);
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

  static inline double GetRadiusInMMToMomentumInGeVConstant() {return 0.299792458; /* GeV/(m*T) */ };
  // static inline double GetPerpMomentumInGeVFromRadiusInMM(double radius) {return GetRadiusInMMToMomentumInGeVConstant() * radius * SANDTrackerKFGeoManager::GetMagneticField(); };
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
  // static TString PrintStateVector(const SANDTrackerKFStateVector& v); 
  // static const SANDTrackerCluster* GetClusterPointer(int clusterID, const std::vector<SANDTrackerCluster>& clusters);
  static TVector3 GetCartesianCoordinateFromCylindrical(double radius, double angle, double x);

  friend class SANDTrackerStrawTubeTracker;
};
#endif