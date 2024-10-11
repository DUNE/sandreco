#ifndef STTKFGEOMANAGER_H
#define STTKFGEOMANAGER_H

#include "TGeoManager.h"

class STTKFGeoManager {
  private:
    static const double kEdepSimDensityToGCM3;
    // https://github.com/ClarkMcGrew/edep-sim/blob/master/README.md#reading-the-output
    // Be aware that in the saved TGeoManager object, the masses and densities are 
    // also in CLHEP units, so that 1 kilogram equals 6.24x10^24^ MeV ns^2^ mm^-2^, 
    // and densities are in units of 6.24x^24^ MeV ns^2^ mm^-5^.

    static TGeoManager* fGeo; // o affini

    static double GetDensityInGCM3() {return fGeo->GetCurrentNode()->GetVolume()->GetMaterial()->GetDensity()/kEdepSimDensityToGCM3; };
    static double GetPathLengthInCM() {return fGeo->GetStep() * 0.1; };

  public:
    // STTKFGeoManager() {fSTT.Init(); fGeo = STTUtils::GetGeoManager(); };
    STTKFGeoManager() {};
    static void Init();

    // static double GetNextZ(const STTPlaneID& plane) { return fSTT.Next(plane, 1)->GetZ(); };
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

};

#endif