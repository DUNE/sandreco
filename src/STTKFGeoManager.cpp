#include "STTKFGeoManager.h"
#include "SANDTrackerUtils.h"

TGeoManager* STTKFGeoManager::fGeo = 0;
const double STTKFGeoManager::kEdepSimDensityToGCM3 = 6.24E18;
// https://github.com/ClarkMcGrew/edep-sim/blob/master/README.md#reading-the-output
// Be aware that in the saved TGeoManager object, the masses and densities are 
// also in CLHEP units, so that 1 kilogram equals 6.24x10^24^ MeV ns^2^ mm^-2^, 
// and densities are in units of 6.24x^24^ MeV ns^2^ mm^-5^.

void STTKFGeoManager::Init() {fGeo = SANDTrackerUtils::GetGeoManager(); }

double STTKFGeoManager::GetCrossedMaterialInGCM2(double z, 
                                           double px, double py, double pz,
                                           double sx, double sy, double sz) {
  // sx, sy, sz represent the directional cosine 
  // and should be correctly normalized:
  // sx^2 + sy^2 + sz^2 = 1
  // The function return the ammount of
  // crossed material in g/cm2
  fGeo->InitTrack(px, py, pz, sx, sy, sz);
  const double* lastPosition = 0;

  double crossedMaterial = 0.;
  // std::cout << __LINE__ << std::endl;
  // std::cout << sx << " " << sy << " " << sz << std::endl;
  // std::cout << fGeo->GetCurrentPoint()[2] << " " << z << std::endl;
  while((lastPosition = fGeo->GetCurrentPoint()) && lastPosition[2] > z) {
    auto density = GetDensityInGCM3();
    auto pathLength = GetPathLengthInCM();
    crossedMaterial += density * pathLength;
    // std::cout << density << " " << pathLength << std::endl;
    fGeo->Step();
  }
  return crossedMaterial;
}

double STTKFGeoManager::GetPathLengthInX0(double z, 
                                           double px, double py, double pz,
                                           double sx, double sy, double sz) {

  fGeo->InitTrack(px, py, pz, sx, sy, sz);
  const double* lastPosition = 0;

  double pathLengthInX0 = 0.;
  int count = 0;
  while((lastPosition = fGeo->GetCurrentPoint()) && lastPosition[2] > z) {
    auto Z = static_cast<int>(fGeo->GetCurrentNode()->GetVolume()->GetMaterial()->GetZ());
    auto A = static_cast<int>(fGeo->GetCurrentNode()->GetVolume()->GetMaterial()->GetA());
    auto name = static_cast<std::string>(fGeo->GetCurrentNode()->GetVolume()->GetMaterial()->GetName());
    count++;
    auto X0 = SANDTrackerUtils::GetX0(Z, A);
    auto density = GetDensityInGCM3();
    auto pathLength = GetPathLengthInCM();
    pathLengthInX0 += pathLength * density / X0;
    fGeo->Step();
  }
  return pathLengthInX0;
}

double STTKFGeoManager::GetPathLengthInCM(double z, 
                                           double px, double py, double pz,
                                           double sx, double sy, double sz) {
  // sx, sy, sz represent the directional cosine 
  // and should be correctly normalized:
  // sx^2 + sy^2 + sz^2 = 1
  // The function return the ammount of
  // crossed material in g/cm2
  fGeo->InitTrack(px, py, pz, sx, sy, sz);
  const double* lastPosition = 0;

  double pathLengthInCM = 0.;
  while((lastPosition = fGeo->GetCurrentPoint()) && lastPosition[2] > z) {
    auto pathLength = GetPathLengthInCM();
    pathLengthInCM += pathLength;
    fGeo->Step();
  }
  return pathLengthInCM;
}


double STTKFGeoManager::GetDE(double z, 
                              double px, double py, double pz,
                              double sx, double sy, double sz,
                              double beta, double mass, int charge) {

  double K = 0.307075; // MeV mol−1 cm2
  double m_e = 0.5109989461; // MeV  
  double gamma = 1 / sqrt(1 - beta*beta);
  double W_max = (2 * m_e * beta * beta * gamma* gamma) /
                 (1 + 2 * gamma * m_e / mass + pow(m_e / mass, 2));

  fGeo->InitTrack(px, py, pz, sx, sy, sz);
  const double* lastPosition = 0;

  double dE = 0.;
  while((lastPosition = fGeo->GetCurrentPoint()) && lastPosition[2] > z) {
    fGeo->FindNextBoundary();
    auto density = GetDensityInGCM3();
    auto pathLength = GetPathLengthInCM();

    auto current_material = fGeo->GetCurrentNode()->GetVolume()->GetMaterial();
    auto name = static_cast<std::string>(current_material->GetName());
    // std::cout << " name " << name << " lastPosition[2] " << lastPosition[2] << std::endl;
    // std::cout << " name " << fGeo->GetCurrentNode()->GetVolume()->GetName() << std::endl;

    int n_elements = current_material->GetNelements();
    double A = 0;
    double Z = 0;
    double w = 0;
    double material_I = 0;
    double num = 0;
    double den = 0;
    for (int i = 0; i < n_elements; i++) {
      current_material->GetElementProp(A, Z, w, i);
      int I = 0;
      if (Z < 13) {
        I = 12 * Z + 7;
      } else {
        I = 9.76 * Z + 58.8 * pow(Z, -1.19) * Z;
      }

      num += w * Z / A * log(I);
      den += w * Z / A;

    }

    material_I = exp(num / den) * 1E-6;
    
    // std::cout << "density " << density << std::endl;
    // std::cout << "pathLength " << pathLength << std::endl;
    // std::cout << "material_I " << material_I << std::endl;

    auto matZ = static_cast<double>(current_material->GetZ());
    auto matA = static_cast<double>(current_material->GetA());

    double plasma_energy = sqrt(density*fabs(matZ/matA))*28.816E-6;
    double delta = log(plasma_energy/material_I)+log(beta*gamma)-0.5;     
    double coeff = K * charge * charge * matZ / matA / (beta * beta);
    double BB_log = log(2 * m_e * beta * beta * gamma * gamma * W_max / (material_I * material_I));
    double dEdx = coeff * (0.5 * BB_log - beta * beta - 0.5*delta);

    dE += dEdx * density * pathLength;
    fGeo->Step(true, true);
  }
  return dE;
}
