#include <ConstField.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterRefTrack.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackPoint.h>

#include <MaterialEffects.h>
#include <RKTrackRep.h>
#include <TGeoMaterialInterface.h>

#include <EventDisplay.h>

#include <PlanarMeasurement.h>

#include <TEveManager.h>
#include <TGeoManager.h>
#include <TVector3.h>
#include <vector>
#include <iomanip>

#include "TDatabasePDG.h"
#include <TMath.h>


int main(int argc, char* argv[]) {

  // init geometry and mag. field
  // new TGeoManager("Geometry", "Geane geometry");
  // TGeoManager::Import("genfitGeom.root");
  TGeoManager *mng = new TGeoManager("world", "the simplest geometry");
  TGeoMaterial *mat = new TGeoMaterial("Vacuum",0,0,0);
  TGeoMedium *med = new TGeoMedium("Vacuum",1,mat);
  TGeoVolume *top = gGeoManager->MakeBox("Top",med,1.e7,1.e7,1.e7);
  gGeoManager->SetTopVolume(top);   // accessing the world
  gGeoManager->CloseGeometry();     // via a global symbol

  const int pdg = atoi(argv[4]);

  const double px = atof(argv[1]);
  const double py = atof(argv[2]);
  const double pz = atof(argv[3]);

  TDatabasePDG pdgdb;
  pdgdb.ReadPDGTable();

  const double mass = pdgdb.GetParticle(pdg)->Mass();
  const int Z = pdgdb.GetParticle(pdg)->Charge()/3;
  const int sign = Z/TMath::Abs(Z);

  std::cout << "pdg   :" << std::setw(25) << pdg << std::endl;
  std::cout << "mass  :" << std::setw(25) << mass << std::endl;
  std::cout << "charge:" << std::setw(25) << Z << std::endl;

  const double x0 = 0.;
  const double y0 = 0.;
  const double z0 = 0.;
  const double t0 = 0.;

  std::cout << "x0    :" << std::setw(25) << x0 << std::endl;
  std::cout << "y0    :" << std::setw(25) << y0 << std::endl;
  std::cout << "z0    :" << std::setw(25) << z0 << std::endl;
  std::cout << "t0    :" << std::setw(25) << t0 << std::endl;

  const double Bx = 1.; // T

  std::cout << "Bx    :" << std::setw(25) << Bx << std::endl;

  std::cout << "px    :" << std::setw(25) << px << std::endl;
  std::cout << "py    :" << std::setw(25) << py << std::endl;
  std::cout << "pz    :" << std::setw(25) << pz << std::endl;

  const double pt = TMath::Sqrt(py*py+pz*pz);

  const double p = TMath::Sqrt(px*px+pt*pt);

  const double e = TMath::Sqrt(p*p + mass*mass);

  const double gamma = e/mass;

  const double beta = p/e;

  std::cout << "pt    :" << std::setw(25) << pt << std::endl;
  std::cout << "p     :" << std::setw(25) << p << std::endl;
  std::cout << "e     :" << std::setw(25) << e << std::endl;
  std::cout << "gamma :" << std::setw(25) << gamma << std::endl;
  std::cout << "beta  :" << std::setw(25) << beta << std::endl;

  const double bx = beta * px / p;
  const double by = beta * py / p;
  const double bz = beta * pz / p;

  const double bt = beta * pt / p;

  std::cout << "bx    :" << std::setw(25) << bx << std::endl;
  std::cout << "by    :" << std::setw(25) << by << std::endl;
  std::cout << "bz    :" << std::setw(25) << bz << std::endl;
  std::cout << "bt    :" << std::setw(25) << bt << std::endl;

  const double c = 299792458.;

  const double r = pt/(c/1E9*Bx*TMath::Abs(Z));

  const double period = 2*TMath::Pi()*r/(bt*c);

  const double yC = y0 + pz/pt * r * sign;
  const double zC = z0 - py/pt * r * sign;

  std::cout << "r     :" << std::setw(25) << r << std::endl;
  std::cout << "T     :" << std::setw(25) << period << std::endl;
  std::cout << "yC    :" << std::setw(25) << yC << std::endl;
  std::cout << "zC    :" << std::setw(25) << zC << std::endl;

  std::cout << "q/|p| :" << std::setw(25) << Z/p << std::endl;
  std::cout << "px/pz :" << std::setw(25) << px/pz << std::endl;
  std::cout << "py/pz :" << std::setw(25) << py/pz << std::endl;
  std::cout << "x0    :" << std::setw(25) << x0 << std::endl;
  std::cout << "y0    :" << std::setw(25) << y0 << std::endl;

  const int sampl_n = 83;
  const double sampl_step = 0.2;
  const double sampl_zmin = 0.0;

  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(Bx * 10, 0., 0.)); // 1 T


  // init event display
  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();


  // init fitter
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();

  // start values for the fit, e.g. from pattern recognition
  TVector3 pos(x0, y0, z0);
  TVector3 mom(px, py, pz);

  // trackrep
  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

  // create track
  genfit::Track fitTrack(rep, pos, mom);


  const int detId(0); // detector ID
  int planeId(0); // detector plane ID
  int hitId(0); // hit ID

  double detectorResolution(0.1); // resolution of planar detectors
  TMatrixDSym hitCov(2);
  hitCov.UnitMatrix();
  hitCov *= detectorResolution*detectorResolution;


  // add some planar hits to track with coordinates I just made up
  TVectorD hitCoords(2);

  std::vector<double> x = {x0};
  std::vector<double> y = {y0};
  std::vector<double> z = {z0};
  std::vector<double> t = {t0};
  std::vector<double> phi = {TMath::ATan2(y.back() - yC, z.back() - zC)};

  double phi1, phi2, dt;

  for(int j = 1; j < sampl_n; j++)
  {
    z.push_back(sampl_zmin + j * sampl_step);
    y.push_back(yC - sign * TMath::Sqrt(r*r - (z.back() - zC)*(z.back() - zC)));
    phi1 = phi.back();
    phi.push_back(TMath::ATan2(y.back() - yC, z.back() - zC));
    phi2 = phi.back();
    dt = sign * (phi2 - phi1) / (2 * TMath::Pi()) * period;
    t.push_back(t.back() + dt);
    x.push_back(x.back() + bx * c * dt);

    /*
    std::cout << std::setw(25) << x.back() << 
      std::setw(25) << y.back() << 
      std::setw(25) << z.back() << 
      std::setw(25) << t.back() << std::endl;
    */
  } 

  /*
  std::vector<double> x = {0.00E+00, 4.00E+00, 8.00E+00, 1.20E+01, 1.60E+01, 2.00E+01, 2.40E+01, 2.80E+01, 3.20E+01, 3.61E+01, 4.01E+01, 4.41E+01, 4.82E+01, 5.22E+01, 5.63E+01, 6.03E+01, 6.44E+01, 6.85E+01, 7.26E+01, 7.67E+01, 8.08E+01, 8.49E+01, 8.91E+01, 9.32E+01, 9.74E+01, 1.02E+02, 1.06E+02, 1.10E+02, 1.14E+02, 1.18E+02, 1.23E+02, 1.27E+02, 1.31E+02, 1.36E+02, 1.40E+02, 1.44E+02, 1.49E+02, 1.53E+02, 1.58E+02, 1.62E+02, 1.67E+02, 1.71E+02, 1.76E+02, 1.81E+02, 1.85E+02, 1.90E+02, 1.95E+02, 2.00E+02, 2.05E+02, 2.10E+02, 2.14E+02, 2.20E+02, 2.25E+02, 2.30E+02, 2.35E+02, 2.40E+02, 2.46E+02, 2.51E+02, 2.57E+02, 2.62E+02, 2.68E+02, 2.74E+02, 2.80E+02, 2.86E+02, 2.92E+02, 2.98E+02, 3.05E+02, 3.11E+02, 3.18E+02, 3.25E+02, 3.32E+02, 3.40E+02, 3.48E+02, 3.56E+02, 3.64E+02, 3.73E+02, 3.82E+02, 3.93E+02, 4.03E+02, 4.15E+02, 4.29E+02, 4.44E+02, 4.63E+02, 5.24E+02};
  std::vector<double> y = {0.00E+00, -1.20E-01, -4.80E-01, -1.08E+00, -1.92E+00, -3.00E+00, -4.32E+00, -5.89E+00, -7.69E+00, -9.74E+00, -1.20E+01, -1.46E+01, -1.74E+01, -2.04E+01, -2.37E+01, -2.72E+01, -3.10E+01, -3.50E+01, -3.93E+01, -4.39E+01, -4.87E+01, -5.37E+01, -5.91E+01, -6.47E+01, -7.06E+01, -7.67E+01, -8.31E+01, -8.98E+01, -9.68E+01, -1.04E+02, -1.12E+02, -1.20E+02, -1.28E+02, -1.36E+02, -1.45E+02, -1.54E+02, -1.63E+02, -1.73E+02, -1.83E+02, -1.94E+02, -2.04E+02, -2.16E+02, -2.27E+02, -2.39E+02, -2.51E+02, -2.64E+02, -2.77E+02, -2.90E+02, -3.04E+02, -3.18E+02, -3.33E+02, -3.48E+02, -3.64E+02, -3.80E+02, -3.97E+02, -4.14E+02, -4.32E+02, -4.50E+02, -4.69E+02, -4.89E+02, -5.10E+02, -5.31E+02, -5.52E+02, -5.75E+02, -5.99E+02, -6.23E+02, -6.48E+02, -6.75E+02, -7.02E+02, -7.31E+02, -7.61E+02, -7.93E+02, -8.26E+02, -8.62E+02, -8.99E+02, -9.39E+02, -9.81E+02, -1.03E+03, -1.08E+03, -1.13E+03, -1.20E+03, -1.27E+03, -1.36E+03, -1.67E+03};
  std::vector<double> z = {0.00E+00, 2.00E+01, 4.00E+01, 6.00E+01, 8.00E+01, 1.00E+02, 1.20E+02, 1.40E+02, 1.60E+02, 1.80E+02, 2.00E+02, 2.20E+02, 2.40E+02, 2.60E+02, 2.80E+02, 3.00E+02, 3.20E+02, 3.40E+02, 3.60E+02, 3.80E+02, 4.00E+02, 4.20E+02, 4.40E+02, 4.60E+02, 4.80E+02, 5.00E+02, 5.20E+02, 5.40E+02, 5.60E+02, 5.80E+02, 6.00E+02, 6.20E+02, 6.40E+02, 6.60E+02, 6.80E+02, 7.00E+02, 7.20E+02, 7.40E+02, 7.60E+02, 7.80E+02, 8.00E+02, 8.20E+02, 8.40E+02, 8.60E+02, 8.80E+02, 9.00E+02, 9.20E+02, 9.40E+02, 9.60E+02, 9.80E+02, 1.00E+03, 1.02E+03, 1.04E+03, 1.06E+03, 1.08E+03, 1.10E+03, 1.12E+03, 1.14E+03, 1.16E+03, 1.18E+03, 1.20E+03, 1.22E+03, 1.24E+03, 1.26E+03, 1.28E+03, 1.30E+03, 1.32E+03, 1.34E+03, 1.36E+03, 1.38E+03, 1.40E+03, 1.42E+03, 1.44E+03, 1.46E+03, 1.48E+03, 1.50E+03, 1.52E+03, 1.54E+03, 1.56E+03, 1.58E+03, 1.60E+03, 1.62E+03, 1.64E+03, 1.67E+03};
  */

  for(unsigned int j = 0; j < sampl_n; j++)
  { 
    hitCoords[0] = x[j]*100;
    hitCoords[1] = y[j]*100;
    genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, nullptr);
    measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,z[j]*100), TVector3(1,0,0), TVector3(0,1,0))), ++planeId);
    fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));
  } 

  //check
  fitTrack.checkConsistency();

  // do the fit
  fitter->processTrack(&fitTrack);

  // print fit result
  std::cout << "****** fit result *******" << std::endl;
  std::cout << "         " << std::setw(25) << "true" << std::setw(25) << "reco" << std::setw(25) << "residual" << std::endl;
  std::cout << "q/|p|   :" << std::setw(25) << Z/p << std::setw(25) << fitTrack.getFittedState().getQop() << std::setw(25) << Z/p - fitTrack.getFittedState().getQop() << std::endl;
  std::cout << "px/pz   :" << std::setw(25) << px/pz << std::setw(25) << fitTrack.getFittedState().getMom().X()/fitTrack.getFittedState().getMom().Z() << std::setw(25) << px/pz - fitTrack.getFittedState().getMom().X()/fitTrack.getFittedState().getMom().Z() << std::endl;
  std::cout << "py/pz   :" << std::setw(25) << py/pz << std::setw(25) << fitTrack.getFittedState().getMom().Y()/fitTrack.getFittedState().getMom().Z() << std::setw(25) << py/pz - fitTrack.getFittedState().getMom().Y()/fitTrack.getFittedState().getMom().Z() << std::endl;
  std::cout << "x0      :" << std::setw(25) << x0 << std::setw(25) << fitTrack.getFittedState().getPos().X() << std::setw(25) << x0 - fitTrack.getFittedState().getPos().X() << std::endl;
  std::cout << "y0      :" << std::setw(25) << y0 << std::setw(25) << fitTrack.getFittedState().getPos().Y() << std::setw(25) << y0 - fitTrack.getFittedState().getPos().Y() << std::endl;
  //fitTrack.getFittedState().Print();
  //fitTrack.Print();
  std::cout << "*************************" << std::endl;

  //check
  fitTrack.checkConsistency();


  display->addEvent(&fitTrack);


  delete fitter;

  // open event display
  display->open();

}


