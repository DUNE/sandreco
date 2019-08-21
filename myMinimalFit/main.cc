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
#include <TApplication.h>
#include <TGeoManager.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TTree.h>
#include <TDatabasePDG.h>
#include <TMath.h>

#include <vector>
#include <iomanip>
#include <cmath>

#include "struct.h"


int main(int argc, char* argv[]) {


  //obtained with the following commands:
  // cd /wd/dune-it/enurec/analysis/kloe-simu
  // $$ .> myMinimalFit/ev_312_mu.txt
  // $$ tEvent->Scan("particles.pdg:Entry$:particles.tr.@digits.size():particles.tr.digits.x:particles.tr.digits.y:particles.tr.digits.z:particles.tr.digits.t","particles.pdg==13&&particles.tr.@digits.size()>200&&particles.primary==1&&z>10500&&Entry$==312")
  // $$ .>
  // x=$(tail -n +5 ev_312_mu.txt | head -n -5 | sed 's/*/ /g' | awk -F' ' '{print $6}'); for i in ${x[@]}; do echo -n "$i, "; done; echo ""
  // px = 751.35641
  // py = -669.0994
  // pz = 1104.7786
  // E  = 1497.9753
  // mass = 105.658
  // charge = -1
  
  // energy in GeV
  // distance in cm
  // magnetic field in kG
  
  // magnetic field is not read from geometry but should be provided explicitly
  // geometry should be in cm and g/cm^3 
  
  
  TFile f("/home/dune-it/data/reco/numu_geoV12_1000.0.reco.root");

  TTree* tev = (TTree*) f.Get("tEvent");
  
  event* evt = new event;
  
  tev->SetBranchAddress("event",&evt);
  
  // init geometry and mag. field
  TGeoManager::Import("geo_v12.root");
  //TApplication* gApp = new TApplication("myMinimalFit",&argc,argv);
  //TEveManager* gEve = new TEveManager(1000, 500);
  
  const int dbgLvl = 0;
  
  const double mm2cm = 0.1;
  const double kG2T = 0.1;
  const double m2cm = 100.;
  const double cm2m = 1./m2cm;
  const double GeV2MeV = 1000.;
  const double MeV2GeV = 1./GeV2MeV;
  
  const double Bx = 6.; // kGauss

  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(Bx, 0., 0.)); // kGauss
  
  genfit::MaterialEffects::getInstance()->setDebugLvl(dbgLvl);

  // init event display
  //genfit::EventDisplay* display = genfit::EventDisplay::getInstance();
  
  TDatabasePDG pdgdb;
  pdgdb.ReadPDGTable();
      
  TRandom3 rand(0);
  
  const double xsigma = 0.02;
  const double psigma = 0.05;
  
  const int nev = tev->GetEntries();
  
  int pdg = 0;
  double px, py, pz;
  double mass;
  int Z, sign;
  double x0,y0,z0,t0;
  double pt, p, e, gamma, beta;
  double bx, by, bz, bt;
  double r, period, yC, zC;
  double r_reco, yC_reco, zC_reco;
  double x0_guess, y0_guess, z0_guess;
  double px_guess, py_guess, pz_guess;
  double x0_reco, y0_reco, z0_reco;
  double px_reco, py_reco, pz_reco;
  double x0_true, y0_true, z0_true;
  double px_true, py_true, pz_true;
  double x0_old, y0_old, z0_old;
  double px_old, py_old, pz_old;
  double x0_new, y0_new, z0_new;
  double px_new, py_new, pz_new;
  double p_reco, pt_reco;
  double ey, ez;
  
  const double c = 299792458.;
  
  const bool debug = false;
  
  //TFile fout("numu_geoV12_100000.0.check.root","RECREATE");
  TFile fout(argv[3],"RECREATE");
  TTree tout("tCheck","tCheck");
  
  tout.Branch("px_old",&px_old,"px_old/D");
  tout.Branch("py_old",&py_old,"py_old/D");
  tout.Branch("pz_old",&pz_old,"pz_old/D");
  
  tout.Branch("px_new",&px_new,"px_new/D");
  tout.Branch("py_new",&py_new,"py_new/D");
  tout.Branch("pz_new",&pz_new,"pz_new/D");
  
  tout.Branch("px_true",&px_true,"px_true/D");
  tout.Branch("py_true",&py_true,"py_true/D");
  tout.Branch("pz_true",&pz_true,"pz_true/D");
  
  tout.Branch("x0_old",&x0_old,"x0_old/D");
  tout.Branch("y0_old",&y0_old,"y0_old/D");
  tout.Branch("z0_old",&z0_old,"z0_old/D");
  
  tout.Branch("x0_new",&x0_new,"x0_new/D");
  tout.Branch("y0_new",&y0_new,"y0_new/D");
  tout.Branch("z0_new",&z0_new,"z0_new/D");
  
  tout.Branch("x0_true",&x0_true,"x0_true/D");
  tout.Branch("y0_true",&y0_true,"y0_true/D");
  tout.Branch("z0_true",&z0_true,"z0_true/D");
  
  int istart = atoi(argv[1]);
  int istop  = atoi(argv[2]); 
  
  for(int i = istart; i < istop; i++)
  {
    if(debug)
      std::cout << "Event: " << i << std::endl;
  
    tev->GetEntry(i);      
    
    for(int k = 0; k < evt->particles.size(); k++)
    {    
      if(evt->particles[k].primary == 1 && evt->particles[k].pdg == 13 && evt->particles[k].has_track == true && evt->particles[k].tr.digits.size() > 0)
      {
        pdg = evt->particles[k].pdg;
    
        px = evt->particles[k].pxtrue*MeV2GeV;
        py = evt->particles[k].pytrue*MeV2GeV;
        pz = evt->particles[k].pztrue*MeV2GeV;
        
        mass = pdgdb.GetParticle(pdg)->Mass();
        Z = pdgdb.GetParticle(pdg)->Charge()/3;
        sign = Z/TMath::Abs(Z);
        
        if(debug)
        {
          std::cout << "pdg   :" << std::setw(25) << pdg << std::endl;
          std::cout << "mass  :" << std::setw(25) << mass << " (GeV)" << std::endl;
          std::cout << "charge:" << std::setw(25) << Z << std::endl;
        }
      
        x0 = evt->particles[k].tr.x0*mm2cm;
        y0 = evt->particles[k].tr.y0*mm2cm;
        z0 = evt->particles[k].tr.z0*mm2cm;
        t0 = evt->particles[k].tr.t0;
        
        px_old = evt->particles[k].pxreco*MeV2GeV;
        py_old = evt->particles[k].pyreco*MeV2GeV;
        pz_old = evt->particles[k].pzreco*MeV2GeV;
        
        x0_old = evt->particles[k].xreco*mm2cm;
        y0_old = evt->particles[k].yreco*mm2cm;
        z0_old = evt->particles[k].zreco*mm2cm;
        
        px_true = px;
        py_true = py;
        pz_true = pz;
        
        x0_true = evt->x*mm2cm;
        y0_true = evt->y*mm2cm;
        z0_true = evt->z*mm2cm;
      
        // problematic events. Check why
        if(isnan(x0) || t0 < 1.)
          continue;
      
        if(debug)
        {
          std::cout << "x0    :" << std::setw(25) << x0 << " (cm)" << std::endl;
          std::cout << "y0    :" << std::setw(25) << y0 << " (cm)" << std::endl;
          std::cout << "z0    :" << std::setw(25) << z0 << " (cm)" << std::endl;
          std::cout << "t0    :" << std::setw(25) << t0 << " (ns)" << std::endl;
        
          std::cout << "Bx    :" << std::setw(25) << Bx << " (kG)" << std::endl;
        
          std::cout << "px    :" << std::setw(25) << px << " (GeV)" << std::endl;
          std::cout << "py    :" << std::setw(25) << py << " (GeV)" << std::endl;
          std::cout << "pz    :" << std::setw(25) << pz << " (GeV)" << std::endl;
        }
      
        pt = TMath::Sqrt(py*py+pz*pz);
      
        p = TMath::Sqrt(px*px+pt*pt);
      
        e = TMath::Sqrt(p*p + mass*mass);
      
        gamma = e/mass;
      
        beta = p/e;
        
        if(debug)
        {
          std::cout << "pt    :" << std::setw(25) << pt << " (GeV)" << std::endl;
          std::cout << "p     :" << std::setw(25) << p << " (GeV)" << std::endl;
          std::cout << "e     :" << std::setw(25) << e << " (GeV)" << std::endl;
          std::cout << "gamma :" << std::setw(25) << gamma << std::endl;
          std::cout << "beta  :" << std::setw(25) << beta << std::endl;
        }
      
        bx = beta * px / p;
        by = beta * py / p;
        bz = beta * pz / p;
      
        bt = beta * pt / p;
        
        if(debug)
        {
          std::cout << "bx    :" << std::setw(25) << bx << std::endl;
          std::cout << "by    :" << std::setw(25) << by << std::endl;
          std::cout << "bz    :" << std::setw(25) << bz << std::endl;
          std::cout << "bt    :" << std::setw(25) << bt << std::endl;
        }
      
        r = pt/(c/1E9*Bx*kG2T*TMath::Abs(Z))*m2cm;
      
        period = 2*TMath::Pi()*r*cm2m/(bt*c);
      
        yC = y0 + pz/pt * r * sign;
        zC = z0 - py/pt * r * sign;
      
        if(debug)
        {
          std::cout << "r     :" << std::setw(25) << r << " (cm)" << std::endl;
          std::cout << "T     :" << std::setw(25) << period << " (s)" << std::endl;
          std::cout << "yC    :" << std::setw(25) << yC << " (cm)" << std::endl;
          std::cout << "zC    :" << std::setw(25) << zC << " (cm)" << std::endl;
        
          std::cout << "q/|p| :" << std::setw(25) << Z/p << std::endl;
          std::cout << "px/pz :" << std::setw(25) << px/pz << std::endl;
          std::cout << "py/pz :" << std::setw(25) << py/pz << std::endl;
          std::cout << "x0    :" << std::setw(25) << x0 << " (cm)" << std::endl;
          std::cout << "y0    :" << std::setw(25) << y0 << " (cm)" << std::endl;
        }
      
        // init fitter
        genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();
        fitter->setDebugLvl(dbgLvl);
        
        // start values for the fit, e.g. from pattern recognition
        x0_guess = x0+xsigma*rand.Gaus();
        y0_guess = y0+xsigma*rand.Gaus();
        z0_guess = z0+xsigma*rand.Gaus();
        px_guess = px*(1+psigma*rand.Gaus());
        py_guess = py*(1+psigma*rand.Gaus());
        pz_guess = pz*(1+psigma*rand.Gaus());
        TVector3 pos(x0_guess, y0_guess, z0_guess);
        TVector3 mom(px_guess, py_guess, pz_guess);
      
        // trackrep
        genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);
        rep->setDebugLvl(dbgLvl);
      
        // create track
        genfit::Track fitTrack(rep, pos, mom);
      
        const int detId(0); // detector ID
        int planeId(0); // detector plane ID
        int hitId(0); // hit ID
      
        double detectorResolution(0.02); // resolution of planar detectors
        TMatrixDSym hitCov(2);
        hitCov.UnitMatrix();
        hitCov *= detectorResolution*detectorResolution;
      
      
        // add some planar hits to track with coordinates I just made up
        TVectorD hitCoords(2);
        
        int n = evt->particles[k].tr.digits.size();
      
        for(unsigned int j = 0; j < n; j++)
        { 
          hitCoords[0] = evt->particles[k].tr.digits[j].x*mm2cm;
          hitCoords[1] = evt->particles[k].tr.digits[j].y*mm2cm;
          genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, nullptr);
          measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,evt->particles[k].tr.digits[j].z*mm2cm), TVector3(1,0,0), TVector3(0,1,0))), ++planeId);
          fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));
        } 
      
        //check
        fitTrack.checkConsistency();
      
        // do the fit
        try
        {
          fitter->processTrack(&fitTrack);
        }
        catch(char* a)
        {
          continue;
        }
      
        // print fit result
        p_reco = fitTrack.getFittedState().getMom().Mag();
        px_reco = fitTrack.getFittedState().getMom().X();
        py_reco = fitTrack.getFittedState().getMom().Y();
        pz_reco = fitTrack.getFittedState().getMom().Z();
        pt_reco = TMath::Sqrt(px_reco*px_reco+py_reco*py_reco+pz_reco*pz_reco);
        x0_reco = fitTrack.getFittedState().getPos().X();
        y0_reco = fitTrack.getFittedState().getPos().Y();
        z0_reco = fitTrack.getFittedState().getPos().Z();
        
        x0_new = x0_reco;
        y0_new = y0_reco;
        z0_new = z0_reco;
        px_new = px_reco;
        py_new = py_reco;
        pz_new = pz_reco;
        
        if(debug)
        {
          std::cout << "****** fit result *******" << std::endl;
          std::cout << "         " << std::setw(25) << "true" << std::setw(25) << "reco" << std::setw(25) << "residual(%)" << std::endl;
          std::cout << "q/|p|   :" << std::setw(25) << Z/p << std::setw(25) << fitTrack.getFittedState().getQop() << std::setw(25) << 100*TMath::Abs(1. - fitTrack.getFittedState().getQop()/(Z/p)) << std::endl;
          std::cout << "px/pz   :" << std::setw(25) << px/pz << std::setw(25) << px_reco/pz_reco << std::setw(25) << 100*TMath::Abs(1. - (px_reco/pz_reco)/(px/pz)) << std::endl;
          std::cout << "py/pz   :" << std::setw(25) << py/pz << std::setw(25) << py_reco/pz_reco << std::setw(25) << 100*TMath::Abs(1. - (py_reco/pz_reco)/(py/pz)) << std::endl;
          std::cout << "x0 (cm) :" << std::setw(25) << x0 << std::setw(25) << x0_reco << std::setw(25) << 100*TMath::Abs(1. - x0_reco/x0) << std::endl;
          std::cout << "y0 (cm) :" << std::setw(25) << y0 << std::setw(25) << y0_reco << std::setw(25) << 100*TMath::Abs(1. - y0_reco/y0) << std::endl;
          //fitTrack.getFittedState().Print();
          //fitTrack.Print();
          std::cout << "*************************" << std::endl;
        }
        
        if(debug)
        {
          // print particle parameters
          std::cout << "****** particle momentum *******" << std::endl;
          std::cout << "         " << std::setw(25) << "true" << std::setw(25) << "reco" << std::setw(25) << "residual(%)" << std::endl;
          std::cout << "p  (GeV):" << std::setw(25) << p << std::setw(25) << p_reco << std::setw(25) << 100*TMath::Abs(1. - p_reco/p) << std::endl;
          std::cout << "px (GeV):" << std::setw(25) << px << std::setw(25) << px_reco << std::setw(25) << 100*TMath::Abs(1. - px_reco/px) << std::endl;
          std::cout << "py (GeV):" << std::setw(25) << py << std::setw(25) << py_reco << std::setw(25) << 100*TMath::Abs(1. - py_reco/py) << std::endl;
          std::cout << "pz (GeV):" << std::setw(25) << pz << std::setw(25) << pz_reco << std::setw(25) << 100*TMath::Abs(1. - pz_reco/pz) << std::endl;
          std::cout << "*************************" << std::endl;
        }
      
        // print track parameters
        r_reco = -m2cm/(fitTrack.getFittedState().getQop()*c/1E9*Bx*kG2T);
        ez = py_reco/pt_reco;
        ey = -pz_reco/pt_reco;
        yC_reco = y0_reco + ey * r_reco;
        zC_reco = z0 + ez * r_reco;
        
        if(debug)
        {
          std::cout << "****** track parameters *******" << std::endl;
          std::cout << "         " << std::setw(25) << "true" << std::setw(25) << "reco" << std::setw(25) << "residual" << std::endl;
          std::cout << "R  (cm) :" << std::setw(25) << r << std::setw(25) << r_reco << std::setw(25) << 100*TMath::Abs(1. - r_reco/r) << std::endl;
          std::cout << "yC (cm) :" << std::setw(25) << yC << std::setw(25) << yC_reco << std::setw(25) << 100*TMath::Abs(1. - yC_reco/yC) << std::endl;
          std::cout << "zC (cm) :" << std::setw(25) << zC << std::setw(25) << zC_reco << std::setw(25) << 100*TMath::Abs(1. - zC_reco/zC) << std::endl;
          std::cout << "*************************" << std::endl;
        }
      
        //display->addEvent(&fitTrack);
        
        tout.Fill();
      
        delete fitter;
      }
    }
  }
  
  fout.cd();
  tout.Write();
  fout.Close();

  // open event display
  //display->open();

}


