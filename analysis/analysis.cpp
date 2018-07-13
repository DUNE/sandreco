#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TGeoManager.h>
#include <TLorentzVector.h>
#include <TGeoNavigator.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TGeoTube.h>
#include <TChain.h>

#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4Event.h"
#include "/mnt/nas01/users/mtenti/sw/edep-sim/edep-sim-bin/include/edep-sim/TG4HitSegment.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

const double angres = 0.00;
const double angresint = 0.00;

std::vector<std::string> volInAC {"volArgonCubeDetector",
    "volLArCryo",
    "volArgonCube",
    "volArgonCubeWall",
    "volLArModule",
    "volLArBot",
    "volLArTop",
    "volLArMod",
    "volPixelPlaneL",
    "volPixelPlaneR",
    "volLArCathode",
    "volLArSubModuleLeft",
    "volLArSubModuleRight",
    "volResistiveFieldTop",
    "volResistiveFieldSide",
    "volResistivePlaneTop",
    "volActiveTopWall",
    "volResistiveWireTop",
    "volResistiveWireSide",
    "volLArLight",
    "volLightUPlane",
    "volLightDPlane",
    "volLArActive",
    "volCylinderLBox",
    "volCylinderTBox",
    "volCylinderRBox",
    "volCylinder",
    "volCylinderT"};

std::vector<std::string> volInKLOE {"KLOEYokeBarrel",
    "KLOEYokeEndcapAL",
    "KLOEYokeEndcapAR",
    "KLOEYokeEndcapBL",
    "KLOEYokeEndcapBR",
    "KLOEYokeEndcapCL",
    "KLOEYokeEndcapCR",
    "KLOEYokeEndcapDL",
    "KLOEYokeEndcapDR",
    "KLOESolenoidCryostatEndcapL",
    "KLOESolenoidCryostatEndcapR",
    "KLOESolenoidCryostatInnerWall",
    "KLOESolenoidCryostatOuterWall",
    "KLOESolenoidCoilShell",
    "KLOESolenoidCoil",
    "KLOEBarrelECAL",
    "KLOEEndcapECALL",
    "KLOEEndcapECALR",
    "KLOETrackingRegion",
    "KLOEEXTTRK",
    "KLOEINTTRK",
    "PassiveLayer",
    "volKLOESTTFULLNEWCONF_STTModuleFULL",
    "volSTTPlane1FULL",
    "volSTTPlane2FULL",
    "STTModule_Radiator_1A",
    "STTModule_Radiator_1B",
    "STTModule_Radiator_2A",
    "STTModule_Radiator_2B",
    "StrawTubeFULL_for_STTPlane1FULL",
    "StrawTubeFULL_for_STTPlane2FULL",
    "StrawTubeFULL_for_STTPlane1FULL_0_0_wire",
    "StrawTubeFULL_for_STTPlane1FULL_0_0_straw",
    "StrawTubeFULL_for_STTPlane2FULL_0_0_wire",
    "StrawTubeFULL_for_STTPlane2FULL_0_0_straw"};

bool isInAc(std::string vol)
{
    for(std::vector<std::string>::iterator it=volInAC.begin(); it != volInAC.end(); ++it)
    {
        if(vol.find(*it)!=std::string::npos)
        {
            return true;
        }
    }
    return false;
}

bool isInKLOE(std::string vol)
{
    for(std::vector<std::string>::iterator it=volInKLOE.begin(); it != volInKLOE.end(); ++it)
    {
        if(vol.find(*it)!=std::string::npos)
        {
            return true;
        }
    }
    return false;
}

double mindist(double s1x, double s1y, double s1z,
               double s2x, double s2y, double s2z,
               double px, double py, double pz)
{
    double segmod = (s1x - s2x)*(s1x - s2x)
                    +(s1y - s2y)*(s1y - s2y)
                    +(s1z - s2z)*(s1z - s2z);

    double prod = (px - s1x)*(s2x - s1x)
                  +(py - s1y)*(s2y - s1y)
                  +(pz - s1z)*(s2z - s1z);

    double t = std::min(std::max(prod/segmod,0.),1.);

    double s3x = s1x + (s2x - s1x) * t;
    double s3y = s1y + (s2y - s1y) * t;
    double s3z = s1z + (s2z - s1z) * t;
    /*
    	std::cout << s1x << " " <<
    				s1y << " " <<
    				s1z << " " <<
    				s2x << " " <<
    				s2y << " " <<
    				s2z << " " <<
    				px << " " <<
    				py << " " <<
    				pz << " " <<
    				s3x << " " <<
    				s3y << " " <<
    				s3z << " " <<
    				sqrt((px - s3x)*(px - s3x)+(py - s3y)*(py - s3y)+(pz - s3z)*(pz - s3z)) << " " << std::endl;
    */
    return sqrt((px - s3x)*(px - s3x)
                +(py - s3y)*(py - s3y)
                +(pz - s3z)*(pz - s3z));
}

double angle(double x1, double y1, double z1,
             double x2, double y2, double z2)
{
    double prod = x1*x2+y1*y2+z1*z2;
    double mag1 = sqrt(x1*x1+y1*y1+z1*z1);
    double mag2 = sqrt(x2*x2+y2*y2+z2*z2);

    return TMath::ACos(prod/(mag1*mag2));
}

bool ishitok2(TG4Event* ev, int trackid,
              double x, double y, double z,
              double tol = 5.)
{
    std::vector<double> vec;
    for(unsigned int jj = 0; jj < ev->Trajectories.size(); jj++)
    {
        if(ev->Trajectories[jj].TrackId == trackid)
        {
            for(unsigned int kk = 0; kk < ev->Trajectories[jj].Points.size()-1; kk++)
            {
                vec.push_back(mindist(ev->Trajectories[jj].Points[kk].Position.X(),
                                      ev->Trajectories[jj].Points[kk].Position.Y(),
                                      ev->Trajectories[jj].Points[kk].Position.Z(),
                                      ev->Trajectories[jj].Points[kk+1].Position.X(),
                                      ev->Trajectories[jj].Points[kk+1].Position.Y(),
                                      ev->Trajectories[jj].Points[kk+1].Position.Z(),
                                      x,y,z));
            }
        }
    }
    if(*std::min_element(vec.begin(),vec.end()) > tol)
        return false;
    else
        return true;
}

bool ishitok(TG4Event* ev, int trackid, TG4HitSegment hit,
             double postol = 5., double angtol = 0.3)
{
    double x = 0.5*(hit.Start.X()+hit.Stop.X());
    double y = 0.5*(hit.Start.Y()+hit.Stop.Y());
    double z = 0.5*(hit.Start.Z()+hit.Stop.Z());

    std::vector<double> dpos;
    std::vector<double> dang;
    for(unsigned int jj = 0; jj < ev->Trajectories.size(); jj++)
    {
        if(ev->Trajectories[jj].TrackId == trackid)
        {
            for(unsigned int kk = 0; kk < ev->Trajectories[jj].Points.size()-1; kk++)
            {
                dpos.push_back(mindist(ev->Trajectories[jj].Points[kk].Position.X(),
                                       ev->Trajectories[jj].Points[kk].Position.Y(),
                                       ev->Trajectories[jj].Points[kk].Position.Z(),
                                       ev->Trajectories[jj].Points[kk+1].Position.X(),
                                       ev->Trajectories[jj].Points[kk+1].Position.Y(),
                                       ev->Trajectories[jj].Points[kk+1].Position.Z(),
                                       x,y,z));
                dang.push_back(angle(ev->Trajectories[jj].Points[kk+1].Position.X()-ev->Trajectories[jj].Points[kk].Position.X(),
                                     ev->Trajectories[jj].Points[kk+1].Position.Y()-ev->Trajectories[jj].Points[kk].Position.Y(),
                                     ev->Trajectories[jj].Points[kk+1].Position.Z()-ev->Trajectories[jj].Points[kk].Position.Z(),
                                     hit.Stop.X()-hit.Start.X(),
                                     hit.Stop.Y()-hit.Start.Y(),
                                     hit.Stop.Z()-hit.Start.Z()));
            }
        }
    }
    int index = std::distance(dpos.begin(), std::min_element(dpos.begin(),dpos.end()));

    //std::cout << dpos[index] << " " << dang[index] << std::endl;

    if(dpos[index] > postol || dang[index] > angtol)
        return false;
    else
        return true;
}

bool strcompare (std::string a, std::string b) {
    return (a.compare(b)==0);
}

double getpar(double x1, double x2, double x0)
{
  return (x0 - x1)/(x2 - x1);
}

double getX(double t, double x1, double x2)
{
  return x1 + (x2 - x1) * t;
}

bool isPointOk1(double max1, double max2, double min1, double min2, double x1, double x2)
{
  return (x1 < max1 && x1 > min1) && (x2 < max2 && x2 > min2);
}

bool isPointOk2(double x1, double x2, double x0)
{
  return (x1 <= x0 && x0 <= x2)||(x2 <= x0 && x0 <= x1);
}

bool intersec(double x1, double x2, 
              double y1, double y2, 
              double z1, double z2, 
              double x0, double& y0, double& z0,
              double ymin, double zmin, 
              double ymax, double zmax)
{
  double t = getpar(x1, x2, x0);
  y0 = getX(t, y1, y2);
  z0 = getX(t, z1, z2);
  
  //std::cout << t << " " << x0 << " " << x1 << " " << x2 << " " << y0 << " " << y1 << " " << y2 << " " << z0 << " " << z1 << " " << z2 << " " << ymin << " " << ymax << " " << zmin << " " << zmax << std::endl;
  
  return isPointOk1(ymax, zmax, ymin, zmin, y0, z0) && isPointOk2(x1, x2, x0);
}

double getL(double x1, double y1, double z1,
            double x2, double y2, double z2,
            double yc, double zc, double rmin, double rmax)
{
  double a0 = - (y2 - y1);
  double b0 = - (z1 - z2);
  double c0 = z1*y2 - z2*y1;
  
  double a1 = -2*zc;
  double b1 = -2*yc;
  double c1min = pow(zc,2) + pow(yc,2) - pow(rmin,2);
  double c1max = pow(zc,2) + pow(yc,2) - pow(rmax,2);
  
  double a = pow(a0,2) + pow(b0,2);
  double b = 2*c0*a0 + a1*pow(b0,2) - a0*b1*b0;
  double cmin = pow(c0,2) - c0*b1*b0 + c1min*pow(b0,2);
  double cmax = pow(c0,2) - c0*b1*b0 + c1max*pow(b0,2);
  
  double pz1 = (-b - sqrt(pow(b,2)-4*a*cmin))/(2*a);
  double pz2 = (-b - sqrt(pow(b,2)-4*a*cmax))/(2*a);
  double py1 = (-a0*pz1 - c0)/b0;
  double py2 = (-a0*pz2 - c0)/b0;
  double px1 = ((x2 - x1)*pz1 - (z1*x2 - z2*x1))/(z2 - z1);
  double px2 = ((x2 - x1)*pz2 - (z1*x2 - z2*x1))/(z2 - z1);
  
  return sqrt(pow(px1-px2,2)+pow(py1-py2,2)+pow(pz1-pz2,2));
}

void classify(const char* finname, const char* foutname)
{
    //TString fname = "/mnt/nas01/users/mtenti/wd/prod/edep-sim/ArgonCube-RT5-petti.0.001-60gev-1E7/muon-ArgonCube-RT5-petti.0.001-60gev-1E7-CC.1.edep-sim.root";

    //TString fname = "/mnt/nas01/users/mtenti/wd/prod/edep-sim/ArgonCube-RT3.5-petti.1E7/muon-ArgonCube-RT3.5-petti.1E7-CC.1.edep-sim.root";

    TString fname = finname;

    //TFile f(fname);
	TChain* t = new TChain("EDepSimEvents","EDepSimEvents");
	t->Add(finname);
	TFile f(t->GetListOfFiles()->At(0)->GetTitle());
    //TTree* t = (TTree*) f.Get("EDepSimEvents");
    TGeoManager* geo = (TGeoManager*) f.Get("EDepSimGeometry");
    TGeoNavigator* nav = geo->GetCurrentNavigator();
    nav->cd("volWorld_PV/volDetEnclosure_PV_0/volKLOEFULLECALSENSITIVE_EXTTRK_NEWGAP_PV_0/KLOETrackingRegion_volume_PV_0");
    
    
    TChain* tgenie = new TChain("DetSimPassThru/gRooTracker");
	  tgenie->Add(finname);

    double kloe_center[] = {0, 0, 0};
    double origin[] = {0, 0, 0};

    nav->LocalToMaster(origin,kloe_center);

    std::cout << "kloe center: " << kloe_center[0] << " " << kloe_center[1] << " " << kloe_center[2] << std::endl;

    TGeoTube* tub = (TGeoTube*) geo->GetTopVolume()->GetNode(0)->GetVolume()->GetNode(1)->GetVolume()->GetNode(42)->GetVolume()->GetShape();
    
    const double XL = 0. - 3042.;
    const double XR = 0. + 3042.;
    const double YU = 4330.96 - 1468. + 2457.6;
    const double YD = 4330.96 - 1468. - 2457.6;
    const double ZB = 7396.1625 - 4256. - 3545.5;
    const double ZF = 7396.1625 - 4256. + 3545.5;

    double rmin = tub->GetRmin();
    double rmax = tub->GetRmax();
    double diron = 50.;
    double dair = 50.;
    double dyoke = 370.;
    
    double radiimax[] = {rmax - dair, rmax - 2*dair - diron, rmax - 3*dair - 2*diron};
    double radiimin[] = {rmax - dair - diron, rmax - 2*dair - 2*diron, rmax - 3*dair - 2*diron - dyoke };

    double sampling_radius[] = {rmax + 5.,
                                rmax - 1 * dair - 0.5 * diron,
                                rmax - 2 * dair - 1.5 * diron,
                                rmin - 5.
                               };

    double RsamplExtTrk = 3590.;
    double RsamplTrkReg = 2000.;
    double RsamplCalo = 2230.;
    double RsamplRmaxYoke = 3300.;
    double RsamplRminYoke = 2930.;

    TG4Event* ev = new TG4Event;

    t->SetBranchAddress("Event",&ev);

    double StdHepP4[2000][4];
    double EvtWght;
    TObjString* EvtCode = new TObjString;

    tgenie->SetBranchAddress("StdHepP4",StdHepP4);
    tgenie->SetBranchAddress("EvtWght",&EvtWght);
    tgenie->SetBranchAddress("EvtCode",&EvtCode);

    int eventid;
    int muonid = -999;
    std::vector<std::string> vols;
    std::string vname;

    std::vector<double> xhitxtr;
    std::vector<double> yhitxtr;
    std::vector<double> zhitxtr;
    std::vector<double> thitxtr;
    std::vector<double> Lhitxtr;

    std::vector<double> x1hitxtr;
    std::vector<double> y1hitxtr;
    std::vector<double> z1hitxtr;
    std::vector<double> x2hitxtr;
    std::vector<double> y2hitxtr;
    std::vector<double> z2hitxtr;

    double xhititr;
    double yhititr;
    double zhititr;
    double thititr;

    double x1hititr;
    double y1hititr;
    double z1hititr;
    double x2hititr;
    double y2hititr;
    double z2hititr;

    double dextL;
    double dexttheta;
    double muextptreco;
    double dintL;
    double dinttheta;
    double muintptreco;
    double muairX1;
    double muairY1;
    double muairZ1;
    double muairPX1;
    double muairPY1;
    double muairPZ1;
    double muairX2;
    double muairY2;
    double muairZ2;
    double muairPX2;
    double muairPY2;
    double muairPZ2;
    double mutrkregX1;
    double mutrkregY1;
    double mutrkregZ1;
    double mutrkregPX1;
    double mutrkregPY1;
    double mutrkregPZ1;
    double mutrkregX2;
    double mutrkregY2;
    double mutrkregZ2;
    double mutrkregPX2;
    double mutrkregPY2;
    double mutrkregPZ2;
    double muRmaxYokeX1;
    double muRmaxYokeY1;
    double muRmaxYokeZ1;
    double muRmaxYokePX1;
    double muRmaxYokePY1;
    double muRmaxYokePZ1;
    double muRmaxYokeX2;
    double muRmaxYokeY2;
    double muRmaxYokeZ2;
    double muRmaxYokePX2;
    double muRmaxYokePY2;
    double muRmaxYokePZ2;
    double muRminYokeX1;
    double muRminYokeY1;
    double muRminYokeZ1;
    double muRminYokePX1;
    double muRminYokePY1;
    double muRminYokePZ1;
    double muRminYokeX2;
    double muRminYokeY2;
    double muRminYokeZ2;
    double muRminYokePX2;
    double muRminYokePY2;
    double muRminYokePZ2;
    double muCaloX1;
    double muCaloY1;
    double muCaloZ1;
    double muCaloPX1;
    double muCaloPY1;
    double muCaloPZ1;
    double muCaloX2;
    double muCaloY2;
    double muCaloZ2;
    double muCaloPX2;
    double muCaloPY2;
    double muCaloPZ2;
    
    double muLArDE;
    double muLArDE2;
    double muLArL;
    double muLArSX;
    double muLArSY;
    double muLArSZ;
    
    std::vector<double>  muLArHitX1;
    std::vector<double>  muLArHitY1;
    std::vector<double>  muLArHitZ1;
    std::vector<double>  muLArHitT1;
    std::vector<double>  muLArHitX2;
    std::vector<double>  muLArHitY2;
    std::vector<double>  muLArHitZ2;
    std::vector<double>  muLArHitT2;
    std::vector<double>  muLArHitSX;
    std::vector<double>  muLArHitSY;
    std::vector<double>  muLArHitSZ;
    std::vector<double>  muLArHitDE;
    std::vector<double>  muLArHitDE2;
    std::vector<double>  muLArHitL;    

    TRandom3 rand;

    bool exited;
    bool entered;
    bool hitcal;
    bool hitxtr;
    bool hititr;
    bool hitstt;
    bool isCC;

    int nhitcal;
    int nhitxtr;
    int nhititr;
    int nhitstt;
    
    double dummyX, dummyY, dummyZ;

    double mupx, mupy, mupz;
    double vx, vy, vz;
    double weight;
    double nuEne;

    bool samplExtTrk = false;
    bool samplTrkReg = false;
    bool samplRminYoke = false;
    bool samplRmaxYoke = false;
    bool samplCalo = false;
    
    bool exitedUp = false;
    bool exitedDw = false;
    bool exitedLf = false;
    bool exitedRg = false;
    bool exitedFr = false;
    bool exitedBk = false;
    
    double exitX = -999.;
    double exitY = -999.;
    double exitZ = -999.;
    double exitPX = -999.;
    double exitPY = -999.;
    double exitPZ = -999.;
    
    bool found = false;

    TFile fout(foutname,"RECREATE");
    TTree tout("classes","classes");
    tout.Branch("eventid",&eventid,"eventid/I");
    tout.Branch("exitedUp",&exitedUp,"exitedUp/O");
    tout.Branch("exitedDw",&exitedDw,"exitedDw/O");
    tout.Branch("exitedLf",&exitedLf,"exitedLf/O");
    tout.Branch("exitedRg",&exitedRg,"exitedRg/O");
    tout.Branch("exitedFr",&exitedFr,"exitedFr/O");
    tout.Branch("exitedBk",&exitedBk,"exitedBk/O");
    tout.Branch("exited",&exited,"exited/O");
    tout.Branch("entered",&entered,"entered/O");
    tout.Branch("hitcal",&hitcal,"hitcal/O");
    tout.Branch("hitxtr",&hitxtr,"hitxtr/O");
    tout.Branch("hititr",&hititr,"hititr/O");
    tout.Branch("hitstt",&hitstt,"hitstt/O");
    tout.Branch("nhitstt",&nhitstt,"nhitstt/I");
    tout.Branch("nhitxtr",&nhitxtr,"nhitxtr/I");
    tout.Branch("nhititr",&nhititr,"nhititr/I");
    tout.Branch("nhitcal",&nhitcal,"nhitcal/I");
    tout.Branch("xhitxtr",&xhitxtr);
    tout.Branch("yhitxtr",&yhitxtr);
    tout.Branch("zhitxtr",&zhitxtr);
    tout.Branch("thitxtr",&thitxtr);
    tout.Branch("Lhitxtr",&Lhitxtr);
    tout.Branch("exitX",&exitX,"exitX/D");
    tout.Branch("exitY",&exitY,"exitY/D");
    tout.Branch("exitZ",&exitZ,"exitZ/D");
    tout.Branch("exitPX",&exitPX,"exitPX/D");
    tout.Branch("exitPY",&exitPY,"exitPY/D");
    tout.Branch("exitPZ",&exitPZ,"exitPZ/D");
    tout.Branch("xhititr",&xhititr,"xhititr/D");
    tout.Branch("yhititr",&yhititr,"yhititr/D");
    tout.Branch("zhititr",&zhititr,"zhititr/D");
    tout.Branch("thititr",&thititr,"thititr/D");
    tout.Branch("x1hititr",&x1hititr,"x1hititr/D");
    tout.Branch("y1hititr",&y1hititr,"y1hititr/D");
    tout.Branch("x2hititr",&x2hititr,"x2hititr/D");
    tout.Branch("y2hititr",&y2hititr,"y2hititr/D");
    tout.Branch("z1hititr",&z1hititr,"z1hititr/D");
    tout.Branch("z2hititr",&z2hititr,"z2hititr/D");
    tout.Branch("x1hitxtr",&x1hitxtr);
    tout.Branch("y1hitxtr",&y1hitxtr);
    tout.Branch("z1hitxtr",&z1hitxtr);
    tout.Branch("x2hitxtr",&x2hitxtr);
    tout.Branch("y2hitxtr",&y2hitxtr);
    tout.Branch("z2hitxtr",&z2hitxtr);
    tout.Branch("dextL",&dextL,"dextL/D");
    tout.Branch("dexttheta",&dexttheta,"dexttheta/D");
    tout.Branch("muextptreco",&muextptreco,"muextptreco/D");
    tout.Branch("dintL",&dintL,"dintL/D");
    tout.Branch("dinttheta",&dinttheta,"dinttheta/D");
    tout.Branch("muintptreco",&muintptreco,"muintptreco/D");
    tout.Branch("muairX1",&muairX1,"muairX1/D");
    tout.Branch("muairY1",&muairY1,"muairY1/D");
    tout.Branch("muairZ1",&muairZ1,"muairZ1/D");
    tout.Branch("muairPX1",&muairPX1,"muairPX1/D");
    tout.Branch("muairPY1",&muairPY1,"muairPY1/D");
    tout.Branch("muairPZ1",&muairPZ1,"muairPZ1/D");
    tout.Branch("muairX2",&muairX2,"muairX2/D");
    tout.Branch("muairY2",&muairY2,"muairY2/D");
    tout.Branch("muairZ2",&muairZ2,"muairZ2/D");
    tout.Branch("muairPX2",&muairPX2,"muairPX2/D");
    tout.Branch("muairPY2",&muairPY2,"muairPY2/D");
    tout.Branch("muairPZ2",&muairPZ2,"muairPZ2/D");
    tout.Branch("mutrkregX1",&mutrkregX1,"mutrkregX1/D");
    tout.Branch("mutrkregY1",&mutrkregY1,"mutrkregY1/D");
    tout.Branch("mutrkregZ1",&mutrkregZ1,"mutrkregZ1/D");
    tout.Branch("mutrkregPX1",&mutrkregPX1,"mutrkregPX1/D");
    tout.Branch("mutrkregPY1",&mutrkregPY1,"mutrkregPY1/D");
    tout.Branch("mutrkregPZ1",&mutrkregPZ1,"mutrkregPZ1/D");
    tout.Branch("mutrkregX2",&mutrkregX2,"mutrkregX2/D");
    tout.Branch("mutrkregY2",&mutrkregY2,"mutrkregY2/D");
    tout.Branch("mutrkregZ2",&mutrkregZ2,"mutrkregZ2/D");
    tout.Branch("mutrkregPX2",&mutrkregPX2,"mutrkregPX2/D");
    tout.Branch("mutrkregPY2",&mutrkregPY2,"mutrkregPY2/D");
    tout.Branch("mutrkregPZ2",&mutrkregPZ2,"mutrkregPZ2/D");
    tout.Branch("muRmaxYokeX1",&muRmaxYokeX1,"muRmaxYokeX1/D");
    tout.Branch("muRmaxYokeY1",&muRmaxYokeY1,"muRmaxYokeY1/D");
    tout.Branch("muRmaxYokeZ1",&muRmaxYokeZ1,"muRmaxYokeZ1/D");
    tout.Branch("muRmaxYokePX1",&muRmaxYokePX1,"muRmaxYokePX1/D");
    tout.Branch("muRmaxYokePY1",&muRmaxYokePY1,"muRmaxYokePY1/D");
    tout.Branch("muRmaxYokePZ1",&muRmaxYokePZ1,"muRmaxYokePZ1/D");
    tout.Branch("muRmaxYokeX2",&muRmaxYokeX2,"muRmaxYokeX2/D");
    tout.Branch("muRmaxYokeY2",&muRmaxYokeY2,"muRmaxYokeY2/D");
    tout.Branch("muRmaxYokeZ2",&muRmaxYokeZ2,"muRmaxYokeZ2/D");
    tout.Branch("muRmaxYokePX2",&muRmaxYokePX2,"muRmaxYokePX2/D");
    tout.Branch("muRmaxYokePY2",&muRmaxYokePY2,"muRmaxYokePY2/D");
    tout.Branch("muRmaxYokePZ2",&muRmaxYokePZ2,"muRmaxYokePZ2/D");
    tout.Branch("muRminYokeX1",&muRminYokeX1,"muRminYokeX1/D");
    tout.Branch("muRminYokeY1",&muRminYokeY1,"muRminYokeY1/D");
    tout.Branch("muRminYokeZ1",&muRminYokeZ1,"muRminYokeZ1/D");
    tout.Branch("muRminYokePX1",&muRminYokePX1,"muRminYokePX1/D");
    tout.Branch("muRminYokePY1",&muRminYokePY1,"muRminYokePY1/D");
    tout.Branch("muRminYokePZ1",&muRminYokePZ1,"muRminYokePZ1/D");
    tout.Branch("muRminYokeX2",&muRminYokeX2,"muRminYokeX2/D");
    tout.Branch("muRminYokeY2",&muRminYokeY2,"muRminYokeY2/D");
    tout.Branch("muRminYokeZ2",&muRminYokeZ2,"muRminYokeZ2/D");
    tout.Branch("muRminYokePX2",&muRminYokePX2,"muRminYokePX2/D");
    tout.Branch("muRminYokePY2",&muRminYokePY2,"muRminYokePY2/D");
    tout.Branch("muRminYokePZ2",&muRminYokePZ2,"muRminYokePZ2/D"); 
    tout.Branch("muCaloX1",&muCaloX1,"muCaloX1/D");
    tout.Branch("muCaloY1",&muCaloY1,"muCaloY1/D");
    tout.Branch("muCaloZ1",&muCaloZ1,"muCaloZ1/D");
    tout.Branch("muCaloPX1",&muCaloPX1,"muCaloPX1/D");
    tout.Branch("muCaloPY1",&muCaloPY1,"muCaloPY1/D");
    tout.Branch("muCaloPZ1",&muCaloPZ1,"muCaloPZ1/D");
    tout.Branch("muCaloX2",&muCaloX2,"muCaloX2/D");
    tout.Branch("muCaloY2",&muCaloY2,"muCaloY2/D");
    tout.Branch("muCaloZ2",&muCaloZ2,"muCaloZ2/D");
    tout.Branch("muCaloPX2",&muCaloPX2,"muCaloPX2/D");
    tout.Branch("muCaloPY2",&muCaloPY2,"muCaloPY2/D");
    tout.Branch("muCaloPZ2",&muCaloPZ2,"muCaloPZ2/D");  
    tout.Branch("muLArDE",&muLArDE,"muLArDE/D");
    tout.Branch("muLArDE2",&muLArDE2,"muLArDE2/D");
    tout.Branch("muLArL",&muLArL,"muLArL/D");
    tout.Branch("muLArSX",&muLArSX,"muLArSX/D");
    tout.Branch("muLArSY",&muLArSY,"muLArSY/D");
    tout.Branch("muLArSZ",&muLArSZ,"muLArSZ/D");
    tout.Branch("muLArHitX1",&muLArHitX1);
    tout.Branch("muLArHitY1",&muLArHitY1);
    tout.Branch("muLArHitZ1",&muLArHitZ1);
    tout.Branch("muLArHitT1",&muLArHitT1);
    tout.Branch("muLArHitX2",&muLArHitX2);
    tout.Branch("muLArHitY2",&muLArHitY2);
    tout.Branch("muLArHitZ2",&muLArHitZ2);
    tout.Branch("muLArHitT2",&muLArHitT2);
    tout.Branch("muLArHitSX",&muLArHitSX);
    tout.Branch("muLArHitSY",&muLArHitSY);
    tout.Branch("muLArHitSZ",&muLArHitSZ);
    tout.Branch("muLArHitDE",&muLArHitDE);
    tout.Branch("muLArHitDE2",&muLArHitDE2);
    tout.Branch("muLArHitL",&muLArHitL);
    tout.Branch("vx",&vx,"vx/D");
    tout.Branch("vy",&vy,"vy/D");
    tout.Branch("vz",&vz,"vz/D");
    tout.Branch("mupx",&mupx,"mupx/D");
    tout.Branch("mupy",&mupy,"mupy/D");
    tout.Branch("mupz",&mupz,"mupz/D");
    tout.Branch("weight",&weight,"weight/D");
    tout.Branch("nuEne",&nuEne,"nuEne/D");
    tout.Branch("isCC",&isCC,"isCC/O");
    
    const int nev = t->GetEntries();
    
    std::cout << "Events: " << nev << std::endl;
    std::cout << std::setw(4) << int(0) << "%" << std::flush;
    
    t->AddFriend(tgenie);

    for(int i = 0; i < nev; i++)
    {
        std::cout << "\b\b\b\b\b" << std::setw(4) << int((i+1.)/nev*100) << "%" << std::flush;
        
        //if(i % 1000 == 0)
        //    std::cout << std::setw(7) << i << std::endl;
    
    
        t->GetEntry(i);

        eventid = ev->EventId;
        muonid = -999;
        vols.clear();
        
        found = false;

        xhitxtr.clear();
        yhitxtr.clear();
        zhitxtr.clear();
        thitxtr.clear();
        Lhitxtr.clear();

        x1hitxtr.clear();
        y1hitxtr.clear();
        z1hitxtr.clear();
        x2hitxtr.clear();
        y2hitxtr.clear();
        z2hitxtr.clear();

        xhititr = -999.;
        yhititr = -999.;
        zhititr = -999.;
        thititr = -999.;

        x1hititr = -999.;
        x2hititr = -999.;
        y1hititr = -999.;
        y2hititr = -999.;
        z1hititr = -999.;
        z2hititr = -999.;

        exited = false;
        exitedUp = false;
        exitedDw = false;
        exitedLf = false;
        exitedRg = false;
        exitedFr = false;
        exitedBk = false;
        entered = false;
        
        exitX = -999;
        exitY = -999;
        exitZ = -999;
        exitPX = -999;
        exitPY = -999;
        exitPZ = -999;

        hitcal = false;
        hitxtr = false;
        hititr = false;
        hitstt = false;

        mupx = -999.;
        mupy = -999.;
        mupz = -999.;

        muairX1 = -999.;
        muairY1 = -999.;
        muairZ1 = -999.;

        muairX2 = -999.;
        muairY2 = -999.;
        muairZ2 = -999.;

        muairPX1 = -999.;
        muairPY1 = -999.;
        muairPZ1 = -999.;

        muairPX2 = -999.;
        muairPY2 = -999.;
        muairPZ2 = -999.;

        mutrkregX1 = -999.;
        mutrkregY1 = -999.;
        mutrkregZ1 = -999.;

        mutrkregX2 = -999.;
        mutrkregY2 = -999.;
        mutrkregZ2 = -999.;

        mutrkregPX1 = -999.;
        mutrkregPY1 = -999.;
        mutrkregPZ1 = -999.;

        mutrkregPX2 = -999.;
        mutrkregPY2 = -999.;
        mutrkregPZ2 = -999.;

        muRmaxYokeX1 = -999.;
        muRmaxYokeY1 = -999.;
        muRmaxYokeZ1 = -999.;

        muRmaxYokeX2 = -999.;
        muRmaxYokeY2 = -999.;
        muRmaxYokeZ2 = -999.;

        muRmaxYokePX1 = -999.;
        muRmaxYokePY1 = -999.;
        muRmaxYokePZ1 = -999.;

        muRmaxYokePX2 = -999.;
        muRmaxYokePY2 = -999.;
        muRmaxYokePZ2 = -999.;

        muRminYokeX1 = -999.;
        muRminYokeY1 = -999.;
        muRminYokeZ1 = -999.;

        muRminYokeX2 = -999.;
        muRminYokeY2 = -999.;
        muRminYokeZ2 = -999.;

        muRminYokePX1 = -999.;
        muRminYokePY1 = -999.;
        muRminYokePZ1 = -999.;

        muRminYokePX2 = -999.;
        muRminYokePY2 = -999.;
        muRminYokePZ2 = -999.;
        
        muCaloX1 = -999.;
        muCaloY1 = -999.;
        muCaloZ1 = -999.;

        muCaloX2 = -999.;
        muCaloY2 = -999.;
        muCaloZ2 = -999.;

        muCaloPX1 = -999.;
        muCaloPY1 = -999.;
        muCaloPZ1 = -999.;

        muCaloPX2 = -999.;
        muCaloPY2 = -999.;
        muCaloPZ2 = -999.;

        samplExtTrk = false;
        samplTrkReg = false;
        samplRminYoke = false;
        samplRmaxYoke = false;
        samplCalo = false;

        vx = ev->Primaries[0].Position.X();
        vy = ev->Primaries[0].Position.Y();
        vz = ev->Primaries[0].Position.Z();

        nhitcal = 0;
        nhitxtr = 0;
        nhititr = 0;
        nhitstt = 0;

        dextL = -999.;
        dexttheta = -999.;
        dintL = -999.;
        dinttheta = -999.;
        muextptreco = -999.;
        muintptreco = -999.;
        
        muLArDE = 0.;
        muLArDE2 = 0.;
        muLArL = 0.;
        muLArSX = 0.;
        muLArSY = 0.;
        muLArSZ = 0.;
    
        muLArHitX1.clear();
        muLArHitY1.clear();
        muLArHitZ1.clear();
        muLArHitT1.clear();
        muLArHitX2.clear();
        muLArHitY2.clear();
        muLArHitZ2.clear();
        muLArHitT2.clear();
        muLArHitSX.clear();
        muLArHitSY.clear();
        muLArHitSZ.clear();
        muLArHitDE.clear();
        muLArHitDE2.clear();
        muLArHitL.clear();

        std::vector<double> hits_r[3];
        std::vector<TG4HitSegment> hits[3];

        std::vector<double> hitsint_r;
        std::vector<TG4HitSegment> hitsint;

        weight = EvtWght;
        nuEne = StdHepP4[0][3];
        if(EvtCode->GetString().Contains("[CC]") == true)
        {
            isCC = true;
        }
        else
        {
            isCC = false;
        }

        for(unsigned int j = 0; j < ev->Primaries[0].Particles.size(); j++)
        {
            if(ev->Primaries[0].Particles[j].PDGCode == 13)
            {
                muonid = ev->Primaries[0].Particles[j].TrackId;
                mupx = ev->Primaries[0].Particles[j].Momentum.X();
                mupy = ev->Primaries[0].Particles[j].Momentum.Y();
                mupz = ev->Primaries[0].Particles[j].Momentum.Z();
                break;
            }
        }

        for(unsigned int j = 0; j < ev->Trajectories.size(); j++)
        {
            if(ev->Trajectories[j].TrackId == muonid)
            {
                for(unsigned int k = 0; k < ev->Trajectories[j].Points.size(); k++)
                {
                    TLorentzVector pos = ev->Trajectories[j].Points[k].Position;
                    TGeoNode* node = nav->FindNode(pos.X(),pos.Y(),pos.Z());
                    if(node == 0)
                    {
                        vname = "no volume";
                    }
                    else
                    {
                        vname = node->GetName();
                    }

                    if(k > 0 )
                    {

                        double r1 = sqrt((kloe_center[1] - ev->Trajectories[j].Points[k-1].Position.Y())*
                                         (kloe_center[1] - ev->Trajectories[j].Points[k-1].Position.Y())+
                                         (kloe_center[2] - ev->Trajectories[j].Points[k-1].Position.Z())*
                                         (kloe_center[2] - ev->Trajectories[j].Points[k-1].Position.Z()));

                        double r2 = sqrt((kloe_center[1] - ev->Trajectories[j].Points[k].Position.Y())*
                                         (kloe_center[1] - ev->Trajectories[j].Points[k].Position.Y())+
                                         (kloe_center[2] - ev->Trajectories[j].Points[k].Position.Z())*
                                         (kloe_center[2] - ev->Trajectories[j].Points[k].Position.Z()));

                        //std::cout << r1 << " " << r2 << " " << RsamplExtTrk << " " << RsamplTrkReg << std::endl;
                        /*
                        if(eventid == 20 || eventid == 58 || eventid == 70 || eventid == 75 || eventid == 82 || eventid == 99)
                        {
                            std::cout << eventid << " " << ev->Trajectories[j].Points[k-1].Position.X() << " " << ev->Trajectories[j].Points[k].Position.X() << " " << ev->Trajectories[j].Points[k-1].Position.Y() << " " << ev->Trajectories[j].Points[k].Position.Y() << " " << ev->Trajectories[j].Points[k-1].Position.Z() << " " << ev->Trajectories[j].Points[k].Position.Z() << std::endl;
                        }
                        */
                        if(found == false && intersec(ev->Trajectories[j].Points[k-1].Position.X(), ev->Trajectories[j].Points[k].Position.X(), 
                                    ev->Trajectories[j].Points[k-1].Position.Y(), ev->Trajectories[j].Points[k].Position.Y(),
                                    ev->Trajectories[j].Points[k-1].Position.Z(), ev->Trajectories[j].Points[k].Position.Z(),
                                    XL, dummyY, dummyZ,
                                    YD, ZB, YU, ZF))
                        {
                            exitedLf = true;
                            exitX = XL;
                            exitY = dummyY;
                            exitZ = dummyZ;
                            exitPX = ev->Trajectories[j].Points[k].Momentum.X();
                            exitPY = ev->Trajectories[j].Points[k].Momentum.Y();
                            exitPZ = ev->Trajectories[j].Points[k].Momentum.Z();
                            found = true;
                        }
                        
                        if(found == false && intersec(ev->Trajectories[j].Points[k-1].Position.X(), ev->Trajectories[j].Points[k].Position.X(), 
                                    ev->Trajectories[j].Points[k-1].Position.Y(), ev->Trajectories[j].Points[k].Position.Y(),
                                    ev->Trajectories[j].Points[k-1].Position.Z(), ev->Trajectories[j].Points[k].Position.Z(),
                                    XR, dummyY, dummyZ,
                                    YD, ZB, YU, ZF))
                        {
                            exitedRg = true;
                            exitX = XR;
                            exitY = dummyY;
                            exitZ = dummyZ;
                            exitPX = ev->Trajectories[j].Points[k].Momentum.X();
                            exitPY = ev->Trajectories[j].Points[k].Momentum.Y();
                            exitPZ = ev->Trajectories[j].Points[k].Momentum.Z();
                            found = true;
                        }
                        
                        if(found == false && intersec(ev->Trajectories[j].Points[k-1].Position.Y(), ev->Trajectories[j].Points[k].Position.Y(), 
                                    ev->Trajectories[j].Points[k-1].Position.X(), ev->Trajectories[j].Points[k].Position.X(),
                                    ev->Trajectories[j].Points[k-1].Position.Z(), ev->Trajectories[j].Points[k].Position.Z(),
                                    YD, dummyX, dummyZ,
                                    XL, ZB, XR, ZF))
                        {
                            exitedDw = true;
                            exitX = dummyX;
                            exitY = YD;
                            exitZ = dummyZ;
                            exitPX = ev->Trajectories[j].Points[k].Momentum.X();
                            exitPY = ev->Trajectories[j].Points[k].Momentum.Y();
                            exitPZ = ev->Trajectories[j].Points[k].Momentum.Z();
                            found = true;
                        }
                        
                        if(found == false && intersec(ev->Trajectories[j].Points[k-1].Position.Y(), ev->Trajectories[j].Points[k].Position.Y(), 
                                    ev->Trajectories[j].Points[k-1].Position.X(), ev->Trajectories[j].Points[k].Position.X(),
                                    ev->Trajectories[j].Points[k-1].Position.Z(), ev->Trajectories[j].Points[k].Position.Z(),
                                    YU, dummyX, dummyZ,
                                    XL, ZB, XR, ZF))
                        {
                            exitedUp = true;
                            exitX = dummyX;
                            exitY = YU;
                            exitZ = dummyZ;
                            exitPX = ev->Trajectories[j].Points[k].Momentum.X();
                            exitPY = ev->Trajectories[j].Points[k].Momentum.Y();
                            exitPZ = ev->Trajectories[j].Points[k].Momentum.Z();
                            found = true;
                        }
                        
                        if(found == false && intersec(ev->Trajectories[j].Points[k-1].Position.Z(), ev->Trajectories[j].Points[k].Position.Z(), 
                                    ev->Trajectories[j].Points[k-1].Position.X(), ev->Trajectories[j].Points[k].Position.X(),
                                    ev->Trajectories[j].Points[k-1].Position.Y(), ev->Trajectories[j].Points[k].Position.Y(),
                                    ZB, dummyX, dummyY,
                                    XL, YD, XR, YU))
                        {
                            exitedBk = true;
                            exitX = dummyX;
                            exitY = dummyY;
                            exitZ = ZB;
                            exitPX = ev->Trajectories[j].Points[k].Momentum.X();
                            exitPY = ev->Trajectories[j].Points[k].Momentum.Y();
                            exitPZ = ev->Trajectories[j].Points[k].Momentum.Z();
                            found = true;
                        }
                        
                        if(found == false && intersec(ev->Trajectories[j].Points[k-1].Position.Z(), ev->Trajectories[j].Points[k].Position.Z(), 
                                    ev->Trajectories[j].Points[k-1].Position.X(), ev->Trajectories[j].Points[k].Position.X(),
                                    ev->Trajectories[j].Points[k-1].Position.Y(), ev->Trajectories[j].Points[k].Position.Y(),
                                    ZF, dummyX, dummyY,
                                    XL, YD, XR, YU))
                        {
                            exitedFr = true;
                            exitX = dummyX;
                            exitY = dummyY;
                            exitZ = ZF;
                            exitPX = ev->Trajectories[j].Points[k].Momentum.X();
                            exitPY = ev->Trajectories[j].Points[k].Momentum.Y();
                            exitPZ = ev->Trajectories[j].Points[k].Momentum.Z();
                            found = true;
                        }

                        if(!samplExtTrk && r1 > RsamplExtTrk && r2 < RsamplExtTrk)
                        {
                            muairPX1 = ev->Trajectories[j].Points[k-1].Momentum.X();
                            muairPY1 = ev->Trajectories[j].Points[k-1].Momentum.Y();
                            muairPZ1 = ev->Trajectories[j].Points[k-1].Momentum.Z();
                            muairPX2 = ev->Trajectories[j].Points[k].Momentum.X();
                            muairPY2 = ev->Trajectories[j].Points[k].Momentum.Y();
                            muairPZ2 = ev->Trajectories[j].Points[k].Momentum.Z();
                            muairX1 = ev->Trajectories[j].Points[k-1].Position.X();
                            muairY1 = ev->Trajectories[j].Points[k-1].Position.Y();
                            muairZ1 = ev->Trajectories[j].Points[k-1].Position.Z();
                            muairX2 = ev->Trajectories[j].Points[k].Position.X();
                            muairY2 = ev->Trajectories[j].Points[k].Position.Y();
                            muairZ2 = ev->Trajectories[j].Points[k].Position.Z();

                            samplExtTrk = true;
                        }

                        if(!samplTrkReg && r1 > RsamplTrkReg && r2 < RsamplTrkReg)
                        {
                            mutrkregPX1 = ev->Trajectories[j].Points[k-1].Momentum.X();
                            mutrkregPY1 = ev->Trajectories[j].Points[k-1].Momentum.Y();
                            mutrkregPZ1 = ev->Trajectories[j].Points[k-1].Momentum.Z();
                            mutrkregPX2 = ev->Trajectories[j].Points[k].Momentum.X();
                            mutrkregPY2 = ev->Trajectories[j].Points[k].Momentum.Y();
                            mutrkregPZ2 = ev->Trajectories[j].Points[k].Momentum.Z();
                            mutrkregX1 = ev->Trajectories[j].Points[k-1].Position.X();
                            mutrkregY1 = ev->Trajectories[j].Points[k-1].Position.Y();
                            mutrkregZ1 = ev->Trajectories[j].Points[k-1].Position.Z();
                            mutrkregX2 = ev->Trajectories[j].Points[k].Position.X();
                            mutrkregY2 = ev->Trajectories[j].Points[k].Position.Y();
                            mutrkregZ2 = ev->Trajectories[j].Points[k].Position.Z();

                            samplTrkReg = true;
                        }

                        if(!samplRmaxYoke && r1 > RsamplRmaxYoke && r2 < RsamplRmaxYoke)
                        {
                            muRmaxYokePX1 = ev->Trajectories[j].Points[k-1].Momentum.X();
                            muRmaxYokePY1 = ev->Trajectories[j].Points[k-1].Momentum.Y();
                            muRmaxYokePZ1 = ev->Trajectories[j].Points[k-1].Momentum.Z();
                            muRmaxYokePX2 = ev->Trajectories[j].Points[k].Momentum.X();
                            muRmaxYokePY2 = ev->Trajectories[j].Points[k].Momentum.Y();
                            muRmaxYokePZ2 = ev->Trajectories[j].Points[k].Momentum.Z();
                            muRmaxYokeX1 = ev->Trajectories[j].Points[k-1].Position.X();
                            muRmaxYokeY1 = ev->Trajectories[j].Points[k-1].Position.Y();
                            muRmaxYokeZ1 = ev->Trajectories[j].Points[k-1].Position.Z();
                            muRmaxYokeX2 = ev->Trajectories[j].Points[k].Position.X();
                            muRmaxYokeY2 = ev->Trajectories[j].Points[k].Position.Y();
                            muRmaxYokeZ2 = ev->Trajectories[j].Points[k].Position.Z();

                            samplRmaxYoke = true;
                        }

                        if(!samplRminYoke && r1 > RsamplRminYoke && r2 < RsamplRminYoke)
                        {
                            muRminYokePX1 = ev->Trajectories[j].Points[k-1].Momentum.X();
                            muRminYokePY1 = ev->Trajectories[j].Points[k-1].Momentum.Y();
                            muRminYokePZ1 = ev->Trajectories[j].Points[k-1].Momentum.Z();
                            muRminYokePX2 = ev->Trajectories[j].Points[k].Momentum.X();
                            muRminYokePY2 = ev->Trajectories[j].Points[k].Momentum.Y();
                            muRminYokePZ2 = ev->Trajectories[j].Points[k].Momentum.Z();
                            muRminYokeX1 = ev->Trajectories[j].Points[k-1].Position.X();
                            muRminYokeY1 = ev->Trajectories[j].Points[k-1].Position.Y();
                            muRminYokeZ1 = ev->Trajectories[j].Points[k-1].Position.Z();
                            muRminYokeX2 = ev->Trajectories[j].Points[k].Position.X();
                            muRminYokeY2 = ev->Trajectories[j].Points[k].Position.Y();
                            muRminYokeZ2 = ev->Trajectories[j].Points[k].Position.Z();

                            samplRminYoke = true;
                        }

                        if(!samplCalo && r1 > RsamplCalo && r2 < RsamplCalo)
                        {
                            muCaloPX1 = ev->Trajectories[j].Points[k-1].Momentum.X();
                            muCaloPY1 = ev->Trajectories[j].Points[k-1].Momentum.Y();
                            muCaloPZ1 = ev->Trajectories[j].Points[k-1].Momentum.Z();
                            muCaloPX2 = ev->Trajectories[j].Points[k].Momentum.X();
                            muCaloPY2 = ev->Trajectories[j].Points[k].Momentum.Y();
                            muCaloPZ2 = ev->Trajectories[j].Points[k].Momentum.Z();
                            muCaloX1 = ev->Trajectories[j].Points[k-1].Position.X();
                            muCaloY1 = ev->Trajectories[j].Points[k-1].Position.Y();
                            muCaloZ1 = ev->Trajectories[j].Points[k-1].Position.Z();
                            muCaloX2 = ev->Trajectories[j].Points[k].Position.X();
                            muCaloY2 = ev->Trajectories[j].Points[k].Position.Y();
                            muCaloZ2 = ev->Trajectories[j].Points[k].Position.Z();

                            samplCalo = true;
                        }

                    }
                    //std::cout << vname.c_str() << std::endl;
                    vols.push_back(vname);
                }

                std::vector<std::string>::iterator it = std::unique (vols.begin(), vols.end(),strcompare);
                vols.resize(std::distance(vols.begin(),it));

                //std::cout << "********************" << std::endl;
                //for(unsigned int k = 0; k < vols.size(); k++)
                //	std::cout << vols[k].c_str() << std::endl;

                for(unsigned int k = 0; k < vols.size() - 1; k++)
                {
                    if(isInAc(vols[k]) && !isInAc(vols[k+1]))
                    {
                        exited = true;
                        //std::cout << i << "\t" << vols[k].c_str() << " " << vols[k+1].c_str() << std::endl;
                    }
                    if(isInKLOE(vols[k]) && !entered)
                    {
                        entered = true;
                        //std::cout << i << "\t" << vols[k].c_str() << std::endl;
                    }
                }
            }
        }

        for (std::map<std::string,std::vector<TG4HitSegment> >::iterator it=ev->SegmentDetectors.begin();
                it!=ev->SegmentDetectors.end(); ++it)
        {
            if(it->first == "StrawTracker")
                for(unsigned int j = 0; j < it->second.size(); j++)
                {
                    for(unsigned int k = 0; k < it->second[j].Contrib.size(); k++)
                    {
                        if(it->second[j].Contrib[k] == muonid)
                        {   /*
                            if(ishitok2(ev,muonid,
                            		0.5*(it->second[j].Start.X()+it->second[j].Stop.X()),
                            		0.5*(it->second[j].Start.Y()+it->second[j].Stop.Y()),
                            		0.5*(it->second[j].Start.Z()+it->second[j].Stop.Z())))*/
                            if(ishitok(ev,muonid,it->second[j]))
                            {
                                hitstt = true;
                                nhitstt++;
                                break;
                            }
                        }
                    }
                }

            if(it->first == "InternalTracker")
                for(unsigned int j = 0; j < it->second.size(); j++)
                {
                    for(unsigned int k = 0; k < it->second[j].Contrib.size(); k++)
                    {
                        if(it->second[j].Contrib[k] == muonid)
                        {
                            if(ishitok(ev,muonid,it->second[j]))
                            {
                                hititr = true;
                                double z = 0.5*(it->second[j].Start.Z()+it->second[j].Stop.Z());
                                double y = 0.5*(it->second[j].Start.Y()+it->second[j].Stop.Y());
                                double r = sqrt((y - kloe_center[1])*(y - kloe_center[1])+(z - kloe_center[2])*(z - kloe_center[2]));

                                hitsint_r.push_back(r);
                                hitsint.push_back(it->second[j]);

                                break;
                            }
                        }
                    }
                }

            if(it->first == "ExternalTracker")
                for(unsigned int j = 0; j < it->second.size(); j++)
                {
                    for(unsigned int k = 0; k < it->second[j].Contrib.size(); k++)
                    {
                        if(it->second[j].Contrib[k] == muonid)
                        {   /*
                            if(ishitok2(ev,muonid,
                                       0.5*(it->second[j].Start.X()+it->second[j].Stop.X()),
                                       0.5*(it->second[j].Start.Y()+it->second[j].Stop.Y()),
                                       0.5*(it->second[j].Start.Z()+it->second[j].Stop.Z())))*/
                            if(ishitok(ev,muonid,it->second[j]))
                            {
                                hitxtr = true;
                                double z = 0.5*(it->second[j].Start.Z()+it->second[j].Stop.Z());
                                double y = 0.5*(it->second[j].Start.Y()+it->second[j].Stop.Y());
                                double r = sqrt((y - kloe_center[1])*(y - kloe_center[1])+(z - kloe_center[2])*(z - kloe_center[2]));

                                if(r < sampling_radius[0] && r > sampling_radius[1])
                                {
                                    hits_r[0].push_back(r);
                                    hits[0].push_back(it->second[j]);
                                }
                                else if(r < sampling_radius[1] && r > sampling_radius[2])
                                {
                                    hits_r[1].push_back(r);
                                    hits[1].push_back(it->second[j]);
                                }
                                else if(r < sampling_radius[2] && r > sampling_radius[3])
                                {
                                    hits_r[2].push_back(r);
                                    hits[2].push_back(it->second[j]);
                                }
                                //std::cout << i << " " << z << " " << y << " " << theta << std::endl;
                                break;
                            }
                        }
                    }
                }

            if(it->first == "EMCal")
                for(unsigned int j = 0; j < it->second.size(); j++)
                {
                    for(unsigned int k = 0; k < it->second[j].Contrib.size(); k++)
                    {
                        if(it->second[j].Contrib[k] == muonid)
                        {   /*
                            if(ishitok2(ev,muonid,
                                       0.5*(it->second[j].Start.X()+it->second[j].Stop.X()),
                                       0.5*(it->second[j].Start.Y()+it->second[j].Stop.Y()),
                                       0.5*(it->second[j].Start.Z()+it->second[j].Stop.Z())))*/
                            if(ishitok(ev,muonid,it->second[j]))
                            {
                                hitcal = true;
                                nhitcal++;
                                break;
                            }
                        }
                    }
                }

            if(it->first == "ArgonCube")
                for(unsigned int j = 0; j < it->second.size(); j++)
                {
                    for(unsigned int k = 0; k < it->second[j].Contrib.size(); k++)
                    {
                        if(it->second[j].Contrib[k] == muonid)
                        { 
                            if(ishitok(ev,muonid,it->second[j]))
                            {
                                muLArDE += it->second[j].EnergyDeposit;
                                muLArDE2 += it->second[j].SecondaryDeposit;
                                muLArL += it->second[j].TrackLength;
                                
                                muLArHitX1.push_back(it->second[j].Start.X());
                                muLArHitY1.push_back(it->second[j].Start.Y());
                                muLArHitZ1.push_back(it->second[j].Start.Z());
                                muLArHitT1.push_back(it->second[j].Start.T());
                                
                                muLArHitX2.push_back(it->second[j].Stop.X());
                                muLArHitY2.push_back(it->second[j].Stop.Y());
                                muLArHitZ2.push_back(it->second[j].Stop.Z());
                                muLArHitT2.push_back(it->second[j].Stop.T());
                                
                                double mag = sqrt(pow(it->second[j].Stop.X()-it->second[j].Start.X(),2)+
                                                  pow(it->second[j].Stop.Y()-it->second[j].Start.Y(),2)+
                                                  pow(it->second[j].Stop.Z()-it->second[j].Start.Z(),2));
                                
                                muLArHitSX.push_back((it->second[j].Stop.X()-it->second[j].Start.X())/mag);
                                muLArHitSY.push_back((it->second[j].Stop.Y()-it->second[j].Start.Y())/mag);
                                muLArHitSZ.push_back((it->second[j].Stop.Z()-it->second[j].Start.Z())/mag);
                                
                                muLArHitDE.push_back(it->second[j].EnergyDeposit);
                                muLArHitDE2.push_back(it->second[j].SecondaryDeposit);
                                muLArHitL.push_back(it->second[j].TrackLength);
                                break;
                            }
                        }
                    }
                }

        }
        
        int nn = 0;
        int nmax = 10;
        
        for(int i = muLArHitSX.size() - 1; i >= 0; i--)
        {
            if(nn > nmax)
              break;
            
            muLArSX += muLArHitSX[i];
            muLArSY += muLArHitSY[i];
            muLArSZ += muLArHitSZ[i];
            
            nn++;
        }
        
        muLArSX /= nn;
        muLArSY /= nn;
        muLArSZ /= nn;

        if(hitsint_r.size() > 0)
        {
            int minindex = std::distance(hitsint_r.begin(), std::min_element(hitsint_r.begin(),hitsint_r.end()));
            int maxindex = std::distance(hitsint_r.begin(), std::max_element(hitsint_r.begin(),hitsint_r.end()));

            x1hititr = hitsint[maxindex].Start.X();
            y1hititr = hitsint[maxindex].Start.Y();
            z1hititr = hitsint[maxindex].Start.Z();
            x2hititr = hitsint[minindex].Stop.X();
            y2hititr = hitsint[minindex].Stop.Y();
            z2hititr = hitsint[minindex].Stop.Z();

            double theta = TMath::ATan((hitsint[maxindex].Start.Y()-hitsint[minindex].Stop.Y())/(hitsint[maxindex].Start.Z()-hitsint[minindex].Stop.Z()));

            theta += rand.Gaus(0,angresint);

            xhititr = 0.5*(hitsint[maxindex].Start.X()+hitsint[minindex].Stop.X());
            zhititr = 0.5*(hitsint[maxindex].Start.Z()+hitsint[minindex].Stop.Z());
            yhititr = 0.5*(hitsint[maxindex].Start.Y()+hitsint[minindex].Stop.Y());

            thititr = theta;

            nhititr = 1;

        }

        for(int l = 0; l < 3; l++)
        {
            if(hits_r[l].size() == 0)
                break;

            int minindex = std::distance(hits_r[l].begin(), std::min_element(hits_r[l].begin(),hits_r[l].end()));
            int maxindex = std::distance(hits_r[l].begin(), std::max_element(hits_r[l].begin(),hits_r[l].end()));

            x1hitxtr.push_back(hits[l][maxindex].Start.X());
            y1hitxtr.push_back(hits[l][maxindex].Start.Y());
            z1hitxtr.push_back(hits[l][maxindex].Start.Z());
            x2hitxtr.push_back(hits[l][minindex].Stop.X());
            y2hitxtr.push_back(hits[l][minindex].Stop.Y());
            z2hitxtr.push_back(hits[l][minindex].Stop.Z());

            double theta = TMath::ATan((hits[l][maxindex].Start.Y()-hits[l][minindex].Stop.Y())/
                                       (hits[l][maxindex].Start.Z()-hits[l][minindex].Stop.Z()));

            theta += rand.Gaus(0,angres);

            double x = 0.5*(hits[l][maxindex].Start.X()+hits[l][minindex].Stop.X());
            double z = 0.5*(hits[l][maxindex].Start.Z()+hits[l][minindex].Stop.Z());
            double y = 0.5*(hits[l][maxindex].Start.Y()+hits[l][minindex].Stop.Y());
            
            double L = getL(hits[l][maxindex].Start.X(), hits[l][maxindex].Start.Y(), hits[l][maxindex].Start.Z(), 
                          hits[l][minindex].Stop.X(), hits[l][minindex].Stop.Y(), hits[l][minindex].Stop.Z(), 
                          kloe_center[1], kloe_center[2], radiimin[l], radiimax[l]);

            xhitxtr.push_back(x);
            yhitxtr.push_back(y);
            zhitxtr.push_back(z);
            thitxtr.push_back(theta);
            Lhitxtr.push_back(L);
            
        }

        nhitxtr = y1hitxtr.size();

        if(nhitxtr>1)
        {       
            
            dextL = sqrt((yhitxtr[0]-yhitxtr[yhitxtr.size()-1])*(yhitxtr[0]-yhitxtr[yhitxtr.size()-1])+
                      (zhitxtr[0]-zhitxtr[zhitxtr.size()-1])*(zhitxtr[0]-zhitxtr[zhitxtr.size()-1]));

            dexttheta = thitxtr[0]-thitxtr[thitxtr.size()-1];

            muextptreco = -1.746067*dextL*(1./6.)/dexttheta;
            
            /*
            int index = 0;
            for(int i = 1; i < nhitxtr; i++)
            {
              if(thitxtr[0] < thitxtr[i])
              {
                index++;
              }
              else
              {
                break;
              }
            }
            
            if(index != 0)
            { 
              dextL = sqrt((yhitxtr[0]-yhitxtr[index])*(yhitxtr[0]-yhitxtr[index])+
                        (zhitxtr[0]-zhitxtr[index])*(zhitxtr[0]-zhitxtr[index]));
  
              dexttheta = thitxtr[0]-thitxtr[index];
  
              muextptreco = -1.746067*dextL*(1./6.)/dexttheta;
            }
            */
        }
  
        if(nhititr>0&&nhitxtr>0)
        {
      		dintL = sqrt((yhititr-yhitxtr[yhitxtr.size()-1])*(yhititr-yhitxtr[yhitxtr.size()-1])+
                  (zhititr-zhitxtr[zhitxtr.size()-1])*(zhititr-zhitxtr[zhitxtr.size()-1]));

        	dinttheta = thitxtr[thitxtr.size()-1] - thititr;
      
      		muintptreco = -1.746067*dintL/dinttheta;
        }

        tout.Fill();
    }
    std::cout << std::endl;
    fout.cd();
    tout.Write();
    fout.Close();
    f.Close();
}
