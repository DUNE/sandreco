#include <TGeoManager.h>
#include <TString.h>
#include <TGeoNode.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TGeoTrd2.h>
#include <TGeoTube.h>
#include <TChain.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TSystem.h>

#include "TG4Event.h"
#include "TG4HitSegment.h"

#include <vector>
#include <map>
#include <iostream>
#include <iomanip>

#include "struct.h"
#include "utils.h"
#include "transf.h"

// Energy MeV
// Distance mm
// Time ns





double Attenuation(double d, int planeID)
{
    /*
       dE/dx attenuation - Ea=p1*exp(-d/atl1)+(1.-p1)*exp(-d/atl2)
       d    distance from photocatode - 2 cells/cell; d1 and d2
       atl1  50. cm
       atl2  430 cm planes 1-2    innermost plane is 1
       380 cm plane 3
       330 cm planes 4-5
       p1   0.35
       */
    const double p1 = 0.35;
    const double alt1 = 500.;
    double alt2 = 0.0;

    switch (planeID) {
        case 0:
        case 1:
            alt2 = 4300.0;
            break;

        case 2:
            alt2 = 3800.0;
            break;

        case 3:
        case 4:
            alt2 = 3300.0;
            break;

        default:
            // std::cout << "planeID out if range" << std::endl;
            alt2 = -999.0;
            break;
    }

    if (ns_Digit::debug) {
        std::cout << "planeID = " << planeID << std::endl;
        std::cout << "\tp1   = " << p1 << std::endl;
        std::cout << "\talt1 = " << alt1 << std::endl;
        std::cout << "\talt2 = " << alt2 << std::endl;
        std::cout << "\tatt  = "
            << p1* TMath::Exp(-d / alt1) + (1. - p1) * TMath::Exp(-d / alt2)
            << std::endl;
    }

    return p1 * TMath::Exp(-d / alt1) + (1. - p1) * TMath::Exp(-d / alt2);
}





double E2PE(double E)
{
    // Average number of photoelectrons = 25*Ea(MeV)
    const double e2p2 = 25.;

    if (ns_Digit::debug)
        std::cout << "E = " << E << " -> p.e. = " << e2p2* E << std::endl;

    return e2p2 * E;
}





double petime(double t0, double d)
{
    /*
       - For each photoelectron: Time for TDC simulation obtained from

       C  PHOTOELECTRON TIME :  Particle TIME in the cell
       C                      + SCINTILLATION DECAY TIME +
       C                      + signal propagation to the cell
       C                      + 1ns  uncertainty

       TPHE = Part_time+TSDEC+DPM1*VLFB+Gauss(1ns)

       VLFB = 5.85 ns/m
       !!!! Input-TDC Scintillation time -
       TSDEC = TSCIN*(1./RNDMPH(1)-1)**TSCEX  (ns)

       TSCIN  3.08  ns
       TSCEX  0.588
       */

    TRandom3 r(0);

    double mm_to_m = 1E-3;

    double tdec =
        ns_Digit::tscin * TMath::Power(1. / r.Uniform() - 1., ns_Digit::tscex);

    double time = t0 + tdec + ns_Digit::vlfb * d * mm_to_m + r.Gaus();

    if (ns_Digit::debug) {
        std::cout << "time : " << time << std::endl;
        std::cout << "t0   : " << t0 << std::endl;
        std::cout << "scint: " << tdec << std::endl;
        std::cout << "prop : " << ns_Digit::vlfb* d* mm_to_m << std::endl;
    }

    return time;
}





bool ProcessHitFluka(const TG4HitSegment& hit, int& modID,
        int& planeID, int& cellID, double& d1, double& d2, double& t,
        double& de)
{
    if (ns_Digit::debug) {
        std::cout << "ProcessHit FLUKA" << std::endl;
    }

    modID = -999;
    planeID = -999;
    cellID = -999;
    d1 = -999;
    d2 = -999;
    t = -999;

    double x = 0.5 * (hit.Start.X() + hit.Stop.X());
    double y = 0.5 * (hit.Start.Y() + hit.Stop.Y());
    double z = 0.5 * (hit.Start.Z() + hit.Stop.Z());

    t = 0.5 * (hit.Start.T() + hit.Stop.T());
    de = hit.EnergyDeposit;
    
    // Global to local coordinates
    //
    TLorentzVector globalPos(x, y, z, 0);
    TLorentzVector localPos = GlobalToLocalCoordinates(globalPos);
    x = localPos.X();
    y = localPos.Y();
    z = localPos.Z();
    double r = sqrt(y*y + z*z);
    if (ns_Digit::debug) std::cout << "x: " << x << "\ty: " << y << "\tz: " << z << "\tr: " << r;

    // hitAngle, cellAngle, modAngle
    // 
    double modDeltaAngle = 2.0 * TMath::Pi() / 24;             // 24 modules in a ring
    double cellDeltaAngle = 2.0 * TMath::Pi() / (24*12);       // 24 modules * 12 cells/module = number of cells in a ring
    double hitAngle;
    // This is the angle w.r.t. the y-axis. In Ideal2RealCal I used the angle w.r.t. the z-axis!
    if (z!=0) {
        if (z<0) hitAngle = 2 * atan( -z / ( y + sqrt( y*y + z*z)));                                  
        if (z>0) hitAngle = 2 * atan( -z / ( y + sqrt( y*y + z*z))) + 2 * TMath::Pi();
    }
    else if (z==0) {
        if (y<0) hitAngle = TMath::Pi();
        if (y>0) hitAngle = 0;
        if (y==0) hitAngle = 0;
    }
    double cellAngle = int(hitAngle / cellDeltaAngle) * cellDeltaAngle + cellDeltaAngle/2;
    double modAngle = int((hitAngle + 0.5 * modDeltaAngle) / modDeltaAngle) * modDeltaAngle;

    // Coordinates rotation and volume finding
    //
    TString str = "";
    double rotated_z = z * cos(-modAngle) - y * sin(-modAngle);
    double rotated_y = z * sin(-modAngle) + y * cos(-modAngle);
    // TODO: check units
    if ( (rotated_y > ns_Draw::kloe_int_R) && (rotated_y < ns_Draw::kloe_int_R + 2 * ns_Digit::ec_dzf) && (abs(x) < ns_Digit::lCalBarrel / 2 * 1000) )   str = "volECAL";        // ECAL barrel
    else if ( (rotated_y < ns_Digit::ec_rf) && (abs(x) > ns_Draw::kloe_int_dx) && (abs(x) < ns_Draw::kloe_int_dx + 2 * ns_Digit::ec_dzf) )      str = "endvolECAL";     // ECAL endcaps
    else if ( (rotated_y < ns_Digit::ec_rf) && (abs(x) < ns_Draw::kloe_int_dx) )                                                                str = "tracker";
    else                                                                                                                                        str = "outside";        // outside
    if (ns_Digit::debug) std::cout << "\tVol: " << str;

    // modID, planeID, cellID, d1, d2
    //
    double cellD = 0;
    if (str=="volECAL") {
        // modID
        modID = int((hitAngle + 0.5 * modDeltaAngle) / modDeltaAngle) % 24;
        // planeID
        planeID = int((rotated_y - ns_Draw::kloe_int_R) / 44);
        if (planeID > 4) planeID = 4;
        // cellID
        cellID = int((hitAngle + 0.5 * modDeltaAngle) / cellDeltaAngle) % 12;   // TODO: check ordering of cells (clockwise or anticlockwise?)
        // d1 distance from left end (x<0)
        d1 = 2150 + x;   // TODO: check units
        // d2 distance from right end (x>0)
        d2 = 2150 - x;   // TODO: check units
        // cellCoord
        cellD = ns_Draw::kloe_int_R + ns_Digit::dzlay[0] / 2;
        for (int planeindex=1; planeindex<planeID+1; planeindex++) cellD += ns_Digit::dzlay[planeindex-1] / 2 + ns_Digit::dzlay[planeindex] / 2;
        ns_Digit::cellCoordBarrel[modID][planeID][cellID][0] = 0;
        ns_Digit::cellCoordBarrel[modID][planeID][cellID][2] = + cellD * sin(-modAngle) - cellD * tan(cellAngle - modAngle) * cos(-modAngle); // TODO: fix tan argument
        ns_Digit::cellCoordBarrel[modID][planeID][cellID][1] = + cellD * cos(-modAngle) + cellD * tan(cellAngle - modAngle) * sin(-modAngle);
    } else if (str=="endvolECAL") {
        // modID
        if (x<0)        modID = 40;
        else if (x>0)   modID = 30;
        // planeID
        planeID = int((abs(x) - ns_Draw::kloe_int_dx) / 44);         // TODO: check units
        if (planeID > 4) planeID = 4;
        // cellID
        cellID = int((z + ns_Digit::ec_rf) / 44);                        
        // d1 distance from top (y>0)
        d1 =  sqrt(ns_Digit::ec_rf * ns_Digit::ec_rf - z * z) - y;
        // d2 distance from bottom (y<0)
        d2 =  sqrt(ns_Digit::ec_rf * ns_Digit::ec_rf - z * z) + y;
        // cellCoord
        cellD = TMath::Sign(1.0, x) * (ns_Draw::kloe_int_dx + ns_Digit::dzlay[0] / 2);
        for (int planeindex=1; planeindex<planeID+1; planeindex++) cellD += TMath::Sign(1.0, x) * (ns_Digit::dzlay[planeindex-1] / 2 + ns_Digit::dzlay[planeindex] / 2);
        ns_Digit::cellCoordEndcap[int(modID/10)][planeID][cellID][0] = cellD; 
        ns_Digit::cellCoordEndcap[int(modID/10)][planeID][cellID][1] = 0; 
        ns_Digit::cellCoordEndcap[int(modID/10)][planeID][cellID][2] = 44 / 2 + cellID * 44 - ns_Digit::ec_rf;
    } else if (str=="tracker" || str=="outside") {
        if (ns_Digit::debug) std::cout << std::endl;
        return false; 
    }
    if (ns_Digit::debug) std::cout << "\tmod " << modID << "\tplane: " << planeID << "\tcell: " << cellID << "\td1: " << d1 << "\td2: " << d2 << std::endl;

    return true;
}





bool ProcessHitEdepSim(TGeoManager* g, const TG4HitSegment& hit, int& modID,
        int& planeID, int& cellID, double& d1, double& d2, double& t,
        double& de)
{
    if (ns_Digit::debug) {
        std::cout << "ProcessHit EDEPSIM" << std::endl;
    }

    modID = -999;
    planeID = -999;
    cellID = -999;
    d1 = -999;
    d2 = -999;
    t = -999;

    double x = 0.5 * (hit.Start.X() + hit.Stop.X());
    double y = 0.5 * (hit.Start.Y() + hit.Stop.Y());
    double z = 0.5 * (hit.Start.Z() + hit.Stop.Z());

    t = 0.5 * (hit.Start.T() + hit.Stop.T());
    de = hit.EnergyDeposit;

    TGeoNode* node = g->FindNode(x, y, z);
    if (node == 0) return false;
    TString str = node->GetName();
    TString str2 = g->GetPath();
    TObjArray* obj = str2.Tokenize("/");

    int size = obj->GetEntries();
    if (size < 6) {
        return false;
    };

    str2 = ((TObjString*)obj->At(5))->GetString();
    delete obj;


    if (ns_Digit::debug) {
        std::cout << "node name: " << str.Data() << std::endl;
    }
    //if node is in ecal in active and in barrel
    if (str.Contains("volECAL") == true && str.Contains("Active") == true && str.Contains("end") == false) { 
        TObjArray* obja = str.Tokenize("_");
        TObjArray* obja2 = str2.Tokenize("_");

        int slabID;
        //top module => modID == 0
        //increasing modID counterclockwise as seen from positive x 
        //(i.e. z(modID==1) < z(modID==0) & z(modID==0) < z(modID==23)) 
        modID = ((TObjString*)obja2->At(3))->GetString().Atoi();
        slabID = ((TObjString*)obja->At(1))->GetString().Atoi();

        delete obja;
        delete obja2;

        // planeID==0 -> smallest slab -> internal
        // planeID==208 -> biggest slab -> external 
        planeID = slabID / 40;

        if (planeID > 4) planeID = 4;

        double Pmaster[3];
        double Plocal[3];
        Pmaster[0] = x;
        Pmaster[1] = y;
        Pmaster[2] = z;

        g->GetCurrentNavigator()->MasterToLocal(Pmaster, Plocal);

        TGeoTrd2* trd = (TGeoTrd2*)node->GetVolume()->GetShape();

        if (ns_Digit::debug) {
            std::cout << "pointer: " << trd << std::endl;
        }

        double dx1 = trd->GetDx1();
        double dx2 = trd->GetDx2();
        double dz = trd->GetDz();
        double dy1 = trd->GetDy1();
        double dy2 = trd->GetDy2();

        // d1 distanza da estremo left (x<0)
        // d2 distanza da estremo right (x>0)
        d1 = dy1 + Plocal[1];
        d2 = dy1 - Plocal[1];

        // http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/html/Detector/Geometry/geomSolids.html
        // if z = -dz -> dx = 2*dx1
        // if z =  dz -> dx = 2*dx2
        // semilarghezza della slab di scintillatore alla quota Plocal[2]
        double dx = 0.5 * Plocal[2]/dz * (dx2 - dx1) + 0.5 * (dx2 + dx1);

        // Cell width at z = Plocal[2]
        double cellw = 2. * dx / 12.;

        // cellID = distanza dall'estremo diviso larghezza cella
        cellID = (Plocal[0] + dx) / cellw;

        if (ns_Digit::debug) {
            std::cout << "hit: " << str.Data() << std::endl;
            std::cout << "\t[x,y,z]                " << x << " " << y << " " << z
                << std::endl;
            std::cout << "\t[modID,planeID,cellID] " << modID << " " << planeID << " "
                << cellID << std::endl;
            std::cout << "\t[d1,d2,t,de]           " << d1 << " " << d2 << " " << t
                << " " << de << std::endl;
        }

        return true;
    }
    // end cap modules
    else if (str.Contains("endvolECAL") == true &&
            str.Contains("Active") == true) {
        TObjArray* obja = str.Tokenize("_");
        TObjArray* obja2 = str2.Tokenize("_");

        int slabID;
        modID = ((TObjString*)obja2->At(4))->GetString().Atoi();
        slabID = ((TObjString*)obja->At(1))->GetString().Atoi();

        // mod == 40 -> left
        // mod == 30 -> right
        if (modID == 0)
            modID = 40;
        else if (modID == 1)
            modID = 30;

        delete obja;
        delete obja2;

        // planeID==0 -> internal
        // planeID==208 -> external
        planeID = slabID / 40;

        if (planeID > 4) planeID = 4;

        double Pmaster[3];
        double Plocal[3];
        Pmaster[0] = x;
        Pmaster[1] = y;
        Pmaster[2] = z;

        g->GetCurrentNavigator()->MasterToLocal(Pmaster, Plocal);

        TGeoTube* tub = (TGeoTube*)node->GetVolume()->GetShape();

        if (ns_Digit::debug) {
            std::cout << "pointer: " << tub << std::endl;
        }

        double rmin = tub->GetRmin();
        double rmax = tub->GetRmax();
        double dz = tub->GetDz();

        // d1 distanza da estremo up (y>0)
        // d2 distanza da estremo down (y<0)
        d1 = rmax * TMath::Sin(TMath::ACos(Plocal[0] / rmax)) - Plocal[1];
        d2 = rmax * TMath::Sin(TMath::ACos(Plocal[0] / rmax)) + Plocal[1];

        cellID = int((Plocal[0] / rmax + 1.) * 45);

        if (ns_Digit::debug) {
            std::cout << "hit: " << str.Data() << std::endl;
            std::cout << "\t[x,y,z]                " << x << " " << y << " " << z
                << std::endl;
            std::cout << "\t[modID,planeID,cellID] " << modID << " " << planeID << " "
                << cellID << std::endl;
        }

        return true;
    } else {
        return false;
    }
}





void SimulatePEFluka(TG4Event* ev, 
        std::map<int, std::vector<double> >& time_pe,
        std::map<int, std::vector<int> >& id_hit,
        std::map<int, double>& L)
{
    int modID, planeID, cellID, id;
    double d1, d2, t0, de;

    TRandom3 r(0);

    //The geometry is needed for distinghish between endcap and barrel

    for (std::map<std::string, std::vector<TG4HitSegment> >::iterator it =
            ev->SegmentDetectors.begin();
            it != ev->SegmentDetectors.end(); ++it) {
        if (it->first == "EMCalSci") {
            for (unsigned int j = 0; j < it->second.size(); j++) {
                if (ProcessHitFluka(it->second[j], modID, planeID, cellID, d1, d2, t0, de) == true )  {
                    double en1 = de * Attenuation(d1, planeID);
                    double en2 = de * Attenuation(d2, planeID);

                    double ave_pe1 = E2PE(en1);
                    double ave_pe2 = E2PE(en2);

                    int pe1 = r.Poisson(ave_pe1);
                    int pe2 = r.Poisson(ave_pe2);

                    id = cellID + 100 * planeID + 1000 * modID;

                    if (ns_Digit::debug) {
                        std::cout << "cell ID: " << id << std::endl;
                        std::cout << "\t" << de << " " << en1 << " " << en2 << std::endl;
                        std::cout << "\t" << ave_pe1 << " " << ave_pe2 << std::endl;
                        std::cout << "\t" << pe1 << " " << pe2 << std::endl;
                    }

                    // cellend 1 -> x < 0 -> ID > 0 -> left
                    // cellend 2 -> x > 0 -> ID < 0 -> right

                    for (int i = 0; i < pe1; i++) {
                        time_pe[id].push_back(petime(t0, d1));
                        id_hit[id].push_back(j);
                        L[id] = d1 + d2;
                    }

                    for (int i = 0; i < pe2; i++) {
                        time_pe[-1 * id].push_back(petime(t0, d2));
                        id_hit[-1 * id].push_back(j);
                        L[-1 * id] = d1 + d2;
                    }
                }
            }
        }
    }
}





void SimulatePEEdepSim(TG4Event* ev, TGeoManager* g,
        std::map<int, std::vector<double> >& time_pe,
        std::map<int, std::vector<int> >& id_hit,
        std::map<int, double>& L)
{
    int modID, planeID, cellID, id;
    double d1, d2, t0, de;

    TRandom3 r(0);

    //The geometry is needed for distinghish between endcap and barrel

    for (std::map<std::string, std::vector<TG4HitSegment> >::iterator it =
            ev->SegmentDetectors.begin();
            it != ev->SegmentDetectors.end(); ++it) {
        if (it->first == "EMCalSci") {
            for (unsigned int j = 0; j < it->second.size(); j++) {
                if (ProcessHitEdepSim(g, it->second[j], modID, planeID, cellID, d1, d2, t0,de) == true )  {
                    double en1 = de * Attenuation(d1, planeID);
                    double en2 = de * Attenuation(d2, planeID);

                    double ave_pe1 = E2PE(en1);
                    double ave_pe2 = E2PE(en2);

                    int pe1 = r.Poisson(ave_pe1);
                    int pe2 = r.Poisson(ave_pe2);

                    id = cellID + 100 * planeID + 1000 * modID;

                    if (ns_Digit::debug) {
                        std::cout << "cell ID: " << id << std::endl;
                        std::cout << "\t" << de << " " << en1 << " " << en2 << std::endl;
                        std::cout << "\t" << ave_pe1 << " " << ave_pe2 << std::endl;
                        std::cout << "\t" << pe1 << " " << pe2 << std::endl;
                    }

                    // cellend 1 -> x < 0 -> ID > 0 -> left
                    // cellend 2 -> x > 0 -> ID < 0 -> right

                    for (int i = 0; i < pe1; i++) {
                        time_pe[id].push_back(petime(t0, d1));
                        id_hit[id].push_back(j);
                        L[id] = d1 + d2;
                    }

                    for (int i = 0; i < pe2; i++) {
                        time_pe[-1 * id].push_back(petime(t0, d2));
                        id_hit[-1 * id].push_back(j);
                        L[-1 * id] = d1 + d2;
                    }
                }
            }
        }
    }
}





void TimeAndSignal(std::map<int, std::vector<double> >& time_pe,
        std::map<int, double>& adc, std::map<int, double>& tdc)
{
    /*
       -  ADC - Proportional to NPHE
       -  TDC - Constant fraction - simulated
       TPHE(1...NPHE) in increasing time order
       IND_SEL= 0.15*NPHE
       TDC_cell = TPHE(IND_SEL)
       */

    const double pe2ADC = 1.0;

    for (std::map<int, std::vector<double> >::iterator it = time_pe.begin();
            it != time_pe.end(); ++it) {
        adc[it->first] = pe2ADC * it->second.size();
        std::sort(it->second.begin(), it->second.end());
        int index = 0.15 * it->second.size();
        tdc[it->first] = it->second[index];
    }
}





void CollectSignal(TGeoManager* geo,
        std::map<int, std::vector<double> >& time_pe,
        std::map<int, double>& adc, std::map<int, double>& tdc,
        std::map<int, double>& L,
        std::map<int, std::vector<int> >& id_hit,
        std::vector<cell>& vec_cell)
{
    for (std::map<int, std::vector<double> >::iterator it = time_pe.begin();
            it != time_pe.end(); ++it) {
        if (it->first < 0) continue;

        cell c;                       
        c.id = it->first;
        c.adc1 = adc[it->first];
        c.tdc1 = tdc[it->first];
        c.adc2 = adc[-1 * it->first];
        c.tdc2 = tdc[-1 * it->first];
        c.pe_time1 = time_pe[it->first];
        c.pe_time2 = time_pe[-1 * it->first];
        c.hindex1 = id_hit[it->first];
        c.hindex2 = id_hit[-1 * it->first];
        c.l = L[it->first];
        c.x = 0;
        c.y = 0;
        c.z = 0;

        c.mod = c.id / 1000;
        c.lay = (c.id - c.mod * 1000) / 100;
        c.cel = c.id - c.mod * 1000 - c.lay * 100;

        double dummyLoc[3];
        double dummyMas[3];

        if(ns_Digit::flukatype==false){

            if (c.mod < 24) {
                dummyLoc[0] = ns_Digit::cxlay[c.lay][c.cel];  //coordinate nel local del modulo in base al numero del modulo e al numero della cella  
                dummyLoc[1] = 0.;
                dummyLoc[2] = ns_Digit::czlay[c.lay];

                geo->cd(TString::Format(ns_Digit::path_barrel_template, c.mod).Data());
                geo->LocalToMaster(dummyLoc, dummyMas);

                c.x = dummyMas[0];   //coordinate della cella nel Global
                c.y = dummyMas[1];
                c.z = dummyMas[2];
            } else if (c.mod == 30 || c.mod == 40)
                // right x > 0 : c.mod = 30
                // left  x < 0 : c.mod = 40
            {
                dummyLoc[0] = ns_Digit::ec_r / 45. * (0.5 + c.cel) - ns_Digit::ec_r;
                dummyLoc[1] = 0.;
                dummyLoc[2] = ns_Digit::czlay[c.lay];

                if (c.mod == 30) {
                    geo->cd(ns_Digit::path_endcapR_template);
                } else if (c.mod == 40) {
                    geo->cd(ns_Digit::path_endcapL_template);
                }

                geo->LocalToMaster(dummyLoc, dummyMas);

                c.x = dummyMas[0];
                c.y = dummyMas[1];
                c.z = dummyMas[2];
            }

        } else if (ns_Digit::flukatype==true) {

            if (c.mod < 24)                                                 // Barrel
            {

                // Local coordinates calculation
                dummyLoc[0] = ns_Digit::cellCoordBarrel[c.mod][c.lay][c.cel][0];  
                dummyLoc[1] = ns_Digit::cellCoordBarrel[c.mod][c.lay][c.cel][1];
                dummyLoc[2] = ns_Digit::cellCoordBarrel[c.mod][c.lay][c.cel][2];

                // Transformation to global coordinates
                dummyMas[0] = LocalToGlobalCoordinates(dummyLoc).X();
                dummyMas[1] = LocalToGlobalCoordinates(dummyLoc).Y();
                dummyMas[2] = LocalToGlobalCoordinates(dummyLoc).Z();

                // Global coordinates are assigned to the cell position
                c.x = dummyMas[0];
                c.y = dummyMas[1];
                c.z = dummyMas[2];

            } else if (c.mod == 30 || c.mod == 40)                          // Endcaps
            {   

                // right x > 0 : c.mod = 30
                // left  x < 0 : c.mod = 40

                // Local coordinates calculation
                dummyLoc[0] = ns_Digit::cellCoordEndcap[int(c.mod/10)][c.lay][c.cel][0];
                dummyLoc[1] = ns_Digit::cellCoordEndcap[int(c.mod/10)][c.lay][c.cel][1];
                dummyLoc[2] = ns_Digit::cellCoordEndcap[int(c.mod/10)][c.lay][c.cel][2];

                // Transformation to global coordinates
                dummyMas[0] = LocalToGlobalCoordinates(dummyLoc).X();
                dummyMas[1] = LocalToGlobalCoordinates(dummyLoc).Y();
                dummyMas[2] = LocalToGlobalCoordinates(dummyLoc).Z();

                // Global coordinates are assigned to the cell position
                c.x = dummyMas[0];
                c.y = dummyMas[1];
                c.z = dummyMas[2];
            }

        }
        vec_cell.push_back(c);
    }
}





void init(TGeoManager* geo)  //initialize the calorimeter dimensions
{

    TGeoTrd2* mod = (TGeoTrd2*)geo->FindVolumeFast("ECAL_lv_PV")->GetShape(); //read from geometry the useful dimensions

    // https://root.cern.ch/root/htmldoc/guides/users-guide/Geometry.html#shapes
    // GetDx1() half length in x at -Dz  =262.55
    // GetDx2() half length in x at +Dz  =292.85
    // Dx1 < Dx2 => -Dz corresponds to minor width => internal side 
    //dz half of module height:           =115

    double xmin; 
    double xmax; 
    double dz; 

    if(ns_Digit::flukatype==false){
        xmin = mod->GetDx1();
        xmax = mod->GetDx2();
        dz = mod->GetDz();

        if((abs(xmin-ns_Digit::xmin_f)>0.2) || (abs(xmax-ns_Digit::xmax_f)>0.2) || (abs(dz-ns_Digit::dz_f)>0.2)) {
            std::cout<<"ERROR ON ECAL GEOMETRY: xmin= "<<xmin<<" instead of "<<ns_Digit::xmin_f<<std::endl; 
            std::cout<<"ERROR ON ECAL GEOMETRY: xmax= "<<xmax<<" instead of "<<ns_Digit::xmax_f<<std::endl; 
            std::cout<<"ERROR ON ECAL GEOMETRY: dz= "<<dz<<" instead of "<<ns_Digit::dz_f<<std::endl; 
            exit(1);
        }

    } else{
        xmin = ns_Digit::xmin_f;
        xmax = ns_Digit::xmax_f;
        dz = ns_Digit::dz_f;
    }    


    double m = 0.5*(xmax - xmin)/dz;
    double q = 0.5*(xmax + xmin);

    // z edge of the cells
    double zlevel[ns_Digit::nLay+1];
    zlevel[0] = -dz;

    for(int i = 0; i < ns_Digit::nLay; i++)
    {
        zlevel[i+1] = zlevel[i] + ns_Digit::dzlay[i];
    }

    // z position of the center of the cells
    for (int i = 0; i < ns_Digit::nLay; i++) {
        ns_Digit::czlay[i] = 0.5* (zlevel[i]+zlevel[i+1]);

        // total module width at the z position of the center of the cell 
        double xwidth = 2*(m*ns_Digit::czlay[i]+q);

        // cell width at the z position of the center of the cell 
        double dx = xwidth/ns_Digit::nCel;

        // x position of the center of the cells
        for (int j = 0; j < ns_Digit::nCel; j++) {
            ns_Digit::cxlay[i][j] = dx * (j+0.5) - xwidth*0.5;
            if(j==11) std::cout<<"lay "<<i<<"cx czlay "<<ns_Digit::cxlay[i][j]<<" "<<ns_Digit::czlay[i]<<std::endl;
        }
    }


    //dimensions of endcaps     
    if(ns_Digit::flukatype==false){
        TGeoTube* ec = (TGeoTube*)geo->FindVolumeFast("ECAL_end_lv_PV")->GetShape();

        ns_Digit::ec_r = ec->GetRmax();   //Maximum radius    =2000
        ns_Digit::ec_dz = ec->GetDz();    //half of thickness =115
        if(abs(ns_Digit::ec_r-ns_Digit::ec_rf)>0.2 || (abs(ns_Digit::ec_dz-ns_Digit::ec_dzf))) {
            std::cout<<"ERROR ON ECAL ENDCAP GEOMETRY: R= "<<ns_Digit::ec_r<<" instead of "<<ns_Digit::ec_rf<<std::endl; 
            std::cout<<"ERROR ON ECAL ENDCAP GEOMETRY: Thickness= "<<ns_Digit::ec_dz<<" instead of "<<ns_Digit::ec_dzf<<std::endl; 
            exit(1);
        }
    } else{
        ns_Digit::ec_r=ns_Digit::ec_rf;
        ns_Digit::ec_dz=ns_Digit::ec_dzf;
    }
}





void DigitizeCal(TG4Event* ev, TGeoManager* geo, std::vector<cell>& vec_cell)
{
    std::map<int, std::vector<double> > time_pe;
    std::map<int, std::vector<int> > id_hit;
    std::map<int, double> adc;
    std::map<int, double> tdc;
    std::map<int, double> L;

    vec_cell.clear();

    if (ns_Digit::debug) {
        std::cout << "SimulatePE" << std::endl;
    }
    if(ns_Digit::flukatype==false) SimulatePEEdepSim(ev, geo, time_pe, id_hit, L);
    if(ns_Digit::flukatype==true) SimulatePEFluka(ev, time_pe, id_hit, L);

    if (ns_Digit::debug) {
        std::cout << "TimeAndSignal" << std::endl;
    }
    TimeAndSignal(time_pe, adc, tdc);
    if (ns_Digit::debug) {
        std::cout << "CollectSignal" << std::endl;
    }
    CollectSignal(geo, time_pe, adc, tdc, L, id_hit, vec_cell); //TODO: provare a lasciarla unica??
}





void Cluster(TG4Event* ev, TGeoManager* geo,
        std::map<std::string, std::vector<hit> >& cluster_map)
{
    cluster_map.clear();

    for (unsigned int j = 0; j < ev->SegmentDetectors["Straw"].size(); j++) {
        const TG4HitSegment& hseg = ev->SegmentDetectors["Straw"].at(j);

        double x = 0.5 * (hseg.Start.X() + hseg.Stop.X());
        double y = 0.5 * (hseg.Start.Y() + hseg.Stop.Y());
        double z = 0.5 * (hseg.Start.Z() + hseg.Stop.Z());

        std::string sttname = geo->FindNode(x, y, z)->GetName();

        hit h;
        h.det = sttname; 
        h.x1 = hseg.Start.X();
        h.y1 = hseg.Start.Y();
        h.z1 = hseg.Start.Z();
        h.t1 = hseg.Start.T();
        h.x2 = hseg.Stop.X();
        h.y2 = hseg.Stop.Y();
        h.z2 = hseg.Stop.Z();
        h.t2 = hseg.Stop.T();
        h.de = hseg.EnergyDeposit;
        h.pid = hseg.PrimaryId;
        h.index = j;

        std::string cluster_name(sttname);
        cluster_name += "_" + std::to_string(hseg.PrimaryId);

        cluster_map[cluster_name].push_back(h);
    }
}





void Cluster2Digit(std::map<std::string, std::vector<hit> >& cluster_map,
        std::vector<digit>& digit_vec)
{
    for (std::map<std::string, std::vector<hit> >::iterator it =
            cluster_map.begin();
            it != cluster_map.end(); ++it) {
        digit d;
        d.de = 0;
        d.det = it->second[0].det;

        for (unsigned int k = 0; k < it->second.size(); k++) {
            d.hindex.push_back(it->second[k].index);
            d.de += it->second[k].de;
        }

        d.hor = (d.det.find("hor") != std::string::npos) ? false : true;

        std::sort(it->second.begin(), it->second.end(), isHitBefore);
        d.t = it->second.at(0).t1;
        d.x = 0.5 * (it->second.front().x1 + it->second.back().x2);
        d.y = 0.5 * (it->second.front().y1 + it->second.back().y2);
        d.z = 0.5 * (it->second.front().z1 + it->second.back().z2);

        digit_vec.push_back(d);
    }
}





void DigitizeStt(TG4Event* ev, TGeoManager* geo, std::vector<digit>& digit_vec)
{
    std::map<std::string, std::vector<hit> > cluster_map;
    digit_vec.clear();

    Cluster(ev, geo, cluster_map);
    Cluster2Digit(cluster_map, digit_vec);
}





void Digitize(const char* finname, const char* foutname)
{
    //TChain* t = new TChain("EDepSimEvents","EDepSimEvents");
    //t->Add(finname);
    //TFile f(t->GetListOfFiles()->At(0)->GetTitle());
    TFile f(finname,"READ");

    //TString (finname.c_str();
    if(TString(finname).Contains("FLUKA") == true) {
        ns_Digit::flukatype=true; //dobbiamo leggere GeneratorName ..cambiare quando fatto
        std::cout<<"This is a FLUKA SIMULATION"<<std::endl; 
    }
    TTree* t = (TTree*) f.Get("EDepSimEvents");

    TGeoManager* geo = 0; 
    TTree* gRooTracker = 0; 
    TTree* InputKinem = 0; 
    TTree* InputFiles = 0; 	

    if(ns_Digit::flukatype==false){
        //initialize geometry for edepsim file
        geo = (TGeoManager*) f.Get("EDepSimGeometry"); 
        gRooTracker = (TTree*) f.Get("DetSimPassThru/gRooTracker");   //FIXME dobbiamo metterlo anche nei file di FLUKA...togliere da questo if quando fatto
        InputKinem = (TTree*) f.Get("DetSimPassThru/InputKinem");
        InputFiles = (TTree*) f.Get("DetSimPassThru/InputFiles");

        init(geo);
    } 

    TG4Event* ev = new TG4Event;
    t->SetBranchAddress("Event",&ev);

    std::vector<digit> digit_vec;    
    std::vector<cell> vec_cell;

    TFile fout(foutname,"RECREATE");
    TTree tout("tDigit","Digitization");
    tout.Branch("cell","std::vector<cell>",&vec_cell);
    tout.Branch("Stt","std::vector<digit>",&digit_vec);

    const int nev = t->GetEntries();

    std::cout << "Events: " << nev << " [";
    std::cout << std::setw(3) << int(0) << "%]" << std::flush;

    for(int i = 0; i < nev; i++)
    {
        t->GetEntry(i);

        std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;

        // FOR CHANGING COORDIN FIXME
        for (auto i = vec_cell.begin(); i != vec_cell.end(); ++i){

            cell cel = *i;
            TLorentzVector prova(cel.x,cel.y,cel.z,0);
            TLorentzVector v = GlobalToLocalCoordinates(prova);
            //(*i).x=v.X()/10;
            //(*i).y=v.Y()/10;
            //(*i).z=v.Z()/10;

            float R=sqrt(v.Y()*v.Y() + v.Z()*v.Z());
            if(R>2230) {std::cout << "ERROR " << std::endl; 
                std::cout << "vec " << cel.x << " " << cel.y << " " << cel.z << std::endl;
                std::cout << "local coord " << v.X() << " " << v.Y() << " " << v.Z() <<std::endl;
                std::cout << "mod lay cell " << cel.mod << " " << cel.lay << " " << cel.cel << std::endl;

            }
        }
        /////////////////////////


        DigitizeCal(ev, geo, vec_cell);
        if(ns_Digit::flukatype==false) DigitizeStt(ev, geo, digit_vec);
        //}else{
        //	std::cout<<"Non ancora implementato in fluka"<<std::endl;
        //	return;
        //}

        //std::cout<<"Prima del salvataggio "<<std::endl;
        //USEFULE FUNCTION FOR CHANGING COORDINATE TO SAND CENTER
        /*
           for (auto i = vec_cell.begin(); i != vec_cell.end(); ++i){
           cell cel=*i;
           TLorentzVector prova(cel.x,cel.y,cel.z,0);
           TLorentzVector v=GlobalToLocalCoordinates(prova);
        //(*i).x=v.X()/10;
        //(*i).y=v.Y()/10;
        //(*i).z=v.Z()/10;

        float R=sqrt(v.Y()*v.Y()+v.Z()*v.Z());
        }
        */
        tout.Fill();
}

std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
std::cout << std::endl;

fout.cd();
tout.Write();
if(ns_Digit::flukatype==false)   geo->Write();
t->CloneTree()->Write();
if(ns_Digit::flukatype==false)  gRooTracker->CloneTree()->Write();
if(ns_Digit::flukatype==false)  InputKinem->CloneTree()->Write();
if(ns_Digit::flukatype==false)  InputFiles->CloneTree()->Write();
fout.Close();

f.Close();
}





void help_digit()
{
    std::cout << "Digitize <input file> <output file>" << std::endl;
    std::cout << "input file name could contain wild card" << std::endl;
}





int main(int argc, char* argv[])
{
    if (argc != 3)
        help_digit();
    else
        Digitize(argv[1], argv[2]);
}
