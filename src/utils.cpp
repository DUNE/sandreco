#include "struct.h"
#include "utils.h"

#include "TG4Event.h"
#include "TG4HitSegment.h"

#include <TChain.h>
#include <TGeoManager.h>
#include <TCanvas.h>
#include <TGeoTrd2.h>
#include <TGeoTube.h>
#include <TRandom3.h>

#include <iostream>
#include "transf.h"

namespace kloe_simu
{
  double cellCoordBarrel[nMod][nLay][nCel][3];
  double cellCoordEndcap[5][nLay][90][3];

  double czlay[nLay];
  double cxlay[nLay][nCel];
  double ec_r;
  double ec_dz;

  double stt_center[3];

  TPRegexp* rST;
  TPRegexp* rSTplane;

  std::map<int, std::map<double, int> > stX;
  std::map<int, double> stL;
  std::map<int, std::map<int, TVector2> > stPos;
  std::map<int, TVector2> tubePos;
  std::map<int, double> t0;
} 

bool kloe_simu::isPeBefore(const pe& p1, const pe& p2)
{
  return p1.time < p2.time;
} 

// value of parameter of segment (y1,z1,y2,z2)
// corresponding to minimal distance to point (y,z)
double kloe_simu::getT(double y1, double y2, double y, double z1, double z2,
                       double z)
{
  double t = 0;
  if (y1 != y2 || z1 != z2) {
    t = -((y1 - y) * (y2 - y1) + (z1 - z) * (z2 - z1)) /
        ((y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
  }
  if (t < 0)
    return 0;
  else if (t > 1)
    return 1;
  else
    return t;
}

int kloe_simu::encodeSTID(int planeid, int tubeid)
{
  return tubeid * 1000 + planeid;
}

void kloe_simu::decodeSTID(int id, int& planeid, int& tubeid)
{
  tubeid = id / 1000;
  planeid = id % 1000;
}

int kloe_simu::encodePlaneID(int moduleid, int type)
{
  return moduleid * 10 + type;
}

void kloe_simu::decodePlaneID(int id, int& moduleid, int& type)
{
  moduleid = id / 10;
  type = id % 10;
}

bool kloe_simu::isST(TString name)
{
  return name.Contains(*kloe_simu::rST);
}

bool kloe_simu::isSTPlane(TString name)
{
  return name.Contains(*kloe_simu::rSTplane);
}

// get local id of tube
int kloe_simu::getSTId(TString name)
{
  int id = -999;

  TObjArray* obj = name.Tokenize("_");

  if (obj->GetEntries() != 10 && obj->GetEntries() != 9) {
    std::cout << "Error: tokenizing " << name.Data() << std::endl;
  } else {
    TString sid = ((TObjString*)obj->At(obj->GetEntries() - 1))->GetString();

    id = sid.Atoi();
  }
  delete obj;

  return id;
}

// get plane id
int kloe_simu::getPlaneID(TString name)
{
  int mod = 0;
  int type = 0;

  TObjArray* obj = name.Tokenize("_");

  if (obj->GetEntries() != 7 && obj->GetEntries() != 6 &&
      obj->GetEntries() != 11 && obj->GetEntries() != 10) {
    std::cout << "Error: tokenizing " << name.Data() << std::endl;
  } else {
    int ipos = (obj->GetEntries() == 10 || obj->GetEntries() == 6) ? 2 : 3;
    TString stype = ((TObjString*)obj->At(ipos))->GetString();
    TString smod = ((TObjString*)obj->At(1))->GetString();

    mod = smod.Atoi();

    if (stype.Contains("hor2"))
      type = 4;
    else if (stype.Contains("ver"))
      type = 1;
    else if (stype.Contains("hor"))
      type = 2;
    else
      std::cout << "Error evaluating type for: " << name.Data() << std::endl;
  }

  delete obj;

  return encodePlaneID(mod, type);
}

// get position of the center of the tube for a plane
void kloe_simu::getSTinfo(TGeoNode* nod, TGeoHMatrix mat, int pid,
                          std::map<double, int>& stX,
                          std::map<int, double>& stL,
                          std::map<int, TVector2>& stPos)
{
  int type;
  int mod;
  decodePlaneID(pid, mod, type);
  int ic = 1 - (type % 2);

  if (ic != 0 && ic != 1)
    std::cout << "Error: ic expected 0 or 1 -> " << ic << std::endl;

  for (int i = 0; i < nod->GetNdaughters(); i++) {
    TGeoNode* dau = nod->GetDaughter(i);
    TGeoTube* tub = (TGeoTube*)dau->GetVolume()->GetShape();
    double lenght = 2 * tub->GetDz();
    TString name = dau->GetName();

    if (!isST(name))
      std::cout << "Error: expected ST but not -> " << name.Data() << std::endl;

    TGeoMatrix* thismat = nod->GetDaughter(i)->GetMatrix();
    TGeoHMatrix mymat = mat * (*thismat);

    int id = getSTId(name);

    TVector2 v;
    v.SetX(mymat.GetTranslation()[2]);
    v.SetY(mymat.GetTranslation()[ic]);

    stX[v.Y()] = id;
    stL[id] = lenght;
    stPos[id] = v;
  }
}

// get position of the center of the tube for each plane
void kloe_simu::getSTPlaneinfo(TGeoNode* nod, TGeoHMatrix mat,
                               std::map<int, std::map<double, int> >& stX,
                               std::map<int, double>& stL,
                               std::map<int, std::map<int, TVector2> >& stPos)
{
  TString name = nod->GetName();
  TGeoMatrix* thismat = nod->GetMatrix();
  TGeoHMatrix mymat = mat * (*thismat);

  int pid = 0;
  double x = 0;
  double z = 0;

  if (isSTPlane(name)) {
    pid = getPlaneID(name);

    std::map<double, int> mstX;
    std::map<int, TVector2> mstPos;

    getSTinfo(nod, mymat, pid, mstX, stL, mstPos);

    stX[pid] = mstX;
    stPos[pid] = mstPos;
  } else {

    for (int i = 0; i < nod->GetNdaughters(); i++) {
      getSTPlaneinfo(nod->GetDaughter(i), mymat, stX, stL, stPos);
    }
  }
}

// get unique id of the tube corresponding to (x,y,z)
int kloe_simu::getSTUniqID(TGeoManager* g, double x, double y, double z)
{
  TString sttname = g->FindNode(x, y, z)->GetName();

  int sid = -999;
  int pid = getPlaneID(sttname);

  if (pid == 0) return -999;

  double xt = 0.;

  if (pid % 2 == 1)
    xt = x;
  else
    xt = y;

  std::map<double, int>::iterator it = kloe_simu::stX.at(pid).lower_bound(xt);

  if (it == kloe_simu::stX.at(pid).begin()) {
    sid = kloe_simu::stX.at(pid).begin()->second;
  } else if (it == kloe_simu::stX.at(pid).end()) {
    sid = kloe_simu::stX.at(pid).rbegin()->second;
  } else {
    TVector2 v1 = kloe_simu::stPos.at(pid).at(it->second);
    TVector2 v2 = kloe_simu::stPos.at(pid).at(std::prev(it)->second);

    TVector2 v(z, xt);

    if ((v - v1).Mod() > (v - v2).Mod()) {
      if ((v - v2).Mod() > 5)
        std::cout << "Error: distance grater than ST radius" << std::endl;

      sid = std::prev(it)->second;
    } else {
      if ((v - v1).Mod() > 5)
        std::cout << "Error: distance grater than ST radius" << std::endl;

      sid = it->second;
    }
  }

  return encodeSTID(pid, sid);
}

bool kloe_simu::isCluBigger(const std::vector<dg_tube>& v1,
                            const std::vector<dg_tube>& v2)
{
  return v1.size() > v2.size();
}

bool kloe_simu::isDigUpstream(const dg_tube& d1, const dg_tube& d2)
{
  return d1.z < d2.z;
}

bool kloe_simu::isHitBefore(hit h1, hit h2)
{
  return h1.t1 < h2.t1;
}

bool kloe_simu::isDigBefore(dg_tube d1, dg_tube d2)
{
  return d1.tdc < d2.tdc;
}

bool kloe_simu::isCellBefore(dg_cell c1, dg_cell c2)
{
  if (c1.ps1.adc.size() == 0 || c1.ps2.adc.size() == 0)
    return false;
  else if (c2.ps1.adc.size() == 0 || c2.ps2.adc.size() == 0)
    return true;
  else
    return ((c1.ps1.tdc.at(0) + c1.ps2.tdc.at(0)) <
            (c2.ps1.tdc.at(0) + c2.ps2.tdc.at(0)));
}

bool kloe_simu::isAfter(particle p1, particle p2)
{
  return p1.tid > p2.tid;
}

bool kloe_simu::isBarrel(TString& str)
{
  // something like: volECALActiveSlab_21_PV_0
  return str.Contains("volECAL") == true && str.Contains("Active") == true &&
         str.Contains("end") == false;
}

bool kloe_simu::isEndCap(TString& str)
{
  // something like: endvolECALActiveSlab_0_PV_0
  return str.Contains("endvolECAL") == true && str.Contains("Active") == true;
}

// find module id and plane id for barrel
void kloe_simu::BarrelModuleAndLayer(TString& str, TString& str2, int& modID,
                                     int& planeID)
{
  TObjArray* obja = str.Tokenize("_");    // BARERL => volECALActiveSlab_21_PV_0
  TObjArray* obja2 = str2.Tokenize("_");  // BARREL => ECAL_lv_PV_18

  int slabID;
  // top module => modID == 0
  // increasing modID counterclockwise as seen from positive x
  //(i.e. z(modID==1) < z(modID==0) & z(modID==0) < z(modID==23))
  modID = ((TObjString*)obja2->At(3))->GetString().Atoi();
  slabID = ((TObjString*)obja->At(1))->GetString().Atoi();

  delete obja;
  delete obja2;

  // planeID==0 -> smallest slab -> internal
  // planeID==208 -> biggest slab -> external
  planeID = slabID / 40;

  if (planeID > 4) planeID = 4;
}

// find module id and plane id for endcaps
void kloe_simu::EndCapModuleAndLayer(TString& str, TString& str2, int& modID,
                                     int& planeID)
{
  TObjArray* obja = str.Tokenize("_");  // ENDCAP => endvolECALActiveSlab_0_PV_0
  TObjArray* obja2 = str2.Tokenize("_");  // ENDCAP => ECAL_end_lv_PV_0

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
}

// find barrel cell from (x,y,z)
void kloe_simu::BarrelCell(double x, double y, double z, TGeoManager* g,
                           TGeoNode* node, int& cellID, double& d1, double& d2)
{
  double Pmaster[3];
  double Plocal[3];
  Pmaster[0] = x;
  Pmaster[1] = y;
  Pmaster[2] = z;

  g->GetCurrentNavigator()->MasterToLocal(Pmaster, Plocal);

  TGeoTrd2* trd = (TGeoTrd2*)node->GetVolume()->GetShape();

  if (kloe_simu::debug) {
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
  double dx = 0.5 * Plocal[2] / dz * (dx2 - dx1) + 0.5 * (dx2 + dx1);

  // Cell width at z = Plocal[2]
  double cellw = 2. * dx / kloe_simu::nCel;

  // cellID = distanza dall'estremo diviso larghezza cella
  cellID = (Plocal[0] + dx) / cellw;
}

// find endcap cell from (x,y,z)
void kloe_simu::EndCapCell(double x, double y, double z, TGeoManager* g,
                           TGeoNode* node, int& cellID, double& d1, double& d2)
{
  double Pmaster[3];
  double Plocal[3];
  Pmaster[0] = x;
  Pmaster[1] = y;
  Pmaster[2] = z;

  g->GetCurrentNavigator()->MasterToLocal(Pmaster, Plocal);

  TGeoTube* tub = (TGeoTube*)node->GetVolume()->GetShape();

  if (kloe_simu::debug) {
    std::cout << "pointer: " << tub << std::endl;
  }

  double rmin = tub->GetRmin();
  double rmax = tub->GetRmax();
  double dz = tub->GetDz();

  // d1 distanza da estremo up (y>0)
  // d2 distanza da estremo down (y<0)
  d1 = rmax * TMath::Sin(TMath::ACos(Plocal[0] / rmax)) - Plocal[1];
  d2 = rmax * TMath::Sin(TMath::ACos(Plocal[0] / rmax)) + Plocal[1];

  cellID = int((Plocal[0] / rmax + 1.) * kloe_simu::nCel_ec * 0.5);
}

// check if the geometry path is ok
bool kloe_simu::CheckAndProcessPath(TString& str2)
{
  // ENDCAP ==> something like:
  // "/volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_lv_PV_18/volECALActiveSlab_21_PV_0"
  // BARREL ==> something like:
  // "/volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_end_lv_PV_0/endvolECALActiveSlab_0_PV_0"
  TObjArray* obj = str2.Tokenize("/");

  int size = obj->GetEntries();
  if (size < 8) {
    return false;
  };

  // BARREL => ECAL_lv_PV_18
  // ENDCAP => ECAL_end_lv_PV_0
  str2 = ((TObjString*)obj->At(6))->GetString();
  delete obj;

  return true;
}

// get cell center from module id, layer id and cell id
void kloe_simu::CellPosition(TGeoManager* geo, int mod, int lay, int cel,
                             double& x, double& y, double& z)
{
  x = 0;
  y = 0;
  z = 0;

  double dummyLoc[3];
  double dummyMas[3];

  if (mod < 24) {

    if (kloe_simu::flukatype == false) {
      dummyLoc[0] = kloe_simu::cxlay[lay][cel];
      dummyLoc[1] = 0.;
      dummyLoc[2] = kloe_simu::czlay[lay];

      geo->cd(TString::Format(kloe_simu::path_barrel_template, mod).Data());
    } else {
      // Local coordinates calculation
      dummyLoc[0] = kloe_simu::cellCoordBarrel[mod][lay][cel][0];
      dummyLoc[1] = kloe_simu::cellCoordBarrel[mod][lay][cel][1];
      dummyLoc[2] = kloe_simu::cellCoordBarrel[mod][lay][cel][2];

      // Transformation to global coordinates
      dummyMas[0] = LocalToGlobalCoordinates(dummyLoc).X();
      dummyMas[1] = LocalToGlobalCoordinates(dummyLoc).Y();
      dummyMas[2] = LocalToGlobalCoordinates(dummyLoc).Z();
    }
  } else if (mod == 30 || mod == 40)
      // right x > 0 : c->mod = 30
      // left  x < 0 : c->mod = 40
  {

    if (kloe_simu::flukatype == false) {

      dummyLoc[0] = 2 * kloe_simu::ec_r / kloe_simu::nCel_ec * (0.5 + cel) -
                    kloe_simu::ec_r;
      dummyLoc[1] = 0.;
      dummyLoc[2] = kloe_simu::czlay[lay];

      if (mod == 30)
        geo->cd(kloe_simu::path_endcapR_template);
      else if (mod == 40)
        geo->cd(kloe_simu::path_endcapL_template);

    } else {
      // Local coordinates calculation
      dummyLoc[0] = kloe_simu::cellCoordEndcap[int(mod / 10)][lay][cel][0];
      dummyLoc[1] = kloe_simu::cellCoordEndcap[int(mod / 10)][lay][cel][1];
      dummyLoc[2] = kloe_simu::cellCoordEndcap[int(mod / 10)][lay][cel][2];

      // Transformation to global coordinates
      dummyMas[0] = LocalToGlobalCoordinates(dummyLoc).X();
      dummyMas[1] = LocalToGlobalCoordinates(dummyLoc).Y();
      dummyMas[2] = LocalToGlobalCoordinates(dummyLoc).Z();
    }
  }

  if (kloe_simu::flukatype == false) geo->LocalToMaster(dummyLoc, dummyMas);

  x = dummyMas[0];
  y = dummyMas[1];
  z = dummyMas[2];
}

// init geometry
// - costruct calo cells
// - find straw tube center
void kloe_simu::init(TGeoManager* geo)
{
  // https://root.cern.ch/root/htmldoc/guides/users-guide/Geometry.html#shapes
  // GetDx1() half length in x at -Dz
  // GetDx2() half length in x at +Dz
  // Dx1 < Dx2 => -Dz corresponds to minor width => internal side
  double xmin;
  double xmax;
  double dz;

  if (kloe_simu::flukatype == false) {
    TGeoTrd2* mod = (TGeoTrd2*)geo->FindVolumeFast("ECAL_lv_PV")->GetShape();
    xmin = mod->GetDx1();
    xmax = mod->GetDx2();
    dz = mod->GetDz();

    if ((abs(xmin - kloe_simu::xmin_f) > 0.2) ||
        (abs(xmax - kloe_simu::xmax_f) > 0.2) ||
        (abs(dz - kloe_simu::dz_f) > 0.2)) {
      std::cout << "ERROR ON ECAL GEOMETRY: xmin= " << xmin
                << " instead of what is expected in Fluka" << kloe_simu::xmin_f
                << std::endl;
      std::cout << "ERROR ON ECAL GEOMETRY: xmax= " << xmax
                << " instead of what is expected in Fluka" << kloe_simu::xmax_f
                << std::endl;
      std::cout << "ERROR ON ECAL GEOMETRY: dz= " << dz
                << " instead of what is expected in Fluka" << kloe_simu::dz_f
                << std::endl;
      // exit(1);
    }

  } else {
    xmin = kloe_simu::xmin_f;
    xmax = kloe_simu::xmax_f;
    dz = kloe_simu::dz_f;
  }

  double m = 0.5 * (xmax - xmin) / dz;
  double q = 0.5 * (xmax + xmin);

  // z edge of the cells
  double zlevel[kloe_simu::nLay + 1];
  zlevel[0] = -dz;

  for (int i = 0; i < kloe_simu::nLay; i++) {
    zlevel[i + 1] = zlevel[i] + kloe_simu::dzlay[i];
  }

  // z position of the center of the cells
  for (int i = 0; i < kloe_simu::nLay; i++) {
    kloe_simu::czlay[i] = 0.5 * (zlevel[i] + zlevel[i + 1]);

    // total module width at the z position of the center of the cell
    double xwidth = 2 * (m * kloe_simu::czlay[i] + q);

    // cell width at the z position of the center of the cell
    double dx = xwidth / kloe_simu::nCel;

    // x position of the center of the cells
    for (int j = 0; j < kloe_simu::nCel; j++) {
      kloe_simu::cxlay[i][j] = dx * (j + 0.5) - xwidth * 0.5;
    }
  }

  if (kloe_simu::flukatype == false) {
    TGeoTube* ec = (TGeoTube*)geo->FindVolumeFast("ECAL_end_lv_PV")->GetShape();
    kloe_simu::ec_r = ec->GetRmax();  // Maximum radius = 2000
    kloe_simu::ec_dz = ec->GetDz();   // half of thickness = 115

    TGeoHMatrix mat = *gGeoIdentity;

    kloe_simu::rST = new TPRegexp(rST_string);
    kloe_simu::rSTplane = new TPRegexp(rSTplane_string);

    getSTPlaneinfo(gGeoManager->GetTopVolume()->GetNode(0), mat, kloe_simu::stX,
                   kloe_simu::stL, kloe_simu::stPos);

    for (std::map<int, std::map<int, TVector2> >::iterator it =
             kloe_simu::stPos.begin();
         it != kloe_simu::stPos.end(); it++) {
      int mod_id = it->first;
      for (std::map<int, TVector2>::iterator ite = it->second.begin();
           ite != it->second.end(); ite++) {
        int tub_id = ite->first;
        int id = encodeSTID(mod_id, tub_id);
        tubePos[id] = ite->second;
      }
    }

    geo->cd(kloe_simu::path_internal_volume);
    double master[3];
    geo->LocalToMaster(kloe_simu::stt_center, master);

    if (abs(kloe_simu::ec_r - kloe_simu::ec_rf) > 0.2 ||
        (abs(kloe_simu::ec_dz - kloe_simu::ec_dzf))) {
      std::cout << "ERROR ON ECAL ENDCAP GEOMETRY: R= " << kloe_simu::ec_r
                << " instead of what is expected in Fluka" << kloe_simu::ec_rf
                << std::endl;
      std::cout << "ERROR ON ECAL ENDCAP GEOMETRY: Thickness= "
                << kloe_simu::ec_dz << " instead of what is expected in Fluka"
                << kloe_simu::ec_dzf << std::endl;
      //  exit(1);
    }
  }
}

// evaluated t0 for each tube assuming the beam bucket that
// produced the neutrino is known
void kloe_simu::initT0(TG4Event* ev)
{
  TRandom3 r(0);
  kloe_simu::t0.clear();

  double t0_beam = ev->Primaries[0].Position.T() -
                   ev->Primaries[0].Position.Z() / kloe_simu::c +
                   r.Gaus(0, kloe_simu::bucket_rms);

  for (std::map<int, std::map<int, TVector2> >::iterator it = stPos.begin();
       it != stPos.end(); it++) {
    kloe_simu::t0[it->first] =
        t0_beam + it->second.begin()->second.X() / kloe_simu::c;
  }
}

int kloe_simu::EncodeID(int mod, int lay, int cel)
{
  return cel + 100 * lay + 1000 * mod;
}

void kloe_simu::DecodeID(int id, int& mod, int& lay, int& cel)
{
  mod = id / 1000;
  lay = (id - mod * 1000) / 100;
  cel = id - mod * 1000 - lay * 100;
}

// evaluate minimum distance between segment (s1x,s1y,s1z) -> (s2x,s2y,s2z)
// and point (px,py,pz)
double kloe_simu::mindist(double s1x, double s1y, double s1z, double s2x,
                          double s2y, double s2z, double px, double py,
                          double pz)
{
  double segmod = (s1x - s2x) * (s1x - s2x) + (s1y - s2y) * (s1y - s2y) +
                  (s1z - s2z) * (s1z - s2z);

  double prod = (px - s1x) * (s2x - s1x) + (py - s1y) * (s2y - s1y) +
                (pz - s1z) * (s2z - s1z);

  double t = std::min(std::max(prod / segmod, 0.), 1.);

  double s3x = s1x + (s2x - s1x) * t;
  double s3y = s1y + (s2y - s1y) * t;
  double s3z = s1z + (s2z - s1z) * t;

  return sqrt((px - s3x) * (px - s3x) + (py - s3y) * (py - s3y) +
              (pz - s3z) * (pz - s3z));
}

// evaluate angle between (x1,y1,z1) and (x2,y2,z2)
double kloe_simu::angle(double x1, double y1, double z1, double x2, double y2,
                        double z2)
{
  double prod = x1 * x2 + y1 * y2 + z1 * z2;
  double mag1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
  double mag2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);

  return TMath::ACos(prod / (mag1 * mag2));
}

// get fiber attenuation factor.
// It depends on distance from pmt (d)
// and planeID (planes have different fibers)
double kloe_simu::AttenuationFactor(double d, int planeID)
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
  double atl2 = 0.0;

  switch (planeID) {
    case 0:
    case 1:
      atl2 = kloe_simu::atl2_01;
      break;

    case 2:
      atl2 = kloe_simu::atl2_2;
      break;

    case 3:
    case 4:
      atl2 = kloe_simu::atl2_34;
      break;

    default:
      // std::cout << "planeID out if range" << std::endl;
      atl2 = -999.0;
      break;
  }

  if (kloe_simu::debug) {
    std::cout << "planeID = " << planeID << std::endl;
    std::cout << "\tp1   = " << kloe_simu::p1 << std::endl;
    std::cout << "\talt1 = " << kloe_simu::atl1 << std::endl;
    std::cout << "\talt2 = " << atl2 << std::endl;
    std::cout << "\tatt  = "
              << kloe_simu::p1* TMath::Exp(-d / kloe_simu::atl1) +
                     (1. - kloe_simu::p1) * TMath::Exp(-d / atl2) << std::endl;
  }

  return kloe_simu::p1 * TMath::Exp(-d / kloe_simu::atl1) +
         (1. - kloe_simu::p1) * TMath::Exp(-d / atl2);
}

// reconstruct t of the hit from tdc1 and tdc2
double kloe_simu::TfromTDC(double t1, double t2, double L)
{
  return 0.5 * (t1 + t2 - kloe_simu::vlfb * L / kloe_simu::m_to_mm);
}

// reconstruct longitudinal coordinate of the hit from tdc1 and tdc2
double kloe_simu::XfromTDC(double t1, double t2)
{
  return 0.5 * (t1 - t2) / kloe_simu::vlfb * kloe_simu::m_to_mm;
}

// energy deposit of the hit from adc1 and adc2 and
// reconstructed longidutinal coordinate
double kloe_simu::EfromADC(double adc1, double adc2, double d1, double d2,
                           int planeID)
{
  double f1 = AttenuationFactor(d1, planeID);
  double f2 = AttenuationFactor(d2, planeID);

  return 0.5 * (adc1 / f1 + adc2 / f2) * kloe_simu::adc2MeV;
}

// reconstruct hit position, time and energy of the cell
void kloe_simu::CellXYZTE(dg_cell c, double& x, double& y, double& z, double& t,
                          double& e)
{
  if (c.id < 25000)  // Barrel
  {
    x = c.x + XfromTDC(c.ps1.tdc.at(0), c.ps2.tdc.at(0));
    y = c.y;
  } else {
    x = c.x;
    y = c.y - XfromTDC(c.ps1.tdc.at(0), c.ps2.tdc.at(0));
  }
  double d1 = 0.5 * c.l + XfromTDC(c.ps1.tdc.at(0), c.ps2.tdc.at(0));
  double d2 = 0.5 * c.l - XfromTDC(c.ps1.tdc.at(0), c.ps2.tdc.at(0));
  z = c.z;
  t = TfromTDC(c.ps1.tdc.at(0), c.ps2.tdc.at(0), c.l);
  e = EfromADC(c.ps1.adc.at(0), c.ps2.adc.at(0), d1, d2, c.lay);
}
