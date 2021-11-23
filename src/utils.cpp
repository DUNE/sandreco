#include "utils.h"
#include "struct.h"

#include "TG4Event.h"
#include "TG4HitSegment.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TGeoManager.h>
#include <TGeoTrd2.h>
#include <TGeoTube.h>
#include <TObjString.h>
#include <TRandom3.h>

#include "transf.h"
#include <iostream>

namespace sand_reco
{
bool flukatype = false;

double cellCoordBarrel[nMod][nLay][nCel][3];
double cellCoordEndcap[5][nLay][90][3];

double czlay[nLay];
double cxlay[nLay][nCel];
double ec_r;
double ec_dz;

double stt_center[3];

TPRegexp* rST;
TPRegexp* r2ST;
TPRegexp* rSTplane;
TPRegexp* rSTmod;

std::map<int, std::map<double, int> > stX;
std::map<int, double> stL;
std::map<int, std::map<int, TVector2> > stPos;
std::map<int, TVector2> tubePos;
std::map<int, double> t0;
}  // namespace sand_reco

bool sand_reco::isPeBefore(const pe& p1, const pe& p2)
{
  return p1.time < p2.time;
}

// value of parameter of segment (y1,z1,y2,z2)
// corresponding to minimal distance to point (y,z)
double sand_reco::getT(double y1, double y2, double y, double z1, double z2,
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

int sand_reco::encodeSTID(int planeid, int tubeid)
{
  return tubeid * 100000 + planeid;
}

void sand_reco::decodeSTID(int id, int& planeid, int& tubeid)
{
  tubeid = id / 100000;
  planeid = id % 100000;  // global id
}

int sand_reco::encodePlaneID(int moduleid, int planelocid, int type)
{
  return moduleid * 100 + planelocid * 10 + type;
}

void sand_reco::decodePlaneID(int id, int& moduleid, int& planelodid, int& type)
{
  moduleid = id / 100;
  planelodid = (id - moduleid * 100) / 10;
  type = id % 10;
}

bool sand_reco::isST(TString name) { return name.Contains(*sand_reco::rST); }

bool sand_reco::isSTPlane(TString name)
{
  return name.Contains(*sand_reco::rSTplane);
}

// get local id of tube
// int sand_reco::getSTId(TString name)
// {
//   int id = -999;

//   TObjArray* obj = name.Tokenize("_");

//   if (obj->GetEntries() != 10 && obj->GetEntries() != 9) {
//     std::cout << "Error: tokenizing " << name.Data() << std::endl;
//   } else {
//     TString sid = ((TObjString*)obj->At(obj->GetEntries() - 1))->GetString();

//     id = sid.Atoi();
//   }
//   delete obj;

//   return id;
// }

// get plane id
int sand_reco::getPlaneID(TString path)
{
  auto obja = sand_reco::rSTplane->MatchS(path);
  auto obja2 = sand_reco::rSTmod->MatchS(path);

  if(obja->GetEntries() == 0)
  {
    delete obja;
    delete obja2;
    return 0;
  }

  int mod = (reinterpret_cast<TObjString*>(obja->At(1)))->GetString().Atoi();
  int icopy = (reinterpret_cast<TObjString*>(obja->At(5)))->GetString().Atoi();
  int type =
      ((reinterpret_cast<TObjString*>(obja->At(4)))->GetString().EqualTo("hh")
           ? 2
           : 1) +
      2 * icopy;
  int icopymod =
      (reinterpret_cast<TObjString*>(obja2->At(3)))->GetString().Atoi();

  delete obja;
  delete obja2;

  // int mod = obja;
  // int type = 0;

  // TObjArray* obj = name.Tokenize("_");

  // if (obj->GetEntries() != 7 && obj->GetEntries() != 6 &&
  //     obj->GetEntries() != 11 && obj->GetEntries() != 10) {
  //   std::cout << "Error: tokenizing " << name.Data() << std::endl;
  // } else {
  //   int ipos = (obj->GetEntries() == 10 || obj->GetEntries() == 6) ? 2 : 3;
  //   TString stype = ((TObjString*)obj->At(ipos))->GetString();
  //   TString smod = ((TObjString*)obj->At(1))->GetString();

  //   mod = smod.Atoi();

  //   if (stype.Contains("hor2"))
  //     type = 4;
  //   else if (stype.Contains("ver"))
  //     type = 1;
  //   else if (stype.Contains("hor"))
  //     type = 2;
  //   else
  //     std::cout << "Error evaluating type for: " << name.Data() << std::endl;
  // }

  // delete obj;

  return encodePlaneID(mod * 10 + icopymod, icopy, type);
}

// get position of the center of the tube for a plane
void sand_reco::getSTinfo(TGeoNode* nod, TGeoHMatrix mat, int pid,
                          std::map<double, int>& stX,
                          std::map<int, double>& stL,
                          std::map<int, TVector2>& stPos)
{
  int type;
  int mod;
  int plloc;
  decodePlaneID(pid, mod, plloc, type);
  int ic = 1 - (type % 2);

  if (ic != 0 && ic != 1)
    std::cout << "Error: ic expected 0 or 1 -> " << ic << std::endl;

  for (int i = 0; i < nod->GetNdaughters(); i++) {
    auto n2straw = nod->GetDaughter(i);
    auto obja = sand_reco::r2ST->MatchS(n2straw->GetName());

    int n2straw_id =
        (reinterpret_cast<TObjString*>(obja->At(5)))->GetString().Atoi();
    delete obja;

    TGeoMatrix* n2strawmat = n2straw->GetMatrix();
    TGeoHMatrix n2strawhmat = mat * (*n2strawmat);

    for (int j = 0; j < n2straw->GetNdaughters(); j++) {
      TGeoNode* dau = n2straw->GetDaughter(j);
      TGeoTube* tub = (TGeoTube*)dau->GetVolume()->GetShape();
      double lenght = 2 * tub->GetDz();
      TString name = dau->GetName();

      if (!isST(name))
        std::cout << "Error: expected ST but not -> " << name.Data()
                  << std::endl;

      TGeoMatrix* thismat = n2straw->GetDaughter(j)->GetMatrix();
      TGeoHMatrix mymat = n2strawhmat * (*thismat);

      auto obja2 = sand_reco::rST->MatchS(name);

      int tid =
          (reinterpret_cast<TObjString*>(obja2->At(obja2->GetEntries() - 2)))
              ->GetString()
              .Atoi();
      delete obja2;

      int id = n2straw_id * 2 + tid;

      TVector2 v;
      v.SetX(mymat.GetTranslation()[2]);
      v.SetY(mymat.GetTranslation()[ic]);

      stX[v.Y()] = id;
      stL[id] = lenght;
      stPos[id] = v;
    }
  }
}

// get position of the center of the tube for each plane
void sand_reco::getSTPlaneinfo(TGeoHMatrix mat,
                               std::map<int, std::map<double, int> >& stX,
                               std::map<int, double>& stL,
                               std::map<int, std::map<int, TVector2> >& stPos)
{
  TGeoNode* nod = gGeoManager->GetCurrentNode();
  TString path = gGeoManager->GetPath();
  TString name = nod->GetName();
  TGeoMatrix* thismat = nod->GetMatrix();
  TGeoHMatrix mymat = mat * (*thismat);

  int pid = 0;
  double x = 0;
  double z = 0;

  if (isSTPlane(name)) {
    pid = getPlaneID(path);

    std::map<double, int> mstX;
    std::map<int, TVector2> mstPos;

    getSTinfo(nod, mymat, pid, mstX, stL, mstPos);

    stX[pid] = mstX;
    stPos[pid] = mstPos;

  } else {

    for (int i = 0; i < nod->GetNdaughters(); i++) {
      gGeoManager->CdDown(i);
      getSTPlaneinfo(mymat, stX, stL, stPos);
      gGeoManager->CdUp();
    }
  }
}

// get unique id of the tube corresponding to (x,y,z)
int sand_reco::getSTUniqID(TGeoManager* g, double x, double y, double z)
{
  int sid = -999;
  int pid = getPlaneID(g->GetPath());

  if (pid == 0) return -999;

  double xt = 0.;

  if (pid % 2 == 1)
    xt = x;
  else
    xt = y;

  std::map<double, int>::iterator it = sand_reco::stX.at(pid).lower_bound(xt);

  if (it == sand_reco::stX.at(pid).begin()) {
    sid = sand_reco::stX.at(pid).begin()->second;
  } else if (it == sand_reco::stX.at(pid).end()) {
    sid = sand_reco::stX.at(pid).rbegin()->second;
  } else {
    TVector2 v1 = sand_reco::stPos.at(pid).at(it->second);
    TVector2 v2 = sand_reco::stPos.at(pid).at(std::prev(it)->second);

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

bool sand_reco::isCluBigger(const std::vector<dg_tube>& v1,
                            const std::vector<dg_tube>& v2)
{
  return v1.size() > v2.size();
}

bool sand_reco::isDigUpstream(const dg_tube& d1, const dg_tube& d2)
{
  return d1.z < d2.z;
}

bool sand_reco::isHitBefore(hit h1, hit h2) { return h1.t1 < h2.t1; }

bool sand_reco::isDigBefore(dg_tube d1, dg_tube d2) { return d1.tdc < d2.tdc; }

bool sand_reco::isCellBefore(dg_cell c1, dg_cell c2)
{
  if (c1.ps1.size() == 0 || c1.ps2.size() == 0)
    return false;
  else if (c2.ps1.size() == 0 || c2.ps2.size() == 0)
    return true;
  else
    return ((c1.ps1.at(0).tdc + c1.ps2.at(0).tdc) <
            (c2.ps1.at(0).tdc + c2.ps2.at(0).tdc));
}

bool sand_reco::isAfter(particle p1, particle p2) { return p1.tid > p2.tid; }

bool sand_reco::isBarrel(TString& str)
{
  // something like: volECALActiveSlab_21_PV_0
  return str.Contains("volECAL") == true && str.Contains("Active") == true &&
         str.Contains("end") == false;
}

bool sand_reco::isEndCap(TString& str)
{
  // something like: endvolECALActiveSlab_0_PV_0
  return str.Contains("endvolECAL") == true && str.Contains("Active") == true;
}

// find module id and plane id for barrel
void sand_reco::BarrelModuleAndLayer(TString& str, TString& str2, int& modID,
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
void sand_reco::EndCapModuleAndLayer(TString& str, TString& str2, int& modID,
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
void sand_reco::BarrelCell(double x, double y, double z, TGeoManager* g,
                           TGeoNode* node, int& cellID, double& d1, double& d2)
{
  double Pmaster[3];
  double Plocal[3];
  Pmaster[0] = x;
  Pmaster[1] = y;
  Pmaster[2] = z;

  g->GetCurrentNavigator()->MasterToLocal(Pmaster, Plocal);

  TGeoTrd2* trd = (TGeoTrd2*)node->GetVolume()->GetShape();

  if (sand_reco::debug) {
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
  double cellw = 2. * dx / sand_reco::nCel;

  // cellID = distanza dall'estremo diviso larghezza cella
  cellID = (Plocal[0] + dx) / cellw;
}

// find endcap cell from (x,y,z)
void sand_reco::EndCapCell(double x, double y, double z, TGeoManager* g,
                           TGeoNode* node, int& cellID, double& d1, double& d2)
{
  double Pmaster[3];
  double Plocal[3];
  Pmaster[0] = x;
  Pmaster[1] = y;
  Pmaster[2] = z;

  g->GetCurrentNavigator()->MasterToLocal(Pmaster, Plocal);

  TGeoTube* tub = (TGeoTube*)node->GetVolume()->GetShape();

  if (sand_reco::debug) {
    std::cout << "pointer: " << tub << std::endl;
  }

  double rmin = tub->GetRmin();
  double rmax = tub->GetRmax();
  double dz = tub->GetDz();

  // d1 distanza da estremo up (y>0)
  // d2 distanza da estremo down (y<0)
  d1 = rmax * TMath::Sin(TMath::ACos(Plocal[0] / rmax)) - Plocal[1];
  d2 = rmax * TMath::Sin(TMath::ACos(Plocal[0] / rmax)) + Plocal[1];

  cellID = int((Plocal[0] / rmax + 1.) * sand_reco::nCel_ec * 0.5);
}

// check if the geometry path is ok
bool sand_reco::CheckAndProcessPath(TString& str2)
{
  // ENDCAP ==> something like:
  // "/volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_lv_PV_18/volECALActiveSlab_21_PV_0"
  // BARREL ==> something like:
  // "/volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_end_lv_PV_0/endvolECALActiveSlab_0_PV_0"
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
void sand_reco::CellPosition(TGeoManager* geo, int mod, int lay, int cel,
                             double& x, double& y, double& z)
{
  x = 0;
  y = 0;
  z = 0;

  double dummyLoc[3];
  double dummyMas[3];

  if (mod < 24) {

    if (sand_reco::flukatype == false) {
      dummyLoc[0] = sand_reco::cxlay[lay][cel];
      dummyLoc[1] = 0.;
      dummyLoc[2] = sand_reco::czlay[lay];

      geo->cd(TString::Format(sand_reco::path_barrel_template, mod).Data());
    } else {
      // Local coordinates calculation
      dummyLoc[0] = sand_reco::cellCoordBarrel[mod][lay][cel][0];
      dummyLoc[1] = sand_reco::cellCoordBarrel[mod][lay][cel][1];
      dummyLoc[2] = sand_reco::cellCoordBarrel[mod][lay][cel][2];

      // Transformation to global coordinates
      dummyMas[0] = LocalToGlobalCoordinates(dummyLoc).X();
      dummyMas[1] = LocalToGlobalCoordinates(dummyLoc).Y();
      dummyMas[2] = LocalToGlobalCoordinates(dummyLoc).Z();
    }
  } else if (mod == 30 || mod == 40)
  // right x > 0 : c->mod = 30
  // left  x < 0 : c->mod = 40
  {

    if (sand_reco::flukatype == false) {

      dummyLoc[0] = 2 * sand_reco::ec_r / sand_reco::nCel_ec * (0.5 + cel) -
                    sand_reco::ec_r;
      dummyLoc[1] = 0.;
      dummyLoc[2] = sand_reco::czlay[lay];

      if (mod == 30)
        geo->cd(sand_reco::path_endcapR_template);
      else if (mod == 40)
        geo->cd(sand_reco::path_endcapL_template);

    } else {
      // Local coordinates calculation
      dummyLoc[0] = sand_reco::cellCoordEndcap[int(mod / 10)][lay][cel][0];
      dummyLoc[1] = sand_reco::cellCoordEndcap[int(mod / 10)][lay][cel][1];
      dummyLoc[2] = sand_reco::cellCoordEndcap[int(mod / 10)][lay][cel][2];

      // Transformation to global coordinates
      dummyMas[0] = LocalToGlobalCoordinates(dummyLoc).X();
      dummyMas[1] = LocalToGlobalCoordinates(dummyLoc).Y();
      dummyMas[2] = LocalToGlobalCoordinates(dummyLoc).Z();
    }
  }

  if (sand_reco::flukatype == false) geo->LocalToMaster(dummyLoc, dummyMas);

  x = dummyMas[0];
  y = dummyMas[1];
  z = dummyMas[2];
}

// init geometry
// - costruct calo cells
// - find straw tube center
void sand_reco::init(TGeoManager* geo)
{
  // https://root.cern.ch/root/htmldoc/guides/users-guide/Geometry.html#shapes
  // GetDx1() half length in x at -Dz
  // GetDx2() half length in x at +Dz
  // Dx1 < Dx2 => -Dz corresponds to minor width => internal side
  double xmin;
  double xmax;
  double dz;

  if (sand_reco::flukatype == false) {
    TGeoTrd2* mod = (TGeoTrd2*)geo->FindVolumeFast("ECAL_lv_PV")->GetShape();
    xmin = mod->GetDx1();
    xmax = mod->GetDx2();
    dz = mod->GetDz();

    if ((abs(xmin - sand_reco::xmin_f) > 0.2) ||
        (abs(xmax - sand_reco::xmax_f) > 0.2) ||
        (abs(dz - sand_reco::dz_f) > 0.2)) {
      std::cout << "ERROR ON ECAL GEOMETRY: xmin= " << xmin
                << " instead of what is expected in Fluka" << sand_reco::xmin_f
                << std::endl;
      std::cout << "ERROR ON ECAL GEOMETRY: xmax= " << xmax
                << " instead of what is expected in Fluka" << sand_reco::xmax_f
                << std::endl;
      std::cout << "ERROR ON ECAL GEOMETRY: dz= " << dz
                << " instead of what is expected in Fluka" << sand_reco::dz_f
                << std::endl;
      // exit(1);
    }

  } else {
    xmin = sand_reco::xmin_f;
    xmax = sand_reco::xmax_f;
    dz = sand_reco::dz_f;
  }

  double m = 0.5 * (xmax - xmin) / dz;
  double q = 0.5 * (xmax + xmin);

  // z edge of the cells
  double zlevel[sand_reco::nLay + 1];
  zlevel[0] = -dz;

  for (int i = 0; i < sand_reco::nLay; i++) {
    zlevel[i + 1] = zlevel[i] + sand_reco::dzlay[i];
  }

  // z position of the center of the cells
  for (int i = 0; i < sand_reco::nLay; i++) {
    sand_reco::czlay[i] = 0.5 * (zlevel[i] + zlevel[i + 1]);

    // total module width at the z position of the center of the cell
    double xwidth = 2 * (m * sand_reco::czlay[i] + q);

    // cell width at the z position of the center of the cell
    double dx = xwidth / sand_reco::nCel;

    // x position of the center of the cells
    for (int j = 0; j < sand_reco::nCel; j++) {
      sand_reco::cxlay[i][j] = dx * (j + 0.5) - xwidth * 0.5;
    }
  }

  if (sand_reco::flukatype == false) {
    TGeoTube* ec = (TGeoTube*)geo->FindVolumeFast("ECAL_end_lv_PV")->GetShape();
    sand_reco::ec_r = ec->GetRmax();  // Maximum radius = 2000
    sand_reco::ec_dz = ec->GetDz();   // half of thickness = 115

    TGeoHMatrix mat = *gGeoIdentity;

    sand_reco::rST = new TPRegexp(rST_string);
    sand_reco::r2ST = new TPRegexp(r2ST_string);
    sand_reco::rSTplane = new TPRegexp(rSTplane_string);
    sand_reco::rSTmod = new TPRegexp(rSTmod_string);

    gGeoManager->CdDown(0);

    getSTPlaneinfo(mat, sand_reco::stX, sand_reco::stL, sand_reco::stPos);

    for (std::map<int, std::map<int, TVector2> >::iterator it =
             sand_reco::stPos.begin();
         it != sand_reco::stPos.end(); it++) {
      int mod_id = it->first;
      for (std::map<int, TVector2>::iterator ite = it->second.begin();
           ite != it->second.end(); ite++) {
        int tub_id = ite->first;
        int id = encodeSTID(mod_id, tub_id);
        tubePos[id] = ite->second;
      }
    }

    geo->cd(sand_reco::path_internal_volume);
    double master[3];
    geo->LocalToMaster(sand_reco::stt_center, master);

    if (abs(sand_reco::ec_r - sand_reco::ec_rf) > 0.2 ||
        (abs(sand_reco::ec_dz - sand_reco::ec_dzf))) {
      std::cout << "ERROR ON ECAL ENDCAP GEOMETRY: R= " << sand_reco::ec_r
                << " instead of what is expected in Fluka" << sand_reco::ec_rf
                << std::endl;
      std::cout << "ERROR ON ECAL ENDCAP GEOMETRY: Thickness= "
                << sand_reco::ec_dz << " instead of what is expected in Fluka"
                << sand_reco::ec_dzf << std::endl;
      //  exit(1);
    }
  }
}

// evaluated t0 for each tube assuming the beam bucket that
// produced the neutrino is known
void sand_reco::initT0(TG4Event* ev)
{
  TRandom3 r(0);
  sand_reco::t0.clear();

  double t0_beam = ev->Primaries[0].Position.T() -
                   ev->Primaries[0].Position.Z() / sand_reco::c +
                   r.Gaus(0, sand_reco::bucket_rms);

  for (std::map<int, std::map<int, TVector2> >::iterator it = stPos.begin();
       it != stPos.end(); it++) {
    sand_reco::t0[it->first] =
        t0_beam + it->second.begin()->second.X() / sand_reco::c;
  }
}

int sand_reco::EncodeID(int mod, int lay, int cel)
{
  return cel + 100 * lay + 1000 * mod;
}

void sand_reco::DecodeID(int id, int& mod, int& lay, int& cel)
{
  mod = id / 1000;
  lay = (id - mod * 1000) / 100;
  cel = id - mod * 1000 - lay * 100;
}

// evaluate minimum distance between segment (s1x,s1y,s1z) -> (s2x,s2y,s2z)
// and point (px,py,pz)
double sand_reco::mindist(double s1x, double s1y, double s1z, double s2x,
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
double sand_reco::angle(double x1, double y1, double z1, double x2, double y2,
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
double sand_reco::AttenuationFactor(double d, int planeID)
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
      atl2 = sand_reco::atl2_01;
      break;

    case 2:
      atl2 = sand_reco::atl2_2;
      break;

    case 3:
    case 4:
      atl2 = sand_reco::atl2_34;
      break;

    default:
      // std::cout << "planeID out if range" << std::endl;
      atl2 = -999.0;
      break;
  }

  if (sand_reco::debug) {
    std::cout << "planeID = " << planeID << std::endl;
    std::cout << "\tp1   = " << sand_reco::p1 << std::endl;
    std::cout << "\talt1 = " << sand_reco::atl1 << std::endl;
    std::cout << "\talt2 = " << atl2 << std::endl;
    std::cout << "\tatt  = "
              << sand_reco::p1 * TMath::Exp(-d / sand_reco::atl1) +
                     (1. - sand_reco::p1) * TMath::Exp(-d / atl2)
              << std::endl;
  }

  return sand_reco::p1 * TMath::Exp(-d / sand_reco::atl1) +
         (1. - sand_reco::p1) * TMath::Exp(-d / atl2);
}

// reconstruct t of the hit from tdc1 and tdc2
double sand_reco::TfromTDC(double t1, double t2, double L)
{
  return 0.5 * (t1 + t2 - sand_reco::vlfb * L / sand_reco::m_to_mm);
}

// reconstruct longitudinal coordinate of the hit from tdc1 and tdc2
double sand_reco::XfromTDC(double t1, double t2)
{
  return 0.5 * (t1 - t2) / sand_reco::vlfb * sand_reco::m_to_mm;
}

// energy deposit of the hit from adc1 and adc2 and
// reconstructed longidutinal coordinate
double sand_reco::EfromADC(double adc1, double adc2, double d1, double d2,
                           int planeID)
{
  double f1 = AttenuationFactor(d1, planeID);
  double f2 = AttenuationFactor(d2, planeID);

  return 0.5 * (adc1 / f1 + adc2 / f2) * sand_reco::adc2MeV;
}

// reconstruct hit position, time and energy of the cell
void sand_reco::CellXYZTE(dg_cell c, double& x, double& y, double& z, double& t,
                          double& e)
{
  if (c.id < 25000)  // Barrel
  {
    x = c.x + XfromTDC(c.ps1.at(0).tdc, c.ps2.at(0).tdc);
    y = c.y;
  } else {
    x = c.x;
    y = c.y - XfromTDC(c.ps1.at(0).tdc, c.ps2.at(0).tdc);
  }
  double d1 = 0.5 * c.l + XfromTDC(c.ps1.at(0).tdc, c.ps2.at(0).tdc);
  double d2 = 0.5 * c.l - XfromTDC(c.ps1.at(0).tdc, c.ps2.at(0).tdc);
  z = c.z;
  t = TfromTDC(c.ps1.at(0).tdc, c.ps2.at(0).tdc, c.l);
  e = EfromADC(c.ps1.at(0).adc, c.ps2.at(0).adc, d1, d2, c.lay);
}
