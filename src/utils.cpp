#include "utils.h"
#include "struct.h"
#include "transf.h"

#include "TG4Event.h"
#include "TG4HitSegment.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TGeoManager.h>
#include <TGeoTrd2.h>
#include <TGeoTube.h>
#include <TObjString.h>
#include <TRandom3.h>

#include <iostream>

namespace sand_reco
{
std::map<int, double> t0;
namespace ecal
{
double acquisition::fixed_thresh_pe = 3.;

namespace endcap
{
double ec_r;
double ec_dz;
}  // namespace endcap

}  // namespace ecal

namespace stt
{

double stt_center[3];

}  // namespace stt

namespace fluka
{
namespace ecal
{
double czlay[nLay];
double cxlay[nLay][nCel];
double cellCoordBarrel[nMod][nLay][nCel][3];
double cellCoordEndcap[5][nLay][90][3];

}  // namespace ecal

namespace stt
{

std::map<int, std::map<double, int> > stX;
std::map<int, double> stL;
std::map<int, std::map<int, TVector2> > stPos;
std::map<int, TVector2> tubePos;

int encodeSTID(int planeid, int tubeid)
{
  return tubeid * 100000 + planeid;
}

int encodePlaneID(int moduleid, int planelocid, int type)
{
  return moduleid * 100 + planelocid * 10 + type;
}

void decodeSTID(int id, int& planeid, int& tubeid)
{
  tubeid = id / 100000;
  planeid = id % 100000;  // global id
}

void decodePlaneID(int id, int& moduleid, int& planelodid, int& type)
{
  moduleid = id / 100;
  planelodid = (id - moduleid * 100) / 10;
  type = id % 10;
}

}  // namespace stt

}  // namespace fluka

}  // namespace sand_reco

bool sand_reco::ecal::isPeBefore(const pe& p1, const pe& p2)
{
  return p1.time < p2.time;
}

// value of parameter of segment (y1,z1,y2,z2)
// corresponding to minimal distance to point (y,z)
double sand_reco::stt::getT(double y1, double y2, double y, double z1,
                            double z2, double z)
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

bool sand_reco::ecal::isCluBigger(const std::vector<dg_wire>& v1,
                                  const std::vector<dg_wire>& v2)
{
  return v1.size() > v2.size();
}

bool sand_reco::stt::isDigUpstream(const dg_wire& d1, const dg_wire& d2)
{
  return d1.z < d2.z;
}

bool sand_reco::ecal::isHitBefore(hit h1, hit h2)
{
  return h1.t1 < h2.t1;
}

bool sand_reco::stt::isDigBefore(dg_wire d1, dg_wire d2)
{
  return d1.tdc < d2.tdc;
}

bool sand_reco::ecal::isCellBefore(dg_cell c1, dg_cell c2)
{
  if (c1.ps1.size() == 0 || c1.ps2.size() == 0)
    return false;
  else if (c2.ps1.size() == 0 || c2.ps2.size() == 0)
    return true;
  else
    return ((c1.ps1.at(0).tdc + c1.ps2.at(0).tdc) <
            (c2.ps1.at(0).tdc + c2.ps2.at(0).tdc));
}

bool sand_reco::isAfter(particle p1, particle p2)
{
  return p1.tid > p2.tid;
}

// get cell center from module id, layer id and cell id
void sand_reco::fluka::ecal::CellPosition(TGeoManager* geo, int det, int mod,
                                          int lay, int cel, double& x,
                                          double& y, double& z)
{
  double dummyMas[3] = {0., 0., 0.};

  if (mod < 24) {

    double dummyLoc[3];

    // Local coordinates calculation
    dummyLoc[0] = cellCoordBarrel[mod][lay][cel][0];
    dummyLoc[1] = cellCoordBarrel[mod][lay][cel][1];
    dummyLoc[2] = cellCoordBarrel[mod][lay][cel][2];

    // Transformation to global coordinates
    dummyMas[0] = LocalToGlobalCoordinates(dummyLoc).X();
    dummyMas[1] = LocalToGlobalCoordinates(dummyLoc).Y();
    dummyMas[2] = LocalToGlobalCoordinates(dummyLoc).Z();

  } else if (mod == 30 || mod == 40)
      // right x > 0 : c->mod = 30
      // left  x < 0 : c->mod = 40
  {

    double dummyLoc[3];

    // Local coordinates calculation
    dummyLoc[0] = cellCoordEndcap[int(mod / 10)][lay][cel][0];
    dummyLoc[1] = cellCoordEndcap[int(mod / 10)][lay][cel][1];
    dummyLoc[2] = cellCoordEndcap[int(mod / 10)][lay][cel][2];

    // Transformation to global coordinates
    dummyMas[0] = LocalToGlobalCoordinates(dummyLoc).X();
    dummyMas[1] = LocalToGlobalCoordinates(dummyLoc).Y();
    dummyMas[2] = LocalToGlobalCoordinates(dummyLoc).Z();
  }

  x = dummyMas[0];
  y = dummyMas[1];
  z = dummyMas[2];
}

// init geometry
// - costruct calo cells
// - find straw tube center
void sand_reco::fluka::init(TGeoManager* geo)
{
  // https://root.cern.ch/root/htmldoc/guides/users-guide/Geometry.html#shapes
  // GetDx1() half length in x at -Dz
  // GetDx2() half length in x at +Dz
  // Dx1 < Dx2 => -Dz corresponds to minor width => internal side
  double xmin;
  double xmax;
  double dz;

  xmin = ecal::xmin_f;
  xmax = ecal::xmax_f;
  dz = ecal::dz_f;

  double m = 0.5 * (xmax - xmin) / dz;
  double q = 0.5 * (xmax + xmin);

  // z edge of the cells
  double zlevel[sand_reco::fluka::ecal::nLay + 1];
  zlevel[0] = -dz;

  for (int i = 0; i < sand_reco::fluka::ecal::nLay; i++) {
    zlevel[i + 1] = zlevel[i] + sand_reco::fluka::ecal::dzlay[i];
  }

  // z position of the center of the cells
  for (int i = 0; i < sand_reco::fluka::ecal::nLay; i++) {
    sand_reco::fluka::ecal::czlay[i] = 0.5 * (zlevel[i] + zlevel[i + 1]);

    // total module width at the z position of the center of the cell
    double xwidth = 2 * (m * sand_reco::fluka::ecal::czlay[i] + q);

    // cell width at the z position of the center of the cell
    double dx = xwidth / sand_reco::fluka::ecal::nCel;

    // x position of the center of the cells
    for (int j = 0; j < sand_reco::fluka::ecal::nCel; j++) {
      sand_reco::fluka::ecal::cxlay[i][j] = dx * (j + 0.5) - xwidth * 0.5;
    }
  }
}

int sand_reco::ecal::decoder::EncodeID(int det, int mod, int lay, int cel)
{
  return cel + 100 * lay + 1000 * mod + det * 1e7;
}

void sand_reco::ecal::decoder::DecodeID(int id, int& det, int& mod, int& lay,
                                        int& cel)
{
  det = id / 1e7;
  id -= det * 1e7;

  mod = id / 1000;
  id -= mod * 1000;

  lay = id / 100;
  id -= lay * 100;

  cel = id;
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
double sand_reco::ecal::attenuation::AttenuationFactor(double d, int planeID)
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
      atl2 = atl2_01;
      break;

    case 2:
      atl2 = atl2_2;
      break;

    case 3:
    case 4:
      atl2 = atl2_34;
      break;

    default:
      // std::cout << "planeID out if range" << std::endl;
      atl2 = -999.0;
      break;
  }

  if (sand_reco::debug) {
    std::cout << "planeID = " << planeID << std::endl;
    std::cout << "\tp1   = " << p1 << std::endl;
    std::cout << "\talt1 = " << atl1 << std::endl;
    std::cout << "\talt2 = " << atl2 << std::endl;
    std::cout << "\tatt  = "
              << p1* TMath::Exp(-d / atl1) + (1. - p1) * TMath::Exp(-d / atl2)
              << std::endl;
  }

  return p1 * TMath::Exp(-d / atl1) + (1. - p1) * TMath::Exp(-d / atl2);
}

// reconstruct t of the hit from tdc1 and tdc2
double sand_reco::ecal::reco::TfromTDC(double t1, double t2, double L)
{
  return 0.5 * (t1 + t2 - scintillation::vlfb * L / conversion::m_to_mm);
}

// reconstruct longitudinal coordinate of the hit from tdc1 and tdc2
double sand_reco::ecal::reco::XfromTDC(double t1, double t2)
{
  return 0.5 * (t1 - t2) / scintillation::vlfb * conversion::m_to_mm;
}

// energy deposit of the hit from adc1 and adc2 and
// reconstructed longidutinal coordinate
double sand_reco::ecal::reco::EfromADC(double adc1, double adc2, double d1,
                                       double d2, int planeID)
{
  double f1 = attenuation::AttenuationFactor(d1, planeID);
  double f2 = attenuation::AttenuationFactor(d2, planeID);

  //return 0.5 * (adc1 / f1 + adc2 / f2) * energy_calibration::adc2MeV;
  return 0.5 * (adc1 / f1 + adc2 / f2) / (energy_calibration::attpassratio * acquisition::pe2ADC * photo_sensor::e2pe);
}

// energy deposit of the hit from a single adc and 
// reconstructed longidutinal coordinate of a cluster
double sand_reco::ecal::reco::EfromADCsingle(double adc, double f)
{
  return adc / (f * energy_calibration::attpassratio * sand_reco::ecal::acquisition::pe2ADC *
                sand_reco::ecal::photo_sensor::e2pe);
}

// reconstruct hit position, time and energy of the cell
void sand_reco::ecal::reco::CellXYZTE(dg_cell c, double& x, double& y,
                                      double& z, double& t, double& e)
{
  if (c.id < 25000)  // Barrel
  {
    x = c.x - XfromTDC(c.ps1.at(0).tdc, c.ps2.at(0).tdc);
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
