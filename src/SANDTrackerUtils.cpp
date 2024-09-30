#include "SANDTrackerUtils.h"

TGeoManager* SANDTrackerUtils::fGeo = 0;
const double SANDTrackerUtils::kMagneticFieldInT = 0.6; /*T*/

const double SANDTrackerUtils::k = 0.299792458;
const double SANDTrackerUtils::c = SANDTrackerUtils::k * 1E3;  // mm/ns

void SANDTrackerUtils::Clear()
{
  // sand_reco::stt::stL.clear();
  // sand_reco::stt::stX.clear();
  // sand_reco::stt::stPos.clear();
  // sand_reco::stt::t0.clear();
  // sand_reco::stt::tubePos.clear();
}

// check if tubes are adjacent
// using tube id
bool SANDTrackerUtils::AreAdjacent(const SANDTrackerCellID &tub1, const SANDTrackerCellID &tub2)
{
  // To Do: this doesn't work for staggered stt tubes
  return (abs(long(tub1()) - long(tub2())) <= 1);
}

const SANDTrackerCellID SANDTrackerUtils::GetTubeID(const SANDTrackerDigitID &id)
{
  SANDTrackerPlaneID pln;
  SANDTrackerCellID tub;
  std::tie(pln, tub) = GetPlaneAndTubeIds(id);
  return tub;
}

const SANDTrackerPlaneID SANDTrackerUtils::GetPlaneID(const SANDTrackerDigitID &id)
{
  SANDTrackerPlaneID pln;
  SANDTrackerCellID tub;
  std::tie(pln, tub) = GetPlaneAndTubeIds(id);
  return pln;
}

const SANDTrackerPlaneLocalID SANDTrackerUtils::GetPlaneLocalID(const SANDTrackerPlaneID &id)
{
  SANDTrackerPlaneLocalID pln;
  SANDTrackerModuleID mid;
  std::tie(mid, pln) = GetModuleIdAndPlaneLocalId(id);
  return pln;
}

const SANDTrackerModuleID SANDTrackerUtils::GetModelID(const SANDTrackerPlaneID &id)
{
  SANDTrackerPlaneLocalID pln;
  SANDTrackerModuleID mid;
  std::tie(mid, pln) = GetModuleIdAndPlaneLocalId(id);
  return mid;
}

std::tuple<SANDTrackerPlaneID, SANDTrackerCellID> SANDTrackerUtils::GetPlaneAndTubeIds(
    const SANDTrackerDigitID &id)
{
  int pid, tid;
  // sand_reco::stt::decodeSTID(id(), pid, tid);
  // return std::make_tuple(SANDTrackerPlaneID(pid), SANDTrackerCellID(tid));
}

std::tuple<SANDTrackerModuleID, SANDTrackerPlaneLocalID> SANDTrackerUtils::GetModuleIdAndPlaneLocalId(
    const SANDTrackerPlaneID &id)
{
  int mid, lid, typ;
  // sand_reco::stt::decodePlaneID(id(), mid, lid, typ);
  // return std::make_tuple(SANDTrackerModuleID(mid), SANDTrackerPlaneLocalID(lid));
}

// init pattern reco
void SANDTrackerUtils::Init(TGeoManager *geo)
{
  InitGeo(geo);
  // SANDTrackerStrawTubeTracker::Init();
}

double SANDTrackerUtils::GetPerpMomentumInGeVFromRadiusInMM(double radius) {
  return GetRadiusInMMToMomentumInGeVConstant() * radius * GetMagneticField(); 
}

double SANDTrackerUtils::GetRadiusInMMFromPerpMomentumInGeV(double perpMom) {
  return perpMom / (GetRadiusInMMToMomentumInGeVConstant() * GetMagneticField()); 
}

TString SANDTrackerUtils::PrintMatrix(const TMatrixD& m) {
  TString str = "\n";
  for(int i = 0; i < m.GetNrows(); i++)
  {
    for(int j = 0; j < m.GetNcols(); j++)
    {
      str += TString::Format("%20.5f  ",m[i][j]);
    }
    str += "\n";
  }
  str += "\n";

  return str;
}

// TString SANDTrackerUtils::PrintStateVector(const SANDTrackerKFStateVector& v) {
//   TString str = "\n";
//   str += TString::Format("%20.5f\n",v.X());
//   str += TString::Format("%20.5f\n",v.Y());
//   str += TString::Format("%20.5f\n",v.SignedInverseRadius());
//   str += TString::Format("%20.5f\n",v.TanLambda());
//   str += TString::Format("%20.5f\n",v.Phi());
//   return str;
// }

// const SANDTrackerCluster* SANDTrackerUtils::GetClusterPointer(int clusterID, const std::vector<SANDTrackerCluster>& clusters)
// {
//   for(auto& cl: clusters)
//     if(cl.GetId() == clusterID)
//       return &cl;
  
//   return 0;
// }

TVector3 SANDTrackerUtils::GetCartesianCoordinateFromCylindrical(double radius, double angle, double x)
{
  auto sandCenter = SANDTrackerUtils::GetSANDInnerVolumeCenterPosition();
  auto xSandCenter = sandCenter[0];
  auto ySandCenter = sandCenter[1];
  auto zSandCenter = sandCenter[2];
  
  auto y = ySandCenter + radius * sin(angle);
  auto z = zSandCenter + radius * cos(angle);

  return {x,y,z};
}