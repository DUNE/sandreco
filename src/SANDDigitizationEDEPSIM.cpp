#include "SANDDigitizationEDEPSIM.h"
#include "SANDDigitization.h"

#include <iomanip>
#include <iostream>

#include <iomanip>

#include "TFile.h"
#include "TTree.h"

#include "utils.h"

using namespace sand_reco;
using namespace sand_geometry;

namespace digitization
{

namespace edep_sim
{
namespace tracker
{
  std::vector<TLorentzVector> WireHitClosestPoints(hit& h, SANDWireInfo& wire)
{
  std::vector<TLorentzVector> closestPoints;

  TVector3 start = {h.x1, h.y1, h.z1};  // hit start
  TVector3 stop = {h.x2, h.y2, h.z2};   // hit end
  TVector3 s = stop - start;            // hit direction

  TVector3 leftend;
  TVector3 rightend;
  if (wire.readout_end() == SANDWireInfo::ReadoutEnd::kFirst) {
    leftend  = wire.getFirstPoint();
    rightend = wire.getSecondPoint();
  } else {
    leftend  = wire.getSecondPoint();
    rightend = wire.getFirstPoint();
  }
  

  TVector3 r(rightend - leftend); // wire direction

  TVector3 d = start - leftend;
  double A = s.Dot(s);    // s . s
  double B = s.Dot(r);    // s . r
  double C = r.Dot(r);    // r . r
  double D = s.Dot(d);    // s . (start - leftend)
  double E = r.Dot(d);    // r . (start - leftend)

  double denominator = A * C - B * B;
  if (denominator != 0) {
    double t = (B * E - C * D) / denominator;
    double t_prime = (A * E - B * D) / denominator;

    t = std::max(0.0, std::min(1.0, t));
    t_prime = std::max(0.0, std::min(1.0, t_prime));

    TVector3 closest_point_hit = start + t * s;
    
    if (t == 0 || t == 1) {
      TVector3 AP = closest_point_hit - leftend; 
      t_prime = AP.Dot(r) / r.Mag2();
      t_prime = std::max(0.0, std::min(1.0, t_prime));
    }

    TVector3 closest_point_wire = leftend + t_prime * r;

    TLorentzVector closest_point_hit_l;
    double fraction = (closest_point_hit - start).Mag() / s.Mag();
    closest_point_hit_l.SetXYZT(closest_point_hit.X(), closest_point_hit.Y(), closest_point_hit.Z(), h.t1 + fraction * (h.t2 - h.t1));

    TLorentzVector closest_point_wire_l;
    closest_point_wire_l.SetXYZT(closest_point_wire.X(), closest_point_wire.Y(), closest_point_wire.Z(), 
                                 closest_point_hit_l.T() + 
                                 ((closest_point_hit - closest_point_wire).Mag() - sand_reco::stt::wire_radius) 
                                 / sand_reco::stt::v_drift);
    closestPoints.push_back(closest_point_hit_l);
    closestPoints.push_back(closest_point_wire_l);

  } else {
    std::cout << "Wire and hit are parallel?" << std::endl;
  }
 
  return closestPoints;

}

double GetMinWireTime(TLorentzVector point, SANDWireInfo& wire)
{
  TVector3 wire_point = wire.getReadoutPoint();

  return point.T() +
         (point.Vect() - wire_point).Mag() / sand_reco::stt::v_signal_inwire;
}

void create_digits_from_hits(const SANDGeoManager& geo,
                             std::map<SANDTrackerCellID, std::vector<hit> >& hits2cell,
                             std::vector<dg_wire>& wire_digits)
{
  wire_digits.clear();

  for (std::map<SANDTrackerCellID, std::vector<hit> >::iterator it = hits2cell.begin();
       it != hits2cell.end(); ++it)  // run over wires
  {
    long did = it->first();  // wire unique id
    auto wire_info = geo.get_cell_info(it->first())->second.wire();
    double wire_time = 999.;
    double drift_time = 999.;
    double signal_time = 999.;
    double t_hit = 999.;

    dg_wire d;
    d.det = it->second[0].det;
    d.did = did;
    d.de = 0;
    // To Do: what point do we want to save? 
    // Center or one of the attachment points?
    d.x = wire_info.center().X();
    d.y = wire_info.center().Y();
    d.z = wire_info.center().Z();
    for (unsigned int i = 0; i < it->second.size();
         i++) {  // run over hits of given wire
      auto running_hit = it->second[i];
      // find hit closest point to wire
      std::vector<TLorentzVector> ClosestPoints =
          digitization::edep_sim::tracker::WireHitClosestPoints(running_hit,
                                                                wire_info);
      if (ClosestPoints.size() == 0) {
        continue;
      }
      TLorentzVector closest_point_hit_l = ClosestPoints[0];
      // find wire closest point to hit : time of closest_point_wire_l  = drift time +
      // hit time
      TLorentzVector closest_point_wire_l = ClosestPoints[1];

      // total time = time 2 signal propagation + drift time + hit time
      double hit_smallest_time =
          digitization::edep_sim::tracker::GetMinWireTime(closest_point_wire_l,
                                                          wire_info);

      if (hit_smallest_time < wire_time) {
        wire_time = hit_smallest_time;
        t_hit = closest_point_hit_l.T();
        drift_time = closest_point_wire_l.T() - t_hit;
        signal_time = hit_smallest_time - closest_point_wire_l.T();
      }
      d.de += running_hit.de;
      d.hindex.push_back(running_hit.index);
    }
    d.tdc = wire_time + rand.Gaus(0, sand_reco::stt::tm_stt_smearing);
    d.t_hit = t_hit;
    d.drift_time = drift_time;
    d.signal_time = signal_time;
    d.adc = d.de;

    wire_digits.push_back(d);
  }
}
}

namespace ecal
{
// convert deposited energy into
// into the mean number of pe
double energy_to_photo_electrons(double E)
{
  if (debug)
    std::cout << "E = " << E
              << " -> p.e. = " << sand_reco::ecal::photo_sensor::e2pe* E
              << std::endl;

  return sand_reco::ecal::photo_sensor::e2pe * E;
}

// process calorimeter hits
// return:
//  - module id
//  - plane id
//  - cell id
//  - d1: distance from pmt1
//  - d2: distance from pmt2
//  - t: time of the hit
//  - de: energy deposit
bool process_hit(const SANDGeoManager& g, const TG4HitSegment& hit, int& detID,
                 int& modID, int& planeID, int& cellID, double& d1, double& d2,
                 double& t, double& de)
{
  if (debug) {
    std::cout << "ProcessHit" << std::endl;
  }

  detID = -999;
  modID = -999;
  planeID = -999;
  cellID = -999;
  d1 = -999;
  d2 = -999;

  double x = 0.5 * (hit.Start.X() + hit.Stop.X());
  double y = 0.5 * (hit.Start.Y() + hit.Stop.Y());
  double z = 0.5 * (hit.Start.Z() + hit.Stop.Z());

  t = 0.5 * (hit.Start.T() + hit.Stop.T());
  de = hit.EnergyDeposit;

  auto hit_direction = hit.Stop - hit.Start;

  g.SetGeoCurrentPoint(x, y, z);
  g.SetGeoCurrentDirection(hit_direction.X(), hit_direction.Y(),
                           hit_direction.Z());

  volume running_volume;
  g.InitVolume(running_volume);

  if (!running_volume.IsActive) return false;

  auto cell_global_id = g.get_ecal_cell_id(x, y, z);

  if (cell_global_id == 999 || cell_global_id == -999) return false;

  g.decode_ecal_cell_id(cell_global_id, detID, modID, planeID, cellID);

  if (cellID == -99) {
    std::cout << "\n";
    std::cout << __FILE__ << " " << __LINE__ << " \n";
    std::cout << "cell_global_id : " << cell_global_id << "\n";
    std::cout << "detID   : " << detID << "\n";
    std::cout << "modID   : " << modID << "\n";
    std::cout << "planeID : " << planeID << "\n";
    std::cout << "cellID  : " << cellID << "\n";
    throw "";
  }
  return true;

  // /////
  // TGeoNode* node = g->FindNode(x, y, z);

  // if (node == 0) return false;

  // TString str = node->GetName();
  // TString str2 = g->GetPath();

  // if (debug) {
  //   std::cout << "node name: " << str.Data() << std::endl;
  // }

  // if (sand_reco::ecal::geometry::CheckAndProcessPath(str2) == false)
  //   return false;
  // //////

  // // barrel modules
  // if (sand_reco::ecal::geometry::isBarrel(str)) {

  //   sand_reco::ecal::geometry::BarrelModuleAndLayer(str, str2, detID, modID,
  //                                                   planeID);

  //   sand_reco::ecal::geometry::BarrelCell(x, y, z, g, node, cellID, d1, d2);

  //   if (debug) {
  //     std::cout << "hit: " << str.Data() << std::endl;
  //     std::cout << "\t[x,y,z]                " << x << " " << y << " " << z
  //               << std::endl;
  //     std::cout << "\t[detID,modID,planeID,cellID] " << detID << " " << modID
  //               << " " << planeID << " " << cellID << std::endl;
  //     std::cout << "\t[d1,d2,t,de]           " << d1 << " " << d2 << " " << t
  //               << " " << de << std::endl;
  //   }

  //   return true;
  // }
  // // end cap modules
  // else if (sand_reco::ecal::geometry::isEndCap(str)) {

  //   sand_reco::ecal::geometry::EndCapModuleAndLayer(str, str2, detID, modID,
  //                                                   planeID);

  //   sand_reco::ecal::geometry::EndCapCell(x, y, z, g, node, cellID, d1, d2);

  //   if (debug) {
  //     std::cout << "hit: " << str.Data() << std::endl;
  //     std::cout << "\t[x,y,z]                " << x << " " << y << " " << z
  //               << std::endl;
  //     std::cout << "\t[detID,modID,planeID,cellID] " << detID << " " << modID
  //               << " " << planeID << " " << cellID << std::endl;
  //   }
  //   return true;
  // } else {
  //   return false;
  // }
}

void simulate_photo_electrons(TG4Event* ev, const SANDGeoManager& g,
                              std::map<int, std::vector<pe> >& photo_el,
                              std::map<int, double>& L)
{
  int detID, modID, planeID, cellID, uniqID;
  double d1, d2, t0, de;

  for (std::map<std::string, std::vector<TG4HitSegment> >::iterator it =
           ev->SegmentDetectors.begin();
       it != ev->SegmentDetectors.end(); ++it) {
    if (it->first == "EMCalSci") {
      for (unsigned int j = 0; j < it->second.size(); j++) {

        if (digitization::edep_sim::ecal::process_hit(g, it->second[j], detID,
                                                      modID, planeID, cellID,
                                                      d1, d2, t0, de) == true) {
          double en1 =
              de * sand_reco::ecal::attenuation::AttenuationFactor(d1, planeID);
          double en2 =
              de * sand_reco::ecal::attenuation::AttenuationFactor(d2, planeID);

          double ave_pe1 =
              digitization::edep_sim::ecal::energy_to_photo_electrons(en1);
          double ave_pe2 =
              digitization::edep_sim::ecal::energy_to_photo_electrons(en2);

          int pe1 = digitization::rand.Poisson(ave_pe1);
          int pe2 = digitization::rand.Poisson(ave_pe2);

          uniqID =
              sand_reco::ecal::decoder::EncodeID(detID, modID, planeID, cellID);

          // cellend 1 -> x < 0 -> ID > 0 -> left
          // cellend 2 -> x > 0 -> ID < 0 -> right

          for (int i = 0; i < pe1; i++) {
            pe this_pe;
            this_pe.time =
                digitization::ecal::photo_electron_time_to_pmt_arrival_time(t0,
                                                                            d1);
            this_pe.h_index = j;
            photo_el[uniqID].push_back(this_pe);
            L[uniqID] = d1 + d2;
          }

          for (int i = 0; i < pe2; i++) {
            pe this_pe;
            this_pe.time =
                digitization::ecal::photo_electron_time_to_pmt_arrival_time(t0,
                                                                            d2);
            this_pe.h_index = j;
            photo_el[-1 * uniqID].push_back(this_pe);
            L[-1 * uniqID] = d1 + d2;
          }
        }
      }
    }
  }
}

// construct calo digit and collect them in a vector
void group_pmts_in_cells(const SANDGeoManager& geo,
                         std::map<int, std::vector<dg_ps> >& ps,
                         std::map<int, double>& L,
                         std::vector<dg_cell>& vec_cell)
{
  std::map<int, dg_cell> map_cell;
  dg_cell* c;
  for (std::map<int, std::vector<dg_ps> >::iterator it = ps.begin();
       it != ps.end(); ++it) {
    int id = abs(it->first);
    // if(id==999) continue;
    c = &(map_cell[id]);

    c->id = id;
    sand_reco::ecal::decoder::DecodeID(c->id, c->det, c->mod, c->lay, c->cel);
    c->l = L[it->first];

    if (it->first >= 0) {
      c->ps1 = it->second;
    } else {
      c->ps2 = it->second;
    }
    if (c->id == 214294) continue;
    auto cell_info = geo.get_ecal_cell_info(c->id);
    c->x = cell_info.x();
    c->y = cell_info.y();
    c->z = cell_info.z();
  }

  for (std::map<int, dg_cell>::iterator it = map_cell.begin();
       it != map_cell.end(); ++it) {
    vec_cell.push_back(it->second);
  }
}

// simulate calorimeter responce for whole event
void digitize_ecal(TG4Event* ev, const SANDGeoManager& geo,
                   std::vector<dg_cell>& vec_cell,
                   ECAL_digi_mode ecal_digi_mode)
{
  std::map<int, std::vector<pe> > photo_el;
  std::map<int, std::vector<dg_ps> > ps;
  std::map<int, double> L;

  vec_cell.clear();

  if (debug) {
    std::cout << "SimulatePE" << std::endl;
  }

  digitization::edep_sim::ecal::simulate_photo_electrons(ev, geo, photo_el, L);
  if (debug) {
    std::cout << "TimeAndSignal" << std::endl;
  }
  digitization::ecal::eval_adc_and_tdc_from_photo_electrons(photo_el, ps,
                                                            ecal_digi_mode);
  if (debug) {
    std::cout << "CollectSignal" << std::endl;
  }
  digitization::edep_sim::ecal::group_pmts_in_cells(geo, ps, L, vec_cell);
}

}  // namespace ecal

namespace stt
{
// Group hits into tube
void group_hits_by_tube(TG4Event* ev, const SANDGeoManager& geo,
                        std::map<SANDTrackerCellID, std::vector<hit> >& hits2Tube)
{
  hits2Tube.clear();

  int skipped_hit = 0;
  int all_hit = ev->SegmentDetectors["Straw"].size();

  for (unsigned int j = 0; j < ev->SegmentDetectors["Straw"].size(); j++) {
    const TG4HitSegment& hseg = ev->SegmentDetectors["Straw"].at(j);

    double x = 0.5 * (hseg.Start.X() + hseg.Stop.X());
    double y = 0.5 * (hseg.Start.Y() + hseg.Stop.Y());
    double z = 0.5 * (hseg.Start.Z() + hseg.Stop.Z());

    SANDTrackerCellID stid = geo.get_stt_tube_id(x, y, z);

    if (stid == -999) {
      // std::cout << std::setprecision(12) << x << " " << y << " " << z <<
      // std::endl;
      // geo.print_stt_tube_id(x,y,z);
      skipped_hit++;
      continue;
    };

    // std::string sttname = "NULL";
    // int stid = -999;  // should be implemented for FLUKA

    // sttname = geo->FindNode(x, y, z)->GetName();

    // stid = sand_reco::stt::getSTUniqID(geo, x, y, z);
    // if (stid == -999) continue;

    hit h;
    h.det = "Straw";
    h.did = stid();
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

    hits2Tube[stid].push_back(h);
  }
  if (skipped_hit != 0) {
    std::cout << "WARNING: " << skipped_hit << " out of " << all_hit
              << " hits skipped due to unexpected volume path!!" << std::endl;
  }
  if (skipped_hit != 0) {
    std::cout << "WARNING: " << skipped_hit << " out of " << all_hit
              << " hits skipped due to unexpected volume path!!" << std::endl;
  }
}

// for each tube simulate tdc and adc
// tdc is the time of closest point to wire + drift time
// adc is the sum of energy deposit within integration time window
void create_digits_from_hits(const SANDGeoManager& geo,
                             std::map<SANDTrackerCellID, std::vector<hit> >& hits2Tube,
                             std::vector<dg_wire>& wire_digits)
{
  wire_digits.clear();

  for (std::map<SANDTrackerCellID, std::vector<hit> >::iterator it = hits2Tube.begin();
       it != hits2Tube.end(); ++it) {
    double min_time_tub   = 1E9;  // mm
    double min_drift_time = 1E9;  // mm
    SANDTrackerCellID cell_global_id = it->first;

    auto cell_info = geo.get_cell_info(cell_global_id)->second;

    SANDTrackerModuleID module_unique_id;
    SANDTrackerPlaneID plane_global_id, plane_local_id, plane_type;
    SANDTrackerCellID cell_local_id;

    SANDGeoManager::decode_cell_id(cell_global_id, plane_global_id, cell_local_id);
    SANDGeoManager::decode_plane_id(plane_global_id, module_unique_id, 
                                    plane_local_id, plane_type);

    dg_wire d;
    d.det = it->second[0].det;
    d.did =  cell_global_id();
    d.de = 0;
    d.hor = (plane_type() % 2 == 0);
    d.t0 = sand_reco::t0[plane_global_id()];
    TVector2 wire;
    if (d.hor == true) {
      d.x = sand_reco::stt::stt_center[0];
      d.y = cell_info.wire().y();
      d.z = cell_info.wire().z();
      wire.SetX(cell_info.wire().z());
      wire.SetY(cell_info.wire().y());
    } else {
      d.x = cell_info.wire().x();
      d.y = sand_reco::stt::stt_center[1];
      d.z = cell_info.wire().z();
      wire.SetX(cell_info.wire().z());
      wire.SetY(cell_info.wire().x());
    }

    for (unsigned int i = 0; i < it->second.size(); i++) {
      double x1 = it->second[i].z1;
      double x2 = it->second[i].z2;
      double t1 = it->second[i].t1;
      double t2 = it->second[i].t2;

      double y1, y2;
      double z1, z2, z;
      double l, dwire;

      if (plane_type == 2) {
        y1 = it->second[i].y1;
        y2 = it->second[i].y2;
        z1 = it->second[i].x1;
        z2 = it->second[i].x2;
        l = sand_reco::stt::getT(y1, y2, cell_info.wire().y(), x1, x2, cell_info.wire().z());
        z = z1 + (z2 - z1) * l;
        dwire = cell_info.wire().x() + cell_info.wire().length() - z;
      } else {
        y1 = it->second[i].x1;
        y2 = it->second[i].x2;
        z1 = it->second[i].y1;
        z2 = it->second[i].y2;
        l = sand_reco::stt::getT(y1, y2, cell_info.wire().x(), x1, x2, cell_info.wire().z());
        z = z1 + (z2 - z1) * l;
        dwire = cell_info.wire().y() + cell_info.wire().length() - z;
      }

      double x = x1 + (x2 - x1) * l;
      double y = y1 + (y2 - y1) * l;
      double t = t1 + (t2 - t1) * l;

      TVector2 min_dist_point(x, y);
      double min_dist_hit = (min_dist_point - wire).Mod();
      double min_drift_hit = (min_dist_hit - sand_reco::stt::wire_radius) /
                                    sand_reco::stt::v_drift;
      double min_time_hit = t + min_drift_hit +
                            dwire / sand_reco::stt::v_signal_inwire;

      if (min_time_hit < min_time_tub) {
        min_time_tub = min_time_hit;
        min_drift_time = min_drift_hit;
      }

      if (t - d.t0 < sand_reco::stt::stt_int_time) d.de += it->second[i].de;

      d.hindex.push_back(it->second[i].index);
    }

    d.tdc = min_time_tub + rand.Gaus(0, sand_reco::stt::tm_stt_smearing);
    d.drift_time = min_drift_time;
    d.adc = d.de;

    wire_digits.push_back(d);
  }
}

// simulate stt responce for whole event
void digitize_stt(TG4Event* ev, const SANDGeoManager& geo,
                  std::vector<dg_wire>& wire_digits)
{
  std::map<SANDTrackerCellID, std::vector<hit> > hits2Tube;
  wire_digits.clear();

  group_hits_by_tube(ev, geo, hits2Tube);
  digitization::edep_sim::tracker::create_digits_from_hits(geo, hits2Tube,
                                                       wire_digits);
}
}  // namespace stt

namespace chamber
{

TVector3 IntersectHitPlane(const TG4HitSegment& hseg, double plane_coordinate,
                           SANDWireInfo::Orient plane_orientation)
{
  // find intersect between hit segment and wire plane
  auto delta = hseg.Stop - hseg.Start;
  double t = -999.;

  TVector3 crossing_point = {-999., -999., -999.};

  if (plane_orientation == SANDWireInfo::Orient::kHorizontal) {
    t = (plane_coordinate - hseg.Start.Y()) / delta.Y();
    crossing_point.SetXYZ(hseg.Start.X() + delta.X() * t, plane_coordinate,
                          hseg.Start.Z() + delta.Z() * t);

  } else {
    t = (plane_coordinate - hseg.Start.X()) / delta.X();
    crossing_point.SetXYZ(plane_coordinate, hseg.Start.Y() + delta.Y() * t,
                          hseg.Start.Z() + delta.Z() * t);
  }
  return crossing_point;
}

void group_hits_by_cell(TG4Event* ev, const SANDGeoManager& geo,
                        std::map<SANDTrackerCellID, std::vector<hit> >& hits2cell)
{
  hits2cell.clear();

  for (unsigned int j = 0; j < ev->SegmentDetectors["DriftVolume"].size(); j++) {
    const TG4HitSegment& hseg = ev->SegmentDetectors["DriftVolume"].at(j);

    int pdg = ev->Trajectories[hseg.GetPrimaryId()].GetPDGCode();
    std::vector<SANDTrackerCellID> ids = geo.get_segment_ids(hseg);
    SANDTrackerCellID id1 = ids[0];
    SANDTrackerCellID id2 = ids[1];

    if (id1 == -999) {
      // std::cout<<"skipping this hit\n";
      continue;
    }

    SANDTrackerPlaneID plane_global_id1;
    SANDTrackerPlaneID plane_global_id2;
    SANDTrackerCellID cell_local_id1;
    SANDTrackerCellID cell_local_id2;
    geo.decode_cell_id(id1, plane_global_id1, cell_local_id1);
    geo.decode_cell_id(id2, plane_global_id2, cell_local_id2);

    // std::cout << id1() << " "  << id2() << " " << std::endl;
    // std::cout << plane_global_id1() << " " << plane_global_id2() << std::endl;
    // std::cout << cell_local_id1() << " " << cell_local_id2() << std::endl;

    if (plane_global_id1 != plane_global_id2) {
      std::cout << "WIRE ID CORRESPONDING TO 2 DIFFERENT DIRFT PLANES"
                << std::endl;
      break;
    }

    long start_id = 999;
    long stop_id = 999;
    if (id2 > id1) {
      start_id = id1();
      stop_id = id2();
    } else if (id2 < id1) {
      start_id = id2();
      stop_id = id1();
    } else  // hit in 1 cell
    {
      hit h;
      h.det = "DriftVolume";
      h.did = id1();
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
      hits2cell[id1].push_back(h);
      continue;
    }

    TVector3 start = {hseg.Start.X(), hseg.Start.Y(), hseg.Start.Z()};
    double hseg_length = (hseg.Stop - hseg.Start).Mag();
    double hseg_dt = (hseg.Stop - hseg.Start).T();
    double hseg_start_t = hseg.Start.T();

    auto& plane = *geo.get_plane_info(SANDTrackerCellID(start_id));

    TVector2 rotated_hit_start_2d_position = geo.GlobalToRotated(TVector2(hseg.Start.X(), hseg.Start.Y()), plane);
    TVector2 rotated_hit_stop_2d_position  = geo.GlobalToRotated(TVector2(hseg.Stop.X(), hseg.Stop.Y())  , plane);

    auto rotated_delta_x = rotated_hit_stop_2d_position.X() - rotated_hit_start_2d_position.X();
    auto rotated_delta_y = rotated_hit_stop_2d_position.Y() - rotated_hit_start_2d_position.Y();
    auto rotated_delta_z = hseg.Stop.Z() - hseg.Start.Z();

    for (auto i = start_id; i <= stop_id; i++) {
      auto cell1 = geo.get_cell_info(i);
      auto cell2 = geo.get_cell_info(i + 1);

      
      TVector2 rotated_start_2d_position = geo.GlobalToRotated(TVector2(start.X(), start.Y()), plane);
      double transverse_coord_start = rotated_start_2d_position.Y();

      double step_coordinate;

      if (cell2 != plane.getIdToCellMapEnd()) {

        SANDWireInfo wire1 = cell1->second.wire();
        SANDWireInfo wire2 = cell2->second.wire();
        
        TVector2 rotated_wire_center1_2d_position = geo.GlobalToRotated(TVector2(wire1.center().X(), wire1.center().Y()), plane);
        double transverse_coord1 = rotated_wire_center1_2d_position.Y();
 
        TVector2 rotated_wire_center2_2d_position = geo.GlobalToRotated(TVector2(wire2.center().X(), wire2.center().Y()), plane);
        double transverse_coord2 = rotated_wire_center2_2d_position.Y();

        double plane_coordinate = (transverse_coord1 + transverse_coord2) * 0.5;
        
        if (fabs(plane_coordinate - transverse_coord_start) < 
            fabs(rotated_hit_stop_2d_position.Y() - transverse_coord_start)) {
          step_coordinate = plane_coordinate;
        } else {
          step_coordinate = rotated_hit_stop_2d_position.Y();
        }
      } else {
        step_coordinate = rotated_hit_stop_2d_position.Y();
      }
      double t = fabs((step_coordinate - transverse_coord_start) / rotated_delta_y);
      
      TVector2 rotated_crossing_point(rotated_start_2d_position.X() + rotated_delta_x * t, 
                                      rotated_start_2d_position.Y() + rotated_delta_y * t);

      TVector2 global_crossing_point = geo.RotatedToGlobal(TVector2(rotated_crossing_point.X(), rotated_crossing_point.Y()), plane);

      TVector3 stop(global_crossing_point.X(), 
                    global_crossing_point.Y(), 
                    start.Z() + rotated_delta_z * t);

      double portion = (start - stop).Mag() / hseg_length;

      TVector3 center = (start + stop) * 0.5;

      // start.Print();
      // start_2d_position.Print();
      // crossing_point.Print();
      // crossing_point_2d.Print();
      // stop.Print();
      // hseg.Stop.Print();
      SANDTrackerCellID cell_id = geo.GetClosestCellToHit(center, plane, false);

      hit h;
      h.det = "DriftVolume";
      h.did = cell_id();
      h.x1 = start.X();
      h.y1 = start.Y();
      h.z1 = start.Z();
      h.t1 = hseg_start_t;
      h.x2 = stop.X();
      h.y2 = stop.Y();
      h.z2 = stop.Z();
      h.t2 = hseg_start_t + t * hseg_dt;
      h.de = hseg.EnergyDeposit * t;
      h.pid = hseg.PrimaryId;
      h.index = j;

      hits2cell[cell_id].push_back(h);

      start = stop;
      hseg_start_t += t * hseg_dt;

      if ((start - hseg.Stop.Vect()).Mag() < 1E-6) {
        break;
      }

    }
  }
}

bool isInWire(SANDWireInfo& wire, TVector3& point)
{
  TVector3 wire3 = {wire.center().x(), wire.center().y(), wire.center().z()};
  return ((wire3 - point).Mag() <= wire.length() * 0.5);
}

bool isInHit(hit& h, TVector3& point)
{
  TVector3 middle = {(h.x1 + h.x2) / 2., (h.y1 + h.y2) / 2.,
                     (h.z1 + h.z2) / 2.};
  TVector3 start = {h.x1, h.y1, h.z1};
  return ((point - middle).Mag() <= (start - middle).Mag());
}



// simulate wire responce for whole event
void digitize_drift(TG4Event* ev, const SANDGeoManager& geo,
                    std::vector<dg_wire>& wire_digits)
{
  std::map<SANDTrackerCellID, std::vector<hit> > hits2cell;
  wire_digits.clear();

  group_hits_by_cell(ev, geo, hits2cell);
  digitization::edep_sim::tracker::create_digits_from_hits(geo, hits2cell,
                                                           wire_digits);
}

}  // namespace chamber

// digitize event
void digitize(const char* finname, const char* foutname,
              ECAL_digi_mode ecal_digi_mode)
{
  TFile f(finname, "READ");

  std::cout << "This is a standard Geant4-edepsim SIMULATION" << std::endl;

  // MC info tree
  TTree* t = (TTree*)f.Get("EDepSimEvents");

  // Event
  TG4Event* ev = new TG4Event;
  t->SetBranchAddress("Event", &ev);

  // TGeoManager pointer used in case of geant4-based simulation
  TGeoManager* geo = 0;

  // Get TGeoManager or additional Tree depending on the simulation chain
  geo = (TGeoManager*)f.Get("EDepSimGeometry");

  if (debug) std::cout << "Inizializzo la geometria" << std::endl;

  // Initialization of detector-geometry-related
  // usefull variables defined in utils.h
  // Specifically:
  // - calo cell center and size
  //   - double czlay[nLay]: "radial" position of the cells
  //   - double cxlay[nLay][nCel]: transverse position of the cells
  //   - double ec_r: ECAL endcap radius
  //   - double ec_dz: ECAL endcap thickness
  // - stt center and length
  //   - double stt_center[3]: SAND internal volume center
  //   - std::map<int, double> sand_reco::stL
  //       - KEY  : id of the tube
  //       - VALUE: length of the tubes
  //   - std::map<int, std::map<double, int>> sand_reco::stX
  //       - KEY  : id of the plane
  //       - VALUE: map with
  //         - KEY  : beam transversal coordinate of the tubes
  //         - VALUE: id of the tube
  //   - std::map<int, std::map<int, TVector2>> sand_reco::stPos
  //       - KEY  : id of the plane
  //       - VALUE: map with
  //         - KEY  : local id of the tube
  //         - VALUE: position of the center of the tube
  //                  (X in TVector2 is Z in the global system)
  //                  (Y in TVector2 is X or Y in the global system)
  //   - std::map<int, TVector2> sand_reco::tubePos
  //       - KEY  : id of the tube
  //       - VALUE: map with
  // sand_reco::init(geo);
  SANDGeoManager sand_geo;
  sand_geo.init(geo);

  // vector of ECAL and STT digits
  std::vector<dg_cell> vec_cell;
  std::vector<dg_wire> wire_digits;

  // output
  TFile fout(foutname, "RECREATE");
  TTree tout("tDigit", "Digitization");

  tout.Branch("dg_cell", "std::vector<dg_cell>", &vec_cell);

  if (geo->FindVolumeFast("STTtracker_PV")) {
    std::cout << "\n--- Digitize STT based simulation ---\n";
  } else {
    std::cout << "\n--- Digitize Drift based simulation ---\n";
  }
  tout.Branch("dg_wire", "std::vector<dg_wire>", &wire_digits);

  // number of events
  const int nev = t->GetEntries();

  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;

  // loop on all input events
  for (int i = 0; i < nev; i++) {
    t->GetEntry(i);

    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i) / nev * 100)
              << "%]" << std::flush;

    // define the T0 for this event
    // for each straw tubs:
    // std::map<int, double> sand_reco::t0
    sand_reco::stt::initT0(ev, sand_geo);
    digitization::edep_sim::ecal::digitize_ecal(ev, sand_geo, vec_cell,
                                                ecal_digi_mode);

    if (geo->FindVolumeFast("STTtracker_PV")) {
      digitization::edep_sim::stt::digitize_stt(ev, sand_geo, wire_digits);
    } else {
      digitization::edep_sim::chamber::digitize_drift(ev, sand_geo,
                                                      wire_digits);
    }

    tout.Fill();
  }
  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;

  sand_geo.PrintCounter();

  // write output
  fout.cd();
  tout.Write();
  geo->Write();
  fout.Close();

  f.Close();

  // cleaning
  // sand_reco::stt::stL.clear();
  // sand_reco::stt::stX.clear();
  // sand_reco::stt::stPos.clear();
  // sand_reco::stt::tubePos.clear();
  sand_reco::t0.clear();
}

}  // namespace edep_sim

}  // namespace digitization