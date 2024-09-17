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
                        std::map<long, std::vector<hit> >& hits2Tube)
{
  hits2Tube.clear();

  int skipped_hit = 0;
  int all_hit = ev->SegmentDetectors["Straw"].size();

  for (unsigned int j = 0; j < ev->SegmentDetectors["Straw"].size(); j++) {
    const TG4HitSegment& hseg = ev->SegmentDetectors["Straw"].at(j);

    double x = 0.5 * (hseg.Start.X() + hseg.Stop.X());
    double y = 0.5 * (hseg.Start.Y() + hseg.Stop.Y());
    double z = 0.5 * (hseg.Start.Z() + hseg.Stop.Z());

    long stid = geo.get_stt_tube_id(x, y, z);

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
    h.det = "STT";
    h.did = stid;
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
                             std::map<long, std::vector<hit> >& hits2Tube,
                             std::vector<dg_wire>& wire_digits)
{
  wire_digits.clear();

  for (std::map<long, std::vector<hit> >::iterator it = hits2Tube.begin();
       it != hits2Tube.end(); ++it) {
    double min_time_tub = 1E9;  // mm
    long did = it->first;

    long supmod, mod, tub, type, pla, plloc;

    SANDGeoManager::decode_wire_id(did, pla, tub);
    SANDGeoManager::decode_plane_id(supmod, pla, mod, plloc, type);

    auto stt_info = geo.get_wire_info(did);

    dg_wire d;
    d.det = it->second[0].det;
    d.did = did;
    d.de = 0;
    d.hor = (type % 2 == 0);
    d.t0 = sand_reco::t0[pla];
    TVector2 wire;
    if (d.hor == true) {
      d.x = sand_reco::stt::stt_center[0];
      d.y = stt_info.y();
      d.z = stt_info.z();
      wire.SetX(stt_info.z());
      wire.SetY(stt_info.y());
    } else {
      d.x = stt_info.x();
      d.y = sand_reco::stt::stt_center[1];
      d.z = stt_info.z();
      wire.SetX(stt_info.z());
      wire.SetY(stt_info.x());
    }

    for (unsigned int i = 0; i < it->second.size(); i++) {
      double x1 = it->second[i].z1;
      double x2 = it->second[i].z2;
      double t1 = it->second[i].t1;
      double t2 = it->second[i].t2;

      double y1, y2;
      double z1, z2, z;
      double l, dwire;

      if (type == 2) {
        y1 = it->second[i].y1;
        y2 = it->second[i].y2;
        z1 = it->second[i].x1;
        z2 = it->second[i].x2;
        l = sand_reco::stt::getT(y1, y2, stt_info.y(), x1, x2, stt_info.z());
        z = z1 + (z2 - z1) * l;
        dwire = stt_info.x() + stt_info.length() - z;
      } else {
        y1 = it->second[i].x1;
        y2 = it->second[i].x2;
        z1 = it->second[i].y1;
        z2 = it->second[i].y2;
        l = sand_reco::stt::getT(y1, y2, stt_info.x(), x1, x2, stt_info.z());
        z = z1 + (z2 - z1) * l;
        dwire = stt_info.y() + stt_info.length() - z;
      }

      double x = x1 + (x2 - x1) * l;
      double y = y1 + (y2 - y1) * l;
      double t = t1 + (t2 - t1) * l;

      TVector2 min_dist_point(x, y);
      double min_dist_hit = (min_dist_point - wire).Mod();
      double min_time_hit = t + (min_dist_hit - sand_reco::stt::wire_radius) /
                                    sand_reco::stt::v_drift +
                            dwire / sand_reco::stt::v_signal_inwire;

      if (min_time_hit < min_time_tub) min_time_tub = min_time_hit;

      if (t - d.t0 < sand_reco::stt::stt_int_time) d.de += it->second[i].de;

      d.hindex.push_back(it->second[i].index);
    }

    d.tdc = min_time_tub + rand.Gaus(0, sand_reco::stt::tm_stt_smearing);
    d.adc = d.de;

    wire_digits.push_back(d);
  }
}

// simulate stt responce for whole event
void digitize_stt(TG4Event* ev, const SANDGeoManager& geo,
                  std::vector<dg_wire>& wire_digits)
{
  std::map<long, std::vector<hit> > hits2Tube;
  wire_digits.clear();

  group_hits_by_tube(ev, geo, hits2Tube);
  digitization::edep_sim::stt::create_digits_from_hits(geo, hits2Tube,
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
                        std::map<long, std::vector<hit> >& hits2cell)
{
  hits2cell.clear();

  for (unsigned int j = 0; j < ev->SegmentDetectors["DriftVolume"].size();
       j++) {
    const TG4HitSegment& hseg = ev->SegmentDetectors["DriftVolume"].at(j);

    int pdg = ev->Trajectories[hseg.GetPrimaryId()].GetPDGCode();
    // std::cout<<"flag,"<<pdg<<","<<hseg.Start.X()<<","<<hseg.Start.Y()<<","<<hseg.Start.Z()<<","
    //                                 <<hseg.Stop.X()<<","<<hseg.Stop.Y()<<","<<hseg.Stop.Z()<<"\n";

    // if(pdg!=13) continue; //->test digit only for muons

    std::vector<long> ids = geo.get_segment_ids(hseg);
    long id1 = ids[0];
    long id2 = ids[1];

    if (id1 == -999) {
      // std::cout<<"skipping this hit\n";
      continue;
    }

    long plane_global_id1;
    long plane_global_id2;
    long wire_local_id1;
    long wire_local_id2;
    geo.decode_wire_id(id1, plane_global_id1, wire_local_id1);
    geo.decode_wire_id(id2, plane_global_id2, wire_local_id2);

    // std::cout << id1 << " "  << id2 << " " << std::endl;
    // std::cout << plane_global_id1 << " " << plane_global_id2 << std::endl;
    // std::cout << wire_local_id1 << " " << wire_local_id2 << std::endl;

    if (plane_global_id1 != plane_global_id2) {
      std::cout << "WIRE ID CORRESPONDING TO 2 DIFFERENT DIRFT PLANES"
                << std::endl;
      break;
    }

    long start_id = 999;
    long stop_id = 999;
    if (id2 > id1) {
      start_id = id1;
      stop_id = id2;
    } else if (id2 < id1) {
      start_id = id2;
      stop_id = id1;
    } else  // hit in 1 cell
    {
      hit h;
      h.det = "DriftChamber";
      h.did = id1;
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

    for (auto i = start_id; i < stop_id; i++) {
      SANDWireInfo wire1 = geo.get_wire_info(i);
      SANDWireInfo wire2 = geo.get_wire_info(i + 1);

      double plane_coordinate =
          (wire1.orientation() == SANDWireInfo::Orient::kHorizontal)
              ? (wire2.y() + wire1.y()) / 2.
              : (wire2.x() + wire1.x()) / 2.;

      TVector3 stop = digitization::edep_sim::chamber::IntersectHitPlane(
          hseg, plane_coordinate, wire1.orientation());

      double portion = (start - stop).Mag() / hseg_length;

      hit h;
      h.det = "DriftChamber";
      h.did = i;
      h.x1 = start.X();
      h.y1 = start.Y();
      h.z1 = start.Z();
      h.t1 = hseg_start_t;
      h.x2 = stop.X();
      h.y2 = stop.Y();
      h.z2 = stop.Z();
      h.t2 = hseg_start_t + portion * hseg_dt;
      h.de = hseg.EnergyDeposit * portion;
      h.pid = hseg.PrimaryId;
      h.index = j;

      start = stop;
      hseg_start_t += portion * hseg_dt;

      hits2cell[i].push_back(h);
    }
  }
}

bool isInWire(SANDWireInfo& wire, TVector3& point)
{
  TVector3 wire3 = {wire.x(), wire.y(), wire.z()};
  return ((wire3 - point).Mag() <= wire.length() * 0.5);
}

bool isInHit(hit& h, TVector3& point)
{
  TVector3 middle = {(h.x1 + h.x2) / 2., (h.y1 + h.y2) / 2.,
                     (h.z1 + h.z2) / 2.};
  TVector3 start = {h.x1, h.y1, h.z1};
  return ((point - middle).Mag() <= (start - middle).Mag());
}

std::vector<TLorentzVector> WireHitClosestPoints(hit& h, SANDWireInfo& arg_wire)
{
  // return hit closest point to wire between hit start and stop and viceversa
  // -----------------------------------------------------------
  TVector3 start = {h.x1, h.y1, h.z1};  // hit start
  TVector3 stop = {h.x2, h.y2, h.z2};   // hit end
  TVector3 s = stop - start;            // hit direction

  TVector3 wire = {arg_wire.x(), arg_wire.y(), arg_wire.z()};
  TVector3 r0 = {-999., -999., -999.};  // point on the wire
  TVector3 r = {-999., -999., -999.};   // wire direction
  double angle = 0.;                    // wire.angle();
  double m = -999.;                     // angular coeffcient
  TVector3 leftend = {arg_wire.x() - arg_wire.length() * std::cos(angle),
                      arg_wire.y() - arg_wire.length() * std::sin(angle),
                      arg_wire.z()};  // wire left/bottom end

  TVector3 rightend = {arg_wire.x() + arg_wire.length() * std::cos(angle),
                       arg_wire.y() + arg_wire.length() * std::sin(angle),
                       arg_wire.z()};  // wire right/top

  double t = -999.;  // real number to find the closest hit point to the wire
  double t_prime =
      -999.;  // real number to find the closest wire point to the hit
  // 1 check wire orientation: horizontal, vertical or inclined
  if (arg_wire.orientation() == SANDWireInfo::Orient::kVertical) {
    m = 0.;
    r = {0., 1., 0.};
    r0 = {wire.x(), 0., wire.z()};
  } else {
    m = std::tan(angle);
    r = {1., m, 0.};
    r0 = {0., wire.y() - m * wire.x(), wire.z()};
  }
  // 2 find closest point
  t = (r0.Dot(s) + start.Dot(r) * s.Dot(r) / (r.Mag() * r.Mag()) -
       start.Dot(s) - r0.Dot(r) * r.Dot(s) / (r.Mag() * r.Mag())) /
      ((s.Mag() * s.Mag()) - (r.Dot(s) * r.Dot(s)) / (r.Mag() * r.Mag()));

  t_prime = (start.Dot(r) + t * s.Dot(r) - r0.Dot(r)) / (r.Mag() * r.Mag());
  // check if it is between start and stop, otherwise set the closest to be
  // either start or stop
  TVector3 v_closest2Wire = {start.X() + t * s.X(), start.Y() + t * s.Y(),
                             start.Z() + t * s.Z()};

  TVector3 v_closest2Hit = {r0.X() + t_prime * r.X(), r0.Y() + t_prime * r.Y(),
                            r0.Z() + t_prime * r.Z()};

  TLorentzVector closest2Wire, closest2Hit;

  if (digitization::edep_sim::chamber::isInWire(arg_wire, v_closest2Hit)) {
    closest2Hit.SetXYZT(v_closest2Hit.X(), v_closest2Hit.Y(), v_closest2Hit.Z(),
                        0.);
  } else {
    // closest is one of the 2 extremes
    if ((leftend - v_closest2Wire).Mag() < (rightend - v_closest2Wire).Mag()) {
      closest2Hit.SetXYZT(leftend.X(), leftend.Y(), leftend.Z(), 0.);
    } else {
      closest2Hit.SetXYZT(rightend.X(), rightend.Y(), rightend.Z(), 0.);
    }
  }

  if (digitization::edep_sim::chamber::isInHit(h, v_closest2Wire)) {
    double fraction = (v_closest2Wire - start).Mag() / s.Mag();
    closest2Wire.SetXYZT(v_closest2Wire.X(), v_closest2Wire.Y(),
                         v_closest2Wire.Z(), h.t1 + fraction * h.t1);
  } else {
    // closest is either start or stop
    if ((start - v_closest2Hit).Mag() < (stop - v_closest2Hit).Mag()) {
      closest2Wire.SetXYZT(h.x1, h.y1, h.z1, h.t1);
    } else {
      closest2Wire.SetXYZT(h.x2, h.y2, h.z2, h.t2);
    }
  }

  closest2Hit.SetT(closest2Wire.T() + ((v_closest2Wire - v_closest2Hit).Mag() -
                                       sand_reco::stt::wire_radius) /
                                          sand_reco::stt::v_drift);

  std::vector<TLorentzVector> closestPoints = {closest2Wire, closest2Hit};

  return closestPoints;
}

double GetMinWireTime(TLorentzVector point, SANDWireInfo& arg_wire)
{
  TVector3 wire = {arg_wire.x(), arg_wire.y(), arg_wire.z()};
  double angle = 0.;  // wire.angle();
  double m = -999.;   // angular coeffcient
  TVector3 leftend = {arg_wire.x() - arg_wire.length() * std::cos(angle),
                      arg_wire.y() - arg_wire.length() * std::sin(angle),
                      arg_wire.z()};  // wire left/bottom end

  return point.T() +
         (point.Vect() - leftend).Mag() / sand_reco::stt::v_signal_inwire;
}

void create_digits_from_hits(const SANDGeoManager& geo,
                             std::map<long, std::vector<hit> >& hits2cell,
                             std::vector<dg_wire>& wire_digits)
{
  wire_digits.clear();

  for (std::map<long, std::vector<hit> >::iterator it = hits2cell.begin();
       it != hits2cell.end(); ++it)  // run over wires
  {
    long did = it->first;  // wire unique id
    auto wire_info = geo.get_wire_info(did);
    double wire_time = 999.;
    double drift_time = 999.;
    double signal_time = 999.;
    double t_hit = 999.;

    dg_wire d;
    d.det = it->second[0].det;
    d.did = did;
    d.de = 0;
    d.hor = (wire_info.orientation() == SANDWireInfo::Orient::kHorizontal);
    d.x = wire_info.x();
    d.y = wire_info.y();
    d.z = wire_info.z();
    for (unsigned int i = 0; i < it->second.size();
         i++) {  // run over hits of given wire
      auto running_hit = it->second[i];
      // find hit closest point to wire
      std::vector<TLorentzVector> ClosestPoints =
          digitization::edep_sim::chamber::WireHitClosestPoints(running_hit,
                                                                wire_info);

      TLorentzVector closest2wire = ClosestPoints[0];
      // find wire closest point to hit : time of closest2hit  = drift time +
      // hit time
      TLorentzVector closest2hit = ClosestPoints[1];

      // total time = time 2 signal propagation + drift time + hit time
      double hit_smallest_time =
          digitization::edep_sim::chamber::GetMinWireTime(closest2hit,
                                                          wire_info);

      if (hit_smallest_time < wire_time) {
        wire_time = hit_smallest_time;
        t_hit = closest2wire.T();
        drift_time = closest2hit.T() - t_hit;
        signal_time = hit_smallest_time - closest2hit.T();
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

// simulate wire responce for whole event
void digitize_drift(TG4Event* ev, const SANDGeoManager& geo,
                    std::vector<dg_wire>& wire_digits)
{
  std::map<long, std::vector<hit> > hits2cell;
  wire_digits.clear();

  group_hits_by_cell(ev, geo, hits2cell);
  digitization::edep_sim::chamber::create_digits_from_hits(geo, hits2cell,
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