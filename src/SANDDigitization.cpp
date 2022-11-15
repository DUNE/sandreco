#include "SANDDigitization.h"

#include <iostream>

using namespace sand_reco;

namespace digitization
{

TRandom3 rand(0);

namespace ecal
{
// simulate pe arrival time to pmt
double photo_electron_time_to_pmt_arrival_time(double t0, double d)
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

  double tdec = sand_reco::ecal::scintillation::tscin *
                TMath::Power(1. / rand.Uniform() - 1.,
                             sand_reco::ecal::scintillation::tscex);

  double time = t0 + tdec +
                sand_reco::ecal::scintillation::vlfb * d * conversion::mm_to_m +
                rand.Gaus();

  if (debug) {
    std::cout << "time : " << time << std::endl;
    std::cout << "t0   : " << t0 << std::endl;
    std::cout << "scint: " << tdec << std::endl;
    std::cout << "prop : "
              << sand_reco::ecal::scintillation::vlfb * d * conversion::mm_to_m
              << std::endl;
  }

  return time;
}

// from simulated pe produce adc e tdc of calo cell
void eval_adc_and_tdc_from_photo_electrons(
    std::map<int, std::vector<pe> >& photo_el,
    std::map<int, std::vector<dg_ps> >& map_pmt)
{
  /*
    -  ADC - Proportional to NPHE
    -  TDC - Constant fraction - simulated
             TPHE(1...NPHE) in increasing time order
             IND_SEL= 0.15*NPHE
             TDC_cell = TPHE(IND_SEL)
  */

  // https://www-sciencedirect-com.ezproxy.cern.ch/science/article/pii/S0168900297013491

  double int_start;
  int pe_count;
  int start_index;
  int index;

  std::vector<pe> photo_el_digit;

  for (std::map<int, std::vector<pe> >::iterator it = photo_el.begin();
       it != photo_el.end(); ++it) {
    // order by arrival time
    std::sort(it->second.begin(), it->second.end(),
              sand_reco::ecal::isPeBefore);

    auto side = 2 * (it->first > 0) - 1;

    photo_el_digit.clear();

    int_start = it->second.front().time;
    pe_count = 0;
    start_index = 0;
    index = 0;

    for (std::vector<pe>::iterator this_pe = it->second.begin();
         this_pe != it->second.end(); ++this_pe) {
      // integrate for int_time
      if (this_pe->time < int_start + sand_reco::ecal::acquisition::int_time) {
        pe_count++;
        photo_el_digit.push_back(*this_pe);
      } else if (this_pe->time > int_start +
                                     sand_reco::ecal::acquisition::int_time +
                                     sand_reco::ecal::acquisition::dead_time) {
        // above threshold -> digit
        if (pe_count > sand_reco::ecal::acquisition::pe_threshold) {
          dg_ps signal;
          signal.side = side;
          signal.adc = sand_reco::ecal::acquisition::pe2ADC * pe_count;
          index =
              int(sand_reco::ecal::acquisition::costant_fraction * pe_count) +
              start_index;
          signal.tdc = it->second[index].time;
          signal.photo_el = photo_el_digit;
          map_pmt[it->first].push_back(signal);
        }
        // get ready for next digiit
        pe_count = 1;
        photo_el_digit.clear();
        int_start = this_pe->time;
        start_index = this_pe - it->second.begin();
      }
    }

    if (pe_count > sand_reco::ecal::acquisition::pe_threshold) {
      dg_ps signal;
      signal.side = side;
      signal.adc = sand_reco::ecal::acquisition::pe2ADC * pe_count;
      index = int(sand_reco::ecal::acquisition::costant_fraction * pe_count) +
              start_index;
      signal.tdc = it->second[index].time;
      signal.photo_el = photo_el_digit;
      map_pmt[it->first].push_back(signal);
    }
  }
}

// construct calo digit and collect them in a vector
void group_pmts_in_cells(TGeoManager* geo,
                         std::map<int, std::vector<dg_ps> >& ps,
                         std::map<int, double>& L,
                         std::vector<dg_cell>& vec_cell)
{
  std::map<int, dg_cell> map_cell;
  dg_cell* c;

  for (std::map<int, std::vector<dg_ps> >::iterator it = ps.begin();
       it != ps.end(); ++it) {
    int id = abs(it->first);

    c = &(map_cell[id]);

    c->id = id;
    sand_reco::ecal::decoder::DecodeID(c->id, c->det, c->mod, c->lay, c->cel);
    c->l = L[it->first];

    if (it->first >= 0) {
      c->ps1 = it->second;
    } else {
      c->ps2 = it->second;
    }
    sand_reco::ecal::geometry::CellPosition(geo, c->det, c->mod, c->lay, c->cel,
                                            c->x, c->y,
                                            c->z);  // ok per fluka e geant4
  }

  for (std::map<int, dg_cell>::iterator it = map_cell.begin();
       it != map_cell.end(); ++it) {
    vec_cell.push_back(it->second);
  }
}
}  // namespace ecal

namespace stt
{

// for each tube simulate tdc and adc
// tdc is the time of closest point to wire + drift time
// adc is the sum of energy deposit within integration time window
void create_digits_from_hits(std::map<int, std::vector<hit> >& hits2Tube,
                             std::vector<dg_tube>& digit_vec)
{
  digit_vec.clear();

  for (std::map<int, std::vector<hit> >::iterator it = hits2Tube.begin();
       it != hits2Tube.end(); ++it) {
    double min_time_tub = 1E9;  // mm
    int did = it->first;

    int mod, tub, type, pla, plloc;
    double dwire = 0.;

    sand_reco::stt::decodeSTID(did, pla, tub);
    sand_reco::stt::decodePlaneID(pla, mod, plloc, type);

    TVector2 wire = sand_reco::stt::tubePos[did];

    dg_tube d;
    d.det = it->second[0].det;
    d.did = did;
    d.de = 0;
    d.hor = (type % 2 == 0);
    d.t0 = sand_reco::stt::t0[pla];

    if (d.hor == true) {
      d.x = sand_reco::stt::stt_center[0];
      d.y = wire.Y();
      d.z = wire.X();
      dwire = d.x - 0.5 * sand_reco::stt::stL[did];
    } else {
      d.x = wire.Y();
      d.y = sand_reco::stt::stt_center[1];
      d.z = wire.X();
      dwire = d.y - 0.5 * sand_reco::stt::stL[did];
    }

    for (unsigned int i = 0; i < it->second.size(); i++) {
      double x1 = it->second[i].z1;
      double x2 = it->second[i].z2;
      double t1 = it->second[i].t1;
      double t2 = it->second[i].t2;

      double y1, y2;

      if (type == 2) {
        y1 = it->second[i].y1;
        y2 = it->second[i].y2;
      } else {
        y1 = it->second[i].x1;
        y2 = it->second[i].x2;
      }

      double l = sand_reco::stt::getT(y1, y2, wire.Y(), x1, x2, wire.X());
      double x = x1 + (x2 - x1) * l;
      double y = y1 + (y2 - y1) * l;
      double t = t1 + (t2 - t1) * l;

      TVector2 min_dist_point(x, y);
      double min_dist_hit = (min_dist_point - wire).Mod();
      double min_time_hit = t +
                            (min_dist_hit - sand_reco::stt::wire_radius) /
                                sand_reco::stt::v_drift +
                            dwire / sand_reco::stt::v_signal_inwire;

      if (min_time_hit < min_time_tub) min_time_tub = min_time_hit;

      if (t - d.t0 < sand_reco::stt::stt_int_time) d.de += it->second[i].de;

      d.hindex.push_back(it->second[i].index);
    }

    d.tdc = min_time_tub + rand.Gaus(0, sand_reco::stt::tm_stt_smearing);
    d.adc = d.de;

    digit_vec.push_back(d);
  }
}
}  // namespace stt

}  // namespace digitization