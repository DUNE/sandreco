#include "SANDDigitization.h"

#include "utils.h"

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
              << sand_reco::ecal::scintillation::vlfb* d* conversion::mm_to_m
              << std::endl;
  }

  return time;
}

// from simulated pe produce adc e tdc of calo cell
void eval_adc_and_tdc_from_photo_electrons(
    std::map<int, std::vector<pe> >& photo_el,
    std::map<int, std::vector<dg_ps> >& map_pmt, ECAL_digi_mode ecal_digi_mode)
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
          switch (ecal_digi_mode) {
            case ECAL_digi_mode::const_fract:
              index = int(sand_reco::ecal::acquisition::costant_fraction *
                          pe_count) +
                      start_index;
              if (debug) std::cout << " Const. Fract. " << index << std::endl;
              break;
            case ECAL_digi_mode::fixed_thresh:
              double tdc_thresh =
                  (sand_reco::ecal::acquisition::fixed_thresh_pe >
                   sand_reco::ecal::acquisition::pe_threshold)
                      ? sand_reco::ecal::acquisition::fixed_thresh_pe
                      : sand_reco::ecal::acquisition::pe_threshold;
              index = TMath::Ceil(tdc_thresh) + start_index;
              if (debug) std::cout << " Fix. Thresh. " << index << std::endl;
              break;
          }
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
}  // namespace ecal

}  // namespace digitization