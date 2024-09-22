#include <TRandom3.h>

#include "utils.h"

#ifndef SANDDIGITIZATION
#define SANDDIGITIZATION

using namespace sand_reco;

namespace digitization
{

enum class DETSIM_TYPE {
  kEdepsim,
  kFluka
};
enum class ECAL_digi_mode {
  const_fract,
  fixed_thresh
};

extern TRandom3 rand;

namespace ecal
{
double photo_electron_time_to_pmt_arrival_time(double t0, double d);

void eval_adc_and_tdc_from_photo_electrons(
    std::map<int, std::vector<pe> >& photo_el,
    std::map<int, std::vector<dg_ps> >& map_pmt, ECAL_digi_mode ecal_digi_mode);
}  // namespace ecal

}  // namespace digitization

#endif