#include <TRandom3.h>

#include "utils.h"

#ifndef SANDDIGITIZATION
#define SANDDIGITIZATION

using namespace sand_reco;

namespace digitization
{

extern TRandom3 rand;

namespace ecal
{
double photo_electron_time_to_pmt_arrival_time(double t0, double d);

void eval_adc_and_tdc_from_photo_electrons(
    std::map<int, std::vector<pe> >& photo_el,
    std::map<int, std::vector<dg_ps> >& map_pmt);

void group_pmts_in_cells(TGeoManager* geo,
                         std::map<int, std::vector<dg_ps> >& ps,
                         std::map<int, double>& L,
                         std::vector<dg_cell>& vec_cell);
}  // namespace ecal

namespace stt
{
void create_digits_from_hits(std::map<int, std::vector<hit> >& hits2Tube,
                             std::vector<dg_tube>& digit_vec);
}  // namespace stt

}  // namespace digitization

#endif