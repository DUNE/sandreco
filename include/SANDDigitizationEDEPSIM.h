#include "TGeoManager.h"

#include "TG4Event.h"
#include "TG4HitSegment.h"

#include "struct.h"

#ifndef SANDDIGITIZATIONEDEPSIM
#define SANDDIGITIZATIONEDEPSIM

namespace digitization
{

namespace edep_sim
{

namespace ecal
{

double energy_to_photo_electrons(double E);

bool process_hit(TGeoManager* g, const TG4HitSegment& hit, int& detID,
                 int& modID, int& planeID, int& cellID, double& d1, double& d2,
                 double& t, double& de);

void simulate_photo_electrons(TG4Event* ev, TGeoManager* g,
                              std::map<int, std::vector<pe> >& photo_el,
                              std::map<int, double>& L);

void digitize_ecal(TG4Event* ev, TGeoManager* geo,
                   std::vector<dg_cell>& vec_cell);

}  // namespace ecal

namespace stt
{
void group_hits_by_tube(TG4Event* ev, TGeoManager* geo,
                        std::map<int, std::vector<hit> >& hits2Tube);

void digitize_stt(TG4Event* ev, TGeoManager* geo,
                  std::vector<dg_tube>& digit_vec);
}  // namespace stt

// digitize event
void digitize(const char* finname, const char* foutname);

}  // namespace edep_sim

}  // namespace digitization

#endif