#include "TGeoManager.h"

#include "TG4Event.h"
#include "TG4HitSegment.h"

#include "SANDGeoManager.h"
#include "struct.h"

#ifndef SANDDIGITIZATIONEDEPSIM
#define SANDDIGITIZATIONEDEPSIM

namespace digitization
{

enum class ECAL_digi_mode ;

namespace edep_sim
{

namespace ecal
{

double energy_to_photo_electrons(double E);

bool process_hit(const SANDGeoManager& g, const TG4HitSegment& hit, int& detID,
                 int& modID, int& planeID, int& cellID, double& d1, double& d2,
                 double& t, double& de);

void simulate_photo_electrons(TG4Event* ev, const SANDGeoManager& g,
                              std::map<int, std::vector<pe> >& photo_el,
                              std::map<int, double>& L);

void group_pmts_in_cells(const SANDGeoManager& geo,
                         std::map<int, std::vector<dg_ps> >& ps,
                         std::map<int, double>& L,
                         std::vector<dg_cell>& vec_cell);

void digitize_ecal(TG4Event* ev, const SANDGeoManager& geo,
                   std::vector<dg_cell>& vec_cell,
                   ECAL_digi_mode  ecal_digi_mode);

}  // namespace ecal

namespace tracker
{
std::vector<TLorentzVector> wireHitClosestPoints(hit& h, const sand_geometry::tracker::WireInfo& wire);

double getMinWireTime(const TLorentzVector& point, const sand_geometry::tracker::WireInfo& wire);

void createDigitsFromHits(const SANDGeoManager& geo,
                             const std::map<sand_geometry::tracker::CellID, std::vector<hit>>& hits2cell,
                             std::vector<dg_wire>& wire_digits);

namespace stt
{
void groupHitsByTube(const TG4Event& ev, const SANDGeoManager& geo,
                        std::map<int, std::vector<hit> >& hits2Tube);

void digitizeStt(const TG4Event& ev, const SANDGeoManager& geo,
                  std::vector<dg_wire>& wire_digits);
}  // namespace stt

namespace chamber
{
void groupHitsByCell(const TG4Event& ev, const SANDGeoManager& geo,
                        std::map<sand_geometry::tracker::CellID, std::vector<hit>>& hits2cell);

void digitizeDrift(const TG4Event& ev, const SANDGeoManager& geo,
                    std::vector<dg_wire>& wire_digits);
}  // namespace chamber
} // namespace tracker

// digitize event
void digitize(const char* finname, const char* foutname,
              ECAL_digi_mode  ecal_digi_mode);

}  // namespace edep_sim

}  // namespace digitization

#endif