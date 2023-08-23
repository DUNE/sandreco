#include "TGeoManager.h"

#include "TG4Event.h"
#include "TG4HitSegment.h"

#include "struct.h"

#ifndef SANDDIGITIZATIONFLUKA
#define SANDDIGITIZATIONFLUKA

namespace digitization
{

enum class ECAL_digi_mode;

namespace fluka
{

namespace ecal
{
double energy_to_photo_electrons(double E);

bool process_hit(const TG4HitSegment& hit, int& detID, int& modID, int& planeID,
                 int& cellID, double& d1, double& d2, double& t, double& de);

void simulate_photo_electrons(TG4Event* ev, TGeoManager* g,
                              std::map<int, std::vector<pe> >& photo_el,
                              std::map<int, double>& L);

void group_pmts_in_cells(TGeoManager* geo,
                         std::map<int, std::vector<dg_ps> >& ps,
                         std::map<int, double>& L,
                         std::vector<dg_cell>& vec_cell);

void digitize_ecal(TG4Event* ev, TGeoManager* geo,
                   std::vector<dg_cell>& vec_cell,
                   ECAL_digi_mode ecal_digi_mode);

}  // namespace ecal

namespace stt
{
void group_hits_by_tube(TG4Event* ev, TGeoManager* geo, int NHits,
                        Int_t DetType[10000], Float_t xPos[10000],
                        Float_t yPos[10000], Float_t zPos[10000],
                        std::map<int, std::vector<hit> >& hits2Tube);

void create_digits_from_hits(std::map<int, std::vector<hit> >& hits2Tube,
                             std::vector<dg_tube>& digit_vec);

void digitize_stt(TG4Event* ev, TGeoManager* geo, int NHits,
                  Int_t DetType[10000], Float_t xPos[10000],
                  Float_t yPos[10000], Float_t zPos[10000],
                  std::vector<dg_tube>& digit_vec);
}  // namespace stt

void digitize(const char* finname, const char* foutname,
              ECAL_digi_mode ecal_digi_mode);

}  // namespace fluka

}  // namespace digitization

#endif