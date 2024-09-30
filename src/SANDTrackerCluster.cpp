#include "SANDTrackerCluster.h"
// #include "STTClusterTrackReco.h"
#include "SANDTrackerUtils.h"

#include "Math/Functor.h"
#include <TMath.h>

#include <numeric>

// evaluate distance between two points
inline double SANDTrackerClusterTools::distance(const Point &p1, const Point &p2)
{
  return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

// get minimal distance of line from point
inline double SANDTrackerClusterTools::getExpR(const Point &p, const Line &l)
{
  return fabs((-p.y + p.x * l.m + l.q) / sqrt(1 + l.m * l.m));
}

// get minimal distance of line from tube
inline double SANDTrackerClusterTools::getExpR(const Tube &t, const Line &l)
{
  return getExpR(t.center, l);
}

// get expected tdc, given a point, a line and a t0
double SANDTrackerClusterTools::getExpTdc(const Point &p, const Line &l, double t0)
{
  auto r = getExpR(p, l);
  // return t0 + r / SANDTrackerUtils::GetSTTElectronDriftVelocity();
}

// get measured tdc
inline double SANDTrackerClusterTools::getMeasTdc(const Tube &t) { return t.tdc; }

// translate tdc to reconstructed radius given t0
inline double SANDTrackerClusterTools::getMeasR(const Tube &t, double t0)
{
  // return (t.tdc - t0) * SANDTrackerUtils::GetSTTElectronDriftVelocity();
}

// evaluate likelihood for a single tube
inline double SANDTrackerClusterTools::evalLike(double tobs, double texp, double sigma)
{
  return pow((tobs - texp) / sigma, 2);
}

// complete likelihood for
// cluster with hit and not-hit tubes
double SANDTrackerClusterTools::likelihood(const double *xx)
{
  double like = 0.;
  // auto par = reinterpret_cast<const InputParams *>(xx);

  // for (Int_t i = 0; i < par->cls.pars.ntub; i++) {
  //   auto like_contrib = 0.;

  //   const Double_t tdc_obs = par->cls.tubs.array[i].tdc;
  //   const Double_t tdc_exp =
  //       getExpTdc(par->cls.tubs.array[i].center, par->trk, par->cls.pars.t0);

  //   if (tdc_obs == -1.) {
  //     const auto drift_time = tdc_exp - par->cls.pars.t0;
  //     if (drift_time <= SANDTrackerUtils::GetTubeMaxDriftTime()) {
  //       like_contrib =
  //           pow((drift_time - SANDTrackerUtils::GetTubeMaxDriftTime()) / 4., 2);
  //     }
  //   } else {
  //     like_contrib = pow((tdc_obs - tdc_exp) / 4., 2);
  //   }

  //   like += like_contrib;
  // }

  return like;
}

// find most distant tubes
std::pair<double, std::pair<int, int>>
    SANDTrackerClusterTools::get_most_distant_tubes_and_dist(const Cluster &cls)
{
  std::map<double, std::pair<int, int>> dists;
  for (int i = 0; i < cls.pars.ntub - 1; i++)
    for (int j = i + 1; j < cls.pars.ntub; j++)
      if (cls.tubs.array[i].tdc != -1. && cls.tubs.array[j].tdc != -1.)
        dists[distance(cls.tubs.array[i].center, cls.tubs.array[j].center)] =
            std::make_pair(i, j);

  return *dists.rbegin();
}

// guess slope +
double SANDTrackerClusterTools::getGuessSlopePlus(double dx, double dy, double dr,
                                          double d)
{
  return (-dx * dy + dr * sqrt(d * d - dr * dr)) / (dr * dr - dx * dx);
}

// guess slope -
double SANDTrackerClusterTools::getGuessSlopeMinus(double dx, double dy, double dr,
                                           double d)
{
  return (-dx * dy - dr * sqrt(d * d - dr * dr)) / (dr * dr - dx * dx);
}

// get sign to be used in the formula
double SANDTrackerClusterTools::getSign(double dx, double dy, double dr, double phi)
{
  return (dx * sin(phi) - dy * cos(phi)) / dr;
}

// get angle: from line angle to angle of the minimal distance point
double SANDTrackerClusterTools::getAlpha(double phi, double sign)
{
  return phi + sign * M_PI * 0.5;
}

// get point in a circle with center in p
// radius r and angle alpha
Point SANDTrackerClusterTools::getPointOnCirle(const Point &p, double r, double alpha)
{
  return {p.x + r * cos(alpha), p.y + r * sin(alpha)};
}

// find intercept of a line given slope and a point
void SANDTrackerClusterTools::setQ(Line &l, const Point &p) { l.q = p.y - l.m * p.x; }

// find all guess lines for a cluster
std::vector<Line> SANDTrackerClusterTools::getGuessLines(const Cluster &cls)
{
  std::vector<Line> l(4);

  auto tubes_and_dist = get_most_distant_tubes_and_dist(cls);
  
  for(int k = 0; k < cls.pars.ntub; k++)
    auto& tb = cls.tubs.array[k];

  auto d = tubes_and_dist.first;
  const Tube t[2] = {cls.tubs.array[tubes_and_dist.second.first],
                     cls.tubs.array[tubes_and_dist.second.second]};
  auto dx = t[1].center.x - t[0].center.x;
  auto dy = t[1].center.y - t[0].center.y;
  auto r1 = getMeasR(t[0], cls.pars.t0);
  auto r2 = getMeasR(t[1], cls.pars.t0);
  auto dr_plus = r2 + r1;
  auto dr_minus = r2 - r1;

  double dradius[] = {dr_minus, dr_minus, dr_plus, dr_plus};

  l[0].m = getGuessSlopePlus(dx, dy, dradius[0], d);
  l[1].m = getGuessSlopeMinus(dx, dy, dradius[1], d);
  l[2].m = getGuessSlopePlus(dx, dy, dradius[2], d);
  l[3].m = getGuessSlopeMinus(dx, dy, dradius[3], d);

  double phi[] = {atan(l[0].m), atan(l[1].m), atan(l[2].m), atan(l[3].m)};

  setQ(l[0],
       getPointOnCirle(t[0].center, r1,
                       getAlpha(phi[0], getSign(dx, dy, dradius[0], phi[0]))));
  setQ(l[1],
       getPointOnCirle(t[0].center, r1,
                       getAlpha(phi[1], getSign(dx, dy, dradius[1], phi[1]))));
  setQ(l[2],
       getPointOnCirle(t[0].center, r1,
                       getAlpha(phi[2], -getSign(dx, dy, dradius[2], phi[2]))));
  setQ(l[3],
       getPointOnCirle(t[0].center, r1,
                       getAlpha(phi[3], -getSign(dx, dy, dradius[3], phi[3]))));

  return l;
}

// find all solution -> optimize them all -> choose best solution
// !!! NO FILTER !!!
const RecoParams GetMinuitBasedOptimizedPar(const Cluster &cls)
{
  // auto guesses = getGuessLines(cls);

  // std::vector<RecoParams> reco_pars;

  // for (auto const &g : guesses) {
  //   InputParams ipar;
  //   ipar.trk = g;
  //   ipar.cls = cls;
  //   reco_pars.push_back(STTClusterTrackReco::Reconstruct(ipar));
  // }

  // std::sort(reco_pars.begin(), reco_pars.end(),
  //           [](const RecoParams &p1, const RecoParams &p2) {
  //             return p1.quality < p2.quality;
  //           });

  // return reco_pars.front();
}

const RecoParams GetBestLikelihoodPar(const Cluster &cls)
{
  auto guesses = getGuessLines(cls);

  std::vector<RecoParams> reco_pars;

  for (auto const &g : guesses) {

    InputParams ipar;
    ipar.trk = g;
    ipar.cls = cls;

    RecoParams opar;
    opar.trk = g;
    opar.t0 = cls.pars.t0;
    opar.quality =
        SANDTrackerClusterTools::likelihood(reinterpret_cast<double *>(&ipar));

    reco_pars.push_back(opar);
  }

  std::sort(reco_pars.begin(), reco_pars.end(),
            [](const RecoParams &p1, const RecoParams &p2) {
              return p1.quality < p2.quality;
            });

  return reco_pars.front();
}

int SANDTrackerCluster::fCounter = 0;
const RecoParams (*SANDTrackerCluster::fGetRecoFunc)(
    const SANDTrackerClusterTools::Cluster &cls) = GetBestLikelihoodPar;

SANDTrackerCluster::SANDTrackerCluster(std::vector<SANDTrackerDigitID> &digits, const SANDTrackerPlaneID &id)
    : fId(fCounter++), fDigits(digits)
{
  // fPlane = &STTStrawTubeTracker::GetPlane(id);
  // fMeanX =
  //     std::accumulate(fDigits.begin(), fDigits.end(), 0.,
  //                     [this](double x, const SANDTrackerDigitID &d) {
  //                       return x + this->GetDigitCoordinate(
  //                                      SANDTrackerDigitCollection::GetDigit(d));
  //                     }) /
  //     fDigits.size();

  // if (fDigits.size() > 1) Reconstruct();
}

const Cluster SANDTrackerCluster::GetExtendedCluster() const
{
  Cluster cls;

  // auto count = 0;

  // std::vector<SANDTrackerCellID> ids;

  // for (auto i = 0u; i < fDigits.size(); i++) {
  //   auto d = SANDTrackerDigitCollection::GetDigit(fDigits.at(i));
  //   cls.tubs.array[i].center.x = d.z;
  //   cls.tubs.array[i].center.y = (d.hor) ? d.y : d.x;
  //   cls.tubs.array[i].tdc = d.tdc;

  //   ids.push_back(SANDTrackerUtils::GetTubeID(d.did));

  //   count++;
  // }

  // SANDTrackerPlaneID planeid = SANDTrackerUtils::GetPlaneID(fDigits.front());

  // std::sort(ids.begin(), ids.end());

  // auto plane = STTStrawTubeTracker::GetPlane(planeid);

  // auto n_tubes_in_module = plane.GetN();
  // auto const id_max = std::min(ids.back()() + 10, n_tubes_in_module - 1);
  // auto const id_min = std::max(ids.front()() - 10, 0u);

  // for (auto this_id = id_min; this_id <= id_max; this_id++) {

  //   if (std::find(ids.begin(), ids.end(), this_id) == ids.end()) {
  //     auto tb = plane.Get(this_id);
  //     cls.tubs.array[count].center.x = tb.GetCenter().Z();
  //     cls.tubs.array[count].center.y = (plane.GetOrientation() == STTPlane::EOrientation::kHorizontal) ? tb.GetCenter().Y() : tb.GetCenter().X();
  //     cls.tubs.array[count].tdc = -1.;
  //     count++;
  //   }
  // }

  // cls.pars.ntub = count;

  return cls;
}

void SANDTrackerCluster::SetRecoAlgo(const RecoAlgo &algo)
{
  switch (algo) {
    case RecoAlgo::ELikelihood:
      fGetRecoFunc = GetBestLikelihoodPar;
      break;
    case RecoAlgo::EMinuit:
      fGetRecoFunc = GetMinuitBasedOptimizedPar;
      break;
  }
}

void SANDTrackerCluster::Reconstruct()
{
  auto cls = GetExtendedCluster();

  auto opar = fGetRecoFunc(cls);

  RecoParams par;
  par.trk = opar.trk;
  par.quality = opar.quality;
  par.t0 = opar.t0;

  fParameters.clear();
  fParameters.push_back(par);
}