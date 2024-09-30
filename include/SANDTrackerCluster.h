#ifndef SANDTrackerCLUSTER_H
#define SANDTrackerCLUSTER_H

#include <vector>

#include "SANDTrackerDigitCollection.h"
#include "SANDGeoManager.h"

namespace SANDTrackerClusterTools
{

// point
struct Point {
  double x;
  double y;
};

// tube
struct Tube {
  Point center;  // position of the center
  double tdc;    // measured tdc
};

// Tube Collection
struct TubeCollection {
  Tube array[3000];  // tubes belonging to the cluster
};

// line
struct Line {
  double m;  // slope
  double q;  // intercept
};

// cluster parameters
struct ClusterParameters {
  double t0;    // t0
  double ntub;  // number of tubs of the cluster
};

// cluster
struct Cluster {
  ClusterParameters pars;
  TubeCollection tubs;
};

// reconstructed parameters
struct RecoParams {
  Line trk;
  double t0;       // t0
  double quality;  // likelihood
};

// input parameter
struct InputParams {
  Line trk;     // line describing the track
  Cluster cls;  // cluster
};

// evaluate distance between two points
inline double distance(const Point &p1, const Point &p2);

// get minimal distance of line from point
inline double getExpR(const Point &p, const Line &l);

// get minimal distance of line from tube
inline double getExpR(const Tube &t, const Line &l);

// get expected tdc, given a point, a line and a t0
double getExpTdc(const Point &p, const Line &l, double t0);

// get measured tdc
inline double getMeasTdc(const Tube &t);

// translate tdc to reconstructed radius given t0
inline double getMeasR(const Tube &t, double t0);

// evaluate likelihood for a single tube
inline double evalLike(double tobs, double texp, double sigma);

// complete likelihood for
// cluster with hit and not-hit tubes
double likelihood(const double *xx);

// find most distant tubes
std::pair<double, std::pair<int, int>> get_most_distant_tubes_and_dist(
    const Cluster &cls);

// guess slope +
double getGuessSlopePlus(double dx, double dy, double dr, double d);

// guess slope -
double getGuessSlopeMinus(double dx, double dy, double dr, double d);

// get sign to be used in the formula
double getSign(double dx, double dy, double dr, double phi);

// get angle: from line angle to angle of the minimal distance point
double getAlpha(double phi, double sign);

// get point in a circle with center in p
// radius r and angle alpha
Point getPointOnCirle(const Point &p, double r, double alpha);

// find intercept of a line given slope and a point
void setQ(Line &l, const Point &p);

// find all guess lines for a cluster
std::vector<Line> getGuessLines(const Cluster &cls);

}  // namespace SANDTrackerClusterTools

using namespace SANDTrackerClusterTools;

class SANDTrackerClusterID : public SingleElStruct<unsigned long>
{
 public:
  SANDTrackerClusterID(unsigned long id) : SingleElStruct<unsigned long>(id){};
  SANDTrackerClusterID() : SingleElStruct<unsigned long>(){};
};

class SANDTrackerCluster
{
 private:
  int fId;
  const SANDTrackerPlane *fPlane;
  double fMeanX;  // to remove
  std::vector<RecoParams> fParameters;
  std::vector<SANDTrackerDigitID> fDigits;

  static int fCounter;
  static const RecoParams (*fGetRecoFunc)(const SANDTrackerClusterTools::Cluster &cls);

  SANDTrackerCluster(std::vector<SANDTrackerDigitID> &digits, const SANDTrackerPlaneID &id);
  void Reconstruct();
  const Cluster GetExtendedCluster() const;

 public:
  SANDTrackerCluster() = default;
  enum class RecoAlgo { ELikelihood, EMinuit };
  inline int GetId() const { return fId; };
  // inline const SANDTrackerPlaneIndex GetPlaneIndex() const
  // {
    // return fPlane->GetIndex();
  // };
  inline SANDTrackerPlaneID GetPlaneId() const 
  { 
    // return fPlane->GetId(); 
  };
  // inline SANDTrackerPlane::EOrientation GetOrientation() const
  // {
    // return fPlane->GetOrientation();
  // };
  inline double GetMeanX() const { return fMeanX; };
  inline double GetZ() const 
  { 
    // return fPlane->GetZ(); 
  };
  inline const std::vector<SANDTrackerDigitID> &GetDigits() const { return fDigits; };
  static void ResetCounter() { fCounter = 0; };
  inline const std::vector<RecoParams> &GetRecoParameters() const
  {
    return fParameters;
  };
  double GetDigitCoordinate(const SANDTrackerDigit &d) const
  {
    // return (GetOrientation() == SANDTrackerPlane::EOrientation::kHorizontal) ? d.y
                                                                    //  : d.x;
  };
  static void SetRecoAlgo(const RecoAlgo &algo);

  friend class SANDTrackerClustersInPlane;
};

#ifdef __MAKECINT__
#pragma link C++ class SANDTrackerClusterTools::Point + ;
#pragma link C++ class SANDTrackerClusterTools::Tube + ;
#pragma link C++ class SANDTrackerClusterTools::TubeCollection + ;
#pragma link C++ class SANDTrackerClusterTools::Line + ;
#pragma link C++ class std::vector < SANDTrackerClusterTools::Line> + ;
#pragma link C++ class SANDTrackerClusterTools::ClusterParameters + ;
#pragma link C++ class SANDTrackerClusterTools::Cluster + ;
#pragma link C++ class SANDTrackerClusterTools::RecoParams + ;
#pragma link C++ class std::vector < SANDTrackerClusterTools::RecoParams> + ;
#pragma link C++ class SANDTrackerClusterTools::InputParams + ;
#pragma link C++ class SANDTrackerPlane + ;
#pragma link C++ class SingleElStruct<unsigned int> + ;
#pragma link C++ class SANDTrackerDigitID + ;
#pragma link C++ class SANDTrackerPlaneIndex + ;
#pragma link C++ class std::vector < SANDTrackerDigitID> + ;
#pragma link C++ class SANDTrackerTubeID + ;
#pragma link C++ class SANDTrackerTube + ;
#pragma link C++ class std::map<SANDTrackerTubeID,SANDTrackerTube> + ;
#pragma link C++ class SANDTrackerPlaneID + ;
#pragma link C++ class SANDTrackerDigit + ;
#pragma link C++ class SANDTrackerCluster + ;
#pragma link C++ class std::vector < SANDTrackerCluster> + ;
#endif

#endif