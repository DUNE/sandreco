#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>

#include <TFile.h>
#include <TTree.h>

#include "TG4Event.h"
#include "TStyle.h"
#include "SANDGeoManager.h"
#include "SANDTrackerDigitCollection.h"
#include "SANDTrackerModule.h"
#include "TrackBuilder.h"
#include "SANDTrackerCluster.h"
#include "SANDTrackerClusterCollection.h"
#include "SANDTrackerUtils.h"
#include "SANDProcessTracklets.h"


struct TripletsInfo{
  int run= 0;
  int event= 0;
  int moduleID = 0;
  double x_true=0.0, y_true=0.0, z_true=0.0;
  double x=0.0, y=0.0, z=0.0;
  double score = 0.0;
  double D =0.0;

  double dx=0.0, dy=0.0;

  int nU = 0;
  int nV = 0;
  int nY = 0;
};

struct TripletsTree{
  TFile* file = nullptr;
  TTree* tree = nullptr;
  TripletsInfo triplets_info;
};

void create_triplets_tree(const std::string& filename, TripletsTree& tripTree)
{
  tripTree.file = TFile::Open(filename.c_str(), "RECREATE");
  if (!tripTree.file || tripTree.file->IsZombie()) {
    throw std::runtime_error("Cannot create ROOT file: " + filename);
  }

  tripTree.tree = new TTree("triplets", "Triplets per event");
  auto& ti = tripTree.triplets_info;

  tripTree.tree->Branch("run",      &ti.run);
  tripTree.tree->Branch("event",    &ti.event);
  tripTree.tree->Branch("moduleID", &ti.moduleID);

  tripTree.tree->Branch("x_true",   &ti.x_true);
  tripTree.tree->Branch("y_true",   &ti.y_true);
  tripTree.tree->Branch("z_true",   &ti.z_true);

  tripTree.tree->Branch("z",        &ti.z);
  tripTree.tree->Branch("x",        &ti.x);
  tripTree.tree->Branch("y",        &ti.y);
  
  tripTree.tree->Branch("score",    &ti.score);
  tripTree.tree->Branch("D",        &ti.D);
  
  tripTree.tree->Branch("dx",       &ti.dx);  
  tripTree.tree->Branch("dy",       &ti.dy);

  tripTree.tree->Branch("nU",       &ti.nU);
  tripTree.tree->Branch("nV",       &ti.nV);
  tripTree.tree->Branch("nY",       &ti.nY);
}

void close_triplets_tree(TripletsTree& tripTree)
{
  if (!tripTree.file) return;

  tripTree.file->cd();
  
  if (tripTree.tree) tripTree.tree->Write();
  tripTree.file->Write();
  tripTree.file->Close();

  tripTree.tree = nullptr;
  tripTree.file = nullptr;
}


namespace sand_reco
{
namespace kf
{
namespace trackbuilding
{

// inline sand_geometry::tracker::plane_iterator getPlane() const { return
// plane_; }; inline sand_geometry::tracker::PlaneID getPlaneId() const { return
// plane_->uId(); }; inline double getRotation() const { return
// plane_->getRotation(); }; //U.V.Y inline double getZ() const { return
// plane_->getPosition().Z(); };

// Da capire: come avere rotazioni del piano corrente, esiste gia?
enum class PlaneOrientation {
  kY,
  kV,
  kU,
  kBho
};

inline PlaneOrientation
    getPlaneOrientation(double rotation, double stereo_angle)
{
  const double theta_threshold = 2.0 * M_PI / 180.0; //2 gradi
  //
  if (std::fabs(rotation) < theta_threshold) {
    return PlaneOrientation::kY;
  }
  if (std::fabs(rotation - stereo_angle) < theta_threshold) {
    return PlaneOrientation::kV;
  }
  if (std::fabs(rotation + stereo_angle) < theta_threshold) {
    return PlaneOrientation::kU;
  }
  return PlaneOrientation::kBho;
}

/// Calcola un id del modulo a partire da un Cluste usando direttamente la
/// geometria: Cluster -> plane_iterator -> Plane::getModule() -> Module::Id()
// moduli forse ce un modo piu intelligente chiedi. Module contine egia una
// mappa std::map<PlaneID, Plane> non so se riesco ad arrivare al module
// aprtendo dal cluster guarda TrackerModule bho 
inline int getModuleID(const tracker::Cluster& cluster)
{
 auto planeIt = cluster.getPlane();
 const auto& plane = *planeIt;

  // PlaneID ultimi due digit del modulo identificano l'orientamento U,V,Y 
  const auto planeId = plane.uId();
  const unsigned long unsigned_planeId = planeId();
  const int moduleID = static_cast<int>(unsigned_planeId / 100);

    // std::cout << "[TrackBuilder] planeId=" << raw
    //         << " -> moduleID=" << moduleID << "\n";

  return moduleID;
}

// costruiamo un punto (x,y) per ogni cluster in goni modulo . prendo baricentro
// dei centri dei fili. penso eissta gia tra l'altro
std::optional<std::pair<double, double>> getClusterCentroidXY(
    const sand_reco::tracker::Cluster& cluster, const SANDGeoManager* sand_geo)
{
  // auto planeIt = cluster.getPlane();
  const auto& digitsID = cluster.getDigits();
  if (digitsID.empty()) {
    return std::nullopt;
  }

  double sumX = 0.0;
  double sumY = 0.0;
  int n = 0;

  // loop sui digit -> cellid-> centro del filo
  for (const auto& digId : digitsID) {
    // auto cellIt =
    //     sand_geo->getCellInfo(sand_geometry::tracker::CellID(digId()));
    // const auto& cell = cellIt->second;

    // TVector3 center = cell.getWire().getCenter();
    // sumX += center.X();
    // sumY += center.Y();
    const auto& digit = sand_reco::tracker::DigitCollection::getDigit(digId);
    sumX += digit.x;
    sumY += digit.y;
    ++n;
  }

  if (n == 0) {
    return std::nullopt;
  }

  return std::make_pair(sumX / n, sumY / n);
}


//Get Truth for each clsuter hit segment detector, una Truth è pos_,dir_ e mom_ da ProcessTracklets
std::optional<Truth> getTruth(const Triplet& triplet){

  const sand_reco::tracker::Cluster* cluster = nullptr;
  if(triplet.clusterU) cluster = triplet.clusterU;
  else if (triplet.clusterV) cluster = triplet.clusterV;
  else if (triplet.clusterY) cluster = triplet.clusterY;

  if(!cluster){ 
    return std::nullopt; 
  }

  TVector3 first_point;
  TVector3 last_point;
  TVector3 first_point_mom;
  TVector3 last_point_mom;

  double z_min = 1e10;
  double z_max = -1e10;
   
  for (const auto& digitId : cluster->getDigits()){
    auto digit = sand_reco::tracker::DigitCollection::getDigit(digitId);

    if(digit.z > z_max){
      z_max = digit.z;
      last_point = TVector3(digit.x, digit.y, digit.z);
      last_point_mom = TVector3(digit.px, digit.py, digit.pz);
    }
    if(digit.z < z_min){
      z_min = digit.z;
      first_point = TVector3(digit.x, digit.y, digit.z);
      first_point_mom = TVector3(digit.px, digit.py, digit.pz);
    }
  }
  //la key della mappa è la z del triplet
  double z_to_triplet = triplet.triplet_pos_.Z();
  
  Truth truth = getTrueTrackletOfCluster(
      first_point, last_point, first_point_mom, last_point_mom,
      z_to_triplet);

  return truth;
}




// buil triplet map from cluster for each module. per ogni cluster un triplet è
// costruito come un punto misurato nello spazio (x,y,z) dove x,y preso dal
// baricentro dei fili e z dalla posizione del cluster cluster.GetZ() e
// Plane::GetModule()->Id(). chiave z valore triplet etc.
//  cambio di coordinate da locali a modulo U,V,Y (getRotation() rtc..?)
// associo anche uno score calcolo le combinazioni e seleziono in base ad uno
// score come in NOMAD

// chiave z valore triplet(x,y,z,score) (score lo calcolo dopo)
TripletMap buildTripletFromClusters(const SANDGeoManager* sand_geo,
                                    const std::vector<sand_reco::tracker::Cluster>& clusters,
                                    double stereo_angle){
  TripletMap z_to_triplet;
  if (!sand_geo) {
    return z_to_triplet;
  }

  // separo i cluster per module view
  std::map<int, std::vector<const tracker::Cluster*>> cluster_u;
  std::map<int, std::vector<const tracker::Cluster*>> cluster_v;
  std::map<int, std::vector<const tracker::Cluster*>> cluster_y;

  for (const auto& cluster : clusters) {
    auto centroidXY = getClusterCentroidXY(cluster, sand_geo);
    if (!centroidXY) continue;

    int moduleID = getModuleID(cluster);
    double rotation = cluster.getRotation();

    PlaneOrientation orientation = getPlaneOrientation(rotation, stereo_angle);

    //--------------------------DEBUG----------------------
                    std::string orientationName;
                switch (orientation) {
                  case PlaneOrientation::kU: orientationName = "U"; break;
                  case PlaneOrientation::kV: orientationName = "V"; break;
                  case PlaneOrientation::kY: orientationName = "Y"; break;
                  default:       orientationName = "Unknown"; break;
                }

                // std::cout << "[TrackBuilder] cluster planeId=" << cluster.getPlaneId()()
                //           << " moduleID=" << moduleID
                //           << " rotation=" << rotation
                //           << " -> view=" << orientationName
                //           << std::endl;
    //----------------------------------------------------------
    switch (orientation) {
      case PlaneOrientation::kU:
        cluster_u[moduleID].push_back(&cluster);
        break;

      case PlaneOrientation::kV:
        cluster_v[moduleID].push_back(&cluster);
        break;

      case PlaneOrientation::kY:
        cluster_y[moduleID].push_back(&cluster);
        break;

      case PlaneOrientation::kBho:  
      default:
        std::cout << "[TrackBuilder] WARNING cluster in module :" << moduleID << "has invalid orientation= " << rotation << "\n";
        break;
    }
  }

  // per ogni modulo costruisco le triplet UVY

  // for (const auto& [moduleID, u_clusters] : cluster_u) {
  //   auto v_clusters_it = cluster_v.find(moduleID);
  //   auto y_clusters_it = cluster_y.find(moduleID);
  //       std::size_t nU = u_clusters.size();
  //   std::size_t nV = (v_clusters_it != cluster_v.end()) ? v_clusters_it->second.size() : 0;
  //   std::size_t nY = (y_clusters_it != cluster_y.end()) ? y_clusters_it->second.size() : 0;

  //   std::cout << "[TrackBuilder] Module " << moduleID
  //             << " has nU=" << nU
  //             << " nV=" << nV
  //             << " nY=" << nY << "\n";

  //   if (v_clusters_it == cluster_v.end() || y_clusters_it == cluster_y.end()) {
  //     continue;  // non trova v o y corrispondenti
  //   }
  //   const auto& v_clusters = v_clusters_it->second;
  //   const auto& y_clusters = y_clusters_it->second;

  std::set<int> allModules;
  for (const auto& kv : cluster_u) allModules.insert(kv.first);
  for (const auto& kv : cluster_v) allModules.insert(kv.first);
  for (const auto& kv : cluster_y) allModules.insert(kv.first);

  std::cout << "[TrackBuilder] nModules with at least one cluster = " << allModules.size() << "\n";

  for (int moduleID : allModules) {
    const auto itU = cluster_u.find(moduleID);
    const auto itV = cluster_v.find(moduleID);
    const auto itY = cluster_y.find(moduleID);

    const auto& u_clusters = (itU != cluster_u.end()) ? itU->second
                                                      : std::vector<const tracker::Cluster*>{};
    const auto& v_clusters = (itV != cluster_v.end()) ? itV->second
                                                      : std::vector<const tracker::Cluster*>{};
    const auto& y_clusters = (itY != cluster_y.end()) ? itY->second
                                                      : std::vector<const tracker::Cluster*>{};

    std::size_t nU = u_clusters.size();
    std::size_t nV = v_clusters.size();
    std::size_t nY = y_clusters.size();

    // std::cout << "[TrackBuilder] Module " << moduleID
    //           << " has nU=" << nU
    //           << " nV=" << nV
    //           << " nY=" << nY << "\n";

    if (nU == 0 || nV == 0 || nY == 0) {
      continue;
    }


    // SANDTrackerTRACKRECO_LOG("INFO",
    // "Module " << moduleID
    //         << " nU=" << u_clusters.size()
    //         << " nV=" << v_clusters.size()
    //         << " nY=" << y_clusters.size());


    double best_score = 1e30;
    bool find_best = false;
    Triplet best_triplet;

    // loop su tutte le combinazioni
    for (const auto* u_cluster : u_clusters) {
      auto u_centroid = getClusterCentroidXY(*u_cluster, sand_geo);
      if (!u_centroid) continue;

      const auto [xU, yU] = *u_centroid;
      double zU = u_cluster->getZ();

      for (const auto* v_cluster : v_clusters) {
        auto v_centroid = getClusterCentroidXY(*v_cluster, sand_geo);
        if (!v_centroid) continue;

        const auto [xV, yV] = *v_centroid;
        double zV = v_cluster->getZ();

        for (const auto* y_cluster : y_clusters) {
          auto y_centroid = getClusterCentroidXY(*y_cluster, sand_geo);
          if (!y_centroid) continue;

          const auto [xY, yY] = *y_centroid;
          double zY = y_cluster->getZ();

          // cambio coordinate
          const double U = yU * cos(stereo_angle) - xU * sin(stereo_angle);
          const double V = yV * cos(stereo_angle) + xV * sin(stereo_angle);
          const double Y = yY;

          double D = U + V - 2.0 * Y * cos(stereo_angle);

          if (std::fabs(D) < best_score) {
          best_score = fabs(D);
          find_best = true;

          const double x_triplet = (V - U) / (2.0 * sin(stereo_angle));
          const double y_triplet = (U + V) / (2.0 * cos(stereo_angle));
          const double z_triplet = (zU + zY + zV) / 3.0;  // grande dubbio, prendo la media in z o direttamente z centrale? prima usavamo quella piu vicina al mctruth

          Triplet triplet;
          triplet.moduleID = moduleID;
          triplet.triplet_pos_ = TVector3(x_triplet, y_triplet, z_triplet);
          triplet.score = best_score;
          triplet.D = D;
          triplet.clusterU = u_cluster;
          triplet.clusterV = v_cluster;
          triplet.clusterY = y_cluster;

          best_triplet = triplet;

          }  //if 
        }  // loop Y
      }  // loop v
    }  // loop u

    if (find_best) {
      double z_module = best_triplet.triplet_pos_.Z();
      z_to_triplet[z_module] = best_triplet;

    }
  }

  std::cout << " builTripletsFromCluster done, found n triplets: "
  << z_to_triplet.size() << "\n";

  return z_to_triplet;
}

}  // namespace trackbuilding
}  // namespace kf
}  // namespace sand_reco

//-----------------------------------------------------------------------------
//                               PROCESS EVENTS
//-----------------------------------------------------------------------------

void processEventWithTriplets(const SANDGeoManager* sand_geo,
                              std::vector<dg_wire>* digits,
                              int run_number,
                              int event_number,
                              TripletsTree& triplets_tree)
{
  using namespace sand_reco;
  using namespace sand_reco::tracker;

  std::cout << "\n[TrackBuilder] === Event " << event_number << " ===\n";

  if (!digits || digits->empty()) {
    std::cout << "[TrackBuilder] Event " << event_number << ": no digits, skipping.\n";
    return;
  }

  //FILL DIGIT MAP
  sand_reco::tracker::DigitCollection::fillMap(digits);
  auto digit_vec = DigitCollection::getDigits();
  if (digit_vec.empty()) {
    std::cout << "[TrackBuilder] Event " << event_number
              << ": DigitCollection empty, skipping.\n";
    return;
  }

  //CLUSTERS
  sand_reco::tracker::ClusterCollection clusters(
      sand_geo,
      digit_vec,
      sand_reco::tracker::ClusterCollection::ClusteringMethod::kCellAdjacency);

  std::vector<Cluster> allClusters;
  for (const auto& container : clusters.getContainers()) {
    for (const auto& c : container->getClusters()) {
      allClusters.push_back(c);
    }
  }

  std::cout << "[TrackBuilder] Event " << event_number
            << ": nClusters=" << allClusters.size() << std::endl;
          // --- DEBUG---

          // std::set<double> rotations_rad;
          // for (const auto& c : allClusters) {
          //   rotations_rad.insert(c.getRotation());
          // }

          // std::cout << " stereo_angle=" << sand_reco::kf::trackbuilding::stereoAngleInRad << " rad"
          //           << " -> unique plane rotations in clusters (rad, deg): ";

          // for (double rot : rotations_rad) {
          //   double deg = rot * 180.0 / M_PI;
          //   std::cout << " [ " << rot << " , " << deg << "° ]";
          // }
          // std::cout << std::endl;

  //BUILD TRIPLETS MAP
  auto tripletMap =
      sand_reco::kf::trackbuilding::buildTripletFromClusters(
          sand_geo,
          allClusters,
          sand_reco::kf::trackbuilding::stereoAngleInRad);

  std::cout << "[TrackBuilder] Event " << event_number
          << ": nTriplets=" << tripletMap.size()
          << "\n";

  for (const auto& [z, triplet] : tripletMap) {
    std::cout << "  Triplet: z=" << z
              << "  module=" << triplet.moduleID
              << "  x=" << triplet.triplet_pos_.X()
              << "  y=" << triplet.triplet_pos_.Y()
              << "  score=" << triplet.score
              << "\n";
  }


  //FILL TREE INFO
  for (const auto& [z, trip] : tripletMap){
    auto& ti = triplets_tree.triplets_info;

    ti.run      = run_number;
    ti.event    = event_number;
    ti.moduleID = trip.moduleID;
    ti.z        = trip.triplet_pos_.Z();
    ti.x        = trip.triplet_pos_.X();
    ti.y        = trip.triplet_pos_.Y();
    ti.score    = trip.score;
  
    ti.nU = 0.0;
    ti.nV = 0.0;
    ti.nY = 0.0;

    double x_true = 0.0;
    double y_true = 0.0;

    auto truth_opt = getTruth(trip);
    if(truth_opt){
      const auto& truth = *truth_opt;
      ti.x_true = truth.pos_.X();
      ti.y_true = truth.pos_.Y();
      ti.dx = ti.x - ti.x_true;
      ti.dy = ti.y - ti.y_true;
    } else {
      ti.x_true = -999999.;
      ti.y_true = -999999.;
      ti.dx = -999999.;
      ti.dy = -999999.;
    }



    triplets_tree.tree->Fill();
  }

}


//-----------------------------------------------------------------------------
//                                    MAIN
//-----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  gStyle->SetOptStat(0);

  if (argc < 3) {
    std::cerr << "Usage: TrackBuilder edepsim.root digits.root" << std::endl;
    return 1;
  }

  //edepsim file
  TFile f(argv[1], "READ");
  if (f.IsZombie()) {
    std::cerr << "[TrackBuilder] Error opening " << argv[1] << std::endl;
    return 1;
  }

  TGeoManager* geo = (TGeoManager*)f.Get("EDepSimGeometry");
  if (!geo) {
    std::cerr << "[TrackBuilder] Cannot find TGeoManager 'EDepSimGeometry' in "
              << argv[1] << std::endl;
    return 1;
  }

  //MC info tree
  TTree* t_h = (TTree*)f.Get("EDepSimEvents");
  TG4Event* ev = nullptr;
  if (t_h) {
    ev = new TG4Event;
    t_h->SetBranchAddress("Event", &ev);
  }

  //Digit file
  TFile f_d(argv[2], "READ");
  if (f_d.IsZombie()) {
    std::cerr << "[TrackBuilder] Error opennig " << argv[2] << std::endl;
    return 1;
  }

  TTree* t = (TTree*)f_d.Get("tDigit");
  if (!t) {
    std::cerr << "[TrackBuilder] Cannot find TTree 'tDigit' in "
              << argv[2] << std::endl;
    return 1;
  }

  std::vector<dg_wire>* digits = nullptr;
  t->SetBranchAddress("dg_wire", &digits);

  //Geomanager
  SANDGeoManager sand_geo;
  sand_geo.init(geo);

  std::string geometry;
  if (geo->FindVolumeFast("STTtracker_PV")) {
    geometry = "STT";
  } else if (geo->FindVolumeFast("SANDtracker_PV")) {
    geometry = "DRIFT";
  } else {
    std::cout << "[TrackBuilder] WARNING: Unknown geometry, default to STT\n";
    geometry = "STT";
  }
  sand_geo.fillAdjacentCells(geometry);

  TripletsTree triplets_tree;
  create_triplets_tree("triplets.root", triplets_tree);
 

  int nev = t->GetEntries();
  for (int i =10; i < 20; i++) {
    if (t_h) t_h->GetEntry(i);
    t->GetEntry(i);

    int run_number   = 0;
    int event_number = i;
    processEventWithTriplets(&sand_geo, digits, run_number, event_number, triplets_tree);
  }

  close_triplets_tree(triplets_tree);
  return 0;
}
