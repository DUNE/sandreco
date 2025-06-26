#pragma once

#include "TMatrixD.h"
#include "TVector3.h"
#include "TVectorD.h"
#include <TDecompSVD.h>
#include "TH3D.h"
#include "TH2D.h"
#include <TDecompLU.h>
#include <TStyle.h>
#include "TEllipse.h"
#include "TBox.h"


#include <Math/Factory.h>
#include <Math/Minimizer.h>
#include <Math/Functor.h>
#include <TCanvas.h>
#include <TPolyLine3D.h>
#include <TLine.h>

#include <SANDTrackerCell.h>
#include "SANDTrackerCluster.h"
#include "SANDTrackerDigitCollection.h"

#include <Line3D.h>
#include "utils.h"

#include <algorithm>
#include "iostream"
#include "math.h"

class TrackletFinder {
  public:
    TrackletFinder() {};
    ~TrackletFinder() {};

    void setVolumeParameters(int* p) {volume_parameters_ = p;};
    void setSigmaPosition(double sp) {sigma_pos_ = sp;};
    void setSigmaAngle(double sa)    {sigma_ang_ = sa;};
    void setCells(const sand_reco::tracker::Cluster& cluster) {cluster_ = cluster;};
    void setTrajectory(TVector3 tp, TVector3 td)     {trajectory_ = Line3D(tp, td);};
    void setDigitCollection(sand_reco::tracker::DigitCollection* digit_collection) {digit_collection_ = digit_collection;};

    bool checkParallel(TVector3 d1, TVector3 d2);
    void linesParallelToWire(Line3D w, double distance, std::vector<Line3D>& lines);
    void computeCellsIntersections();
    void computeCellsBands();
    void getScanningAreaVertices();
    void computeDriftTime();

    const std::map<sand_reco::tracker::DigitID, double>& getDigitToDriftTimeMap() const {return digitId_to_drift_time_;};

    std::vector<TVectorD> findTracklets();

    void clear();

    void draw3D();
    void draw3DWires();
    void draw2DWires();
    void draw2DDistance(TFile* h);
    void draw2DDigits();

  private:
    sand_reco::tracker::Cluster cluster_;
    sand_reco::tracker::DigitCollection* digit_collection_;
    std::map<sand_reco::tracker::DigitID, double> digitId_to_drift_time_;
    Line3D trajectory_;

    std::vector<TVector3> cells_intersections_;
    std::vector<Line3D>  cells_bands_;

    double sigma_pos_; // mm
    double sigma_ang_; // rad

    int* volume_parameters_;
    TVector3 mean_point_3d_;

    TCanvas* c2_ = nullptr;
    TCanvas* c3_ = nullptr;

};

