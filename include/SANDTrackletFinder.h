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

#include <CLine3D.h>
#include "utils.h"

#include <algorithm>
#include "iostream"
#include "math.h"

class TrackletFinder {
  public:
    TrackletFinder() {};
    ~TrackletFinder() {};

    void SetVolumeParameters(int* p) {_volume_parameters = p;};
    void SetSigmaPosition(double sp) {_sigma_pos = sp;};
    void SetSigmaAngle(double sa)    {_sigma_ang = sa;};
    void SetCells(const SANDTrackerCluster& cluster) {_cluster = cluster;};
    void SetTrajectory(TVector3 tp, TVector3 td)     {_trajectory = CLine3D(tp, td);};
    void SetDigitCollection(SANDTrackerDigitCollection* digit_collection) {_digit_collection = digit_collection;};

    bool CheckParallel(TVector3 d1, TVector3 d2);
    void LinesParallelToWire(CLine3D w, double distance, std::vector<CLine3D>& lines);
    void ComputeCellsIntersections();
    void ComputeCellsBands();
    void GetScanningAreaVertices();
    void ComputeDriftTime();

    const std::map<SANDTrackerDigitID, double>& GetDigitToDriftTimeMap() const {return _digitId_to_drift_time;};

    std::vector<TVectorD> FindTracklets();

    void Clear();

    void Draw3D();
    void Draw3DWires();
    void Draw2DWires();
    void Draw2DDistance();
    void Draw2DDigits();

  private:
    SANDTrackerCluster _cluster;
    SANDTrackerDigitCollection* _digit_collection;
    std::map<SANDTrackerDigitID, double> _digitId_to_drift_time;
    CLine3D _trajectory;

    std::vector<TVector3> _cells_intersections;
    std::vector<CLine3D>  _cells_bands;

    double _sigma_pos; // mm
    double _sigma_ang; // rad

    int* _volume_parameters;

    TCanvas* _c2 = nullptr;
    TCanvas* _c3 = nullptr;

};

