#pragma once

#include "TMatrixD.h"
#include "TVector3.h"
#include "TVectorD.h"
#include <TDecompSVD.h>
#include "TH3D.h"
#include "TH2D.h"
#include <TDecompLU.h>
#include <TStyle.h>


#include <Math/Factory.h>
#include <Math/Minimizer.h>
#include <Math/Functor.h>
#include <TCanvas.h>
#include <TPolyLine3D.h>
#include <TLine.h>

#include <SANDTrackerCell.h>
#include <CLine3D.h>
#include "utils.h"

#include "iostream"
#include "math.h"

class TrackletFinder {
  public:
    TrackletFinder() {};
    ~TrackletFinder() {};

    void SetVolumeParameters(int* p) {_volume_parameters = p;};
    void SetSigmaPosition(double sp) {_sigma_pos = sp;};
    void SetSigmaAngle(double sa)    {_sigma_ang = sa;};
    void SetCells(std::map<double, SANDTrackerCell>* cells) {_fired_cells = cells;};
    void SetTrajectory(TVector3 tp, TVector3 td)       {_trajectory = CLine3D(tp, td);};

    bool CheckParallel(TVector3 d1, TVector3 d2);
    void LinesParallelToWire(CLine3D w, double distance, std::vector<CLine3D>& lines);
    void ComputeCellsIntersections();
    void ComputeCellsBands();

    std::vector<TVectorD> FindTracklets();

    void Clear();

    void Draw3D();
    void Draw3DWires();
    void Draw2DWires();
    void Draw2DDistance();

  private:
    const std::map<double, SANDTrackerCell>* _fired_cells = nullptr;
    CLine3D _trajectory;

    std::vector<TVector3> _cells_intersections;
    std::vector<CLine3D>  _cells_bands;

    double _sigma_pos; // mm
    double _sigma_ang; // rad

    int* _volume_parameters;

    TCanvas* _c2 = nullptr;
    TCanvas* _c3 = nullptr;

};

