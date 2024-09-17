#pragma once
#include "SANDTrackerCell.h"
#include <vector>
#include <map>

class SANDTrackerPlane {
private:
  long _id;
  double _rotation;
  TVector3 _center_position;
  TVector3 _dimensions;
  std::map<double, long> _coord_to_id_map;
  std::map<long, SANDTrackerCell> _id_to_wire_map;
  
  std::vector<SANDTrackerCell> _vCells;
  std::vector<long> _vFiredCellsIndex;

public:
  SANDTrackerPlane(long idPlane) { _id = idPlane; }
  SANDTrackerPlane(long idPlane, std::vector<SANDTrackerCell> cells) {
    _id = idPlane;
    _vCells = cells;
  }

  long id() const { return _id; }
  void addCell(const SANDTrackerCell cell) { _vCells.push_back(cell); }
  int nCells() const { return _vCells.size(); }
  SANDTrackerCell &getCell(int);

  std::vector<long> firedCellsIndex() const { return _vFiredCellsIndex; }
  void propagateTrack(const CLine3D track);
  void clear();

};
