#include "struct.h"

#ifndef TRACK3D_H
#define TRACK3D_H

struct track3D
{
  int tid;
  double px;
  double py;
  double pz;
  double yc;
  double zc;
  double r;
  double a;
  double b;
  double h;
  double ysig;
  double x0;
  double y0;
  double z0;
  double t0;
  int ret_ln;
  double chi2_ln;
  int ret_cr;
  double chi2_cr;
  std::vector<digit> clX;
  std::vector<digit> clY;
};

#ifdef __MAKECINT__
#pragma link C++ class std::vector < track3D > +;
#pragma link C++ class track3D + ;
#endif
#endif