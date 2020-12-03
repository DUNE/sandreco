// File struct.h
#include <vector>
#include <map>

#ifndef STRUCT_H
#define STRUCT_H

struct cell
{
  int id;
  double z;
  double y;
  double x;
  double l;
  double adc1;
  double tdc1;
  double adc2;
  double tdc2;
  int mod;
  int lay;
  int cel;
  std::vector<double> pe_time1;
  std::vector<int> hindex1;
  std::vector<double> pe_time2;
  std::vector<int> hindex2;
};

struct cluster
{
  int tid;
  double x;
  double y;
  double z;
  double t;
  double e;
  double sx;
  double sy;
  double sz;
  double varx;
  double vary;
  double varz;
  std::vector<cell> cells;
};

struct hit
{
  std::string det;
  double x1;
  double y1;
  double z1;
  double t1;
  double x2;
  double y2;
  double z2;
  double t2;
  double de;
  int pid;
  int index;
};

struct digit
{
  std::string det;
  double x;
  double y;
  double z;
  double t;
  double de;
  bool hor;
  std::vector<int> hindex;
};

struct track
{
  int tid;
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
  std::vector<digit> digits;
  //std::vector<digit> clX;
  //std::vector<digit> clY;
};

struct particle
{
  int primary;
  int pdg;
  int tid;
  int parent_tid;
  double charge;
  double mass;
  double pxtrue;
  double pytrue;
  double pztrue;
  double Etrue;
  double xtrue;
  double ytrue;
  double ztrue;
  double ttrue;

  double pxreco;
  double pyreco;
  double pzreco;
  double Ereco;
  double xreco;
  double yreco;
  double zreco;
  double treco;
  bool kalman_ok;
  bool has_track;
  double charge_reco;
  track tr;

  bool has_cluster;
  cluster cl;

  bool has_daughter;
  std::vector<particle> daughters;
};

struct event
{
  double x;
  double y;
  double z;
  double t;
  double Enu;
  double pxnu;
  double pynu;
  double pznu;
  double Enureco;
  double pxnureco;
  double pynureco;
  double pznureco;
  std::vector<particle> particles;
};

struct gcell
{
  int id;
  double Z[4];
  double Y[4];
  double adc;
  double tdc;
};

#ifdef __MAKECINT__
#pragma link C++ class std::map < int, std::vector < double >> +;
#pragma link C++ class std::map < int, std::vector < int >> +;
#pragma link C++ class std::map < int, double > +;
#pragma link C++ class std::vector < cell > +;
#pragma link C++ class std::map < std::string, std::vector < hit >> +;
#pragma link C++ class std::vector < digit > +;
#pragma link C++ class std::vector < track > +;
#pragma link C++ class std::vector < cluster > +;
#pragma link C++ class std::vector < particle > +;
#pragma link C++ class digit + ;
#pragma link C++ class cell + ;
#pragma link C++ class cluster + ;
#pragma link C++ class track + ;
#pragma link C++ class particle + ;
#pragma link C++ class event + ;
#endif
#endif
