// File struct.h
#include <map>
#include <string>
#include <vector>

#ifndef STRUCT_H
#define STRUCT_H

struct pe {
  double time;
  int h_index;
};
struct cluster_generator{
  int pdg_code; 
  int parent_pdg_code;
  int track_id;
  int parent_track_id; 
  double dep_energy;
  double initial_energy;
  double initial_momentum; 
  double initial_x; 
  double initial_y;
  double initial_z;  
};

struct truecluster{
    int tid;
    double x;
    double y;
    double z;
    double t; 
    double e;
    double sx;
    double sy;
    double sz;
    int ntot_cell;
    int cell_l0;
    int cell_l1; 
    int cell_l2;
    int cell_l3;
    int cell_l4;
    double energy_l0;
    double energy_l1;
    double energy_l2;
    double energy_l3;
    double energy_l4;
    double lay0_maxE;
    double lay1_maxE;
    double lay2_maxE;
    double lay3_maxE;
    double lay4_maxE;
    double asymmetry; 
    double Eoverp;
    bool moregens= false;
    std::vector<cluster_generator> vec_generator;
};

struct hit {
  std::string det;
  int did;
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

// photo-signal
struct dg_ps {
  int side;
  double adc;
  double tdc;
  std::vector<pe> photo_el;
};

struct dg_cell {
  int id;
  double z;
  double y;
  double x;
  double l;
  int mod;
  int lay;
  int cel;
  int det;
  std::vector<dg_ps> ps1;
  std::vector<dg_ps> ps2;
};

struct reco_cell {
  int id;
  double z;
  double y;
  double x;
  double l;
  int mod;
  int lay;
  double e;
  double t; 
  dg_ps ps1;
  dg_ps ps2;
  int fired_pmt;
};

struct dg_tube {
  std::string det;
  int did;
  double x;
  double y;
  double z;
  double t0;
  double de;
  double adc;
  double tdc;
  bool hor;
  std::vector<int> hindex;
};

struct cluster {
  int tid;
  double x;
  double y;
  double z;
  double t;
  double e;
  double ax;
  double ay;
  double az;
  double sx;
  double sy;
  double sz;
  double varx;
  double vary;
  double varz;
  int type;
  std::vector<reco_cell> reco_cells; 
};

struct track {
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
  std::vector<dg_tube> clX;
  std::vector<dg_tube> clY;
};

struct particle {
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

struct event {
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

struct gcell {
  int id;
  double Z[4];
  double Y[4];
  double adc;
  double tdc;
};

#endif
