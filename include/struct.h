#include <map>
#include <string>
#include <vector>
#include <cstdint>
#include <cstddef>
#include <TString.h>
#include <TGeoManager.h>

#ifndef STRUCT_H
#define STRUCT_H

struct pe
{
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
    double vis_e;
    int n_traj;
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
struct dg_ps
{
  int side;
  double adc;
  double tdc;
  std::vector<pe> photo_el;
};

struct dg_cell
{
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

/**
 * @struct reco_cell
 * @brief Represents a reconstructed detector cell with spatial, energy, and timing information.
 * 
 * This structure holds information about a detector cell, including its position, 
 * energy deposit, time, and associated photodetector signals.
 */
struct reco_cell {
  int id;      /**< Unique identifier of the cell */
  double z;    /**< Z-coordinate of the cell position */
  double y;    /**< Y-coordinate of the cell position */
  double x;    /**< X-coordinate of the cell position */
  double l;    /**< Reconstructed cell length */
  int mod;     /**< Module number where the cell is located */
  int lay;     /**< Layer number of the cell */
  double e;    /**< Energy deposited in the cell */
  double t;    /**< Time of the recorded signal */

  dg_ps ps1;   /**< Photodetector signal from the first side */
  dg_ps ps2;   /**< Photodetector signal from the second side */

  int fired_pmt; /**< Indicates which photodetectors were triggered:
                  *   - 1 if only ps1 is set
                  *   - 2 if only ps2 is set
                  *   - 3 if both ps1 and ps2 are set
                  */
};

struct dg_wire
{
  std::string det;
  long did;
  double x;
  double y;
  double z;
  double t0;
  double de;
  double adc;
  double tdc = 1e9;
  bool hor;
  double wire_length;
  std::vector<int> hindex;
  /*
    ADDENDUM
    tdc = drift_time + signal_time + t_hit
    added to check validity of track fitting
    reconstruction method for drift chamber
  */
  // true quantities
  double t_hit = 1e9;
  double signal_time = 1e9;
  double drift_time = 1e9;
  // measured quantities
  double t_hit_measured = 1e9;        // via global trigger
  double signal_time_measured = 1e9;  // exploit different wire orientation
  double drift_time_measured =
      1e9;  // tdc - signal_time_measured - t_hit_measured

  double missing_coordinate = 1e9;
};

struct cluster
{
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
  int type;   // type 1 barrel
              // type 2 endcap
              // type 3 mixed
  std::vector<reco_cell> reco_cells; 
};

struct track
{
  int tid = -1;
  double yc = NAN;
  double zc = NAN;
  double r = NAN;
  double a = NAN;
  double b = NAN;
  double h = NAN;
  double ysig = NAN;
  double x0 = NAN;
  double y0 = NAN;
  double z0 = NAN;
  double t0 = NAN;
  int ret_ln = -1;
  double chi2_ln = NAN;
  int ret_cr = -1;
  double chi2_cr = NAN;
  std::vector<dg_wire> clX;
  std::vector<dg_wire> clY;
};

struct track_grain_lens {
  double xstart, ystart, zstart;
  double xend, yend, zend;
  int track_ids;
  int PDG_reco;
	double energy_reco;
};

struct vertex_grain_lens {
  double xvtx, yvtx, zvtx;
};

struct particle
{
  int primary;
  int pdg;
  int pdg_true;
  int pdg_reco;
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

  //track_grain_lens trackReco; // traccia ricostruita dal TTree

  double xend= 0;
  double yend = 0;
  double zend = 0; 
  bool filled = false;
  bool isinGRAIN = false;

};

struct vertex 
{
  int id;
  double x;
  double y;
  double z;
  std::vector<int> track_ids;
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

struct volume
{
  TGeoVolume* geo_volume;
  TString volume_path;
  bool IsActive;
};

struct wf_grain
{
   std::vector<pe> photo_el;
   std::vector<float> samples;
   float t0;
};

// grain detector response

struct dg_grain
{
   uint16_t channel_id;
   double time_rising_edge;
   double time_over_threshold;
   double charge;
   std::vector<int> h_indices; 
};

using grain_sparse_image = std::vector<dg_grain>;

// grain spill slicer

struct grain_dense_image 
{
    std::vector<float> charge;
    std::vector<float> time;
    uint16_t camera_id;  
    std::vector<std::vector<int>> h_indices_img;   
};

// volumereco

struct voxel_grain
{
    size_t vox_dims_x;
    size_t vox_dims_y;
    size_t vox_dims_z;
    
    std::vector<float> voxels;
    float at(size_t x_idx, size_t y_idx, size_t z_idx) const
    {
        return voxels.at(z_idx + y_idx * vox_dims_z + x_idx * vox_dims_y * vox_dims_z); 
    }

    float& at(size_t x_idx, size_t y_idx, size_t z_idx)
    {
        return voxels.at(z_idx + y_idx * vox_dims_z + x_idx * vox_dims_y * vox_dims_z); 
    }

};

// volumereco analyis or lens analyis

struct tracklet_grain
{
    float x; 
    float y; 
    float z;
    float px;
    float py;
    float pz;
    float energy;
    float thickness;
    std::vector<int> h_indices;     
};

struct cluster_grain
{
    std::vector<tracklet_grain> tracklets;
    float cx;
    float cy;
    float cz;
};

struct lar_point {
  double x;
  double y;
  double z;
  double time;
  double energy;
  double xsize;
  double ysize;
  double zsize;
};

struct lar_track {
  int tid;
  std::string method;
  int visibility; 

  double xstart, ystart, zstart;
  double xmcstart, ymcstart, zmcstart;

  double xend, yend, zend;
  double xmcend, ymcend, zmcend;
 
  double recoE, trueE;
  double time;
 
  double xdir_zx, ydir_zx, zdir_zx;
  double xdir_zy, ydir_zy, zdir_zy;
 
  bool hasSTTdigits;
  int contained;
  int STTmatch, STTmatch2;
  
  std::vector<lar_point> points;
};

struct grain_event {
  int evnum;
  double xvtx, yvtx, zvtx;
  double xmcvtx, ymcvtx, zmcvtx;
  int n_ttrack;
  int n_rtrack;
  double layer;
  double tolerance;
  std::vector<lar_track> tracks;
};

#endif