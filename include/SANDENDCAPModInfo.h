#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoTube.h>
#include <TMath.h>
#include <TObject.h>

#include <iostream>
#include <map>

#ifndef SANDENDCAPModInfo_H
#define SANDENDCAPModInfo_H

// class for storing geometric info of the SAND ECAL cells
class SANDENDCAPModInfo : public TObject
{
 public:
  enum class Orient { kHorizontal, kVertical };

 private:
  int id_;    // id of the cell
  double x_;  // x position of the center of the module
  double y_;  // y position of the center of the module
  double z_;  // z position of the center of the module
  double l_hor_;
  double l_vert_;
  double r_max_;
  double r_min_;
  int n_sec_;     // number of sections in the module (either 4 or 5)
  double lmin_;   // min (outer) cell path length
  double lmax_;   // max (inner) cell path length
  double width_;  // width of the module (locally along the x direction)
  double mod_dz_; // total depth of the module (including Al thickness)
  double Al_dz_;  // half_thickness of the Al plate
  TString path_;
  TGeoNode* mod_node_;       // corresponding node of the module
  TGeoHMatrix mod_hmatrix_;  // h_matrix of the module
  // std::map<int, TGeoNode> cell_sections_;  // map of the sections within the
  //                                          // module cell (indicized by the
  //                                          // section id) (could be something
  //                                          // else apart from TGeoNode)

  // double length_;       // length of the cell?
  // Orient orientation_;  // orientation of the cell?
  // private Al_thick setter from mod_node_
  void get_Al_dz();
  void compute_min_max_l();

 public:
  SANDENDCAPModInfo();  // Default constructor
  // SANDENDCAPModInfo(int id, double x, double y, double z, double length,
  //                   Orient orientation);  // parametric constructor
  SANDENDCAPModInfo(int arg_id, TGeoNode* arg_mod_node,
                    const TGeoHMatrix& arg_hmatrix);
  // Setter methods for the attributes
  void id(int arg_id);
  void x(double arg_x);
  void y(double arg_y);
  void z(double arg_z);
  void width(double arg_width);
  void Al_dz(double arg_Al_thick);
  // void orientation(Orient arg_orientation);
  // Getter methods for the attributes
  int id() const;
  double x() const;
  double y() const;
  double z() const;
  double width() const;
  double mod_dz() const;
  double l_hor() const;
  double l_vert() const;
  double rmin() const;
  double rmax() const;
  double Al_dz() const;
  TString path() const;
  TGeoNode* mod_node() const;
  TGeoHMatrix mod_hmatrix() const;
  // PMT pos. and path length computation
  // void get_ecal_endcap_cell_local_id(double x, double y, double z,
  //                                    int& cell_id) const;+
  
  // compute the arc length along the curved sections at a given depth
  double get_curv_arc_len(double depth) const;
  // compute the total cell path length given the depth along the module (w.r.t. the outer layer)
  double get_cell_tot_len(double depth) const;

  ClassDef(SANDENDCAPModInfo, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class SANDENDCAPModInfo + ;
#endif

#endif
  