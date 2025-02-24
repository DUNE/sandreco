#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoTube.h>
#include <TMath.h>
#include <TObject.h>

#include <iostream>
#include <map>

#ifndef ENDCAPModInfo_H
#define ENDCAPModInfo_H

namespace sand_geometry
{

namespace ecal
{

// class for storing geometric info of the SAND ECAL cells
class ENDCAPModInfo : public TObject
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
  double al_dz_;  // half_thickness of the Al plate
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
  void computeAlDz();
  void computeMinMaxL();

 public:
  ENDCAPModInfo();  // Default constructor
  // ENDCAPModInfo(int id, double x, double y, double z, double length,
  //                   Orient orientation);  // parametric constructor
  ENDCAPModInfo(int arg_id, TGeoNode* arg_mod_node,
                    const TGeoHMatrix& arg_hmatrix);
  // Setter methods for the attributes
  void setId(int arg_id);
  void setX(double arg_x);
  void setY(double arg_y);
  void setZ(double arg_z);
  void setWidth(double arg_width);
  void setAlDz(double arg_al_thick);
  // void orientation(Orient arg_orientation);
  // Getter methods for the attributes
  int getId() const;
  int getNSections() const;
  double getX() const;
  double getY() const;
  double getZ() const;
  double getWidth() const;
  double getModDz() const;
  double getLHor() const;
  double getLVert() const;
  double getRMin() const;
  double getRMax() const;
  double getAlDz() const;
  TString getPath() const;
  TGeoNode* getModNode() const;
  TGeoHMatrix getModHMatrix() const;
  // PMT pos. and path length computation
  // void get_ecal_endcap_cell_local_id(double x, double y, double z,
  //                                    int& cell_id) const;+
  
  // compute the arc length along the curved sections at a given depth
  double getCurvatureArcLength(double depth) const;
  // compute the total cell path length given the depth along the module (w.r.t. the outer layer)
  double getCellTotalLength(double depth) const;

  ClassDef(ENDCAPModInfo, 1);
};
} // namespace ecal
} // namesoace sand_geometry

#ifdef __MAKECINT__
#pragma link C++ class sand_geometry::ecal::ENDCAPModInfo + ;
#endif

#endif