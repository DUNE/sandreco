#include <TGeoManager.h>
#include<TGeoBBox.h>
#include <TGeoNode.h>
#include <TObject.h>
#include <map>

#ifndef SANDENDCAPModInfo_H
#define SANDENDCAPModInfo_H

// class for storing geometric info of the SAND ECAL cells
class SANDENDCAPModInfo : public TObject
{
 public:
  enum class Orient { kHorizontal, kVertical };

 private:
  int id_;              // id of the cell
  double x_;            // x position of the center of the module
  double y_;            // y position of the center of the module
  double z_;            // z position of the center of the module
  double width_;        // width of the module (locally along the x direction)
  double Al_thick_;      // thickness of the Al plate
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
  void get_Al_thick();

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
  void Al_thick(double arg_Al_thick);
  // void orientation(Orient arg_orientation);
  // Getter methods for the attributes
  int id();
  double x();
  double y();
  double z();
  double width();
  double Al_thick();
  TString path();
  TGeoNode* mod_node();
  TGeoHMatrix mod_hmatrix();
  // PMT pos. and path length computation
  void get_ecal_endcap_cell_local_id(double x, double y, double z,
                                     int& cell_id) const;

  ClassDef(SANDENDCAPModInfo, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class SANDENDCAPModInfo + ;
#endif

#endif
