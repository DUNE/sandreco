/*Implementation of the SANDECALCellInfo class,
for storing geometric info of the SAND ECAL cells*/

#include "SANDENDCAPModInfo.h"

// Default constructor
SANDENDCAPModInfo::SANDENDCAPModInfo() {}

// Parametric constructor
SANDENDCAPModInfo::SANDENDCAPModInfo(int arg_id, TGeoNode* arg_mod_node, const TGeoHMatrix &arg_hmatrix 
                                    /*double arg_x, double arg_y,
                                    double arg_z, double arg_width ,
                                    Orient arg_orientation */)
    : id_(arg_id),
      mod_node_(arg_mod_node),
      mod_hmatrix_(arg_hmatrix)
// orientation_(arg_orientation)
{
  // compute and fill the module translations
  x_ = mod_hmatrix_.GetTranslation()[0];
  y_ = mod_hmatrix_.GetTranslation()[1];
  z_ = mod_hmatrix_.GetTranslation()[2];
  // set the module width
  width_ = 2 * ((TGeoBBox*)mod_node_->GetVolume()->GetShape())->GetDX();
  // set the Al_thick
  get_Al_thick();
  //set the node path
  path_ = gGeoManager->GetPath();
}

// Setter methods for the attributes
void SANDENDCAPModInfo::id(int arg_id) { id_ = arg_id; }
void SANDENDCAPModInfo::x(double arg_x) { x_ = arg_x; }
void SANDENDCAPModInfo::y(double arg_y) { y_ = arg_y; }
void SANDENDCAPModInfo::z(double arg_z) { z_ = arg_z; }
void SANDENDCAPModInfo::width(double arg_width) { width_ = arg_width; }
// private setter that computes the Al layer thickness from a module section
void SANDENDCAPModInfo::Al_thick(double arg_Al_thick){
  Al_thick_ = arg_Al_thick;
}
void SANDENDCAPModInfo::get_Al_thick()
{
  Al_thick_ = 0;
  auto temp_d_node = mod_node_->GetDaughter(0);
  for (int i = 0; i < temp_d_node->GetNdaughters(); i++) {
    if (((TString)temp_d_node->GetDaughter(i)->GetName()).Contains("Alplate"))
      Al_thick_ =
          2 * ((TGeoBBox*)temp_d_node->GetDaughter(i)->GetVolume()->GetShape())
                  ->GetDZ();
  }
}
// void SANDENDCAPModInfo::orientation(Orient arg_orientation)
// {
//   orientation_ = arg_orientation;
// }

// Getter methods for the attributes
int SANDENDCAPModInfo::id() { return id_; }
double SANDENDCAPModInfo::x() { return x_; }
double SANDENDCAPModInfo::y() { return y_; }
double SANDENDCAPModInfo::z() { return z_; }
double SANDENDCAPModInfo::width() { return width_; }
double SANDENDCAPModInfo::Al_thick() { return Al_thick_; }
TString SANDENDCAPModInfo::path() { return path_; }
TGeoNode* SANDENDCAPModInfo::mod_node() { return mod_node_; }
TGeoHMatrix SANDENDCAPModInfo::mod_hmatrix() { return mod_hmatrix_; }
// SANDENDCAPModInfo::Orient SANDENDCAPModInfo::orientation()
// {
//   return orientation_;
// }

// PMT local_id and path length computation (could be split into two either)
// this should be equivalent to the older get_ecal_endcap_cell_local_id function
// void get_ecal_endcap_cell_local_id(double x, double y, double z, const TGeoNode* const slab_node,int& cell_local_id/* ,
//                           std::pair<double, double> path_lens */) const
// {
//   // convert (x,y,z) to local module coordinates
//   double master[3];
//   double local[3];
//   master[0] = x;
//   master[1] = y;
//   master[2] = z;

//   mod_node_->MasterToLocal(master, local);

//   for (int i = 0; i < mod_node_->GetNdaughters(); i++) {
//     auto sec_node = mod_node_->GetDaughter(i);
//     // this won't work: both curved and horizontal sections have the same name
//     if (slab_node->GetPath().Contains(sec_node->GetName())) {
//     }
//   }
// }