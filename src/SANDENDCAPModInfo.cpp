/*Implementation of the SANDENDCAPModInfo class,
for storing geometric info of the composite endcap modules*/

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
  // set the number of sections
  n_sec_ = mod_node_->GetNdaughters();
  // set the module width
  width_ = 2 * ((TGeoBBox*)mod_node_->GetVolume()->GetShape())->GetDX();
  // set the Al_dz
  get_Al_dz();
  // compute the min and max path lengths (and intialize other parameters)
  compute_min_max_l();
  // set the node path
  path_ = gGeoManager->GetPath();

  // if (mod_node_->GetNdaughters() > 4) {
  //   std::cout << "hor: ("
  //             << (mod_hmatrix_ * (*mod_node_->GetDaughter(4)->GetMatrix()))
  //                    .GetTranslation()[0]
  //             << ", "
  //             << (mod_hmatrix_ * (*mod_node_->GetDaughter(4)->GetMatrix()))
  //                    .GetTranslation()[1]
  //             << ", "
  //             << (mod_hmatrix_ * (*mod_node_->GetDaughter(4)->GetMatrix()))
  //                    .GetTranslation()[2]
  //             << ")\n";
  // }
}

// Setter methods for the attributes
void SANDENDCAPModInfo::id(int arg_id) { id_ = arg_id; }
void SANDENDCAPModInfo::x(double arg_x) { x_ = arg_x; }
void SANDENDCAPModInfo::y(double arg_y) { y_ = arg_y; }
void SANDENDCAPModInfo::z(double arg_z) { z_ = arg_z; }
void SANDENDCAPModInfo::width(double arg_width) { width_ = arg_width; }
// private setter that computes the Al layer dzness from a module section
void SANDENDCAPModInfo::Al_dz(double arg_Al_dz) { Al_dz_ = arg_Al_dz; }
void SANDENDCAPModInfo::get_Al_dz()
{
  Al_dz_ = 0;
  auto temp_d_node = mod_node_->GetDaughter(0);
  for (int i = 0; i < temp_d_node->GetNdaughters(); i++) {
    if (((TString)temp_d_node->GetDaughter(i)->GetName()).Contains("Alplate"))
      Al_dz_ = ((TGeoBBox*)temp_d_node->GetDaughter(i)->GetVolume()->GetShape())
                   ->GetDZ();
  }
}

void SANDENDCAPModInfo::compute_min_max_l()
{
  for (int i = 0; i < n_sec_; i++) {
    if (((TString)mod_node_->GetDaughter(i)->GetName()).Contains("vert")) {
      l_vert_ =
          2 * ((TGeoBBox*)mod_node_->GetDaughter(i)->GetVolume()->GetShape())
                  ->GetDY();
      mod_dz_ = ((TGeoBBox*)mod_node_->GetDaughter(i)->GetVolume()->GetShape())
                    ->GetDZ();
      break;
    }
  }
  for (int i = 0; i < n_sec_; i++) {
    if (((TString)mod_node_->GetDaughter(i)->GetName()).Contains("hor")) {
      l_hor_ =
          2 * ((TGeoBBox*)mod_node_->GetDaughter(i)->GetVolume()->GetShape())
                  ->GetDY();
      break;
    }
  }
  for (int i = 0; i < n_sec_; i++) {
    if (((TString)mod_node_->GetDaughter(i)->GetName()).Contains("curv")) {
      r_max_ =
          ((TGeoTubeSeg*)mod_node_->GetDaughter(i)->GetVolume()->GetShape())
              ->GetRmax();
      // correct the min radius for the Al plate thickness
      r_min_ =
          ((TGeoTubeSeg*)mod_node_->GetDaughter(i)->GetVolume()->GetShape())
              ->GetRmin() +
          2 * Al_dz_;
      break;
    }
  }

  // compute lmax_ and lmin_ based on the number of sections
  lmax_ =
      l_vert_ + 2 * (0.5 * M_PI * r_max_) + ((n_sec_ == 5) ? 2 : 1) * l_hor_;
  lmin_ =
      l_vert_ + 2 * (0.5 * M_PI * r_min_) + ((n_sec_ == 5) ? 2 : 1) * l_hor_;
}

// void SANDENDCAPModInfo::orientation(Orient arg_orientation)
// {
//   orientation_ = arg_orientation;
// }

// Getter methods for the attributes
int SANDENDCAPModInfo::id() const { return id_; }
int SANDENDCAPModInfo::n_sections() const { return n_sec_; }
double SANDENDCAPModInfo::x() const { return x_; }
double SANDENDCAPModInfo::y() const { return y_; }
double SANDENDCAPModInfo::z() const { return z_; }
double SANDENDCAPModInfo::width() const { return width_; }
double SANDENDCAPModInfo::mod_dz() const { return mod_dz_; }
double SANDENDCAPModInfo::l_hor() const { return l_hor_; }
double SANDENDCAPModInfo::l_vert() const { return l_vert_; }
double SANDENDCAPModInfo::rmin() const { return r_min_; }
double SANDENDCAPModInfo::rmax() const { return r_max_; }
double SANDENDCAPModInfo::Al_dz() const { return Al_dz_; }
TString SANDENDCAPModInfo::path() const { return path_; }
TGeoNode* SANDENDCAPModInfo::mod_node() const { return mod_node_; }
TGeoHMatrix SANDENDCAPModInfo::mod_hmatrix() const { return mod_hmatrix_; }

double SANDENDCAPModInfo::get_curv_arc_len(double depth) const
{
  return 0.5 * M_PI * (r_max_ - depth);
}
// compute the total cell path length given the depth along the module (w.r.t.
// the inner layer)
double SANDENDCAPModInfo::get_cell_tot_len(double depth) const
{
  return lmax_ - M_PI * depth;
}
