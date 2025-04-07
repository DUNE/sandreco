/*Implementation of the SANDENDCAPModInfo class,
for storing geometric info of the composite endcap modules*/

#include "SANDENDCAPModInfo.h"

namespace sand_geometry
{

namespace ecal
{

// Default constructor
ENDCAPModInfo::ENDCAPModInfo() {}

// Parametric constructor
ENDCAPModInfo::ENDCAPModInfo(int arg_id, TGeoNode* arg_mod_node, const TGeoHMatrix &arg_hmatrix 
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
  computeAlDz();
  // compute the min and max path lengths (and intialize other parameters)
  computeMinMaxL();
  // set the node path
  path_ = gGeoManager->GetPath();
}

// Setter methods for the attributes
void ENDCAPModInfo::setId(int arg_id) { id_ = arg_id; }
void ENDCAPModInfo::setX(double arg_x) { x_ = arg_x; }
void ENDCAPModInfo::setY(double arg_y) { y_ = arg_y; }
void ENDCAPModInfo::setZ(double arg_z) { z_ = arg_z; }
void ENDCAPModInfo::setWidth(double arg_width) { width_ = arg_width; }
// private setter that computes the Al layer dzness from a module section
void ENDCAPModInfo::setAlDz(double arg_al_dz) { al_dz_ = arg_al_dz; }
void ENDCAPModInfo::computeAlDz()
{
  al_dz_ = 0;
  auto temp_d_node = mod_node_->GetDaughter(0);
  for (int i = 0; i < temp_d_node->GetNdaughters(); i++) {
    if (((TString)temp_d_node->GetDaughter(i)->GetName()).Contains("Alplate"))
      al_dz_ = ((TGeoBBox*)temp_d_node->GetDaughter(i)->GetVolume()->GetShape())
                   ->GetDZ();
  }
}

void ENDCAPModInfo::computeMinMaxL()
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
          2 * al_dz_;
      break;
    }
  }

  // compute lmax_ and lmin_ based on the number of sections
  lmax_ =
      l_vert_ + 2 * (0.5 * M_PI * r_max_) + ((n_sec_ == 5) ? 2 : 1) * l_hor_;
  lmin_ =
      l_vert_ + 2 * (0.5 * M_PI * r_min_) + ((n_sec_ == 5) ? 2 : 1) * l_hor_;
}

// Getter methods for the attributes
int ENDCAPModInfo::getId() const { return id_; }
int ENDCAPModInfo::getNSections() const { return n_sec_; }
double ENDCAPModInfo::getX() const { return x_; }
double ENDCAPModInfo::getY() const { return y_; }
double ENDCAPModInfo::getZ() const { return z_; }
double ENDCAPModInfo::getWidth() const { return width_; }
double ENDCAPModInfo::getModDz() const { return mod_dz_; }
double ENDCAPModInfo::getLHor() const { return l_hor_; }
double ENDCAPModInfo::getLVert() const { return l_vert_; }
double ENDCAPModInfo::getRMin() const { return r_min_; }
double ENDCAPModInfo::getRMax() const { return r_max_; }
double ENDCAPModInfo::getAlDz() const { return al_dz_; }
TString ENDCAPModInfo::getPath() const { return path_; }
TGeoNode* ENDCAPModInfo::getModNode() const { return mod_node_; }
TGeoHMatrix ENDCAPModInfo::getModHMatrix() const { return mod_hmatrix_; }

double ENDCAPModInfo::getCurvatureArcLength(double depth) const
{
  return 0.5 * M_PI * (r_max_ - depth);
}
// compute the total cell path length given the depth along the module (w.r.t.
// the inner layer)
double ENDCAPModInfo::getCellTotalLength(double depth) const
{
  return lmax_ - M_PI * depth;
}
} // namespace ecal
} // namesoace sand_geometry