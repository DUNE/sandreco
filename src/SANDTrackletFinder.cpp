#include "SANDTrackletFinder.h"

double MinimizingFunction(const double* params, const sand_reco::tracker::Cluster& cluster, const std::map<sand_reco::tracker::DigitID, double>& digitId_to_drift_time)
{
  double dx = cos(params[2]);
  double dy = sin(params[3]);
  double dz = sin(params[2]);
  TVector3 pos(params[0], params[1], cluster.getZ());
  TVector3 dir(dx, dy, dz);
  dir = dir * (1. / dir.Mag());

  double sum = 0.0;
  auto sand_geo = cluster.getSandGeoManager();
  for (const auto& digitId_and_time : digitId_to_drift_time) {
    auto cell = sand_geo->getCellInfo(sand_geometry::tracker::CellID(digitId_and_time.first()))->second;
    TVector3 n = cell.getWire().getDirection().Cross(dir);
    double d = fabs(n.Dot(cell.getWire().getCenter() - pos)) / n.Mag();
    sum += (d - cell.getDriftVelocity() * digitId_and_time.second) 
            * (d - cell.getDriftVelocity() * digitId_and_time.second); 
  }
  return sum;
}

bool TrackletFinder::checkParallel(TVector3 d1, TVector3 d2)
{
  return d1.Dot(d2) == 1;
}

void TrackletFinder::linesParallelToWire(Line3D w, double distance, std::vector<Line3D>& lines)
{
    double perp_x_dir = -w.getDirection().Y();
    double perp_y_dir =  w.getDirection().X();

    double x_perp = w.getPoint().X() + perp_x_dir * distance / sqrt(perp_x_dir*perp_x_dir + perp_y_dir*perp_y_dir);
    double y_perp = w.getPoint().Y() + perp_y_dir * distance / sqrt(perp_x_dir*perp_x_dir + perp_y_dir*perp_y_dir);
    TVector3 direction_perp(perp_x_dir, perp_y_dir, 0);
    TVector3 point_perp(x_perp, y_perp, 0);

    lines.push_back(Line3D(point_perp, w.getDirection()));
    
    perp_x_dir =  w.getDirection().Y();
    perp_y_dir = -w.getDirection().X();

    x_perp = w.getPoint().X() + perp_x_dir * distance / sqrt(perp_x_dir*perp_x_dir + perp_y_dir*perp_y_dir);
    y_perp = w.getPoint().Y() + perp_y_dir * distance / sqrt(perp_x_dir*perp_x_dir + perp_y_dir*perp_y_dir);
    direction_perp = TVector3(perp_x_dir, perp_y_dir, 0);
    TVector3 point_perp2(x_perp, y_perp, 0);

    // point_perp2.Print();
    // w.getDirection().Print();
    lines.push_back(Line3D(point_perp2, w.getDirection()));
}

void TrackletFinder::computeCellsIntersections()
{
  std::vector<Line3D> lines;

  for (auto digitID:cluster_.getDigits()) {
    // std::cout << digitID() << std::endl;
    auto cell = cluster_.getSandGeoManager()->getCellInfo(sand_geometry::tracker::CellID(digitID()))->second;
    auto cell_size = cell.getSize();
    double h = cell_size.h;
    double w = cell_size.w;
    Line3D line_from_wire(cell.getWire().getCenter(), cell.getWire().getDirection());
    linesParallelToWire(line_from_wire, h, lines);
  }
  
  for (uint i = 0; i < lines.size() - 1; i++) {
    for (uint j = i + 1; j < lines.size(); j++) {
      // std::cout << i << " " << j << std::endl;
      bool is_parallel = checkParallel(lines[i].getDirection(), lines[j].getDirection());
      if (!is_parallel) {
        double t = (lines[i].getPoint() - lines[j].getPoint()).Dot(lines[i].getDirection() - lines[j].getDirection() * lines[i].getDirection().Dot(lines[j].getDirection())) /
                     (pow(lines[i].getDirection().Dot(lines[j].getDirection()), 2) - 1);
        
        TVector3 intersection = lines[i].getPoint() + t * lines[i].getDirection();
        cells_intersections_.push_back(intersection);
        // intersection.Print();
      }
    }
  }
}

void TrackletFinder::getScanningAreaVertices()
{
  double x_min = 1E9;
  double y_min = 1E9;
  double x_max = -1E9;
  double y_max = -1E9;
  for (auto p:cells_intersections_) {
    if (p.X() < x_min) {
      x_min = p.X();
    }
    if (p.Y() < y_min) {
      y_min = p.Y();
    }
    if (p.X() > x_max) {
      x_max = p.X();
    }
    if (p.Y() > y_max) {
      y_max = p.Y();
    }
  }
  const auto plane_half_dimension = cluster_.getPlane()->getDimension() * 0.5;
  const auto plane_position  = cluster_.getPlane()->getPosition();

  if (x_min <  plane_position.X() - plane_half_dimension.X()) {
    x_min = plane_position.X() - plane_half_dimension.X();
  }
  if (y_min <  plane_position.Y() - plane_half_dimension.Y()) {
    y_min = plane_position.Y() - plane_half_dimension.Y();
  }
  if (x_max >  plane_position.X() + plane_half_dimension.X()) {
    x_max = plane_position.X()  + plane_half_dimension.X();
  }
  if (y_max >  plane_position.Y() + plane_half_dimension.Y()) {
    y_max = plane_position.Y() + plane_half_dimension.Y();
  }

  cells_intersections_.clear();
  cells_intersections_.push_back(TVector3(x_min, y_min, 0));
  cells_intersections_.push_back(TVector3(x_max, y_max, 0));
}

void TrackletFinder::computeDriftTime()
{
  TVector2 mean_point_2d(0, 0);
  for (const auto& point:cells_intersections_) {
    mean_point_2d = mean_point_2d + TVector2(point.X(), point.Y()); 
  }
  mean_point_2d = mean_point_2d * ( 1. / cells_intersections_.size());

  for (const auto& digit_id:cluster_.getDigits()) {
    auto cell = cluster_.getSandGeoManager()->getCellInfo(sand_geometry::tracker::CellID(digit_id()))->second;
    
    // To Do: this part is used multiple times. It's point-line distance.
    //        Write it once in sandgeomangaer or utils.

    TVector3 leftend = cell.getWire().getReadoutPoint();
    TVector3 rightend = cell.getWire().getOppositePointToReadout();

    TVector3 r = cell.getWire().getDirection();
    mean_point_3d_ = TVector3(mean_point_2d.X(), 
                              mean_point_2d.Y(),
                              cell.getWire().getCenter().Z());

    TVector3 AP = mean_point_3d_ - leftend;
    double t = AP.Dot(r) / r.Mag2();
    t = std::max(0.0, std::min(1.0, t));

    TVector3 closest_point = leftend + t * r;
    double wire_time = (closest_point - leftend).Mag() / sand_reco::stt::v_signal_inwire;
    digitId_to_drift_time_[digit_id] = digit_collection_->getDigit(digit_id).tdc - digit_collection_->getDigit(digit_id).t_hit - wire_time;
  }
}

// To Do: find a minimizier able to escape local minima
std::vector<TVectorD> TrackletFinder::findTracklets()
{
  std::vector<TVectorD> minima;

  auto minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "");
  minimizer->SetMaxFunctionCalls(100000); 
  minimizer->SetMaxIterations(100000);
  minimizer->SetTolerance(1e-9);

  computeCellsIntersections();
  computeDriftTime();

  if (cells_intersections_.size() == 0) {
    return minima;
  }

  getScanningAreaVertices();
  setTrajectory(mean_point_3d_, TVector3(0,0,1));

  auto digitId_to_drift_time = digitId_to_drift_time_;
  auto cluster = cluster_;
  ROOT::Math::Functor functor_cells([&digitId_to_drift_time, &cluster](const double* params) { return MinimizingFunction(params, cluster, digitId_to_drift_time); }, 4);
  
  minimizer->SetFunction(functor_cells);
    
  double theta_xz = atan(trajectory_.getDirection().Z() / trajectory_.getDirection().X());
  double theta_yz = atan(trajectory_.getDirection().Y() / trajectory_.getDirection().Z());
  double theta_width = M_PI_4;

  // std::cout << cells_intersections_[1].X() << " " << cells_intersections_[0].X() << std::endl;
  // std::cout << cells_intersections_[1].Y() << " " << cells_intersections_[0].Y() << std::endl;

  double x_width = cells_intersections_[1].X() - cells_intersections_[0].X();
  double y_width = cells_intersections_[1].Y() - cells_intersections_[0].Y();
  
  // To Do: should be a config parameter
  int subdivisions = 2;
  double x_sub_width = x_width / subdivisions;
  double y_sub_width = y_width / subdivisions;
  double theta_sub_width = theta_width / subdivisions;

  std::vector<TVectorD> sampling_points;
  for (int i = 0; i <= subdivisions; i++) {
    for (int j = 0; j <= subdivisions; j++) {
      for (int k = 0; k <= subdivisions; k++) {
        for (int z = 0; z <= subdivisions; z++) {
          TVectorD point(4);

          point[0] = cells_intersections_[0].X() + i * x_sub_width;
          point[1] = cells_intersections_[0].Y() + j * y_sub_width;
          point[2] = theta_xz - theta_width + k * theta_sub_width;
          point[3] = theta_yz - theta_width + z * theta_sub_width;

          sampling_points.push_back(point);
        }
      }
    }  
  }
  // std::cout << "SAMPLING SIZE: " << sampling_points.size() << std::endl;
  for (uint i = 0; i < sampling_points.size(); i++) {
    double starting_point[4] = {sampling_points[i][0], sampling_points[i][1], sampling_points[i][2], sampling_points[i][3]};
    minimizer->SetLimitedVariable(0, "px", starting_point[0], 1, cells_intersections_[0].X() - 200, cells_intersections_[1].X() + 200);
    minimizer->SetLimitedVariable(1, "py", starting_point[1], 1, cells_intersections_[0].Y() - 200, cells_intersections_[1].Y() + 200);
    minimizer->SetLimitedVariable(2, "dx", starting_point[2], 0.1, theta_xz - theta_width, theta_xz + theta_width);
    minimizer->SetLimitedVariable(3, "dy", starting_point[3], 0.1, theta_yz - theta_width, theta_yz + theta_width);

    minimizer->Minimize();
    const double *xs = minimizer->X();
    TVectorD min(5);
    min[0] = xs[0];
    min[1] = xs[1];
    min[2] = xs[2];
    min[3] = xs[3];
    min[4] = minimizer->MinValue();
    minima.push_back(min);
  }
  return minima;
}

void TrackletFinder::clear()
{
  digit_collection_ = nullptr;
  cells_intersections_.clear();
  digitId_to_drift_time_.clear();
}

void TrackletFinder::draw3DWires() {
  gStyle->SetOptStat(0);
  c3_ = new TCanvas("c3D","c3D",1500,1500);
  auto plane_half_dimension = cluster_.getPlane()->getDimension() * 0.5;
  auto plane_position  = cluster_.getPlane()->getPosition();
  TH3D h("", "", 2 * plane_half_dimension.X(), plane_position.X() - plane_half_dimension.X(), plane_position.X() + plane_half_dimension.X(),
                 2 * plane_half_dimension.Y(), plane_position.Y() - plane_half_dimension.Y(), plane_position.Y() + plane_half_dimension.Y(),
                 2 * plane_half_dimension.Z(), plane_position.Z() - plane_half_dimension.Z(), plane_position.Z() + plane_half_dimension.Z());
  h.SetTitle(";x;y;z");
  h.Draw();

  for (auto digitID:cluster_.getDigits()) {
    auto cell = cluster_.getSandGeoManager()->getCellInfo(sand_geometry::tracker::CellID(digitID()))->second;
    TPolyLine3D* pl2 = new TPolyLine3D(2);
    TVector3 l_start = cell.getWire().getFirstPoint();
    pl2->SetPoint(0, l_start.X(), l_start.Y(), l_start.Z());
    TVector3 l_end = cell.getWire().getSecondPoint();
    pl2->SetPoint(1, l_end.X(), l_end.Y(), l_end.Z());
    pl2->Draw("same");
  }

  c3_->SaveAs("./c3D.C");
}

TVector3 getCylinderCoordinates(Line3D w, double t, double radius, double theta) {
  double x = w.getPoint().X() + t*w.getDirection().X() + radius*(cos(theta)*w.getU().X() + sin(theta)*w.getV().X());
  double y = w.getPoint().Y() + t*w.getDirection().Y() + radius*(cos(theta)*w.getU().Y() + sin(theta)*w.getV().Y());
  double z = w.getPoint().Z() + t*w.getDirection().Z() + radius*(cos(theta)*w.getU().Z() + sin(theta)*w.getV().Z());

  return TVector3(x, y, z);
}

void TrackletFinder::draw3D()
{
  gStyle->SetOptStat(0);
  if (!c3_) {
    c3_ = new TCanvas("c3D","c3D",1500,1500);
  }

  auto sand_geo = cluster_.getSandGeoManager();

  int ccc = 3;
  for (const auto& digitId_and_time : digitId_to_drift_time_) {
    auto cell = sand_geo->getCellInfo(sand_geometry::tracker::CellID(digitId_and_time.first()))->second;
    
  auto plane_half_dimension = cluster_.getPlane()->getDimension() * 0.5;
  auto plane_position  = cluster_.getPlane()->getPosition();
  TH3D* h3 = new TH3D("", "", 2 * plane_half_dimension.X(), plane_position.X() - plane_half_dimension.X(), plane_position.X() + plane_half_dimension.X(),
                              2 * plane_half_dimension.Y(), plane_position.Y() - plane_half_dimension.Y(), plane_position.Y() + plane_half_dimension.Y(),
                              2 * plane_half_dimension.Z(), plane_position.Z() - plane_half_dimension.Z(), plane_position.Z() + plane_half_dimension.Z());
    h3->SetTitle(";x;y;z");
    h3->SetMarkerColor(ccc);
    Line3D line_from_wire(cell.getWire().getCenter(), cell.getWire().getDirection());
    for (int t = -100; t < 100; t++) {
      for (int theta = 0; theta < 628; theta++) {
        TVector3 p = getCylinderCoordinates(line_from_wire, (double)t, cell.getDriftVelocity() * digitId_and_time.second, theta/100.);
        h3->Fill(p.X(), p.Y(), p.Z());
      }
    }
    h3->Draw("same scat");
    ccc++;
  }

  c3_->SaveAs("./c3D.png");
  
}

void TrackletFinder::draw2DWires()
{
  gStyle->SetOptStat(0);
  if (!c2_) {
    c2_ = new TCanvas("c2D","c2D",1500,1500);
    TH2D* h2 = new TH2D("h","h", cells_intersections_[1].X()  - cells_intersections_[0].X(), cells_intersections_[0].X(), cells_intersections_[1].X(),
                                 cells_intersections_[1].Y()  - cells_intersections_[0].Y(), cells_intersections_[0].Y(), cells_intersections_[1].Y());

    h2->Draw();
  }

  // for (auto l:_cells_bands) {
  //   TVector3 second_point;
  //   second_point = l.getPoint() + 100*l.getDirection();
  //   TLine* tl = new TLine(l.getPoint().X(), l.getPoint().Y(), second_point.X(), second_point.Y());
  //   tl->SetLineStyle(2);
  //   tl->Draw("same");
  // }
  
  for (auto digitID:cluster_.getDigits()) {
    auto cell = cluster_.getSandGeoManager()->getCellInfo(sand_geometry::tracker::CellID(digitID()))->second;
    TVector3 l_start = cell.getWire().getFirstPoint();
    TVector3 l_end   = cell.getWire().getSecondPoint();

    TLine* tl = new TLine(l_start.X(), l_start.Y(), l_end.X(), l_end.Y());
    tl->Draw("same");
  }

  c2_->SaveAs("./c2D.png");

}

void TrackletFinder::draw2DDistance()
{
  gStyle->SetOptStat(0);
  TCanvas c("c2DMinimization","c2DMinimization",1500,1500);
  if (!c2_) {
    c2_ = new TCanvas("c2D","c2D",1500,1500);

  }
  
  TH2D* h2 = new TH2D("h","h", cells_intersections_[1].X()  - cells_intersections_[0].X(), 
                               cells_intersections_[0].X(), cells_intersections_[1].X(),
                               cells_intersections_[1].Y()  - cells_intersections_[0].Y(), 
                               cells_intersections_[0].Y(), cells_intersections_[1].Y());

  double min;
  auto digitId_to_drift_time = digitId_to_drift_time_;
  auto sand_geo = cluster_.getSandGeoManager();
  
  double theta_xz = atan(trajectory_.getDirection().Z() / trajectory_.getDirection().X()) * 1000;
  double theta_yz = atan(trajectory_.getDirection().Y() / trajectory_.getDirection().Z()) * 1000;
  int count = 0;
  for (int px = cells_intersections_[0].X(); px <= cells_intersections_[1].X(); px++) {
    for (int py = cells_intersections_[0].Y(); py <= cells_intersections_[1].Y(); py++) {
      min = 1E9;
      for (int angle_xz = theta_xz - 400; angle_xz < theta_xz + 400 ; angle_xz+=10) {
        for (int angle_yz = theta_yz - 400; angle_yz < theta_yz + 400; angle_yz+=10) {
          double p[4] = {px / 1., py / 1., angle_xz / 1000., angle_yz / 1000.};
          double tmp_min = MinimizingFunction(p, cluster_, digitId_to_drift_time);
          if (tmp_min < min) {
            min = tmp_min;
          }
        }
      }
      h2->Fill(px, py, min);
      count++;
    }
  }

  h2->GetZaxis()->SetRangeUser(0, 10e-1);
  h2->Draw("colz");

  c2_->SaveAs("./c2D.png");

}