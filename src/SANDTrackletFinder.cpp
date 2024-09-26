#include "SANDTrackletFinder.h"

double MinimizingFunction(const double* params, const  std::map<std::pair<double, double>, SANDTrackerCell>* cells)
{
  double dx = cos(params[2]);
  double dy = sin(params[3]);
  double dz = sin(params[2]);
  // TVector3 pos(params[0], params[1], (*cells)[0].wire().center().Z());
  TVector3 pos(params[0], params[1], 2391);
  TVector3 dir(dx, dy, dz);
  dir = dir * (1. / dir.Mag());

  double sum = 0.0;
  for (const auto& c : (*cells)) {
    TVector3 n = c.second.wire().getDirection().Cross(dir);
    double d = fabs(n.Dot(c.second.wire().center() - pos)) / n.Mag();
    sum += (d - c.second.driftVelocity() * c.first.first) 
            * (d - c.second.driftVelocity() * c.first.first); 
  }
  return sum;
}

bool TrackletFinder::CheckParallel(TVector3 d1, TVector3 d2)
{
  return d1.Dot(d2) == 1;
}

void TrackletFinder::LinesParallelToWire(CLine3D w, double distance, std::vector<CLine3D>& lines)
{
    double perp_x_dir = -w.getDirection().Y();
    double perp_y_dir =  w.getDirection().X();

    double x_perp = w.getPoint().X() + perp_x_dir * distance / sqrt(perp_x_dir*perp_x_dir + perp_y_dir*perp_y_dir);
    double y_perp = w.getPoint().Y() + perp_y_dir * distance / sqrt(perp_x_dir*perp_x_dir + perp_y_dir*perp_y_dir);
    TVector3 direction_perp(perp_x_dir, perp_y_dir, 0);
    TVector3 point_perp(x_perp, y_perp, 0);

    lines.push_back(CLine3D(point_perp, w.getDirection()));
    
    perp_x_dir =  w.getDirection().Y();
    perp_y_dir = -w.getDirection().X();

    x_perp = w.getPoint().X() + perp_x_dir * distance / sqrt(perp_x_dir*perp_x_dir + perp_y_dir*perp_y_dir);
    y_perp = w.getPoint().Y() + perp_y_dir * distance / sqrt(perp_x_dir*perp_x_dir + perp_y_dir*perp_y_dir);
    direction_perp = TVector3(perp_x_dir, perp_y_dir, 0);
    TVector3 point_perp2(x_perp, y_perp, 0);

    lines.push_back(CLine3D(point_perp2, w.getDirection()));
}

void TrackletFinder::ComputeCellsIntersections()
{
  std::vector<CLine3D> lines;

  if (_fired_cells) {
    for (auto c:(*_fired_cells)) {
      double h, w;
      c.second.size(h, w);
      CLine3D line_from_wire(c.second.wire().center(), c.second.wire().getDirection());
      LinesParallelToWire(line_from_wire, w, lines);
    }
  }
  
  for (uint i = 0; i < lines.size() - 1; i++) {
    for (uint j = i + 1; j < lines.size(); j++) {

      bool is_parallel = CheckParallel(lines[i].getDirection(), lines[j].getDirection());
      double len = 0;
      if (!is_parallel) {
        len = (lines[i].getPoint() - lines[j].getPoint()).Dot(lines[i].getDirection() - lines[j].getDirection() * lines[i].getDirection().Dot(lines[j].getDirection())) /
                     (pow(lines[i].getDirection().Dot(lines[j].getDirection()), 2) - 1);
        
        // double len = ((lines[i].getPoint().Y() - lines[j].getPoint().Y()) - lines[j].getDirection().Y() / lines[j].getDirection().X() * (lines[i].getPoint().X() - lines[j].getPoint().X())) /
        //              (lines[i].getDirection().Y() - lines[j].getDirection().Y() / lines[j].getDirection().X() * lines[i].getDirection().X());
      
        TVector3 intersection = lines[i].getPoint() + len * lines[i].getDirection();
        _cells_intersections.push_back(intersection);
      }
    }
  }

  double x_min = 1E9;
  double y_min = 1E9;
  double x_max = -1E9;
  double y_max = -1E9;
  for (auto p:_cells_intersections) {
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
  if (x_min <  _volume_parameters[1]) {
    x_min = _volume_parameters[1];
  }
  if (y_min <  _volume_parameters[4]) {
    y_min = _volume_parameters[4];
  }
  if (x_max >  _volume_parameters[2]) {
    x_max = _volume_parameters[2];
  }
  if (y_max >  _volume_parameters[5]) {
    y_max = _volume_parameters[5];
  }

  _cells_intersections.clear();
  _cells_intersections.push_back(TVector3(x_min, y_min, 0));
  _cells_intersections.push_back(TVector3(x_max, y_max, 0));
}

void TrackletFinder::ComputeCellsBands()
{
  if (_fired_cells) {
    for (auto c:(*_fired_cells)) {
      double h, w;
      c.second.size(h, w);
      CLine3D line_from_wire(c.second.wire().center(), c.second.wire().getDirection());
      LinesParallelToWire(line_from_wire, w, _cells_bands);
    }
  }
}

std::vector<TVectorD> TrackletFinder::FindTracklets()
{
  auto minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "");
  minimizer->SetMaxFunctionCalls(100000); 
  minimizer->SetMaxIterations(100000);
  minimizer->SetTolerance(1e-9);

  ComputeCellsIntersections();
  ComputeCellsBands();
  auto cells = _fired_cells;
  ROOT::Math::Functor functor_cells([&cells](const double* params) { return MinimizingFunction(params, cells); }, 4);
  
  if (_fired_cells) minimizer->SetFunction(functor_cells);
    
  double theta_xz = atan(_trajectory.getDirection().Z() / _trajectory.getDirection().X());
  double theta_yz = atan(_trajectory.getDirection().Y() / _trajectory.getDirection().Z());
  double theta_width = 0.4;

  std::cout << _cells_intersections[1].X() << " " << _cells_intersections[0].X() << std::endl;
  std::cout << _cells_intersections[1].Y() << " " << _cells_intersections[0].Y() << std::endl;

  double x_width = _cells_intersections[1].X() - _cells_intersections[0].X();
  double y_width = _cells_intersections[1].Y() - _cells_intersections[0].Y();
  
  int subdivisions = 3;
  double x_sub_width = x_width / subdivisions;
  double y_sub_width = y_width / subdivisions;
  double theta_sub_width = theta_width / subdivisions;

  std::vector<TVectorD> sampling_points;
  for (int i = 0; i <= subdivisions; i++) {
    for (int j = 0; j <= subdivisions; j++) {
      for (int k = 0; k <= subdivisions; k++) {
        for (int z = 0; z <= subdivisions; z++) {
          TVectorD point(4);

          point[0] = _cells_intersections[0].X() + i * x_sub_width;
          point[1] = _cells_intersections[0].Y() + j * y_sub_width;
          point[2] = theta_xz - theta_width + k * theta_sub_width;
          point[3] = theta_yz - theta_width + z * theta_sub_width;

          sampling_points.push_back(point);
        }
      }
    }  
  }

  std::vector<TVectorD> minima;
  for (uint i = 0; i <sampling_points.size(); i++) {
    double starting_point[4] = {sampling_points[i][0], sampling_points[i][1], sampling_points[i][2], sampling_points[i][3]};
    minimizer->SetLimitedVariable(0, "px", starting_point[0], 0.01, _cells_intersections[0].X(), _cells_intersections[1].X());
    minimizer->SetLimitedVariable(1, "py", starting_point[1], 0.01, _cells_intersections[0].Y(), _cells_intersections[1].Y());
    minimizer->SetLimitedVariable(2, "dx", starting_point[2], 0.001, theta_xz - theta_width, theta_xz + theta_width);
    minimizer->SetLimitedVariable(3, "dy", starting_point[3], 0.001, theta_yz - theta_width, theta_yz + theta_width);

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

  for (uint i = 0; i < minima.size() - 1; i++) {
    for (uint j = i + 1; j < minima.size(); j++) {
      if (fabs(minima[i][0] - minima[j][0]) < 2 * _sigma_pos &&
          fabs(minima[i][1] - minima[j][1]) < 2 * _sigma_pos &&
          fabs(minima[i][2] - minima[j][2]) < 2 * _sigma_ang &&
          fabs(minima[i][3] - minima[j][3]) < 2 * _sigma_ang) {
            minima.erase(minima.begin() + j);
            j--;
      }
    }
  }
  for (uint i = 0; i < minima.size(); i++) {
    std::cout << minima[i][0] << " " << minima[i][1] << " " << cos(minima[i][2]) << " " << sin(minima[i][3]) << " " << sin(minima[i][2]) << " " << minima[i][4]<< std::endl;
  }
  std::cout << "DONE" << std::endl;            
  
  return minima;
}

void TrackletFinder::Clear()
{
  _fired_cells = nullptr;
  _cells_intersections.clear();
  _cells_bands.clear();
}

void TrackletFinder::Draw3DWires() {
  gStyle->SetOptStat(0);
  _c3 = new TCanvas("c3D","c3D",1500,1500);
  TH3D h("", "", _volume_parameters[0],_volume_parameters[1], _volume_parameters[2],
                  _volume_parameters[3],_volume_parameters[4], _volume_parameters[5],
                  _volume_parameters[6],_volume_parameters[7], _volume_parameters[8]);
  h.SetTitle(";x;y;z");
  h.Draw();

  for (auto c:(*_fired_cells)) {  
    TPolyLine3D* pl2 = new TPolyLine3D(2);
    TVector3 l_start = c.second.wire().getPoints()[0];
    pl2->SetPoint(0, l_start.X(), l_start.Y(), l_start.Z());
    TVector3 l_end = c.second.wire().getPoints()[1];
    pl2->SetPoint(1, l_end.X(), l_end.Y(), l_end.Z());
    pl2->Draw("same");
  }

  _c3->SaveAs("./c3D.C");
}

TVector3 GetCylinderCoordinates(CLine3D w, double t, double radius, double theta) {
  double x = w.getPoint().X() + t*w.getDirection().X() + radius*(cos(theta)*w.getU().X() + sin(theta)*w.getV().X());
  double y = w.getPoint().Y() + t*w.getDirection().Y() + radius*(cos(theta)*w.getU().Y() + sin(theta)*w.getV().Y());
  double z = w.getPoint().Z() + t*w.getDirection().Z() + radius*(cos(theta)*w.getU().Z() + sin(theta)*w.getV().Z());

  return TVector3(x, y, z);
}

void TrackletFinder::Draw3D()
{
  gStyle->SetOptStat(0);
  if (!_c3) {
    _c3 = new TCanvas("c3D","c3D",1500,1500);
  }

  int ccc = 3;
  for (auto c:(*_fired_cells)) {
    TH3D* h2 = new TH3D("", "", _volume_parameters[0],_volume_parameters[1], _volume_parameters[2],
                                _volume_parameters[3],_volume_parameters[4], _volume_parameters[5],
                                _volume_parameters[6],_volume_parameters[7], _volume_parameters[8]);
    h2->SetTitle(";x;y;z");
    h2->SetMarkerColor(ccc);
    CLine3D line_from_wire(c.second.wire().center(), c.second.wire().getDirection());
    for (int t = -100; t < 100; t++) {
      for (int theta = 0; theta < 628; theta++) {
        TVector3 p = GetCylinderCoordinates(line_from_wire, (double)t, c.second.driftVelocity() * c.first.first, theta/100.);
        h2->Fill(p.X(), p.Y(), p.Z());
      }
    }
    h2->Draw("same scat");
    ccc++;
  }

  _c3->SaveAs("./c3D.png");
  
}

void TrackletFinder::Draw2DWires()
{
  gStyle->SetOptStat(0);
  if (!_c2) {
    _c2 = new TCanvas("c2D","c2D",1500,1500);
    TH2D* h2 = new TH2D("h","h", _cells_intersections[1].X()  - _cells_intersections[0].X(), _cells_intersections[0].X(), _cells_intersections[1].X(),
                                 _cells_intersections[1].Y()  - _cells_intersections[0].Y(), _cells_intersections[0].Y(), _cells_intersections[1].Y());

    h2->Draw();
  }

  for (auto l:_cells_bands) {
    TVector3 second_point;
    second_point = l.getPoint() + 100*l.getDirection();
    TLine* tl = new TLine(l.getPoint().X(), l.getPoint().Y(), second_point.X(), second_point.Y());
    tl->SetLineStyle(2);
    tl->Draw("same");
  }
  
  for (auto c:(*_fired_cells)) {
    TVector3 l_start = c.second.wire().getPoints()[0];
    TVector3 l_end   = c.second.wire().getPoints()[1];

    TLine* tl = new TLine(l_start.X(), l_start.Y(), l_end.X(), l_end.Y());
    tl->Draw("same");
  }

  _c2->SaveAs("./c2D.png");

}

void TrackletFinder::Draw2DDistance()
{
  gStyle->SetOptStat(0);
  TCanvas c("c2DMinimization","c2DMinimization",1500,1500);
  if (!_c2) {
    _c2 = new TCanvas("c2D","c2D",1500,1500);

  }
  
  TH2D* h2 = new TH2D("h","h", _cells_intersections[1].X()  - _cells_intersections[0].X(), 
                               _cells_intersections[0].X(), _cells_intersections[1].X(),
                               _cells_intersections[1].Y()  - _cells_intersections[0].Y(), 
                               _cells_intersections[0].Y(), _cells_intersections[1].Y());

  double min;
  auto cells = _fired_cells;

  double theta_xz = atan(_trajectory.getDirection().Z() / _trajectory.getDirection().X()) * 1000;
  double theta_yz = atan(_trajectory.getDirection().Y() / _trajectory.getDirection().Z()) * 1000;
  int count = 0;
  for (int px = _cells_intersections[0].X(); px <= _cells_intersections[1].X(); px++) {
    for (int py = _cells_intersections[0].Y(); py <= _cells_intersections[1].Y(); py++) {
      min = 1E9;
      for (int angle_xz = theta_xz - 400; angle_xz < theta_xz + 400 ; angle_xz+=10) {
        for (int angle_yz = theta_yz - 400; angle_yz < theta_yz + 400; angle_yz+=10) {
          double p[4] = {px / 1., py / 1., angle_xz / 1000., angle_yz / 1000.};
          double tmp_min = MinimizingFunction(p, cells);
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

  _c2->SaveAs("./c2D.png");

}