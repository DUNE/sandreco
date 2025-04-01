///////////////////////////////////////////
//      Original author: M. Pozzato      //
//      Ported by: V.pia                 //
///////////////////////////////////////////

#include "SANDTrackerVertexing.h"

void TrackerVertexing::clearAll() {
  tracks_.clear();
  vertices_.clear();
  vertices_2_prong_.clear();
  vertices_multi_prong_.clear();
  vertices_list_.clear();
}

int TrackerVertexing::vertexEstimate(Track tr1, int tr1Index, Track tr2, int tr2Index,
                   Vertex &vtx) {
  // first line
  double a, b, c, d, e, f, t1;
  // second line
  double a1, b1, c1, d1, e1, f1, t2;
  
  a = tr1.x - tr1.tx * tr1.z;
  b = tr1.tx;
  c = tr1.y - tr1.ty * tr1.z;
  d = tr1.ty;
  e = 0;
  f = 1;

  a1 = tr2.x - tr2.tx * tr2.z;
  b1 = tr2.tx;
  c1 = tr2.y - tr2.ty * tr2.z;
  d1 = tr2.ty;
  e1 = 0;
  f1 = 1;

  double k1 = (-b * b - d * d - f * f);
  double k2 = (b1 * b + d1 * d + f1 * f);
  double k3 = b * (a1 - a) + d * (c1 - c) + f * (e1 - e);
  double k11 = (-d * d1 - b * b1 - f * f1);
  double k21 = (b1 * b1 + d1 * d1 + f1 * f1);
  double k31 = b1 * (a1 - a) + d1 * (c1 - c) + f1 * (e1 - e);
  t2 = (k11 * k3 / k1 - k31) / (k21 - k11 * k2 / k1);
  t1 = (-k2 * t2 - k3) / k1;
  double px = a + b * t1;
  double py = c + d * t1;
  double pz = e + f * t1;

  double p1x = a1 + b1 * t2;
  double p1y = c1 + d1 * t2;
  double p1z = e1 + f1 * t2;

  double vx = (px + p1x) / 2.;
  double vy = (py + p1y) / 2.;
  double vz = (pz + p1z) / 2.;

  double dist = sqrt((px - p1x) * (px - p1x) + (py - p1y) * (py - p1y) +
                     (pz - p1z) * (pz - p1z));

  vtx.clean();
  vtx.x = vx;
  vtx.y = vy;
  vtx.z = vz;
  vtx.id = id_vert_;
  id_vert_++;
  vtx.nProng = 2;
  vtx.indexTrack.push_back(tr1Index);
  vtx.indexTrack.push_back(tr2Index);
  if (dist == 0) dist = 1E-15;
  vtx.impact.push_back(dist / 2.);
  vtx.impact.push_back(dist / 2.);
  return 0;
}

int TrackerVertexing::doVertex() {
  int nTracks = static_cast<int>(tracks_.size());
  Vertex vtx;
  for (int i = 0; i < nTracks - 1; i++) {
    for (int j = i + 1; j < nTracks; j++) {
      vertexEstimate(tracks_.at(i), i, tracks_.at(j), j, vtx);

      if (fabs(vtx.z - tracks_.at(i).z) < dz_ && fabs(vtx.z - tracks_.at(j).z) < dz_ &&
          vtx.impact.at(0) < ip_ && vtx.impact.at(1) < ip_) {
        vertices_.push_back(vtx);
      }
    }
  }
  return static_cast<int>(vertices_.size());
}

void TrackerVertexing::dumpVertex(const char *outputFile) {
  if (outputFile) {
    std::ofstream outfile(outputFile);

    int nVertex = static_cast<int>(vertices_.size());
    std::cout << "nVertex: " << nVertex << std::endl;
    for (int i = 0; i < nVertex; i++) {
    //  if (vertices_.at(i).flag == 2) continue; // remove from output file reversed vertex
      outfile << vertices_.at(i).id << "\t" << vertices_.at(i).nProng << "\t"
              << vertices_.at(i).x << "\t" << vertices_.at(i).y << "\t"
              << vertices_.at(i).z << "\t" << vertices_.at(i).flag << std::endl;
      for (int j = 0; j < static_cast<int>(vertices_.at(i).indexTrack.size());
           j++) {
        outfile << tracks_.at(vertices_.at(i).indexTrack.at(j)).id << "\t"
                << tracks_.at(vertices_.at(i).indexTrack.at(j)).x << "\t"
                << tracks_.at(vertices_.at(i).indexTrack.at(j)).y << "\t"
                << tracks_.at(vertices_.at(i).indexTrack.at(j)).z << "\t"
                << tracks_.at(vertices_.at(i).indexTrack.at(j)).tx << "\t"
                << tracks_.at(vertices_.at(i).indexTrack.at(j)).ty << "\t"
                << vertices_.at(i).impact.at(j) << std::endl;
      }
    }
  }
}

int TrackerVertexing::selectVertex() {
  int nVertex = static_cast<int>(vertices_.size());
  std::vector<Vertex>::iterator iter;
  std::vector<Vertex> tmp;
  std::vector<int> index;
  double x = vertices_.at(0).x;
  double y = vertices_.at(0).y;
  double z = vertices_.at(0).z;
  index.clear();
  tmp.clear();
  index.push_back(0);
  int merged = 0;
  int indexPos = 0;

  double sumx = x;
  double sumy = y;
  double sumz = z;

  for (iter = vertices_.begin() + 1; iter != vertices_.end(); iter++) {
    double dist = sqrt((x - (*iter).x) * (x - (*iter).x) +
                       (y - (*iter).y) * (y - (*iter).y) +
                       (z - (*iter).z) * (z - (*iter).z));
    indexPos++;
    if (dist < merging_radius_) {
      index.push_back(indexPos);
      tmp.push_back((*iter));
      sumx += (*iter).x;
      sumy += (*iter).y;
      sumz += (*iter).z;
      merged++;
      x = sumx / (merged + 1);
      y = sumy / (merged + 1);
      z = sumz / (merged + 1);
    }
  }
  if (merged == 0) {
    vertices_2_prong_.push_back(vertices_.at(0));
    vertices_.erase(vertices_.begin());
  } else {
    tmp.push_back(vertices_.at(0));
    vertices_list_.push_back(tmp);
    for (int j = static_cast<int>(index.size()) - 1; j >= 0; j--) {
      iter = vertices_.begin();
      vertices_.erase(iter + index.at(j));
    }
  }
  return 0;
}

double pointTo3DlineDistance(double x0, double y0, double z0, Track tr) {
  double dist = 0;
  double qp_x = tr.x - x0;
  double qp_y = tr.y - y0;
  double qp_z = tr.z - z0;
  double vx = tr.tx;
  double vy = tr.ty;
  double vz = 1;
  double qr =
      (qp_x * vx + qp_y * vy + qp_z * vz) / sqrt(vx * vx + vy * vy + vz * vz);
  double qp = sqrt(qp_x * qp_x + qp_y * qp_y + qp_z * qp_z);
  dist = sqrt(qp * qp - qr * qr);
  return dist;
}

int TrackerVertexing::mergeVertex() {
  int nVertex = static_cast<int>(vertices_list_.size());
  if (nVertex == 0) {
    std::cout << "No vertex to merge" << std::endl;
    return 0;
  }
  for (int i = 0; i < nVertex; i++) {
    int nVtx = static_cast<int>(vertices_list_.at(i).size());
    double x = 0;
    double y = 0;
    double z = 0;
    double w = 0;
    double wSum = 0;

    for (int j = 0; j < nVtx; j++) {
      w = 1. / (2 * (vertices_list_.at(i).at(j).impact.at(0)));
      x += vertices_list_.at(i).at(j).x * w;
      y += vertices_list_.at(i).at(j).y * w;
      z += vertices_list_.at(i).at(j).z * w;
      wSum += w;
    }
    x /= wSum;
    y /= wSum;
    z /= wSum;
    Vertex vt;
    vt.x = x;
    vt.y = y;
    vt.z = z;
    vt.id = i;
    std::vector<int> trkIndex;

    for (int j = 0; j < nVtx; j++) {
      int nProng = vertices_list_.at(i).at(j).nProng;
      for (int k = 0; k < nProng; k++) {
        trkIndex.push_back(vertices_list_.at(i).at(j).indexTrack.at(k));
      }
    }

    std::sort(trkIndex.begin(), trkIndex.end());

    int oldIndex = -999;
    int nIndex = static_cast<int>(trkIndex.size());
    int prongs = 0;
    for (int k = 0; k < nIndex; k++) {
      if (trkIndex.at(k) != oldIndex) {
        oldIndex = trkIndex.at(k);
        double dist = pointTo3DlineDistance(x, y, z, tracks_.at(oldIndex));
        vt.indexTrack.push_back(oldIndex);
        vt.impact.push_back(dist);
        prongs++;
      }
    }
    vt.nProng = prongs;

    vertices_multi_prong_.push_back(vt);
  }

  return nVertex;
}

void TrackerVertexing::refineVertexPosition(double stepSize, int nSteps) {
    int nVertex = static_cast<int>(vertices_multi_prong_.size());
    double length = stepSize * nSteps;
    
    for (int i = 0; i < nVertex; i++) {
        int nProng = vertices_multi_prong_.at(i).nProng;
        double sumIP = 1e3;
        double newSumIP = 0;
        double vx = vertices_multi_prong_.at(i).x;
        double vy = vertices_multi_prong_.at(i).y;
        double vz = vertices_multi_prong_.at(i).z;
        
        double startX = vx - 0.5 * length;
        double startY = vy - 0.5 * length;
        double startZ = vz - 0.5 * length;
        int best_j = 0, best_k = 0, best_l = 0;;

        for (int j = 0; j < nSteps; j++) {
            vx = startX + j * stepSize;
            for (int k = 0; k < nSteps; k++) {
                vy = startY + k * stepSize;
                for (int l = 0; l < nSteps; l++) {
                    newSumIP = 0;
                    vz = startZ + l * stepSize;
                    for (int n = 0; n < nProng; n++) {
                        double dist = pointTo3DlineDistance(
                            vx, vy, vz, tracks_.at(vertices_multi_prong_.at(i).indexTrack.at(n)));
                        newSumIP += dist;
                    }
                    if (newSumIP < sumIP)
                    {
                        sumIP = newSumIP;
                        best_j = j;
                        best_k = k;
                        best_l = l;
                    }
                }
            }
                      
        }

        vx = startX + best_j * stepSize;
        vy = startY + best_k * stepSize;
        vz = startZ + best_l * stepSize;
        vertices_multi_prong_.at(i).x = vx;
        vertices_multi_prong_.at(i).y = vy;
        vertices_multi_prong_.at(i).z = vz;
        for (int k = 0; k < nProng; k++) {
            double dist = pointTo3DlineDistance(
                vx, vy, vz, tracks_.at(vertices_multi_prong_.at(i).indexTrack.at(k)));
            vertices_multi_prong_.at(i).impact.at(k) = dist;
        }
    }
}

// Notice: this was done for OPERA. It is not updated to work with SAND.
void TrackerVertexing::flagVertex() {
  int nVertex = static_cast<int>(vertices_.size());
  bool isFirst;
  for (int i = 0; i < nVertex; i++) {
    isFirst = true;
    int flag = -1;  // 1 = neutral forword; 2 = neutral backword; 3 = charged;
    for (int j = 0; j < static_cast<int>(vertices_.at(i).indexTrack.size()); j++) {
      int trkIndex = vertices_.at(i).indexTrack.at(j);
      // tracks_.at(trkIndex).indexVertex.push_back(vertices_.at(i).id);
      tracks_.at(trkIndex).indexVertex.push_back(i);

      if (isFirst) {
        isFirst = false;
        if ((tracks_.at(trkIndex).z - vertices_.at(i).z) > 0)
          flag = 1;
        else
          flag = 2;
      } else {
        if ((tracks_.at(trkIndex).z - vertices_.at(i).z) > 0 && flag == 2)
          flag = 3;
        else if ((tracks_.at(trkIndex).z - vertices_.at(i).z) < 0 && flag == 1)
          flag = 3;
      }
    }
    vertices_.at(i).flag = flag;
  }
}

void TrackerVertexing::setParameters(double dz, double ip, double merging_radius, std::vector<Track> tracks) {
  dz_ = dz;
  ip_ = ip;
  merging_radius_ = merging_radius;
  tracks_ = tracks;

  std::cout << tracks_.size() << " tracks read." << std::endl;
  std::cout << "DZ: " << dz_ << "\nIP: " << ip_
            << "\nMerging Radius: " << merging_radius << std::endl;
}

int TrackerVertexing::run() {
  int pairVertexes = doVertex();
  dumpVertex("2-Prong.txt");
  std::cout << "Start selecting neighboor vertexes (" << pairVertexes
            << " - 2 Prong )" << std::endl;
  while (static_cast<int>(vertices_.size()) > 0) selectVertex();

  std::cout << "Start merging..." << std::endl;
  mergeVertex();

  refineVertexPosition();

  for (int i = 0; i < static_cast<int>(vertices_multi_prong_.size()); i++) {
    vertices_.push_back(vertices_multi_prong_.at(i));
  }
  int size = static_cast<int>(vertices_multi_prong_.size());
  for (int j = 0; j < static_cast<int>(vertices_2_prong_.size()); j++) {
    vertices_2_prong_.at(j).id = j + size;
    vertices_.push_back(vertices_2_prong_.at(j));
  }

  std::cout << "Flag vertexes\n";
  flagVertex();

  std::cout << "\n\n======== Results ========" << std::endl;
  std::cout << "2-Prong: " << vertices_2_prong_.size() << std::endl;
  std::cout << "Multi-Prong: " << vertices_multi_prong_.size() << std::endl;
  std::cout << "=========================" << std::endl;
  std::cout << "\nDump vertexes\n";
  dumpVertex("vertices.txt");

  clearAll();
  return 0;
}
