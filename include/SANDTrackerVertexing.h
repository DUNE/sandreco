//////////////////////////////////
//      Author: M. Pozzato      //
//////////////////////////////////

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "struct.h"

class Track {
  public:
   int id;
   double x, y, z, tx, ty;
   std::vector<int> indexVertex;
   Track() {
     x = y = z = tx = ty = 0.;
     id = 0;
   }
   ~Track(){};
 };

class Vertex {
  public:
   int id, flag;  // 1 = neutral forward; 2 = neutral backward; 3 = charged;
   double x, y, z;
   int nProng;
   std::vector<double> impact;
   std::vector<int> indexTrack;
   Vertex() {
     x = y = z = 0.;
     id = nProng = 0;
   }
   ~Vertex(){};
   void clean() {
     x = y = z = 0.;
     id = nProng = 0;
     impact.clear();
     indexTrack.clear();
   }
 };

class TrackerVertexing {

  public:

  private:

    int run();
    void setParameters(double dz, double ip, double merging_radius, std::vector<Track> tracks);
    int doVertex();
    int vertexEstimate(Track tr1, int tr1Index, Track tr2, int tr2Index, Vertex &vtx);
    void flagVertex();
    void dumpVertex(const char *outputFile);
    int selectVertex();
    void clearAll();
    int mergeVertex();
    void refineVertexPosition(double stepSize = 0.001, int nSteps = 50);

    int id_vert_;
    double dz_;
    double ip_;
    double merging_radius_;
    std::vector<Track> tracks_;
    std::vector<Vertex> vertices_;
    std::vector<Vertex> vertices_2_prong_;
    std::vector<Vertex> vertices_multi_prong_;
    std::vector<std::vector<Vertex>> vertices_list;
  
};