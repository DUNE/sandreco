#include <TVector3.h>
#include <TTree.h>
#include <TFile.h>
#include <TMarker.h>
#include <TArrow.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <random>

#include "SANDTrackerVertexing.h"
#include "struct.h"
#include "utils.h"

int main(int argc, char* argv[])
{
  gStyle->SetOptStat(0);  

  if (argc != 7) {
    std::cout << "TestVertexing <number of slices> <tracks position sigma> " 
              << "<track angle sigma> <dz> <impact parameter> <merging radius> (all mm)" << std::endl;
    return 1;
  }
  
  TRandom3 r(0);

  int n_slice = std::stoi(argv[1]);

  
  for (int i = 0; i < n_slice; i++) {
    int n_vertices = r.Uniform(1,9);
    std::cout << "Generating slice with " << n_vertices << " vertices" << std::endl;
    
    std::vector<Track> tracks;
    for (int v = 0; v < n_vertices; v++) {
      std::cout << "Generating vertex " << v << std::endl;
      
      double vtx_x = r.Uniform(0,100);
      double vtx_y = r.Uniform(0,100);
      double vtx_z = r.Uniform(0,100);
      
      std::cout << "  > position = (" << vtx_x << ", " << vtx_y << ", " << vtx_z << ")" << std::endl;
      
      int n_tracks = r.Uniform(1, 9);
      std::cout << "  > n tracks = " << n_tracks << std::endl;
      
      for (int t = 0; t < n_tracks; t++) {
        std::cout << "  > Generating track " << t << std::endl;
        Track track;
        track.id = n_slice * 100 + v * 10 + t;
        track.x  = vtx_x + r.Gaus(0, std::stoi(argv[2]));
        track.y  = vtx_y + r.Gaus(0, std::stoi(argv[2]));
        track.z  = vtx_z + r.Gaus(0, std::stoi(argv[2]));
        
        track.tx = tan(r.Uniform(-M_PI_2, M_PI_2)) + r.Gaus(0, std::stoi(argv[3]));
        track.ty = tan(r.Uniform(-M_PI_2, M_PI_2)) + r.Gaus(0, std::stoi(argv[3]));
        std::cout << "    > tan angle x = " << track.tx << std::endl;
        std::cout << "    > tan angle y = " << track.ty << std::endl;

        tracks.push_back(track);
      }
      
    }
    std::cout << "Total tracks: " << tracks.size() << std::endl;

    TrackerVertexing vertex_finder;
    vertex_finder.setParameters(std::stoi(argv[4]), std::stoi(argv[5]), std::stoi(argv[6]), tracks);
    vertex_finder.run();

  }

  return 0;
}