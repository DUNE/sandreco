// File loader.C 
#include <vector> 
#include <map>

#ifndef LOADER_C
#define LOADER_C

struct cell {
  int id;
  double z;
  double y;
  double adc1;
  double tdc1;
  double adc2;
  double tdc2;
  int mod;
  int lay;
  int cel;
  std::vector<double> pe_time1;
  std::vector<int> hindex1;
  std::vector<double> pe_time2;
  std::vector<int> hindex2;
};

struct hit {
  std::string det;
  double x1;
  double y1;
  double z1;
  double t1;
  double x2;
  double y2;
  double z2;
  double t2;
  double de;
  int pid;
  int index;
};

struct digit {
  std::string det;
  double x;
  double y;
  double z;
  double t;
  double de;
  bool hor;
  std::vector<int> hindex;
};

struct track {
  int tid;
  double yc;
  double zc;
  double r;
  double a;
  double b;
  double h;
  double x0;
  double y0;
  double z0;
  double t0;
  std::vector<digit> digits;
};

class gcell {
  public:
    int id;
    double Z[4];
    double Y[4];
    double adc;
    double tdc;
};

bool isHitBefore(hit h1, hit h2)
{
  return h1.t1 < h2.t1;
}

bool isDigBefore(digit d1, digit d2)
{
  return d1.t < d2.t;
}

#endif

#ifdef __MAKECINT__ 
#pragma link C++ class std::map<int,std::vector<double> >+; 
#pragma link C++ class std::map<int,std::vector<int> >+;
#pragma link C++ class std::map<int,double>+;
#pragma link C++ class std::vector<cell>+; 
#pragma link C++ class std::map<std::string,std::vector<hit> >+; 
#pragma link C++ class std::vector<digit>+; 
#pragma link C++ class std::vector<track>+;
#endif
