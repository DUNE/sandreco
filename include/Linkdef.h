#ifdef __CINT__
#include "struct.h"

#pragma link C++ class std::map < int, std::vector < double>> + ;
#pragma link C++ class std::map < int, std::vector < int>> + ;
#pragma link C++ class std::map < int, double> + ;
#pragma link C++ class std::vector < pe> + ;
#pragma link C++ class std::vector < dg_ps> + ;
#pragma link C++ class std::vector < dg_cell> + ;
#pragma link C++ class std::map < std::string, std::vector < hit>> + ;
#pragma link C++ class std::vector < dg_tube> + ;
#pragma link C++ class std::vector < track> + ;
#pragma link C++ class std::vector < cluster> + ;
#pragma link C++ class std::vector < particle> + ;
#pragma link C++ class pe + ;
#pragma link C++ class dg_ps + ;
#pragma link C++ class dg_tube + ;
#pragma link C++ class dg_cell + ;
#pragma link C++ class cluster + ;
#pragma link C++ class track + ;
#pragma link C++ class particle + ;
#pragma link C++ class event + ;
#pragma link C++ class Helix + ;
#pragma link C++ class Line + ;
#endif