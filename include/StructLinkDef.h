#ifdef __CINT__
#include "struct.h"

//#pragma link C++ class std::map < int, std::vector < double>> + ;
//#pragma link C++ class std::map < int, std::vector < int>> + ;
//#pragma link C++ class std::map < int, double> + ;
#pragma link C++ class std::vector < pe> + ;
#pragma link C++ class std::vector < truecluster> + ;
#pragma link C++ class std::vector < cluster_generator> + ;
#pragma link C++ class std::vector < dg_ps> + ;
#pragma link C++ class std::vector < dg_cell> + ;
#pragma link C++ class std::vector < reco_cell > +;
#pragma link C++ class std::vector < vertex > +;
//#pragma link C++ class std::map < std::string, std::vector < hit>> + ;
#pragma link C++ class std::vector < dg_wire > +;
#pragma link C++ class std::vector < track > +;
#pragma link C++ class std::vector < cluster > +;
#pragma link C++ class particle+;
#pragma link C++ class std::vector < particle> + ;
#pragma link C++ class pe + ;
#pragma link C++ class cluster_generator + ;
#pragma link C++ class truecluster + ;
#pragma link C++ class dg_ps + ;
#pragma link C++ class dg_wire + ;
#pragma link C++ class dg_cell + ;
#pragma link C++ class reco_cell + ;
// #pragma link C++ class incomplete_cell + ;
#pragma link C++ class cluster + ;
#pragma link C++ class track + ;
// #pragma link C++ class particle + ;
#pragma link C++ class event + ;
#pragma link C++ class vertex + ;

// NEW --SILVIA     
#pragma link C++ class grain_event + ;
#pragma link C++ class lar_track + ;
#pragma link C++ class lar_point + ;
#pragma link C++ class track_grain_lens + ;
#pragma link C++ class vertex_grain_lens + ;
#pragma link C++ class std::vector<track_grain_lens>+;
#pragma link C++ class std::vector<vertex_grain_lens>+;

#endif
