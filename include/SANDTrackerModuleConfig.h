#include <math.h>

namespace TrackerModuleConfiguration 
{
  std::map<std::string, double> _id_to_angle =
  {
    {"0", 0},
    {"1", 0},
    {"2", M_PI_4}
  };
  
  std::map<std::string, double> _id_to_offset =
  {
    {"0", 10},
    {"1", 5},
    {"2", 10}
  };
  
  std::map<std::string, double> _id_to_spacing =
  {
    {"0", 10},
    {"1", 10},
    {"2", 10}
  };
};