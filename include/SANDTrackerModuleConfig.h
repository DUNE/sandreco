#include <math.h>

namespace TrackerModuleConfiguration 
{
  namespace Drift
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
      {"2", 20}
    };
    
    std::map<std::string, double> _id_to_length =
    {
      {"0", 500},
      {"1", 500},
      {"2", 500}
    };
  };
  
  namespace STT
  {
    std::map<std::string, double> _id_to_angle =
    {
      {"1", M_PI_2},
      {"2", 0}
    };
  }
};