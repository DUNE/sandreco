#include <math.h>

namespace tracker_module_configuration 
{
  namespace drift
  {
    std::map<std::string, double> id_to_angle =
    {
      {"0", 0},
      {"1", -M_PI / 36.},
      {"2", M_PI / 36.}
    };
    
    std::map<std::string, double> id_to_offset =
    {
      {"0", 10},
      {"1", 10},
      {"2", 10}
    };
    
    std::map<std::string, double> id_to_spacing =
    {
      {"0", 10},
      {"1", 10},
      {"2", 10}
    };
    
    std::map<std::string, double> id_to_length =
    {
      {"0", 0},
      {"1", 0},
      {"2", 0}
    };

    std::map<std::string, double> id_to_velocity =
    {
      {"0", 0.05},
      {"1", 0.05},
      {"2", 0.05}
    };
  }
  
  namespace stt
  {
    std::map<std::string, double> id_to_angle =
    {
      {"1", M_PI_2},
      {"2", 0}
    };

    std::map<std::string, double> id_to_velocity =
    {
      {"1", 0.05},
      {"2", 0.05},
    };
  }
}