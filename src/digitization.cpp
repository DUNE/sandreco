#include "SANDDigitization.h"
#include "SANDDigitizationEDEPSIM.h"
#include "SANDDigitizationFLUKA.h"

#include <iostream>

enum class DETSIM_TYPE {kEdepsim, kFluka};

void help_digit()
{
  std::cout << "Digitize <MC file> <digit file>" << std::endl;
  std::cout << "MC file name could contain wild card" << std::endl;
  
  std::cout
      << "usage: Digitize <MC file> <digit file> [detsim_type]\n";
  std::cout << "    - detsim_type: 'detsim_type::edepsim' (default) \n";
  std::cout << "                   'detsim_type::fluka' \n";
}

int main(int argc, char* argv[])
{
  if (argc < 3 || argc > 4) {
    help_digit();
    return -1;
  }

  if (argc > 3 && strcmp(argv[3], "detsim_type::fluka") == 0) {
    std::cout << "DETSIM_TYPE: FLUKA\n";
    digitization::fluka::digitize(argv[1], argv[2]);
  } else {
    std::cout << "DETSIM_TYPE: EDEPSIM\n";
    digitization::edep_sim::digitize(argv[1], argv[2]);
  }
}
