#include "SANDDigitization.h"
#include "SANDDigitizationEDEPSIM.h"
#include "SANDDigitizationFLUKA.h"

#include <iostream>

void help_digit()
{
  std::cout << "usage: Digitize <MC file> <digit file> [detsim_type] "
               "[ecal_digi_mode]\n";
  std::cout << "    - detsim_type: 'detsim_type::edepsim' (default) \n";
  std::cout << "                   'detsim_type::fluka' \n";
  std::cout
      << "    - ecal_digi_mode: 'ecal_digi_mode::const_fract' (default) \n";
  std::cout << "                      'ecal_digi_mode::fixed_thresh' \n";
}

int main(int argc, char* argv[])
{
  if (argc < 3 || argc > 6) {
    help_digit();
    return -1;
  }

  auto detsim_type = digitization::DetSimType::kEdepsim;
  auto ecal_digi_mode = digitization::EcalDigiMode::const_fract;

  for (int i = 3; i < argc; i++) {
    if (strcmp(argv[i], "detsim_type::fluka") == 0) {
      detsim_type = digitization::DetSimType::kFluka;
    } else if (strcmp(argv[i], "ecal_digi_mode::fixed_thresh") == 0) {
      ecal_digi_mode = digitization::EcalDigiMode::fixed_thresh;
      sand_reco::ecal::acquisition::fixed_thresh_pe = atof(argv[++i]);
    }
  }

  std::cout << (detsim_type == digitization::DetSimType::kEdepsim
                    ? "DetSimType: EDEPSIM\n"
                    : "DetSimType: FLUKA\n");
  std::cout << (ecal_digi_mode == digitization::EcalDigiMode::const_fract
                    ? "EcalDigiMode: constant fraction\n"
                    : "EcalDigiMode: fixed threshold\n");

  if (detsim_type == digitization::DetSimType::kEdepsim) {
    digitization::edep_sim::digitize(argv[1], argv[2], ecal_digi_mode);
  } else {
    digitization::fluka::digitize(argv[1], argv[2], ecal_digi_mode);
  }
}
