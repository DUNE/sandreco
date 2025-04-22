#include "struct.h"

#include <vector>

#ifndef SANDGRAININFO_H
#define SANDGRAININFO_H


int CountPhotons(const std::vector<grain_dense_image>& raw_sipm_data);

float CountCharge(const std::vector<grain_dense_image>& raw_sipm_data);

#endif
