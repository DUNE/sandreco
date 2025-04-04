#include "struct.h"

#include <vector>
#include <numeric>

int CountPhotons(const std::vector<grain_dense_image>& raw_sipm_data) {
    int total_photons = std::accumulate(raw_sipm_data.begin(), raw_sipm_data.end(), 0, 
        [](int sum, const grain_dense_image& s) {return sum + s.h_indices_img.size();});
    return total_photons;
}

float CountCharge(const std::vector<grain_dense_image>& raw_sipm_data) {
    float total_charge = std::accumulate(raw_sipm_data.begin(), raw_sipm_data.end(), 0.0f, 
        [](float sum, const grain_dense_image& s) {return sum + std::accumulate(s.charge.begin(), s.charge.end(), 0.0f);});
    return total_charge;
}