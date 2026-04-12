#ifndef AREA_WEIGHTING_H
#define AREA_WEIGHTING_H
#include <vector>
#include <array>
#include <iostream>
std::vector<std::vector<std::pair<int, double>>> get_area_weights(
    std::vector<float>& lat_c, std::vector<float>& lon_c, std::vector<float>& lat_b,
    std::vector<float>& lon_b, size_t lines_c, size_t pixels_c, size_t lines_b, size_t pixels_b,
    std::vector<int>& index_search_results, bool compute_l1b_areas = false);
double get_l1b_area(int index);
#endif