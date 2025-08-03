#ifndef EXPAND_3D_HPP
#define EXPAND_3D_HPP

#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>
#include <netcdf>
#include <boost/lexical_cast.hpp>

/**
 * @brief Parse the requested string of wavelength
 *
 * @param wv_list user supplied wv3d list
 * @param wv3d_list  wavelength and their indexes;
 * @param wv_requested_vals user requested wavelength list
 * @param wv_requested_indexes slices indexes of the 3D variable
 * @param exit_on_not_found hard exit if wavelength not found
 */
void parse_wv_list(const std::string &wv_list, const std::unordered_map<int, int> &wv3d_list,
                   std::vector<int> &wv_requested_vals, std::vector<int> &wv_requested_indexes,
                   bool exit_on_not_found = true);
#endif