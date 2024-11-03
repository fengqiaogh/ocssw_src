#ifndef EXPAND_3D_HPP
#define EXPAND_3D_HPP

#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>
#include <netcdf>
#include "genutils.h"
#include "readL2scan.h"
#include <boost/lexical_cast.hpp>
#define EXIT_LOG(...)                                                                     \
    {                                                                                     \
        __VA_ARGS__;                                                                      \
        std::cerr << "Exiting. See  " << __FILE__ << ":" << __LINE__ << std::endl; \
        exit(EXIT_FAILURE);                                                               \
    }
void set_mapper(const std::string &input_file_name, std::string &products_requested, std::vector<std::string> &products_requested_separated, const std::string &requested_wavelengths);
void set_prodname_3d_to_l2(const std::vector<std::string> &prodparam, l2_prod &l2_str, std::vector<std::string> &l2_prodname, std::vector<std::string> &l3_prodname, std::vector<int32_t> &thirdDimId, std::vector<float> &min_value, std::vector<float> &max_value);
bool set_l2_flags_use(const std::string &flagsuse);
void parse_wv_list(const std::string &wv_list,
                   const std::unordered_map<int, int> &wv3d_list, std::vector<int> & wv_requested_vals, 
                  std::vector<int>  & wv_requested_indexes, bool exit_on_not_found = true);
float l3_attr(const std::string &inp3);
#endif