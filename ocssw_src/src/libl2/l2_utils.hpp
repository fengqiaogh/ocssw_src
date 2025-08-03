#ifndef L2_UTILS_HPP
#define L2_UTILS_HPP
#include "scaled_nc_var.hpp"
#include <unordered_set>
#include <find_variable.hpp>
#include <cmath>
#include <algorithm>

const std::unordered_set<std::string> lat_possible_names = {
    "latitude", "lan", "Latitude"};  // should be made static and moved out to the cpp.
const std::unordered_set<std::string> lon_possible_names = {"longitude", "lon", "Longitude"};
const std::vector<std::string> possible_sozlen_names = {"solzen",
                                                        "solar_zenith_angle",
                                                        "solz",
                                                        "Solar Zenith Angle",
                                                        "csol_z",
                                                        "Center Solar Zenith Angle",
                                                        "center_solar_zenith_angle"};

const int brake_scan_fill_value{-9999};
const uint8_t scan_within_the_group{1};
const uint8_t scan_not_within_the_group{255};
const size_t bin_allocator_normalization{50000000};
const size_t max_per_bin_allocation{20};
const size_t min_per_bin_allocation{2};
const std::string min_max_token = "=;";
const std::string min_max_range_delimiter = ":";
const std::string product_delimiter = ", ";
#define EXIT_LOG(...)                                                              \
    {                                                                              \
        __VA_ARGS__;                                                               \
        std::cerr << "Exiting. See  " << __FILE__ << ":" << __LINE__ << std::endl; \
        exit(EXIT_FAILURE);                                                        \
    }

/**
 * @brief Find latitude and longitude variables in a netCDF file
 * @tparam File_Class Type of netCDF file class
 * @tparam Var_Class Type of netCDF variable class
 * @param file_nc The netCDF file object
 * @param lat_nc Output parameter for latitude variable
 * @param lon_nc Output parameter for longitude variable
 * @param message Additional error info to pass (such as filename)
 * @throws Exits if lat/lon variables not found or not 2D
 */
template <typename File_Class, typename Var_Class>
void find_lat_lon(File_Class &file_nc, Var_Class &lat_nc, Var_Class &lon_nc, const std::string &message = "") {
    auto possible_vars = find_all_variables<File_Class, Var_Class>(file_nc, lat_possible_names);
    if (possible_vars.empty()) {
        EXIT_LOG(std::cerr << "No valid latitude is found. " + message <<  std::endl)
    }
    // pick the first match
    lat_nc = possible_vars.begin()->second;
    possible_vars = find_all_variables<File_Class, Var_Class>(file_nc, lon_possible_names);
    if (possible_vars.empty()) {
        EXIT_LOG(std::cerr << "No valid longitude is found. " + message << std::endl)
    }
    // pick the first match
    lon_nc = possible_vars.begin()->second;
    if (lat_nc.getDimCount() != 2 || lon_nc.getDimCount() != 2) {
        EXIT_LOG(std::cerr << "Latitude and Longitude are not two dimensional variables. " + message << std::endl)
    }
};

/**
 * @brief Check if value is valid latitude
 * @tparam T Numeric type (float, double, etc)
 * @param lat Latitude value to check
 * @return true if latitude is valid (between -90 and 90 degrees and finite)
 */
template <typename T>
bool valid_lat(const T &lat) {
    return std::abs(lat) <= 90.0 && std::isfinite(lat);
}

/**
 * @brief Check if value is valid longitude
 * @tparam T Numeric type (float, double, etc)
 * @param lon Longitude value to check
 * @return true if longitude is valid (between -360 and 360 degrees and finite)
 */
template <typename T>
bool valid_lon(const T &lon) {
    return std::abs(lon) <= 360.0 && std::isfinite(lon);
}

/**
 * @brief Allocate and zero-initialize memory for an array
 * @tparam T Type of array elements
 * @param number_of_bins Number of elements to allocate
 * @return Pointer to allocated memory, or nullptr if allocation fails
 */
template <typename T>
T *calloc(size_t number_of_bins) {
    return static_cast<T *>(calloc(number_of_bins, sizeof(T)));
}

/**
 * @brief Allocate and zero-initialize memory for an array via pointer
 * @tparam T Type of array elements
 * @param address Pointer to store allocated memory address
 * @param number_of_bins Number of elements to allocate
 */
template <typename T>
void calloc(T **address, size_t number_of_bins) {
    *address = static_cast<T *>(calloc(number_of_bins, sizeof(T)));
}

/**
 * @brief Reallocate memory for an array
 * @tparam T Type of array elements
 * @param address Pointer to memory to reallocate
 * @param number_of_bins New number of elements
 * @return Pointer to reallocated memory
 */
template <typename T>
T *realloc(T **address, size_t number_of_bins) {
    return static_cast<T *>(realloc(*address, number_of_bins * sizeof(T)));
}

/**
 * @brief Resize netCDF buffer based on count dimensions
 * @tparam T Type of buffer elements
 * @param data Vector to resize
 * @param count Vector containing dimension sizes
 */
template <class T>
void resize_nc_buffer(std::vector<T> &data, std::vector<size_t> &count) {
    size_t total_size = 1;
    for (const auto &dim : count)
        total_size *= dim;
    data.resize(total_size);
}

/**
 * @brief Clear vector memory by swapping with empty vector
 * @tparam T Type of vector elements
 * @param vec Vector to clear
 */
template <class T>
void free_vector(std::vector<T>& vec) {
    std::vector<T>().swap(vec);
}

/**
 * @brief Get attribute value as float
 * @tparam T Source type of attribute
 * @tparam V Type of netCDF attribute (defaults to netCDF::NcVarAtt)
 * @param attr NetCDF attribute to read
 * @return Attribute value converted to float
 */
template <class T, class V = netCDF::NcVarAtt>
float get_attr(V &attr) {
    T val;
    attr.getValues(&val);
    return static_cast<float>(val);
}

/**
 * @brief Round vector elements to new type
 * @tparam T Source type of vector elements
 * @tparam F Destination type for rounded elements
 * @param inp Input vector
 * @return Vector with elements rounded and converted to type F
 */
template <class T, class F>
std::vector<F> vector_rounder(const std::vector<T> &inp) {
    std::vector<F> result(inp.size());
    std::transform(inp.begin(), inp.end(), result.begin(),
                   [](const auto &c) { return static_cast<F>(std::round(c)); });
    return result;
}

/**
 * @brief Convert std::string to char*
 * @param data String to convert
 * @return Pointer to converted string in static buffer
 */
char* convert_string(const std::string& data);

#endif
