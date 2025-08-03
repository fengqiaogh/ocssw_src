#ifndef GEO_SPATIAL_H
#define GEO_SPATIAL_H
#include "l2_utils.hpp"
#include <iostream>
#include <string>
#include <netcdf>
#include <vector>
#include <unordered_map>
#include <cmath>

class Geospatialbounds {
   private:
    std::unordered_map<std::string, std::vector<float>> scan_line_vars;
    float min_lon = NAN;
    float max_lon = NAN;
    std::string geospatial_bounds_wkt{};
    std::string platform{};
    std::string instrument{};
    std::string day_night_flag{};
    std::string time_coverage_start{},time_coverage_end{};
    std::vector<float> gringpointlatitude{}, gringpointlongitude{};
    std::vector<float> lat, lon;

    void get_lat_lon(const netCDF::NcFile &nc_file);

    size_t number_of_lines{}, pixels_per_line{};
    std::string geospatial_lon_min;

    std::string get_day_night_flag(const netCDF::NcFile &nc_file);

   public:
    Geospatialbounds();

    explicit Geospatialbounds(const std::string &path_nc);

    explicit Geospatialbounds(const netCDF::NcFile &nc);

    /**
     * @brief Get pointer to ending latitude values
     * @return Pointer to array of ending latitude values
     */
    const float *get_elat();

    /**
     * @brief Get pointer to starting latitude values
     * @return Pointer to array of starting latitude values
     */
    const float *get_slat();

    /**
     * @brief Get map of scan line boundary values
     * @return Reference to map containing scan line boundary vectors
     */
    const std::unordered_map<std::string, std::vector<float>> &get_bounds();

    /**
     * @brief Get geographic ring point coordinates
     * @return Pair of vectors containing latitude and longitude coordinates
     */
    std::pair<std::vector<float>, std::vector<float>> get_gring();

    /**
     * @brief Get Well-Known Text representation of geographic bounds
     * @return WKT string describing geographic bounds
     */
    std::string get_bounds_wkt();

    /**
     * @brief Get minimum longitude value
     * @return Minimum longitude in degrees
     */
    float get_min_lon();

    /**
     * @brief Get maximum longitude value
     * @return Maximum longitude in degrees
     */
    float get_max_lon();

    /**
     * @brief Get platform name
     * @return Reference to platform name string
     */
    const std::string &get_platform() const {
        return platform;
    };

    /**
     * @brief Get instrument name
     * @return Reference to instrument name string
     */
    const std::string &get_instrument() const {
        return instrument;
    };

    /**
     * @brief Get day/night flag
     * @return Reference to day/night flag string
     */
    const std::string &get_day_night_flag() const {
        return day_night_flag;
    };

    /**
     * @brief Get coverage start time
     * @return Reference to time coverage start string
     */
    const std::string &get_time_coverage_start() const {
        return time_coverage_start;
    };

    /**
     * @brief Get coverage end time
     * @return Reference to time coverage end string
     */
    const std::string &get_time_coverage_end() const {
        return time_coverage_end;
    };

    /**
     * @brief Allocate and populate vector from netCDF variable
     * @param[out] Vector to allocate and populate
     * @param[in] var NetCDF variable to read from
     */
    static void allocate_var(std::vector<float> &, const netCDF::NcVar &var);

    /**
     * @brief Allocate and populate vector from netCDF attribute
     * @param[out] Vector to allocate and populate
     * @param[in] att NetCDF attribute to read from
     */
    static void allocate_attr(std::vector<float> &, const netCDF::NcAtt &att);

    /**
     * @brief Calculate geographic ring points from input vectors
     * @param[in] Five vectors containing coordinate data
     * @return Pair of vectors with calculated latitude and longitude points
     */
    static std::pair<std::vector<float>, std::vector<float>> calc_gring(const std::vector<float> &,
                                                                        const std::vector<float> &,
                                                                        const std::vector<float> &,
                                                                        const std::vector<float> &,
                                                                        const std::vector<float> &);

    /**
     * @brief Convert geographic ring points to WKT string
     * @param[in] Pair of vectors containing latitude and longitude points
     * @return WKT string representation of ring points
     */
    static std::string calc_gring(const std::pair<std::vector<float>, std::vector<float>> &);

    /**
     * @brief Calculate scan line boundaries
     * @param[in] First vector of coordinate data
     * @param[in] Second vector of coordinate data
     * @param[in] Number of lines
     * @param[in] Pixels per line
     * @param[out] Map to store calculated boundary vectors
     */
    static void calc_scan_line(const std::vector<float> &, const std::vector<float> &, size_t, size_t,
                               std::unordered_map<std::string, std::vector<float>> &);
};

#endif
