#ifndef GEO_SPATIAL_H
#define GEO_SPATIAL_H

#include <iostream>
#include <string>
#include <netcdf>
#include <vector>
#include <unordered_map>
#include <cmath>

class Geospatialbounds {
private:
    std::unordered_map<std::string, std::vector<float>> scan_line_vars;
    std::vector<std::string> possible_sozlen_names = {"solzen", "solar_zenith_angle", "solz", "Solar Zenith Angle",
                                                      "csol_z", "Center Solar Zenith Angle", "center_solar_zenith_angle"};
    float min_lon = NAN;
    float max_lon = NAN;
    std::string geospatial_bounds_wkt{};
    std::string platform{};
    std::string instrument{};
    std::string day_night_flag{};
    std::string time_coverage_start{};
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

    const float *get_elat();

    const float *get_slat();

    const std::unordered_map<std::string, std::vector<float>> &get_bounds();

    std::pair<std::vector<float>, std::vector<float>> get_gring();

    std::string get_bounds_wkt();

    float get_min_lon();

    float get_max_lon();

    const std::string &get_platform() const {
        return platform;
    };

    const std::string &get_instrument() const {
        return instrument;
    };

    const std::string &get_day_night_flag() const {
        return day_night_flag;
    };

    const std::string &get_time_coverage_start() const {
        return time_coverage_start;
    };

    static void allocate_var(std::vector<float> &, const netCDF::NcVar &var);

    static void allocate_attr(std::vector<float> &, const netCDF::NcAtt &att);

    static std::pair<std::vector<float>, std::vector<float>> calc_gring(const std::vector<float> &,
                                                                        const std::vector<float> &,
                                                                        const std::vector<float> &,
                                                                        const std::vector<float> &,
                                                                        const std::vector<float> &);

    static std::string calc_gring(const std::pair<std::vector<float>, std::vector<float>> &);

    static void calc_scan_line(const std::vector<float> &, const std::vector<float> &, size_t, size_t,
                               std::unordered_map<std::string, std::vector<float>> &);
};

#endif
