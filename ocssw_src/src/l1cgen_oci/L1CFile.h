#ifndef OCSSW_L1CFILE_H
#define OCSSW_L1CFILE_H
#include "scaled_nc_file.hpp"
#include <filesystem>
#include "bin_l1b.h"
#include "utils.h"
#include "algorithm"
#include <algorithm>
#include "L1BFile.h"
constexpr float fillValue = BAD_FLT;
constexpr short fillValueInt = BAD_INT;
constexpr float scale_factor = 0.01f;
constexpr float add_offset = 0.0f;
constexpr float max_i = 999.f;
constexpr float min_i = 0.f;
constexpr float max_i_stdev = 800.f;
constexpr float min_i_stdev = 0.f;
constexpr short max_azimuth = 180 * 100;
constexpr short min_azimuth = -max_azimuth;
constexpr short max_zenith = 180 * 100;
constexpr short min_zenith = 0;
constexpr short max_scattering = 180 * 100;
constexpr short min_scattering = 0;
constexpr short max_number_of_observations = 999;
constexpr short min_number_of_observations = 0;
constexpr double offset_time_max = 200;
constexpr double offset_time_min = -200;
constexpr double nadir_time_max = 172800.;
constexpr double nadir_time_min = 0;
constexpr u_char min_uchar = 0;
constexpr u_char max_uchar = 10;
constexpr u_char fill_value_char = 255;
constexpr short max_height = 10'000;
constexpr short min_height = -10'000;
constexpr short min_height_stdev = 0;
constexpr float min_wavelength = 0;
constexpr float max_wavelength = 10000;
constexpr float min_f0 = 0;
constexpr float max_f0 = 4000;
constexpr float min_width = 0;
constexpr float max_width = 100;
constexpr float oci_bandpass = 5.0;
constexpr float min_tilt = -90.0;
constexpr float max_tilt = 90.0;
class L1CFile {
   private:
    netCDF::NcFile file;
    std::vector<float> i_data;
    std::vector<float> i_stdev_data;
    size_t n_lines{};
    size_t n_pixels{};
    std::vector<float> latitude_data;
    std::vector<float> longitude_data;
    std::vector<short> height_data;
    std::vector<float> height_stdev{};
    std::vector<double> nadir_view_time_data;
    std::vector<float> sensor_zenith_angles_data, solar_zenith_angles_data, sensor_azimuth_angles_data,
        solar_azimuth_angles_data, scattering_angles_data;
    netCDF::NcVar lat_var, lon_var, height_var, height_stdev_var, sensor_zenith_angles_var,
        solar_zenith_angles_var, sensor_azimuth_angles_var, solar_azimuth_angles_var, scattering_angles_var;
    netCDF::NcVar aggregate_height_var, aggregate_height_stdev_var;
    netCDF::NcVar i_var;
    netCDF::NcVar i_stdev_var;
    netCDF::NcVar number_of_observations_var;
    netCDF::NcVar nadir_view_time_var;
    netCDF::NcVar view_time_offsets_var;
    std::vector<short> number_of_observations;
    std::vector<short> number_of_observations_band;
    size_t n_bands{};
    std::string o_file;
    std::vector<double> view_time_offsets_data;
    std::string units_since{};
    netCDF::NcVar qc;
    std::vector<u_char> quality_flags{};
    netCDF::NcVar sensor_view_angle_var;
    netCDF::NcVar scan_quality_var;
    double offset_time{};
    std::vector<float> tilt;
    std::vector<u_char> scan_quality_flags;
    std::vector<size_t> number_of_observations_line;
    std::vector<float> area_total{};
    std::vector<float> area_total_per_band{};
    std::vector<float> area_second_momentum{};
    bool cloud_correction = false;
    std::vector<float> aggregated_height{};
    std::vector<float> aggregated_height_stdev{};

   public:
    L1CFile() = default;
    friend int bin_l1b_file(L1BFile& l1bfile, L1CFile& l1cfile);
    explicit L1CFile(const std::string& output_file_name,
                     const std::multimap<std::string, netCDF::NcGroupAtt>& l1c_global_attributes,
                     const L1BFile& l1bfile, size_t n_lines, size_t n_pixels,
                     const std::vector<float>& latitude, const std::vector<float>& longitude,
                     const std::vector<short>& height, const std::vector<double>& nadir_view_time,
                     const InputAttributes& input_attributes);
    void close();
    virtual ~L1CFile() {
        close();
    }
};

#endif  // OCSSW_L1CFILE_H
