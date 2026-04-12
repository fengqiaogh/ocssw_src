#ifndef OCSSW_UTILS_H
#define OCSSW_UTILS_H
#include <vector>
#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <cmath>
#include "genutils.h"
#include <array>
#include <unordered_map>
#include <variant>
constexpr double EarthRadius = 6371.0;
constexpr double semi_axis_a = 6378.137;
constexpr double f = 1.0 / 298.257223563;
constexpr double semi_axis_b = (1 - f) * semi_axis_a;
constexpr double e2 = 2 * f - f * f;
size_t count_valid(const std::vector<int>& indexes);
typedef std::unordered_map<std::string, std::variant<std::string, float, double, int>> InputAttributes;
#define NCPP_ERROR(...)                                                                \
    try {                                                                              \
        __VA_ARGS__;                                                                   \
    } catch (const netCDF::exceptions::NcException& e) {                               \
        fprintf(stderr, "-E- : %s:%d. NetCDF error %s", __FILE__, __LINE__, e.what()); \
        exit(EXIT_FAILURE);                                                            \
    }
#define NC_ERROR(e)                                                                                  \
    {                                                                                                \
        int ret = e;                                                                                 \
        if (ret != NC_NOERR) {                                                                       \
            fprintf(stderr, "-E-: %s:%d. NetCDF error: %s\n", __FILE__, __LINE__, nc_strerror(ret)); \
            exit(EXIT_FAILURE);                                                                      \
        }                                                                                            \
    }
class StdoutRedirector {
   private:
    int saved_stdout;
    int devnull;
    bool is_redirected;

   public:
    StdoutRedirector();
    void redirect();

    void restore();

    ~StdoutRedirector() {
        restore();  // Ensure cleanup happens
    }
};

template <class T>
std::string to_string(T&& obj) {
    if constexpr (std::is_same_v<std::string, std::decay_t<T>>) {
        return obj;
    } else {
        return std::to_string(std::forward<T>(obj));
    }
}

bool is_point_inside_cgal(std::array<std::pair<float, float>, 4>& quad, std::pair<float, float>& point);

struct Visitor {
    template <class T>
    std::string operator()(T&& arg) const {
        return to_string(arg);
    }
};

/**
 * @brief Calculate the ECEF 3D position at a given geodetic height above the ellipsoid surface.
 *
 * This function computes the 3D Cartesian coordinates of a point at a specified geodetic height (elevation)
 * above the Earth's ellipsoid. The point is located along the zenith direction (normal to the ellipsoid)
 * from the surface position at the given latitude/longitude.
 *
 * For a height of 0, this returns the position on the ellipsoid surface. For positive heights,
 * it returns a point elevated above the surface, useful for representing terrain elevation,
 * cloud top heights, or other altitude-dependent features.
 *
 * @param latitude  Geodetic latitude in degrees
 * @param longitude Geodetic longitude in degrees
 * @param height    Height above the ellipsoid surface in kilometers
 * @return A Point representing the 3D ECEF position in kilometers relative to Earth's center
 */
template <typename T>
T geodetic_to_ECEF(double latitude, double longitude, double height) {
    // precompute cosine and sine of the latitude
    double sin_lat = sin(latitude * OEL_DEGRAD);
    double cos_lat = cos(latitude * OEL_DEGRAD);
    // compute prime vertical radius of curvature
    double N = semi_axis_a / sqrt(1 - e2 * sin_lat * sin_lat);
    // compute x, y and z
    double x = (N + height) * cos(longitude * OEL_DEGRAD) * cos_lat;
    double y = (N + height) * sin(longitude * OEL_DEGRAD) * cos_lat;
    double z = (N * (1 - e2) + height) * sin_lat;
    return {x, y, z};
}

std::vector<std::tuple<std::string, double, double>> get_time_coverage(const std::vector<std::string>& files);

double dot_product(const std::array<double, 3>& v1, const std::array<double, 3>& v2);

std::array<double, 3> cross_product(const std::array<double, 3>& a, const std::array<double, 3>& b);
std::array<double, 3> add_vectors(const std::array<double, 3>& a, const std::array<double, 3>& b);

std::array<double, 3> operator+(const std::array<double, 3>& a, const std::array<double, 3>& b);
std::array<double, 3> operator*(double alpha, const std::array<double, 3>& a);
double operator*(const std::array<double, 3>& a, const std::array<double, 3>& b);
bool operator==(const std::array<double, 3>& a, const std::array<double, 3>& b);
std::array<double, 3> operator%(const std::array<double, 3>& a, const std::array<double, 3>& b);

std::array<double, 3> multiply_vector(double alpha, const std::array<double, 3>& b);
std::array<double, 3> rotate_around_vector(double cos_theta, double sin_theta,
                                           const std::array<double, 3>& axis,
                                           const std::array<double, 3>& vector);

std::array<double, 3> rotation_alignment(const std::array<double, 3>& normal_vector,
                                        const std::array<double, 3>& vector);

std::string expand_env_variable_path(const std::string& file_env_path);   
bool valid_coordinates(float lat, float lon);
bool invert3x3_row_major(double inp[3][3], double inv[3][3]);                                     
#endif  // OCSSW_UTILS_H
