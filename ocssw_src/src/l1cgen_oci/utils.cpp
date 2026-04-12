#include "utils.h"
#include <cstdlib>
#include <array>
#include <list>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <regex>
#include <iostream>
#include <filesystem>

typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Ray_3 Ray;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;
typedef std::list<Triangle>::const_iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

size_t count_valid(const std::vector<int>& indexes) {
    size_t count_valid = 0;
    for (const auto& idx : indexes) {
        if (idx != -1) {
            count_valid++;
        }
    }
    return count_valid;
}

size_t count_valid_scan(const int* indexes, size_t n_pixels) {
    size_t count_valid = 0;
    for (size_t ip = 0; ip < n_pixels; ip++) {
        int idx = indexes[ip];
        if (idx != -1) {
            count_valid++;
        }
    }
    return count_valid;
}

// find last valid_scan
size_t count_valid_scans(const std::vector<int>& indexes, size_t n_pixels) {
    size_t nscan = indexes.size() / n_pixels;
    for (size_t iscan = 0; iscan < nscan; iscan++) {
        size_t count_valid = count_valid_scan(indexes.data() + iscan * n_pixels, n_pixels);
        if (count_valid == 0) {
            return iscan + 1;
        }
    }
    return nscan;
}
std::array<double, 3> latlon_to_xyz(const std::pair<float, float>& latlon) {
    return geodetic_to_ECEF<std::array<double, 3>>(latlon.first, latlon.second, 0.0);
};
/*
 * @brief Check if a point is inside a quadrilateral defined by four latitude/longitude pairs.
 * @param quad - array of four latitude/longitude pairs defining the quadrilateral
 * @param point - latitude/longitude pair representing the point to check
 * @return true if the point is inside the quadrilateral, false otherwise
 */
bool is_point_inside_cgal(std::array<std::pair<float, float>, 4>& quad, std::pair<float, float>& point) {
    auto p_xyz = latlon_to_xyz(point);
    std::array<std::array<double, 3>, 4> quad_xyz{};
    for (size_t i = 0; i < 4; i++) {
        quad_xyz[i] = latlon_to_xyz(quad[i]);
    }
    Point p(p_xyz[0], p_xyz[1], p_xyz[2]);
    Point origin(0.0, 0.0, 0.0);

    Point a(quad_xyz[0][0], quad_xyz[0][1], quad_xyz[0][2]);
    Point b(quad_xyz[1][0], quad_xyz[1][1], quad_xyz[1][2]);
    Point c(quad_xyz[2][0], quad_xyz[2][1], quad_xyz[2][2]);
    Point d(quad_xyz[3][0], quad_xyz[3][1], quad_xyz[3][2]);

    std::list<Triangle> triangles;
    triangles.emplace_back(a, b, c);
    triangles.emplace_back(a, b, d);
    triangles.emplace_back(a, d, c);
    // constructs AABB tree
    Tree tree(triangles.begin(), triangles.end());
    // counts #intersections
    Ray ray_query(p, origin);
    bool has_intersection = tree.do_intersect(ray_query);

    return has_intersection;
}

double dot_product(const std::array<double, 3>& v1, const std::array<double, 3>& v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

std::array<double, 3> cross_product(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
};

std::array<double, 3> add_vectors(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

std::array<double, 3> multiply_vector(double alpha, const std::array<double, 3>& a) {
    return {a[0] * alpha, a[1] * alpha, a[2] * alpha};
}

std::array<double, 3> operator+(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return add_vectors(a, b);
}

std::array<double, 3> operator*(double alpha, const std::array<double, 3>& a) {
    return multiply_vector(alpha, a);
}

double operator*(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return dot_product(a, b);
}

std::array<double, 3> operator%(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return cross_product(a, b);
}

std::array<double, 3> rotate_around_vector(double cos_theta, double sin_theta,
                                           const std::array<double, 3>& axis,
                                           const std::array<double, 3>& vector) {
    return (cos_theta * vector + (1.0 - cos_theta) * (axis * vector) * axis) + (sin_theta * axis % vector);
}

bool operator==(const std::array<double, 3>& v1, const std::array<double, 3>& v2) {
    if (std::abs(v1[0] - v2[0]) > 1e-12)
        return false;
    if (std::abs(v1[1] - v2[1]) > 1e-12)
        return false;
    if (std::abs(v1[2] - v2[2]) > 1e-12)
        return false;
    return true;
}

std::array<double, 3> rotation_alignment(const std::array<double, 3>& normal_vector,
                                         const std::array<double, 3>& vector) {
    // rotate by theta' = theta - pi/2
    double norm = std::sqrt(normal_vector * normal_vector);
    double cos_theta = normal_vector[2] / norm;  // cos(theta') =  cos(theta - pi/2) = sin(theta) = z/r
    double sin_theta = -std::sqrt(normal_vector[0] * normal_vector[0] + normal_vector[1] * normal_vector[1]) /
                       norm;  // sin(theta') =  sin(theta - pi/2) = -cos(theta)
    double norm2 = std::sqrt(normal_vector[0] * normal_vector[0] + normal_vector[1] * normal_vector[1]);
    std::array<double, 3> axis = {-normal_vector[1], normal_vector[0], 0};
    std::array<double, 3> ans = rotate_around_vector(cos_theta, sin_theta, (1 / norm2) * axis, vector);
    return ans;
}

#ifdef WINPATH
#define SEP "\\"
#define REGEX ".*%(\\w+)%.*"
#else
#define SEP "/"
#define REGEX ".*\\$(\\w+).*"
#endif

std::string expand_env_variable_path(const std::string& file_env_path) {
    std::vector<std::string> path_sep;
    boost::split(path_sep, file_env_path, boost::is_any_of(SEP));
    std::filesystem::path out_path{};
    for (const auto& el : path_sep) {
        std::string s = el;
        std::regex rgx(REGEX);
        std::smatch match;
        if (std::regex_search(s, match, rgx)) {
            std::string env_ = match[1];
            out_path /= std::getenv(env_.c_str());
        } else {
            out_path /= s;
        }
    }
    return out_path;
}

bool valid_coordinates(float lat, float lon) {
    if (std::abs(lat) <= 90 && std::abs(lon) <= 180)
        return true;
    else
        return false;
}

StdoutRedirector::StdoutRedirector() : saved_stdout(-1), devnull(-1), is_redirected(false) {
}

void StdoutRedirector::redirect() {
    if (is_redirected) {
        return;  // Already redirected
    }

    std::cout.flush();  // Flush C++ stream before redirecting

    saved_stdout = dup(STDOUT_FILENO);
    if (saved_stdout == -1) {
        std::cerr << "Failed to duplicate stdout" << std::endl;
        return;
    }

    devnull = open("/dev/null", O_WRONLY);
    if (devnull == -1) {
        std::cerr << "Failed to open /dev/null" << std::endl;
        close(saved_stdout);
        saved_stdout = -1;
        return;
    }

    dup2(devnull, STDOUT_FILENO);
    is_redirected = true;
}

void StdoutRedirector::restore() {
    if (!is_redirected) {
        return;  // Nothing to restore
    }

    if (saved_stdout != -1) {
        dup2(saved_stdout, STDOUT_FILENO);
        close(saved_stdout);
        saved_stdout = -1;
    }

    if (devnull != -1) {
        close(devnull);
        devnull = -1;
    }

    is_redirected = false;
}

bool invert3x3_row_major(double inp[3][3], double inv[3][3]) {
    // cofactors
    const double inp_00 = inp[0][0], inp_01 = inp[0][1], inp_02 = inp[0][2];
    const double inp_10 = inp[1][0], inp_11 = inp[1][1], inp_12 = inp[1][2];
    const double inp_20 = inp[2][0], inp_21 = inp[2][1], inp_22 = inp[2][2];
    const double c00 = (inp_11 * inp_22 - inp_12 * inp_21);
    const double c01 = -(inp_10 * inp_22 - inp_12 * inp_20);
    const double c02 = (inp_10 * inp_21 - inp_11 * inp_20);
    const double c10 = -(inp_01 * inp_22 - inp_02 * inp_21);
    const double c11 = (inp_00 * inp_22 - inp_02 * inp_20);
    const double c12 = -(inp_00 * inp_21 - inp_01 * inp_20);
    const double c20 = (inp_01 * inp_12 - inp_02 * inp_11);
    const double c21 = -(inp_00 * inp_12 - inp_02 * inp_10);
    const double c22 = (inp_00 * inp_11 - inp_01 * inp_10);
    const double det = inp_00 * c00 + inp_01 * c01 + inp_02 * c02;
    if (std::abs(det) < 1e-12)
        return false;
    const double invDet = double(1) / det;
    inv[0][0] = c00 * invDet;
    inv[0][1] = c10 * invDet;
    inv[0][2] = c20 * invDet;
    inv[1][0] = c01 * invDet;
    inv[1][1] = c11 * invDet;
    inv[1][2] = c21 * invDet;
    inv[2][0] = c02 * invDet;
    inv[2][1] = c12 * invDet;
    inv[2][2] = c22 * invDet;
    return true;
}