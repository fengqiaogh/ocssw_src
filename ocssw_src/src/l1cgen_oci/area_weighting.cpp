#include "area_weighting.h"
#include "utils.h"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

namespace bg = boost::geometry;

typedef bg::model::d2::point_xy<double> point_t;
typedef bg::model::polygon<point_t> polygon_t;

double find_intersection_area(polygon_t& poly1, polygon_t& poly2) {
    if (bg::disjoint(poly1, poly2)) {
        return 0.0;
    }

    std::vector<polygon_t> intersection_result;

    bg::intersection(poly1, poly2, intersection_result);

    double area = 0.0;
    for (const auto& poly : intersection_result) {
        area += bg::area(poly);
    }
    return area;
}


template <size_t N>
std::vector<int> create_vertices(std::vector<float>& lat, std::vector<float>& lon,
                                 const std::array<std::pair<int, int>, N>& coords, int ipixel, int iline,
                                 int pixels, int index) {
    int last_pixel = -1;
    int last_line = -1;
    int first_pixel = -1;
    int first_line = -1;
    std::vector<int> indexes_to_insert;
    // we start from the left middle, pixel - 1, line
    for (const auto& [p, l] : coords) {
        int index_w = l * pixels + p;
        if (valid_coordinates(lat[index_w], lon[index_w])) {
            if (last_pixel != -1) {
                int dy_2 = l - iline;
                int dx_2 = p - ipixel;
                int dy_1 = last_line - iline;
                int dx_1 = last_pixel - ipixel;
                int cross = dy_2 * dx_1 - dy_1 * dx_2;
                if (cross > 0) {
                    indexes_to_insert.emplace_back(index);
                }
            } else {
                first_pixel = ipixel;
                first_line = iline;
            }
            last_pixel = ipixel;
            last_line = iline;
            indexes_to_insert.emplace_back(index_w);
        }
    }
    if (first_line != -1 && indexes_to_insert.size() > 1) {
        int dy_2 = first_line - iline;
        int dx_2 = first_pixel - ipixel;
        int dy_1 = last_line - iline;
        int dx_1 = last_pixel - ipixel;
        int cross = dy_2 * dx_1 - dy_1 * dx_2;
        if (cross > 0) {
            indexes_to_insert.emplace_back(index);
        }
    }

    return indexes_to_insert;
}

std::vector<std::vector<int>> build_l1b_grid(std::vector<float>& lat, std::vector<float>& lon, int lines,
                                             int pixels) {
    std::vector<std::vector<int>> results(lines * pixels);
    assert(lat.size() == (size_t)(pixels * lines));
    assert(lon.size() == (size_t)(pixels * lines));
    for (int iline = 0; iline < lines; iline++)
        for (int ipixel = 0; ipixel < pixels; ipixel++) {
            int index = iline * pixels + ipixel;
            if (!valid_coordinates(lat[index], lon[index])) {
                continue;
            }
            std::vector<int> indexes_to_insert;
            if (iline > 0 && iline < lines - 1)
                if (ipixel > 0 && ipixel < pixels - 1) {
                    std::array<std::pair<int, int>, 8> coords = {
                        std::pair<int, int>{ipixel - 1, iline - 1},
                        {ipixel - 1, iline},
                        {ipixel - 1, iline + 1},
                        {ipixel, iline + 1},
                        {ipixel + 1, iline + 1},
                        {ipixel + 1, iline},
                        {ipixel + 1, iline - 1},
                        {ipixel, iline - 1},
                    };
                    indexes_to_insert = create_vertices(lat, lon, coords, ipixel, iline, pixels, index);
                }
            if (iline == 0 && ipixel > 0 && ipixel < pixels - 1) {
                std::array<std::pair<int, int>, 5> coords = {
                    std::pair<int, int>{ipixel - 1, iline},
                    {ipixel - 1, iline + 1},
                    {ipixel, iline + 1},
                    {ipixel + 1, iline + 1},
                    {ipixel + 1, iline},

                };
                indexes_to_insert = create_vertices(lat, lon, coords, ipixel, iline, pixels, index);
            }

            if (iline == lines - 1 && ipixel > 0 && ipixel < pixels - 1) {
                std::array<std::pair<int, int>, 5> coords = {
                    std::pair<int, int>{ipixel - 1, iline - 1},
                    {ipixel - 1, iline},
                    {ipixel + 1, iline},
                    {ipixel + 1, iline - 1},
                    {ipixel, iline - 1},
                };
                indexes_to_insert = create_vertices(lat, lon, coords, ipixel, iline, pixels, index);
            }
            if (iline > 0 && iline < lines - 1 && ipixel == 0) {
                std::array<std::pair<int, int>, 5> coords = {
                    std::pair<int, int>{ipixel, iline + 1},
                    {ipixel + 1, iline + 1},
                    {ipixel + 1, iline},
                    {ipixel + 1, iline - 1},
                    {ipixel, iline - 1},
                };
                indexes_to_insert = create_vertices(lat, lon, coords, ipixel, iline, pixels, index);
            }
            if (iline > 0 && iline < lines - 1 && ipixel == pixels - 1) {
                std::array<std::pair<int, int>, 5> coords = {
                    std::pair<int, int>{ipixel - 1, iline - 1},
                    {ipixel - 1, iline},
                    {ipixel - 1, iline + 1},
                    {ipixel, iline + 1},
                    {ipixel, iline - 1},
                };
                indexes_to_insert = create_vertices(lat, lon, coords, ipixel, iline, pixels, index);
            }
            if (iline == 0 && ipixel == 0) {
                std::array<std::pair<int, int>, 4> coords = {
                    std::pair<int, int>{ipixel, iline},
                    {ipixel, iline + 1},
                    {ipixel + 1, iline + 1},
                    {ipixel + 1, iline},

                };
                indexes_to_insert = create_vertices(lat, lon, coords, ipixel, iline, pixels, index);
            }

            if (iline == 0 && ipixel == pixels - 1) {
                std::array<std::pair<int, int>, 4> coords = {
                    std::pair<int, int>{ipixel - 1, iline},
                    {ipixel - 1, iline + 1},
                    {ipixel, iline + 1},
                    {ipixel, iline},
                };
                indexes_to_insert = create_vertices(lat, lon, coords, ipixel, iline, pixels, index);
            }

            if (iline == lines - 1 && ipixel == 0) {
                std::array<std::pair<int, int>, 4> coords = {
                    std::pair<int, int>{ipixel, iline},
                    {ipixel + 1, iline},
                    {ipixel + 1, iline - 1},
                    {ipixel, iline - 1},
                };
                indexes_to_insert = create_vertices(lat, lon, coords, ipixel, iline, pixels, index);
            }

            if (iline == lines - 1 && ipixel == pixels - 1) {
                std::array<std::pair<int, int>, 4> coords = {
                    std::pair<int, int>{ipixel - 1, iline - 1},
                    {ipixel - 1, iline},
                    {ipixel, iline},
                    {ipixel, iline - 1},
                };
                indexes_to_insert = create_vertices(lat, lon, coords, ipixel, iline, pixels, index);
            }
            // more cases left
            results[index] = indexes_to_insert;
        }

    return results;
}

std::vector<std::array<int, 4>> build_l1c_grid(std::vector<float>& lat, std::vector<float>& lon, int lines,
                                               int pixels) {
    std::vector<std::array<int, 4>> results(lines * pixels);
    assert(lat.size() == (size_t)(pixels * lines));
    assert(lon.size() == (size_t)(pixels * lines));
    for (int iline = 0; iline < lines; iline++)
        for (int ipixel = 0; ipixel < pixels; ipixel++) {
            int index = iline * pixels + ipixel;
            int iline_up = std::min(iline + 1, lines - 1);
            int iline_down = std::max(iline - 1, 0);
            int ipixel_left = std::max(ipixel - 1, 0);
            int ipixel_right = std::min(ipixel + 1, pixels - 1);
            int index_left_down = iline_down * pixels + ipixel_left;
            int index_left_up = iline_up * pixels + ipixel_left;
            int index_right_up = iline_up * pixels + ipixel_right;
            int index_right_down = iline_down * pixels + ipixel_right;
            results[index] = {index_left_down, index_left_up, index_right_up, index_right_down};
        }

    return results;
}

namespace {
std::vector<double> l1b_areas;
}

std::vector<std::vector<std::pair<int, double>>> get_area_weights(
    std::vector<float>& lat_c, std::vector<float>& lon_c, std::vector<float>& lat_b,
    std::vector<float>& lon_b, size_t lines_c, size_t pixels_c, size_t lines_b, size_t pixels_b,
    std::vector<int>& index_search_results, bool compute_l1b_areas) {
    std::vector<std::vector<std::pair<int, double>>> results(lines_b * pixels_b);
    if (compute_l1b_areas) {
        l1b_areas.resize(lines_b * pixels_b);
    }
    auto grid_l1c = build_l1c_grid(lat_c, lon_c, lines_c, pixels_c);
    auto grid_l1b = build_l1b_grid(lat_b, lon_b, lines_b, pixels_b);
    for (size_t ip = 0; ip < lines_b * pixels_b; ip++) {
        if (compute_l1b_areas) {
            std::array<double, 3> center = geodetic_to_ECEF<std::array<double, 3>>(lat_b[ip], lon_b[ip], 0);

            auto vertices = grid_l1b[ip];
            std::vector<point_t> poly_vector;
            polygon_t poly;
            for (int iv : vertices) {
                std::array<double, 3> vert =
                    geodetic_to_ECEF<std::array<double, 3>>(lat_b[iv], lon_b[iv], 0.0);
                vert = 0.5 * (vert + center);
                auto aligned = rotation_alignment(center, vert);
                poly_vector.emplace_back(point_t(aligned[0], aligned[1]));
            }
            bg::assign_points(poly, poly_vector);
            bg::correct(poly);
            l1b_areas[ip] = bg::area(poly);
        }

        int index_l1c = index_search_results[ip];
        if (index_l1c == -1)
            continue;
        std::vector<std::pair<int, double>> area_index;
        int i_l1c_line = index_l1c / pixels_c;
        int i_l1c_pixel = index_l1c % pixels_c;
        int i_min_line = std::max(0, i_l1c_line - 1);
        int i_max_line = std::min((int)lines_c - 1, i_l1c_line + 1);
        int i_min_pixel = std::max(0, i_l1c_pixel - 1);
        int i_max_pixel = std::min((int)pixels_c - 1, i_l1c_pixel + 1);
        for (int ii_l = i_min_line; ii_l <= i_max_line; ii_l++)
            for (int ii_p = i_min_pixel; ii_p <= i_max_pixel; ii_p++) {
                int index_l1c_window = ii_l * pixels_c + ii_p;
                std::array<double, 3> center_l1b =
                    geodetic_to_ECEF<std::array<double, 3>>(lat_b[ip], lon_b[ip], 0);
                std::array<double, 3> center_l1c = geodetic_to_ECEF<std::array<double, 3>>(
                    lat_c[index_l1c_window], lon_c[index_l1c_window], 0);
                auto vertices = grid_l1b[ip];
                std::vector<point_t> poly_vector;
                polygon_t poly;
                for (int iv : vertices) {
                    std::array<double, 3> vert =
                        geodetic_to_ECEF<std::array<double, 3>>(lat_b[iv], lon_b[iv], 0.0);
                    vert = 0.5 * (vert + center_l1b);
                    auto aligned = rotation_alignment(center_l1c, vert);
                    poly_vector.emplace_back(point_t(aligned[0], aligned[1]));
                }
                bg::assign_points(poly, poly_vector);
                bg::correct(poly);

                auto vertices_l1c = grid_l1c[index_l1c_window];
                std::vector<point_t> poly_vector_l1c;
                polygon_t poly_l1c;
                for (int iv : vertices_l1c) {
                    std::array<double, 3> vert =
                        geodetic_to_ECEF<std::array<double, 3>>(lat_c[iv], lon_c[iv], 0.0);
                    vert = 0.5 * (vert + center_l1c);
                    auto aligned = rotation_alignment(center_l1c, vert);
                    poly_vector_l1c.emplace_back(point_t(aligned[0], aligned[1]));
                }

                bg::assign_points(poly_l1c, poly_vector_l1c);
                bg::correct(poly_l1c);
                double area_intersection = find_intersection_area(poly_l1c, poly);
                if (area_intersection != 0) {
                    area_index.emplace_back(std::pair<int, double>{index_l1c_window, area_intersection});
                }
                if (index_l1c_window == index_l1c && area_intersection == 0) {
                    index_search_results[ip] = -1;
                }
            }
        results[ip] = area_index;
    }
    return results;
}

double get_l1b_area(int index) {
    if (l1b_areas.empty()) {
        exit(1);
    }
    return l1b_areas[index];
}