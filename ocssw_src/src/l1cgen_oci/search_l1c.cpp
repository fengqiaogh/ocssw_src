#include "search_l1c.h"
#include "utils.h"
using Point = bg::model::d2::point_xy<double>;
using Polygon = bg::model::polygon<Point>;
using LineString = bg::model::linestring<Point>;

/**
 * @brief Convert geodetic coordinates to a unit 3D ECEF vector.
 *
 *
 * @param lat Latitude in degrees 
 * @param lon Longitude in degrees
 *
 * @return A 3D unit vector in ECEF coordinates representing the surface position
 *         at the given latitude and longitude.
 *
 */
point_3d get_norm_vector(double lat, double lon) {
    return geodetic_to_ECEF<point_3d>(lat, lon, 0.0);
}

/**
 * @brief Construct a SearchL1C spatial index for efficient L1C grid point queries.
 *
 * @param lat_data_c Vector of latitude values for the L1C grid (in degrees).
 * @param lon_data_c Vector of longitude values for the L1C grid (in degrees).
 *        Both vectors must have the same size (n_lines × n_pixels elements).
 * @param n_lines Number of scan lines in the L1C grid.
 * @param n_pixels Number of pixels per scan line in the L1C grid.

 */
SearchL1C::SearchL1C(const std::vector<float>& lat_data_c, const std::vector<float>& lon_data_c,
                     size_t n_lines, size_t n_pixels)
    : n_lines(n_lines), n_pixels(n_pixels) {
    std::vector<value_t> setl1c;
    for (std::size_t i = 0; i < lat_data_c.size(); ++i) {
        setl1c.emplace_back(get_norm_vector(lat_data_c[i], lon_data_c[i]), i);
    }
    lats_c = std::vector<double>(lat_data_c.begin(), lat_data_c.end());
    lons_c = std::vector<double>(lon_data_c.begin(), lon_data_c.end());
    rtree_setl1c = bgi::rtree<value_t, bgi::quadratic<16>>(setl1c.begin(), setl1c.end());
    // find max distance
    for (size_t i_line = 0; i_line < n_lines - 1; ++i_line) {
        for (size_t i_pixel = 0; i_pixel < n_pixels - 1; ++i_pixel) {
            size_t index = i_line * n_pixels + i_pixel;
            size_t index_right = index + 1;
            size_t index_down = index + n_pixels;
            size_t index_down_right = index_down + 1;
            double dist_01 =
                haversineDistance(lats_c[index], lons_c[index], lats_c[index_right], lons_c[index_right]);
            double dist_02 =
                haversineDistance(lats_c[index], lons_c[index], lats_c[index_down], lons_c[index_down]);
            double dist_03 = haversineDistance(lats_c[index], lons_c[index], lats_c[index_down_right],
                                               lons_c[index_down_right]);
            double dist_12 = haversineDistance(lats_c[index_down], lons_c[index_down], lats_c[index_right],
                                               lons_c[index_right]);
            double dist_13 = haversineDistance(lats_c[index_down], lons_c[index_down],
                                               lats_c[index_down_right], lons_c[index_down_right]);
            double dist_23 = haversineDistance(lats_c[index_right], lons_c[index_right],
                                               lats_c[index_down_right], lons_c[index_down_right]);
            max_dist = std::max(std::max(dist_01, dist_02), std::max(dist_12, dist_13));
            max_dist = std::max(max_dist, std::max(dist_23, dist_03));
        }
    }
}

/**
 * @brief Query the spatial index to find nearest L1C grid points for given lat/lon coordinates.
 *
 *
 * @param lat A vector of latitude values 
 * @param lon A vector of longitude values 
 *
 * @return A vector of indices into the L1C grid. Each index corresponds to the nearest
 *         grid point for the corresponding input coordinate. Returns -1 for invalid
 *         coordinates (out of range) or for points that lie outside the L1C grid bounds.
 */
std::vector<int> SearchL1C::query_lat_lon(const std::vector<float>& lat,
                                          const std::vector<float>& lon) const {
    std::vector<int> nearest_indices;

    nearest_indices.reserve(lat.size());
    for (size_t i = 0; i < lat.size(); ++i) {
        // check if lat[i], lon[i] is valid
        if (std::abs(lat[i]) > 90.0f || std::abs(lon[i]) > 180.0f) {
            nearest_indices.emplace_back(-1);
            continue;
        }
        point_3d q = get_norm_vector(lat[i], lon[i]);
        std::vector<value_t> result;
        rtree_setl1c.query(bgi::nearest(q, 1), std::back_inserter(result));
        nearest_indices.emplace_back(result.front().second);
    }

    // go over all indees
    for (size_t i = 0; i < nearest_indices.size(); ++i) {
        // get line and pixel in L1C
        int index_l1c = nearest_indices[i];
        if (index_l1c == -1) {
            continue;
        }
        size_t line = index_l1c / n_pixels;
        size_t pixel = index_l1c % n_pixels;
        // NN algorithm doesn't tell if a point inside l1c grid.
        // if line == 0 or line == n_lines -1 or pixel == 0 or pixel == n_pixels -1, check that
        // haversineDistance between the l1b and l1c point is less than max_dist, in km.
        if (line == 0 || line == n_lines - 1 || pixel == 0 || pixel == n_pixels - 1) {
            double dist = haversineDistance(lat[i], lon[i], lats_c[index_l1c], lons_c[index_l1c]);
            if (dist > max_dist) {
                nearest_indices[i] = -1;
                continue;
            }
            // exclude pixels  that outside of the grid. the haversine distance check implemented above is
            // fast and eliminates most of the pixels outside for points close to the border, check if they
            // are actually inside build 4 points move left, move up
            size_t line_up = std::max(0, (int)line - 1);
            size_t pixel_left = std::max(0, (int)pixel - 1);
            size_t index_up_left = line_up * n_pixels + pixel_left;
            // move right, move down
            size_t line_down = std::min((int)n_lines - 1, (int)line + 1);
            size_t pixel_right = std::min((int)n_pixels - 1, (int)pixel + 1);
            size_t index_down_right = line_down * n_pixels + pixel_right;
            size_t index_up_right = line_up * n_pixels + pixel_right;
            size_t index_down_left = line_down * n_pixels + pixel_left;
            std::array<std::pair<float, float>, 4> quad = {
                std::pair{lats_c[index_up_left], lons_c[index_up_left]},
                {lats_c[index_down_right], lons_c[index_down_right]},
                {lats_c[index_down_left], lons_c[index_down_left]},
                {lats_c[index_up_right], lons_c[index_up_right]}};
            std::pair<float, float> point = {lat[i], lon[i]};
            if (!is_point_inside_cgal(quad, point)) {
                nearest_indices[i] = -1;
            }
        }
    }
    return nearest_indices;
}

/**
 * @brief Calculate the great-circle distance between two points on a sphere using the Haversine formula.
 *
 * @param lat1 Latitude of the first point in degrees 
 * @param lon1 Longitude of the first point in degrees 
 * @param lat2 Latitude of the second point in degrees 
 * @param lon2 Longitude of the second point in degrees 
 * @param radius Radius of the sphere in the desired distance units
 */
double haversineDistance(double lat1, double lon1, double lat2, double lon2, double radius) {
    // Convert degrees to radians
    double lat1_rad = lat1 * OEL_DEGRAD;
    double lon1_rad = lon1 * OEL_DEGRAD;
    double lat2_rad = lat2 * OEL_DEGRAD;
    double lon2_rad = lon2 * OEL_DEGRAD;

    // Difference in coordinates
    double dlat = lat2_rad - lat1_rad;
    double dlon = lon2_rad - lon1_rad;

    // Haversine formula
    double a = sin(dlat / 2) * sin(dlat / 2) + cos(lat1_rad) * cos(lat2_rad) * sin(dlon / 2) * sin(dlon / 2);

    // Prevent sqrt domain error from floating point inaccuracies
    if (a > 1.0)
        a = 1.0;

    double c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));

    return radius * c;
}