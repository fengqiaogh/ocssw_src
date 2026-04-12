#include "surface_intersector.h"
#include "surface_correction.h"
namespace srf {
/**
 * @brief Correct  latitude, longitude, and height for surface intersections.
 *
 *
 * @param latitude    latitude values in degrees
 *                             Updated with corrected values where intersections with the surface occur.
 * @param longitude   longitude values in degrees
 *                             Updated with corrected values where intersection swith the surface occur.
 * @param height      elevation in kilometers. Updated with corrected
 *                             values where intersections with the surface occur.
 * @param sensor_zenith_angle   Sensor zenith angles
 * @param sensor_azimuth_angle  Sensor azimuth angles
 * @param latitude_grid   The reference surface grid latitude value
 * @param longitude_grid  The reference surface grid longitude values
 * @param height_grid     The reference surface grid elevation in kilometers.
 * @param mask_grid       Binary mask for valid reference surface grid points.
 * @param lines        Number of lines in lat/lon grid
 * @param pixels       Number of pixels in lat/lon grid
 * @param lines_grid   Number of lines in the reference surface grid.
 * @param pixels_grid  Number of pixels per line in the reference surface grid.

 */
void surface_correction(std::vector<float>& latitude, std::vector<float>& longitude,
                        std::vector<float>& height, std::vector<float>& sensor_zenith_angle,
                        std::vector<float>& sensor_azimuth_angle, const std::vector<float>& latitude_grid,
                        const std::vector<float>& longitude_grid, const std::vector<float>& height_grid,
                        const std::vector<int8_t>& mask_grid, size_t lines, size_t pixels, size_t lines_grid,
                        size_t pixels_grid) { // height must be in kilometers
    assert(latitude.size() == lines * pixels);
    assert(longitude.size() == lines * pixels);
    assert(height.size() == lines * pixels);
    assert(sensor_zenith_angle.size() == lines * pixels);
    assert(sensor_azimuth_angle.size() == lines * pixels);
    assert(latitude_grid.size() == lines_grid * pixels_grid);
    assert(longitude_grid.size() == lines_grid * pixels_grid);
    assert(height_grid.size() == lines_grid * pixels_grid);
    assert(mask_grid.size() == lines_grid * pixels_grid);

    std::vector<Point> points_grid(lines_grid * pixels_grid);
    for (size_t i = 0; i < lines_grid * pixels_grid; i++) {
        points_grid[i] = elevated_position(latitude_grid[i], longitude_grid[i], height_grid[i]);
    }
    SurfaceIntersector intersector(points_grid, lines_grid, pixels_grid, mask_grid);
    std::vector<Point> points_surface(lines * pixels);
    std::vector<Point> points_satellite(lines * pixels);
    for (size_t i = 0; i < lines * pixels; i++) {
        points_surface[i] = elevated_position(latitude[i], longitude[i], height[i]);
        points_satellite[i] = sensor_position(sensor_zenith_angle[i], sensor_azimuth_angle[i], latitude[i],
                                              longitude[i], height[i]);
    }

    for (size_t i = 0; i < lines * pixels; i++) {
        if (is_bad(points_surface[i]) || is_bad(points_satellite[i]))
            continue;
        if (intersector.doesIntersect(points_surface[i], points_satellite[i])) {
            std::vector<Point> intersection_points =
                intersector.findAllIntersections(points_surface[i], points_satellite[i]);
            // find the point that closest to the satellite
            double min_dist = std::numeric_limits<double>::max();
            for (const auto& intersection : intersection_points) {
                double dist = (points_surface[i] - intersection).squared_length();
                if (dist < min_dist) {
                    min_dist = dist;
                    points_surface[i] = intersection;
                }
            }
            // updated latitude, longitude and height
            const auto [latitude_corrected, longitude_corrected, height_corrected] =
                get_lat_lon_from_xyz(points_surface[i].x(), points_surface[i].y(), points_surface[i].z());
            float sensor_zenith_angle_corrected = sensor_zenith_angle[i];
            float sensor_azimuth_angle_corrected = sensor_azimuth_angle[i];
            recalculate_sensor_angles(sensor_zenith_angle[i], sensor_azimuth_angle[i], latitude[i],
                                      longitude[i], latitude_corrected, longitude_corrected,
                                      sensor_zenith_angle_corrected, sensor_azimuth_angle_corrected);
            // santiy check here
            // if (std::abs(sensor_zenith_angle[i] - sensor_zenith_angle_corrected) > 1e-2 ||
            //     std::abs(sensor_azimuth_angle[i] - sensor_azimuth_angle_corrected) > 1e-2) {
            //     fprintf(stderr, "-E-: Fatal error: %f vs %f, %f vs %f \n", sensor_zenith_angle[i],
            //             sensor_zenith_angle_corrected, sensor_azimuth_angle[i],
            //             sensor_azimuth_angle_corrected);
            //     exit(EXIT_FAILURE);
            // }

            latitude[i] = latitude_corrected;
            longitude[i] = longitude_corrected;
            height[i] = height_corrected;
            sensor_zenith_angle[i] = sensor_zenith_angle_corrected;
            sensor_azimuth_angle[i] = sensor_azimuth_angle_corrected;
            
        }
    }
}

/**
 * @brief Recompute sensor viewing angles when the ground location shifts.
 *
 * This function takes the original sensor zenith/azimuth angles defined at an
 * original (latitude_old, longitude_old) ground point and transforms them to
 * the corresponding angles at a nearby corrected ground point
 * (latitude_new, longitude_new).
 *
 * @param sensor_zenith_angle_old  Sensor zenith angle (degrees) at the original
 *                                 ground location.
 * @param sensor_azimuth_angle_old Sensor azimuth angle (degrees) at the original
 *                                 ground location.
 * @param latitude_old            Original ground latitude (degrees).
 * @param longitude_old           Original ground longitude (degrees).
 * @param latitude_new            Corrected ground latitude (degrees).
 * @param longitude_new           Corrected ground longitude (degrees).
 * @param[out] sensor_zenith_angle_new  Computed sensor zenith angle (degrees) at the
 *                                     corrected ground location.
 * @param[out] sensor_azimuth_angle_new Computed sensor azimuth angle (degrees) at the
 *                                     corrected ground location.
 */
void recalculate_sensor_angles(float sensor_zenith_angle_old, float sensor_azimuth_angle_old,
                               float latitude_old, float longitude_old, float latitude_new,
                               float longitude_new, float& sensor_zenith_angle_new,
                               float& sensor_azimuth_angle_new) {
    const auto [vector_zenith_old, vector_azimuth_old, vector_east_old] =
        get_local_vectors(latitude_old, longitude_old);
    const auto [vector_zenith_new, vector_azimuth_new, vector_east_new] =
        get_local_vectors(latitude_new, longitude_new);
    double cos_az_old = cos(sensor_azimuth_angle_old * OEL_DEGRAD);
    double sin_az_old = sin(sensor_azimuth_angle_old * OEL_DEGRAD);
    double cos_zen_old = cos(sensor_zenith_angle_old * OEL_DEGRAD);
    double sin_zen_old = sin(sensor_zenith_angle_old * OEL_DEGRAD);
    std::array<double, 3> vector_ray{};
    double matrix[3][3], inverse[3][3];
    for (size_t i = 0; i < 3; i++) {
        vector_ray[i] = vector_zenith_old[i] * cos_zen_old +
                        vector_azimuth_old[i] * sin_zen_old * cos_az_old +
                        vector_east_old[i] * sin_zen_old * sin_az_old;
        matrix[i][0] = vector_zenith_new[i];
        matrix[i][1] = vector_azimuth_new[i];
        matrix[i][2] = vector_east_new[i];
    }
    bool status = invert3x3_row_major(matrix, inverse);
    double out[3];
    if (status) {
        for (size_t i = 0; i < 3; i++) {
            out[i] = 0;
            for (size_t j = 0; j < 3; j++) {
                out[i] += inverse[i][j] * vector_ray[j];
            }
        }
        double cos_sen_new = out[0];
        double az_new = atan2(out[2], out[1]);
        if (std::abs(cos_sen_new) <= 1) {
            sensor_zenith_angle_new = acos(cos_sen_new) * OEL_RADEG;
            sensor_azimuth_angle_new = az_new * OEL_RADEG;
        }
    }
}
}  // namespace srf