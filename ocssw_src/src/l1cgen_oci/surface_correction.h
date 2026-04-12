#ifndef OCSSW_SURFACE_CORRECTION_H
#define OCSSW_SURFACE_CORRECTION_H
namespace srf {
void surface_correction(std::vector<float>& latitude, std::vector<float>& longitude,
                        std::vector<float>& height, std::vector<float>& sensor_zenith_angle,
                        std::vector<float>& sensor_azimuth_angle, const std::vector<float>& latitude_grid,
                        const std::vector<float>& longitude_grid, const std::vector<float>& height_grid,
                        const std::vector<int8_t>& mask_grid, size_t lines, size_t pixels, size_t lines_grid,
                        size_t pixels_grid);
void recalculate_sensor_angles(float sensor_zenith_angle_old, float sensor_azimuth_angle_old,
                               float latitude_old, float longitude_old, float latitude_new,
                               float longitude_new, float& sensor_zenith_angle_new,
                               float& sensor_azimuth_angle_new);
}
#endif  // OCSSW_SURFACE_CORRECTION_H
