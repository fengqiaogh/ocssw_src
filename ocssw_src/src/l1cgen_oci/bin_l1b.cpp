#include "bin_l1b.h"
#include "L1BFile.h"
#include "L1CFile.h"
#include <cassert>
#include <optional>

/**
 * @brief Bins L1B files into L1C file.
 * *
 * @param[in] l1bfile L1BFile object containing the input L1B data and geometry information.
 * @param[in,out] l1cfile L1CFile object that accumulates the binned data and statistics.
 *
 * @return 0 on successful completion, non-zero error code if L1B file reading fails.
 *
 */
int bin_l1b_file(L1BFile& l1bfile, L1CFile& l1cfile) {
    size_t n_scans = l1bfile.get_l1file().nscan;
    size_t n_pixels = l1bfile.get_l1file().npix;
    size_t n_bands = l1cfile.n_bands;
    const l1str& l1rec = l1bfile.get_l1rec();
    const std::vector<int>& bin_indexes = l1bfile.bin_indexes;
    assert(bin_indexes.size() == n_pixels * n_scans);
    size_t number_of_valid = count_valid(bin_indexes);
    if (number_of_valid == 0) {
        std::cout << "No valid pixels found in L1B file " << l1bfile.ifile << std::endl;
        return 0;
    }
    // read L1B file line by line ( 50 lines per iteration
    for (int i_line = 0; i_line < (int)n_scans; ++i_line) {
        if (i_line % 50 == 0) {
            // write processing info
            std::cout << "Processing line: " << i_line << std::endl;
        }
        double scantime = l1rec.scantime - l1cfile.offset_time;
        int status = l1bfile.read_line(i_line);
        if (status != 0)
            return status;
        for (int i_pixel = 0; i_pixel < (int)n_pixels; ++i_pixel) {
            int index_l1c = bin_indexes[i_line * n_pixels + i_pixel];
            if (index_l1c == -1)
                continue;
            // bin tilt
            // compute tilt line
            int line_l1c = index_l1c / l1cfile.n_pixels;
            float& tilt = l1cfile.tilt[line_l1c];
            u_char& scan_quality_flag = l1cfile.scan_quality_flags[line_l1c];
            size_t& n_obs_line = l1cfile.number_of_observations_line[line_l1c];
            if (l1rec.tilt != BAD_FLT) {
                if (n_obs_line == 0) {
                    tilt = l1rec.tilt;
                    scan_quality_flag = l1rec.navwarn[i_pixel];
                    n_obs_line = 1;

                } else {
                    tilt = (tilt * n_obs_line + l1rec.tilt) / (n_obs_line + 1.0f);
                    n_obs_line += 1;
                    scan_quality_flag |= l1rec.navwarn[i_pixel];
                }
            }
            // bin geometry
            float senz = l1rec.senz[i_pixel];
            float sena = l1rec.sena[i_pixel];
            float sola = l1rec.sola[i_pixel];
            float solz = l1rec.solz[i_pixel];
            float height_l1b = l1rec.height[i_pixel];
            if (l1cfile.cloud_correction) {
                senz = l1bfile.senz_corrected[i_line * n_pixels + i_pixel];
                sena = l1bfile.sena_corrected[i_line * n_pixels + i_pixel];
            }
            // check if geometry is valid:
            if (senz == BAD_FLT || sena == BAD_FLT || sola == BAD_FLT || solz == BAD_FLT)
                continue;
            // check if between valid max/min ranges
            if (senz > max_zenith * scale_factor || senz < min_zenith * scale_factor)
                continue;
            if (solz > max_zenith * scale_factor || solz < min_zenith * scale_factor)
                continue;
            if (sola > max_azimuth * scale_factor || sola < min_azimuth * scale_factor)
                continue;
            if (sena > max_azimuth * scale_factor || sena < min_azimuth * scale_factor)
                continue;
            // bin quality for l1b
            if (l1bfile.format.type == FT_OCIL1B) {
                for (size_t i_band = 0; i_band < n_bands; ++i_band)
                    l1cfile.quality_flags[index_l1c * n_bands + i_band] |=
                        l1bfile.quality_flags[i_pixel * n_bands + i_band];
            }

            if (l1bfile.area_weights.size() == 0) {  // no area weigting
                short& n_observations = l1cfile.number_of_observations[index_l1c];
                float& sensor_azimuth_angle = l1cfile.sensor_azimuth_angles_data[index_l1c];
                float& sensor_zenith_angle = l1cfile.sensor_zenith_angles_data[index_l1c];
                float& solar_azimuth_angle = l1cfile.solar_azimuth_angles_data[index_l1c];
                float& solar_zenith_angle = l1cfile.solar_zenith_angles_data[index_l1c];
                float& scattering_angle = l1cfile.scattering_angles_data[index_l1c];
                double& view_time_offsets = l1cfile.view_time_offsets_data[index_l1c];
                float& height_stdev = l1cfile.height_stdev[index_l1c];
                float cose = cos((senz + 180) * OEL_DEGRAD);
                float cosu = cos(solz * OEL_DEGRAD);
                float term1 = cose * cosu;
                float term2 = sqrt((1 - cose * cose)) * sqrt((1 - cosu * cosu));
                float term3 = cos((sena + 180 - sola) * OEL_DEGRAD);
                float sca = acos(term1 + term2 * term3) * OEL_RADEG;
                if (n_observations == 0) {
                    sensor_azimuth_angle = sena;
                    sensor_zenith_angle = senz;
                    solar_azimuth_angle = sola;
                    solar_zenith_angle = solz;
                    scattering_angle = sca;
                    view_time_offsets = scantime;
                    height_stdev = 0;
                    if (l1cfile.cloud_correction) {
                        float& aggregated_height = l1cfile.aggregated_height[index_l1c];
                        float& aggregated_height_stdev = l1cfile.aggregated_height_stdev[index_l1c];
                        float height = l1bfile.height_data_corrected[i_line * n_pixels + i_pixel] * 1000;
                        aggregated_height = height;
                        aggregated_height_stdev = 0;
                    }
                    n_observations = 1;
                } else {
                    sensor_azimuth_angle =
                        (sensor_azimuth_angle * n_observations + sena) / (n_observations + 1.0f);
                    solar_azimuth_angle =
                        (solar_azimuth_angle * n_observations + sola) / (n_observations + 1.0f);
                    sensor_zenith_angle =
                        (sensor_zenith_angle * n_observations + senz) / (n_observations + 1.0f);
                    solar_zenith_angle =
                        (solar_zenith_angle * n_observations + solz) / (n_observations + 1.0f);
                    scattering_angle = (scattering_angle * n_observations + sca) / (n_observations + 1.0f);
                    view_time_offsets =
                        (view_time_offsets * n_observations + scantime) / (n_observations + 1.0f);
                    height_stdev =
                        (height_stdev * n_observations + height_l1b * height_l1b) / (n_observations + 1.0f);
                    if (l1cfile.cloud_correction) {
                        float& mean = l1cfile.aggregated_height[index_l1c];
                        float& stdev = l1cfile.aggregated_height_stdev[index_l1c];
                        float value = l1bfile.height_data_corrected[i_line * n_pixels + i_pixel] * 1000;
                        float diff_prev = value - mean;
                        mean = mean + diff_prev / ((float)n_observations + 1.0f);
                        float diff_curr = value - mean;
                        stdev = (stdev * (float)n_observations + diff_curr * diff_prev) /
                                ((float)n_observations + 1.0f);
                    }
                    n_observations += 1;
                }

                // bin intensity
                for (size_t i_band = 0; i_band < n_bands; ++i_band) {
                    float Lt = l1rec.Lt[n_bands * i_pixel + i_band];
                    if (Lt == BAD_FLT)
                        continue;
                    float value = Lt * 10;
                    if (value > max_i || value < min_i)
                        continue;
                    short& count = l1cfile.number_of_observations_band[index_l1c * n_bands + i_band];
                    float& mean = l1cfile.i_data[index_l1c * n_bands + i_band];
                    float& stdev = l1cfile.i_stdev_data[index_l1c * n_bands + i_band];
                    if (count == 0) {
                        mean = value;
                        stdev = 0;
                        count = 1;
                    } else {
                        float diff_prev = value - mean;
                        mean = mean + diff_prev / ((float)count + 1.0f);
                        float diff_curr = value - mean;
                        stdev = (stdev * (float)count + diff_curr * diff_prev) / ((float)count + 1.0f);
                        count += 1;
                    }
                }
            } else {
                std::vector<std::pair<int, double>> area_weights =
                    l1bfile.area_weights[i_line * n_pixels + i_pixel];
                for (const auto& [index_l1c_area, area] : area_weights) {
                    short& n_observations = l1cfile.number_of_observations[index_l1c_area];
                    float& sensor_azimuth_angle = l1cfile.sensor_azimuth_angles_data[index_l1c_area];
                    float& sensor_zenith_angle = l1cfile.sensor_zenith_angles_data[index_l1c_area];
                    float& solar_azimuth_angle = l1cfile.solar_azimuth_angles_data[index_l1c_area];
                    float& solar_zenith_angle = l1cfile.solar_zenith_angles_data[index_l1c_area];
                    float& scattering_angle = l1cfile.scattering_angles_data[index_l1c_area];
                    double& view_time_offsets = l1cfile.view_time_offsets_data[index_l1c_area];
                    float & height_stdev = l1cfile.height_stdev[index_l1c_area];
                    float& area_total = l1cfile.area_total[index_l1c_area];

                    float cose = cos((senz + 180) * OEL_DEGRAD);
                    float cosu = cos(solz * OEL_DEGRAD);
                    float term1 = cose * cosu;
                    float term2 = sqrt((1 - cose * cose)) * sqrt((1 - cosu * cosu));
                    float term3 = cos((sena + 180 - sola) * OEL_DEGRAD);
                    float sca = acos(term1 + term2 * term3) * OEL_RADEG;
                    if (n_observations == 0) {
                        sensor_azimuth_angle = sena;
                        sensor_zenith_angle = senz;
                        solar_azimuth_angle = sola;
                        solar_zenith_angle = solz;
                        scattering_angle = sca;
                        view_time_offsets = scantime;
                        height_stdev = 0;
                        if (l1cfile.cloud_correction) {
                            float& aggregated_height = l1cfile.aggregated_height[index_l1c_area];
                            float& aggregated_height_stdev = l1cfile.aggregated_height_stdev[index_l1c_area];
                            float height = l1bfile.height_data_corrected[i_line * n_pixels + i_pixel] * 1000;
                            aggregated_height = height;
                            aggregated_height_stdev = 0;
                        }
                        n_observations = 1;
                        area_total = area;
                    } else {
                        sensor_azimuth_angle =
                            (sensor_azimuth_angle * area_total + sena * area) / (area_total + area);
                        solar_azimuth_angle =
                            (solar_azimuth_angle * area_total + sola * area) / (area_total + area);
                        sensor_zenith_angle =
                            (sensor_zenith_angle * area_total + senz * area) / (area_total + area);
                        solar_zenith_angle =
                            (solar_zenith_angle * area_total + solz * area) / (area_total + area);
                        scattering_angle = (scattering_angle * area_total + sca * area) / (area_total + area);
                        view_time_offsets =
                            (view_time_offsets * area_total + scantime * area) / (area_total + area);
                        height_stdev = (height_stdev * area_total + height_l1b * height_l1b * area) / (area_total + area);
                        if (l1cfile.cloud_correction) {
                            float& aggregated_height = l1cfile.aggregated_height[index_l1c_area];
                            float& aggregated_height_stdev = l1cfile.aggregated_height_stdev[index_l1c_area];
                            float height = l1bfile.height_data_corrected[i_line * n_pixels + i_pixel] * 1000;
                            aggregated_height_stdev =
                                ((aggregated_height_stdev + aggregated_height * aggregated_height) *
                                     area_total +
                                 height * height * area) /
                                (area_total + area);
                            aggregated_height =
                                (aggregated_height * area_total + height * area) / (area_total + area);
                            aggregated_height_stdev =
                                aggregated_height_stdev - aggregated_height * aggregated_height;
                        }
                        area_total += area;
                        n_observations += 1;
                        for (size_t i_band = 0; i_band < n_bands; ++i_band) {
                            float Lt = l1rec.Lt[n_bands * i_pixel + i_band];
                            if (Lt == BAD_FLT)
                                continue;
                            float value = Lt * 10;
                            if (value > max_i || value < min_i)
                                continue;
                            short& count =
                                l1cfile.number_of_observations_band[index_l1c_area * n_bands + i_band];
                            float& area_band = l1cfile.area_total_per_band[index_l1c_area * n_bands + i_band];
                            float& mean = l1cfile.i_data[index_l1c_area * n_bands + i_band];
                            float& stdev = l1cfile.i_stdev_data[index_l1c_area * n_bands + i_band];
                            float& area_second_momentum =
                                l1cfile.area_second_momentum[index_l1c_area * n_bands + i_band];
                            if (count == 0) {
                                mean = value;
                                stdev = 0;
                                area_second_momentum = value * value;
                                area_band = area;
                                count = 1;
                            } else {
                                mean = (mean * area_band + value * area) / (area_band + area);
                                area_second_momentum =
                                    (area_second_momentum * area_band + value * value * area) /
                                    (area_band + area);
                                stdev = area_second_momentum - mean * mean;
                                area_band += area;
                                count++;
                            }
                        }
                    }
                }
            }
        }
    }
    // binning here
    return 0;
}