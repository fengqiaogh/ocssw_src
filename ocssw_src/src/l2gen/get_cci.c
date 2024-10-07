/**
 *
 * @name get_cci.c
 * @authors Jakob Lindo
 * @date Feb 2024
 * @details Takes in a level 1 structure and outputs the chlorophyll-carotenoid index for each pixel
 * @cite Gamon, J. A., Huemmrich, K. F., Wong, C. Y., Ensminger, I., Garrity, S., Hollinger, D. Y., Noormets,
 * A., and Peñuelas, J., A remotely sensed pigment index reveals photosynthetic phenology in evergreen
 * conifers. Proceedings of the National Academy of Sciences, 113(46), 13087-13092,
 * doi.org/10.1073/pnas.1606162113, (2016)
 *
 */

#include "l12_proto.h"
#include "vegetation_indices.h"

static int32_t num_bands;  // Number of bands the instrument used has
static uint idx_modis_band_11 = -1;
static uint idx_modis_band_1 = -1;
static int idx_red_min = BAD_INT;  // A band index for red
static int idx_red_max = BAD_INT;  // A band index for red
static int idx_grn_min = BAD_INT;  // A band index for green
static int idx_grn_max = BAD_INT;  // A band index for green

/**
 * @name get_cci
 * @brief Calculate a chlorophyll-carotenoid index at in red and green spectra.
 * @cite Gamon, J. A., Huemmrich, K. F., Wong, C. Y., Ensminger, I., Garrity, S., Hollinger, D. Y., Noormets,
 * A., and Peñuelas, J., A remotely sensed pigment index reveals photosynthetic phenology in evergreen
 * conifers. Proceedings of the National Academy of Sciences, 113(46), 13087-13092,
 * doi.org/10.1073/pnas.1606162113, (2016)
 * @param l1rec A level one file
 * @param cci An array into which the result of calculation will be stored
 */
void calculate_cci(l1str *l1rec, float cci[]) {
    double rhos_grn = -0.0;
    double rhos_red = 0.0;
    int32_t num_pixels = l1rec->npix;

    for (int32_t idx_pixel = 0; idx_pixel < num_pixels; idx_pixel++) {
        cci[idx_pixel] = UNDEFINED;
        int32_t idx_band = l1rec->l1file->nbands * idx_pixel;

        if (instrument_is_hyperspectral(num_bands)) {  // In future, check for hyperspectral ins., not just OCI
            uint grn_width = idx_grn_max - idx_grn_min + 1;  // Difference == width - 1; add 1 back
            uint red_width = idx_red_max - idx_red_min + 1;  // Difference == width - 1; add 1 back
            float *start_of_rhos_grn = &l1rec->rhos[idx_band + idx_grn_min];
            float *start_of_rhos_red = &l1rec->rhos[idx_band + idx_red_min];

            rhos_grn = average_rhos_values(start_of_rhos_grn, grn_width);
            rhos_red = average_rhos_values(start_of_rhos_red, red_width);
        } else {  // Multi band instrument
            rhos_grn = l1rec->rhos[idx_band + idx_modis_band_11];
            rhos_red = l1rec->rhos[idx_band + idx_modis_band_1];
        }

        double rhos_values[] = {rhos_grn, rhos_red};
        double pixel_elevation = l1rec->dem[idx_pixel];
        double landmask = l1rec->flags[idx_pixel] & LAND_MASK;
        size_t len_rhos_vals = sizeof(rhos_values) / sizeof(rhos_values[0]);
        if (invalid_pixel(pixel_elevation, landmask, rhos_values, len_rhos_vals)) {
            l1rec->flags[idx_pixel] |= PRODFAIL;
            continue;
        }

        double cci_value = (rhos_grn - rhos_red) / (rhos_grn + rhos_red);
        cci[idx_pixel] = clamp(cci_value, VI_MINVAL, VI_MAXVAL);
    }

    return;
}

void get_cci(l1str *l1rec, float cci[]) {
    num_bands = l1rec->l1file->nbands;

    static bool first_call = true;
    if (first_call) {
        first_call = false;

        if (instrument_is_hyperspectral(num_bands)) {
            idx_grn_min = bindex_get(modis_b11_min);
            idx_grn_max = bindex_get(modis_b11_max);
            idx_red_min = bindex_get(red_min);
            idx_red_max = bindex_get(red_max);

            if (idx_grn_min == DEFAULT_BAND_INDEX || idx_grn_max == DEFAULT_BAND_INDEX ||
                idx_red_min == DEFAULT_BAND_INDEX || idx_red_max == DEFAULT_BAND_INDEX) {
                printf("CCI requires bands near 555 and 645 nm\n");
                exit(EXIT_FAILURE);
            }
        } else {
            idx_modis_band_11 = bindex_get(555);
            idx_modis_band_1 = bindex_get(645);

            if (idx_modis_band_11 < 0 || idx_modis_band_1 < 0) {
                printf("CCI requires bands near 555 and 645 nm\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    calculate_cci(l1rec, cci);

    return;
}
