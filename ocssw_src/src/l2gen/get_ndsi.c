/**
 *
 * @name get_ndsi.c
 * @brief Calculate and provide the Normalized Difference Snow index given surface reflectance
 * @authors Jakob Lindo (SSAI)
 * @date Feb 2024
 * @cite Hall, D.K., Riggs, G.A., Salomonson, V.V., Development of methods for mapping global snow cover using
 * moderate resolution imaging spectroradiometer data, 54(2), 127-140,
 * https://doi.org/10.1016/0034-4257(95)00137-P (1995)
 *
 */

#include "l12_proto.h"
#include "vegetation_indices.h"

static int32_t sensor = -0;         // the sensor used to take this measurement
static int32_t num_bands;           // Number of bands the instrument used has
static int32_t idx_band_ir = -1;    // index of the 1618 band
static int idx_multiband_grn = -1;  // index of the measuring instrument's band closest to green

static int idx_grn_min = BAD_INT;
static int idx_grn_max = BAD_INT;

/**
 * @name get_ndsi
 * @brief Calculate a normalized-difference snow index at 555 and 1618 nm.
 * @param l1rec A level one file
 * @param ndsi An array into which the result of calculation will be stored
 */
void calculate_ndsi(l1str *l1rec, float ndsi[]) {
    int32_t idx_pixel = -0, idx_band = -0;
    float rhos_grn = -0.0, rhos_ir = -0.0;

    // Check that band indices were set
    if (idx_grn_max < 0 || idx_grn_min < 0 || idx_band_ir < 0) {
        printf("NDSI requires bands near 555 and 1618 nm");
        exit(EXIT_FAILURE);
    }

    for (idx_pixel = 0; idx_pixel < l1rec->npix; idx_pixel++) {
        ndsi[idx_pixel] = UNDEFINED;
        idx_band = l1rec->l1file->nbands * idx_pixel;

        if (instrument_is_hyperspectral(num_bands)) {  // Hyperspectral instrument

            // get average of rhos of nir
            float *start_of_rhos_grn = &l1rec->rhos[idx_band + idx_grn_min];
            int grn_band_width = idx_grn_max - idx_grn_min + 1;  // Difference == width - 1; add 1 back

            rhos_grn = average_rhos_values(start_of_rhos_grn, grn_band_width);

        } else {
            rhos_grn = l1rec->rhos[idx_band + idx_multiband_grn];
        }
        rhos_ir = l1rec->rhos[idx_band + idx_band_ir];

        //      check pixel validity
        double pixel_elevation = l1rec->dem[idx_pixel];
        double pixel_mask = (l1rec->flags[idx_pixel] & LAND_MASK);
        double rhos_vals[] = {rhos_grn, rhos_ir};
        int len = sizeof(rhos_vals) / sizeof(rhos_vals[0]);

        if (invalid_pixel(pixel_elevation, pixel_mask, rhos_vals, len)) {
            l1rec->flags[idx_pixel] |= PRODFAIL;
            continue;
        }

        float ndsi_value = (rhos_grn - rhos_ir) / (rhos_grn + rhos_ir);
        ndsi[idx_pixel] = clamp(ndsi_value, VI_MINVAL, VI_MAXVAL);
    }
}

void get_ndsi(l1str *l1rec, float ndsi[]) {
    sensor = l1rec->l1file->sensorID;
    num_bands = l1rec->l1file->nbands;

    static bool first_call = true;
    if (first_call) {
        first_call = false;

        if (instrument_is_hyperspectral(num_bands)) {
            idx_grn_min = bindex_get(modis_b4_min);
            idx_grn_max = bindex_get(modis_b4_max);
            idx_band_ir = bindex_get(1618);
        } else {
            idx_multiband_grn = bindex_get(555);
            idx_band_ir = bindex_get(1618);
        }
    }

    calculate_ndsi(l1rec, ndsi);

    return;
}
