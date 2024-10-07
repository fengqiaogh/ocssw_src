/**
 * @name get_ndii.c - Normalized Difference Infrared Index
 * @brief Calculate and provide the Normalized Difference Infrared Index given surface reflectance
 * @cite Hardisky M.A., Klemas V., Smart R.M., The influence of soil-salinity, growth form, and leaf moisture
 * on the spectral radiance of Spartina alterniflora canopies. Photogrammetric Engineering and Remote Sensing
 * 49, 77–83 (1983).
 * @authors Jakob Lindo (SSAI)
 * @date Feb 2024
 */

#include "l12_proto.h"
#include "vegetation_indices.h"

static int32_t idx_band_ir = -1;  // index of the 1618 band
// index indicating location of the rhos values of NIR in a multiband instrument
static int idx_multiband_nir;

static int idx_nir_min = BAD_INT;  // A band index for NIR
static int idx_nir_max = BAD_INT;  // A band index for NIR

static int32_t num_bands;  // Number of bands the instrument used has

/**
 * @name get_ndii
 * @brief Calculate a normalized difference infrared index
 * @cite Hardisky M.A., Klemas V., Smart R.M., The influence of soil-salinity, growth form, and leaf moisture
 * on the spectral radiance of Spartina alterniflora canopies. Photogrammetric Engineering and Remote Sensing
 * 49, 77–83 (1983).
 * @param l1rec A level one file
 * @param ndii An array into which the result of calculation will be stored
 */
void calculate_ndii(l1str *l1rec, float ndii[]) {
    int32_t idx_pixel = -0, idx_band = -0;
    float rhos_nir = -0.0, rhos_ir = -0.0;

    // Check that band indices were set
    if (idx_nir_max < 0 || idx_nir_min < 0 || idx_band_ir < 0) {
        printf("NDII requires NIR bands between 670 and 875, and at 1618 nm");
        exit(EXIT_FAILURE);
    }

    for (idx_pixel = 0; idx_pixel < l1rec->npix; idx_pixel++) {
        ndii[idx_pixel] = UNDEFINED;
        idx_band = l1rec->l1file->nbands * idx_pixel;

        if (instrument_is_hyperspectral(num_bands)) {  // Hyperspectral instrument

            // get average of rhos of nir
            float *start_of_rhos_nir = &l1rec->rhos[idx_band + idx_nir_min];
            int nir_band_width = idx_nir_max - idx_nir_min + 1;

            rhos_nir = average_rhos_values(start_of_rhos_nir, nir_band_width);
            rhos_ir = l1rec->rhos[idx_band + idx_band_ir];

        } else {
            rhos_nir = l1rec->rhos[idx_band + idx_multiband_nir];
            rhos_ir = l1rec->rhos[idx_band + idx_band_ir];
        }

        //      check pixel validity
        double pixel_elevation = l1rec->dem[idx_pixel];
        double pixel_mask = (l1rec->flags[idx_pixel] & LAND_MASK);
        double rhos_vals[] = {rhos_nir, rhos_ir};
        int len = sizeof(rhos_vals) / sizeof(rhos_vals[0]);

        if (invalid_pixel(pixel_elevation, pixel_mask, rhos_vals, len)) {
            l1rec->flags[idx_pixel] |= PRODFAIL;
            continue;
        }

        float ndii_value = (rhos_nir - rhos_ir) / (rhos_nir + rhos_ir);
        ndii[idx_pixel] = clamp(ndii_value, VI_MINVAL, VI_MAXVAL);
    }
}

/**
 * @brief Main entry point for getting ndii
 * @param l1rec A level 1 record
 * @param prod A caller-provided array that will contain the result of calculation
 */
void get_ndii(l1str *l1rec, float prod[]) {
    num_bands = l1rec->l1file->nbands;

    static bool first_call = true;
    if (first_call) {
        first_call = false;

        if (instrument_is_hyperspectral(num_bands)) {  // TODO: Get rid of magic numbers
            idx_nir_min = bindex_get(nir_min);
            idx_nir_max = bindex_get(nir_max);
            idx_band_ir = bindex_get(1618);
        } else {
            idx_multiband_nir = bindex_get(859);
            idx_band_ir = bindex_get(1618);
        }

        // Should I check for MODIS(A|T)?
    }

    calculate_ndii(l1rec, prod);  // Produce ndii product
}