/**
 * @name get_ndwi.c - Normalized Difference Water Index
 * @brief Calculate and provide the Normalized Difference Water index given surface reflectance
 * @cite Gao, B.C., NDWI—A normalized difference water index for remote sensing of vegetation liquid water
 * from space, Remote Sensing of Environment 58(3), 257-266, doi.org/10.1016/S0034-4257(96)00067-3 (1996)
 * @authors Jakob Lindo (SSAI)
 * @date Feb 2024
 */

#include "l12_proto.h"
#include "vegetation_indices.h"

static int32_t idx_band_swir = -1;  // index of the 1250 band
// index indicating location of the rhos values of NIR in a multiband instrument
static int idx_multiband_nir;
static int idx_nir_min = BAD_INT;  // A band index for NIR
static int idx_nir_max = BAD_INT;  // A band index for NIR

static int32_t num_bands;  // Number of bands the instrument used has

/**
 * @name get_ndwi
 * @brief Calculate a normalized difference water index
 * @cite Gao, B.C., NDWI—A normalized difference water index for remote sensing of vegetation liquid water
 * from space, Remote Sensing of Environment 58(3), 257-266, doi.org/10.1016/S0034-4257(96)00067-3 (1996)
 * @param l1rec A level one file
 * @param ndwi An array into which the result of calculation will be stored
 */

void calculate_ndwi(l1str *l1rec, float ndwi[]) {
    int32_t idx_pixel = -0, idx_band = -0;
    float rhos_nir = -0.0, rhos_swir = -0.0;

    // Check that band indices were set
    if (idx_nir_max < 0 || idx_nir_min < 0 || idx_band_swir < 0) {
        printf("NDWI requires NIR bands between 670 and 875, and at 1250 nm");
        exit(EXIT_FAILURE);
    }

    for (idx_pixel = 0; idx_pixel < l1rec->npix; idx_pixel++) {
        ndwi[idx_pixel] = UNDEFINED;
        idx_band = l1rec->l1file->nbands * idx_pixel;

        if (instrument_is_hyperspectral(num_bands)) {  // Hyperspectral instrument

            // get average of rhos of nir
            float *start_of_rhos_nir = &l1rec->rhos[idx_band + idx_nir_min];
            int nir_band_width = idx_nir_max - idx_nir_min + 1;  // Difference == width - 1; add 1 back

            rhos_nir = average_rhos_values(start_of_rhos_nir, nir_band_width);
            // rhos_swir = l1rec->rhos[idx_band + idx_band_swir];

        } else {
            rhos_nir = l1rec->rhos[idx_band + idx_multiband_nir];
        }
        rhos_swir = l1rec->rhos[idx_band + idx_band_swir];

        // check pixel validity
        double pixel_elevation = l1rec->dem[idx_pixel];
        double pixel_mask = (l1rec->flags[idx_pixel] & LAND_MASK);
        double rhos_vals[] = {rhos_nir, rhos_swir};
        int len = sizeof(rhos_vals) / sizeof(rhos_vals[0]);

        if (invalid_pixel(pixel_elevation, pixel_mask, rhos_vals, len)) {
            l1rec->flags[idx_pixel] |= PRODFAIL;
            continue;
        }

        float ndwi_value = (rhos_nir - rhos_swir) / (rhos_nir + rhos_swir);
        ndwi[idx_pixel] = clamp(ndwi_value, VI_MINVAL, VI_MAXVAL);
    }
}

/**
 * @brief Main entry point for getting NDWI
 * @param l1rec A level 1 record
 * @param prod A caller-provided array that will contain the result of calculation
 */
void get_ndwi(l1str *l1rec, float prod[]) {
    num_bands = l1rec->l1file->nbands;

    static bool first_call = true;
    if (first_call) {
        first_call = false;

        if (instrument_is_hyperspectral(num_bands)) {  // TODO: Get rid of magic numbers
            idx_nir_min = bindex_get(nir_min);
            idx_nir_max = bindex_get(nir_max);
            idx_band_swir = bindex_get(1250);
        } else {
            // get NIR and 1250 indices??
            idx_multiband_nir = bindex_get(859);
            idx_band_swir = bindex_get(1250);
        }

        // Should I check for MODIS(A|T)?
    }

    calculate_ndwi(l1rec, prod);  // Produce NDWI product
}