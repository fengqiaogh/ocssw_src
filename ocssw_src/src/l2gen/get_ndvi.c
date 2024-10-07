/*---------------------------------------------------------------------*/
/* get_ndvi.c - vegetation index classification for MSl12.             */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/* Outputs:                                                            */
/*     ndvi  - vegetation index for land, 1 value per pixel.           */
/*                                                                     */
/* Written by: Bryan Franz, SAIC-GSC, February 2000                    */
/*             Jakob Lindo, SSAI, February 2024                        */
/*                                                                     */
/*---------------------------------------------------------------------*/

#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"
#include "vegetation_indices.h"

static int idx_red_min = BAD_INT;
static int idx_red_max = BAD_INT;
static int idx_nir_min = BAD_INT;
static int idx_nir_max = BAD_INT;
static int idx_blu_min = BAD_INT;
static int idx_blu_max = BAD_INT;

// Multi-band constants
static float undef = BAD_FLT;
static int ibblue = -1;
static int ibred = -1;
static int ibnir = -1;

// as per C. Tucker (11/2014), no cut-off at -2
static float minval = -1000.0;
static float maxval = 1000.0;
static int32_t mask = LAND;

static int32_t sensor = -1;  // An enum indicating which sensor was used
static int32_t num_bands;    // Number of bands the instrument used has

/**
 * @short Calculate Normalized Difference Vegetation Index CAT_ix 32
 * @param l1rec Level 1 record
 * @param ndvi Array to place resultant values into
 */
void get_ndvi(l1str *l1rec, float ndvi[]) {
    int32_t idx_pixel = -0;  // The index of the current pixel
    int32_t idx_band = -0;   // Offset from l1rec->l1file->nbands; start of this pixel's bands
    float red = -0.0;
    float nir = -0.0;

    if ((!instrument_is_hyperspectral(num_bands) && (ibred < 0 || ibnir < 0)) ||
        (instrument_is_hyperspectral(num_bands) &&
         (idx_red_min < 0 || idx_nir_min < 0 || idx_red_max < 0 || idx_nir_max < 0))) {
        printf("NDVI requires bands near 670 and 865nm\n");
        exit(1);
    }

    for (idx_pixel = 0; idx_pixel < l1rec->npix; idx_pixel++) {
        ndvi[idx_pixel] = undef;
        idx_band = l1rec->l1file->nbands * idx_pixel;

        if (instrument_is_hyperspectral(num_bands)) {
            int red_width = idx_red_max - idx_red_min + 1;  // Difference == width - 1; add 1 back
            int nir_width = idx_nir_max - idx_nir_min + 1;  // Difference == width - 1; add 1 back
            float *start_of_rhos_red = &l1rec->rhos[idx_band + idx_red_min];
            float *start_of_rhos_nir = &l1rec->rhos[idx_band + idx_nir_min];

            red = average_rhos_values(start_of_rhos_red, red_width);  // Average rho_s b/n 620 and 670
            nir = average_rhos_values(start_of_rhos_nir, nir_width);  // Average rho_s b/n 841 and 876
        } else {
            red = l1rec->rhos[idx_band + ibred];
            nir = l1rec->rhos[idx_band + ibnir];
        }

        double rhos_values[] = {red, nir};
        if (invalid_pixel(l1rec->dem[idx_pixel], l1rec->flags[idx_pixel] & mask, rhos_values, 2)) {
            l1rec->flags[idx_pixel] |= PRODFAIL;
            continue;
        }
        float ndvi_value = (nir - red) / (nir + red);
        ndvi[idx_pixel] = clamp(ndvi_value, (int)minval, (int)maxval);
    }
}

/**
 * @short Calculate Enhanced Vegetation Index CAT_ix 38
 * @param l1rec Level 1 record
 * @param evi Array to place resultant values into
 */
void get_evi(l1str *l1rec, float evi[]) {
    static float L = 1.0, c1 = 6.0, c2 = 7.5;

    int32_t idx_pixel = -0;  // The index of the current pixel
    int32_t idx_band = -0;   // Offset from l1rec->l1file->nbands; start of this pixel's bands
    float blu = -0.0;
    float red = -0.0;
    float nir = -0.0;
    double val = -0.0;

    if ((!instrument_is_hyperspectral(num_bands) && (ibblue < 0 || ibred < 0 || ibnir < 0)) ||
        (instrument_is_hyperspectral(num_bands) && (idx_red_min < 0 || idx_red_max < 0 || idx_nir_min < 0 ||
                                                    idx_nir_max < 0 || idx_blu_max < 0 || idx_blu_max < 0))) {
        printf("EVI requires bands near 412, 670, and 865nm\n");
        exit(1);
    }

    for (idx_pixel = 0; idx_pixel < l1rec->npix; idx_pixel++) {
        evi[idx_pixel] = undef;
        idx_band = l1rec->l1file->nbands * idx_pixel;

        if (instrument_is_hyperspectral(num_bands)) {
            int blu_width = idx_blu_max - idx_blu_min + 1;  // Difference == width - 1; add 1 back
            int red_width = idx_red_max - idx_red_min + 1;  // Difference == width - 1; add 1 back
            int nir_width = idx_nir_max - idx_nir_min + 1;  // Difference == width - 1; add 1 back
            float *start_of_rhos_blu = &l1rec->rhos[idx_band + idx_blu_max];
            float *start_of_rhos_red = &l1rec->rhos[idx_band + idx_red_min];
            float *start_of_rhos_nir = &l1rec->rhos[idx_band + idx_nir_min];

            blu = average_rhos_values(start_of_rhos_blu, blu_width);  // Average rho_s b/n 455 and 480
            red = average_rhos_values(start_of_rhos_red, red_width);  // Average rho_s b/n 565 and 670
            nir = average_rhos_values(start_of_rhos_nir, nir_width);  // Average rho_s b/n 840 and 875
        } else {
            blu = l1rec->rhos[idx_band + ibblue];
            red = l1rec->rhos[idx_band + ibred];
            nir = l1rec->rhos[idx_band + ibnir];
        }

        double rhos_values[] = {red, nir, blu};
        if (invalid_pixel(l1rec->dem[idx_pixel], l1rec->flags[idx_pixel] & mask, rhos_values, 3)) {
            l1rec->flags[idx_pixel] |= PRODFAIL;
            continue;
        }
        double evi_val = 0xdeadbeef;
        if (blu > 0.0 && (blu <= red || red <= nir)) {
            /* Most cases - EVI formula */

            if ((val = L + nir + c1 * red - c2 * blu) == 0)
                continue;
            else
                evi_val = 2.5 * (nir - red) / val;

        } else {
            /* Backup - SAVI formula */

            if ((val = 0.5 + nir + red) == 0)
                continue;
            else
                evi_val = 1.5 * (nir - red) / val;
        }
        evi[idx_pixel] = clamp(evi_val, minval, maxval);
    }
}

// From Compton Tucker 07/30/2016
// EVI3=2.5*(nir-red)/(nir+6*red-7.5*blue+1)
// EVI2=2.5*(NIR-Red)/(NIR+2.4*Red+1)

/**
 * @short Calculate Enhanced Vegetation Index - EVI2 CAT_ix 63
 * @param l1rec Level 1 record
 * @param evi2 Array to place resultant values into
 */
void get_evi2(l1str *l1rec, float evi2[]) {
    int32_t idx_pixel = -0;  // The index of the current pixel
    int32_t idx_band = -0;   // Offset from l1rec->l1file->nbands; start of this pixel's bands
    float red;
    float nir;

    if ((!instrument_is_hyperspectral(num_bands) && (ibred < 0 || ibnir < 0)) ||
        (instrument_is_hyperspectral(num_bands) &&
         (idx_red_min < 0 || idx_nir_min < 0 || idx_red_max < 0 || idx_nir_max < 0))) {
        printf("EVI2 requires bands near 670, and 865nm\n");
        exit(1);
    }

    for (idx_pixel = 0; idx_pixel < l1rec->npix; idx_pixel++) {
        evi2[idx_pixel] = undef;
        idx_band = l1rec->l1file->nbands * idx_pixel;

        if (instrument_is_hyperspectral(num_bands)) {
            int red_width = idx_red_max - idx_red_min + 1;  // Difference == width - 1; add 1 back
            int nir_width = idx_nir_max - idx_nir_min + 1;  // Difference == width - 1; add 1 back
            float *start_of_rhos_red = &l1rec->rhos[idx_band + idx_red_min];
            float *start_of_rhos_nir = &l1rec->rhos[idx_band + idx_nir_min];

            red = average_rhos_values(start_of_rhos_red, red_width);  // Average rho_s b/n 565 and 670
            nir = average_rhos_values(start_of_rhos_nir, nir_width);  // Average rho_s b/n 840 and 875
        } else {
            red = l1rec->rhos[idx_band + ibred];
            nir = l1rec->rhos[idx_band + ibnir];
        }

        double rhos_values[] = {red, nir};
        if (invalid_pixel(l1rec->dem[idx_pixel], l1rec->flags[idx_pixel] & mask, rhos_values, 2)) {
            l1rec->flags[idx_pixel] |= PRODFAIL;
            continue;
        }
        double evi2_val = -0.0;
        if (red <= nir) {  // Wouldn't this always be the case?
            evi2_val = 2.5 * (nir - red) / (nir + 2.4 * red + 1);
        }
        evi2[idx_pixel] = clamp(evi2_val, minval, maxval);
    }
}

/**
 * @short Calculate Enhanced Vegetation Index - EVI3 CAT_ix 64
 * @param l1rec Level 1 record
 * @param evi3 Array to place resultant values into
 */
void get_evi3(l1str *l1rec, float evi3[]) {
    static float L = 1.0, c1 = 6.0, c2 = 7.5;

    int32_t idx_pixel = -0;  // The index of the current pixel
    int32_t idx_band = -0;   // Offset from l1rec->l1file->nbands; start of this pixel's bands
    float blu;
    float red;
    float nir;
    double val;

    if ((!instrument_is_hyperspectral(num_bands) && (ibblue < 0 || ibred < 0 || ibnir < 0)) ||
        (instrument_is_hyperspectral(num_bands) && (idx_red_min < 0 || idx_red_max < 0 || idx_nir_min < 0 ||
                                                    idx_nir_max < 0 || idx_blu_max < 0 || idx_blu_max < 0))) {
        printf("EVI requires bands near 412, 670, and 865nm\n");
        exit(1);
    }

    for (idx_pixel = 0; idx_pixel < l1rec->npix; idx_pixel++) {
        evi3[idx_pixel] = undef;
        idx_band = l1rec->l1file->nbands * idx_pixel;

        if (instrument_is_hyperspectral(num_bands)) {
            int blu_width = idx_blu_max - idx_blu_min + 1;  // Difference == width - 1; add 1 back
            int red_width = idx_red_max - idx_red_min + 1;  // Difference == width - 1; add 1 back
            int nir_width = idx_nir_max - idx_nir_min + 1;  // Difference == width - 1; add 1 back
            float *start_of_rhos_blu = &l1rec->rhos[idx_band + idx_blu_max];
            float *start_of_rhos_red = &l1rec->rhos[idx_band + idx_red_min];
            float *start_of_rhos_nir = &l1rec->rhos[idx_band + idx_nir_min];

            blu = average_rhos_values(start_of_rhos_blu, blu_width);  // Average rho_s b/n 455 and 480
            red = average_rhos_values(start_of_rhos_red, red_width);  // Average rho_s b/n 565 and 670
            nir = average_rhos_values(start_of_rhos_nir, nir_width);  // Average rho_s b/n 840 and 875
        } else {
            blu = l1rec->rhos[idx_band + ibblue];
            red = l1rec->rhos[idx_band + ibred];
            nir = l1rec->rhos[idx_band + ibnir];
        }

        double rhos_values[] = {red, nir, blu};
        if (invalid_pixel(l1rec->dem[idx_pixel], l1rec->flags[idx_pixel] & mask, rhos_values, 3)) {
            l1rec->flags[idx_pixel] |= PRODFAIL;
            continue;
        }

        double evi3_val = -0.0;
        if (blu > 0.0 && (blu <= red || red <= nir)) {
            /* Most cases - EVI formula */

            if ((val = L + nir + c1 * red - c2 * blu) == 0)
                continue;
            else
                evi3_val = 2.5 * (nir - red) / val;
        }
        evi3[idx_pixel] = clamp(evi3_val, minval, maxval);
    }
}

void get_ndvi_evi(l1str *l1rec, int prodnum, float prod[]) {
    sensor = l1rec->l1file->sensorID;
    num_bands = l1rec->l1file->nbands;

    static bool firstCall = true;
    if (firstCall) {
        firstCall = false;
        if (instrument_is_hyperspectral(num_bands)) {
            idx_blu_min = bindex_get(blu_min);
            idx_blu_max = bindex_get(blu_max);
            idx_red_min = bindex_get(red_min);
            idx_red_max = bindex_get(red_max);
            idx_nir_min = bindex_get(nir_min);
            idx_nir_max = bindex_get(nir_max);
        } else {
            ibblue = bindex_get(412);
            ibred = bindex_get(670);
            ibnir = bindex_get(865);
        }

        if (sensor == MODISA || sensor == MODIST) {
            ibred = bindex_get(645);
            ibnir = bindex_get(859);
        }
    }

    switch (prodnum) {
        case CAT_ndvi:
            get_ndvi(l1rec, prod);
            break;
        case CAT_evi:
            get_evi(l1rec, prod);
            break;
        case CAT_evi2:
            get_evi2(l1rec, prod);
            break;
        case CAT_evi3:
            get_evi3(l1rec, prod);
            break;
        default:
            printf("Error: %s : Unknown product specifier: %d\n", __FILE__, prodnum);
            exit(FATAL_ERROR);
    }
}