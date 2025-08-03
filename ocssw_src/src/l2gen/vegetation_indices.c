/**
 * @name vegetation_indices.c - Utility functions for vegetation index calculations
 * @brief Functions useful to the implementation of vegetation index algorithms
 * @authors Jakob Lindo (SSAI)
 * @date Feb 2024
 */

#include "vegetation_indices.h"
#include <math.h>

#define NUM_BANDS_TO_BE_HYPERSPECTRAL 50

// Band constants for use with hyperspectral instuments
const int32_t blu_min = 459;        // blue band boundary for hyperspectral instruments
const int32_t blu_max = 479;        // blue band boundary for hyperspectral instruments
const int32_t red_min = 620;        // red band boundary for hyperspectral instruments
const int32_t red_max = 670;        // red band boundary for hyperspectral instruments
const int32_t nir_min = 841;        // near infrared band boundary for hyperspectral instruments
const int32_t nir_max = 876;        // near infrared band boundary for hyperspectral instruments
const int32_t grn_min = 480;        // green band boundary for hyperspectral instruments
const int32_t grn_max = 537;        // green band boundary for hyperspectral instruments
const int32_t modis_b4_min = 545;   // MODIS Green
const int32_t modis_b4_max = 565;   // MODIS Green
const int32_t modis_b11_min = 526;  // MODIS Green
const int32_t modis_b11_max = 536;  // MODIS Green

float average_rhos_values(const float rhos_values[], size_t length) {
    double sum = 0;
    for (int i = 0; i < length; i++) {
        sum += rhos_values[i];
    }
    return (sum / length);
}

bool invalid_pixel(double pixel_elevation, double pixel_mask, const double rhos_values[],
                   int len_rhos_values) {
    if (pixel_elevation <= MIN_PIXEL_ELEVATION || pixel_mask == 0)
        return true;

    for (int idx_rhos = 0; idx_rhos < len_rhos_values; idx_rhos++) {
        double rho_s = rhos_values[idx_rhos];
        if (rho_s < 0.0 || 1.0 < rho_s)
            return true;
    }

    // Nothing wrong found
    return false;
}

float clamp(float pixel_value, int minimum, int maximum) {
    return fmin(fmax(pixel_value, minimum), maximum);
}

bool instrument_is_hyperspectral(int32_t num_bands) {
    return num_bands > NUM_BANDS_TO_BE_HYPERSPECTRAL; 
}
