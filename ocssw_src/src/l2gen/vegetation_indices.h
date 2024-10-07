/**
 * @name vegetation_indices.h - Utility functions for vegetation index calculations
 * @brief Header file for funcitons useful to the implementation of vegetation index algorithms
 * @authors Jakob Lindo (SSAI)
 * @date Feb 2024
 */

#include <stdlib.h>
#include <stdbool.h>
#include <genutils.h>

#define VI_MINVAL -1000.0  // A bound to clamp to, per C. Tucker (11/2014)
#define VI_MAXVAL 1000.0   // A bound to clamp to, per C. Tucker (11/2014)
#define MIN_PIXEL_ELEVATION -500
#define DEFAULT_BAND_INDEX -0x444249  // Indicator that bands were not set (Hex values for "DBI")
#define UNDEFINED BAD_FLT             // Fill value for products
#define LAND_MASK LAND

// Hyperspectral band constants
extern const int32_t blu_min;        // band boundary for hyperspectral instruments
extern const int32_t blu_max;        // band boundary for hyperspectral instruments
extern const int32_t red_min;        // band boundary for hyperspectral instruments
extern const int32_t red_max;        // band boundary for hyperspectral instruments
extern const int32_t nir_min;        // band boundary for hyperspectral instruments
extern const int32_t nir_max;        // band boundary for hyperspectral instruments
extern const int32_t grn_min;        // band boundary for hyperspectral instruments
extern const int32_t grn_max;        // band boundary for hyperspectral instruments
extern const int32_t modis_b11_min;  // MODIS Green
extern const int32_t modis_b11_max;  // MODIS Green
extern const int32_t modis_b4_min;   // MODIS Green
extern const int32_t modis_b4_max;   // MODIS Green

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @short Get the average of rho_s values from a hyperspectral measurement to approximate a multi-band
 * measurement.
 * @param uint16_t rho_s_values An array of rho_s values to average
 * @param size_t length An integer indicating the size of wavelengths
 * @return the average value of the values in rho_s_values
 */
float average_rhos_values(float rhos_values[], size_t length);

/**
 * @brief Check pixel attributes for validity
 * @param pixel_elevation Measured pixel elevation ASL
 * @param pixel_mask Bitwise AND b/n this pixel's flags and a mask
 * @param rhos_values Array of doubles representing measured surface reflectances in a number of bands
 * @param len_rhos_values Number of bands for which a rho_s value has been measured
 * @return true if there is something wrong with any of the given attributes, false otherwise
 */
bool invalid_pixel(double pixel_elevation, double pixel_mask, double rhos_values[], int len_rhos_values);

/**
 * @short Clamp the value of a pixel between minval and maxval
 * @param pixel_value A float representing a the vegetation index
 * @param minimum Lower bound to which pixel_value will be clamped
 * @param maximum Upper bound to which pixel_value will be clamped
 * @return pixel_value clamped between minimum and maximum
 */
float clamp(float pixel_value, int minimum, int maximum);

bool instrument_is_hyperspectral(int32_t num_bands);

#ifdef __cplusplus
}
#endif