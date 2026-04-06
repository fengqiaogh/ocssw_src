#ifndef L1_OCI_PRIVATE_H
#define L1_OCI_PRIVATE_H

/*
    Defines private data for OCI that other records may not need for each line.
    Attach this structure to the private_data pointer of l1str 
*/

// m12 and m13 always have 3 coefs
static const int POLARIZATION_COEFS = 3;
static const int HAM_SIDES = 2;

static const size_t expected_num_blue_bands = 119;
static const size_t expected_num_red_bands = 163;
static const size_t expected_num_SWIR_bands = 9;

static const size_t skipped_blue_bands = 3; // Overlap b/n blue and red, prefer red
static const size_t skipped_SWIR_bands = 2; // Two bands (1250 & 1615) have low-gain and high-gain channels

/**
 * 
 * @note
 * Some blue bands overlap (3 of them at the end of the blue bands) with red. The blue bands in that overlap
 * are less than useful, so we skip those.
 * 
 * SWIR bands 2 and 5 (1250 and 1615) have low-gain and high-gain channels. 1250 low is 2 while the high is 3.
 * 1615 low is 5 while the high is 6.
 *
 * (119 - 3) + 163 + (9 - 2)
 */
static const size_t tot_num_bands = (expected_num_blue_bands - skipped_blue_bands) + expected_num_red_bands +
                                    (expected_num_SWIR_bands - skipped_SWIR_bands);

typedef struct PolcorOciData {
    int ncid_L1B;
    float *blueM12Coefs;
    float *blueM13Coefs;
    float *redM12Coefs;
    float *redM13Coefs;
    float *swirM12Coefs;
    float *swirM13Coefs;
    float *ccdScanAngle;    // size == numPix
    int num_blue_bands;
    int num_red_bands;
    int num_swir_bands;
} PolcorOciData;

#endif