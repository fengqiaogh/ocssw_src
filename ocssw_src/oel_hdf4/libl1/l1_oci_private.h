#ifndef L1_OCI_PRIVATE_H
#define L1_OCI_PRIVATE_H

/*
    Defines private data for OCI that other records may not need for each line.
    Attach this structure to the private_data pointer of l1str 
*/

// m12 and m13 always have 3 coefs
const int POLARIZATION_COEFS = 3;
const int HAM_SIDES = 2;

typedef struct PolcorOciData {
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