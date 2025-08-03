// Polarization Correction for OCI
// File is the same as polcor_hawkeye.cpp, but with some modifications
#include <math.h>
#include <allocate2d.h>
#include "l12_proto.h"
#include "polcor_oci.h"
#include "l1_oci_private.h"     // access polcor data set by l1 reader in l1rec->private_data

using namespace std;

/**
 * Calculates the polarization correction for the current pixel and sets it to index.
 * 
 * Index is the relative position of the current pixel to the polcor array, which 
 * has size numBands * numPixels for the current scan.
 * 
 * 
 * Calculate the value of m12 and m13 using the Second Degree Polynomial Function:
 *  Ax^2 + Bx + C, where:
 *      A = coef[2], 
 *      B = coef[1],
 *      C = coef[0],
 *      x = ccd scan angle
 * 
 */
void calculateAndSetPolCor(l1str* &l1rec, int currPix, int index, vector<float> m12Coefs, vector<float> m13Coefs, float scanAngle) {

    double alpha = l1rec->alpha[currPix] / OEL_RADEG;
    double L_x = l1rec->Lt[index] / l1rec->tg_sol[index] / l1rec->tg_sen[index];
    double L_qp = l1rec->L_q[index] * cos(2 * alpha) + l1rec->L_u[index] * sin(2 * alpha);
    double L_up = l1rec->L_u[index] * cos(2 * alpha) - l1rec->L_q[index] * sin(2 * alpha);

    /*
        If L_qp and L_up are NaN, skip. polcor == 1.0.
        --
        l1rec->L_u and L_q are NaN when airmass is negative because log(negative) is undefined.
        l1rec->L_u and L_q are calculated in rayleigh.c and it uses the function ray_press_wang to
        calculate a fraction. If airmass is negative, that fraction is NaN since it uses log of airmass. 
        
        This instance happens when solar_zenith (solz) is >= 90, resulting in a negative cos after
        it is converted into radian. This value carries over when calculating airmass.
    */
    if (isnan(L_qp) || isnan(L_up)) return;
    
    // calculate m12
    float A = m12Coefs[2];
    float B = m12Coefs[1];
    float C = m12Coefs[0];
    double m12 = (A * pow(scanAngle, 2.0)) + (B * scanAngle) + C;

    // calculate m13
    A = m13Coefs[2];
    B = m13Coefs[1];
    C = m13Coefs[0];
    double m13 = (A * pow(scanAngle, 2.0)) + (B * scanAngle) + C;

    // calculate polcor 
    l1rec->polcor[index] = 1.0 / (1.0 - m12 * L_qp / L_x - m13 * L_up / L_x);
    l1rec->dpol[index] = sqrt(pow(l1rec->L_q[index], 2.0) + pow(l1rec->L_u[index], 2.0)) / L_x;
}

/**
 * For the current scan, grab and merge the blue and red m12 and m13 bands.
 * Also grab the swir bands based on saturation (max 2).
 * Then, for every band, apply the polarization correction for currPix
 * @param l1rec - L1B file -- contains the current scan line this pixel belongs to
 * @param currPix - current pixel to apply polcor
 * 
 */
extern "C" void polcor_oci(l1str *l1rec, int32_t currPix) {

    // Skip if it is land pixel. atmocor_land.c sets polcor to 1 for land.
    int isLandPixel = (int)l1rec->land[currPix];
    if (isLandPixel) {
        return;
    }
    
    // grab the private data used for polarization from the private data pointer
    PolcorOciData* polcorPrivateData = (PolcorOciData*)l1rec->private_data;


    int numBands = l1rec->l1file->nbands;
    int numBlueBands = polcorPrivateData->num_blue_bands;
    int numRedBands = polcorPrivateData->num_red_bands;
    int numSwirBands = polcorPrivateData->num_swir_bands;


    // track the number of bands that were used in calculating polcor because
    // the array for l1rec->polcor, Lt, tg_sol, tg_sen, L_q and L_u are flatten and
    // are of size number of bands * number of pixels. Blue, Red and Swir bands during
    // this calculation step are separate, so we need to track the band count so we 
    // access the correct array. 
    int bandCount = 0; 
    //int currHamSide = hamSide[currScanLine];

    // grab the scan angle for the current pixel from l1rec
    float scanAngle = polcorPrivateData->ccdScanAngle[currPix];

    // Process Blue Bands first, skipping last 3 bands
    // NOTE: **ORDER MATTERS** -- when merging the band coefs, the order is blue, red, swir
    // Skip last 3 bands because it overlaps with red. Use red over blue. m12 and m13 
    // already skips the blue bands, but not the SWIR ones.
    for (int currBand = 0; currBand < numBlueBands-3; currBand++) {

        // relative position for this pixel and this band in the polcor array.
        // polcor = numBands * numpixels for the scan line.
        // to get the index of currPix relative to it, it will be 
        // currpix * numBands + curr band position.
        // not using currBand but bandCount bc some bands are skipped and because
        // this calculation is done for each band color
        int relativeIndex = currPix * numBands + bandCount;

        // grab the blue m12 and m13 coef from the private_data. l1 reader reads in the m12 and m13 blue coefs
        // for the current ham side, so no need to worry about which side it is. 
        // m12 and 13 is flatten, so it is of size: num_blue_bands * polarization coef (3)
        vector<float> m12BlueCoefs = {
            polcorPrivateData->blueM12Coefs[currBand * HAM_SIDES * POLARIZATION_COEFS + l1rec->mside * POLARIZATION_COEFS + 0], // current band's coef 0
            polcorPrivateData->blueM12Coefs[currBand * HAM_SIDES * POLARIZATION_COEFS + l1rec->mside * POLARIZATION_COEFS + 1], // coef 1
            polcorPrivateData->blueM12Coefs[currBand * HAM_SIDES * POLARIZATION_COEFS + l1rec->mside * POLARIZATION_COEFS + 2]  // coef 2
        };

        vector<float> m13BlueCoefs = {
            polcorPrivateData->blueM13Coefs[currBand * HAM_SIDES * POLARIZATION_COEFS + l1rec->mside * POLARIZATION_COEFS + 0], // current band's coef 0
            polcorPrivateData->blueM13Coefs[currBand * HAM_SIDES * POLARIZATION_COEFS + l1rec->mside * POLARIZATION_COEFS + 1], // coef 1
            polcorPrivateData->blueM13Coefs[currBand * HAM_SIDES * POLARIZATION_COEFS + l1rec->mside * POLARIZATION_COEFS + 2]  // coef 2
        };

        // vector<float> m12Coef = m12[currHamSide][bandCount];
        // vector<float> m13Coef = m13[currHamSide][bandCount];

        calculateAndSetPolCor(l1rec, currPix, relativeIndex, m12BlueCoefs, m13BlueCoefs, scanAngle);
        bandCount++;
    }


    // Process Red Bands
    for (int currBand = 0; currBand < numRedBands; currBand++) {

        int relativeIndex = currPix * numBands + bandCount;
        // m12 and m13 are merged, so use the current band count to get the red band
        vector<float> m12RedCoefs = {
            polcorPrivateData->redM12Coefs[currBand * HAM_SIDES * POLARIZATION_COEFS + l1rec->mside * POLARIZATION_COEFS + 0], // current band's coef 0
            polcorPrivateData->redM12Coefs[currBand * HAM_SIDES * POLARIZATION_COEFS + l1rec->mside * POLARIZATION_COEFS + 1], // coef 1
            polcorPrivateData->redM12Coefs[currBand * HAM_SIDES * POLARIZATION_COEFS + l1rec->mside * POLARIZATION_COEFS + 2]  // coef 2
        };

        vector<float> m13RedCoefs = {
            polcorPrivateData->redM13Coefs[currBand * HAM_SIDES * POLARIZATION_COEFS + l1rec->mside * POLARIZATION_COEFS + 0], // current band's coef 0
            polcorPrivateData->redM13Coefs[currBand * HAM_SIDES * POLARIZATION_COEFS + l1rec->mside * POLARIZATION_COEFS + 1], // coef 1
            polcorPrivateData->redM13Coefs[currBand * HAM_SIDES * POLARIZATION_COEFS + l1rec->mside * POLARIZATION_COEFS + 2]  // coef 2
        };
        // vector<float> m12Coef = m12[currHamSide][bandCount];
        // vector<float> m13Coef = m13[currHamSide][bandCount];
        calculateAndSetPolCor(l1rec, currPix, relativeIndex, m12RedCoefs, m13RedCoefs, scanAngle);
        bandCount++;
    }


    // check if the current pixels SWIR high gain bands are satured
    // if they are, we skip bands 3 and 6
    //bool skipHighGain = areHighGainBandsSaturated(currPix);
    bool skipHighGain = (int)l1rec->hilt[currPix]; // cast to int since it is a char or byte

    //this var will track the band for polcor, which will NOT increment when dropping
    // the high or low gain bands. bandCount will increment because it accounts for
    // all SWIR bands and it needs to update to grab the correct m12 and m13 coefs 
    //int polcorBandCount = bandCount;
    
    // Process SWIR Bands
    for (int currBand = 0; currBand < numSwirBands; currBand++) {

        // if skipping high gain, skip bands 3 and 6
        if (skipHighGain && (currBand == 3 || currBand == 6)) {
            continue;
        }

        // if not skipping high gain, throw out the low gain ones
        if (!skipHighGain && (currBand == 2 || currBand == 4)) {
            continue;
        }

        // band count because it is used to save into polcor arr
        int relativeIndex = currPix * numBands + bandCount;
        // vector<float> m12Coef = m12[currHamSide][bandCount];
        // vector<float> m13Coef = m13[currHamSide][bandCount];

        vector<float> m12SwirCoefs = {
            polcorPrivateData->swirM12Coefs[currBand * HAM_SIDES * POLARIZATION_COEFS + l1rec->mside * POLARIZATION_COEFS + 0], // current band's coef 0
            polcorPrivateData->swirM12Coefs[currBand * HAM_SIDES * POLARIZATION_COEFS + l1rec->mside * POLARIZATION_COEFS + 1], // coef 1
            polcorPrivateData->swirM12Coefs[currBand * HAM_SIDES * POLARIZATION_COEFS + l1rec->mside * POLARIZATION_COEFS + 2]  // coef 2
        };

        vector<float> m13SwirCoefs = {
            polcorPrivateData->swirM13Coefs[currBand * HAM_SIDES * POLARIZATION_COEFS + l1rec->mside * POLARIZATION_COEFS + 0], // current band's coef 0
            polcorPrivateData->swirM13Coefs[currBand * HAM_SIDES * POLARIZATION_COEFS + l1rec->mside * POLARIZATION_COEFS + 1], // coef 1
            polcorPrivateData->swirM13Coefs[currBand * HAM_SIDES * POLARIZATION_COEFS + l1rec->mside * POLARIZATION_COEFS + 2]  // coef 2
        };

        calculateAndSetPolCor(l1rec, currPix, relativeIndex, m12SwirCoefs, m13SwirCoefs, scanAngle);
        bandCount++;
        
    }
   
}
