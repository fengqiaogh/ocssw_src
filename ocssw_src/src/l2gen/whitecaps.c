#include "l12_proto.h"

/* -----------------------------------------------------------------
 * whitecaps() - whitecap radiances per band in mW/cm^2/um/sr
 *
 * Inputs:
 * mu0[]     - Cosine of solar zenith angle
 * mu[]      - Cosine of sensor zenith angle
 * Fo[]      - Solar irradiance (mW/cm^2/um/Sr)
 * tau_r[]   - Rayleigh optical thickness per band
 * press     - Atmospheric pressure (millibars)
 * ws        - Wind speed (m/s)
 * ws_max    - Maximum allowed wind speed (m/s)
 *
 * Ouputs:
 * tLf[]     - Whitecap radiance at sensor (mW/cm^2/um/Sr)
 *
 * History:
 *
 * B. A. Franz, SAIC, Ocean Discipline Processing Group, Feb. 1999
 * 
 * - Reduced whitecap reflectance by 75%, as suggested by Gordon
 *   and implemented by the SeaWiFS project. BAF, 26 June 1998.
 *
 * - Added Frouin wavelength dependence. BAF, 13 May 1999.
 *
 * - Raised albedo multiplier from 0.25 to 0.4 and updated to new
 *   Frouin wavelength dependence, BAF, 24 May 1999
 *
 * - Converted to C.  Externalized albedos, BAF, 5 Aug 2004.
 *
 * - remove ozone transmittance, BAF, Aug 2006.
 * 
 * - replace with Stramska model for eval, SWB, Nov 2008.
 *  ----------------------------------------------------------------- */

void whitecap_spectral_shape(int32_t nwave, float *wave, float *awc) {
    int nwc_tab = 8;
    float awc_wav[] = {412, 443, 490, 510, 555, 670, 765, 865};
    float awc_tab[] = {1.0, 1.0, 1.0, 1.0, 1.0, 0.889225, 0.760046, 0.644950};
    int iw;

    for (iw = 0; iw < nwave; iw++) {
        if (wave[iw] < 1000)
            awc[iw] = linterp(awc_wav, awc_tab, nwc_tab, wave[iw]);
        else
            awc[iw] = 0.0;
    }
}

void whitecaps(int32_t sensorID, int32_t evalmask, int32_t nwave,
        float ws, float ws_max, float rhof[]) {
    static int firstCall = 1;
    static float *awhite;

    float rhowc; /* White cap reflectance   */
    int32_t iw;

    // get spectral shape of whitecap reflectance, interpolated to
    // sensor spectral bands

    if (firstCall) {
        float *wave;
        firstCall = 0;
        rdsensorinfo(sensorID, evalmask, "fwave", (void **) &wave);
        if ((awhite = (float *) calloc(nwave, sizeof (float))) == NULL) {
            printf("Unable to allocate space for awhite.\n");
            exit(1);
        }
        whitecap_spectral_shape(nwave, wave, awhite);
    }

    // Stramska and Petelski (2003) wc fractional coverage for underdeveloped seas
    float ws_min = 6.33;
    if (ws > ws_min && ws_max > ws_min) {
        rhowc = 1.925E-5 * pow((MIN(ws, ws_max) - ws_min), 3.00);
        for (iw = 0; iw < nwave; iw++)
            rhof[iw] = awhite[iw] * rhowc;
    } else {
        for (iw = 0; iw < nwave; iw++)
            rhof[iw] = 0.0;
    }

}


