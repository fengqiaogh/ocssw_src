#include "l1.h"

#include <math.h>

#define WAVE_INDEX_NUM 13000
#define WAVE_INDEX_MIN   300
#define WAVE_INDEX_MAX   (WAVE_INDEX_MIN + WAVE_INDEX_NUM)

static int band_index[WAVE_INDEX_NUM + 1];

/* ---------------------------------------------------------------------------------------- */
/* windex_set() - load direct access wavelength index                                       */

/* ---------------------------------------------------------------------------------------- */
void bindex_set(int32_t wave[], int nwave, int dwave_vswir) {
    int iw, iw1, iw2, i;
    int dwave;

    for (i = 0; i < WAVE_INDEX_NUM; i++) {
        band_index[i] = -1;
    }

    for (iw = 0; iw < nwave; iw++) {
        if (wave[iw] > 3000)
            dwave = 1000;
        else
            dwave = dwave_vswir;
        iw1 = MAX(wave[iw] - dwave / 2, WAVE_INDEX_MIN) - WAVE_INDEX_MIN;
        iw2 = MIN(wave[iw] + dwave / 2, WAVE_INDEX_MAX) - WAVE_INDEX_MIN;
        if (band_index[iw1] != -1) {
            /* adjust for overlapping band-passes */
            iw1 = MAX(wave[iw]-(wave[iw] - wave[iw - 1]) / 2, WAVE_INDEX_MIN) - WAVE_INDEX_MIN;
        }
        for (i = iw1; i <= iw2; i++) {
            band_index[i] = iw; // eliminated bindx
        }
    }
}


/* ---------------------------------------------------------------------------------------- */
/* bindex_get() - retrieve direct access wavelength index                                   */

/* ---------------------------------------------------------------------------------------- */
int bindex_get(int32_t wave) {
    if (wave >= WAVE_INDEX_MIN && wave < WAVE_INDEX_MAX)
        return (band_index[wave - WAVE_INDEX_MIN]);
    else
        return (-1);
}


/* ---------------------------------------------------------------------------------------- */
/* bindex_get_555() - retrieve direct access wavelength index of band nearest to 555nm      */

/* ---------------------------------------------------------------------------------------- */
int bindex_get_555(void) {
    int ib = bindex_get(547); /* note: we want 547 for MODIS, not 555, hence the order */

    if (ib < 0) ib = bindex_get(550);
    if (ib < 0) ib = bindex_get(555);
    if (ib < 0) ib = bindex_get(560);
    if (ib < 0) ib = bindex_get(565);

    return (ib);
}


/* ---------------------------------------------------------------------------------------- */
/* windex() - return wavelength index of table which is closest to sensor wavelength        */

/* ---------------------------------------------------------------------------------------- */
int windex(float wave, float twave[], int ntwave) {
    int iw, index;
    float wdiff;
    float wdiffmin = 99999.;

    for (iw = 0; iw < ntwave; iw++) {

        /* break on exact match */
        if (twave[iw] == wave) {
            index = iw;
            break;
        }

        /* look for closest */
        wdiff = fabs(twave[iw] - wave);
        if (wdiff < wdiffmin) {
            wdiffmin = wdiff;
            index = iw;
        }
    }

    return (index);
}

int windex_(float *wave, float twave[], int *ntwave) {
    return (windex(*wave, twave, *ntwave));
}


/* ---------------------------------------------------------------------------------------- */
/* invbindex() - invert the band index, given a wavelength index (oy vey!)                  */

/* ---------------------------------------------------------------------------------------- */
int invbindx(int band, int32_t *bindx, int nbands) {
    int i;

    for (i = 0; i < nbands; i++)
        if (band == bindx[i])
            return (i);

    return (-1);
}
