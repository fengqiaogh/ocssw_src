/**
 *  @file smile.c
 *  @brief MERIS and OLCI smile correction
 *  @author Don Shea
 */

#include "smile.h"

#include <genutils.h>

#include <math.h>
#include <stdlib.h>


// max line length in a bandinfo.txt file
#define BANDINFO_LINE_MAX 255

// this is also used as a flag to see if smile correction has been initialized
static int numBands = 0;

static int numDetectors;

static int *switch_land; // [bands]
static int *lower_land; // [bands]
static int *upper_land; // [bands]
static int *switch_water; // [bands]
static int *lower_water; // [bands]
static int *upper_water; // [bands]
static float *theoWL; // [bands]
static float *theoE0; // [bands]

static float *detectorWL; // [detector][band]
static float *detectorE0; // [detector][band]

// set the detector center wavelengths
//
// note: that the dimensions are detectorWLs[numDetectors][numBands]

/**
 * Setup the smile correction information.
 *
 * @param num_bands number of bands for this sensor
 * @param num_detectors number of detectors (used in detectorWL, detectorE0)
 * @param bandinfo_filename full path to the band_info.txt file
 * @param detectorWLs detector central wavelengths (detectorWLs[numDetectors][numBands])
 * @param detectorE0s detector solar flux (detectorE0s[numDetectors][numBands])
 */
void smile_init(int num_bands, int num_detectors, const char* bandinfo_filename,
        float* detectorWLs, float* detectorE0s) {
    char line[BANDINFO_LINE_MAX];
    FILE *fp = NULL;
    int count;
    int dummy;
    int band;

    if (numBands != 0) {
        printf("-E- %s:%d - Can not re-initialize the smile correction.", __FILE__, __LINE__);
        exit(1);
    }
    if (bandinfo_filename == NULL) {
        printf("-E- %s:%d - bandinfo_filename can not be NULL.", __FILE__, __LINE__);
        exit(1);
    }
    if (detectorWLs == NULL) {
        printf("-E- %s:%d - detectorWLs can not be NULL.", __FILE__, __LINE__);
        exit(1);
    }
    if (detectorE0s == NULL) {
        printf("-E- %s:%d - detectorE0s can not be NULL.", __FILE__, __LINE__);
        exit(1);
    }

    numBands = num_bands;
    numDetectors = num_detectors;
    detectorWL = detectorWLs;
    detectorE0 = detectorE0s;

    switch_land = (int*) allocateMemory(sizeof (int)*numBands, "smile switch_land");
    lower_land = (int*) allocateMemory(sizeof (int)*numBands, "smile lower_land");
    upper_land = (int*) allocateMemory(sizeof (int)*numBands, "smile upper_land");
    switch_water = (int*) allocateMemory(sizeof (int)*numBands, "smile switch_water");
    lower_water = (int*) allocateMemory(sizeof (int)*numBands, "smile lower_water");
    upper_water = (int*) allocateMemory(sizeof (int)*numBands, "smile upper_water");
    theoWL = (float*) allocateMemory(sizeof (float)*numBands, "smile theoWL");
    theoE0 = (float*) allocateMemory(sizeof (float)*numBands, "smile theoE0");

    /*-------------------------------------------------------------*/
    // read the band_info file
    if ((fp = fopen(bandinfo_filename, "r")) == NULL) {
        fprintf(stderr,
                "-E- %s:%d: unable to open %s for reading\n", __FILE__, __LINE__, bandinfo_filename);
        exit(1);
    }

    // discard the first line of column labels
    fgets(line, BANDINFO_LINE_MAX, fp);
    for (band = 0; band < numBands; band++) {
        if (!fgets(line, BANDINFO_LINE_MAX, fp)) {
            fprintf(stderr, "-E- %s line %d: unable to read band %d from file %s\n",
                    __FILE__, __LINE__, band, bandinfo_filename);
            exit(1);
        }
        count = sscanf(line, "%d %d %d %d %d %d %d %f %f", &dummy,
                &switch_water[band], &lower_water[band], &upper_water[band],
                &switch_land[band], &lower_land[band], &upper_land[band],
                &theoWL[band], &theoE0[band]);
        if (count != 9) {
            fprintf(stderr,
                    "-E- %s line %d: unable to read band %d line from file %s, count = %d\n",
                    __FILE__, __LINE__, band, bandinfo_filename, count);
            exit(1);
        }

        // now change the indexes to 0 based array indexes
        lower_land[band]--;
        upper_land[band]--;
        lower_water[band]--;
        upper_water[band]--;

    } // for band
    fclose(fp);

    if (want_verbose)
        printf("\nSmile corrections enabled.\n\n");

}

/* ------------------------------------------------------------------------ */
/*                                                                          */
/* smile_calc_delta()                                                       */
/*   calculates the smile delta for all of the bands of a pixel             */
/*                                                                          */
/* int *shouldCorrect  true if the reflectance correction should be made    */
/* int *indexes1       index of lower band to use for interpolation         */
/* int *indexes2       upper band to use for intrepolation                  */
/* float *radiances    original measurment                                  */
/* float *theoretWLs   theoretical wavelengths for each band                */
/* float *theoretE0s   theoretical sun spectral fluxes for each band        */
/* float *detectorWLs  detector wavelengths for each band                   */
/* float *detectorE0s  sun spectral fluxes for the detector wavelengths     */
/* float *delta        resulting smile delta correction                     */
/*                                                                          */
/* D. Shea, SAIC, Jan 2009.                                                 */

/* ------------------------------------------------------------------------ */
void smile_calc_delta(int *shouldCorrect,
        int *indexes1,
        int *indexes2,
        float *radiances,
        float *theoretWLs,
        float *theoretE0s,
        float fsol,
        float *detectorWLs,
        float *detectorE0s,
        float *delta) {
    double r0, r1, r2, rc, dl, dr;
    int i0, i1, i2;

    for (i0 = 0; i0 < numBands; i0++) {

        // perform irradiance correction
        r0 = radiances[i0] / detectorE0s[i0];
        rc = r0 * theoretE0s[i0] * fsol;
        delta[i0] = rc - radiances[i0];
        if (shouldCorrect[i0]) {
            // perform reflectance correction
            i1 = indexes1[i0];
            i2 = indexes2[i0];
            r1 = radiances[i1] / detectorE0s[i1];
            r2 = radiances[i2] / detectorE0s[i2];
            dl = (theoretWLs[i0] - detectorWLs[i0]) / (detectorWLs[i2] - detectorWLs[i1]);
            dr = (r2 - r1) * dl * theoretE0s[i0] * fsol;
            delta[i0] += dr;
            delta[i0] /= 10.0; // put delta into l2gen units
        } else {
            delta[i0] /= 10.0; // put delta into l2gen units
        }
    } // for bands

}


/* ------------------------------------------------------------------------ */
/* radcor()                                                          */
/*     loads smile calibration deltas into the radcor array of l1rec        */
/*                                                                          */
/* l1rec  level 1 record to set radcor in                                   */
/* ip     pixel number used to calculate the smile correction               */
/* land   0=ocean, else=land                                                   */
/*                                                                          */
/*                                                                          */
/* D. Shea, SAIC, Jan 2009.                                                 */

/* ------------------------------------------------------------------------ */
void radcor(l1str *l1rec, int32_t ip, int32_t land, int32_t escorrected) {

    int ib;
    int ipb;
    int *shouldCorrect;
    int correctionPossible;
    int *index1;
    int *index2;
    int detectorIndex;
    static float *Ltemp;

    // not initialized so do nothing
    if (numBands == 0)
        return;

    if (!Ltemp) {
        Ltemp = (float*) allocateMemory(sizeof (float)*numBands, "smile temp Lt array");
    }

    /*-------------------------------------------------------------*/
    // do the correction

    detectorIndex = l1rec->pixdet[ip];
    correctionPossible = ((detectorIndex >= 0) && (detectorIndex < numDetectors));

    if (correctionPossible) {
        if (land) {
            shouldCorrect = switch_land;
            index1 = lower_land;
            index2 = upper_land;
        } else {
            shouldCorrect = switch_water;
            index1 = lower_water;
            index2 = upper_water;
        }

        // correct all bands for ozone
        for (ib = 0; ib < numBands; ib++) {
            ipb = ip * l1rec->l1file->nbands + ib; // use the real nbands here to index properly

            /* Correct for ozone absorption.  We correct for inbound and outbound here, */
            /* then we put the inbound back when computing Lw.                          */
             Ltemp[ib] = l1rec->Lt[ipb] * 10.0;
        } // for ib
        // OLCI and MERIS 'SAFE' format "solar_flux" data are Earth-Sun distance corrected,
        // MERIS N1 are not account for that with fsol, as appropriate
        float fsol = l1rec->fsol;
        if (!escorrected)
            fsol = 1.0;
        smile_calc_delta(shouldCorrect,
                index1,
                index2,
                Ltemp,
                theoWL,
                theoE0,
                fsol,
                &(detectorWL[detectorIndex * numBands]),
                &(detectorE0[detectorIndex * numBands]),
                &(l1rec->radcor[ip * l1rec->l1file->nbands]));

    } else {
        // if correction not possible, ensure radcor delta is zero
        for (ib = 0; ib < numBands; ib++) {
            ipb = ip * l1rec->l1file->nbands + ib;
            l1rec->radcor[ipb] = 0;
        }
    } // correction possible

}
