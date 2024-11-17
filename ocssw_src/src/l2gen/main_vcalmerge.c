/* ========================================================================
 * vcalmerge - reads original and vicarious data from L2 file or files
 *	       and writes valid cross-calibration data into an HDF/netCDF4 file
 *
 * Synopsis:
 *
 *   vcalmerge par=parfile 
 *
 * Description:
 * 
 * Modification history:
 *
 *     Programmer     Organization      Date         Description of change
 *   --------------   ------------    --------       ---------------------
 *   Ewa Kwiatkowska   SAIC        27 August   2003    Original development
 *   Ewa Kwiatkowska   SAIC         6 November 2006    getl1rec change
 *   Ewa Kwiatkowska   SAIC         9 October  2007    complete redevelopment
 *					                                   to extract xcal data
 *					                                   from L2 files
 *
 *   Joel Gales        Futuretech  24 October  2012    Fix sensor attribute bug
 *   Joel Gales        Futuretech  14 June     2013    Add support for NETCDF4
 *
 *   Sean Bailey       Futuretech  23 July     2014    Modified to use the new
 *                                                     calfile_utils and to use
 *                                                     the clo library
 * ======================================================================== */

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <libgen.h>
#include <math.h>

#include <readL2scan.h>

#include "l12_proto.h"
#include <setupflags.h>
#include "calfile_utils.h"
#include <hdf5.h>

int16_t fileID = 0;
#define NOVALIDPIX 111

float ***aots, ***nlws; //arrays to contain the values needed to avg over for inversion

/* -------------------------------------------------------------------- *
 *                              main                                    *
 * -------------------------------------------------------------------- */
int main(int argc, char *argv[]) {

    int32_t iscan = 0; /* input scan number                  */
    int32_t oscan = 0; /* output scan number                 */
    int32_t spix = 0; /* start pixel for subscene process   */
    int32_t epix = -1; /* end pixel for subscene process     */
    int32_t ipix;
    int32_t dpix;
    int32_t sscan = 0; /* start scan for subscene process    */
    int32_t escan = -1; /* end scan for subscene process      */
    int32_t dscan;
    int32_t opix = 0; /* output pixels number               */
    int32_t npix = 0; /* Number of output pixels per scan   */
    int32_t nscan = 0; /* Number of output scans             */
    int32_t npixs = 0;
    int32_t i, j, l, np, Ltinx, vLtinx;
    int status;
    static int opened = 0;

    l2_prod l2_str; /* input data for L2 calibrator         */
    filehandle l2file; /* input metadata for L2 calibrator     */
    int32_t nbands, nprods;
    int8 *iptr;
    idDS ds_id;

    uint32_t flagusemask;
    uint32_t required_mask;

    calstr** calrecArray; /* calrec for each pixel in a scan    */
    static int sensorID;
    int nvars1d = 2;
    char vars1Dnames[nvars1d][32]; // floating point one dimensional data to include in output - lon, lat

    char l2_flags[1000];

    if (argc == 1) {
        l2gen_usage("vcalmerge");
        return 0;
    }

    want_verbose = 0;

    setvbuf(stdout, NULL, _IOLBF, 0);
    setvbuf(stderr, NULL, _IOLBF, 0);

    /* Initialize file handles */
    // cdata_();
    filehandle_init(&l2file);
    msl12_input_init();

    /* Parse input parameters */
    if (msl12_input(argc, argv, "vcalmerge", &l2file) != 0) {
        printf("-E- %s: Error parsing input parameters.\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    if (input->ifile[0][0] == '\0') {
        printf("-E- %s: No L2 file provided, exiting...\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    if (access(input->ifile[0], F_OK) || access(input->ifile[0], R_OK)) {
        printf("-E- %s: Input file '%s' does not exist or cannot open.\n",
                argv[0], input->ifile[0]);
        exit(EXIT_FAILURE);
    }

    printf("\n");

    /*
     * If output exists, check sensorID
     *
     */
    int existingSensorID;
    sensorID = l2file.sensorID;

    /* Save old HDF5 error handler */
    H5E_auto_t old_func;
    void *old_client_data;
    H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);

    /* Turn off error handling */
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    if (Hishdf(input->ofile[0]) == TRUE || H5Fis_hdf5(input->ofile[0]) == TRUE) {
        ds_id = openDS(input->ofile[0]);
        status = readAttr(ds_id, "sensorID", &existingSensorID);
        if (status) {
            printf("-E- %s: Problem reading output file: %s\n", argv[0], input->ofile[0]);
            exit(EXIT_FAILURE);
        }
        endDS(ds_id);
        if (sensorID != existingSensorID) {
            printf("-E- %s: Mixing L2 data from different sensors, %s and %s.\n",
                    argv[0], sensorId2SensorName(sensorID), sensorId2SensorName(existingSensorID));
            exit(EXIT_FAILURE);
        }

    }

    /* Restore previous HDF5 error handler */
    H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);
    l2_str.nrec = 0;
    if (openL2(input->ifile[0], 0x0, &l2_str)) {
        printf("-E- %s: Error reading L2 data %s.\n", argv[0], input->ifile[0]);
        exit(EXIT_FAILURE);
    }


    /*
     * Check for required data sets: mside, detnum, pixnum
     */
    if (l2_str.mside == NULL) {
        printf("-E- %s: mside is missing from the L2 data file %s.\n", argv[0],
                l2_str.filename);
        exit(EXIT_FAILURE);
    }
    if (l2_str.detnum == NULL) {
        printf("-E- %s: detnum is missing from the L2 data file %s.\n", argv[0],
                l2_str.filename);
        exit(EXIT_FAILURE);
    }
    if (l2_str.pixnum == NULL) {
        printf("-E- %s: pixnum is missing from the L2 data file %s.\n", argv[0],
                l2_str.filename);
        exit(EXIT_FAILURE);
    }

    /*
     * Make sure all L2 products required for the xcalibration are present
     */

    Ltinx = -1;
    vLtinx = -1;
    for (j = 0; j < l2_str.nprod; j++) {
        if (strncmp(l2_str.prodname[j], "l2_flags", 8) == 0)
            continue;
        if (strncmp(l2_str.prodname[j], "Lt_", 3) == 0) {
            if (Ltinx < 0)
                Ltinx = j;
        }
        if (strncmp(l2_str.prodname[j], "vLt_", 4) == 0) {
            if (vLtinx < 0)
                vLtinx = j;
        }
    }
    if (Ltinx == -1) {
        printf(
                "-E- %s: Lt TOA calibration data are missing from the L2 data file %s.\n",
                argv[0], l2_str.filename);
        exit(EXIT_FAILURE);
    }
    if (vLtinx == -1) {
        printf(
                "-E- %s: vLt TOA calibration data are missing from the L2 data file %s.\n",
                argv[0], l2_str.filename);
        exit(EXIT_FAILURE);
    }

    nbands = rdsensorinfo(sensorID, 0, NULL, NULL);
    nprods = l2_str.nprod;

    strcpy(vars1Dnames[0], "lon");
    strcpy(vars1Dnames[1], "lat");

    /*
     * Set up the l2_flag masking
     */
    strcpy(l2_flags, l2_str.flagnames);
    l = strlen(l2_flags);
    setupflags(l2_flags, input->flaguse, &flagusemask, &required_mask, &status,l2_str.l2_bits);
    if (status < 0) {
        printf("-E- %s: Error setting up L2 flags.\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    printf("\n");

    /*
     * Open the output file.  It will be created if necessary
     */

    long runCounter = 0;
    char outprod[L1_MAXPROD][32];
    for (i = 0; i < nprods; i++) {
        strcpy(outprod[i], l2_str.prodname[i]);
    }

    printf("\n%s: Appending file: %s.\n", argv[0], input->ifile[0]);

    /* Set the end pixel if it was not set by command argument	        */
    if (l1_input->epixl == -1 || l1_input->epixl > l2_str.nsamp)
        l1_input->epixl = l2_str.nsamp;
    if (l1_input->eline == -1 || l1_input->eline > l2_str.nrec)
        l1_input->eline = l2_str.nrec;
    if (l1_input->spixl < 1)
        l1_input->spixl = 1;
    if (l1_input->sline < 1)
        l1_input->sline = 1;

    spix = MAX(l1_input->spixl - 1, 0);
    epix = MIN(l1_input->epixl - 1, l2_str.nsamp - 1);
    dpix = MAX(l1_input->dpixl, 1);
    sscan = MAX(l1_input->sline - 1, 0);
    escan = MIN(l1_input->eline - 1, l2_str.nrec - 1);
    dscan = MAX(l1_input->dline, 1);

    if (sscan > escan || spix > epix) {
        printf("-E- %s: scan and pixel limits make no sense.\n", argv[0]);
        printf(" start scan  = %d\n", sscan + 1);
        printf(" end   scan  = %d\n", escan + 1);
        printf(" start pixel = %d\n", spix + 1);
        printf(" end   pixel = %d\n", epix + 1);
        exit(EXIT_FAILURE);
    }

    npix = (epix - spix) / dpix + 1;
    nscan = (escan - sscan) / dscan + 1;

    int quarterscans = floor(nscan / 4);

    if ((iptr = (int8 *) malloc(nscan * sizeof (int8))) == NULL) {
        printf("-E- %s: Error allocating memory to the pixel index.\n",
                argv[0]);
        exit(EXIT_FAILURE);
    }

    npixs = 0;
    oscan = 0;
    for (iscan = sscan; iscan <= escan; iscan += dscan) {

        if (readL2(&l2_str, l2_str.fileindex, iscan, vLtinx, NULL)) {
            printf("%s: Error reading L2 data file %s at scan %d.\n", argv[0],
                    l2_str.filename, iscan);
            continue;
        }
        if ((iscan % quarterscans) == 0)
            printf("Selecting valid pixels for scan %d\n", iscan);

        opix = 0;
        l = 0;
        for (ipix = spix; ipix <= epix; ipix += dpix) {

            if (((l2_str.l2_flags[ipix] & flagusemask) == 0)
                    && (l2_str.l2_data[vLtinx][ipix] > 0.0)) {
                l = 1;
                ++npixs;
            }
            ++opix;
        }

        if (opix != npix) {
            printf(
                    "%s: Error: Incorrect number of pixels %d obtained in the scan %d.\n",
                    argv[0], opix, oscan);
        }

        if (l)
            iptr[oscan] = 1;
        else
            iptr[oscan] = 0;
        ++oscan;
    }

    printf("Number of valid cross-calibration pixels %d in file %s\n", npixs,
            input->ifile[0]);

    /*
     * if any valid pixels were found, write them to the output file
     */
    if (npixs > 0) {
        opened = 1;
        ds_id = calfile_open(input->ofile[0], sensorID, npixs, 1, nprods, nvars1d, outprod,
                vars1Dnames, &runCounter, CROSSCAL);
        /*
         * Allocate calrec structure
         */
        if ((calrecArray = (calstr **) malloc(npixs * sizeof (calstr *)))
                == NULL) {
            printf(
                    "-E- : Error allocating memory for crosscal record structures\n");
            exit(EXIT_FAILURE);
        }
        for (j = 0; j < npixs; j++) {
            calrecArray[j] = alloc_calrec(1, nbands, nprods, nvars1d);
        }

        oscan = 0;
        for (iscan = sscan; iscan <= escan; iscan += dscan) {
            np = 0;
            /*
             * Give a little progress...
             */
            if ((iscan % quarterscans) == 0)
                printf("Storing pixels at scan %d\n", iscan);

            if (iptr[oscan++]) {
                if (readL2(&l2_str, l2_str.fileindex, iscan, -1, NULL)) {
                    printf("%s: Error reading L2 data file %s at scan %d.\n",
                            argv[0], l2_str.filename, iscan);
                    continue;
                }

                opix = 0;
                for (ipix = spix; ipix <= epix; ipix += dpix) {

                    if (((l2_str.l2_flags[ipix] & flagusemask) == 0)
                            && (l2_str.l2_data[Ltinx][ipix] > 0.0)
                            && (l2_str.l2_data[Ltinx][ipix] < 10000.)
                            && (l2_str.l2_data[vLtinx][ipix] > 0.0)) {

                        calrecArray[np]->sensorID = sensorID;
                        calrecArray[np]->year = (int16) l2_str.year;
                        calrecArray[np]->day = (int16) l2_str.day;
                        calrecArray[np]->msec = (int32) l2_str.msec;
                        calrecArray[np]->iscan = (int16) iscan;
                        calrecArray[np]->mside = (int32) (l2_str.mside[iscan]);
                        calrecArray[np]->detnum = (int32) (l2_str.detnum[iscan]);
                        calrecArray[np]->pixnum = (int32) (l2_str.pixnum[ipix]);
                        calrecArray[np]->vars1D[0] = l2_str.longitude[ipix];
                        calrecArray[np]->vars1D[1] = l2_str.latitude[ipix];
                        for (j = 0; j < nbands; j++) {
                            calrecArray[np]->Lt[j][0] =
                                    (float) l2_str.l2_data[Ltinx + j][ipix];
                            calrecArray[np]->vLt[j][0] =
                                    (float) l2_str.l2_data[vLtinx + j][ipix];
                        }
                        // The calrec entries for the other request products are populated here
                        for (j = 0; j < nprods; j++) {
                            if (strcmp(l2_str.prodname[j], "pixnum") != 0
                                    && strcmp(l2_str.prodname[j], "mside") != 0
                                    && strncmp(l2_str.prodname[j], "Lt", 2) != 0
                                    && strncmp(l2_str.prodname[j], "vLt", 3) != 0
                                    && strcmp(l2_str.prodname[j], "detnum") != 0
                                    && strcmp(l2_str.prodname[j], "l2_flags") != 0) {

                                // fill up calrec data
                                calrecArray[np]->data[j][0] = l2_str.l2_data[j][ipix];
                            }
                        }

                        ++np;
                    }

                    ++opix;
                }

                if (opix != npix) {
                    /*
                     * this really should never happen...
                     */
                    printf(
                            "%s: Error: Incorrect number of pixels %d obtained in the scan %d.\n",
                            argv[0], opix, oscan);
                }
                if (np > 0) {
                    // Write out the valid pixels for this scan

                    for (j = 0; j < np; j++) {
                        calfile_write(ds_id, calrecArray[j], runCounter, 1, 1,
                                nprods, nbands, nvars1d, outprod, vars1Dnames, CROSSCAL);
                        runCounter++;
                    }
                }
            }
        }

        // Free up some memory
        for (i = 0; i < np; i++) {
            free_calrec(calrecArray[i], nbands, nprods);
        }
    }

    free(iptr);

    if (opened == 1) {
        if (calfile_close(ds_id) != 0) {
            printf("Trouble closing file %s\n", input->ofile[0]);
            exit(EXIT_FAILURE);
        }
        printf("Completed processing file %s\n", input->ifile[0]);
        exit(SUCCESS);
    } else {
        printf("Completed processing file %s, but no valid pixels were found\n", input->ifile[0]);
        exit(NOVALIDPIX);
    }

}

