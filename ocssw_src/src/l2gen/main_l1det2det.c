/* =====================================================================*/
/*                                                                      */
/* Program: l1det2det - multi-sensor detector-to-detector calibration   */
/*                                                                      */
/*                                                                      */
/* Written By:                                                          */
/*                                                                      */
/*     Ewa J. Kwiatkowska                                               */
/*     SAIC                                                             */
/*     NASA/Ocean Color Project                                         */
/*     June 2007                                                        */
/*                                                                      */
/*     Programmer     Organization      Date       Description of change*/
/*   --------------   ------------    --------     ---------------------*/
/*   Joel Gales       Futuretech      08/01/13     Add support for      */
/*                                                 NETCDF4              */
/*   B. Franz         NASA            09/10/13     Remove references to */
/*                                                 resolution and 1640  */
/*                                                 and 500m MODIS code  */
/*   Sean Bailey      Futuretech      04/15/13     Entirely reworked to */
/*                                                 handle any mission,  */
/*                                                 work with existing   */
/*                                                 product definitions  */
/*                                                 and require complete */
/*                                                 scans                */
/* =====================================================================*/

#include <stdio.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <libgen.h>
#include <setupflags.h>
#include "l12_proto.h"
#include "calfile_utils.h"
#include <hdf5.h>
#include <l2_wrapper.h>

int read9km_mask(char *file, char *mask);
int is_masked(int32_t bin, char *mask, int32_t nrows);
int16_t fileID = 0;

float ***aots, ***nlws; //arrays to contain the values needed to avg over for inversion

/* -------------------------------------------------------------------- */
/*                            main                                      */

/* -------------------------------------------------------------------- */
int main(int argc, char* argv[]) {

    static int32_t nrows = 2160;
    long iscan = 0; /* input scan number                  */
    long npix = 0; /* input number pixels per scan       */
    long spix = 0; /* start pixel for subscene process   */
    long epix = -1; /* end pixel for subscene process     */
    long sscan = 0; /* start scan for subscene process    */
    long escan = -1; /* end scan for subscene process      */
    long nscan;
    long nbands;
    int32_t ndets;
    long nruns; /* number of output detector runs     */
    long ipix;
    long jscan; /* loop variables                     */
    int band;
    long scanidx;
    long runidx;
    int initial_mirror_side;
    int goodscan;
    idDS ds_id;

    l1str l1rec; /* generic level-1b scan structure    */
    l2str l2rec; /* generic level-2  scan structure    */
    tgstr tgrec; /* structure to store target values   */
    aestr aerec; /* structure to store aerosol values  */
    calstr** detrecArray; /* detrec for each pixel in a scan    */

    filehandle l1file; /* input l1 file handle               */

    double start_time;
    long i, j, ib, runCounter, nir_l, nprods;
    int32 *numbin, *basebin;
    int32_t bin, ibin, bin_row, bin_col;
    long **valid_pixel_array; //Valid pixel array 0=invalid, 1=pixel valid
    int numContiguousDetectorRuns; // number of scans for detector runs
    int status;

    char *mask, *maskfile, maskname[1000];
    float lat, lon, latbin;
    int32_t flagusemask;

    if (argc == 1) {
        l2gen_usage("l1det2det");
        return 0;
    }

    want_verbose = 1;

    setvbuf(stdout, NULL, _IOLBF, 0);
    setvbuf(stderr, NULL, _IOLBF, 0);

    /* Initialize file handles */
    //cdata_();
    filehandle_init(&l1file);
    msl12_input_init();

    /* Parse input parameters */
    if (msl12_input(argc, argv, "l1det2det", &l1file) != 0) {
        printf("-E- %s: Error parsing input parameters.\n", argv[0]);
        exit(FATAL_ERROR);
    }

    if (access(input->ifile[0], F_OK) || access(input->ifile[0], R_OK)) {
        printf("-E- %s: Input file '%s' does not exist or cannot open.\n",
                argv[0], input->ifile[0]);
        exit(FATAL_ERROR);
    }

    /*                                                                  */
    /* Open input file and get sensor   */
    /* and scan information from handle */
    /*                                                                  */
    strcpy(l1file.name, input->ifile[0]);
    l1file.geofile = input->geofile;

    if (openl1(&l1file) != 0) {
        printf("-E- %s: Error opening %s for reading.\n",
                argv[0], l1file.name);
        exit(FATAL_ERROR);
    }

    /*
     * Check output file for matching sensorID
     */
    int existingSensorID;

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
        if (l1file.sensorID != existingSensorID) {
            printf("-E- %s: Mixing L2 data from different sensors, %s and %s.\n",
                    argv[0], sensorId2SensorName(l1file.sensorID), sensorId2SensorName(existingSensorID));
            exit(EXIT_FAILURE);
        }

    }
    /* Restore previous HDF5 error handler */
    H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);

    npix = l1file.npix;
    numContiguousDetectorRuns = input->ybox;
    if (numContiguousDetectorRuns < 1)
        numContiguousDetectorRuns = 1;

    if (numContiguousDetectorRuns > 2) {
        printf("ybox > 2, sorry, only letting you use 2 - talk to management if you want more...\n");
        numContiguousDetectorRuns = 2;
    }
    /*                                                                  */
    /* Allocate memory for L1 and L2 scan data                          */
    /*                                                                  */
    if (alloc_l1(&l1file, &l1rec) == 0) {
        printf("-E- %s: Unable to allocate L1 record.\n", argv[0]);
        exit(FATAL_ERROR);
    }
    if (alloc_l2(&l1rec, &l2rec) == 0) {
        printf("-E- %s: Unable to allocate L2 record.\n", argv[0]);
        exit(FATAL_ERROR);
    }

    /* Set the end pixel if it was not set by command argument          */
    if (l1_input->epixl == -1 || l1_input->epixl > l1file.npix)
        l1_input->epixl = l1file.npix;
    if (l1_input->eline == -1 || l1_input->eline > l1file.nscan)
        l1_input->eline = l1file.nscan;
    if (l1_input->spixl < 1)
        l1_input->spixl = 1;
    if (l1_input->sline < 1)
        l1_input->sline = 1;

    spix = MAX(l1_input->spixl - 1, 0);
    epix = MIN(l1_input->epixl - 1, l1file.npix - 1);
    sscan = MAX(l1_input->sline - 1, 0);
    escan = MIN(l1_input->eline - 1, l1file.nscan - 1);

    ndets = l1file.ndets;
    // ensure there are full scans
    if (sscan > 0) {
        while (sscan % ndets)
            sscan++;
    }
    int numContiguousDets = ndets * numContiguousDetectorRuns;
    escan = ((escan + 1) / numContiguousDets) * numContiguousDets - 1;
    //escan = ((escan + 1) / ndets) * ndets - 1;

    if (sscan > escan || spix > epix) {
        printf("-E- %s: scan and pixel limits make no sense.\n", argv[0]);
        printf(" start scan  = %ld\n", sscan + 1);
        printf(" end   scan  = %ld\n", escan + 1);
        printf(" start pixel = %ld\n", spix + 1);
        printf(" end   pixel = %ld\n", epix + 1);
        exit(FATAL_ERROR);
    }

    /* Note: for the L1 file, npix is still the native scan pixel count */
    l1file.spix = spix;
    l1file.epix = epix;

    nscan = (escan - sscan) + 1;
    npix = (epix - spix) + 1;
    nbands = l1file.nbands;

    /*                                                                  */
    /* Transfer any additional header info to the L2 record headers     */
    /*                                                                  */
    l2rec.bindx = l1file.bindx;

    /*                                                                  */
    /* Set up parameters for strict ocean color forward processing      */
    /*                                                                  */
    input->proc_ocean = 1;
    input->proc_land = 0;
    input->proc_sst = 0;
    input->atmocor = 1;
    input->aer_opt = AERRH;
    l1_input->outband_opt = 0;
    input->mode = FORWARD;
    char outprod[L1_MAXPROD][32];
    int nvars1d = 2;
    char vars1Dnames[nvars1d][32]; // floating point one dimensional data to include in output - lon, lat
    nprods = prodlist(l1file.sensorID, 0, (const char *) &(input->l2prod), (const char *) &(input->def_l2prod), outprod);

    /*
     * Allocating memory for the multidimentional arrays of pixel index, nLw and aot
     */
    if ((valid_pixel_array = (long **) malloc(nscan * sizeof (long *))) == NULL) {
        printf("-E- %s: Error allocating memory to the pixel index.\n", argv[0]);
        exit(FATAL_ERROR);
    }

    if ((aots = (float ***) malloc(nscan * sizeof (float *))) == NULL) {
        printf("-E- : Error allocating memory to in situ inputs\n");
        exit(FATAL_ERROR);
    }
    if ((nlws = (float ***) malloc(nscan * sizeof (float *))) == NULL) {
        printf("-E- : Error allocating memory to in situ inputs\n");
        exit(FATAL_ERROR);
    }

    for (i = 0; i < nscan; i++) {
        if ((valid_pixel_array[i] = (long *) calloc(npix, sizeof (long))) == NULL) {
            printf("-E- %s: Error allocating memory to the pixel index.\n", argv[0]);
            exit(FATAL_ERROR);
        }

        if ((aots[i] = (float **) malloc(npix * sizeof (float *))) == NULL) {
            printf("-E- %s: Error allocating memory to the AOT data.\n", argv[0]);
            exit(FATAL_ERROR);
        }
        if ((nlws[i] = (float **) malloc(npix * sizeof (float *))) == NULL) {
            printf("-E- %s: Error allocating memory to the nLw data.\n", argv[0]);
            exit(FATAL_ERROR);
        }
        for (j = 0; j < npix; j++) {
            if ((aots[i][j] = (float *) calloc(nbands, sizeof (float))) == NULL) {
                printf("-E- %s: Error allocating memory to the aot data.\n", argv[0]);
                exit(FATAL_ERROR);
            }
            if ((nlws[i][j] = (float *) calloc(nbands, sizeof (float))) == NULL) {
                printf("-E- %s: Error allocating memory to the nLw data.\n", argv[0]);
                exit(FATAL_ERROR);
            }
        }
    }

    // Allocate memory for and set up deepwater mask
    if ((mask = (char *) malloc(5940424 * sizeof (char))) == NULL) {
        printf("-E- %s: Error allocating memory to the mask data\n", argv[0]);
        exit(FATAL_ERROR);
    }
    if ((maskfile = getenv("OCDATAROOT")) == NULL) {
        printf("-E- %s: OCDATAROOT directory undefined.\n", argv[0]);
        exit(FATAL_ERROR);
    }
    strcpy(maskname, maskfile);
    strcat(maskname, "/common/deep_water_mask_9km.dat\x0");

    if (read9km_mask(maskname, mask) > 0) {
        printf("-E- %s: Failed reading the deep water mask file\n", argv[0]);
        exit(FATAL_ERROR);
    }

    /* ----------------- */
    /* Set up bin arrays */
    /* ----------------- */
    if ((numbin = (int32 *) calloc(nrows, sizeof (int32))) == NULL) {
        printf("-E- %s: Error allocating memory to numbin\n", argv[0]);
        exit(FATAL_ERROR);
    }
    if ((basebin = (int32 *) calloc(nrows + 1, sizeof (int32))) == NULL) {
        printf("-E- %s: Error allocating memory to basebin\n", argv[0]);
        exit(FATAL_ERROR);
    }

    bin = 0;
    for (i = 0; i < nrows; i++) {
        latbin = (i + 0.5) * (180.0 / nrows) - 90.0;
        numbin[i] = (int32) (cos(latbin * OEL_PI / 180.0) * (2.0 * nrows) + 0.5);
        bin += numbin[i];
    }
    basebin[0] = 1;
    for (i = 1; i <= nrows; i++) {
        basebin[i] = basebin[i - 1] + numbin[i - 1];
    }

    // NOTE: this hard-coded string will need to be updated if we ever change the l2gen flag definitions
    char buf[1024];
    uint32 required;
    strcpy(buf, "ATMFAIL,LAND,PRODWARN,HIGLINT,HILT,HISATZEN,COASTZ,SPARE,STRAYLIGHT,CLDICE,COCCOLITH,TURBIDW,HISOLZEN,SPARE,LOWLW,CHLFAIL,NAVWARN,ABSAER,SPARE,MAXAERITER,MODGLINT,CHLWARN,ATMWARN,OPSHAL,SEAICE,NAVFAIL,FILTER,SPARE,BOWTIEDEL,HIPOL,PRODFAIL,SPARE");
    setupflags(buf, input->flaguse, (uint32 *) & flagusemask, &required, &status);
    printf("flagusemask (%d): %d",status,flagusemask);
    printf("input flags: %s",input->flaguse);
    strcpy(vars1Dnames[0], "lon");
    strcpy(vars1Dnames[1], "lat");

    nir_l = bindex_get(input->aer_wave_long);

    printf("\n\nBegin %s Processing\n", "l1det2det\x0");
    printf("\nSensor is %s\n", sensorId2SensorName(l1file.sensorID));
    printf("Sensor ID is %d\n", l1file.sensorID);
    printf("Sensor has %d reflective bands\n", l1file.nbands);
    printf("Sensor has %d emissive bands\n", l1file.nbandsir);
    printf("Sensor resolution %d m\n", l1_input->resolution);
    printf("Number of along-track detectors per band is %d\n", ndets);
    printf("Number of input pixels per scan is %d\n", l1file.npix);
    printf("Processing pixels %ld to %ld\n", spix + 1, epix + 1);
    printf("Processing scans %ld to %ld\n", sscan + 1, escan + 1);

    start_time = now();
    printf("\nBegin l1det2det processing at %s\n\n", ydhmsf(start_time, 'L'));

    /*
     *  Read file scan by scan, convert to L2
     *  determine valid detector "runs" and
     *  populate the nLw and tau arrays needed
     *  for averaging to insert into inverse
     *  processing.
     */
    for (iscan = sscan; iscan <= escan; iscan++) {

        scanidx = iscan - sscan;

        if (getl1rec(iscan, 1, &l1rec) != 0) {
            exit(FATAL_ERROR);
        }
        // First Line!
        if (iscan == sscan) {
            initial_mirror_side = l1rec.mside;
        }

        if ((iscan % 50) == 0)
            printf("Processing scan #%6ld (%ld of %ld) after %6.0f seconds\n",
                iscan, iscan - sscan + 1, escan - sscan + 1, now() - start_time);

        /*                                                              */
        /* Convert the L1B radiances to L2                              */
        /*                                                              */
        convl12(&l1rec, &l2rec, 0, l1rec.npix - 1, NULL);

        /*                                                      */
        /* Loop through each pixel and check exclusion criteria */
        /*  on pixel flags and values                           */
        /*                                                      */
        for (ipix = 0; ipix < l1rec.npix; ipix++) {
            ibin = 1;

            lat = l1rec.lat[ipix];
            lon = l1rec.lon[ipix];
            bin_row = (int32_t) ((90.0 + lat) * ((float64) nrows / (float64) 180.0));
            bin_col = (int32_t) ((float64) numbin[bin_row] * (lon + 180.0) / (float64) 360.0);
            bin = basebin[bin_row] + bin_col;
            ibin = (int32_t) is_masked(bin, mask, nrows);

            if (((l1rec.flags[ipix] & flagusemask) == 0) &&
                    l2rec.chl[ipix] > 0.0 &&
                    l2rec.chl[ipix] <= input->chlthreshold &&
                    l2rec.taua[ipix * nbands + nir_l] > 0.0 &&
                    l2rec.taua[ipix * nbands + nir_l] <= input->aotthreshold && ibin > 0) {

                valid_pixel_array[scanidx][ipix] = 1;
                for (ib = 0; ib < nbands; ib++) {
                    nlws[scanidx][ipix][ib] = l2rec.nLw [ipix * nbands + ib];
                    aots[scanidx][ipix][ib] = l2rec.taua[ipix * nbands + ib];
                }
            }
            /*
             * Once a detector run is complete,
             * check the all lines in the scan to ensure all detectors are valid
             */
            if (l1rec.detnum == ndets - 1) {
                if (numContiguousDetectorRuns == 1 || l1rec.mside != initial_mirror_side) {
                    for (i = 0; i < ndets * numContiguousDetectorRuns; i++) {
                        if (!valid_pixel_array[scanidx - i][ipix]) {
                            for (j = 0; j < ndets * numContiguousDetectorRuns; j++) {
                                valid_pixel_array[scanidx - j][ipix] = 0;
                            }
                            break;
                        }
                    }

                }
            }
        }
    }

    printf("\nEnd of line-by-line processing at %s\n", ydhmsf(now(), 'L'));
    printf("Processing Rate = %f scans/sec\n", nscan / (now() - start_time));

    // Free memory associated with deepwater mask
    free(mask);
    free(numbin);
    free(basebin);

    /*
     * Count the number of valid detector runs (i.e. scans with all detectors valid)
     */
    nruns = 0;
    for (iscan = sscan; iscan <= escan; iscan += ndets) {
        scanidx = iscan - sscan;
        for (ipix = 0; ipix < l1rec.npix; ipix++) {
            if (valid_pixel_array[scanidx][ipix])
                nruns++;
        }
    }

    if (!nruns) {
        printf("-W- %s: No analysis pixel data found in the L2 file %s\n", argv[0], input->ifile[0]);
    } else {

        printf("\n\nNow extracting valid pixel runs\n");
        printf("Along-track pixel runs of %d found: %ld\n", numContiguousDets, nruns);

        start_time = now();

        /*
         * Allocate detrec, aer and tgt structures
         */
        if ((detrecArray = (calstr **) malloc(l1rec.npix * sizeof (calstr *))) == NULL) {
            printf("-E- : Error allocating memory for detector record structures\n");
            exit(EXIT_FAILURE);
        }
        for (i = 0; i < l1rec.npix; i++) {
            detrecArray[i] = alloc_calrec(ndets, nbands, nprods, nvars1d);
        }


        alloc_aer(l1rec.npix, nbands, &aerec);
        alloc_target(l1rec.npix, nbands, &tgrec);
        /*
         * set up the bits for inversion
         */
        input->aer_opt = FIXAOT;
        input->mode = INVERSE_NLW;
        aerec.mode = ON;
        tgrec.mode = ON;
        l2rec.tgrec = &tgrec;

        /*
         * Open the output file.  It will be created if necessary
         */

        runCounter = 0;
        ds_id = calfile_open(input->ofile[0], l1file.sensorID, nruns, ndets, nprods, nvars1d, outprod, vars1Dnames, &runCounter, DET2DET);

        /*
         * Fill up detrec with detector info and grab avg aot/nlws
         * for use in inversion to fill in vLt array in detrec
         *  along with the other
         * detector values...
         */

        runCounter--; //nruns index
        int runNum = 0; //just for nice printing on output...

        for (iscan = sscan; iscan <= escan; iscan += ndets * numContiguousDetectorRuns) {
            goodscan = 0;
            runidx = iscan - sscan;
            for (ipix = 0; ipix < l1rec.npix; ipix++) {
                if (valid_pixel_array[runidx][ipix]) {
                    goodscan = 1;
                    break;
                }
            }
            if (goodscan) {
                // For "good scans", loop through detectors for the valid runs and get the pixel info for detrec
                for (jscan = iscan; jscan < iscan + ndets * numContiguousDetectorRuns; jscan++) {
                    scanidx = jscan - sscan;
                    int detidx = jscan - iscan;
                    if (detidx >= ndets)
                        detidx -= ndets;

                    if (getl1rec(jscan, 1, &l1rec) != 0) {
                        exit(EXIT_FAILURE);
                    }
                    // For the first line in the scan, fill up the per run values in detrec
                    if (detidx == 0) {
                        int16_t year, day;
                        double sec;
                        unix2yds(l1rec.scantime, &year, &day, &sec);

                        for (ipix = 0; ipix < l1rec.npix; ipix++) {
                            if (valid_pixel_array[scanidx][ipix]) {
                                detrecArray[ipix]->sensorID = l1rec.l1file->sensorID;
                                detrecArray[ipix]->year = (int32_t) year;
                                detrecArray[ipix]->day = (int32_t) day;
                                detrecArray[ipix]->msec = (int32_t) (sec * 1.e3);
                                detrecArray[ipix]->mside = l1rec.mside;
                                detrecArray[ipix]->iscan = (int16_t) jscan;
                                detrecArray[ipix]->pixnum = l1rec.pixnum[ipix];
                                detrecArray[ipix]->vars1D[0] = l1rec.lon[ipix];
                                detrecArray[ipix]->vars1D[1] = l1rec.lat[ipix];
                                // average the nLws and aots over the pixels in ndets*numContiguousDetectorRuns
                                // and fills in the aerec and tgrec arrays needed for the inversion step
                                inversion_init(ndets*numContiguousDetectorRuns, runidx, nbands, ipix, &aerec, &tgrec);
                            }
                        }
                    }
                    // The Lt array is filled and the avg nLw and taua needed for the inversion
                    // are set up in this loop
                    for (ipix = 0; ipix < l1rec.npix; ipix++) {
                        if (valid_pixel_array[scanidx][ipix]) {
                            // average the nLws and aots over the pixels in ndets*numContiguousDetectorRuns
                            // and fills in the aerec and tgrec arrays needed for the inversion step
                            //                            inversion_init(ndets*numContiguousDetectorRuns, runidx, l1rec.nbands, ipix, &aerec, &tgrec);
                            for (band = 0; band < nbands; band++) {
                                ib = ipix * nbands + band;
                                detrecArray[ipix]->Lt[band][detidx] = (float) l1rec.Lt[ib];
                            }
                        }
                    }
                    /*                                                              */
                    /* Convert the L1B radiances to L2                              */
                    /*                                                              */
                    convl12(&l1rec, &l2rec, 0, l1rec.npix - 1, &aerec);
                    /*                                                              */
                    /* Convert back to L1B to get vLts                              */
                    /*                                                              */
                    convl21(&l2rec, &tgrec, 0, l1rec.npix - 1, l1rec.Lt, NULL);

                    // set detrec data and vLts obtained via inversion step
                    for (ipix = 0; ipix < l1rec.npix; ipix++) {
                        if (valid_pixel_array[scanidx][ipix]) {
                            // The vLts are filled here
                            for (band = 0; band < nbands; band++) {
                                ib = ipix * nbands + band;
                                detrecArray[ipix]->vLt[band][detidx] = (float) l1rec.Lt[ib];
                            }
                            // The detrec entries for the other request products are populated here
                            for (i = 0; i < nprods; i++) {
                                if (strcmp(outprod[i], "pixnum") != 0 &&
                                        strcmp(outprod[i], "mside") != 0 &&
                                        strncmp(outprod[i], "Lt", 2) != 0 &&
                                        strncmp(outprod[i], "vLt", 3) != 0 &&
                                        strcmp(outprod[i], "detnum") != 0 &&
                                        strcmp(outprod[i], "l2_flags") != 0) {
                                    // get the product index record
                                    l2prodstr *p;

                                    if ((p = get_l2prod_index(outprod[i], l1rec.l1file->sensorID,
                                            nbands + l1file.nbandsir, l1rec.npix, l1rec.l1file->nscan, l1rec.l1file->iwave)) == NULL) {
                                        fprintf(stderr,
                                                "-E- %s line %d: product index failure.\n",
                                                __FILE__, __LINE__);
                                        return (1);
                                    };

                                    // extract the product
                                    float *pbuf;
                                    pbuf = prodgen(p, &l2rec);

                                    // fill up detrec data
                                    detrecArray[ipix]->data[i][detidx] = pbuf[ipix];
                                }
                            }
                        }
                    }

                    // Write out the valid runs for this scan
                    if (detidx == ndets - 1) {
                        for (ipix = 0; ipix < l1rec.npix; ipix++) {
                            if (valid_pixel_array[scanidx][ipix]) {
                                runCounter++;
                                runNum++;
                                if ((runNum % 500) == 0)
                                    printf("Writing run %d of %ld\n", runNum, nruns);
                                calfile_write(ds_id, detrecArray[ipix], runCounter, nruns, ndets, nprods, nbands, nvars1d, outprod, vars1Dnames, DET2DET);
                            }
                        }
                    }
                }
            }
        }
        printf("\nEnd pixel-run extraction at %s, pixel-run number found: %ld\n\n\n", ydhmsf(now(), 'L'), nruns);
        if (calfile_close(ds_id) != 0) {
            printf("Trouble closing file %s\n", input->ofile[0]);
            exit(EXIT_FAILURE);
        }
        // Free up some memory
        for (i = 0; i < l1rec.npix; i++) {
            free_calrec(detrecArray[i], nbands, nprods);
        }

    }

    for (ib = 0; ib < nbands; ib++) {
        free(aots[ib]);
        free(nlws[ib]);
    }
    free(aots);
    free(nlws);

    /*                                                                  */
    /* Close all files                                                  */
    /*                                                                  */
    closel1(&l1file);

    if (nruns > 0) return (EXIT_SUCCESS);
    else return (NOMATCH_ERROR);
}
