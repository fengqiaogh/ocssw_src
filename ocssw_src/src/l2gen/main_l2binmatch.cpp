/* ========================================================================
 *  l2binmatch creates a file with L2 sensor validation data
 * 
 *  L2 files for a given time period from the sensor being validated
 *  L3 daily, 8-day, etc. bin file from the sensor used as validation truth
 *
 * ======================================================================== */
//todo implement product selection via the l2prod parameter
//todo implement bandshift option for handling different L2 and L3 sensors
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "l2binmatch_input.h"
#include <l2_wrapper.h>
#include <L3File.h>
#include "calfile_utils.h"
#include <setupflags.h>

#include <vector>

using namespace std;
using namespace l3;

int16_t fileID = 0;

float ***aots, ***nlws; //arrays to contain the values needed to avg over for inversion
instr *input; /* input parameters structure         */

int main(int argc, char* argv[]) {
    long runCounter;
    int32_t i, j, k;
    int32_t iscan = 0; /* input scan number                       */
    int32_t ipix; /* input pixel number                      */
    int32_t npixs; /* number of matched pixels in the L2 file */
    int32_t npixsOut; /* number of matched pixels to output */

    int32_t spix = 0; /* start pixel for subscene process   */
    int32_t epix = -1; /* end pixel for subscene process     */
    int32_t dpix;
    int32_t sscan = 0; /* start scan for subscene process    */
    int32_t escan = -1; /* end scan for subscene process      */
    int32_t dscan;
    int32_t nprod, nbands, status;

    float lon, lat;

    static l2_prod l2_str; /* input data for l2 calibrator         */
    meta_l2Type *meta_l2; /* input metadata for l2 calibrator     */

    int32_t l2sensorID, l3sensorID;

    uint32 flagusemask;
    uint32 required_mask;
    char buf[5000];
    idDS ds_id;

    //    static instr input; /* input parameters structure         */
    calstr* pixrec; /* output detector-run structure      */
    std::vector<calstr*> pixrecArray; /* array of pixrec to output    */
    int nvars1d = 2; // initialize to 2, since lon and lat are always output
    char vars1Dnames[L1_MAXPROD][32]; // floating point one dimensional data to include in output
    char outprod[L1_MAXPROD][32];

    /* hold all of the command line options */
    clo_optionList_t* list = clo_createList();

    if (argc == 1 || strcmp(argv[1], "-version") == 0
            || strcmp(argv[1], "--version") == 0) {
        want_verbose = 0;
        l2binmatch_init_options(list);
        l2binmatch_read_options(list, argc, argv);
        clo_printUsage(list);
        exit(EXIT_FAILURE);
    }
    want_verbose = 1;

    /* Parse input parameters */
    if (l2binmatch_input(argc, argv, list) != 0) {
        printf("-E- %s: Error parsing input parameters.\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    printf("Opening %s\n", input->ifile[0]);
    if (openL2(input->ifile[0], 0x0, &l2_str)) {
        printf("-E- %s: Error reading L2 data %s.\n", argv[0], input->ifile[0]);
        exit(EXIT_FAILURE);
    }

    if ((meta_l2 = (meta_l2Type *) malloc(sizeof (meta_l2Type))) == NULL) {
        printf("-E- %s: Error allocating memory for L2 metadata.\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    if (readL2meta(meta_l2, 0)) {
        printf("-E- %s: Error reading L2 metadata %s.\n", argv[0],
                input->ifile[0]);
        exit(EXIT_FAILURE);
    }

    /* grab the target sensor name from the L2 meta-data */
    l2sensorID = sensorName2SensorId(meta_l2->sensor_name);

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
            printf("-E- %s: Problem reading output file: %s\n", argv[0],
                    input->ofile[0]);
            exit(EXIT_FAILURE);
        }
        endDS(ds_id);
        if (l2sensorID != existingSensorID) {
            printf(
                    "-E- %s: Mixing L2 data from different sensors, %s and %s.\n",
                    argv[0], sensorId2SensorName(l2sensorID),
                    sensorId2SensorName(existingSensorID));
            exit(EXIT_FAILURE);
        }

    }
    /* Restore previous HDF5 error handler */
    H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);



    // open the target L3 file
    L3File* l3File = new L3File();
    if (!l3File->open(input->tgtfile)) {
        printf("-E- Could not open ifile=\"%s\".\n", input->tgtfile);
        exit(EXIT_FAILURE);
    }
    l3sensorID = sensorName2SensorId(
            l3File->getMetaData()->sensor_name);
    if (l3sensorID == -1) {
        l3sensorID = instrumentPlatform2SensorId(
                l3File->getMetaData()->sensor,
                l3File->getMetaData()->mission);
        if (l3sensorID == -1) {
            printf("-E- Unknown sensor name %s\n",
                    l3File->getMetaData()->sensor_name);
            exit(EXIT_FAILURE);
        }
    }

    // Get L3 product list
    nprod = l3File->getHdfBinObject()->nprod();
    char *fullprodlist = (char *) malloc(l3File->getHdfBinObject()->query());

    l3File->getHdfBinObject()->query(fullprodlist);

    l3File->setActiveProductList(fullprodlist);

    //Setup flag mask
    strcpy(buf, l2_str.flagnames);
    setupflags(buf, input->flaguse, &flagusemask, &required_mask, &status);
    if (status < 0) {
        printf("-E- %s: Error reading L2 flags %s.\n", argv[0], input->ifile[0]);
        exit(EXIT_FAILURE);
    }

    // initialize the depth file e.g., ETOPO1
    dem_init();

    nbands = rdsensorinfo(l2sensorID, 0, NULL, NULL);

    // Set up list of products to match with L3 file
    int32_t *l2prodinx;
    if ((l2prodinx = (int32_t *) malloc(nprod * sizeof (int32_t))) == NULL) {
        printf("-E- %s: Error allocating memory for L2 product index.\n",
                argv[0]);
        exit(EXIT_FAILURE);
    }

    const char *l3prod;
    for (i = 0; i < nprod; i++) {
        l3prod = l3File->getHdfBinObject()->getProdName(i);
        for (j = 0; j < l2_str.nprod; j++) {
            if (strcmp(l2_str.prodname[j], l3prod) == 0) {
                l2prodinx[i] = j;
                strcpy(outprod[i], l3prod);
                break;
            }
        }
    }
    // Set up list of products to output as 1D arrays (i.e. NOT matched to a L3 product)
    strcpy(vars1Dnames[0], "lon"); // lon and lat you get by default
    strcpy(vars1Dnames[1], "lat");

    int32_t oneDprodinx[L1_MAXPROD];
    k = 2; //start at 2 since lon and lat are defacto
    for (j = 0; j < l2_str.nprod; j++) {
        int skip = 0;
        for (i = 0; i < nprod; i++) {
            if (strcmp(l2_str.prodname[j], outprod[i]) == 0) {
                skip = 1;
                break;
            }
        }
        if (skip == 0 && strcmp(l2_str.prodname[j], "pixnum") != 0
                && strcmp(l2_str.prodname[j], "mside") != 0
                && strcmp(l2_str.prodname[j], "detnum") != 0
                && strcmp(l2_str.prodname[j], "l2_flags") != 0
                && strcmp(l2_str.prodname[j], "latitude") != 0
                && strcmp(l2_str.prodname[j], "longitude") != 0) {
            oneDprodinx[k++] = j;
            strcpy(vars1Dnames[nvars1d++], l2_str.prodname[j]);
        }
    }

    /* Set the end pixel if it was not set by command argument          */
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

    for (iscan = sscan; iscan <= escan; iscan += dscan) {

        if ((iscan % 100) == 0)
            printf("Processing scan %d of %d\n", iscan, l2_str.nrec);

        if (readL2(&(l2_str), 0, iscan, -1, NULL)) {
            printf("%s: Error reading L2 data file %s at scan %d.\n", argv[0],
                    l2_str.filename, iscan);
            continue;
        }
        for (ipix = spix; ipix <= epix; ipix += dpix) {
            // Make sure pixels pass the flags defined as masks in flaguse
            if ((l2_str.l2_flags[ipix] & flagusemask) == 0) {

                lat = l2_str.latitude[ipix];
                lon = l2_str.longitude[ipix];

                L3Bin* l3Bin = l3File->getClosestBin(lat, lon);

                if (l3Bin) {
                    /* check:
                     *  depth is greater than input depth
                     *  minimum nscenes is met
                     *  minimum nsamples is met
                     */

                    if (input->vcal_depth > get_dem(lat, lon)
                            && l3Bin->getNscenes() >= input->vcal_min_nscene
                            && l3Bin->getNobs() >= input->vcal_min_nbin) {
                        pixrec = alloc_calrec(2, nbands, nprod, nvars1d);
                        pixrec->sensorID = l2sensorID;
                        pixrec->year = l2_str.year;
                        pixrec->day = l2_str.day;
                        pixrec->msec = l2_str.msec;
                        pixrec->iscan = iscan;
                        pixrec->mside = l2_str.mside[iscan];
                        pixrec->detnum = l2_str.detnum[iscan];
                        if (l2_str.pixnum)
                            pixrec->pixnum = l2_str.pixnum[ipix];
                        else
                            pixrec->pixnum = ipix;
                        pixrec->vars1D[0] = lon;
                        pixrec->vars1D[1] = lat;
                        // populate the "2D" matched pixel arrays
                        for (j = 0; j < nprod; j++) {
                            pixrec->data[j][0] = l2_str.l2_data[l2prodinx[j]][ipix];
                            pixrec->data[j][1] = l3Bin->getMean(j);
                        }
                        // populate the "1D" arrays
                        for (j = 2; j < nvars1d; j++) {
                            pixrec->vars1D[j] = l2_str.l2_data[oneDprodinx[j]][ipix];
                        }

                        pixrecArray.push_back(pixrec);
                    }
                }
            }
        }
    }
    npixs = pixrecArray.size();
    printf("Number of pixels matched: %d\n", npixs);

    if (npixs > 2 * input->subsamp) {
        j = 0;
        npixsOut = npixs / input->subsamp;
        ds_id = calfile_open(input->ofile[0], l2sensorID, npixsOut, 2, nprod, nvars1d,
                outprod, vars1Dnames, &runCounter, BINMATCH);
        for (i = 0; i < npixs; i += input->subsamp) {
            if (j < npixsOut) {
                pixrec = pixrecArray[i];
                calfile_write(ds_id, pixrec, runCounter++, npixsOut, 2, nprod, nbands,
                        nvars1d, outprod, vars1Dnames, BINMATCH);
                free_calrec(pixrec, nbands, nprod);
                j++;
            }
        }

        if (calfile_close(ds_id) != 0) {
            printf("Trouble closing file %s\n", input->ofile[0]);
            exit(EXIT_FAILURE);
        } else {
            printf("Number of pixels output: %d\n", npixsOut);
            printf("Done processing file %s\n", input->ofile[0]);
        }
    } // Free up some memory

    closeL2(&l2_str, 0);
    freeL2meta(meta_l2);

    freeL2(&l2_str);

    exit(EXIT_SUCCESS);

}
