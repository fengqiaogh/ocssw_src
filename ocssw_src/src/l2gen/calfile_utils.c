/* ========================================================================
 *  Procedures create, read, and write detector/mirror side pixel data
 *
 *     Programmer     Organization      Date         Description of change
 *   --------------   ------------    --------       ---------------------
 *   Sean Bailey      Futuretech       23 July 2014    Original development
 *                                                     based on l1det_hdf.c
 * ======================================================================== */

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <libgen.h>
#include <math.h>

#include <netcdf.h>
#include "passthebuck.h"
#include "filehandle.h"

#include "dfutils.h"
#include "l2prod_struc.h"
#include "calfile_utils.h"
#include "l12_proto.h"

/**
 * det_open
 * opens an output file, if it doesn't exist it creates one, if it does exist it appends
 * by copying existing file to a temporary location and bumping the size of the arrays to
 * create a new output file of the same name
 *
 * @param ofile  - filename to open for writing
 * @param sensorID  - sensor ID number
 * @param ydim   - size of y-dimension for output arrays
 * @param xdim   - size of x-dimension for output arrays
 * @param nprods - number of output products
 * @param nvars1d - number of 1D output products
 * @param l2prods - output product names
 * @param vars1Dnames - 1D product names
 * @param numExistingRecords - number of existing records, used for appending...
 * @param ctype - calibration type (e.g. DET2DET, CROSSCAL, BINMATCH)
 * @return
 */
idDS calfile_open(char *ofile, int sensorID, int ydim, int xdim, int nprods,
        int nvars1d, char l2prods[L1_MAXPROD][32], char vars1Dnames[L1_MAXPROD][32],
        long* numExistingRecords, caltype ctype) {
    idDS ds_id;
    int status;
    int i, j;
    *numExistingRecords = 0;
    int32 dm[3];
    int start[3] = {0, 0, 0};
    int stride[3] = {1, 1, 1};
    int count[3] = {1, 1, 0};
    char name[MAX_NC_NAME];

    ds_id = openDS(ofile);
    if (ds_id.fid == FAIL) {
        status = calfile_create(ofile, &ds_id, sensorID, ydim, xdim, nprods,
                nvars1d, l2prods, vars1Dnames, ctype);
        if (status) {
            printf("Whoops! problem creating output file %s\n", ofile);
            exit(EXIT_FAILURE);
        }
    } else {
        //get the number of runs in the existing file
        status = getDimsDS(ds_id, "fileID", dm);

        if (status) {
            printf(
                    "Sorry, Charlie, output file %s already exists, but does not contain fileID array...try a new output file\n",
                    ofile);
            exit(EXIT_FAILURE);
        }
        // Bump the fileID
        start[0] = dm[0] - 1;
        int16_t fileIDold;
        DPTB(readDS(ds_id, "fileID", start, stride, count, &fileIDold));
        fileID = fileIDold + 1;

        *numExistingRecords = dm[0];
        // close the old file so we can move it to ".old"
        status = endDS(ds_id);

        char *oldfile = (char *) malloc(strlen(ofile) + 5);
        strcpy(oldfile, ofile);
        strcat(oldfile, ".old");
        char cmd[2048];
        sprintf(cmd, "mv %s %s", ofile, oldfile);
        system(cmd);
        // bump nruns by numExistingRecords and create new ofile
        ydim += *numExistingRecords;
        status = calfile_create(ofile, &ds_id, sensorID, ydim, xdim, nprods,
                nvars1d, l2prods, vars1Dnames, ctype);
        if (status) {
            printf("Whoops! problem creating output file %s\n", ofile);
            exit(EXIT_FAILURE);
        }

        // open old file
        idDS old_ds_id;
        old_ds_id = openDS(oldfile);
        if (old_ds_id.fid == FAIL) {
            printf("Whoops! problem opening old file %s\n", oldfile);
            exit(EXIT_FAILURE);
        }

        start[0] = 0;
        count[0] = dm[0];
        count[1] = xdim;
        VOIDP data = malloc(count[0] * count[1] * sizeof (float));

        // loop through old file and write records to new ofile

        // copy filename global attributes
        for (i = 0; i < fileIDold + 1; i++) {
            sprintf(name, "filename%03d", i);
            char* val = readAttrStr(old_ds_id, name);
            if (val) {
                DPTB(SetChrGA(ds_id, name, val));
                free(val);
            }
        }

        // copy data
        DPTB(readDS(old_ds_id, "fileID", start, stride, count, data));
        DPTB(writeDS(ds_id, "fileID", data, 0, 0, 0, count[0], count[1],
                count[2]));

        DPTB(readDS(old_ds_id, "year", start, stride, count, data));
        DPTB(writeDS(ds_id, "year", data, 0, 0, 0, count[0], count[1],
                count[2]));

        DPTB(readDS(old_ds_id, "day", start, stride, count, data));
        DPTB(writeDS(ds_id, "day", data, 0, 0, 0, count[0], count[1], count[2]));

        DPTB(readDS(old_ds_id, "msec", start, stride, count, data));
        DPTB(writeDS(ds_id, "msec", data, 0, 0, 0, count[0], count[1],
                count[2]));

        DPTB(readDS(old_ds_id, "iscan", start, stride, count, data));
        DPTB(writeDS(ds_id, "iscan", data, 0, 0, 0, count[0], count[1],
                count[2]));

        for (i = 0; i < nvars1d; i++) {
            DPTB(readDS(old_ds_id, vars1Dnames[i], start, stride, count, data));
            DPTB(writeDS(ds_id, vars1Dnames[i], data, 0, 0, 0, count[0], count[1],
                    count[2]));
        }
        //        DPTB(readDS(old_ds_id, "lon", start, stride, count, data));
        //        DPTB(
        //                writeDS(ds_id, "lon", data, 0, 0, 0, count[0], count[1],
        //                        count[2]));
        //
        //        DPTB(readDS(old_ds_id, "lat", start, stride, count, data));
        //        DPTB(
        //                writeDS(ds_id, "lat", data, 0, 0, 0, count[0], count[1],
        //                        count[2]));

        DPTB(readDS(old_ds_id, "mside", start, stride, count, data));
        DPTB(writeDS(ds_id, "mside", data, 0, 0, 0, count[0], count[1],
                count[2]));

        DPTB(readDS(old_ds_id, "pixnum", start, stride, count, data));
        DPTB(writeDS(ds_id, "pixnum", data, 0, 0, 0, count[0], count[1],
                count[2]));

        if (ctype == CROSSCAL || ctype == BINMATCH) {
            DPTB(readDS(old_ds_id, "detnum", start, stride, count, data));
            DPTB(writeDS(ds_id, "detnum", data, 0, 0, 0, count[0], count[1],
                    count[2]));
        }

        for (i = 0; i < nprods; i++) {
            int skip = 0;
            for (j = 0; j < nvars1d; j++) {
                if (strcmp(l2prods[i], vars1Dnames[j]) == 0) {
                    skip = 1;
                    break;
                }
            }
            if (skip == 0 && strcmp(l2prods[i], "pixnum") != 0
                    && strcmp(l2prods[i], "mside") != 0
                    && strcmp(l2prods[i], "detnum") != 0
                    && strcmp(l2prods[i], "l2_flags") != 0) {
                DPTB(readDS(old_ds_id, l2prods[i], start, stride, count, data));
                DPTB(writeDS(ds_id, l2prods[i], data, 0, 0, 0, count[0],
                        count[1], count[2]));
            }
        }
        // close old file and remove it.
        free(data);
        status = endDS(old_ds_id);
        sprintf(cmd, "rm -f %s", oldfile);
        system(cmd);

    }
    return ds_id;
}

/**
 * calfile_close - close an open output file
 * @param ds_id
 * @return
 */
int calfile_close(idDS ds_id) {
    int status;
    status = endDS(ds_id);
    return status;
}

/**
 * calfile_create - creates a output file
 * @param ofile  - filename to open for writing
 * @param ds_id  - HDF/netCDF dataset id (idDS)
 * @param sensorID  - sensorID
 * @param ydim   - size of y-dimension for output arrays
 * @param xdim   - size of x-dimension for output arrays
 * @param nprods - number of output products
 * @param nvars1d - number of 1D output products
 * @param l2prods - output product names
 * @param vars1Dnames - 1D product names
 * @param ctype - calibration type (e.g. DET2DET, CROSSCAL)
 * @return
 */
int calfile_create(char *ofile, idDS *ds_id, int sensorID, int ydim, int xdim,
        int nprods, int nvars1d, char l2prods[L1_MAXPROD][32],
        char vars1Dnames[L1_MAXPROD][32], caltype ctype) {

    char title[255];
    int i, j;
    int dumdim;
    int status;
    char name[MAX_NC_NAME];
    int32 dm[3];
    const char dm_name[3][80];
    dm[0] = ydim;
    dm[1] = xdim;
    switch (ctype) {

    case DET2DET:
        strcpy((char *) dm_name[0], "runs");
        strcpy((char *) dm_name[1], "detectors");
        sprintf(title,
                "%s Level-1 detector/mirror-side runs of consecutive pixels",
                sensorId2SensorName(sensorID));
        break;
    case CROSSCAL:
        strcpy((char *) dm_name[0], "pixels");
        strcpy((char *) dm_name[1], "xdim");
        sprintf(title, "%s Level-2 cross-calibration pixels",
                sensorId2SensorName(sensorID));
        break;
    case BINMATCH:
        strcpy((char *) dm_name[0], "pixels");
        strcpy((char *) dm_name[1], "L2L3");
        sprintf(title, "%s Level-2/3 matched pixels",
                sensorId2SensorName(sensorID));
        break;
    default:
        fprintf(stderr, "Unknown caltype! \n");
        exit(EXIT_FAILURE);
        break;
    }

    if (strcmp(input->oformat, "netCDF4") == 0) {
        *ds_id = startDS(ofile, DS_NCDF, DS_WRITE, input->deflate);
        for (i = 0; i < 2; i++) {
            status = nc_def_dim((*ds_id).fid, dm_name[i], dm[i], &dumdim);
            if (status != NC_NOERR) {
                fprintf(stderr, "%s\n", nc_strerror(status));
                exit(EXIT_FAILURE);
            }
        }
    } else {
        *ds_id = startDS(ofile, DS_HDF, DS_WRITE, 0);
    }

    /*                                                                  */
    /* Create the SDSes                                                 */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    PTB(SetChrGA(*ds_id, "title", title));
    PTB(SetI32GA(*ds_id, "sensorID", (int32) sensorID));
    sprintf(name, "filename%03d", fileID);
    PTB(SetChrGA(*ds_id, name, basename(input->ifile[0])));

    PTB(createDS(*ds_id, (int) sensorID, "fileID", dm, dm_name));
    PTB(createDS(*ds_id, (int) sensorID, "year", dm, dm_name));
    PTB(createDS(*ds_id, (int) sensorID, "day", dm, dm_name));
    PTB(createDS(*ds_id, (int) sensorID, "msec", dm, dm_name));

    PTB(createDS(*ds_id, (int) sensorID, "iscan", dm, dm_name));
    //    PTB(createDS(*ds_id, (int ) sensorID, "lon", dm, dm_name));
    //    PTB(createDS(*ds_id, (int ) sensorID, "lat", dm, dm_name));
    PTB(createDS(*ds_id, (int) sensorID, "mside", dm, dm_name));
    PTB(createDS(*ds_id, (int) sensorID, "pixnum", dm, dm_name));

    for (i = 0; i < nvars1d; i++) {
        PTB(createDS(*ds_id, (int) sensorID, vars1Dnames[i], dm, dm_name));
    }

    if (ctype == CROSSCAL || ctype == BINMATCH) {
        PTB(createDS(*ds_id, (int) sensorID, "detnum", dm, dm_name));
    }

    for (i = 0; i < nprods; i++) {
        int skip = 0;
        for (j = 0; j < nvars1d; j++) {
            if (strcmp(l2prods[i], vars1Dnames[j]) == 0) {
                skip = 1;
                break;
            }
        }
        if (skip == 0 && strcmp(l2prods[i], "pixnum") != 0
                && strcmp(l2prods[i], "mside") != 0
                && strcmp(l2prods[i], "detnum") != 0
                && strcmp(l2prods[i], "l2_flags") != 0) {
            PTB(
                    createDS(*ds_id, (int) sensorID, l2prods[i], dm,
                    dm_name));
        }
    }

    return (LIFE_IS_GOOD);

}

/**
 * calfile_write() - writes data to output file
 * @param ds_id  - HDF/netCDF dataset id (idDS)
 * @param calrec - structure containing data to write
 * @param recnum - record number for record to write
 * @param ydim   - size of y-dimension for output arrays
 * @param xdim   - size of x-dimension for output arrays
 * @param nprods - number of output products
 * @param nbands - number of wavelength bands
 * @param nvars1d - number of 1 dimensional arrays
 * @param l2prods - output product names
 * @param vars1Dnames - 1D product names
 * @param ctype - calibration type (e.g. DET2DET, CROSSCAL)
 * @return
 */
int calfile_write(idDS ds_id, calstr *calrec, int recnum, int ydim, int xdim,
        int nprods, int nbands, int nvars1d, char l2prods[L1_MAXPROD][32],
        char vars1Dnames[L1_MAXPROD][32], caltype ctype) {

    int32 i, j;
    static int firstWrite = 1;
    static l2prodstr** prodptr;
    VOIDP *pbuf;
    int iwave[nbands + NBANDSIR], *iwave_p;

    static productInfo_t **v1Dp_info;
    static productInfo_t **p_info;

    if (firstWrite) {
        if ((v1Dp_info = (productInfo_t **) malloc(nvars1d * sizeof (struct productInfo_t *))) == NULL) {
            printf("%s -Error: Cannot allocate memory to productInfo data\n",
                    __FILE__);
            exit(FATAL_ERROR);
        }
        if ((p_info = (productInfo_t **) malloc(nprods * sizeof (struct productInfo_t *))) == NULL) {
            printf("%s -Error: Cannot allocate memory to productInfo data\n",
                    __FILE__);
            exit(FATAL_ERROR);
        }
        rdsensorinfo(calrec->sensorID, 0, "iwave", (void **) &iwave_p);
        for (i = 0; i < nbands + NBANDSIR; i++) {
            iwave[i] = iwave_p[i];
        }
        if ((prodptr = (l2prodstr **) malloc(L1_MAXPROD * sizeof (l2prodstr *)))
                == NULL) {
            printf(
                    "-E- : Error allocating memory for product info record structures\n");
            exit(EXIT_FAILURE);
        }
        for (i = 0; i < nprods; i++) {
            p_info[i] = allocateProductInfo();
            if (!findProductInfo(l2prods[i], calrec->sensorID, p_info[i])) {
                printf("%s not found in XML product table\n", l2prods[i]);
                exit(EXIT_FAILURE);
            }
            if ((prodptr[i] = get_l2prod_index(l2prods[i], calrec->sensorID,
                    nbands, xdim, ydim, iwave)) == NULL) {
                fprintf(stderr, "-E- %s line %d: product index failure.\n",
                        __FILE__, __LINE__);
                return (1);
            };
        }
        for (i = 0; i < nvars1d; i++) {
            v1Dp_info[i] = allocateProductInfo();
            if (!findProductInfo(vars1Dnames[i], calrec->sensorID, v1Dp_info[i])) {
                printf("%s not found in XML product table\n", vars1Dnames[i]);
                exit(EXIT_FAILURE);
            }
            j = i + nprods;
            if ((prodptr[j] = get_l2prod_index(vars1Dnames[i], calrec->sensorID,
                    nbands, xdim, ydim, iwave)) == NULL) {
                fprintf(stderr, "-E- %s line %d: product index failure.\n",
                        __FILE__, __LINE__);
                return (1);
            };
        }

        firstWrite = 0;
    }
    PTB(writeDS(ds_id, "fileID", (VOIDP) & fileID, recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "year", (VOIDP) & (calrec->year), recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "day", (VOIDP) & (calrec->day), recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "msec", (VOIDP) & (calrec->msec), recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "iscan", (VOIDP) & (calrec->iscan), recnum, 0, 0, 1, 1, 1));
    for (i = 0; i < nvars1d; i++) {
        j = i + nprods;
        if ((v1Dp_info[i])->scaleFactor == 1.0 && (v1Dp_info[i])->addOffset == 0.0)
            pbuf = (VOIDP) & calrec->vars1D[i];
        else
            pbuf = scale_sds(&calrec->vars1D[i], v1Dp_info[i], 1);
        PTB(writeDS(ds_id, vars1Dnames[i], pbuf, recnum, 0, 0, 1, 1, 1));
    }
    for (i = 0; i < xdim; i++) {
        PTB(writeDS(ds_id, "pixnum", (VOIDP) & (calrec->pixnum), recnum, i, 0, 1, 1, 1));
    }
    PTB(writeDS(ds_id, "mside", (VOIDP) & (calrec->mside), recnum, 0, 0, 1, 1, 1));

    if (ctype == CROSSCAL || ctype == BINMATCH) {
        PTB(writeDS(ds_id, "detnum", (VOIDP) & (calrec->detnum), recnum, 0, 0, 1, 1, 1));
    }
    int Ltix = 0;
    int vLtix = 0;
    for (i = 0; i < nprods; i++) {
        int skip = 0;
        for (j = 0; j < nvars1d; j++) {
            if (strcmp(l2prods[i], vars1Dnames[j]) == 0) {
                skip = 1;
                break;
            }
        }
        if (skip == 0 && strcmp(l2prods[i], "pixnum") != 0
                && strcmp(l2prods[i], "mside") != 0
                && strcmp(l2prods[i], "detnum") != 0
                && strcmp(l2prods[i], "l2_flags") != 0) {

            if (strncmp(l2prods[i], "Lt", 2) == 0) {
                pbuf = (VOIDP) calrec->Lt[Ltix++];
            } else if (strncmp(l2prods[i], "vLt", 3) == 0) {
                pbuf = (VOIDP) calrec->vLt[vLtix++];
            } else {
                if ((p_info[i])->scaleFactor == 1.0 && (p_info[i])->addOffset == 0.0)
                    pbuf = (VOIDP) calrec->data[i];
                else
                    pbuf = scale_sds(calrec->data[i], p_info[i], xdim);
            }
            PTB(writeDS(ds_id, l2prods[i], pbuf, recnum, 0, 0, 1, xdim, 1));
        }
    }

    return (LIFE_IS_GOOD);

}

/**

 */

/**
 * inversion_init
 * Initialize the average nLw and aot arrays for input into inversion
 * @param ndets  - number of contiguous detectors
 * @param iscan  - scan number of l1 file
 * @param nbands - number of wavelength bands
 * @param ipix   - pixel number
 * @param aerec  - aerosol calibration record structure
 * @param tgrec  - nlw calibration record structure
 */
void inversion_init(long ndets, long iscan, int nbands, long ipix, aestr *aerec,
        tgstr *tgrec) {
    int band;
    long kscan;
    float avg_nLw;
    float avg_aot;

    for (band = 0; band < nbands; band++) {
        avg_nLw = 0;
        avg_aot = 0;
        for (kscan = iscan; kscan < iscan + ndets; kscan++) {
            avg_nLw += nlws[kscan][ipix][band];
            avg_aot += aots[kscan][ipix][band];
        }
        avg_nLw /= ndets;
        avg_aot /= ndets;
        if (avg_nLw < 0)
            avg_nLw = 0.0;
        if (avg_aot < 0)
            avg_aot = 0.0;

        tgrec->nLw[ipix * nbands + band] = avg_nLw;
        aerec->taua[ipix * nbands + band] = avg_aot;
    }
}

/**
 * alloc_calrec
 * Allocate space for a detector run record
 *
 * @param ydim - y-dimension of records to allocate
 * @param nbands - number of wavelength bands
 * @param nprods - number of output products (2D)
 * @param nvar1d - number of 1 dimensional arrays in output product
 * @return - allocated calrec structure
 */
calstr* alloc_calrec(int ydim, int nbands, int nprods, int nvar1d) {

    calstr* record = calloc(1, sizeof (calstr));
    int band, prod;

    if ((record->Lt = (float **) malloc(nbands * sizeof (float *))) == NULL) {
        printf("%s -Error: Cannot allocate memory to detector data\n",
                __FILE__);
        exit(FATAL_ERROR);
    }
    if ((record->vLt = (float **) malloc(nbands * sizeof (float *))) == NULL) {
        printf("%s -Error: Cannot allocate memory to detector data\n",
                __FILE__);
        exit(FATAL_ERROR);
    }

    if ((record->data = (float **) malloc(nprods * sizeof (float *))) == NULL) {
        printf("%s -Error: Cannot allocate memory to detector data\n",
                __FILE__);
        exit(FATAL_ERROR);
    }

    if ((record->vars1D = (float *) calloc(nvar1d, sizeof (float))) == NULL) {
        printf("%s -Error: Cannot allocate memory to 1D data array\n",
                __FILE__);
        exit(FATAL_ERROR);
    }


    for (band = 0; band < nbands; band++) {
        if ((record->Lt[band] = (float *) calloc(ydim, sizeof (float)))
                == NULL) {
            printf("%s -Error: Cannot allocate memory to detector data\n",
                    __FILE__);
            exit(FATAL_ERROR);
        }
        if ((record->vLt[band] = (float *) calloc(ydim, sizeof (float)))
                == NULL) {
            printf("%s -Error: Cannot allocate memory to detector data\n",
                    __FILE__);
            exit(FATAL_ERROR);
        }

    }
    for (prod = 0; prod < nprods; prod++) {
        if ((record->data[prod] = (float *) calloc(ydim, sizeof (float)))
                == NULL) {
            printf("%s -Error: Cannot allocate memory to detector data\n",
                    __FILE__);
            exit(FATAL_ERROR);
        }
    }
    return record;
}

/**
 * free_calrec - free allocated calibration structures
 * @param calrec
 * @param nbands - number of wavelength bands
 * @param nprods - number of output products (2D)
 * @param nvars1d - number of 1 dimensional arrays in output product
 */
void free_calrec(calstr *calrec, int nbands, int nprods) {
    int band, prod;
    for (band = 0; band < nbands; band++) {
        free(calrec->Lt[band]);
        free(calrec->vLt[band]);
    }
    for (prod = 0; prod < nprods; prod++) {
        free(calrec->data[prod]);
    }
    free(calrec->Lt);
    free(calrec->vLt);
    free(calrec->data);
    free(calrec->vars1D);

    free(calrec);
}

/**
 * free_prodptr - free array of product info structures
 * @param p
 */
void free_prodptr(l2prodstr **p) {
    int i;
    for (i = 0; i < L1_MAXPROD; i++) {
        free(p[i]);
    }
    free(p);
}

