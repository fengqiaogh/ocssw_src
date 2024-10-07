/* ========================================================================
 *  Procedures create, read, and write mscalibration pixel data
 *
 * int crosscal_append(char *xcalfile, mscalstr calstr) - creates file, if
 *			it does not exist already, and appends data
 * int crosscal_read(char *xcalfile, int32_t subsmpl, mscalstr *calstr) - 
 *			reads subsampled pixel cross-calibration data and
 * 			allocates memory and data to the calstr structure
 * int crosscal_add(char *xcalfile, int32_t subsmpl, mscalstr *calstr, int32_t *npixs, int32_t *ngranuls) 
 *                      reads cross-calibration data from an hdf file,
 *          		subsamples if needed, and adds pixel and granule
 *          		data to an already existing calstr structure
 *
 *     Programmer     Organization      Date       Description of change
 *   --------------   ------------    --------     ---------------------
 *   Ewa Kwiatkowska  SAIC           10 September 2003    Original development
 *   Joel Gales       Futuretech     24 October   2012    Add ds_id.fftype
 *   Joel Gales       Futuretech     14 June      2013    Add support for NETCDF4 output
 *   Joel Gales       Futuretech     31 July      2013    Remove dead code
 *
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

#include "l12_proto.h" // Must be first to find netcdf.h
#include "passthebuck.h"
#include "l2prod_struc.h"
#include "filehandle.h"
#include "hdf_utils.h"
#include "mfhdf.h"
#include "mscal_struc.h"

static int32_t nbands;

//int rdstrideSDS(int32 fileID, char *sdsname, int32 stride[2], VOIDP array_data);
int read_glbl_attr(idDS ds_id, char *name, VOIDP ptr);
int crosscal_create(char *xcalfile, idDS *ds_id, mscalstr calstr, int32_t npixs);
int crosscal_readblocks(char *xcalfile, idDS *ds_id, int32_t *spix, int32_t *totalpixs, mscalstr *calstr);
int crosscal_writeblocks(char *xcalfile, idDS ds_id, int32_t *spix, mscalstr calstr, int32_t nfiles);

int crosscal_create(char *xcalfile, idDS *ds_id, mscalstr calstr, int32_t npixs) {

    char title[255];
    l2prodstr p, *ptr;
    int j;
    int dumdim;
    int32 dm[3];
    const char dm_name[3][80];

    sprintf(title, "%s Level-2 cross-calibration pixels", sensorName[calstr.sensorID]);

    if (strcmp(calstr.oformat, "netCDF4") == 0) {
        *ds_id = startDS(xcalfile, DS_NCDF, DS_WRITE, 0);
        if (nc_def_dim((*ds_id).fid, "Number_of_Pixels", npixs, &dumdim)
                != NC_NOERR) exit(1);
    } else {
        *ds_id = startDS(xcalfile, DS_HDF, DS_WRITE, 0);
    }

    dm[0] = npixs;
    strcpy((char *) dm_name[0], "Number of Pixels");
    /*                                                                  */
    /* Create the pixel SDSes                                           */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    PTB(SetChrGA(*ds_id, "title", title));
    PTB(SetI32GA(*ds_id, "sensorID\0", (int32) calstr.sensorID));
    PTB(SetChrGA(*ds_id, "Input Parameters", calstr.input_parms));

    PTB(createDS(*ds_id, (int) calstr.sensorID, "fileID", dm, dm_name));

    PTB(createDS(*ds_id, (int) calstr.sensorID, "year", dm, dm_name));
    PTB(createDS(*ds_id, (int) calstr.sensorID, "day", dm, dm_name));
    PTB(createDS(*ds_id, (int) calstr.sensorID, "msec", dm, dm_name));
    PTB(createDS(*ds_id, (int) calstr.sensorID, "iscan", dm, dm_name));

    PTB(createDS(*ds_id, (int) calstr.sensorID, "mside", dm, dm_name));
    PTB(createDS(*ds_id, (int) calstr.sensorID, "detnum", dm, dm_name));
    PTB(createDS(*ds_id, (int) calstr.sensorID, "pixnum", dm, dm_name));

    PTB(createDS(*ds_id, (int) calstr.sensorID, "clon", dm, dm_name));
    PTB(createDS(*ds_id, (int) calstr.sensorID, "clat", dm, dm_name));

    /*                                                                  */
    /* Create the geophysical SDSes for the L2 xcalibration products    */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    ptr = &p;
    for (j = 0; j < calstr.nprods; j++) {

        if (strcmp(calstr.l2prods[j], "l2_flags") == 0);
        else
            if (strcmp(calstr.l2prods[j], "mside") == 0);
        else
            if (strcmp(calstr.l2prods[j], "detnum") == 0);
        else
            if (strcmp(calstr.l2prods[j], "pixnum") == 0);
        else {

            if ((ptr = get_l2prod_index((char *) calstr.l2prods[j], (int32) calstr.sensorID,
                    (int32) calstr.nbands, 0, 1, calstr.Lambda)) == NULL) {

                fprintf(stderr,
                        "-E- %s line %d: product index failure.\n",
                        __FILE__, __LINE__);
                return (FATAL_ERROR);
            };

            PTB(createDS(*ds_id, (int) calstr.sensorID, calstr.l2prods[j], dm, dm_name));
        }
    }


    return (LIFE_IS_GOOD);

}

int crosscal_append(char *xcalfile, mscalstr calstr) {

    idDS ds_id_old, ds_id;
    char command[2 * FILENAME_MAX], tempname[FILENAME_MAX];
    int exists;
    int32_t spix_old, spix, total_old, i;
    mscalstr calstr_old;

    nbands = NBANDS; /* calstr.nbands; */

    spix_old = -1;
    total_old = 0;
    calstr_old.data = NULL;
    calstr_old.npixs = calstr.npixs;
    strcpy(calstr_old.oformat, calstr.oformat);
    if (crosscal_readblocks(xcalfile, &ds_id_old, &spix_old, &total_old, &calstr_old) != LIFE_IS_GOOD) {
        printf("-E- %s line %d: Cannot read already existing SDS data %s\n", __FILE__, __LINE__, xcalfile);
        return (HDF_FUNCTION_ERROR);
    }

    if (calstr_old.npixs == 0L) exists = 0;
    else {
        exists = 1;
        if (calstr_old.sensorID != calstr.sensorID) {
            printf("-E- %s line %d: Attempt to write into a file whose data come from a different sensor %s\n", __FILE__, __LINE__, sensorName[calstr_old.sensorID]);
            free_calstr(calstr_old, 1);
            endDS(ds_id_old);
            return (HDF_FUNCTION_ERROR);
        }
    }

    if (exists) {

        /*
                p = (char *)name;
                strcpy(tempname, xcalfile);
                p = dirname(tempname);
                strcpy(tempname, p);
                strcat(tempname, "/temp.hdf");
         */
        strcpy(tempname, xcalfile);
        strcat(tempname, "_tmp");

        if (crosscal_create(tempname, &ds_id, calstr, total_old + calstr.npixs) != LIFE_IS_GOOD) {
            printf("-E- %s line %d: Cannot read already existing SDS data %s\n", __FILE__, __LINE__, tempname);
            free_calstr(calstr_old, 1);
            endDS(ds_id);
            endDS(ds_id_old);
            return (HDF_FUNCTION_ERROR);
        }

        spix = 0L;

        do {

            if (crosscal_writeblocks(tempname, ds_id, &spix, calstr_old, 0L) != LIFE_IS_GOOD) {
                printf("-E- %s line %d: Cannot write existing SDS data into the temp file %s\n", __FILE__, __LINE__, tempname);
                free_calstr(calstr_old, 1);
                endDS(ds_id);
                endDS(ds_id_old);
                return (HDF_FUNCTION_ERROR);
            }

            if (crosscal_readblocks(xcalfile, &ds_id_old, &spix_old, &total_old, &calstr_old) != LIFE_IS_GOOD) {
                printf("-E- %s line %d: Cannot read already existing SDS data %s\n", __FILE__, __LINE__, xcalfile);
                free_calstr(calstr_old, 1);
                endDS(ds_id);
                return (HDF_FUNCTION_ERROR);
            }

        } while (calstr_old.npixs > 0L);


        for (i = 0; i < calstr.npixs; i++) calstr.fileID[i] += (int16) calstr_old.nfiles;

        if (crosscal_writeblocks(tempname, ds_id, &spix, calstr, calstr_old.nfiles) != LIFE_IS_GOOD) {
            printf("-E- %s line %d: Cannot write new SDS data into the temp file %s\n", __FILE__, __LINE__, tempname);
            endDS(ds_id);
            return (HDF_FUNCTION_ERROR);
        }

        if (endDS(ds_id) == FAIL) {
            printf("-E- %s line %d: Could not close HDF file, %s.\n", __FILE__, __LINE__, tempname);
            return (HDF_FUNCTION_ERROR);
        }

        sprintf(command, "mv %s %s", tempname, xcalfile);
        system(command);

    } else {

        if (crosscal_create(xcalfile, &ds_id, calstr, calstr.npixs) != LIFE_IS_GOOD) {
            printf("-E- %s line %d: Could not open HDF file, %s .\n", __FILE__, __LINE__, xcalfile);
            endDS(ds_id);
            return (HDF_FUNCTION_ERROR);
        }

        spix = 0L;

        if (crosscal_writeblocks(xcalfile, ds_id, &spix, calstr, 0L) != LIFE_IS_GOOD) {
            printf("-E- %s line %d: Cannot write new SDS data into the file %s\n", __FILE__, __LINE__, xcalfile);
            endDS(ds_id);
            return (HDF_FUNCTION_ERROR);
        }

        if (endDS(ds_id) == FAIL) {
            printf("-E- %s line %d: Could not close HDF file, %s.\n", __FILE__, __LINE__, xcalfile);
            return (HDF_FUNCTION_ERROR);
        }
    }


    return (LIFE_IS_GOOD);

}



/* ------------------------------------------------------         */
/* crosscal_readblocks() - reads consecutive cross-calibration    */
/*          blocks of data of the size of 400,000 pixels       	  */
/*          calptr comprises the actual data for the 400,0000     */
/*          pixels, memory is allocated with the first use        */

/* ------------------------------------------------------         */

int crosscal_readblocks(char *xcalfile, idDS *ds_id, int32_t *spix, int32_t *totalpixs, mscalstr *calstr) {

    int32 dim_sizes[H4_MAX_VAR_DIMS], sds_id;
    char name[H4_MAX_NC_NAME];
    int32_t nfiles = 0, i, length, l;
    int32 n_datasets, n_file_attr, rank, num_type, attributes;
    char input_parms[16384];
    int *varids;

    if (*spix >= *totalpixs) {
        if (endDS(*ds_id) == FAIL) {
            printf("-E- %s line %d: Could not close HDF file, %s .\n", __FILE__, __LINE__, xcalfile);
            return (HDF_FUNCTION_ERROR);
        }
        free_calstr(*calstr, 1);
        calstr->npixs = 0L;
        return (LIFE_IS_GOOD);
    }


    if (*spix < 0L) {

        calstr->npixs = 0L;
        calstr->nfiles = 0L;

        if (strcmp(calstr->oformat, "netCDF4") == 0) {
            *ds_id = startDS(xcalfile, DS_NCDF, DS_READ, 0);
        } else {
            *ds_id = startDS(xcalfile, DS_HDF, DS_READ, 0);
        }
        if ((*ds_id).fid == FAIL) {
            return (LIFE_IS_GOOD);
        }

        if (strcmp(calstr->oformat, "netCDF4") == 0) {
            nc_inq_varids((*ds_id).fid, &n_datasets, NULL);
            varids = (int *) calloc(n_datasets, sizeof (int));
            nc_inq_varids((*ds_id).fid, &n_datasets, varids);
        } else {
            SDfileinfo((*ds_id).fid, &n_datasets, &n_file_attr);
        }

        strcpy(name, "sensorID\0");
        if (read_glbl_attr(*ds_id, name, (VOIDP) &(calstr->sensorID)) != LIFE_IS_GOOD) {
            printf("-E- %s line %d: Could not read HDF sensor attribute, %s .\n", __FILE__, __LINE__, xcalfile);
            endDS(*ds_id);
            return (HDF_FUNCTION_ERROR);
        }

        strcpy(name, "Input Parameters\0");
        read_glbl_attr(*ds_id, name, (VOIDP) input_parms);
        length = strlen(input_parms);
        if ((calstr->input_parms = (char *) malloc((length + 1) * sizeof (char))) == NULL) {
            printf("-E- %s line %d: Error allocating memory to MScalmerge input parameter text.\n", __FILE__, __LINE__);
            endDS(*ds_id);
            return (HDF_FUNCTION_ERROR);
        }
        strncpy(calstr->input_parms, input_parms, length);
        calstr->input_parms[length] = '\x0';

        nfiles = 0L;
        do {
            sprintf(name, "filename%d", nfiles);

            if (findAttr(*ds_id, name) == FAIL) break;
            else ++nfiles;

        } while (1);


        strcpy(name, "fileID\0");
        if (getDimsDS(*ds_id, name, dim_sizes) != NC_NOERR) {
            printf("-E- %s line %d: Could not read file dimensions, %s .\n", __FILE__, __LINE__, xcalfile);
            endDS(*ds_id);
            return (HDF_FUNCTION_ERROR);
        }
        if (dim_sizes[0] < 0) {
            printf("-E- %s line %d: The rank of the requested parameter (%s) is incorrect\n", __FILE__, __LINE__, name);
            endDS(*ds_id);
            return (HDF_FUNCTION_ERROR);
        }
        *totalpixs = (int32_t) dim_sizes[0];

        if (nfiles == 0L || *totalpixs == 0L) {
            printf("-E- %s line %d: There are no data in the file, %s .\n", __FILE__, __LINE__, name);
            endDS(*ds_id);
            return (LIFE_IS_GOOD);
        }

        calstr->nprods = 1000;
        if ((calstr->l2prods = (prname *) malloc(calstr->nprods * sizeof (prname))) == NULL) {
            printf("-E- %s line %d: Error allocating memory to l2 product names.\n", __FILE__, __LINE__);
            exit(FATAL_ERROR);
        }
        calstr->nprods = 0;

        for (i = 0; i < n_datasets; i++) {

            if (strcmp(calstr->oformat, "netCDF4") == 0) {
                nc_inq_varname((*ds_id).fid, varids[i], name);
            } else {
                sds_id = SDselect((*ds_id).fid, i);
                SDgetinfo(sds_id, name, &rank, dim_sizes, &num_type, &attributes);
            }

            if (strcmp(name, "fileID") == 0);
            else
                if (strcmp(name, "l2_flags") == 0);
            else
                if (strcmp(name, "year") == 0);
            else
                if (strcmp(name, "day") == 0);
            else
                if (strcmp(name, "msec") == 0);
            else
                if (strcmp(name, "iscan") == 0);
            else
                if (strcmp(name, "mside") == 0);
            else
                if (strcmp(name, "detnum") == 0);
            else
                if (strcmp(name, "pixnum") == 0);
            else
                if (strcmp(name, "lon") == 0);
            else
                if (strcmp(name, "lat") == 0);
            else {

                strcpy(calstr->l2prods[calstr->nprods], name);
                calstr->nprods++;
            }
        }
        free(varids);

        if ((calstr->l2prods = (prname *) realloc((void *) calstr->l2prods, calstr->nprods * sizeof (prname))) == NULL) {
            printf("-E- %s line %d: Error reallocating memory to l2 product names.\n", __FILE__, __LINE__);
            return (HDF_FUNCTION_ERROR);
        }
        calstr->Lambda = NULL;

        if (*totalpixs > 1000000)
            length = alloc_calstr(nfiles, 1000000L, calstr);
        else
            length = alloc_calstr(nfiles, *totalpixs, calstr);


        for (i = 0; i < nfiles; i++) {
            sprintf(name, "filename%d", i);

            if (read_glbl_attr(*ds_id, name, (VOIDP) calstr->filenames[i]) != SUCCESS) {
                printf("-E- %s line %d: Error reading filename attributes (%s)\n", __FILE__, __LINE__, name);
                endDS(*ds_id);
                free_calstr(*calstr, 1);
                calstr->npixs = 0L;
                calstr->nfiles = 0L;
                return (HDF_FUNCTION_ERROR);
            }
        }

        *spix = 0L;

    } else {

        if (*spix + calstr->npixs > *totalpixs) calstr->npixs = *totalpixs - *spix;
    }
    int32 start[3] = {*spix, 0, 0};
    int32 count[3] = {calstr->npixs, 0, 0};
    int32 stride[3] = {1, 1, 1};

    PTB(readDS(*ds_id, "fileID\0", start, stride, count, (VOIDP) calstr->fileID));
    PTB(readDS(*ds_id, "year\0", start, stride, count, (VOIDP) calstr->year));
    PTB(readDS(*ds_id, "day\0", start, stride, count, (VOIDP) calstr->day));
    PTB(readDS(*ds_id, "msec\0", start, stride, count, (VOIDP) calstr->msec));
    PTB(readDS(*ds_id, "iscan\0", start, stride, count, (VOIDP) calstr->iscan));
    PTB(readDS(*ds_id, "mside\0", start, stride, count, (VOIDP) calstr->mside));
    PTB(readDS(*ds_id, "detnum\0", start, stride, count, (VOIDP) calstr->detnum));
    PTB(readDS(*ds_id, "pixnum\0", start, stride, count, (VOIDP) calstr->pixnum));
    PTB(readDS(*ds_id, "lon\0", start, stride, count, (VOIDP) calstr->lon));
    PTB(readDS(*ds_id, "lat\0", start, stride, count, (VOIDP) calstr->lat));

    l = 0;
    for (i = 0; i < calstr->nprods; i++) {
        if (strcmp(calstr->l2prods[i], "l2_flags") == 0);
        else
            if (strcmp(calstr->l2prods[i], "mside") == 0);
        else
            if (strcmp(calstr->l2prods[i], "detnum") == 0);
        else
            if (strcmp(calstr->l2prods[i], "pixnum") == 0);
        else {
            PTB(readDS(*ds_id, calstr->l2prods[i], start, stride, count,
                    (VOIDP)&(calstr->ddata[l * calstr->npixs])));
            ++l;
        }
    }


    *spix += calstr->npixs;

    return (LIFE_IS_GOOD);

}


/* ------------------------------------------------------         */
/* crosscal_writeblocks() - writes consecutive cross-calibration  */
/*          blocks of data                              	  */
/*          calptr comprises the actual data for the 400,0000     */
/*          pixels, memory is allocated with the first use        */

/* ------------------------------------------------------         */

int crosscal_writeblocks(char *xcalfile, idDS ds_id, int32_t *spix, mscalstr calstr, int32_t nfiles) {

    char name[H4_MAX_NC_NAME];
    int32_t i, j, l;

    if (nfiles < 0) nfiles = 0L;

    if (*spix <= 0 || nfiles > 0) {
        for (i = 0; i < calstr.nfiles; i++) {

            j = i + nfiles;
            sprintf(name, "filename%d", j);

            PTB(SetChrGA(ds_id, name, (char *) calstr.filenames[i]));
        }
    }

    PTB(writeDS(ds_id, "fileID\0", (VOIDP) calstr.fileID, (int32) * spix, 0, 0, (int32) calstr.npixs, 1, 1));
    PTB(writeDS(ds_id, "year\0", (VOIDP) calstr.year, (int32) * spix, 0, 0, (int32) calstr.npixs, 1, 1));
    PTB(writeDS(ds_id, "day\0", (VOIDP) calstr.day, (int32) * spix, 0, 0, (int32) calstr.npixs, 1, 1));
    PTB(writeDS(ds_id, "msec\0", (VOIDP) calstr.msec, (int32) * spix, 0, 0, (int32) calstr.npixs, 1, 1));
    PTB(writeDS(ds_id, "iscan\0", (VOIDP) calstr.iscan, (int32) * spix, 0, 0, (int32) calstr.npixs, 1, 1));
    PTB(writeDS(ds_id, "mside\0", (VOIDP) calstr.mside, (int32) * spix, 0, 0, (int32) calstr.npixs, 1, 1));
    PTB(writeDS(ds_id, "detnum\0", (VOIDP) calstr.detnum, (int32) * spix, 0, 0, (int32) calstr.npixs, 1, 1));
    PTB(writeDS(ds_id, "pixnum\0", (VOIDP) calstr.pixnum, (int32) * spix, 0, 0, (int32) calstr.npixs, 1, 1));
    PTB(writeDS(ds_id, "lon\0", (VOIDP) calstr.lon, (int32) * spix, 0, 0, (int32) calstr.npixs, 1, 1));
    PTB(writeDS(ds_id, "lat\0", (VOIDP) calstr.lat, (int32) * spix, 0, 0, (int32) calstr.npixs, 1, 1));

    l = 0;
    for (i = 0; i < calstr.nprods; i++) {
        if (strcmp(calstr.l2prods[i], "l2_flags\0") == 0);
        else
            if (strcmp(calstr.l2prods[i], "mside\0") == 0);
        else
            if (strcmp(calstr.l2prods[i], "detnum\0") == 0);
        else
            if (strcmp(calstr.l2prods[i], "pixnum\0") == 0);
        else {
            PTB(writeDS(ds_id, calstr.l2prods[i], (VOIDP) (calstr.ddata + l * calstr.npixs),
                    (int32) * spix, 0, 0, (int32) calstr.npixs, 1, 1));
            ++l;
        }
    }


    *spix += calstr.npixs;

    return (LIFE_IS_GOOD);

}

int read_glbl_attr(idDS ds_id, char *name, VOIDP ptr) {
    if (readAttr(ds_id, name, (VOIDP) ptr)) {
        printf("-E- %s line %d: Could not get global attribute, %s.\n", __FILE__, __LINE__, (name));
        return (HDF_FUNCTION_ERROR);
    }

    return (LIFE_IS_GOOD);
}

int32_t alloc_calstr(int32_t nfiles, int32_t npixs, mscalstr *calstr) {

    int32_t len, l;
    unsigned char *p;


    l = 5 * sizeof (int16) + sizeof (int32) + 2 * sizeof (uint8) + 2 * sizeof (float);
    len = nfiles * sizeof (stname) + npixs*l;

    for (l = 0; l < calstr->nprods; l++) {

        if (strcmp(calstr->l2prods[l], "l2_flags") == 0);
        else
            if (strcmp(calstr->l2prods[l], "mside") == 0);
        else
            if (strcmp(calstr->l2prods[l], "detnum") == 0);
        else
            if (strcmp(calstr->l2prods[l], "pixnum") == 0);
        else
            len += npixs * sizeof (float32);
    }

    if ((p = (unsigned char *) malloc(len)) == NULL) {
        printf("%s -Error: Cannot allocate memory to cross-calibration data\n", __FILE__);
        exit(FATAL_ERROR);
    }
    calstr->nfiles = nfiles;
    calstr->npixs = npixs;
    calstr->data = p;
    calstr->filenames = (stname *) p;
    p += nfiles * sizeof (stname);
    calstr->fileID = (int16 *) p;
    p += npixs * sizeof (int16);
    calstr->year = (int16 *) p;
    p += npixs * sizeof (int16);
    calstr->day = (int16 *) p;
    p += npixs * sizeof (int16);
    calstr->msec = (int32 *) p;
    p += npixs * sizeof (int32);
    calstr->iscan = (int16 *) p;
    p += npixs * sizeof (int16);
    calstr->mside = (uint8 *) p;
    p += npixs * sizeof (uint8);
    calstr->detnum = (uint8 *) p;
    p += npixs * sizeof (uint8);
    calstr->pixnum = (int16 *) p;
    p += npixs * sizeof (int16);
    calstr->lon = (float *) p;
    p += npixs * sizeof (float);
    calstr->lat = (float *) p;
    p += npixs * sizeof (float);
    calstr->ddata = (float32 *) p;

    return (len);

}

void free_calstr(mscalstr calstr, int all) {
    if (calstr.data != NULL) free(calstr.data);
    calstr.data = NULL;

    if (all) {
        if (calstr.Lambda != NULL) free(calstr.Lambda);
        calstr.Lambda = NULL;
        if (calstr.l2prods != NULL) free(calstr.l2prods);
        calstr.l2prods = NULL;
        if (calstr.input_parms != NULL) free(calstr.input_parms);
        calstr.input_parms = NULL;
    }
}


