/* =========================================================== */
/* Module ncdf_utils.c                                         */
/*                                                             */
/* NCDF4 I/O utilities.                                        */
/*                                                             */
/* Written By:                                                 */
/*     Joel Gales, Futurtech                                   */
/*                                                             */
/* Modification History:                                       */
/*     Joel Gales, Futuretech, OBPG Project, Nov 2013.         */
/*           Add support for CF-compliant metadata             */
/* =========================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h> 
#include <time.h>
#include <math.h>
#include <netcdf.h>
#include <nc4utils.h>
#include <ncdfbin_utils.h>

#define MAX_PRODS 2048
#define FILE_IS_32 0
#define FILE_IS_64 1


// had to define some evil file level global variables
static nc_type binListType = NC_NAT;
static int binListDim = -1;
static int binListVarid = -1;
static int binDataDim = -1;
static int numberOfProducts = 0;
static int binDataVarid[MAX_PRODS];
static int binIndexDim = -1;
static int binIndexVarid = -1;
static int qualityDim = -1;
static int qualityVarid = -1;
static int sizeOfFile = -1;

int defineBinList_nc(int32_t deflate, int32_t grpid) {
    // Define BinList DataSet
    int status;

    // make sure we do not bounce between 32 and 64 bit bin file
    if(sizeOfFile == FILE_IS_64) {
        (void) fprintf(stderr, "line %d of %s: defining 32bit BinList when already in 64bit mode\n", __LINE__, __FILE__);
        exit(1);
    }
    sizeOfFile = FILE_IS_32;
    
    if (nc_inq_varid(grpid, "BinList", &binListVarid) == NC_NOERR) {
        (void) fprintf(stderr, "line %d of %s: BinList is already defined\n", __LINE__, __FILE__);
        exit(1);
    }

    status = nc_def_compound(grpid, sizeof (binListStruct_nc), "binListType", &binListType);
    check_err(status, __LINE__, __FILE__);

    status = nc_insert_compound(grpid, binListType, "bin_num",
            NC_COMPOUND_OFFSET(binListStruct_nc, binnum),
            NC_UINT);
    check_err(status, __LINE__, __FILE__);

    status = nc_insert_compound(grpid, binListType, "nobs",
            NC_COMPOUND_OFFSET(binListStruct_nc, nobs),
            NC_SHORT);
    check_err(status, __LINE__, __FILE__);

    status = nc_insert_compound(grpid, binListType, "nscenes",
            NC_COMPOUND_OFFSET(binListStruct_nc, nscenes),
            NC_SHORT);
    check_err(status, __LINE__, __FILE__);

    status = nc_insert_compound(grpid, binListType, "weights",
            NC_COMPOUND_OFFSET(binListStruct_nc, weights),
            NC_FLOAT);
    check_err(status, __LINE__, __FILE__);

    status = nc_insert_compound(grpid, binListType, "time_rec",
            NC_COMPOUND_OFFSET(binListStruct_nc, time_rec),
            NC_FLOAT);
    check_err(status, __LINE__, __FILE__);

    status = nc_def_dim(grpid, "binListDim", NC_UNLIMITED, &binListDim);
    check_err(status, __LINE__, __FILE__);

    status = nc_def_var(grpid, "BinList", binListType, 1, &binListDim, &binListVarid);
    check_err(status, __LINE__, __FILE__);

    /* First set chunking */
    size_t chunksize = 256;
    nc_init_compress(grpid, binListVarid, &binListDim, 1, &chunksize, deflate);

    return 0;
}

int defineBinList64_nc(int32_t deflate, int32_t grpid) {
    // Define BinList DataSet
    int status;

    // make sure we do not bounce between 32 and 64 bit bin file
    if(sizeOfFile == FILE_IS_32) {
        (void) fprintf(stderr, "line %d of %s: defining 64bit BinList when already in 32bit mode\n", __LINE__, __FILE__);
        exit(1);
    }
    sizeOfFile = FILE_IS_64;
    
    if (nc_inq_varid(grpid, "BinList", &binListVarid) == NC_NOERR) {
        (void) fprintf(stderr, "line %d of %s: BinList is already defined\n", __LINE__, __FILE__);
        exit(1);
    }

    //Turns out the size of the structure is 24 even though it's 20 if you count up the bytes
    //size_t offset = sizeof (binListStruct64_nc);      
    status = nc_def_compound(grpid, sizeof (binListStruct64_nc), "binListType", &binListType);
    check_err(status, __LINE__, __FILE__);

    status = nc_insert_compound(grpid, binListType, "bin_num",
            NC_COMPOUND_OFFSET(binListStruct64_nc, binnum),
            NC_UINT64);
    check_err(status, __LINE__, __FILE__);

    status = nc_insert_compound(grpid, binListType, "nobs",
            NC_COMPOUND_OFFSET(binListStruct64_nc, nobs),
            NC_SHORT);
    check_err(status, __LINE__, __FILE__);

    status = nc_insert_compound(grpid, binListType, "nscenes",
            NC_COMPOUND_OFFSET(binListStruct64_nc, nscenes),
            NC_SHORT);
    check_err(status, __LINE__, __FILE__);

    status = nc_insert_compound(grpid, binListType, "weights",
            NC_COMPOUND_OFFSET(binListStruct64_nc, weights),
            NC_FLOAT);
    check_err(status, __LINE__, __FILE__);

    status = nc_insert_compound(grpid, binListType, "time_rec",
            NC_COMPOUND_OFFSET(binListStruct64_nc, time_rec),
            NC_FLOAT);
    check_err(status, __LINE__, __FILE__);

    status = nc_def_dim(grpid, "binListDim", NC_UNLIMITED, &binListDim);
    check_err(status, __LINE__, __FILE__);

    status = nc_def_var(grpid, "BinList", binListType, 1, &binListDim, &binListVarid);
    check_err(status, __LINE__, __FILE__);

    /* First set chunking */
    size_t chunksize = 256;
    nc_init_compress(grpid, binListVarid, &binListDim, 1, &chunksize, deflate);

    return 0;
}

int writeBinList_nc(int32_t grpid, int32_t nbins_to_write, const void *data) {
    // Write BinList DataSet
    int status;

    static size_t startp = 0;
    size_t countp;

    if (binListVarid == -1) {
        (void) fprintf(stderr, "line %d of %s: BinList variable needs to be defined first\n",
                __LINE__, __FILE__);
        exit(1);
    }

    countp = nbins_to_write;
    status = nc_put_vara(grpid, binListVarid, &startp, &countp, data);
    check_err(status, __LINE__, __FILE__);
    startp += countp;

    return 0;
}

int defineBinData_nc(int32_t deflate, int32_t grpid, int32_t nprods, char** prodnames) {
    int status = 0;

    int varid;
    nc_type binDataType = -1;
    
    status = nc_inq_typeid(grpid, "binDataType", &binDataType);

    if (status != NC_NOERR) {
      status = nc_def_compound(grpid, 8, "binDataType", &binDataType);
      check_err(status, __LINE__, __FILE__);

      status = nc_insert_compound(grpid, binDataType, "sum", 0, NC_FLOAT);
      check_err(status, __LINE__, __FILE__);

      status = nc_insert_compound(grpid, binDataType, "sum_squared", 4, NC_FLOAT);
      check_err(status, __LINE__, __FILE__);

      status = nc_def_dim(grpid, "binDataDim", NC_UNLIMITED, &binDataDim);
      check_err(status, __LINE__, __FILE__);
    }
    
    int prod;
    for (prod = 0; prod < nprods; prod++) {
        if (nc_inq_varid(grpid, prodnames[prod], &varid) == NC_NOERR) {
            fprintf(stderr, "line %d of %s: BinData for %s is already defined\n",
                    __LINE__, __FILE__, prodnames[prod]);
            exit(1);
        }

        status = nc_def_var(grpid, prodnames[prod], binDataType, 1, &binDataDim, &varid);
        if (status != NC_NOERR) {
            report_err(status, __LINE__, __FILE__);
            fprintf(stderr, "trying to create binData for product %s\n", prodnames[prod]);
            exit(1);
        }
        binDataVarid[numberOfProducts++] = varid;
        if (numberOfProducts > MAX_PRODS) {
                fprintf(stderr, "line %d of %s: Max number of output products exceeded\n",
                        __LINE__, __FILE__);
                exit(1);
        }

        /* First set chunking */
        size_t chunksize = 256;
        nc_init_compress(grpid, varid, &binDataDim, 1, &chunksize, deflate);

    }

    return 0;
}

int writeBinData_nc(int32_t grpid, int32_t nbins_to_write, int32_t iprod, const void *data) {
    int status = 0;

    static size_t startp[MAX_PRODS];
    size_t countp;

    if (iprod >= numberOfProducts) {
        fprintf(stderr, "line %d of %s: product index %d out of range\n",
                __LINE__, __FILE__, iprod);
        exit(1);
    }

    countp = nbins_to_write;
    status = nc_put_vara(grpid, binDataVarid[iprod], &startp[iprod], &countp, data);
    check_err(status, __LINE__, __FILE__);
    startp[iprod] += countp;

    return 0;
}

int defineBinIndex_nc(int32_t deflate, int32_t grpid) {
    // Define BinIndex DataSet
    int status;

    nc_type binIndexType;
    int varid;

    // make sure we do not bounce between 32 and 64 bit bin file
    if(sizeOfFile == FILE_IS_64) {
        (void) fprintf(stderr, "line %d of %s: defining 32bit BinIndex when already in 64bit mode\n", __LINE__, __FILE__);
        exit(1);
    }
    sizeOfFile = FILE_IS_32;
    
    if (nc_inq_varid(grpid, "BinIndex", &varid) == NC_NOERR) {
        fprintf(stderr, "line %d of %s: BinIndex is already defined\n",
                __LINE__, __FILE__);
        exit(1);
    }

    status = nc_def_compound(grpid, sizeof(binIndexStruct_nc), "binIndexType", &binIndexType);
    check_err(status, __LINE__, __FILE__);

    status = nc_insert_compound(grpid, binIndexType, "start_num",
            NC_COMPOUND_OFFSET(binIndexStruct_nc, start_num), NC_UINT);
    check_err(status, __LINE__, __FILE__);

    status = nc_insert_compound(grpid, binIndexType, "begin",
            NC_COMPOUND_OFFSET(binIndexStruct_nc, begin), NC_UINT);
    check_err(status, __LINE__, __FILE__);

    status = nc_insert_compound(grpid, binIndexType, "extent",
            NC_COMPOUND_OFFSET(binIndexStruct_nc, extent), NC_UINT);
    check_err(status, __LINE__, __FILE__);

    status = nc_insert_compound(grpid, binIndexType, "max",
            NC_COMPOUND_OFFSET(binIndexStruct_nc, max), NC_UINT);
    check_err(status, __LINE__, __FILE__);
    
    status = nc_def_dim(grpid, "binIndexDim", NC_UNLIMITED, &binIndexDim);
    check_err(status, __LINE__, __FILE__);

    status = nc_def_var(grpid, "BinIndex", binIndexType, 1, &binIndexDim,
            &binIndexVarid);
    check_err(status, __LINE__, __FILE__);

    /* First set chunking */
    size_t chunksize = 256;
    nc_init_compress(grpid, binIndexVarid, &binIndexDim, 1, &chunksize, deflate);

    return 0;
}

int defineBinIndex64_nc(int32_t deflate, int32_t grpid) {
    // Define BinIndex DataSet
    int status;

    nc_type binIndexType;
    int varid;

    // make sure we do not bounce between 32 and 64 bit bin file
    if(sizeOfFile == FILE_IS_32) {
        (void) fprintf(stderr, "line %d of %s: defining 64bit BinIndex when already in 32bit mode\n", __LINE__, __FILE__);
        exit(1);
    }
    sizeOfFile = FILE_IS_64;
    
    if (nc_inq_varid(grpid, "BinIndex", &varid) == NC_NOERR) {
        fprintf(stderr, "line %d of %s: BinIndex is already defined\n",
                __LINE__, __FILE__);
        exit(1);
    }

    status = nc_def_compound(grpid, sizeof(binIndexStruct64_nc), "binIndexType", &binIndexType);
    check_err(status, __LINE__, __FILE__);

    status = nc_insert_compound(grpid, binIndexType, "start_num",
            NC_COMPOUND_OFFSET(binIndexStruct64_nc, start_num), NC_UINT64);
    check_err(status, __LINE__, __FILE__);

    status = nc_insert_compound(grpid, binIndexType, "begin",
            NC_COMPOUND_OFFSET(binIndexStruct64_nc, begin), NC_UINT64);
    check_err(status, __LINE__, __FILE__);

    status = nc_insert_compound(grpid, binIndexType, "extent",
            NC_COMPOUND_OFFSET(binIndexStruct64_nc, extent), NC_UINT);
    check_err(status, __LINE__, __FILE__);

    status = nc_insert_compound(grpid, binIndexType, "max",
            NC_COMPOUND_OFFSET(binIndexStruct64_nc, max), NC_UINT64);
    check_err(status, __LINE__, __FILE__);

    status = nc_def_dim(grpid, "binIndexDim", NC_UNLIMITED, &binIndexDim);
    check_err(status, __LINE__, __FILE__);

    status = nc_def_var(grpid, "BinIndex", binIndexType, 1, &binIndexDim,
            &binIndexVarid);
    check_err(status, __LINE__, __FILE__);

    /* First set chunking */
    size_t chunksize = 256;
    nc_init_compress(grpid, binIndexVarid, &binIndexDim, 1, &chunksize, deflate);


    return 0;
}

int writeBinIndex_nc(int32_t grpid, int32_t n_write, const void *data) {
    // Write BinIndex DataSet
    int status;

    static size_t start = 0;
    size_t count = n_write;

    status = nc_put_vara(grpid, binIndexVarid, &start, &count, data);
    check_err(status, __LINE__, __FILE__);
    start += count;

    return 0;
}

int defineQuality_nc(int32_t deflate, int32_t grpid) {
    // Define Quality DataSet
    int status;
    int varid;

    if (nc_inq_varid(grpid, "qual_l3", &varid) == NC_NOERR) {
        fprintf(stderr, "line %d of %s: BinIndex is already defined\n",
                __LINE__, __FILE__);
        exit(1);
    }

    status = nc_def_dim(grpid, "qualityDim", NC_UNLIMITED, &qualityDim);
    check_err(status, __LINE__, __FILE__);

    status = nc_def_var(grpid, "qual_l3", NC_BYTE, 1, &qualityDim, &qualityVarid);
    check_err(status, __LINE__, __FILE__);

    /* First set chunking */
    size_t chunksize = 256;
    nc_init_compress(grpid, qualityVarid, &qualityDim, 1, &chunksize, deflate);

    return 0;
}

int writeQuality_nc(int32_t grpid, int32_t nbins_to_write, const void *data) {
    // Write Quality DataSet
    int status;
    static size_t startp;
    size_t countp;

    countp = nbins_to_write;
    status = nc_put_vara(grpid, qualityVarid, &startp, &countp, data);
    check_err(status, __LINE__, __FILE__);
    startp += countp;

    return 0;
}

