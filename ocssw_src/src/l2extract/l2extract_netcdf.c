/*
 * subroutine to extract a netCDF L2 file
 */

#include <netcdf.h>

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>


#include <dfutils.h>
#include <genutils.h>

#include <l12_parms.h>
#include <l2_flags.h>

#include "l2extract.h"
#include "l2extract_3d_wv.h"

#define MAX_VARIABLES 512

#define NCDIE(function)                                                                                  \
    {                                                                                                    \
        int status = function;                                                                           \
        switch (status) {                                                                                \
            case NC_NOERR:                                                                               \
                break;                                                                                   \
            default:                                                                                     \
                printf("NetCDF error: file %s, line %d, %s\n", __FILE__, __LINE__, nc_strerror(status)); \
                exit(1);                                                                                 \
        }                                                                                                \
    }

void copyGlobalAttributes(int ncid_r, int ncid_w) {
    int numAtts;
    int i;
    char name[NC_MAX_NAME + 1];

    NCDIE(nc_inq_natts(ncid_r, &numAtts));
    for (i = 0; i < numAtts; i++) {
        NCDIE(nc_inq_attname(ncid_r, NC_GLOBAL, i, name));
        NCDIE(nc_copy_att(ncid_r, NC_GLOBAL, name, ncid_w, NC_GLOBAL));
    }
}

void copyVariableAttributes(int ncid_r, int varid_r, int ncid_w, int varid_w) {
    int numAtts;
    int i;
    char name[NC_MAX_NAME + 1];
    // int status;

    NCDIE(nc_inq_varnatts(ncid_r, varid_r, &numAtts));
    for (i = 0; i < numAtts; i++) {
        NCDIE(nc_inq_attname(ncid_r, varid_r, i, name));
        NCDIE(nc_copy_att(ncid_r, varid_r, name, ncid_w, varid_w));
    }
}

int checkIfInProdlist(const char *prodlist, const char *name) {
    char buf[256];

    if (prodlist[0] != 0) {
        strcpy(buf, ",");
        strcat(buf, name);
        strcat(buf, ",");
        char *cptr = strstr(prodlist, buf);
        if (cptr != NULL)
            return 1;
    } else {
        return 1;  // if prodlist is empty always return true
    }

    return 0;
}
/**
 * @brief Copies all the attributes and chunking/compression properties
 *
 * @param name - name of the variable
 * @param ncid_r - id of a nc input group
 * @param ncid_w - id of a nc output group
 * @param varid_r - id of a nc input variable
 * @param varid_w - id of a nc output variable
 * @param dimIds - ids of dimensions
 * @param chunkSize - chunk size
 * @param deflateLevel - deflation level
 * @param shuffle - shuffle flag. it stores the first byte of all of a variable's values in the chunk
 * contiguously, followed by all the second bytes, and so on. If the values are not all wildly different, this
 * can make the data more easily compressible
 * @param numAtts - number of attributs
 * @param chunkNelems  - The number of chunk slots in the raw data chunk cache hash table will be put here
 * @param chunkPreemption - The preemption will be put here. The preemtion value is between 0 and 1 inclusive
 * and indicates how much chunks that have been fully read are favored for preemption. A value of zero means
 * fully read chunks are treated no differently than other chunks (the preemption is strictly LRU) while a
 * value of one means fully read chunks are always preempted before other chunks
 * @param deflate - 	True to turn on deflation for this variable.
 * @param typeSize - type
 * @param numDims - number of dimensions
 * @param xtype - nc type of the variable.
 */
void set_attributes(const char *name, int ncid_r, int ncid_w, int varid_r, int *varid_w, int *dimIds,
                    size_t *chunkSize, int *deflateLevel, int *shuffle, int *numAtts, size_t *chunkNelems,
                    float *chunkPreemption, int *deflate, size_t *typeSize, int *numDims, nc_type *xtype) {
    // int status;
    NCDIE(nc_inq_var(ncid_r, varid_r, NULL, xtype, numDims, NULL, numAtts));
    NCDIE(nc_get_var_chunk_cache(ncid_r, varid_r, chunkSize, chunkNelems, chunkPreemption));
    NCDIE(nc_inq_var_deflate(ncid_r, varid_r, shuffle, deflate, deflateLevel));

    NCDIE(nc_def_var(ncid_w, name, *xtype, *numDims, dimIds, varid_w));
    NCDIE(nc_set_var_chunk_cache(ncid_w, *varid_w, *chunkSize, *chunkNelems, *chunkPreemption));
    NCDIE(nc_def_var_deflate(ncid_w, *varid_w, *shuffle, *deflate, *deflateLevel));

    copyVariableAttributes(ncid_r, varid_r, ncid_w, *varid_w);

    NCDIE(nc_inq_type(ncid_r, *xtype, NULL, typeSize));
}

/**
 * copy a piece of a variable and all of it's attributes to another file
 *
 * @param ncid_r netCDF file or group to read
 * @param name name of the variable
 * @param start location to start copying from
 * @param count how many of each dimension to copy
 * @param dimIds dimension IDs from the destination file to attach to the new variable
 * @param ncid_w netCDF file or group to write the variable to
 * @param data write this data to the new variable or copy from read variable if NULL
 */
void copyVariable(int ncid_r, const char *name, size_t *start, size_t *count, int *dimIds, int ncid_w,
                  void *data) {
    int status;
    int varid_r;
    int varid_w;
    nc_type xtype;
    int numDims;
    int numAtts;
    size_t chunkSize;
    size_t chunkNelems;
    float chunkPreemption;

    size_t typeSize;
    int arraySize;
    int i;

    static int localDataSize = 0;
    static char *localData = NULL;

    int shuffle;
    int deflate;
    int deflateLevel;

    status = nc_inq_varid(ncid_r, name, &varid_r);
    if (status == NC_NOERR) {
        set_attributes(name, ncid_r, ncid_w, varid_r, &varid_w, dimIds, &chunkSize, &deflateLevel, &shuffle,
                       &numAtts, &chunkNelems, &chunkPreemption, &deflate, &typeSize, &numDims, &xtype);
        if (data == NULL) {
            // calc array size in num of bytes
            arraySize = typeSize;
            for (i = 0; i < numDims; i++) {
                arraySize *= count[i];
            }

            // allocate array
            if (arraySize > localDataSize) {
                if (localData)
                    free(localData);
                localDataSize = arraySize;
                localData = (char *)malloc(localDataSize);
                if (localData == NULL) {
                    printf("-E- %s %d: could not allocate data for variable %s\n", __FILE__, __LINE__, name);
                    exit(EXIT_FAILURE);
                }
            }
            NCDIE(nc_get_vara(ncid_r, varid_r, start, count, localData));
            NCDIE(nc_put_var(ncid_w, varid_w, localData));
        } else {
            NCDIE(nc_put_var(ncid_w, varid_w, data));
        }
    }
}

/**
 * @brief
 * copy a piece of a variable  and all of it's attributes to another file. Slicing is performed along a
 * selected dimension defined by a set of indexes along the dimension
 * @param ncid_r netCDF file or group to read
 * @param name name of the variable
 * @param start location to start copying from
 * @param count how many of each dimension to copy
 * @param dimIds dimension IDs from the destination file to attach to the new
 * variable
 * @param ncid_w netCDF file or group to write the variable to
 * @param data write this data to the new variable or copy from read variable if
 * @param dim_of_indexes - index of the dimension along which slicing is performed.
 * @param indexes indexes along the dimension to be selected for output
 * NULL
 */
void copyVariableSelectedIndexes(int ncid_r, const char *name, size_t *start, size_t *count, int *dimIds,
                                 int ncid_w, void *data, int *indexes, int dim_of_indexes) {
    if (dim_of_indexes < 0) {
        printf("Supplied dim_of_indexes is negative %d\n", dim_of_indexes);
        exit(EXIT_FAILURE);
    }
    size_t start_local[3] = {0, 0, 0};
    size_t count_local[3] = {1, 1, 1};
    size_t start_global[3];
    size_t count_global[3];
    for (size_t i = 0; i < 3; i++) {
        start_global[i] = start[i];
        count_global[i] = count[i];
    }
    int status;
    int varid_r;
    int varid_w;
    nc_type xtype;
    int numDims;
    int numAtts;
    size_t chunkSize;
    size_t chunkNelems;
    float chunkPreemption;

    size_t typeSize;
    int arraySize;
    int i;

    static int localDataSize = 0;
    static char *localData = NULL;

    int shuffle;
    int deflate;
    int deflateLevel;

    status = nc_inq_varid(ncid_r, name, &varid_r);
    if (status == NC_NOERR) {
        set_attributes(name, ncid_r, ncid_w, varid_r, &varid_w, dimIds, &chunkSize, &deflateLevel, &shuffle,
                       &numAtts, &chunkNelems, &chunkPreemption, &deflate, &typeSize, &numDims, &xtype);
        if (data == NULL) {
            // calc array size in num of bytes
            arraySize = typeSize;
            if (dim_of_indexes >= numDims) {
                printf(
                    "Supplied dim_of_indexes execedes the number of dimensions "
                    "%d, %d\n",
                    dim_of_indexes, numDims);
                exit(EXIT_FAILURE);
            }
            for (i = 0; i < numDims; i++) {
                if (i == dim_of_indexes)
                    continue;
                arraySize *= count_global[i];
                count_local[i] = count_global[i];
            }

            // allocate array
            if (arraySize > localDataSize) {
                if (localData)
                    free(localData);
                localDataSize = arraySize;
                localData = (char *)malloc(localDataSize);
                if (localData == NULL) {
                    printf("-E- %s %d: could not allocate data for variable %s\n", __FILE__, __LINE__, name);
                    exit(EXIT_FAILURE);
                }
            }
            size_t total_size = count_global[dim_of_indexes];
            count_global[dim_of_indexes] = 1;
            for (int i_index = 0; i_index < total_size; i_index++) {
                start_global[dim_of_indexes] = indexes[i_index];
                start_local[dim_of_indexes] = i_index;
                NCDIE(nc_get_vara(ncid_r, varid_r, start_global, count_global, localData));
                NCDIE(nc_put_vara(ncid_w, varid_w, start_local, count_local, localData));
            }
        } else {
            NCDIE(nc_put_var(ncid_w, varid_w, data));
        }
    }
}

/**
 * extract a L2 netCDF file
 *
 * @param infile input file name
 * @param outfile output file name
 * @param spix start pixel (1 based)
 * @param epix ending pixel (1 based)
 * @param sscan start line (1 based)
 * @param escan end line (1 based)
 * @param prodlist product list, comma separated, empty string outputs all products
 * @return 0 = success
 */
int extractNetCDF(const char *infile, const char *outfile, int spix, int epix, int sscan, int escan,
                  const char *prodlist, const char *wavelist) {
    int status;
    idDS ds_id_r, ds_id_w;  // data set file ids for reading and writing
    int rootGroup_r;        // netCDF group for read file root
    int rootGroup_w;        // netCDF group for write file root
    int group_r;            // netCDF group for read file
    int group_w;            // netCDF group for write file
    int subGroup_r;         // netCDF sub group for read file
    int subGroup_w;         // netCDF sub group for write file
    int dim_id;
    int varid;
    size_t start[3] = {0, 0, 0};
    size_t count[3] = {1, 1, 1};
    int dimIds[3] = {0, 0, 0};

    size_t numLines_r;                    // number of line in read file
    size_t numPixels_r;                   // number of pixels in read file
    size_t numLines_w;                    // number of line in write file
    size_t numPixels_w;                   // number of pixels in write file
    size_t numTotalBands;                 // number of total bands (vis + IR) in file
    size_t numVisibleBands;               // number of visible bands in file
    size_t numBandsPerPixel;              // bands_per_pixel dimension
    size_t numReflectanceLocationValues;  // number_of_reflectance_location_values dimension
    size_t number_of_wv_3d = 0;           // wv 3d
    int reflectanceLocationValuesFound = 0;
    int numLinesDimId;                     // Number_of_Scan_Lines dimension id
    int numPixelsDimId;                    // Number_of_Scan_Lines dimension id
    int numBandsPerPixelDimId;             // bands_per_pixel dimension id
    int numReflectiveLocationValuesDimId;  // number_of_reflectance_locations dimension id
    int numTotalBandsDimId;                // Number_of_Scan_Lines dimension id
    int numVisibleBandsDimId;              // Number_of_Scan_Lines dimension id
    int numwv3dDimId;                      // netcdf id for number_of_wv_3d dim

    int i;
    int numVariables;
    int variableIds[MAX_VARIABLES];
    char name[256];

    // data sets read in
    float *latArray = NULL;
    float *lonArray = NULL;
    int *flagArray = NULL;

    int *wv_values_to_pass = NULL;
    int *wv_indexes_to_pass = NULL;
    int wv_num_to_pass = 0;

    get_wv3_indexes(infile, wavelist, &wv_values_to_pass, &wv_indexes_to_pass, &wv_num_to_pass, prodlist);

    ds_id_r = openDS(infile);
    if (ds_id_r.fid == -1) {
        printf("could not open \"%s\" for reading.\n", infile);
        exit(EXIT_FAILURE);
    }
    if (ds_id_r.fftype != DS_NCDF) {
        printf("could not open \"%s\" is not a netCDF4 file.\n", infile);
        exit(EXIT_FAILURE);
    }
    rootGroup_r = ds_id_r.fid;

    /* Get # of scans and # of pixels */
    /* ------------------------------ */
    status = nc_inq_dimid(rootGroup_r, "number_of_lines", &dim_id);
    if (status != NC_NOERR) {
        nc_inq_dimid(rootGroup_r, "Number_of_Scan_Lines", &dim_id);
    }
    nc_inq_dimlen(rootGroup_r, dim_id, &numLines_r);

    status = nc_inq_dimid(rootGroup_r, "pixels_per_line", &dim_id);
    if (status != NC_NOERR) {
        nc_inq_dimid(rootGroup_r, "Pixels_per_Scan_Line", &dim_id);
    }
    nc_inq_dimlen(rootGroup_r, dim_id, &numPixels_r);

    status = nc_inq_dimid(rootGroup_r, "number_of_bands", &dim_id);
    if (status != NC_NOERR) {
        nc_inq_dimid(rootGroup_r, "total_band_number", &dim_id);
    }
    nc_inq_dimlen(rootGroup_r, dim_id, &numTotalBands);

    status = nc_inq_dimid(rootGroup_r, "number_of_reflective_bands", &dim_id);
    if (status != NC_NOERR) {
        nc_inq_dimid(rootGroup_r, "band_number", &dim_id);
    }
    nc_inq_dimlen(rootGroup_r, dim_id, &numVisibleBands);

    status = nc_inq_dimid(rootGroup_r, "bands_per_pixel", &dim_id);
    if (status == NC_NOERR) {
        nc_inq_dimlen(rootGroup_r, dim_id, &numBandsPerPixel);
    } else {
        numBandsPerPixel = numTotalBands;
    }

    status = nc_inq_dimid(rootGroup_r, "wavelength_3d", &dim_id);
    if (status == NC_NOERR)
        nc_inq_dimlen(rootGroup_r, dim_id, &number_of_wv_3d);

    status = nc_inq_dimid(rootGroup_r, "number_of_reflectance_location_values", &dim_id);
    if (status == NC_NOERR) {
        nc_inq_dimlen(rootGroup_r, dim_id, &numReflectanceLocationValues);
        reflectanceLocationValuesFound = 1;
    }

    char *title = readAttrStr(ds_id_r, "title");
    if (strstr(title, "Level-2") == NULL) {
        printf("\"%s\" is not a L2 file.\n", infile);
        exit(EXIT_FAILURE);
    }
    free(title);

    if (sscan < 1) {
        sscan = 1;
    }
    if (sscan >= numLines_r) {
        printf("sscan needs to be less than number of scans in file.\n");
        exit(BOUNDS_ERROR);
    }
    if (escan < 1 || escan > numLines_r) {
        escan = numLines_r;
    }
    if (escan < sscan) {
        printf("escan needs to be greater than sscan.\n");
        exit(BOUNDS_ERROR);
    }
    numLines_w = escan - sscan + 1;

    if (spix < 1) {
        spix = 1;
    }
    if (spix >= numPixels_r) {
        printf("spix needs to be less than number of pixels in file.\n");
        exit(BOUNDS_ERROR);
    }
    if (epix < 1 || epix > numPixels_r) {
        epix = numPixels_r;
    }
    if (epix < spix) {
        printf("epix needs to be greater than spix.\n");
        exit(BOUNDS_ERROR);
    }

    // --------------------------------------------------------
    // make sure all products requested exist in the file
    // --------------------------------------------------------
    if (prodlist[0] != 0) {
        status = nc_inq_ncid(rootGroup_r, "geophysical_data", &group_r);
        if (status != NC_NOERR) {
            printf("-E- Could not open \"geophysical_data\" group.\n");
            exit(EXIT_FAILURE);
        }
        char *word;
        char *tmpProdList = strdup(prodlist);
        word = strtok(tmpProdList, ",");
        while (word) {
            if (word[0] != 0) {
                status = nc_inq_varid(group_r, word, &varid);
                if (status != NC_NOERR) {
                    printf("-E- Could not find product \"%s\" in L2 file.\n", word);
                    exit(EXIT_FAILURE);
                }
            }
            word = strtok(NULL, ",");
        }
        free(tmpProdList);
    }

    // --------------------------------------------------------
    // set to navigation data group
    // --------------------------------------------------------
    nc_inq_ncid(rootGroup_r, "navigation_data", &group_r);
   
    printf("sscan: %d  escan: %d\n", sscan, escan);
    printf("spixl: %d  epixl: %d\n", spix, epix);

    numPixels_w = epix - spix + 1;

    // read lat and lon
    latArray = malloc(numLines_w * numPixels_w * sizeof(float));
    lonArray = malloc(numLines_w * numPixels_w * sizeof(float));
    if (latArray == NULL || lonArray == NULL) {
        printf("could not allocate memory for lat and lon arrays\n");
        exit(EXIT_FAILURE);
    }
    start[0] = sscan - 1;
    start[1] = spix - 1;
    count[0] = numLines_w;
    count[1] = numPixels_w;
    NCDIE(nc_inq_varid(group_r, "latitude", &varid));
    NCDIE(nc_get_vara_float(group_r, varid, start, count, latArray));
    NCDIE(nc_inq_varid(group_r, "longitude", &varid));
    NCDIE(nc_get_vara_float(group_r, varid, start, count, lonArray));

    // --------------------------------------------------------
    // set to geophysical data group
    // --------------------------------------------------------
    nc_inq_ncid(rootGroup_r, "geophysical_data", &group_r);

    flagArray = malloc(numLines_w * numPixels_w * sizeof(int));
    if (flagArray == NULL) {
        printf("could not allocate memory for flag array\n");
        exit(EXIT_FAILURE);
    }
    start[0] = sscan - 1;
    start[1] = spix - 1;
    count[0] = numLines_w;
    count[1] = numPixels_w;
    NCDIE(nc_inq_varid(group_r, "l2_flags", &varid));
    NCDIE(nc_get_vara_int(group_r, varid, start, count, flagArray));

    /* Calculate new navigation metadata */
    /* --------------------------------- */
#define GEOBOX_INC 20.0

    int startDone = 0;
    int centerDone = 0;
    int endDone = 0;
    int ccol = numPixels_w / 2;
    int line, scol, ecol;
    int lineStart;
    float last_lat = 0.0;
    int lastGoodLine = 0;
    float geobox[4][100];
    int32_t geobox_cnt = 0;
    float gring_fval[100];
    int32_t gring_ival[100];

    float startCenterLat, startCenterLon;
    float endCenterLat, endCenterLon;
    float geoLatMin, geoLatMax;
    float geoLonMin, geoLonMax;

    // loop to find beginning info
    for (line = 0; line < numLines_w; line++) {
        lineStart = numPixels_w * line;

        // find good start pixel
        if (!startDone) {
            for (scol = 0; scol <= ccol; scol++) {
                // check NAVFAIL flag
                if (flagArray[line * numPixels_w + scol] ^ NAVFAIL) {
                    startDone = 1;
                    geobox[0][geobox_cnt] = lonArray[lineStart + scol];
                    geobox[1][geobox_cnt] = latArray[lineStart + scol];
                    break;
                }  // flag good
            }      // for col
        }          // if not start

        // find good center pixel
        if (!centerDone) {
            // check NAVFAIL flag
            if (flagArray[line * numPixels_w + ccol] ^ NAVFAIL) {
                centerDone = 1;
                startCenterLon = lonArray[lineStart + ccol];
                startCenterLat = latArray[lineStart + ccol];
                last_lat = startCenterLat;
            }  // flag good
        }      // if not center

        // find good end pixel
        if (!endDone) {
            for (ecol = numPixels_w - 1; ecol >= ccol; ecol--) {
                // check NAVFAIL flag
                if (flagArray[line * numPixels_w + ecol] ^ NAVFAIL) {
                    endDone = 1;
                    geobox[2][geobox_cnt] = lonArray[lineStart + ecol];
                    geobox[3][geobox_cnt] = latArray[lineStart + ecol];
                    break;
                }  // flag good
            }      // for col
        }          // if not start

        if (startDone && centerDone && endDone)
            break;
    }

    // set the min and max lat lon values
    geoLonMin = geoLonMax = geobox[0][geobox_cnt];
    geoLatMin = geoLatMax = geobox[1][geobox_cnt];
    if (geoLonMin > geobox[2][geobox_cnt])
        geoLonMin = geobox[2][geobox_cnt];
    if (geoLatMin > geobox[3][geobox_cnt])
        geoLatMin = geobox[3][geobox_cnt];
    if (geoLonMax < geobox[2][geobox_cnt])
        geoLonMax = geobox[2][geobox_cnt];
    if (geoLatMax < geobox[3][geobox_cnt])
        geoLatMax = geobox[3][geobox_cnt];
    if (geoLonMin > startCenterLon)
        geoLonMin = startCenterLon;
    if (geoLonMax < startCenterLon)
        geoLonMax = startCenterLon;

    geobox_cnt++;

    // loop through the rest of the lines
    for (; line < numLines_w; line++) {
        lineStart = numPixels_w * line;

        // find first good pixel on line
        for (scol = 0; scol <= ccol; scol++) {
            // check NAVFAIL flag
            if (flagArray[line * numPixels_w + scol] ^ NAVFAIL)
                break;
        }
        if (scol == ccol)  // could not find start col, so skip this line
            continue;

        // find last good pixel
        for (ecol = numPixels_w - 1; ecol >= ccol; ecol--) {
            // check NAVFAIL flag
            if (flagArray[line * numPixels_w + ecol] ^ NAVFAIL)
                break;
        }
        if (ecol < ccol)  // could not find end col, so skip this line
            continue;

        lastGoodLine = line;

        // set min/max for every line
        if (geoLonMax < lonArray[lineStart + scol])
            geoLonMax = lonArray[lineStart + scol];
        if (geoLonMax < lonArray[lineStart + ecol])
            geoLonMax = lonArray[lineStart + ecol];
        if (geoLonMin > lonArray[lineStart + scol])
            geoLonMin = lonArray[lineStart + scol];
        if (geoLonMin > lonArray[lineStart + ecol])
            geoLonMin = lonArray[lineStart + ecol];

        if (geoLatMax < latArray[lineStart + scol])
            geoLatMax = latArray[lineStart + scol];
        if (geoLatMax < latArray[lineStart + ecol])
            geoLatMax = latArray[lineStart + ecol];
        if (geoLatMin > latArray[lineStart + scol])
            geoLatMin = latArray[lineStart + scol];
        if (geoLatMin > latArray[lineStart + ecol])
            geoLatMin = latArray[lineStart + ecol];

        // load up geobox
        if (fabs(last_lat - latArray[lineStart + ccol]) > GEOBOX_INC) {
            geobox[0][geobox_cnt] = lonArray[lineStart + scol];
            geobox[1][geobox_cnt] = latArray[lineStart + scol];
            geobox[2][geobox_cnt] = lonArray[lineStart + ecol];
            geobox[3][geobox_cnt] = latArray[lineStart + ecol];
            last_lat = latArray[lineStart + ccol];
            geobox_cnt++;
        }

    }  // for lines

    // make sure we add the last line
    lineStart = numPixels_w * lastGoodLine;

    // find first good pixel on line
    for (scol = 0; scol < ccol; scol++) {
        // check NAVFAIL flag
        if (flagArray[lastGoodLine * numPixels_w + scol] ^ NAVFAIL)
            break;
    }

    // find last good pixel
    for (ecol = numPixels_w - 1; ecol >= ccol; ecol--) {
        // check NAVFAIL flag
        if (flagArray[lastGoodLine * numPixels_w + ecol] ^ NAVFAIL)
            break;
    }

    endCenterLon = lonArray[lineStart + ccol];
    endCenterLat = latArray[lineStart + ccol];

    geobox[0][geobox_cnt] = lonArray[lineStart + scol];
    geobox[1][geobox_cnt] = latArray[lineStart + scol];
    geobox[2][geobox_cnt] = lonArray[lineStart + ecol];
    geobox[3][geobox_cnt] = latArray[lineStart + ecol];
    geobox_cnt++;

    /* Create output netCDF file (delete if it exists) */
    /* ------------------------- */
    if (access(outfile, F_OK) != -1) {
        if (unlink(outfile)) {
            printf("could not delete existing output file %s\n", outfile);
            exit(EXIT_FAILURE);
        }
    }
    ds_id_w = startDS(outfile, DS_NCDF, DS_WRITE, 5);
    rootGroup_w = ds_id_w.fid;

    /* Create dimensions */
    /* ----------------- */
    NCDIE(nc_def_dim(rootGroup_w, "number_of_lines", numLines_w, &numLinesDimId));
    NCDIE(nc_def_dim(rootGroup_w, "pixels_per_line", numPixels_w, &numPixelsDimId));
    NCDIE(nc_def_dim(rootGroup_w, "bands_per_pixel", numBandsPerPixel, &numBandsPerPixelDimId));
    if (number_of_wv_3d > 0) {
        if (wv_num_to_pass > 0) {
            number_of_wv_3d = wv_num_to_pass;
        }
        NCDIE(nc_def_dim(rootGroup_w, "wavelength_3d", number_of_wv_3d, &numwv3dDimId));
    }

    if (reflectanceLocationValuesFound) {
        NCDIE(nc_def_dim(rootGroup_w, "number_of_reflectance_location_values", numReflectanceLocationValues,
                         &numReflectiveLocationValuesDimId));
    }
    NCDIE(nc_def_dim(rootGroup_w, "number_of_bands", numTotalBands, &numTotalBandsDimId));
    NCDIE(nc_def_dim(rootGroup_w, "number_of_reflective_bands", numVisibleBands, &numVisibleBandsDimId));

    // --------------------------------------------------------
    // set to sensor_band_parameters group
    // --------------------------------------------------------
    nc_inq_ncid(rootGroup_r, "sensor_band_parameters", &group_r);
    nc_def_grp(rootGroup_w, "sensor_band_parameters", &group_w);

    // copy vcal_gain(visibleBands)
    start[0] = 0;
    count[0] = numVisibleBands;
    dimIds[0] = numVisibleBandsDimId;

    // loop through all of the variables
    status = nc_inq_varids(group_r, &numVariables, variableIds);
    if (status == NC_NOERR) {
        if (numVariables > MAX_VARIABLES) {
            printf("-E- %s %d: too many variables in group \"sensor_band_parameters\"\n",
                __FILE__, __LINE__);
            exit(EXIT_SUCCESS);
        }
        for (i = 0; i < numVariables; i++) {
            status = nc_inq_varname(group_r, variableIds[i], name);
            if ((strcmp(name, "wavelength") == 0)) {
                count[0] = numTotalBands;
                dimIds[0] = numTotalBandsDimId;
                copyVariable(group_r, name, start, count, dimIds, group_w, NULL);
            } else {
                if ((strcmp(name, "wavelength_3d") == 0)) {
                    count[0] = number_of_wv_3d;
                    dimIds[0] = numwv3dDimId;
                    if (wv_num_to_pass > 0) {
                        copyVariableSelectedIndexes(group_r, name, start, count, dimIds, group_w, NULL,
                                                    wv_indexes_to_pass, 0);
                        continue;
                    }
                } else {
                    count[0] = numVisibleBands;
                    dimIds[0] = numVisibleBandsDimId;
                }
                copyVariable(group_r, name, start, count, dimIds, group_w, NULL);
            }
        }
    }

    // --------------------------------------------------------
    // set to sensor_band_parameters group
    // --------------------------------------------------------
    nc_inq_ncid(rootGroup_r, "scan_line_attributes", &group_r);
    nc_def_grp(rootGroup_w, "scan_line_attributes", &group_w);

    // set up the copy parameters
    start[0] = sscan - 1;
    count[0] = numLines_w;
    dimIds[0] = numLinesDimId;

    // loop through all of the variables
    status = nc_inq_varids(group_r, &numVariables, variableIds);
    if (status == NC_NOERR) {
        if (numVariables > MAX_VARIABLES) {
            printf("-E- %s %d: too many variables in group \"scan_line_attributes\"\n",
                __FILE__, __LINE__);
            exit(EXIT_SUCCESS);
        }
        for (i = 0; i < numVariables; i++) {
            status = nc_inq_varname(group_r, variableIds[i], name);
            copyVariable(group_r, name, start, count, dimIds, group_w, NULL);
        }
    }

    // --------------------------------------------------------
    // set to geophysical_data group
    // --------------------------------------------------------
    nc_inq_ncid(rootGroup_r, "geophysical_data", &group_r);
    nc_def_grp(rootGroup_w, "geophysical_data", &group_w);

    // set up the copy parameters
    start[0] = sscan - 1;
    start[1] = spix - 1;
    count[0] = numLines_w;
    count[1] = numPixels_w;
    dimIds[0] = numLinesDimId;
    dimIds[1] = numPixelsDimId;

    // loop through all of the variables
    status = nc_inq_varids(group_r, &numVariables, variableIds);
    if (status == NC_NOERR) {
        if (numVariables > MAX_VARIABLES) {
            printf("-E- %s %d: too many variables in group \"geophysical_data\"\n",
                __FILE__, __LINE__);
            exit(EXIT_SUCCESS);
        }
        for (i = 0; i < numVariables; i++) {
            status = nc_inq_varname(group_r, variableIds[i], name);
            {
                int varids_dims[3] = {0, 0, 0};
                char dim_name[50];
                nc_inq_vardimid(group_r, variableIds[i], varids_dims);
                nc_inq_dimname(group_r, varids_dims[2], dim_name);
                if (strcmp(dim_name, "wavelength_3d") == 0) {
                    dimIds[2] = numwv3dDimId;
                    count[2] = number_of_wv_3d;
                    if (wv_indexes_to_pass > 0) {
                        if (checkIfInProdlist(prodlist, name) == 0)
                            continue;
                        copyVariableSelectedIndexes(group_r, name, start, count, dimIds, group_w, NULL,
                                                    wv_indexes_to_pass, 2);
                        continue;
                    }
                } else {
                    dimIds[2] = 0;
                    count[2] = 1;
                }
            }
            if (checkIfInProdlist(prodlist, name) == 0)
                continue;
            copyVariable(group_r, name, start, count, dimIds, group_w, NULL);
        }
    }

    // --------------------------------------------------------
    // set to navigation data group
    // --------------------------------------------------------
    nc_inq_ncid(rootGroup_r, "navigation_data", &group_r);
    nc_def_grp(rootGroup_w, "navigation_data", &group_w);

    start[0] = sscan - 1;
    start[1] = spix - 1;
    count[0] = numLines_w;
    count[1] = numPixels_w;
    dimIds[0] = numLinesDimId;
    dimIds[1] = numPixelsDimId;
    copyVariable(group_r, "longitude", start, count, dimIds, group_w, lonArray);
    copyVariable(group_r, "latitude", start, count, dimIds, group_w, latArray);

    start[0] = sscan - 1;
    count[0] = numLines_w;
    dimIds[0] = numLinesDimId;
    copyVariable(group_r, "tilt", start, count, dimIds, group_w, NULL);

    // fill up gring
    int j = 1;
    gring_fval[0] = geobox[0][0];
    for (i = 0; i < geobox_cnt; i++) {
        gring_fval[j++] = geobox[2][i];
    }
    for (i = 0; i < geobox_cnt - 1; i++) {
        gring_fval[j++] = geobox[0][geobox_cnt - 1 - i];
    }
    NCDIE(nc_put_att_float(group_w, NC_GLOBAL, "gringpointlongitude", NC_FLOAT, j, gring_fval));

    j = 1;
    gring_fval[0] = geobox[1][0];
    gring_ival[0] = j;
    for (i = 0; i < geobox_cnt; i++) {
        gring_ival[j] = j + 1;
        gring_fval[j++] = geobox[3][i];
    }
    for (i = 0; i < geobox_cnt - 1; i++) {
        gring_ival[j] = j + 1;
        gring_fval[j++] = geobox[1][geobox_cnt - 1 - i];
    }
    NCDIE(nc_put_att_float(group_w, NC_GLOBAL, "gringpointlatitude", NC_FLOAT, j, gring_fval));
    NCDIE(nc_put_att_int(group_w, NC_GLOBAL, "gringpointsequence", NC_INT, j, gring_ival));

    // --------------------------------------------------------
    // set to processing_control group
    // --------------------------------------------------------
    nc_inq_ncid(rootGroup_r, "processing_control", &group_r);
    nc_def_grp(rootGroup_w, "processing_control", &group_w);

    copyGlobalAttributes(group_r, group_w);

    // sub group input_parameters
    nc_inq_ncid(group_r, "input_parameters", &subGroup_r);
    nc_def_grp(group_w, "input_parameters", &subGroup_w);

    copyGlobalAttributes(subGroup_r, subGroup_w);

    // sub group flag_percentages
    nc_inq_ncid(group_r, "flag_percentages", &subGroup_r);
    nc_def_grp(group_w, "flag_percentages", &subGroup_w);

    copyGlobalAttributes(subGroup_r, subGroup_w);

    // --------------------------------------------------------
    // copy global attributes
    // --------------------------------------------------------
    copyGlobalAttributes(rootGroup_r, rootGroup_w);

    // write modified global attrbutes
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "start_center_longitude", NC_FLOAT, 1, &startCenterLon));
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "start_center_latitude", NC_FLOAT, 1, &startCenterLat));
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "end_center_longitude", NC_FLOAT, 1, &endCenterLon));
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "end_center_latitude", NC_FLOAT, 1, &endCenterLat));

    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "northernmost_latitude", NC_FLOAT, 1, &geoLatMax));
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "southernmost_latitude", NC_FLOAT, 1, &geoLatMin));
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "easternmost_longitude", NC_FLOAT, 1, &geoLonMax));
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "westernmost_longitude", NC_FLOAT, 1, &geoLonMin));

    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "geospatial_lat_max", NC_FLOAT, 1, &geoLatMax));
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "geospatial_lat_min", NC_FLOAT, 1, &geoLatMin));
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "geospatial_lon_max", NC_FLOAT, 1, &geoLonMax));
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "geospatial_lon_min", NC_FLOAT, 1, &geoLonMin));
    nc_close(rootGroup_r);
    nc_close(rootGroup_w);

    free(latArray);
    free(lonArray);
    free(flagArray);

    return EXIT_SUCCESS;
}
