/*
 * subroutine to extract a netCDF L2 file
 */

#include <netcdf.h>

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

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

void copyGlobalAttributes(int groupIdRead, int groupIdWrite) {
    int nAttributes;
    int i;
    char name[NC_MAX_NAME + 1];

    NCDIE(nc_inq_natts(groupIdRead, &nAttributes));
    for (i = 0; i < nAttributes; i++) {
        NCDIE(nc_inq_attname(groupIdRead, NC_GLOBAL, i, name));
        NCDIE(nc_copy_att(groupIdRead, NC_GLOBAL, name, groupIdWrite, NC_GLOBAL));
    }
}

void copyVariableAttributes(int groupIdRead, int variableRead, int groupIdWrite, int variableWrite) {
    int nAttributes;
    int i;
    char name[NC_MAX_NAME + 1];

    NCDIE(nc_inq_varnatts(groupIdRead, variableRead, &nAttributes));
    for (i = 0; i < nAttributes; i++) {
        NCDIE(nc_inq_attname(groupIdRead, variableRead, i, name));
        NCDIE(nc_copy_att(groupIdRead, variableRead, name, groupIdWrite, variableWrite));
    }
}

int checkIfInProdlist(const char *productList, const char *name) {
    char buf[256];

    if (productList[0] != 0) {
        strcpy(buf, ",");
        strcat(buf, name);
        strcat(buf, ",");
        char *cptr = strstr(productList, buf);
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
 * @param groupIdRead - id of a nc input group
 * @param groupIdWrite - id of a nc output group
 * @param variableRead - id of a nc input variable
 * @param variableWrite - id of a nc output variable
 * @param dimIds - ids of dimensions *
 * @param count - dimension of the extraction
 * @param chunkSizes - array of chunk sizes
 * @param deflateLevel - deflation level
 * @param shuffle - shuffle flag. it stores the first byte of all of a variable's values in the chunk
 * contiguously, followed by all the second bytes, and so on. If the values are not all wildly different, this
 * can make the data more easily compressible
 * @param nAttributes - number of attributs
 * @param deflate - 	True to turn on deflation for this variable.
 * @param typeSize - type
 * @param numDims - number of dimensions
 * @param xtype - nc type of the variable.
 */
void setAttributes(const char *name, int groupIdRead, int groupIdWrite, int variableRead, int *variableWrite,
                   int *dimIds, size_t *count, size_t *chunkSizes, int *deflateLevel, int *shuffle,
                   int *nAttributes, int *deflate, size_t *typeSize, int *numDims, nc_type *xtype) {
    int i;
    int chunkingStorage;

    NCDIE(nc_inq_var(groupIdRead, variableRead, NULL, xtype, numDims, NULL, nAttributes));
    NCDIE(nc_inq_var_chunking(groupIdRead, variableRead, &chunkingStorage, chunkSizes));
    NCDIE(nc_inq_var_deflate(groupIdRead, variableRead, shuffle, deflate, deflateLevel));

    for (i = 0; i < 3; i++) {
        if (count[i] != 1 && count[i] < chunkSizes[i])
            chunkSizes[i] = count[i];
    }

    NCDIE(nc_def_var(groupIdWrite, name, *xtype, *numDims, dimIds, variableWrite));
    NCDIE(nc_def_var_chunking(groupIdWrite, *variableWrite, chunkingStorage, chunkSizes));
    NCDIE(nc_def_var_deflate(groupIdWrite, *variableWrite, *shuffle, *deflate, *deflateLevel));

    copyVariableAttributes(groupIdRead, variableRead, groupIdWrite, *variableWrite);

    NCDIE(nc_inq_type(groupIdRead, *xtype, NULL, typeSize));
}

/**
 * @brief copy a piece of a variable and all of it's attributes to another file
 *
 * @param groupIdRead netCDF file or group to read
 * @param name name of the variable
 * @param start location to start copying from
 * @param count how many of each dimension to copy
 * @param dimIds dimension IDs from the destination file to attach to the new variable
 * @param groupIdWrite netCDF file or group to write the variable to
 * @param data write this data to the new variable or copy from read variable if NULL
 */
void copyVariable(int groupIdRead, const char *name, size_t *start, size_t *count, int *dimIds,
                  int groupIdWrite, void *data) {
    int status;
    int variableRead;
    int variableWrite;
    nc_type xtype;
    int numDims;
    int nAttributes;
    size_t chunkSizes[3]={0,0,0};

    size_t typeSize;
    int arraySize;
    int i;

    static int localDataSize = 0;
    static char *localData = NULL;

    int shuffle;
    int deflate;
    int deflateLevel;

    status = nc_inq_varid(groupIdRead, name, &variableRead);
    if (status == NC_NOERR) {
        setAttributes(name, groupIdRead, groupIdWrite, variableRead, &variableWrite, dimIds, count,
                      chunkSizes, &deflateLevel, &shuffle, &nAttributes, &deflate, &typeSize, &numDims,
                      &xtype);
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
            NCDIE(nc_get_vara(groupIdRead, variableRead, start, count, localData));
            NCDIE(nc_put_var(groupIdWrite, variableWrite, localData));
        } else {
            NCDIE(nc_put_var(groupIdWrite, variableWrite, data));
        }
    }
}

/**
 * @brief
 * copy a piece of a variable  and all of it's attributes to another file. Slicing is performed along a
 * selected dimension defined by a set of indexes along the dimension
 * @param groupIdRead netCDF file or group to read
 * @param name name of the variable
 * @param start location to start copying from
 * @param count how many of each dimension to copy
 * @param dimIds dimension IDs from the destination file to attach to the new
 * variable
 * @param groupIdWrite netCDF file or group to write the variable to
 * @param data write this data to the new variable or copy from read variable if
 * @param dimensionIndexes - index of the dimension along which slicing is performed.
 * @param indexes indexes along the dimension to be selected for output
 * NULL
 */
void copyVariableSelectedIndexes(int groupIdRead, const char *name, size_t *start, size_t *count, int *dimIds,
                                 int groupIdWrite, void *data, int *indexes, int dimensionIndexes) {
    if (dimensionIndexes < 0) {
        printf("Supplied dimensionIndexes is negative %d\n", dimensionIndexes);
        exit(EXIT_FAILURE);
    }
    size_t startLocal[3] = {0, 0, 0};
    size_t countLocal[3] = {1, 1, 1};
    size_t startGlobal[3];
    size_t countGlobal[3];
    size_t chunkSizes[3];
    for (size_t i = 0; i < 3; i++) {
        startGlobal[i] = start[i];
        countGlobal[i] = count[i];
    }
    int status;
    int variableRead;
    int variableWrite;
    nc_type xtype;
    int numDims;
    int nAttributes;

    size_t typeSize;
    int arraySize;
    int i;

    static int localDataSize = 0;
    static char *localData = NULL;

    int shuffle;
    int deflate;
    int deflateLevel;

    status = nc_inq_varid(groupIdRead, name, &variableRead);
    if (status == NC_NOERR) {
        setAttributes(name, groupIdRead, groupIdWrite, variableRead, &variableWrite, dimIds, count,
                      chunkSizes, &deflateLevel, &shuffle, &nAttributes, &deflate, &typeSize, &numDims,
                      &xtype);

        if (data == NULL) {
            // calc array size in num of bytes
            arraySize = typeSize;
            if (dimensionIndexes >= numDims) {
                printf(
                    "Supplied dimensionIndexes execedes the number of dimensions "
                    "%d, %d\n",
                    dimensionIndexes, numDims);
                exit(EXIT_FAILURE);
            }
            for (i = 0; i < numDims; i++) {
                if (i == dimensionIndexes)
                    continue;
                arraySize *= countGlobal[i];
                countLocal[i] = countGlobal[i];
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
            size_t total_size = countGlobal[dimensionIndexes];
            countGlobal[dimensionIndexes] = 1;
            for (int i_index = 0; i_index < total_size; i_index++) {
                startGlobal[dimensionIndexes] = indexes[i_index];
                startLocal[dimensionIndexes] = i_index;
                NCDIE(nc_get_vara(groupIdRead, variableRead, startGlobal, countGlobal, localData));
                NCDIE(nc_put_vara(groupIdWrite, variableWrite, startLocal, countLocal, localData));
            }
        } else {
            NCDIE(nc_put_var(groupIdWrite, variableWrite, data));
        }
    }
}

/**
 * @brief extract a L2 netCDF file
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
int extractNetCDF(const char *inFile, const char *outFile, int sPixel, int ePixel, int sScan, int eScan,
                  const char *productList, const char *waveList) {
    int status;
    int dsIdRead, dsIdWrite;  // data set file ids for reading and writing
    int groupRead;            // netCDF group for read file
    int groupWrite;           // netCDF group for write file
    int subGroupRead;         // netCDF sub group for read file
    int subGroupWrite;        // netCDF sub group for write file
    int dimId;
    int varId;
    size_t start[3] = {0, 0, 0};
    size_t count[3] = {1, 1, 1};
    int dimIds[3] = {0, 0, 0};

    size_t numLinesRead;                  // number of line in read file
    size_t numPixelsRead;                 // number of pixels in read file
    size_t numLinesWrite;                 // number of line in write file
    size_t numPixelsWrite;                // number of pixels in write file
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
    char name[256], title[256];

    // data sets read in
    float *latArray = NULL;
    float *lonArray = NULL;
    int *flagArray = NULL;

    int *wvValuesPass = NULL;
    int *wvIndexesPass = NULL;
    int wvNumPass = 0;

    getWave3Indexes(inFile, waveList, &wvValuesPass, &wvIndexesPass, &wvNumPass, productList);

    status = nc_open(inFile, NC_NOWRITE, &dsIdRead);
    if (status != NC_NOERR) {
        printf("could not open \"%s\" for reading.\n", inFile);
        exit(EXIT_FAILURE);
    }

    /* Get # of scans and # of pixels */
    /* ------------------------------ */
    status = nc_inq_dimid(dsIdRead, "number_of_lines", &dimId);
    if (status != NC_NOERR) {
        nc_inq_dimid(dsIdRead, "Number_of_Scan_Lines", &dimId);
    }
    nc_inq_dimlen(dsIdRead, dimId, &numLinesRead);

    status = nc_inq_dimid(dsIdRead, "pixels_per_line", &dimId);
    if (status != NC_NOERR) {
        nc_inq_dimid(dsIdRead, "Pixels_per_Scan_Line", &dimId);
    }
    nc_inq_dimlen(dsIdRead, dimId, &numPixelsRead);

    status = nc_inq_dimid(dsIdRead, "number_of_bands", &dimId);
    if (status != NC_NOERR) {
        nc_inq_dimid(dsIdRead, "total_band_number", &dimId);
    }
    nc_inq_dimlen(dsIdRead, dimId, &numTotalBands);

    status = nc_inq_dimid(dsIdRead, "number_of_reflective_bands", &dimId);
    if (status != NC_NOERR) {
        nc_inq_dimid(dsIdRead, "band_number", &dimId);
    }
    nc_inq_dimlen(dsIdRead, dimId, &numVisibleBands);

    status = nc_inq_dimid(dsIdRead, "bands_per_pixel", &dimId);
    if (status == NC_NOERR) {
        nc_inq_dimlen(dsIdRead, dimId, &numBandsPerPixel);
    } else {
        numBandsPerPixel = numTotalBands;
    }

    status = nc_inq_dimid(dsIdRead, "wavelength_3d", &dimId);
    if (status == NC_NOERR)
        nc_inq_dimlen(dsIdRead, dimId, &number_of_wv_3d);

    status = nc_inq_dimid(dsIdRead, "number_of_reflectance_location_values", &dimId);
    if (status == NC_NOERR) {
        nc_inq_dimlen(dsIdRead, dimId, &numReflectanceLocationValues);
        reflectanceLocationValuesFound = 1;
    }

    status = nc_get_att_text(dsIdRead, NC_GLOBAL, "title", title);
    if (strstr(title, "Level-2") == NULL) {
        printf("\"%s\" is not a L2 file.\n", inFile);
        exit(EXIT_FAILURE);
    }

    if (sScan < 1) {
        sScan = 1;
    }
    if (sScan >= numLinesRead) {
        printf("sscan needs to be less than number of scans in file.\n");
        exit(BOUNDS_ERROR);
    }
    if (eScan < 1 || eScan > numLinesRead) {
        eScan = numLinesRead;
    }
    if (eScan < sScan) {
        printf("escan needs to be greater than sscan.\n");
        exit(BOUNDS_ERROR);
    }
    numLinesWrite = eScan - sScan + 1;

    if (sPixel < 1) {
        sPixel = 1;
    }
    if (sPixel >= numPixelsRead) {
        printf("spix needs to be less than number of pixels in file.\n");
        exit(BOUNDS_ERROR);
    }
    if (ePixel < 1 || ePixel > numPixelsRead) {
        ePixel = numPixelsRead;
    }
    if (ePixel < sPixel) {
        printf("ePixel needs to be greater than spix.\n");
        exit(BOUNDS_ERROR);
    }

    // --------------------------------------------------------
    // make sure all products requested exist in the file
    // --------------------------------------------------------
    if (productList[0] != 0) {
        status = nc_inq_ncid(dsIdRead, "geophysical_data", &groupRead);
        if (status != NC_NOERR) {
            printf("-E- Could not open \"geophysical_data\" group.\n");
            exit(EXIT_FAILURE);
        }
        char *word;
        char *tmpProdList = strdup(productList);
        word = strtok(tmpProdList, ",");
        while (word) {
            if (word[0] != 0) {
                status = nc_inq_varid(groupRead, word, &varId);
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
    nc_inq_ncid(dsIdRead, "navigation_data", &groupRead);

    printf("sscan: %d  escan: %d\n", sScan, eScan);
    printf("spixl: %d  epixl: %d\n", sPixel, ePixel);

    numPixelsWrite = ePixel - sPixel + 1;

    // read lat and lon
    latArray = malloc(numLinesWrite * numPixelsWrite * sizeof(float));
    lonArray = malloc(numLinesWrite * numPixelsWrite * sizeof(float));
    if (latArray == NULL || lonArray == NULL) {
        printf("could not allocate memory for lat and lon arrays\n");
        exit(EXIT_FAILURE);
    }
    start[0] = sScan - 1;
    start[1] = sPixel - 1;
    count[0] = numLinesWrite;
    count[1] = numPixelsWrite;
    NCDIE(nc_inq_varid(groupRead, "latitude", &varId));
    NCDIE(nc_get_vara_float(groupRead, varId, start, count, latArray));
    NCDIE(nc_inq_varid(groupRead, "longitude", &varId));
    NCDIE(nc_get_vara_float(groupRead, varId, start, count, lonArray));

    // --------------------------------------------------------
    // set to geophysical data group
    // --------------------------------------------------------
    nc_inq_ncid(dsIdRead, "geophysical_data", &groupRead);

    flagArray = malloc(numLinesWrite * numPixelsWrite * sizeof(int));
    if (flagArray == NULL) {
        printf("could not allocate memory for flag array\n");
        exit(EXIT_FAILURE);
    }
    start[0] = sScan - 1;
    start[1] = sPixel - 1;
    count[0] = numLinesWrite;
    count[1] = numPixelsWrite;
    NCDIE(nc_inq_varid(groupRead, "l2_flags", &varId));
    NCDIE(nc_get_vara_int(groupRead, varId, start, count, flagArray));

    /* Calculate new navigation metadata */
    /* --------------------------------- */
#define GEOBOX_INC 20.0

    int startDone = 0;
    int centerDone = 0;
    int endDone = 0;
    int centerCol = numPixelsWrite / 2;
    int line, startCol, endCol;
    int lineStart;
    float lastLat = 0.0;
    int lastGoodLine = 0;
    float geoBox[4][100];
    int32_t geoBoxCount = 0;
    float gringFloatValue[100];
    int32_t gringIntValue[100];

    float startCenterLat, startCenterLon;
    float endCenterLat, endCenterLon;
    float geoLatMin, geoLatMax;
    float geoLonMin, geoLonMax;

    // loop to find beginning info
    for (line = 0; line < numLinesWrite; line++) {
        lineStart = numPixelsWrite * line;

        // find good start pixel
        if (!startDone) {
            for (startCol = 0; startCol <= centerCol; startCol++) {
                // check NAVFAIL flag
                if (flagArray[line * numPixelsWrite + startCol] ^ NAVFAIL) {
                    startDone = 1;
                    geoBox[0][geoBoxCount] = lonArray[lineStart + startCol];
                    geoBox[1][geoBoxCount] = latArray[lineStart + startCol];
                    break;
                }  // flag good
            }  // for col
        }  // if not start

        // find good center pixel
        if (!centerDone) {
            // check NAVFAIL flag
            if (flagArray[line * numPixelsWrite + centerCol] ^ NAVFAIL) {
                centerDone = 1;
                startCenterLon = lonArray[lineStart + centerCol];
                startCenterLat = latArray[lineStart + centerCol];
                lastLat = startCenterLat;
            }  // flag good
        }  // if not center

        // find good end pixel
        if (!endDone) {
            for (endCol = numPixelsWrite - 1; endCol >= centerCol; endCol--) {
                // check NAVFAIL flag
                if (flagArray[line * numPixelsWrite + endCol] ^ NAVFAIL) {
                    endDone = 1;
                    geoBox[2][geoBoxCount] = lonArray[lineStart + endCol];
                    geoBox[3][geoBoxCount] = latArray[lineStart + endCol];
                    break;
                }  // flag good
            }  // for col
        }  // if not start

        if (startDone && centerDone && endDone)
            break;
    }

    // set the min and max lat lon values
    geoLonMin = geoLonMax = geoBox[0][geoBoxCount];
    geoLatMin = geoLatMax = geoBox[1][geoBoxCount];
    if (geoLonMin > geoBox[2][geoBoxCount])
        geoLonMin = geoBox[2][geoBoxCount];
    if (geoLatMin > geoBox[3][geoBoxCount])
        geoLatMin = geoBox[3][geoBoxCount];
    if (geoLonMax < geoBox[2][geoBoxCount])
        geoLonMax = geoBox[2][geoBoxCount];
    if (geoLatMax < geoBox[3][geoBoxCount])
        geoLatMax = geoBox[3][geoBoxCount];
    if (geoLonMin > startCenterLon)
        geoLonMin = startCenterLon;
    if (geoLonMax < startCenterLon)
        geoLonMax = startCenterLon;

    geoBoxCount++;

    // loop through the rest of the lines
    for (; line < numLinesWrite; line++) {
        lineStart = numPixelsWrite * line;

        // find first good pixel on line
        for (startCol = 0; startCol <= centerCol; startCol++) {
            // check NAVFAIL flag
            if (flagArray[line * numPixelsWrite + startCol] ^ NAVFAIL)
                break;
        }
        if (startCol == centerCol)  // could not find start col, so skip this line
            continue;

        // find last good pixel
        for (endCol = numPixelsWrite - 1; endCol >= centerCol; endCol--) {
            // check NAVFAIL flag
            if (flagArray[line * numPixelsWrite + endCol] ^ NAVFAIL)
                break;
        }
        if (endCol < centerCol)  // could not find end col, so skip this line
            continue;

        lastGoodLine = line;

        // set min/max for every line
        if (geoLonMax < lonArray[lineStart + startCol])
            geoLonMax = lonArray[lineStart + startCol];
        if (geoLonMax < lonArray[lineStart + endCol])
            geoLonMax = lonArray[lineStart + endCol];
        if (geoLonMin > lonArray[lineStart + startCol])
            geoLonMin = lonArray[lineStart + startCol];
        if (geoLonMin > lonArray[lineStart + endCol])
            geoLonMin = lonArray[lineStart + endCol];

        if (geoLatMax < latArray[lineStart + startCol])
            geoLatMax = latArray[lineStart + startCol];
        if (geoLatMax < latArray[lineStart + endCol])
            geoLatMax = latArray[lineStart + endCol];
        if (geoLatMin > latArray[lineStart + startCol])
            geoLatMin = latArray[lineStart + startCol];
        if (geoLatMin > latArray[lineStart + endCol])
            geoLatMin = latArray[lineStart + endCol];

        // load up geobox
        if (fabs(lastLat - latArray[lineStart + centerCol]) > GEOBOX_INC) {
            geoBox[0][geoBoxCount] = lonArray[lineStart + startCol];
            geoBox[1][geoBoxCount] = latArray[lineStart + startCol];
            geoBox[2][geoBoxCount] = lonArray[lineStart + endCol];
            geoBox[3][geoBoxCount] = latArray[lineStart + endCol];
            lastLat = latArray[lineStart + centerCol];
            geoBoxCount++;
        }

    }  // for lines

    // make sure we add the last line
    lineStart = numPixelsWrite * lastGoodLine;

    // find first good pixel on line
    for (startCol = 0; startCol < centerCol; startCol++) {
        // check NAVFAIL flag
        if (flagArray[lastGoodLine * numPixelsWrite + startCol] ^ NAVFAIL)
            break;
    }

    // find last good pixel
    for (endCol = numPixelsWrite - 1; endCol >= centerCol; endCol--) {
        // check NAVFAIL flag
        if (flagArray[lastGoodLine * numPixelsWrite + endCol] ^ NAVFAIL)
            break;
    }

    endCenterLon = lonArray[lineStart + centerCol];
    endCenterLat = latArray[lineStart + centerCol];

    geoBox[0][geoBoxCount] = lonArray[lineStart + startCol];
    geoBox[1][geoBoxCount] = latArray[lineStart + startCol];
    geoBox[2][geoBoxCount] = lonArray[lineStart + endCol];
    geoBox[3][geoBoxCount] = latArray[lineStart + endCol];
    geoBoxCount++;

    /* Create output netCDF file (delete if it exists) */
    /* ------------------------- */
    if (access(outFile, F_OK) != -1) {
        if (unlink(outFile)) {
            printf("could not delete existing output file %s\n", outFile);
            exit(EXIT_FAILURE);
        }
    }

    status = nc_create(outFile, NC_NETCDF4, &dsIdWrite);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: Could not create NCDF4 file, %s .\n", __FILE__, __LINE__, outFile);
        exit(EXIT_FAILURE);
    }

    /* Create dimensions */
    /* ----------------- */
    NCDIE(nc_def_dim(dsIdWrite, "number_of_lines", numLinesWrite, &numLinesDimId));
    NCDIE(nc_def_dim(dsIdWrite, "pixels_per_line", numPixelsWrite, &numPixelsDimId));
    NCDIE(nc_def_dim(dsIdWrite, "bands_per_pixel", numBandsPerPixel, &numBandsPerPixelDimId));
    if (number_of_wv_3d > 0) {
        if (wvNumPass > 0) {
            number_of_wv_3d = wvNumPass;
        }
        NCDIE(nc_def_dim(dsIdWrite, "wavelength_3d", number_of_wv_3d, &numwv3dDimId));
    }

    if (reflectanceLocationValuesFound) {
        NCDIE(nc_def_dim(dsIdWrite, "number_of_reflectance_location_values", numReflectanceLocationValues,
                         &numReflectiveLocationValuesDimId));
    }
    NCDIE(nc_def_dim(dsIdWrite, "number_of_bands", numTotalBands, &numTotalBandsDimId));
    NCDIE(nc_def_dim(dsIdWrite, "number_of_reflective_bands", numVisibleBands, &numVisibleBandsDimId));

    // --------------------------------------------------------
    // set to sensor_band_parameters group
    // --------------------------------------------------------
    nc_inq_ncid(dsIdRead, "sensor_band_parameters", &groupRead);
    nc_def_grp(dsIdWrite, "sensor_band_parameters", &groupWrite);

    // copy vcal_gain(visibleBands)
    start[0] = 0;
    count[0] = numVisibleBands;
    dimIds[0] = numVisibleBandsDimId;

    // loop through all of the variables
    status = nc_inq_varids(groupRead, &numVariables, variableIds);
    if (status == NC_NOERR) {
        if (numVariables > MAX_VARIABLES) {
            printf("-E- %s %d: too many variables in group \"sensor_band_parameters\"\n", __FILE__, __LINE__);
            exit(EXIT_SUCCESS);
        }
        for (i = 0; i < numVariables; i++) {
            status = nc_inq_varname(groupRead, variableIds[i], name);
            if ((strcmp(name, "wavelength") == 0)) {
                count[0] = numTotalBands;
                dimIds[0] = numTotalBandsDimId;
                copyVariable(groupRead, name, start, count, dimIds, groupWrite, NULL);
            } else {
                if ((strcmp(name, "wavelength_3d") == 0)) {
                    count[0] = number_of_wv_3d;
                    dimIds[0] = numwv3dDimId;
                    if (wvNumPass > 0) {
                        copyVariableSelectedIndexes(groupRead, name, start, count, dimIds, groupWrite, NULL,
                                                    wvIndexesPass, 0);
                        continue;
                    }
                } else {
                    count[0] = numVisibleBands;
                    dimIds[0] = numVisibleBandsDimId;
                }
                copyVariable(groupRead, name, start, count, dimIds, groupWrite, NULL);
            }
        }
    }

    // --------------------------------------------------------
    // set to sensor_band_parameters group
    // --------------------------------------------------------
    nc_inq_ncid(dsIdRead, "scan_line_attributes", &groupRead);
    nc_def_grp(dsIdWrite, "scan_line_attributes", &groupWrite);

    // set up the copy parameters
    start[0] = sScan - 1;
    count[0] = numLinesWrite;
    dimIds[0] = numLinesDimId;

    // loop through all of the variables
    status = nc_inq_varids(groupRead, &numVariables, variableIds);
    if (status == NC_NOERR) {
        if (numVariables > MAX_VARIABLES) {
            printf("-E- %s %d: too many variables in group \"scan_line_attributes\"\n", __FILE__, __LINE__);
            exit(EXIT_SUCCESS);
        }
        for (i = 0; i < numVariables; i++) {
            status = nc_inq_varname(groupRead, variableIds[i], name);
            copyVariable(groupRead, name, start, count, dimIds, groupWrite, NULL);
        }
    }

    // --------------------------------------------------------
    // set to geophysical_data group
    // --------------------------------------------------------
    nc_inq_ncid(dsIdRead, "geophysical_data", &groupRead);
    nc_def_grp(dsIdWrite, "geophysical_data", &groupWrite);

    // set up the copy parameters
    start[0] = sScan - 1;
    start[1] = sPixel - 1;
    count[0] = numLinesWrite;
    count[1] = numPixelsWrite;
    dimIds[0] = numLinesDimId;
    dimIds[1] = numPixelsDimId;

    // loop through all of the variables
    status = nc_inq_varids(groupRead, &numVariables, variableIds);
    if (status == NC_NOERR) {
        if (numVariables > MAX_VARIABLES) {
            printf("-E- %s %d: too many variables in group \"geophysical_data\"\n", __FILE__, __LINE__);
            exit(EXIT_SUCCESS);
        }
        for (i = 0; i < numVariables; i++) {
            status = nc_inq_varname(groupRead, variableIds[i], name);
            {
                int varids_dims[3] = {0, 0, 0};
                char dim_name[50];
                nc_inq_vardimid(groupRead, variableIds[i], varids_dims);
                nc_inq_dimname(groupRead, varids_dims[2], dim_name);
                if (strcmp(dim_name, "wavelength_3d") == 0) {
                    dimIds[2] = numwv3dDimId;
                    count[2] = number_of_wv_3d;
                    if (wvIndexesPass > 0) {
                        if (checkIfInProdlist(productList, name) == 0)
                            continue;
                        copyVariableSelectedIndexes(groupRead, name, start, count, dimIds, groupWrite, NULL,
                                                    wvIndexesPass, 2);
                        continue;
                    }
                } else {
                    dimIds[2] = 0;
                    count[2] = 1;
                }
            }
            if (checkIfInProdlist(productList, name) == 0)
                continue;
            copyVariable(groupRead, name, start, count, dimIds, groupWrite, NULL);
        }
    }

    // --------------------------------------------------------
    // set to navigation data group
    // --------------------------------------------------------
    nc_inq_ncid(dsIdRead, "navigation_data", &groupRead);
    nc_def_grp(dsIdWrite, "navigation_data", &groupWrite);

    start[0] = sScan - 1;
    start[1] = sPixel - 1;
    count[0] = numLinesWrite;
    count[1] = numPixelsWrite;
    dimIds[0] = numLinesDimId;
    dimIds[1] = numPixelsDimId;
    copyVariable(groupRead, "longitude", start, count, dimIds, groupWrite, lonArray);
    copyVariable(groupRead, "latitude", start, count, dimIds, groupWrite, latArray);

    start[0] = sScan - 1;
    count[0] = numLinesWrite;
    dimIds[0] = numLinesDimId;
    copyVariable(groupRead, "tilt", start, count, dimIds, groupWrite, NULL);

    // fill up gring
    int j = 1;
    gringFloatValue[0] = geoBox[0][0];
    for (i = 0; i < geoBoxCount; i++) {
        gringFloatValue[j++] = geoBox[2][i];
    }
    for (i = 0; i < geoBoxCount - 1; i++) {
        gringFloatValue[j++] = geoBox[0][geoBoxCount - 1 - i];
    }
    NCDIE(nc_put_att_float(groupWrite, NC_GLOBAL, "gringpointlongitude", NC_FLOAT, j, gringFloatValue));

    j = 1;
    gringFloatValue[0] = geoBox[1][0];
    gringIntValue[0] = j;
    for (i = 0; i < geoBoxCount; i++) {
        gringIntValue[j] = j + 1;
        gringFloatValue[j++] = geoBox[3][i];
    }
    for (i = 0; i < geoBoxCount - 1; i++) {
        gringIntValue[j] = j + 1;
        gringFloatValue[j++] = geoBox[1][geoBoxCount - 1 - i];
    }
    NCDIE(nc_put_att_float(groupWrite, NC_GLOBAL, "gringpointlatitude", NC_FLOAT, j, gringFloatValue));
    NCDIE(nc_put_att_int(groupWrite, NC_GLOBAL, "gringpointsequence", NC_INT, j, gringIntValue));

    // --------------------------------------------------------
    // set to processing_control group
    // --------------------------------------------------------
    nc_inq_ncid(dsIdRead, "processing_control", &groupRead);
    nc_def_grp(dsIdWrite, "processing_control", &groupWrite);

    copyGlobalAttributes(groupRead, groupWrite);

    // sub group input_parameters
    nc_inq_ncid(groupRead, "input_parameters", &subGroupRead);
    nc_def_grp(groupWrite, "input_parameters", &subGroupWrite);

    copyGlobalAttributes(subGroupRead, subGroupWrite);

    // sub group flag_percentages
    nc_inq_ncid(groupRead, "flag_percentages", &subGroupRead);
    nc_def_grp(groupWrite, "flag_percentages", &subGroupWrite);

    copyGlobalAttributes(subGroupRead, subGroupWrite);

    // --------------------------------------------------------
    // copy global attributes
    // --------------------------------------------------------
    copyGlobalAttributes(dsIdRead, dsIdWrite);

    // write modified global attrbutes
    NCDIE(nc_put_att_float(dsIdWrite, NC_GLOBAL, "start_center_longitude", NC_FLOAT, 1, &startCenterLon));
    NCDIE(nc_put_att_float(dsIdWrite, NC_GLOBAL, "start_center_latitude", NC_FLOAT, 1, &startCenterLat));
    NCDIE(nc_put_att_float(dsIdWrite, NC_GLOBAL, "end_center_longitude", NC_FLOAT, 1, &endCenterLon));
    NCDIE(nc_put_att_float(dsIdWrite, NC_GLOBAL, "end_center_latitude", NC_FLOAT, 1, &endCenterLat));

    NCDIE(nc_put_att_float(dsIdWrite, NC_GLOBAL, "northernmost_latitude", NC_FLOAT, 1, &geoLatMax));
    NCDIE(nc_put_att_float(dsIdWrite, NC_GLOBAL, "southernmost_latitude", NC_FLOAT, 1, &geoLatMin));
    NCDIE(nc_put_att_float(dsIdWrite, NC_GLOBAL, "easternmost_longitude", NC_FLOAT, 1, &geoLonMax));
    NCDIE(nc_put_att_float(dsIdWrite, NC_GLOBAL, "westernmost_longitude", NC_FLOAT, 1, &geoLonMin));

    NCDIE(nc_put_att_float(dsIdWrite, NC_GLOBAL, "geospatial_lat_max", NC_FLOAT, 1, &geoLatMax));
    NCDIE(nc_put_att_float(dsIdWrite, NC_GLOBAL, "geospatial_lat_min", NC_FLOAT, 1, &geoLatMin));
    NCDIE(nc_put_att_float(dsIdWrite, NC_GLOBAL, "geospatial_lon_max", NC_FLOAT, 1, &geoLonMax));
    NCDIE(nc_put_att_float(dsIdWrite, NC_GLOBAL, "geospatial_lon_min", NC_FLOAT, 1, &geoLonMin));
    nc_close(dsIdRead);
    nc_close(dsIdWrite);

    free(latArray);
    free(lonArray);
    free(flagArray);

    return EXIT_SUCCESS;
}
