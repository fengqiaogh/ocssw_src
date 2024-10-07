#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include <clo.h>
#include <hdf_utils.h>
#include <genutils.h>

#define PROGRAM "ice2hdf"
#define VERSION "1.2"
#define FMT_L2HDF 2

#define DATA_CITATION "Meier, W., F. Fetterer, K. Knowles, M. Savoie, M. J. Brodzik. 2006, updated quarterly. Sea ice concentrations from Nimbus-7 SMMR and DMSP SSM/I passive microwave data. Boulder, Colorado USA: National Snow and Ice Data Center. Digital media."
#define WEB_SITE "http://nsidc.org/data/docs/daac/nsidc0051_gsfc_seaice.gd.html"
#define FTP_SITE "ftp://sidads.colorado.edu/pub/DATASETS/seaice/polar-stereo/nasateam/near-real-time"

#define HEADER_SIZE 300

#define NORTH_ROWS 448
#define NORTH_COLS 304
#define SOUTH_ROWS 332
#define SOUTH_COLS 316

#define MAX_DATA     250
#define POLAR_HOLE   251
#define MISSING_DATA 255

// data file structure

typedef struct data_file_t {
    char* filename;
    uint8 *dataPtr;
    int32 numRows;
    int32 numCols;
    int hasHeader;
    int32 year;
    int32 doy;
} data_file_t;

/* -------------------------------------------------------------------- */
/*                         allocateDataFile                             */

/* -------------------------------------------------------------------- */
data_file_t* allocateDataFile(char* filename, int32 numRows, int32 numCols) {
    data_file_t* dataFile = (data_file_t*) malloc(sizeof (data_file_t));

    dataFile->filename = strdup(filename);
    dataFile->dataPtr = (uint8*) malloc(numRows * numCols);
    dataFile->numRows = numRows;
    dataFile->numCols = numCols;
    dataFile->hasHeader = 0;
    dataFile->year = 0;
    dataFile->doy = 0;

    return dataFile;
}

/* -------------------------------------------------------------------- */
/*                            readFile                                  */

/* -------------------------------------------------------------------- */
void readFile(data_file_t* dataFile) {
    char header[HEADER_SIZE];
    int result;
    FILE *fin;

    fin = fopen(dataFile->filename, "rb");
    if (fin == NULL) {
        fprintf(stderr, "-E- %s line %d: Could not open file, \"%s\" .\n",
                __FILE__, __LINE__, dataFile->filename);
        exit(EXIT_FAILURE);
    }

    if (dataFile->hasHeader) {
        result = fread(header, 1, HEADER_SIZE, fin);
        if (result != HEADER_SIZE) {
            fprintf(stderr,
                    "-E- %s line %d: Error reading file header, \"%s\" .\n",
                    __FILE__, __LINE__, dataFile->filename);
            exit(EXIT_FAILURE);
        }

        //
        // parse the header
        //
        int32 numRows = atoi(&header[12]);
        int32 numCols = atoi(&header[6]);
        dataFile->year = atoi(&header[102]);
        dataFile->doy = atoi(&header[108]);

        if (numRows != dataFile->numRows) {
            fprintf(stderr,
                    "-E- %s line %d: Bad number of rows(%d) in \"%s\" should be %d\n",
                    __FILE__, __LINE__, numRows, dataFile->filename,
                    dataFile->numRows);
            exit(EXIT_FAILURE);
        }

        if (numCols != dataFile->numCols) {
            fprintf(stderr,
                    "-E- %s line %d: Bad number of columns(%d) in \"%s\" should be %d\n",
                    __FILE__, __LINE__, numCols, dataFile->filename,
                    dataFile->numCols);
            exit(EXIT_FAILURE);
        }

    }

    //
    // read the data
    //
    result = fread(dataFile->dataPtr, 1, dataFile->numRows * dataFile->numCols,
            fin);
    if (result != dataFile->numRows * dataFile->numCols) {
        fprintf(stderr, "-E- %s line %d: Error reading file data, \"%s\" \n",
                __FILE__, __LINE__, dataFile->filename);
        exit(EXIT_FAILURE);
    }

}

/* -------------------------------------------------------------------- */
/* maskData                                                             */
/*   Set the data to 255 anywhere the mask dataset is set to the        */
/*   given values.                                                      */

/* -------------------------------------------------------------------- */
void maskData(data_file_t* dataFile, data_file_t* maskFile, int numVals,
        int* vals) {
    if (dataFile->numRows != maskFile->numRows) {
        fprintf(stderr,
                "-E- %s line %d: data file \"%s\" and mask file \"%s\" have different number of Rows\n",
                __FILE__, __LINE__, dataFile->filename, maskFile->filename);
        exit(EXIT_FAILURE);
    }
    if (dataFile->numCols != maskFile->numCols) {
        fprintf(stderr,
                "-E- %s line %d: data file \"%s\" and mask file \"%s\" have different number of Columns\n",
                __FILE__, __LINE__, dataFile->filename, maskFile->filename);
        exit(EXIT_FAILURE);
    }

    int i, j;
    for (i = 0; i < dataFile->numRows * dataFile->numCols; i++) {
        for (j = 0; j < numVals; j++) {
            if (maskFile->dataPtr[i] == vals[j]) {
                dataFile->dataPtr[i] = MISSING_DATA;
            }
        }
    }
}

/* -------------------------------------------------------------------- */
/* fillHoles                                                            */
/*   Fill the holes in the data by iterativly averaging the surrounding */
/*   pixels.  Assume anything over 250 is a hole.                       */
/*                                                                      */

/* -------------------------------------------------------------------- */
void fillHoles(data_file_t* dataFile) {
    int delta = 1; // size of the box around the pixel

    uint8 * data2;
    int holesFound = 1;
    int row, col;
    int r, c;
    int r1, r2;
    int c1, c2;
    int count;
    int sum;

    data2 = (uint8*) malloc(dataFile->numRows * dataFile->numCols);
    if (data2 == NULL) {
        fprintf(stderr, "-E- %s line %d: Could not malloc.\n", __FILE__,
                __LINE__);
        exit(EXIT_FAILURE);
    }

    while (holesFound) {
        holesFound = 0;

        // copy the data so we don't key off of newly filled values
        memcpy(data2, dataFile->dataPtr, dataFile->numRows * dataFile->numCols);

        for (row = 0; row < dataFile->numRows; row++) {
            for (col = 0; col < dataFile->numCols; col++) {

                if (data2[row * dataFile->numCols + col] > MAX_DATA) {
                    holesFound = 1;

                    // calc limits of the box around the pixel
                    r1 = row - delta;
                    if (r1 < 0)
                        r1 = 0;
                    r2 = row + delta;
                    if (r2 >= dataFile->numRows)
                        r2 = dataFile->numRows - 1;

                    c1 = col - delta;
                    if (c1 < 0)
                        c1 = 0;
                    c2 = col + delta;
                    if (c2 >= dataFile->numCols)
                        c2 = dataFile->numCols - 1;

                    // intialize the average variables
                    count = 0;
                    sum = 0;

                    // sum surrounding pixels
                    for (r = r1; r <= r2; r++) {
                        for (c = c1; c <= c2; c++) {
                            if (data2[r * dataFile->numCols + c] <= MAX_DATA) {
                                count++;
                                sum += data2[r * dataFile->numCols + c];
                            }
                        } // for c1 to c2
                    } // for r1 to r2

                    if (count > 1) {
                        dataFile->dataPtr[row * dataFile->numCols + col] = sum
                                / count;
                    }

                } // pixel is a hole
            } // for cols
        } // for rows

    } // while holes

    free(data2);
}

/* -------------------------------------------------------------------- */
/*                            main                                      */

/* -------------------------------------------------------------------- */
int main(int argc, char* argv[]) {
    clo_optionList_t * list;
    clo_option_t* option;

    char tmpStr[2048];
    int initialNumOptions;
    int numExtraOptions;

    data_file_t *northData;
    data_file_t *southData;
    int northNumRegions;
    int* northRegionList;
    data_file_t *northRegion;
    int southNumRegions;
    int* southRegionList;
    data_file_t *southRegion;

    int32 hdfStat;
    idDS ds_id;

    char northFilename[1024];
    char southFilename[1024];
    char outFilename[1024];
    char northRegionFilename[1024];
    char southRegionFilename[1024];

    sprintf(tmpStr, "ice2hdf %s (%s %s)", VERSION, __DATE__, __TIME__);
    clo_setVersion(tmpStr);

    list = clo_createList();

    strcpy(tmpStr, "Usage: ice2hdf [options] northFile southFile outFile\n\n");
    strcat(tmpStr,
            "Program to convert the NSIDC ice data into an HDF file that l2gen can use.\n");
    strcat(tmpStr, "Documentation - ");
    strcat(tmpStr, WEB_SITE);
    strcat(tmpStr, "\nData - ");
    strcat(tmpStr, FTP_SITE);
    strcat(tmpStr, "\n");
    clo_setHelpStr(tmpStr);

    clo_addOption(list, "fill_pole", CLO_TYPE_BOOL, "true", "Fill the Arctic polar hole with 100% ice");
    strcpy(tmpStr, "Mask the values in the given\n        northern regions\n");
    strcat(tmpStr, "        0: Lakes, extended coast\n");
    strcat(tmpStr, "        1: Non-regional ocean\n");
    strcat(tmpStr, "        2: Sea of Okhotsk and Japan\n");
    strcat(tmpStr, "        3: Bering Sea\n");
    strcat(tmpStr, "        4: Hudson Bay\n");
    strcat(tmpStr, "        5: Baffin Bay/Davis Strait/Labrador Sea\n");
    strcat(tmpStr, "        6: Greenland Sea\n");
    strcat(tmpStr, "        7: Kara and Barents Seas\n");
    strcat(tmpStr, "        8: Arctic Ocean\n");
    strcat(tmpStr, "        9: Canadian Archipelago\n");
    strcat(tmpStr, "       10: Gulf of St. Lawrence\n");
    strcat(tmpStr, "       11: Land\n");
    strcat(tmpStr, "       12: Coast");
    clo_addOption(list, "region_mask_north", CLO_TYPE_INT, "[0,11,12]", tmpStr);
    clo_addOption(list, "region_file_north", CLO_TYPE_IFILE,
            "$OCDATAROOT/common/region_n.msk",
            "\n        Region mask file for the northern hemisphere");

    strcpy(tmpStr, "Mask the values in the given\n        southern regions\n");
    strcat(tmpStr, "        2: Weddell Sea\n");
    strcat(tmpStr, "        3: Indian Ocean\n");
    strcat(tmpStr, "        4: Pacific Ocean\n");
    strcat(tmpStr, "        5: Ross Sea\n");
    strcat(tmpStr, "        6: Bellingshausen Amundsen Sea\n");
    strcat(tmpStr, "       11: Land\n");
    strcat(tmpStr, "       12: Coast");
    clo_addOption(list, "region_mask_south", CLO_TYPE_INT, "[11,12]", tmpStr);
    clo_addOption(list, "region_file_south", CLO_TYPE_IFILE,
            "$OCDATAROOT/common/region_s.msk",
            "\n        Region mask file for the southern hemisphere");

    initialNumOptions = clo_getNumOptions(list);
    clo_readArgs(list, argc, argv);
    numExtraOptions = clo_getNumOptions(list) - initialNumOptions;

    if (numExtraOptions != 3) {
        printf("ERROR - wrong number of parameters given\n\n");
        clo_printUsage(list);
        return 1;
    }

    option = clo_getOption(list, initialNumOptions);
    parse_file_name(option->key, northFilename);

    option = clo_getOption(list, initialNumOptions + 1);
    parse_file_name(option->key, southFilename);

    option = clo_getOption(list, initialNumOptions + 2);
    parse_file_name(option->key, outFilename);

    northRegionList = clo_getInts(list, "region_mask_north", &northNumRegions);
    parse_file_name(clo_getString(list, "region_file_north"),
            northRegionFilename);

    southRegionList = clo_getInts(list, "region_mask_south", &southNumRegions);
    parse_file_name(clo_getString(list, "region_file_south"),
            southRegionFilename);

    time_t currentTime = time(NULL);
    char* currentTimeStr = strdup(asctime(gmtime(&currentTime)));
    currentTimeStr[24] = '\0';

    //
    // read the data files
    //
    northData = allocateDataFile(northFilename, NORTH_ROWS, NORTH_COLS);
    northData->hasHeader = 1;
    readFile(northData);

    if (clo_getBool(list, "fill_pole")) {
        int i;
        for (i = 0; i < NORTH_ROWS * NORTH_COLS; i++) {
            if (northData->dataPtr[i] == POLAR_HOLE) {
                northData->dataPtr[i] = MAX_DATA;
            }
        }
    }

    southData = allocateDataFile(southFilename, SOUTH_ROWS, SOUTH_COLS);
    southData->hasHeader = 1;
    readFile(southData);

    if (northNumRegions > 0) {
        northRegion = allocateDataFile(northRegionFilename, NORTH_ROWS,
                NORTH_COLS);
        northRegion->hasHeader = 1;
        readFile(northRegion);
        maskData(northData, northRegion, northNumRegions, northRegionList);
    }

    if (southNumRegions > 0) {
        southRegion = allocateDataFile(southRegionFilename, SOUTH_ROWS,
                SOUTH_COLS);
        southRegion->hasHeader = 1;
        readFile(southRegion);
        maskData(southData, southRegion, southNumRegions, southRegionList);
    }

    fillHoles(northData);
    fillHoles(southData);

    //
    // open the HDF file
    //
    ds_id = startDS(outFilename, DS_HDF, DS_WRITE, 0);
    if (ds_id.fid == FAIL) {
        fprintf(stderr, "-E- %s line %d: Could not create HDF file, %s .\n",
                __FILE__, __LINE__, outFilename);
        return (HDF_FUNCTION_ERROR);
    }

    /*                                                                  */
    /* Create the north and south data array SDSes                      */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    PTB(createDS(ds_id, /* file id         */
            "north", /* short name      */
            "North Pole Ice Fraction", /* long name       */
            NULL, /* standard name   */
            NULL, /* units           */
            0, 0, /* range           */
            0, 0, /* slope           */
            DFNT_INT8, /* HDF number type */
            2, /* rank            */
            northData->numRows, northData->numCols, 0, /* dimension sizes */
            "north_rows", "north_cols", NULL /* dimension names */
            ));

    PTB(createDS(ds_id, /* file id         */
            "south", /* short name      */
            "South Pole Ice Fraction", /* long name       */
            NULL, /* standard name   */
            NULL, /* units           */
            0, 0, /* range           */
            0, 0, /* slope           */
            DFNT_INT8, /* HDF number type */
            2, /* rank            */
            southData->numRows, southData->numCols, 0, /* dimension sizes */
            "south_rows", "south_cols", NULL /* dimension names */
            ));

    // Modify to use new HDF4/NCDF4 libhdfutils library  JMG  09/28/12
    PTB(SetChrGA(ds_id, "Title", "Sea Ice Concentration Daily"));
    PTB(SetChrGA(ds_id, "Creation Date", currentTimeStr));
    PTB(SetChrGA(ds_id, "Created By", "Don Shea, SAIC, from NSIDC NRT data (http://nsidc.org)"));
    PTB(SetChrGA(ds_id, "Reference1", DATA_CITATION));
    PTB(SetChrGA(ds_id, "Software Name", PROGRAM));
    PTB(SetChrGA(ds_id, "Software Version", VERSION));
    PTB(SetI32GA(ds_id, "year", northData->year));
    PTB(SetI32GA(ds_id, "day_of_year", northData->doy));

    PTB(sd_writedata(ds_id.fid, "north", northData->dataPtr, 0, 0, 0, northData->numRows, northData->numCols, 0));
    PTB(sd_writedata(ds_id.fid, "south", southData->dataPtr, 0, 0, 0, southData->numRows, southData->numCols, 0));

    //
    // close the HDF file
    //
    hdfStat = endDS(ds_id);
    if (hdfStat == FAIL) {
        fprintf(stderr, "-E- %s line %d: Could not close HDF file, %s .\n",
                __FILE__, __LINE__, outFilename);
        return (HDF_FUNCTION_ERROR);
    }

    return (LIFE_IS_GOOD);
}
