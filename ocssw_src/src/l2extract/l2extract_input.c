#include "l2extract.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <genutils.h>
#include <clo.h>
#include <dfutils.h>
#include <sensorInfo.h>
#include <unistd.h>

/** add all of the accepted command line options to list */
int l2extract_init_options(clo_optionList_t* list, const char* softwareVersion) {
    clo_setVersion2("l2extract", softwareVersion);

    clo_setHelpStr(
        "Usage: l2extract argument-list"
        "\n   or: l2extract ifile spix epix sline eline pix_sub sc_sub ofile <prodlist>"
        "\n"
        "\n  This program takes a product (or products if netCDF output) from an L2 file"
        "\n  and does the extraction"
        "\n"
        "\n  The argument list is a set of keyword=value pairs.  Arguments can"
        "\n  be specified on the command line, or put into a parameter file, or the"
        "\n  two methods can be used together, with command line overriding."
        "\n"
        "\nThe list of valid keywords follows:"
        "\n");

    // files
    clo_addOption(list, "ifile", CLO_TYPE_IFILE, NULL, "input L2 filename");
    clo_addOption(list, "ofile", CLO_TYPE_OFILE, "output", "output filename");

    // data, projection, coverage
    clo_addOption(list, "product", CLO_TYPE_STRING, NULL,
                  "comma separated list of products, empty string outputs all\n        all products");
    clo_addOption(list, "spix", CLO_TYPE_INT, "1", "start pixel number (1-based).");
    clo_addOption(list, "epix", CLO_TYPE_INT, "-1", "end pixel number. -1 = last pixel (1-based).");
    clo_addOption(list, "sline", CLO_TYPE_INT, "1", "start line (1-based).");
    clo_addOption(list, "eline", CLO_TYPE_INT, "-1", "end line.  -1 = last line (1-based).");
    // clo_addOption(list, "pix_sub", CLO_TYPE_INT, "1", "pixel subsampling rate.");
    // clo_addOption(list, "sc_sub", CLO_TYPE_INT, "1", "scan line subsampling rate.");
    clo_addOption(list, "suite", CLO_TYPE_STRING, NULL, "suite for default parameters");
    clo_addOption(list, "verbose", CLO_TYPE_BOOL, "false", "more information reporting");
    clo_addOption(list, "wavelist", CLO_TYPE_STRING, NULL,
                  "comma separated list of 3D wavelengths and/or colon\n"
                  "        separated nnn:nnn range of wavelengths, empty string outputs all\n"
                  "        3D wavelengths (i.e. wavelist=353,355,358,360,360:370,450:600,700)");
    return 0;
}

int getSensorId(const char* fileName) {
    idDS dsId;
    int sensorId = -1;

    dsId = openDS(fileName);
    if (dsId.fid == -1) {
        printf("-E- %s: Input file '%s' does not exist or cannot open.  \n", __FILE__, fileName);
        exit(EXIT_FAILURE);
    }

    if (dsId.fftype == DS_NCDF) {
        char* instrumentStr = readAttrStr(dsId, "instrument");
        if (instrumentStr) {
            char* platformStr = readAttrStr(dsId, "platform");
            if (platformStr) {
                sensorId = instrumentPlatform2SensorId(instrumentStr, platformStr);
                free(platformStr);
            }
            free(instrumentStr);
        } else {
            char* sensorStr = readAttrStr(dsId, "Sensor");
            if (sensorStr) {
                sensorId = sensorName2SensorId(sensorStr);
                free(sensorStr);
            }
        }
    } else {
        char* sensorNameStr = readAttrStr(dsId, "Sensor Name");
        if (sensorNameStr) {
            sensorId = sensorName2SensorId(sensorNameStr);
            free(sensorNameStr);
        }
    }

    if (sensorId == -1) {
        printf("Did not find a valid sensor ID - using OCRVC as the sensor ID.\n");
        sensorId = OCRVC;
    }

    endDS(dsId);
    return sensorId;
}

/*
 Read the command line option and all of the default parameter files.

 This is the order for loading the options:
 - load the main program defaults file
 - load the command line (including specified par files)
 - re-load the command line disabling file descending so command
 line arguments will over ride

 */
int l2extract_read_options(clo_optionList_t* list, int argc, char* argv[]) {
    char* dataRoot;
    char tmpStr[FILENAME_MAX];

    assert(list);

    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        fprintf(stderr, "-E- OCDATAROOT environment variable is not defined.  \n");
        return (-1);
    }

    // disable the dump option until we have read all of the files
    clo_setEnableDumpOptions(0);

    // load program defaults
    sprintf(tmpStr, "%s/common/l2extract_defaults.par", dataRoot);
    if (access(tmpStr, R_OK) != -1) {
        if (want_verbose)
            printf("Loading program default parameters from %s\n", tmpStr);
        clo_readFile(list, tmpStr);
    }

    // read all arguments
    clo_readArgs(list, argc, argv);

    // get sensor directory
    char* ifileStr;
    if (clo_getPositionNumOptions(list) == 0) {
        ifileStr = clo_getString(list, "ifile");
    } else {
        ifileStr = clo_getPositionString(optionList, 0);
    }

    int sensorId = getSensorId(ifileStr);
    int subsensorId = sensorId2SubsensorId(sensorId);
    const char* sensorDir = sensorId2SensorDir(sensorId);

    // load the sensor specific defaults file
    sprintf(tmpStr, "%s/%s/l2extract_defaults.par", dataRoot, sensorDir);
    // int defaultLoaded = 0;
    if (access(tmpStr, R_OK) != -1) {
        if (want_verbose)
            printf("Loading default parameters from %s\n", tmpStr);
        clo_readFile(list, tmpStr);
    }

    if (subsensorId != -1) {
        sprintf(tmpStr, "%s/%s/%s/l2extract_defaults.par", dataRoot, sensorDir,
                subsensorId2SubsensorDir(subsensorId));
        if (access(tmpStr, R_OK) != -1) {
            if (want_verbose)
                printf("Loading default parameters from %s\n", tmpStr);
            clo_readFile(list, tmpStr);
        }
    }

    // load the suite specific defaults file
    clo_option_t* option = clo_findOption(list, "suite");
    if (clo_isOptionSet(option)) {
        int suiteLoaded = 0;
        const char* suiteStr = clo_getOptionString(option);

        // common suite
        sprintf(tmpStr, "%s/common/l2extract_defaults_%s.par", dataRoot, suiteStr);
        if (access(tmpStr, R_OK) != -1) {
            if (want_verbose)
                printf("Loading default parameters from %s\n", tmpStr);
            clo_readFile(list, tmpStr);
            suiteLoaded = 1;
        }

        // sensor suite
        sprintf(tmpStr, "%s/%s/l2extract_defaults_%s.par", dataRoot, sensorDir, suiteStr);
        if (access(tmpStr, R_OK) != -1) {
            if (want_verbose)
                printf("Loading default parameters from %s\n", tmpStr);
            clo_readFile(list, tmpStr);
            suiteLoaded = 1;
        }

        // sub-sensor suite
        if (subsensorId != -1) {
            sprintf(tmpStr, "%s/%s/%s/l2extract_defaults_%s.par", dataRoot, sensorDir,
                    subsensorId2SubsensorDir(subsensorId), suiteStr);
            if (access(tmpStr, R_OK) != -1) {
                if (want_verbose)
                    printf("Loading default parameters from %s\n", tmpStr);
                clo_readFile(list, tmpStr);
                suiteLoaded = 1;
            }
        }

        if (!suiteLoaded) {
            printf("-E- Failed to load parameters for suite %s for sensor %s\n", suiteStr,
                   sensorId2SensorName(sensorId));
            exit(EXIT_FAILURE);
        }
    }
    // enable the dump option
    clo_setEnableDumpOptions(1);
    // make the command line over ride
    clo_readArgs(list, argc, argv);

    return 0;
}
