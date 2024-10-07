#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <genutils.h>
#include <clo.h>
#include <sensorInfo.h>

#include "l2brsgen.h"

//-----------------------------------------------------------------------

/** add all of the accepted command line options to list */
int l2brsgen_init_options(clo_optionList_t* list) {
    char tmpStr[2048];

    sprintf(tmpStr, "l2brsgen %s (%s %s)", VERSION, __DATE__, __TIME__);
    clo_setVersion(tmpStr);

    sprintf(tmpStr, "Usage: l2brsgen argument-list\n\n");

    strcat(tmpStr, "  This program takes a product from a L2 file, subsamples the file\n");
    strcat(tmpStr, "  and writes a browse file\n\n");

    strcat(tmpStr, "  The argument-list is a set of keyword=value pairs. The arguments can\n");
    strcat(tmpStr, "  be specified on the commandline, or put into a parameter file, or the\n");
    strcat(tmpStr, "  two methods can be used together, with commandline over-riding.\n\n");
    strcat(tmpStr, "The list of valid keywords follows:\n");
    clo_setHelpStr(tmpStr);

    clo_addOption(list, "ifile", CLO_TYPE_IFILE, NULL, "input L2 file name");
    clo_addOption(list, "ofile", CLO_TYPE_OFILE, "output", "output filename");
    clo_addOption(list, "prod", CLO_TYPE_STRING, "chlor_a", "product name");
    clo_addOption(list, "quality", CLO_TYPE_INT, "999", "highest quality value acceptable");
    clo_addOption(list, "rflag", CLO_TYPE_STRING, "ORIGINAL", "replacement flag");
    clo_addOption(list, "flaguse", CLO_TYPE_STRING, NULL, "Flags used to mask data");
    clo_addOption(list, "chl_flags", CLO_TYPE_STRING, "ATMFAIL,HILT,STRAYLIGHT,CLDICE,LOWLW,CHLWARN,CHLFAIL,NAVWARN,MAXAERITER,NAVFAIL,FILTER,HIGLINT",
            "Flags used to mask data for chl product if flaguse not set");
    clo_addOption(list, "sst_flags", CLO_TYPE_STRING, "SSTFAIL",
            "Flags used to mask data for sst product if flaguse not set");

    clo_addOption(list, "spixl", CLO_TYPE_INT, "1", "start pixel number");
    clo_addOption(list, "epixl", CLO_TYPE_INT, "-1", "end pixel number (-1=the last pixel)");
    clo_addOption(list, "dpixl", CLO_TYPE_INT, "1", "pixel subsampling interval");
    clo_addOption(list, "sline", CLO_TYPE_INT, "1", "start line number");
    clo_addOption(list, "eline", CLO_TYPE_INT, "-1", "end line number (-1=the last line)");
    clo_addOption(list, "dline", CLO_TYPE_INT, "1", "line subsampling interval");

    clo_addOption(list, "apply_pal", CLO_TYPE_BOOL, "no", "apply color palette, false = grayscale");
    clo_addOption(list, "palfile", CLO_TYPE_IFILE, "default", "palette filename.  \"default\" means the\n        palette is chosen using the product table file");
    clo_addOption(list, "palette_dir", CLO_TYPE_IFILE, "$OCDATAROOT/common/palette", "directory\n        containing the palette files");
    clo_addOption(list, "product_table", CLO_TYPE_IFILE, "$OCDATAROOT/common/l2brsgen_product_table.dat", "product table");

    clo_addOption(list, "datamin", CLO_TYPE_FLOAT, "0.0", "minimum value for data scaling\n        (default see product_table)");
    clo_addOption(list, "datamax", CLO_TYPE_FLOAT, "0.0", "maximum value for data scaling\n        (default see product_table)");
    clo_addOption(list, "stype", CLO_TYPE_INT, "0", "scaling type (default see product_table)\n        1: LINEAR\n        2: LOG");

    strcpy(tmpStr, "format of the output file\n");
    strcat(tmpStr, "        hdf4: (1) HDF browse file\n");
    strcat(tmpStr, "        png:  (5) PNG color or grayscale image file\n");
    strcat(tmpStr, "        ppm:  (7) PPM color or PGM grayscale image file\n");
    clo_addOption(list, "oformat", CLO_TYPE_STRING, "HDF4", tmpStr);

    return 0;
}

//-----------------------------------------------------------------------

/**
   Read the command line option and all of the default parameter files.

   This is the order for loading the options:
    - load the main program defaults file
    - load sensor specific defaults file
    - load the command line (including specified par files)
    - re-load the command line disabling file descending so command
       line arguments will over ride
    - opens l2_str
    - loads meta_l2

 */
int l2brsgen_read_options(clo_optionList_t* list, int argc, char* argv[],
        l2_prod *l2_str, meta_l2Type *meta_l2) {
    char *dataRoot;
    char tmpStr[FILENAME_MAX];
    int sensorId;

    assert(list);

    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        fprintf(stderr, "-E- OCDATAROOT environment variable is not defined.\n");
        return (-1);
    }

    // disable the dump option until we have read all of the files
    clo_setEnableDumpOptions(0);
    clo_readArgs(list, argc, argv);

    // open ifile to get sensor info
    parse_file_name(clo_getString(list, "ifile"), tmpStr);
    openL2(tmpStr, 0x0, l2_str);
    readL2meta(meta_l2, l2_str->fileindex);
    sscanf(meta_l2->title, "%s", tmpStr);
    sensorId = sensorName2SensorId(tmpStr);
    if (sensorId == -1) {
        fprintf(stderr, "-E- Could not find sensor %s in sensorName list.\n", tmpStr);
        return (-1);
    }

    // load program defaults
    sprintf(tmpStr, "%s/common/l2brsgen_defaults.par", dataRoot);
    if (want_verbose) {
        fprintf(stderr, "Loading default parameters from %s\n", tmpStr);
    }
    clo_readFile(list, tmpStr);

    // load the sensor specific defaults file
    if (sensorId != -1) {
        sprintf(tmpStr, "%s/%s/l2brsgen_defaults.par", dataRoot, sensorId2SensorDir(sensorId));

        // see if the file exists
        if (!access(tmpStr, R_OK)) {
            if (want_verbose)
                printf("Loading default parameters for %s from %s\n",
                    sensorId2SensorName(sensorId), tmpStr);
            clo_readFile(list, tmpStr);
        }
    }

    clo_setEnableDumpOptions(1);
    clo_readArgs(list, argc, argv);

    return 0;
}
