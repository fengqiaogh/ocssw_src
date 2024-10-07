/*
 * l3bindump_input.c
 *
 *  Created on: Mar 14, 2014
 *      Author: Sean Bailey
 */

#include "l3bindump.h"
#include <string.h>
#include <stdlib.h>

/** add all of the accepted command line options to list */
int l3bindump_init_options(clo_optionList_t* list) {
    char tmpStr[2048];

    clo_setSelectOptionKeys(NULL);

    char softwareVersion[200];
    sprintf(softwareVersion, "%d.%d.%d-%s", VERSION_MAJOR, VERSION_MINOR,
            VERSION_PATCH, GITSHA);

    clo_setVersion2("l3bindump", softwareVersion);
    //    clo_addXmlProgramMetadata("progressRegex", "Processing scan .+?\\((\\d+) of (\\d+)\\)");

    sprintf(tmpStr, "Usage: l3bindump argument-list\n\n");
    strcat(tmpStr, "  The argument-list is a set of keyword=value pairs. The arguments can\n");
    strcat(tmpStr, "  be specified on the command line, or put into a parameter file, or the\n");
    strcat(tmpStr, "  two methods can be used together, with command line overriding.\n\n");
    strcat(tmpStr, "  return value: 0=OK, 1=error, 110=requested bin(s) not found\n");
    strcat(tmpStr, "  file data.\n\n");
    strcat(tmpStr, "  There are 3 use cases:\n");
    strcat(tmpStr, "     1) dump the bin requested by bin number\n");
    strcat(tmpStr, "         use: bin_number=<the number>\n");
    strcat(tmpStr, "     2) region defined by lat, lon, and radius (in km)\n");
    strcat(tmpStr, "         use: lat=<latitude> lon=<longitude> radius=<radius in km>\n");
    strcat(tmpStr, "     3) region defined by north, south, east, west\n");
    strcat(tmpStr, "         use: north=<N> south=<S> east=<E> west=<W>\n");

    strcat(tmpStr, "The list of valid keywords follows:\n");
    clo_setHelpStr(tmpStr);

    clo_addOption(list, "ifile", CLO_TYPE_IFILE, NULL, "input L3 bin file name");
    clo_addOption(list, "ofile", CLO_TYPE_OFILE, NULL, "output file name");

    strcpy(tmpStr, "output file format\n");
    strcat(tmpStr, "        txt:  plain text columnar format\n");
    strcat(tmpStr, "    seabass:  SeaBASS format");
    clo_addOption(list, "oformat", CLO_TYPE_STRING, "txt", tmpStr);

    clo_addOption(list, "l3bprod", CLO_TYPE_STRING, "Unspecified", "binned product to extract");
    clo_addOption(list, "bin_number", CLO_TYPE_INT64, "-1", "bin number");
    clo_addOption(list, "north", CLO_TYPE_FLOAT, "-999", "north boundary");
    clo_addOption(list, "south", CLO_TYPE_FLOAT, "-999", "south boundary");
    clo_addOption(list, "east", CLO_TYPE_FLOAT, "-999", "east boundary");
    clo_addOption(list, "west", CLO_TYPE_FLOAT, "-999", "west boundary");
    clo_addOption(list, "lat", CLO_TYPE_FLOAT, "-999", "latitude");
    clo_addOption(list, "lon", CLO_TYPE_FLOAT, "-999", "longitude");
    clo_addOption(list, "radius", CLO_TYPE_FLOAT, "-999", "radius in km");
    clo_addOption(list, "verbose", CLO_TYPE_BOOL, "0", "verbose output");

    return 0;
}

int l3bindump_load_input(clo_optionList_t *list, instr *input) {
    const char* tmpStr;
    char tmp_file[FILENAME_MAX];
    char *strVal;
    clo_option_t *option;
    int numOptions;
    int optionId;
    char keyword[FILENAME_MAX];
    int count;
    char **strArray;
    int i;
    numOptions = clo_getNumOptions(list);
    for (optionId = 0; optionId < numOptions; optionId++) {
        option = clo_getOption(list, optionId);

        // ignore options of type CLO_TYPE_HELP
        if (option->dataType == CLO_TYPE_HELP)
            continue;

        strcpy(keyword, option->key);

        /* change keyword to lower case */
        strVal = keyword;
        while (*strVal != '\0') {
            *strVal = tolower(*strVal);
            strVal++;
        }

        if (strcmp(keyword, "help") == 0)
            ;
        else if (strcmp(keyword, "version") == 0)
            ;
        else if (strncmp(keyword, "dump_options", 12) == 0)
            ;
        else if (strcmp(keyword, "par") == 0)
            ;
        else if (strcmp(keyword, "ifile") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->ifile, tmp_file);

        } else if (strcmp(keyword, "ofile") == 0) {
            if (clo_isOptionSet(option)) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->ofile, tmp_file);
            }
        } else if (strcmp(keyword, "oformat") == 0) {
            strVal = clo_getOptionString(option);
            tmpStr = getFileFormatName(strVal);
            if (tmpStr == NULL) {
                printf("-E- l2gen_load_input: oformat=%s is not a recognized file format\n",
                        strVal);
                return -1;
            }
            strcpy(input->oformat, tmpStr);

        } else if (strcmp(keyword, "l3bprod") == 0) {
            strArray = clo_getOptionStrings(option, &count);
            input->l3bprod[0] = '\0';
            for (i = 0; i < count; i++) {
                if (i != 0)
                    strcat(input->l3bprod, " ");
                strcat(input->l3bprod, strArray[i]);
            }
        } else if (strcmp(keyword, "west") == 0) {
            input->west = clo_getOptionFloat(option);
        } else if (strcmp(keyword, "east") == 0) {
            input->east = clo_getOptionFloat(option);
        } else if (strcmp(keyword, "north") == 0) {
            input->north = clo_getOptionFloat(option);
        } else if (strcmp(keyword, "south") == 0) {
            input->south = clo_getOptionFloat(option);
        } else if (strcmp(keyword, "lat") == 0) {
            input->lat = clo_getOptionFloat(option);
        } else if (strcmp(keyword, "lon") == 0) {
            input->lon = clo_getOptionFloat(option);
        } else if (strcmp(keyword, "bin_number") == 0) {
            input->bin_number = clo_getOptionInt64(option);
        } else if (strcmp(keyword, "radius") == 0) {
            input->radius = clo_getOptionFloat(option);
        } else if (strcmp(keyword, "verbose") == 0) {
            input->verbose = clo_getOptionBool(option);
        } else {
            printf("-E- Invalid argument \"%s\"\n", keyword);
            clo_dumpOption(option);
            exit(1);
        }

    } // for optionIDs

    /*
     * check for valid use cases
     */
    if (input->bin_number >= 0 && (input->radius != -999 ||
            input->lat != -999 || input->lon != -999 ||
            input->north != -999 || input->south != -999 ||
            input->west != -999 || input->east != -999)) {
        printf("-E- Invalid argument set: bin_number cannot be provided if lon/lat/radius or NSWE are provided \n");
        clo_dumpOption(option);
        exit(1);
    }
    if (input->radius != -999 && (input->lat == -999 || input->lon == -999)) {
        printf("-E- Invalid argument set: radius reqiures lon/lat to be provided\n");
        clo_dumpOption(option);
        exit(1);
    }
    if (input->lat != -999 && input->lon == -999) {
        printf("-E- Invalid argument set: if lat is set, lon must be provided\n");
        clo_dumpOption(option);
        exit(1);
    }
    if (input->lat == -999 && input->lon != -999) {
        printf("-E- Invalid argument set: if lon is set, lat must be provided\n");
        clo_dumpOption(option);
        exit(1);
    }
    if ((input->radius != -999 || input->lat != -999 || input->lon != -999) &&
            (input->north != -999 || input->south != -999 ||
            input->west != -999 || input->east != -999)) {
        printf("-E- Choose either lon/lat/radius or NSWE\n");
        clo_dumpOption(option);
        exit(1);
    }

    return 0;
}


//-----------------------------------------------------------------------

/*
   Read the command line option and all of the default parameter files.

   This is the order for loading the options:
    - read the command line to get the ifile and suite options
    - load the main program defaults file
    - load the command line (including specified par files)
    - re-load the command line so they take precedence

 */
int l3bindump_read_options(clo_optionList_t* list, int argc, char* argv[]) {
    char *dataRoot;
    char tmpStr[FILENAME_MAX];
    char progName[] = "l3bindump";

    assert(list);

    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        printf("-E- OCDATAROOT environment variable is not defined.\n");
        return (-1);
    }

    // disable the dump option until we have read all of the files
    clo_setEnableDumpOptions(0);
    clo_readArgs(list, argc, argv);

    // load program defaults
    sprintf(tmpStr, "%s/common/%s_defaults.par", dataRoot, progName);
    if (want_verbose)
        printf("Loading default parameters from %s\n", tmpStr);
    clo_readFile(list, tmpStr);

    // re-load the command line and par file
    if (want_verbose)
        printf("Loading command line parameters\n\n");
    // enable the dump option the last time through
    clo_setEnableDumpOptions(1);
    clo_readArgs(list, argc, argv);

    return 0;
}

void l3bindump_input_init(instr *input) {
    /*                                                                  */
    /* Set input values to defaults                                     */
    /*                                                                  */

    input->ifile[0] = '\0';
    input->ofile [0] = '\0';
    input->oformat[0] = '\0';
    input->l3bprod[0] = '\0';

    input->west = -999;
    input->east = -999;
    input->north = -999;
    input->south = -999;
    input->lat = -999;
    input->lon = -999;
    input->bin_number = -1;
    input->radius = -999;
    input->verbose = 0;

    return;
}

int l3bindump_usage(char *prog) {
    clo_optionList_t* list;

    list = clo_createList();
    l3bindump_init_options(list);
    clo_printUsage(list);

    return 0;
}
