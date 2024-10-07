#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <genutils.h>
#include <clo.h>
#include "l2mapgen.h"
#include "l2mapgen_input.h"

int input_init(instr *input_str) {
    input_str->ifile[0] = '\0';
    input_str->ofile[0] = '\0';
    input_str->palfile[0] = '\0';
    input_str->palette_dir[0] = '\0';
    input_str->product_table[0] = '\0';
    input_str->flaguse[0] = '\0';
    input_str->parms[0] = '\0';

    input_str->prod [0] = '\0';
    input_str->stype = 0;
    input_str->datamin = 0.0;
    input_str->datamax = 0.0;

    input_str->west = 0.0;
    input_str->east = 0.0;
    input_str->south = 0.0;
    input_str->north = 0.0;
    input_str->width = 800;
    input_str->threshold = 5.;
    input_str->mask = 0;
    input_str->quality = 2;
    input_str->apply_pal = 0;
    input_str->outmode = 1;

    return 0;
}

//-----------------------------------------------------------------------

/** add all of the accepted command line options to list */
int l2mapgen_init_options(clo_optionList_t* list) {
    char tmpStr[2048];

    sprintf(tmpStr, "l2mapgen %s (%s %s)", VERSION, __DATE__, __TIME__);
    clo_setVersion(tmpStr);

    sprintf(tmpStr, "Usage: l2mapgen argument-list\n\n");

    strcat(tmpStr, "  This program takes a product from a L2 file, maps it using a Plate\n");
    strcat(tmpStr, "  Carree cylindrical projection, and produces a gray scale PGM or\n");
    strcat(tmpStr, "  color PPM file.\n\n");

    strcat(tmpStr, "  The argument-list is a set of keyword=value pairs. The arguments can\n");
    strcat(tmpStr, "  be specified on the commandline, or put into a parameter file, or the\n");
    strcat(tmpStr, "  two methods can be used together, with commandline over-riding.\n\n");
    strcat(tmpStr, "The list of valid keywords follows:\n");
    clo_setHelpStr(tmpStr);

    clo_addOption(list, "ifile", CLO_TYPE_IFILE, NULL, "input L2 file name or file with a list of files names");
    clo_addOption(list, "ofile", CLO_TYPE_OFILE, NULL, "output map filename (NULL=STDOUT)");
    clo_addOption(list, "prod", CLO_TYPE_STRING, NULL, "product name");
    clo_addOption(list, "apply_pal", CLO_TYPE_BOOL, "false", "apply color palette, false = grayscale");
    clo_addOption(list, "palfile", CLO_TYPE_IFILE, "default", "palette filename");
    clo_addOption(list, "palette_dir", CLO_TYPE_IFILE, "$OCDATAROOT/common/palette", "palette directory");
    clo_addOption(list, "product_table", CLO_TYPE_IFILE, "$OCDATAROOT/common/smigen_product_table.dat", "product table");
    clo_addOption(list, "flaguse", CLO_TYPE_STRING, NULL, "flags to be masked");
    clo_addOption(list, "quality", CLO_TYPE_INT, "2", "minimum allowable quality level for SST.  Valid\n        only for SST  and only if qual_sst or qual_sst4 SDS exist");
    clo_addOption(list, "mask", CLO_TYPE_BOOL, "no", "apply mask to land, cloud and glint (see below)");
    clo_addOption(list, "datamin", CLO_TYPE_FLOAT, "0.0", "minimum value for data scaling\n        (default see SMI product table)");
    clo_addOption(list, "datamax", CLO_TYPE_FLOAT, "0.0", "maximum value for data scaling\n        (default see SMI product table)");
    clo_addOption(list, "stype", CLO_TYPE_INT, "0", "scaling type (default see SMI product table)\n        1: LINEAR\n        2: LOG");
    clo_addOption(list, "east", CLO_TYPE_FLOAT, "0.0", "Map East longitude\n        (default=scene(s) Easternmost Longitude)");
    clo_addOption(list, "west", CLO_TYPE_FLOAT, "0.0", "Map West longitude\n        (default=scene(s) Westernmost Longitude)");
    clo_addOption(list, "north", CLO_TYPE_FLOAT, "0.0", "Map North latitude\n        (default=scene(s) Northernmost Longitude)");
    clo_addOption(list, "south", CLO_TYPE_FLOAT, "0.0", "Map South latitude\n        (default=scene(s) Southernmost Longitude)");
    clo_addOption(list, "width", CLO_TYPE_INT, "800", "width of the output image");
    clo_addOption(list, "threshold", CLO_TYPE_FLOAT, "5", "minimum percentage of the area of interest\n        that must receive valid pixel data before an image is generated");
    strcpy(tmpStr, "format of the output file\n");
    strcat(tmpStr, "        ppm: PPM or PGM image file (alias 1)\n");
    strcat(tmpStr, "        png: PNG color or grayscale image file (alias 2)\n");
    strcat(tmpStr, "        tiff: TIFF color or grayscale geo tiff image file (alias 3)\n");
    clo_addOption(list, "outmode", CLO_TYPE_STRING, "ppm", tmpStr);

    strcpy(tmpStr, "\n   If the \"mask\" option is set, the output PGM image will be masked for\n");
    strcat(tmpStr, "   flags defined in the flaguse parameter. The \"no data\" pixel value will\n");
    strcat(tmpStr, "   change from 0 to 255, and pixel values 252, 253, and 254 will represent the\n");
    strcat(tmpStr, "   sunglint, land, and all other (e.g. clouds/ice,hilt,atmfail,navfail,chlfail)\n");
    strcat(tmpStr, "   masks, respectively. NOTE: sunglint is NOT masked by default, but if it is\n");
    strcat(tmpStr, "   added to the flaguse parameter, it will be distinguished in the masking as\n");
    strcat(tmpStr, "   medium gray.  If a palette is applied and the mask option is set, the\n");
    strcat(tmpStr, "   palette values will be modified:\n");
    strcat(tmpStr, "                  Value   R       G       B\n");
    strcat(tmpStr, "                  252     128     128     128\n");
    strcat(tmpStr, "                  253     160     82      45\n");
    strcat(tmpStr, "                  254     255     255     255\n");
    strcat(tmpStr, "                  255     0       0       0\n\n");
    strcat(tmpStr, "   By default, this program sends its results to standard output as a\n");
    strcat(tmpStr, "   PGM-formatted binary data stream.  Save it to a file via \">\" or pipe it\n");
    strcat(tmpStr, "   to your favorite image display program.  The output image is rendered in\n");
    strcat(tmpStr, "   a Plate Carree projection.");
    clo_addOption(list, "help1", CLO_TYPE_HELP, NULL, tmpStr);

    return 0;
}

//-----------------------------------------------------------------------

/* 
   Read the command line option and all of the default parameter files.

   This is the order for loading the options:
    - load the main program defaults file
    - load the command line (including specified par files)
    - re-load the command line disabling file decending so command
       line arguments will over ride

 */
int l2mapgen_read_options(clo_optionList_t* list, int argc, char* argv[]) {
    char *dataRoot;
    char tmpStr[FILENAME_MAX];

    assert(list);

    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        fprintf(stderr, "-E- OCDATAROOT environment variable is not defined.\n");
        return (-1);
    }

    // load program defaults
    sprintf(tmpStr, "%s/common/l2mapgen_defaults.par", dataRoot);
    if (want_verbose) {
        fprintf(stderr, "Loading default parameters from %s\n", tmpStr);
    }
    clo_readFile(list, tmpStr);
    clo_readArgs(list, argc, argv);

    return 0;
}

//-----------------------------------------------------------------------------

int l2mapgen_load_input(clo_optionList_t *list, instr *input) {
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
        else if (strcmp(keyword, "dump_options") == 0)
            ;
        else if (strcmp(keyword, "dump_options_paramfile") == 0)
            ;
        else if (strcmp(keyword, "dump_options_xmlfile") == 0)
            ;
        else if (strcmp(keyword, "help1") == 0)
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

        } else if (strcmp(keyword, "palfile") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->palfile, tmp_file);

        } else if (strcmp(keyword, "palette_dir") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->palette_dir, tmp_file);

        } else if (strcmp(keyword, "product_table") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->product_table, tmp_file);

        } else if (strcmp(keyword, "flaguse") == 0) {
            if (clo_isOptionSet(option)) {
                strArray = clo_getOptionStrings(option, &count);
                input->flaguse[0] = '\0';
                for (i = 0; i < count; i++) {
                    if (strcasecmp("default", strArray[i]) == 0) {
                        if (input->flaguse[0] != '\0')
                            strcat(input->flaguse, ",");
                        strcat(input->flaguse, DEF_FLAG);
                    } else {
                        if (input->flaguse[0] != '\0')
                            strcat(input->flaguse, ",");
                        strcat(input->flaguse, strArray[i]);
                    }

                } // for count
            } // if option set

        } else if (strcmp(keyword, "prod") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->prod, tmp_file);

        } else if (strcmp(keyword, "stype") == 0) {
            input->stype = clo_getOptionInt(option);

        } else if (strcmp(keyword, "width") == 0) {
            input->width = clo_getOptionInt(option);

        } else if (strcmp(keyword, "threshold") == 0) {
            input->threshold = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "apply_pal") == 0) {
            input->apply_pal = clo_getOptionBool(option);

        } else if (strcmp(keyword, "mask") == 0) {
            input->mask = clo_getOptionBool(option);

        } else if (strcmp(keyword, "quality") == 0) {
            input->quality = clo_getOptionInt(option);

        } else if (strcmp(keyword, "datamin") == 0) {
            input->datamin = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "datamax") == 0) {
            input->datamax = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "west") == 0) {
            input->west = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "east") == 0) {
            input->east = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "south") == 0) {
            input->south = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "north") == 0) {
            input->north = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "outmode") == 0) {
            char *str = clo_getOptionString(option);
            if (strcmp("1", str) == 0) {
                input->outmode = 1;
            } else if (strcasecmp("ppm", str) == 0) {
                input->outmode = 1;
            } else if (strcmp("2", str) == 0) {
                input->outmode = 2;
            } else if (strcasecmp("png", str) == 0) {
                input->outmode = 2;
            } else if (strcmp("3", str) == 0) {
                input->outmode = 3;
            } else if (strcasecmp("tiff", str) == 0) {
                input->outmode = 3;
            } else {
                fprintf(stderr, "Invalid value for outmode \"%s\"\n", str);
                return -1;
            }

        } else {
            fprintf(stderr, "Invalid argument \"%s\"\n", keyword);
            return -1;
        }
    } // for optionId

    return 0;
}

/*-----------------------------------------------------------------------------
    Function:  msmapl2_input

    Returns:   int (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        Convert the arguments from the command line into a structure input
        variable.

    Parameters: (in calling order)
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       argc          I   number of arguments
        char      **argv        I   list of arguments
        instr     input         O   structure variable for inputs

----------------------------------------------------------------------------*/

int l2mapgen_input(int argc, char **argv, clo_optionList_t* list, instr* input) {
    char str_buf[4096];


    /* initialize the option list with descriptions and default values */
    l2mapgen_init_options(list);

    /*                                                                  */
    /* Set input values to defaults                                     */
    /*                                                                  */
    if (input_init(input) != 0) {
        fprintf(stderr, "-E- %s: Error initializing input structure.\n", __FILE__);
        return (-1);
    }

    /* read the command line options into list */
    if (l2mapgen_read_options(list, argc, argv) != 0) {
        fprintf(stderr, "-E- %s: Error reading program options.\n", __FILE__);
        return (-1);
    }

    /* load options from list into input structure */
    if (l2mapgen_load_input(list, input) != 0) {
        fprintf(stderr, "-E- %s: Error loading options into input structure.\n", __FILE__);
        return (-1);
    }

    /*                                                                  */
    /* Build string of parameters for metadata                          */
    /*                                                                  */
    sprintf(str_buf, "IFILE=%s|", input->ifile);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "OFILE=%s|", input->ofile);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "PALFILE=%s|", input->palfile);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "PALETTE_DIR=%s|", input->palette_dir);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "PRODUCT_TABLE=%s|", input->product_table);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "FLAGUSE=%s|", input->flaguse);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "PROD=%s|", input->prod);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "STYPE=%d|", input->stype);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "DATAMIN=%f|", input->datamin);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "DATAMAX=%f|", input->datamax);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "WEST=%f|", input->west);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "EAST=%f|", input->east);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "SOUTH=%f|", input->south);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "NORTH=%f|", input->north);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "WIDTH=%d|", input->width);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "THRESHOLD=%f|", input->threshold);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "MASK=%d|", input->mask);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "QUALITY=%d|", input->quality);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "APPLY_PAL=%d|", input->apply_pal);
    strcat(input->parms, str_buf);

    return 0;
}
