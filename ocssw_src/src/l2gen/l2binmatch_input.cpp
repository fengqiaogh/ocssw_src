#include "l2binmatch_input.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <genutils.h>
#include <clo.h>
#include "version.h"

// need a place to store the default product lists
static char default_l2prod[L1_PRODSTRLEN];

static char default_flaguse[1024];

int l2binmatch_input_init() {

    input = (instr *) allocateMemory(sizeof (instr), "input structure");

    //    input->ifile[0] = '\0';
    //    input->ofile[0] = '\0';
    strcpy(default_flaguse, "ATMFAIL,LAND,HIGLINT,HILT,HISATZEN,STRAYLIGHT,CLDICE,COCCOLITH,LOWLW,CHLFAIL,NAVWARN,ABSAER,MAXAERITER,ATMWARN,HISOLZEN,NAVFAIL");
    strcpy(input->flaguse, default_flaguse);

    //    input->l2prod [0] = '\0';

    input->vcal_depth = -1000;
    input->vcal_min_nbin = 4;
    input->vcal_min_nscene = 3;
    input->subsamp = 1;
    input->band_shift_opt = 0;
    input->deflate = 0;
    input->demfile[0] = '\0';

    l1_input_init();

    return 0;
}

//-----------------------------------------------------------------------

/** add all of the accepted command line options to list */
int l2binmatch_init_options(clo_optionList_t* list) {
    char tmpStr[2048];

    sprintf(tmpStr, "l2binmatch %d.%d.%d-%s (%s %s)", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH, GITSHA, __DATE__, __TIME__);
    clo_setVersion(tmpStr);

    sprintf(tmpStr, "Usage: l2binmatch argument-list\n\n");

    strcat(tmpStr, "  This program matches pixels from a L2 file with bins from a L3 file\n");

    strcat(tmpStr, "  The argument-list is a set of keyword=value pairs. The arguments can\n");
    strcat(tmpStr, "  be specified on the commandline, or put into a parameter file, or the\n");
    strcat(tmpStr, "  two methods can be used together, with commandline over-riding.\n\n");
    strcat(tmpStr, "The list of valid keywords follows:\n");
    clo_setHelpStr(tmpStr);

    clo_addOption(list, "ifile", CLO_TYPE_IFILE, NULL, "input L2 file name");
    clo_addOption(list, "ofile", CLO_TYPE_OFILE, NULL, "output matchup file");
    clo_addOption(list, "l2prod", CLO_TYPE_STRING, NULL, "products to be included in ofile");
    clo_addOption(list, "flaguse", CLO_TYPE_STRING, default_flaguse, "flags to mask");

    strcpy(tmpStr, "output file format\n");
    strcat(tmpStr, "        netcdf4: output a netCDF version 4 file\n");
    strcat(tmpStr, "        hdf4:    output a HDF version 4 file");
    clo_addOption(list, "oformat", CLO_TYPE_STRING, "netCDF4", tmpStr);

    clo_addOption(list, "deflate", CLO_TYPE_INT, "0", "deflation level");
    clo_addOption(list, "tgtfile", CLO_TYPE_IFILE, NULL, "L3 bin target file");
    clo_addOption(list, "spixl", CLO_TYPE_INT, "1", "start pixel number");
    clo_addOption(list, "epixl", CLO_TYPE_INT, "-1", "end pixel number (-1=the last pixel)");
    clo_addOption(list, "dpixl", CLO_TYPE_INT, "1", "pixel sub-sampling interval");
    clo_addOption(list, "sline", CLO_TYPE_INT, "1", "start line number");
    clo_addOption(list, "eline", CLO_TYPE_INT, "-1", "end line number (-1=the last line)");
    clo_addOption(list, "dline", CLO_TYPE_INT, "1", "line sub-sampling interval");
    clo_addOption(list, "subsamp", CLO_TYPE_INT, "1", "valid pixel sub-sampling");

    strcpy(tmpStr, "bandshifting option \n");
    strcat(tmpStr, "        1: apply bio-optical bandshift\n");
    strcat(tmpStr, "        0: linear interpolation");
    clo_addOption(list, "band_shift_opt", CLO_TYPE_INT, "0", tmpStr);

    strcpy(tmpStr, "depth to use to exclude data from target file\n");
    strcat(tmpStr, "   e.g. -1000 excludes depths less than 1000m\n");
    clo_addOption(list, "vcal_depth", CLO_TYPE_FLOAT, "-1000.0", tmpStr);

    strcpy(tmpStr, "minimum # of samples in a bin for acceptance");
    clo_addOption(list, "vcal_min_nbin", CLO_TYPE_INT, "4", tmpStr);

    strcpy(tmpStr, "minimum # of scenes in a bin for acceptance");
    clo_addOption(list, "vcal_min_nscene", CLO_TYPE_INT, "3", tmpStr);

    clo_addOption(list, "demfile", CLO_TYPE_IFILE, "$OCDATAROOT/common/ETOPO1_ocssw.nc", "global elevation netCDF file");

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
int l2binmatch_read_options(clo_optionList_t* list, int argc, char* argv[]) {
    char *dataRoot;
    char tmpStr[FILENAME_MAX];
    clo_option_t* option;
    assert(list);

    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        fprintf(stderr, "-E- OCDATAROOT environment variable is not defined.\n");
        return (-1);
    }

    // load program defaults
    sprintf(tmpStr, "%s/common/l2binmatch_defaults.par", dataRoot);
    if (want_verbose) {
        fprintf(stderr, "Loading default parameters from %s\n", tmpStr);
    }
    clo_readFile(list, tmpStr);
    // set the default l2prod lists before the command line or par file
    // is loaded
    option = clo_findOption(list, "l2prod");
    if (option && clo_isOptionSet(option))
        strcpy(default_l2prod, option->valStr);

    clo_readArgs(list, argc, argv);

    return 0;
}

//-----------------------------------------------------------------------------

int l2binmatch_load_input(clo_optionList_t *list) {
    char tmp_file[FILENAME_MAX];
    const char* tmpStr;
    char *strVal;
    clo_option_t *option;
    int numOptions;
    int optionId;
    char keyword[FILENAME_MAX];
    int count;
    char **strArray;
    int i;

    // first load up the default_l2prod
    if (default_l2prod[0]) {
        strcpy(input->def_l2prod[0], default_l2prod);
    }

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
        else if (strcmp(keyword, "par") == 0)
            ;
        else if (strcmp(keyword, "ifile") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->ifile[0], tmp_file);

        } else if (strcmp(keyword, "ofile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->ofile[0], tmp_file);
            }

        } else if (strcmp(keyword, "oformat") == 0) {
            strVal = clo_getOptionString(option);
            tmpStr = getFileFormatName(strVal);
            if (tmpStr == NULL) {
                printf("-E- l2binmatch_load_input: oformat=%s is not a recognized file format\n",
                        strVal);
                return -1;
            }
            strcpy(input->oformat, tmpStr);

        } else if (strcmp(keyword, "flaguse") == 0) {
            strArray = clo_getOptionStrings(option, &count);
            input->flaguse[0] = '\0';
            for (i = 0; i < count; i++) {
                if (i != 0)
                    strcat(input->flaguse, ",");
                strcat(input->flaguse, strArray[i]);
            }

        } else if (strcmp(keyword, "l2prod") == 0) {
            strArray = clo_getOptionStrings(option, &count);
            input->l2prod[0][0] = '\0';
            for (i = 0; i < count; i++) {
                if (i != 0)
                    strcat(input->l2prod[0], " ");
                strcat(input->l2prod[0], strArray[i]);
            }

        } else if (strcmp(keyword, "deflate") == 0) {
            input->deflate = clo_getOptionInt(option);

        } else if (strcmp(keyword, "spixl") == 0) {
            l1_input->spixl = clo_getOptionInt(option);

        } else if (strcmp(keyword, "epixl") == 0) {
            l1_input->epixl = clo_getOptionInt(option);

        } else if (strcmp(keyword, "dpixl") == 0) {
            l1_input->dpixl = clo_getOptionInt(option);

        } else if (strcmp(keyword, "sline") == 0) {
            l1_input->sline = clo_getOptionInt(option);

        } else if (strcmp(keyword, "eline") == 0) {
            l1_input->eline = clo_getOptionInt(option);

        } else if (strcmp(keyword, "dline") == 0) {
            l1_input->dline = clo_getOptionInt(option);

        } else if (strcmp(keyword, "subsamp") == 0) {
            input->subsamp = clo_getOptionInt(option);

        } else if (strcmp(keyword, "band_shift_opt") == 0) {
            input->band_shift_opt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "tgtfile") == 0) {
            if (clo_isOptionSet(option)) {
                strVal = clo_getOptionString(option);
                parse_file_name(strVal, tmp_file);
                strcpy(input->tgtfile, tmp_file);
            }

        } else if (strcmp(keyword, "vcal_depth") == 0) {
            input->vcal_depth = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "vcal_min_nbin") == 0) {
            input->vcal_min_nbin = clo_getOptionInt(option);

        } else if (strcmp(keyword, "vcal_min_nscene") == 0) {
            input->vcal_min_nscene = clo_getOptionInt(option);

        } else if (strcmp(keyword, "demfile") == 0) {
            strVal = clo_getOptionString(option);
            parse_file_name(strVal, tmp_file);
            strcpy(input->demfile, tmp_file);

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

int l2binmatch_input(int argc, char **argv, clo_optionList_t* list) {
    char str_buf[4096];


    /* initialize the option list with descriptions and default values */
    l2binmatch_init_options(list);

    /*                                                                  */
    /* Set input values to defaults                                     */
    /*                                                                  */
    if (l2binmatch_input_init() != 0) {
        fprintf(stderr, "-E- %s: Error initializing input structure.\n", __FILE__);
        return (-1);
    }

    /* read the command line options into list */
    if (l2binmatch_read_options(list, argc, argv) != 0) {
        fprintf(stderr, "-E- %s: Error reading program options.\n", __FILE__);
        return (-1);
    }

    /* load options from list into input structure */
    if (l2binmatch_load_input(list) != 0) {
        fprintf(stderr, "-E- %s: Error loading options into input structure.\n", __FILE__);
        return (-1);
    }

    /*                                                                  */
    /* Build string of parameters for metadata                          */
    /*                                                                  */

    strcat(l1_input->input_parms, "\n");
    sprintf(str_buf, "ifile = %s ", input->ifile[0]);
    strcat(l1_input->input_parms, str_buf);
    strcat(l1_input->input_parms, "\n");

    sprintf(str_buf, "ofile = %s", input->ofile[0]);
    strcat(l1_input->input_parms, str_buf);
    strcat(l1_input->input_parms, "\n");

    sprintf(str_buf, "l2prod = %s", input->l2prod[0]);
    strcat(l1_input->input_parms, str_buf);
    strcat(l1_input->input_parms, "\n");

    sprintf(str_buf, "oformat = %s", input->oformat);
    strcat(l1_input->input_parms, str_buf);
    strcat(l1_input->input_parms, "\n");

    sprintf(str_buf, "deflate = %5d", input->deflate);
    strcat(l1_input->input_parms, str_buf);
    strcat(l1_input->input_parms, "\n");

    sprintf(str_buf, "spixl = %5d", l1_input->spixl);
    strcat(l1_input->input_parms, str_buf);
    strcat(l1_input->input_parms, "\n");

    sprintf(str_buf, "epixl = %5d", l1_input->epixl);
    strcat(l1_input->input_parms, str_buf);
    strcat(l1_input->input_parms, "\n");

    sprintf(str_buf, "dpixl = %5d", l1_input->dpixl);
    strcat(l1_input->input_parms, str_buf);
    strcat(l1_input->input_parms, "\n");

    sprintf(str_buf, "sline = %5d", l1_input->sline);
    strcat(l1_input->input_parms, str_buf);
    strcat(l1_input->input_parms, "\n");

    sprintf(str_buf, "eline = %5d", l1_input->eline);
    strcat(l1_input->input_parms, str_buf);
    strcat(l1_input->input_parms, "\n");

    sprintf(str_buf, "dline = %5d", l1_input->dline);
    strcat(l1_input->input_parms, str_buf);
    strcat(l1_input->input_parms, "\n");

    sprintf(str_buf, "subsamp = %5d", input->subsamp);
    strcat(l1_input->input_parms, str_buf);
    strcat(l1_input->input_parms, "\n");

    sprintf(str_buf, "band_shift_opt = %3d", input->band_shift_opt);
    strcat(l1_input->input_parms, str_buf);
    strcat(l1_input->input_parms, "\n");

    sprintf(str_buf, "vcal_depth = %8.4f", input->vcal_depth);
    strcat(l1_input->input_parms, str_buf);
    strcat(l1_input->input_parms, "\n");

    sprintf(str_buf, "vcal_min_nbin = %d", input->vcal_min_nbin);
    strcat(l1_input->input_parms, str_buf);
    strcat(l1_input->input_parms, "\n");

    sprintf(str_buf, "vcal_min_nscene = %d", input->vcal_min_nscene);
    strcat(l1_input->input_parms, str_buf);
    strcat(l1_input->input_parms, "\n");

    sprintf(str_buf, "demfile = %s", input->demfile);
    strcat(l1_input->input_parms, str_buf);
    strcat(l1_input->input_parms, "\n");

    return 0;
}
