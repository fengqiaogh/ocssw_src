#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <mfhdf.h>

#include "genutils.h"
#include "l3bin_input.h"

static char mainProgramName[50];

static char *l3bin_optionKeys[] = {
    "help",
    "version",
    "verbose",
    "dump_options",
    "dump_options_paramfile",
    "dump_options_xmlfile",
    "par",
    "pversion",
    "ifile",
    "ofile",
    "oformat",
    "merged",
    "latnorth",
    "latsouth",
    "loneast",
    "lonwest",
    "sday",
    "eday",
    "deflate",
    "orbit1",
    "orbit2",
    "median",
    "noext",
    "unit_wgt",
    "composite_scheme",
    "composite_prod",
    "reduce_fac",
    "resolve",
    "prod",
    "doi",
    NULL
};

static char *l3binmerge_optionKeys[] = {
    "help",
    "version",
    "verbose",
    "dump_options",
    "dump_options_paramfile",
    "dump_options_xmlfile",
    "par",
    "pversion",
    "eval",
    "ifile",
    "ofile",
    "latnorth",
    "latsouth",
    "loneast",
    "lonwest",
    "noext",
    "u",
    "prod",
    NULL
};

int input_init(instr *input_str) {
    input_str->infile[0] = '\0';
    input_str->ofile[0] = '\0';
    input_str->pfile[0] = '\0';

    strcpy(input_str->out_parm, "DEFAULT");
    strcpy(input_str->pversion, "Unspecified");

    input_str->syear = 9999;
    input_str->sday = 999;

    input_str->eyear = 9999;
    input_str->eday = 999;

    input_str->sorbit = -1;
    input_str->eorbit = -1;

    input_str->reduce_fac = 1;
    input_str->resolve[0] = '\0';

    input_str->noext = 0;

    input_str->merged[0] = '\0';

    input_str->loneast = +180;
    input_str->lonwest = -180;
    input_str->latnorth = +90;
    input_str->latsouth = -90;

    input_str->verbose = 0;
    input_str->unit_wgt = 0;
    input_str->median = 0;
    input_str->union_bins = 0;

    input_str->deflate = 0;
    input_str->oformat[0] = '\0';

    input_str->composite_prod[0] = '\0';
    input_str->composite_scheme[0] = '\0';
    input_str->doi[0] = '\0';

    return 0;
}

int l3bin_init_options(clo_optionList_t* list, const char* prog, const char* version) {
    char tmpStr[2048];
    clo_option_t* option;

    // set the min program name
    strcpy(mainProgramName, prog);

    if (!strcmp(prog, "l3bin")) {
        clo_setSelectOptionKeys(l3bin_optionKeys);
    } else if (!strcmp(prog, "l3binmerge")) {
        clo_setSelectOptionKeys(l3binmerge_optionKeys);
    }
    sprintf(tmpStr, "%s ifile=input-file ofile=output-file prod=prodlist\n\n", prog);
    strcat(tmpStr, "  The input file is a list of L3 binned files.\n");
    strcat(tmpStr, "  The argument-list is a set of keyword=value pairs. The arguments can\n");
    strcat(tmpStr, "  be specified on the commandline, or put into a parameter file, or the\n");
    strcat(tmpStr, "  two methods can be used together, with commandline over-riding.\n\n");
    strcat(tmpStr, "  return value: 0=OK, 1=error, 110=no pixels binned. \n\n");
    strcat(tmpStr, "The list of valid keywords follows:\n");
    clo_setHelpStr(tmpStr);

    // add the parfile alias for backward compatibility
    clo_addAlias(list, "par", "parfile");

    strcpy(tmpStr, "input file name with list of L3 files");
    option = clo_addOption(list, "ifile", CLO_TYPE_IFILE, NULL, tmpStr);
    clo_addOptionAlias(option, "in");
    clo_addOptionAlias(option, "infile");

    option = clo_addOption(list, "ofile", CLO_TYPE_OFILE, "output", "output bin file name");
    clo_addOptionAlias(option, "out");
    option = clo_addOption(list, "merged", CLO_TYPE_OFILE, NULL, "merged file name");

    strcpy(tmpStr, "output file format\n");
    strcat(tmpStr, "           hdf4:    output a HDF4 file\n");
    strcat(tmpStr, "           netCDF4: output a netCDF4 file\n");
    strcat(tmpStr, "           hdf5:    output a HDF5 file\n");
    clo_addOption(list, "oformat", CLO_TYPE_STRING, "netCDF4", tmpStr);

    strcpy(tmpStr, "set to 1 to suppress generation of\n        external files");
    clo_addOption(list, "noext", CLO_TYPE_BOOL, "off", tmpStr);

    clo_addOption(list, "verbose", CLO_TYPE_BOOL, "off", "Allow more verbose screen messages");
    clo_addOption(list, "latnorth", CLO_TYPE_FLOAT, "+90", "northern most latitude");
    clo_addOption(list, "latsouth", CLO_TYPE_FLOAT, "-90", "southern most latitude");
    clo_addOption(list, "loneast", CLO_TYPE_FLOAT, "+180", "eastern most longitude");
    clo_addOption(list, "lonwest", CLO_TYPE_FLOAT, "-180", "western most longitude");
    clo_addOption(list, "reduce_fac", CLO_TYPE_INT, "1", "scale reduction factor (power of 2)");
    strcpy(tmpStr, "bin resolution, overrides reduce_frac if defined\n");
    strcat(tmpStr, "         H: 0.5km\n");
    strcat(tmpStr, "         Q: 250m\n");
    strcat(tmpStr, "        HQ: 100m\n");
    strcat(tmpStr, "        HH: 50m\n");
    strcat(tmpStr, "         1: 1.1km\n");
    strcat(tmpStr, "         2: 2.3km\n");
    strcat(tmpStr, "         4: 4.6km\n");
    strcat(tmpStr, "         9: 9.2km\n");
    strcat(tmpStr, "        18: 18.5km\n");
    strcat(tmpStr, "        36: 36km\n");
    strcat(tmpStr, "        1D: 1 degree\n");
    strcat(tmpStr, "        HD: 0.5 degree\n");
    strcat(tmpStr, "        QD: 0.25 degree");
    option = clo_addOption(list, "resolve", CLO_TYPE_STRING, NULL, tmpStr);
    clo_addOption(list, "pversion", CLO_TYPE_STRING, "unspecified", "production version");
    clo_addOption(list, "sday", CLO_TYPE_INT, "1970001", "start datadate (YYYYDDD) ");
    clo_addOption(list, "eday", CLO_TYPE_INT, "2038018", "end datadate (YYYYDDD)");
    clo_addOption(list, "deflate", CLO_TYPE_INT, "5", "deflate level");

    clo_addOption(list, "orbit1", CLO_TYPE_INT, "-1", "sorbit");
    clo_addOption(list, "orbit2", CLO_TYPE_INT, "-1", "eorbit");
    clo_addOption(list, "median", CLO_TYPE_INT, "0", "median");
    clo_addOption(list, "unit_wgt", CLO_TYPE_INT, "0", "unit_wgt");
    clo_addOption(list, "composite_scheme", CLO_TYPE_STRING, NULL, "composite scheme (min/max)");
    clo_addOption(list, "composite_prod", CLO_TYPE_STRING, NULL, "composite product fieldname");

    strcpy(tmpStr, "bin products\n        [default=all products in L3 file]\n");
    option = clo_addOption(list, "prod", CLO_TYPE_STRING, "DEFAULT", tmpStr);
    clo_addOptionAlias(option, "out_parm");

    option = clo_addOption(list, "union_bins", CLO_TYPE_BOOL, "off", "Output file contains the union of input bins");
    clo_addOptionAlias(option, "u");

    clo_addOption(list, "doi", CLO_TYPE_STRING, NULL, "Digital Object Identifier (DOI) string");
   
    clo_setVersion(version);
    return 0;
}

int l3bin_load_input(clo_optionList_t* list, instr *input) {
    char *tmp_str;
    char keyword [50];
    char tmp_file[FILENAME_MAX];
    int numOptions, optionId;
    clo_option_t *option;

    numOptions = clo_getNumOptions(list);
    for (optionId = 0; optionId < numOptions; optionId++) {
        option = clo_getOption(list, optionId);

        // ignore options of type CLO_TYPE_HELP
        if (option->dataType == CLO_TYPE_HELP)
            continue;

        strcpy(keyword, option->key);

        /* change keyword to lower case */
        tmp_str = keyword;
        while (*tmp_str != '\0') {
            if (isupper(*tmp_str)) *tmp_str = tolower(*tmp_str);
            tmp_str++;
        }
        if (strcmp(keyword, "help") == 0)
            ;
        else if (strcmp(keyword, "version") == 0)
            ;
        else if (strncmp(keyword, "dump_options", 12) == 0)
            ;
        else if (strncmp(keyword, "par", 3) == 0)
            ;
        else if (strcmp(keyword, "ifile") == 0) {
            if (clo_isOptionSet(option)) {
                parse_file_name(clo_getOptionString(option), tmp_file);
                strcpy(input->infile, tmp_file);
            }
        } else if (strcmp(keyword, "ofile") == 0) {
            if (clo_isOptionSet(option)) {
                parse_file_name(clo_getOptionString(option), tmp_file);
                strcpy(input->ofile, tmp_file);
            }
        } else if (strcmp(keyword, "pfile") == 0) {
            if (clo_isOptionSet(option)) {
                parse_file_name(clo_getOptionString(option), tmp_file);
                strcpy(input->pfile, tmp_file);
            }
        } else if (strcmp(keyword, "pversion") == 0) {
            strcpy(input->pversion, clo_getOptionString(option));

        } else if (strcmp(keyword, "syear") == 0) {
            input->syear = clo_getOptionInt(option);

        } else if (strcmp(keyword, "eyear") == 0) {
            input->eyear = clo_getOptionInt(option);

        } else if (strcmp(keyword, "sday") == 0) {
            input->sday = clo_getOptionInt(option);

        } else if (strcmp(keyword, "eday") == 0) {
            input->eday = clo_getOptionInt(option);

        } else if (strcmp(keyword, "orbit1") == 0) {
            input->sorbit = clo_getOptionInt(option);

        } else if (strcmp(keyword, "orbit2") == 0) {
            input->eorbit = clo_getOptionInt(option);

        } else if (strcmp(keyword, "prod") == 0) {
            strcpy(input->out_parm, ":");
            strcat(input->out_parm, clo_getOptionRawString(option));
            strcat(input->out_parm, ":");

        } else if (strcmp(keyword, "reduce_fac") == 0) {
            input->reduce_fac = clo_getOptionInt(option);

       } else if (strcmp(keyword, "resolve") == 0) {
            if (clo_isOptionSet(option)) {
                strcpy(input->resolve, clo_getOptionString(option));
            }

        } else if (strcmp(keyword, "noext") == 0) {
            input->noext = clo_getOptionBool(option);

        } else if (strcmp(keyword, "merged") == 0) {
            if (clo_isOptionSet(option)) {
                parse_file_name(clo_getOptionString(option), tmp_file);
                strcpy(input->merged, tmp_file);
            }
        } else if (strcmp(keyword, "loneast") == 0) {
            input->loneast = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "lonwest") == 0) {
            input->lonwest = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "latnorth") == 0) {
            input->latnorth = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "latsouth") == 0) {
            input->latsouth = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "verbose") == 0) {
            input->verbose = clo_getOptionBool(option);

        } else if (strcmp(keyword, "unit_wgt") == 0) {
            input->unit_wgt = clo_getOptionInt(option);

        } else if (strcmp(keyword, "union_bins") == 0) {
            input->union_bins = clo_getOptionBool(option);

        } else if (strcmp(keyword, "median") == 0) {
            input->median = clo_getOptionInt(option);

        } else if (strcmp(keyword, "deflate") == 0) {
            input->deflate = clo_getOptionInt(option);

        } else if (strcmp(keyword, "oformat") == 0) {
            const char* tmpStr = getFileFormatName(clo_getOptionString(option));
            strcpy(input->oformat, tmpStr);

        } else if (strcmp(keyword, "composite_prod") == 0) {
            if (clo_isOptionSet(option)) {
                parse_file_name(clo_getOptionString(option), tmp_file);
                strcpy(input->composite_prod, tmp_file);
            }
        } else if (strcmp(keyword, "composite_scheme") == 0) {
            if (clo_isOptionSet(option)) {
                parse_file_name(clo_getOptionString(option), tmp_file);
                strcpy(input->composite_scheme, tmp_file);
            }
        } else if (strcmp(keyword, "doi") == 0) {
            if (clo_isOptionSet(option)) {
                strcpy(input->doi, clo_getOptionString(option));
            }
        } else {
            goto Invalid_return;

        }

    }
    return 0;


Invalid_return:
    printf("Invalid argument \"%s\"\n", keyword);
    return -1;
}

/*-----------------------------------------------------------------------------
    Function:  l3bin_input

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

int l3bin_input(int argc, char **argv, instr *input, const char* prog, const char* version) {
    char str_buf[4096];

    /*                                                                  */
    /* Set input values to defaults                                     */
    /*                                                                  */
    if (input_init(input) != 0) {
        printf("-E- %s: Error initializing input structure.\n", __FILE__);
        return (-1);
    }

    /* hold all of the command line options */
    clo_optionList_t* list;

    list = clo_createList();

    /* initialize the option list with descriptions and default values */
    l3bin_init_options(list, prog, version);

    if (argc == 1) {
        clo_printUsage(list);
        exit(1);
    }

    clo_setEnableDumpOptions(1);

    // read command line args to get the ifile parameter
    clo_readArgs(list, argc, argv);

    if (l3bin_load_input(list, input) != 0) {
        printf("-E- %s: Error loading options into input structure.\n", __FILE__);
        clo_deleteList(list);
        return (-1);
    }

    clo_deleteList(list);

    readFileList(input->infile);
    
    /*                                                                  */
    /* Build string of parameters for metadata                          */
    /*                                                                  */
    sprintf(str_buf, "infile = %s\n", input->infile);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "ofile = %s\n", input->ofile);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "pfile = %s\n", input->ofile);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "oformat = %s\n", input->oformat);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "syear = %d\n", input->syear);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "eyear = %d\n", input->eyear);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "sday = %d\n", input->sday);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "eday = %d\n", input->eday);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "sorbit = %d\n", input->sorbit);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "eorbit = %d\n", input->eorbit);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "out_parm = %s\n", input->out_parm);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "processing_version = %s\n", input->pversion);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "reduce_fac = %d\n", input->reduce_fac);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "resolve = %s\n", input->resolve);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "merged = %s\n", input->merged);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "loneast = %f\n", input->loneast);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "lonwest = %f\n", input->lonwest);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "latnorth = %f\n", input->latnorth);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "latsouth = %f\n", input->latsouth);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "verbose = %d\n", input->verbose);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "unit_wgt = %d\n", input->unit_wgt);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "median = %d\n", input->median);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "deflate = %d\n", input->deflate);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "composite_prod = %s\n", input->composite_prod);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "composite_scheme = %s\n", input->composite_scheme);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "doi = %s\n", input->doi);
    strcat(input->parms, str_buf);

    return 0;
}
