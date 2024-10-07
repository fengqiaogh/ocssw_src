#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <netcdf>
#include <unistd.h>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

#include "l2bin_input.h"
#include "genutils.h"
#include "passthebuck.h"
#include "sensorInfo.h"

// store the name of program we are running.
static char mainProgramName[50];

int l2bin_init_options(clo_optionList_t* list, const char* prog, const char* version) {
    char tmpStr[2048];
    clo_option_t* option;

    // set the min program name
    strcpy(mainProgramName, prog);

    sprintf(tmpStr, "Usage: %s argument-list\n\n", prog);
    strcat(tmpStr, "  The argument-list is a set of keyword=value pairs. The arguments can\n");
    strcat(tmpStr, "  be specified on the commandline, or put into a parameter file, or the\n");
    strcat(tmpStr, "  two methods can be used together, with commandline over-riding.\n\n");
    strcat(tmpStr, "  return value: 0=OK, 1=error, 110=north,south,east,west does not intersect\n");
    strcat(tmpStr, "  file data.\n\n");
    strcat(tmpStr, "The list of valid keywords follows:\n");
    clo_setHelpStr(tmpStr);

    // add the parfile alias for backward compatibility
    clo_addAlias(list, "par", "parfile");

    strcpy(tmpStr, "input L2 file name");
    option = clo_addOption(list, "ifile", CLO_TYPE_IFILE, NULL, tmpStr);
    clo_addOptionAlias(option, "infile");

    clo_addOption(list, "ofile", CLO_TYPE_OFILE, "output", "output file name");
    clo_addOption(list, "fileuse", CLO_TYPE_OFILE, NULL, "write filenames of the input files used into this file");

    clo_addOption(list, "suite", CLO_TYPE_STRING, NULL, "suite for default parameters");
    clo_addOption(list, "qual_prod", CLO_TYPE_STRING, NULL, "quality product field name");

    clo_addOption(list, "deflate", CLO_TYPE_INT, "5", "deflation level.  0=off or 1=low through 9=high");

    clo_addOption(list, "verbose", CLO_TYPE_BOOL, "off", "Allow more verbose screen messages");
    clo_addOption(list, "night", CLO_TYPE_BOOL, "off", "set to 1 for SST night processing");
    clo_addOption(list, "qual_max", CLO_TYPE_INT, "2", "maximum acceptable quality");
    clo_addOption(list, "rowgroup", CLO_TYPE_INT, "-1", "# of bin rows to process at once.");
    clo_addOption(list, "sday", CLO_TYPE_INT, "1970001", "start datadate (YYYYDDD) [ignored for \"regional\" prodtype]");
    clo_addOption(list, "eday", CLO_TYPE_INT, "2038018", "end datadate (YYYYDDD) [ignored for \"regional\" prodtype]");
    clo_addOption(list, "latnorth", CLO_TYPE_FLOAT, "90", "northern most latitude");
    clo_addOption(list, "latsouth", CLO_TYPE_FLOAT, "-90", "southern most latitude");
    clo_addOption(list, "loneast", CLO_TYPE_FLOAT, "0", "eastern most longitude");
    clo_addOption(list, "lonwest", CLO_TYPE_FLOAT, "0", "western most longitude");
    clo_addOption(list, "minobs", CLO_TYPE_INT, "0", "required minimum number of observations");

    strcpy(tmpStr, "equator crossing time delta in\n         minutes\n");
    strcat(tmpStr, "         Caveat...if zero, the sensor default equator crossing time will be used\n");
    strcat(tmpStr, "         This is not necessarily noon");
    clo_addOption(list, "delta_crossing_time", CLO_TYPE_FLOAT, "0.0", tmpStr);

    strcpy(tmpStr, "bin resolution\n");
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
    option = clo_addOption(list, "resolution", CLO_TYPE_STRING, "H", tmpStr);
    clo_addOptionAlias(option, "resolve");
    clo_addOption(list, "prodtype", CLO_TYPE_STRING, "day", "product type (Set to \"regional\" to bin all scans.)");
    clo_addOption(list, "pversion", CLO_TYPE_STRING, "unspecified", "production version");

    clo_addOption(list, "composite_scheme", CLO_TYPE_STRING, NULL, "composite scheme (min/max)");
    clo_addOption(list, "composite_prod", CLO_TYPE_STRING, NULL, "composite product fieldname");
    strcpy(tmpStr, "flags masked [see /SENSOR/l2bin_defaults.par]");
    clo_addOption(list, "flaguse", CLO_TYPE_STRING, DEF_FLAG, tmpStr);

    strcpy(tmpStr, "l3bprod = bin products [default=all products]\n");
    strcat(tmpStr, "        Set to \"ALL\" or \"all\" for all L2 products in 1st input file.\n");
    strcat(tmpStr, "        Use ':' or ',' or ' ' as delimiters.\n");
    strcat(tmpStr, "        Use ';' or '=' to delineate minimum values.");
    clo_addOption(list, "l3bprod", CLO_TYPE_STRING, "ALL", tmpStr);

    clo_addOption(list, "area_weighting", CLO_TYPE_INT, "0", "Enable area weighting\n        0: off\n        1: pixel box\n        2: pixel bounding box\n        3: pixel polygon");

    clo_addOption(list, "output_wavelengths", CLO_TYPE_STRING, "ALL", "comma separated list of\n        wavelengths for multi-wavelength products");
    clo_addOption(list, "doi", CLO_TYPE_STRING, NULL, "Digital Object Identifier (DOI) string");
    strcpy(tmpStr,"comma separated list of output L3 product names.\n");
    strcat(tmpStr, "        This option allows the user to specify the output product names which differ from the original l2 product names.\n");
    strcat(tmpStr, "        Usage: original_l2_name:output_l3_name, i.e. oprodname=cloud_flag:cloud_fraction");   
    clo_addOption(list, "oprodname", CLO_TYPE_STRING, NULL, tmpStr);
    clo_setVersion(version);
    return 0;
}

int l2bin_load_input(clo_optionList_t* list, instr *input) {
    char *tmp_str;
    char keyword[50];
    char *parm_str;
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

        if (strcmp(keyword, "help") == 0) {
        }
        else if (strcmp(keyword, "version") == 0) {
        }
        else if (strncmp(keyword, "dump_options", 12) == 0) {
        }
        else if (strncmp(keyword, "par", 3) == 0) {
        }
        else if (strcmp(keyword, "ifile") == 0) {
            if (clo_isOptionSet(option)) {
                parm_str = clo_getOptionString(option);
                parse_file_name(parm_str, tmp_file);
                strcpy(input->infile, tmp_file);
            }
        } else if (strcmp(keyword, "ofile") == 0) {
            if (clo_isOptionSet(option)) {
                parm_str = clo_getOptionString(option);
                parse_file_name(parm_str, tmp_file);
                strcpy(input->ofile, tmp_file);
            }
        } else if (strcmp(keyword, "fileuse") == 0) {
            if (clo_isOptionSet(option)) {
                parm_str = clo_getOptionString(option);
                parse_file_name(parm_str, tmp_file);
                strcpy(input->fileuse, tmp_file);
            }
        } else if (strcmp(keyword, "sday") == 0) {
            if (clo_isOptionSet(option))
                input->sday = clo_getOptionInt(option);

        } else if (strcmp(keyword, "eday") == 0) {
            if (clo_isOptionSet(option))
                input->eday = clo_getOptionInt(option);

        } else if (strcmp(keyword, "resolution") == 0) {
            parm_str = clo_getOptionString(option);
            parse_file_name(parm_str, tmp_file);
            strcpy(input->resolve, tmp_file);

        } else if (strcmp(keyword, "rowgroup") == 0) {
            input->rowgroup = clo_getOptionInt(option);

        } else if (strcmp(keyword, "flaguse") == 0) {
	    string flags = clo_getOptionRawString(option);
	    boost::replace_all(flags, "default", DEF_FLAG);
	    strcpy(input->flaguse, flags.c_str());

        } else if (strcmp(keyword, "l3bprod") == 0) {
            parm_str = clo_getOptionRawString(option);
            strcpy(input->l3bprod, parm_str);

        } else if (strcmp(keyword, "prodtype") == 0) {
            parm_str = clo_getOptionString(option);
            strcpy(input->prodtype, parm_str);

        } else if (strcmp(keyword, "output_wavelengths") == 0) {
            parm_str = clo_getOptionRawString(option);
            strcpy(input->output_wavelengths, parm_str);

        } else if (strcmp(keyword, "pversion") == 0) {
            parm_str = clo_getOptionString(option);
            strcpy(input->pversion, parm_str);

        } else if (strcmp(keyword, "suite") == 0) {
            if (clo_isOptionSet(option)) {
                parm_str = clo_getOptionString(option);
                strcpy(input->suite, parm_str);
            }
        } else if (strcmp(keyword, "latsouth") == 0) {
            input->latsouth = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "latnorth") == 0) {
            input->latnorth = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "lonwest") == 0) {
            input->lonwest = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "loneast") == 0) {
            input->loneast = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "meminfo") == 0) {
            input->meminfo = clo_getOptionInt(option);

        } else if (strcmp(keyword, "dcinfo") == 0) {
            input->dcinfo = clo_getOptionInt(option);

        } else if (strcmp(keyword, "night") == 0) {
            input->night = clo_getOptionBool(option);

        } else if (strcmp(keyword, "verbose") == 0) {
            input->verbose = clo_getOptionBool(option);

        } else if (strcmp(keyword, "minobs") == 0) {
            input->minobs = clo_getOptionInt(option);

        } else if (strcmp(keyword, "delta_crossing_time") == 0) {
            input->deltaeqcross = clo_getOptionFloat(option);

        } else if (strcmp(keyword, "deflate") == 0) {
            input->deflate = clo_getOptionInt(option);

        } else if (strcmp(keyword, "qual_max") == 0) {
            input->qual_max = (uint8_t) clo_getOptionInt(option);

        } else if (strcmp(keyword, "qual_prod") == 0) {
            if (clo_isOptionSet(option)) {
                parm_str = clo_getOptionString(option);
                parse_file_name(parm_str, tmp_file);
                strcpy(input->qual_prod, tmp_file);
            }

        } else if(strcmp(keyword, "oprodname") == 0){
            if (clo_isOptionSet(option)) {
                parm_str = clo_getOptionRawString(option);
                strcpy(input->output_product_names, parm_str);
                }
        } else if (strcmp(keyword, "composite_prod") == 0) {
            if (clo_isOptionSet(option)) {
                parm_str = clo_getOptionString(option);
                parse_file_name(parm_str, tmp_file);
                strcpy(input->composite_prod, tmp_file);
            }
        } else if (strcmp(keyword, "composite_scheme") == 0) {
            if (clo_isOptionSet(option)) {
                parm_str = clo_getOptionString(option);
                parse_file_name(parm_str, tmp_file);
                strcpy(input->composite_scheme, tmp_file);
            }
        } else if (strcmp(keyword, "area_weighting") == 0) {
            if (clo_isOptionSet(option)) {
                input->area_weighting = clo_getOptionInt(option);
            } else {
                input->area_weighting = 0;                
            }
        } else if (strcmp(keyword, "doi") == 0) {
            if (clo_isOptionSet(option)) {
                strcpy(input->doi, clo_getOptionString(option));
            }
        } else {
            printf("-E- Invalid argument \"%s\"\n", keyword);
            exit(EXIT_FAILURE);
        }

    }

    return 0;
}

int input_init(instr *input_str) {
    input_str->infile[0] = '\0';
    input_str->ofile[0] = '\0';
    input_str->pfile[0] = '\0';

    input_str->fileuse[0] = '\0';
    input_str->qual_prod[0] = '\0';
    input_str->composite_prod[0] = '\0';
    input_str->composite_scheme[0] = '\0';
    input_str->output_product_names[0] = '\0';
    strcpy(input_str->pversion, "Unspecified");
    strcpy(input_str->prodtype, "day");

    strcpy(input_str->l3bprod, "ALL");
    strcpy(input_str->output_wavelengths, "ALL");

    input_str->sday = 1970001;
    input_str->eday = 2038018;

    input_str->resolve[0] = '\0';

    input_str->rowgroup = -1;

    input_str->night = 0;
    input_str->verbose = 0;
    input_str->minobs = 0;
    input_str->deltaeqcross = 0.0;

    input_str->meminfo = 0;
    input_str->dcinfo = 0;

    input_str->latsouth = -90.0;
    input_str->latnorth = +90.0;
    input_str->lonwest = 0.0;
    input_str->loneast = 0.0;

    input_str->qual_max = 255;

    input_str->deflate = 0;

    strcpy(input_str->suite, "");
    
    input_str->area_weighting = 0;
    input_str->doi[0] = '\0';

    return 0;
}

/*-----------------------------------------------------------------------------
    Function:  l2bin_input

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

int l2bin_input(int argc, char **argv, instr *input, const char* prog, const char* version) {

    char str_buf[4096];

    char *dataRoot;
    int sensorId;
    int subsensorId = -1;
    char localSuite[FILENAME_MAX];
    char localIfile[FILENAME_MAX];

    /*                                                                  */
    /* Set input values to defaults                                     */
    /*                                                                  */
    if (input_init(input) != 0) {
        printf("-E- %s: Error initializing input structure.\n", __FILE__);
        exit(EXIT_FAILURE);
    }

    /* hold all of the command line options */
    clo_optionList_t* list;

    list = clo_createList();

    /* initialize the option list with descriptions and default values */
    l2bin_init_options(list, prog, version);

    if (argc == 1) {
        clo_printUsage(list);
        exit(EXIT_SUCCESS);
    }

    // disable the dump option until we have read all of the files
    clo_setEnableDumpOptions(0);
    clo_readArgs(list, argc, argv);

    // get list of input files
    strcpy(localIfile, clo_getString(list, "ifile"));
    input->files = readFileList(localIfile);
    if (input->files.size() == 0) {
        printf("-E- No NetCDF input files found in %s.\n", localIfile);
        exit(EXIT_FAILURE);
    }

    // see if suite param was set
    localSuite[0] = '\0';
    if (clo_isSet(list, "suite")) {
        strcpy(localSuite, clo_getString(list, "suite"));
    } // suite option was set

    // find the sensor and sub-sensor ID for first input file
    std::string instrument, platform;
    try {
        NcFile nc_input(input->files[0], NcFile::read);
        nc_input.getAtt("instrument").getValues(instrument);
        nc_input.getAtt("platform").getValues(platform);
        nc_input.close();
    } catch (NcException const & e) {
        printf("-E- L2 file %s does not have instrument/platform attributes.\n", input->files[0].c_str());
        exit(EXIT_FAILURE);
    }
    sensorId = instrumentPlatform2SensorId(instrument.c_str(), platform.c_str());
    subsensorId = sensorId2SubsensorId(sensorId);

    if (sensorId == -1) {
        printf("-E- Can not look up sensor ID for %s.\n", localIfile);
        exit(EXIT_FAILURE);
    }

    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        printf("OCDATAROOT environment variable is not defined.\n");
        exit(EXIT_FAILURE);
    }

    // load l2bin program defaults
    sprintf(str_buf, "%s/common/l2bin_defaults.par", dataRoot);
    if (access(str_buf, R_OK) != -1) {
        if (want_verbose)
            printf("Loading default parameters from %s\n", str_buf);
        clo_readFile(list, str_buf);
    }

    // sensor defaults
    sprintf(str_buf, "%s/%s/l2bin_defaults.par", dataRoot, sensorId2SensorDir(sensorId));
    if (access(str_buf, R_OK) != -1) {
        if (want_verbose)
            printf("Loading default parameters from %s\n", str_buf);
        clo_readFile(list, str_buf);
    }

    // subsensor defaults
    if (subsensorId != -1) {
        sprintf(str_buf, "%s/%s/%s/l2bin_defaults.par", dataRoot,
                sensorId2SensorDir(sensorId), subsensorId2SubsensorDir(subsensorId));
        if (access(str_buf, R_OK) != -1) {
            if (want_verbose)
                printf("Loading default parameters from %s\n", str_buf);
            clo_readFile(list, str_buf);
        }
    }

    // load suite default files
    if (localSuite[0] == 0) {
        if (clo_isSet(list, "suite"))
            strcpy(localSuite, clo_getString(list, "suite"));
    }

    // Check for suite entry
    if (localSuite[0] != 0) {
        int suiteLoaded = 0;

        // load common suite defaults
        sprintf(str_buf, "%s/common/l2bin_defaults_%s.par", dataRoot, localSuite);
        if (access(str_buf, R_OK) != -1) {
            suiteLoaded = 1;
            if (want_verbose)
                printf("Loading default parameters from %s\n", str_buf);
            clo_readFile(list, str_buf);
        }

        // sensor suite defaults
        sprintf(str_buf, "%s/%s/l2bin_defaults_%s.par", dataRoot,
                sensorId2SensorDir(sensorId), localSuite);
        if (access(str_buf, R_OK) != -1) {
            suiteLoaded = 1;
            if (want_verbose)
                printf("Loading default parameters from %s\n", str_buf);
            clo_readFile(list, str_buf);
        }

        // subsensor suite defaults
        if (subsensorId != -1) {
            sprintf(str_buf, "%s/%s/%s/l2bin_defaults_%s.par", dataRoot,
                    sensorId2SensorDir(sensorId), subsensorId2SubsensorDir(subsensorId),
                    localSuite);
            if (access(str_buf, R_OK) != -1) {
                suiteLoaded = 1;
                if (want_verbose)
                    printf("Loading default parameters from %s\n", str_buf);
                clo_readFile(list, str_buf);
            }
        }

        if (!suiteLoaded) {
            printf("-E- Failed to load parameters for suite %s for sensor %s\n", localSuite,
                   sensorId2SensorName(sensorId));
            exit(EXIT_FAILURE);
        }

    }

    // re-load the command line and par file
    if (want_verbose)
        printf("Loading command line parameters\n\n");
    clo_setEnableDumpOptions(1);
    clo_readArgs(list, argc, argv);

    // load input struct with command line arguments
    if (l2bin_load_input(list, input) != 0) {
        printf("-E- %s: Error loading options into input structure.\n", __FILE__);
        exit(EXIT_FAILURE);
    }

    /*                                                                  */
    /* Build string of parameters for metadata                          */
    /*                                                                  */
    sprintf(str_buf, "infile = %s\n", input->infile);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "ofile = %s\n", input->ofile);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "fileuse = %s\n", input->fileuse);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "sday = %d\n", input->sday);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "eday = %d\n", input->eday);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "latnorth = %f\n", input->latnorth);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "latsouth = %f\n", input->latsouth);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "loneast = %f\n", input->loneast);
    strcat(input->parms, str_buf);
    sprintf(str_buf, "lonwest = %f\n", input->lonwest);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "resolve = %s\n", input->resolve);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "rowgroup = %d\n", input->rowgroup);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "flaguse = %s\n", input->flaguse);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "l3bprod = %s\n", input->l3bprod);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "output_wavelengths = %s\n", input->output_wavelengths);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "prodtype = %s\n", input->prodtype);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "pversion = %s\n", input->pversion);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "suite = %s\n", input->suite);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "night = %d\n", input->night);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "verbose = %d\n", input->verbose);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "minobs = %d\n", input->minobs);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "delta_crossing_time = %f\n", input->deltaeqcross);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "deflate = %d\n", input->deflate);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "qual_prod = %s\n", input->qual_prod);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "composite_prod = %s\n", input->composite_prod);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "composite_scheme = %s\n", input->composite_scheme);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "qual_max = %d\n", input->qual_max);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "area_weighting = %d\n", input->area_weighting);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "doi = %s\n", input->doi);
    strcat(input->parms, str_buf);

    sprintf(str_buf, "oprodname = %s\n", input->output_product_names);
    strcat(input->parms, str_buf);
    
    clo_deleteList(list);

    return 0;
}
