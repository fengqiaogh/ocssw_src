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
#include "l2_utils.hpp"
// store the name of program we are running.
static char mainProgramName[50];

int l2bin_init_options(clo_optionList_t* list, const char* prog, const char* version) {
    char tmpStr[2048];
    clo_option_t* option;

    // set the min program name
    strncpy(mainProgramName, prog, sizeof(mainProgramName) - 1);

    snprintf(tmpStr,sizeof(tmpStr), "Usage: %s argument-list\n\n", prog);
    strncat(tmpStr, "  The argument-list is a set of keyword=value pairs. The arguments can\n", sizeof(tmpStr) - strlen(tmpStr) -1);
    strncat(tmpStr, "  be specified on the commandline, or put into a parameter file, or the\n", sizeof(tmpStr) - strlen(tmpStr) -1);
    strncat(tmpStr, "  two methods can be used together, with commandline over-riding.\n\n", sizeof(tmpStr) - strlen(tmpStr) -1);
    strncat(tmpStr, "  return value: 0=OK, 1=error, 110=north,south,east,west does not intersect\n", sizeof(tmpStr) - strlen(tmpStr) -1);
    strncat(tmpStr, "  file data.\n\n", sizeof(tmpStr) - strlen(tmpStr) -1);
    strncat(tmpStr, "The list of valid keywords follows:\n", sizeof(tmpStr) - strlen(tmpStr) -1);
    clo_setHelpStr(tmpStr);

    // add the parfile alias for backward compatibility
    clo_addAlias(list, "par", "parfile");

    strncpy(tmpStr, "input L2 file name", sizeof(tmpStr) - 1);
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

    strncpy(tmpStr, "equator crossing time delta in\n         minutes\n", sizeof(tmpStr) - 1);
    strncat(tmpStr, "         Caveat...if zero, the sensor default equator crossing time will be used\n", sizeof(tmpStr) - strlen(tmpStr) - 1);
    strncat(tmpStr, "         This is not necessarily noon", sizeof(tmpStr) - strlen(tmpStr) - 1);
    clo_addOption(list, "delta_crossing_time", CLO_TYPE_FLOAT, "0.0", tmpStr);

    strncpy(tmpStr, "bin resolution\n", sizeof(tmpStr) - strlen(tmpStr) - 1);
    strncat(tmpStr, "         H: 0.5km\n", sizeof(tmpStr) - strlen(tmpStr) - 1);
    strncat(tmpStr, "         Q: 250m\n", sizeof(tmpStr) - strlen(tmpStr) - 1);
    strncat(tmpStr, "        HQ: 100m\n", sizeof(tmpStr) - strlen(tmpStr) - 1);
    strncat(tmpStr, "        HH: 50m\n", sizeof(tmpStr) - strlen(tmpStr) - 1);
    strncat(tmpStr, "         1: 1.1km\n", sizeof(tmpStr) - strlen(tmpStr) - 1);
    strncat(tmpStr, "         2: 2.3km\n", sizeof(tmpStr) - strlen(tmpStr) - 1);
    strncat(tmpStr, "         4: 4.6km\n", sizeof(tmpStr) - strlen(tmpStr) - 1);
    strncat(tmpStr, "         9: 9.2km\n", sizeof(tmpStr) - strlen(tmpStr) - 1);
    strncat(tmpStr, "        18: 18.5km\n", sizeof(tmpStr) - strlen(tmpStr) - 1);
    strncat(tmpStr, "        36: 36km\n", sizeof(tmpStr) - strlen(tmpStr) - 1);
    strncat(tmpStr, "        1D: 1 degree\n", sizeof(tmpStr) - strlen(tmpStr) - 1);  
    strncat(tmpStr, "        HD: 0.5 degree\n", sizeof(tmpStr) - strlen(tmpStr) - 1);
    strncat(tmpStr, "        QD: 0.25 degree", sizeof(tmpStr) - strlen(tmpStr) - 1); 
    option = clo_addOption(list, "resolution", CLO_TYPE_STRING, "H", tmpStr);
    clo_addOptionAlias(option, "resolve");
    clo_addOption(list, "prodtype", CLO_TYPE_STRING, "day", "product type (Set to \"regional\" to bin all scans.)");
    clo_addOption(list, "pversion", CLO_TYPE_STRING, "unspecified", "production version");

    clo_addOption(list, "composite_scheme", CLO_TYPE_STRING, NULL, "composite scheme (min/max)");
    clo_addOption(list, "composite_prod", CLO_TYPE_STRING, NULL, "composite product fieldname");
    strncpy(tmpStr, "flags masked [see /SENSOR/l2bin_defaults.par]", sizeof(tmpStr) - 1);
    clo_addOption(list, "flaguse", CLO_TYPE_STRING, DEF_FLAG, tmpStr);

    strncpy(tmpStr, "l3bprod = bin products [default=all products]\n", sizeof(tmpStr) - 1);
    strncat(tmpStr, "        Set to \"ALL\" or \"all\" for all L2 products in 1st input file.\n", sizeof(tmpStr) - strlen(tmpStr) - 1);
    strncat(tmpStr, "        Use ',' or ' ' as delimiters between products.\n", sizeof(tmpStr) - strlen(tmpStr) - 1);
    strncat(tmpStr, "        Use '=' and ':' to delineate min/max values. i.e. prod=min:max", sizeof(tmpStr) - strlen(tmpStr) - 1);
    clo_addOption(list, "l3bprod", CLO_TYPE_STRING, "ALL", tmpStr);

    clo_addOption(list, "area_weighting", CLO_TYPE_INT, "0", "Enable area weighting\n        0: off\n        1: pixel box\n        2: pixel bounding box\n        3: pixel polygon");

    clo_addOption(list, "output_wavelengths", CLO_TYPE_STRING, "ALL", "comma separated list of\n        wavelengths for multi-wavelength products.\n        Usage - 'output_wavelengths=wave1=354:2200;wave2=870,1640'");
    clo_addOption(list, "doi", CLO_TYPE_STRING, NULL, "Digital Object Identifier (DOI) string");
    strncpy(tmpStr,"comma separated list of output L3 product names.\n", sizeof(tmpStr) - 1);   
    strncat(tmpStr, "        This option allows the user to specify the output product names which differ from the original l2 product names.\n", sizeof(tmpStr) - strlen(tmpStr) - 1);
    strncat(tmpStr, "        Usage - 'original_l2_name:output_l3_name', i.e. 'oprodname=cloud_flag:cloud_fraction'", sizeof(tmpStr) - strlen(tmpStr) - 1);   
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

        strncpy(keyword, option->key, sizeof(keyword) - 1);

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
                strncpy(input->ofile, tmp_file, sizeof(input->ofile) - 1);
            }
        } else if (strcmp(keyword, "fileuse") == 0) {
            if (clo_isOptionSet(option)) {
                parm_str = clo_getOptionString(option);
                parse_file_name(parm_str, tmp_file);
                strncpy(input->fileuse, tmp_file, sizeof(input->fileuse) - 1);
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
            strncpy(input->resolve, tmp_file, sizeof(input->resolve) - 1);

        } else if (strcmp(keyword, "rowgroup") == 0) {
            input->rowgroup = clo_getOptionInt(option);

        } else if (strcmp(keyword, "flaguse") == 0) {
	    string flags = clo_getOptionRawString(option);
	    boost::replace_all(flags, "default", DEF_FLAG);
	    strncpy(input->flaguse, flags.c_str(), sizeof(input->flaguse) - 1);

        } else if (strcmp(keyword, "l3bprod") == 0) {
            parm_str = clo_getOptionRawString(option);
            strncpy(input->l3bprod, parm_str, sizeof(input->l3bprod) - 1);

        } else if (strcmp(keyword, "prodtype") == 0) {
            parm_str = clo_getOptionString(option);
            strncpy(input->prodtype, parm_str, sizeof(input->prodtype) - 1);

        } else if (strcmp(keyword, "output_wavelengths") == 0) {
            parm_str = clo_getOptionRawString(option);
            strncpy(input->output_wavelengths, parm_str, sizeof(input->output_wavelengths) - 1);

        } else if (strcmp(keyword, "pversion") == 0) {
            parm_str = clo_getOptionString(option);
            strncpy(input->pversion, parm_str, sizeof(input->pversion) - 1);

        } else if (strcmp(keyword, "suite") == 0) {
            if (clo_isOptionSet(option)) {
                parm_str = clo_getOptionString(option);
                strncpy(input->suite, parm_str, sizeof(input->suite) - 1);
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
                strncpy(input->qual_prod, tmp_file, sizeof(input->qual_prod) - 1);
            }

        } else if(strcmp(keyword, "oprodname") == 0){
            if (clo_isOptionSet(option)) {
                parm_str = clo_getOptionRawString(option);
                strncpy(input->output_product_names, parm_str, sizeof(input->output_product_names) - 1);
                }
        } else if (strcmp(keyword, "composite_prod") == 0) {
            if (clo_isOptionSet(option)) {
                parm_str = clo_getOptionString(option);
                parse_file_name(parm_str, tmp_file);
                strncpy(input->composite_prod, tmp_file, sizeof(input->composite_prod) - 1);
            }
        } else if (strcmp(keyword, "composite_scheme") == 0) {
            if (clo_isOptionSet(option)) {
                parm_str = clo_getOptionString(option);
                parse_file_name(parm_str, tmp_file);
                strncpy(input->composite_scheme, tmp_file, sizeof(input->composite_scheme) - 1);
            }
        } else if (strcmp(keyword, "area_weighting") == 0) {
            if (clo_isOptionSet(option)) {
                input->area_weighting = clo_getOptionInt(option);
            } else {
                input->area_weighting = 0;                
            }
        } else if (strcmp(keyword, "doi") == 0) {
            if (clo_isOptionSet(option)) {
                strncpy(input->doi, clo_getOptionString(option), sizeof(input->doi) - 1);
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
    strncpy(input_str->pversion, "Unspecified", sizeof(input_str->pversion) - 1); 
    strncpy(input_str->prodtype, "day", sizeof(input_str->prodtype) - 1);

    strncpy(input_str->l3bprod, "ALL", sizeof(input_str->l3bprod) - 1);
    strncpy(input_str->output_wavelengths, "ALL", sizeof(input_str->output_wavelengths) - 1);

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

    strncpy(input_str->suite, "", sizeof(input_str->suite) - 1);
    
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

    char str_buf[4096 + 32]; // padding due to ofile, ifile sizes (4096)

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
    strncpy(localIfile, clo_getString(list, "ifile"), sizeof(localIfile) - 1);
    input->files = readFileList(localIfile);
    if (input->files.size() == 0) {
        printf("-E- No NetCDF input files found in %s.\n", localIfile);
        exit(EXIT_FAILURE);
    }

    // see if suite param was set
    localSuite[0] = '\0';
    if (clo_isSet(list, "suite")) {
        strncpy(localSuite, clo_getString(list, "suite"), sizeof(localSuite) - 1);
    } // suite option was set

    // find the sensor and sub-sensor ID for first input file
    std::string instrument, platform;
    try {
        NcFile nc_input(input->files[0], NcFile::read);
        nc_input.getAtt("instrument").getValues(instrument);
        nc_input.getAtt("platform").getValues(platform);
        nc_input.close();
    } catch (NcException const & e) {
        printf("-Warning-: L2 file %s does not have instrument/platform attributes.\n", input->files[0].c_str());
    }
    sensorId = instrumentPlatform2SensorId(instrument.c_str(), platform.c_str());
    subsensorId = sensorId2SubsensorId(sensorId);

    if (sensorId == -1) {
        printf("-Warning-: Can not look up sensor ID for %s.\n", localIfile);
    }

    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        printf("OCDATAROOT environment variable is not defined.\n");
        exit(EXIT_FAILURE);
    }

    // load l2bin program defaults
    snprintf(str_buf, sizeof(str_buf), "%s/common/l2bin_defaults.par", dataRoot);
    if (access(str_buf, R_OK) != -1) {
        if (want_verbose)
            printf("Loading default parameters from %s\n", str_buf);
        clo_readFile(list, str_buf);
    }

    // sensor defaults
    snprintf(str_buf, sizeof(str_buf), "%s/%s/l2bin_defaults.par", dataRoot, sensorId2SensorDir(sensorId));
    if (access(str_buf, R_OK) != -1) {
        if (want_verbose)
            printf("Loading default parameters from %s\n", str_buf);
        clo_readFile(list, str_buf);
    }

    // subsensor defaults
    if (subsensorId != -1) {
        snprintf(str_buf, sizeof(str_buf), "%s/%s/%s/l2bin_defaults.par", dataRoot,
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
        snprintf(str_buf, sizeof(str_buf), "%s/common/l2bin_defaults_%s.par", dataRoot, localSuite);
        if (access(str_buf, R_OK) != -1) {
            suiteLoaded = 1;
            if (want_verbose)
                printf("Loading default parameters from %s\n", str_buf);
            clo_readFile(list, str_buf);
        }

        // sensor suite defaults
        snprintf(str_buf, sizeof(str_buf), "%s/%s/l2bin_defaults_%s.par", dataRoot,
                sensorId2SensorDir(sensorId), localSuite);
        if (access(str_buf, R_OK) != -1) {
            suiteLoaded = 1;
            if (want_verbose)
                printf("Loading default parameters from %s\n", str_buf);
            clo_readFile(list, str_buf);
        }

        // subsensor suite defaults
        if (subsensorId != -1) {
            snprintf(str_buf, sizeof(str_buf), "%s/%s/%s/l2bin_defaults_%s.par", dataRoot,
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
    snprintf(str_buf, sizeof(str_buf), "infile = %s\n", input->infile);
    strncpy(input->parms, str_buf, sizeof(input->parms) - 1);
    snprintf(str_buf, sizeof(str_buf), "ofile = %s\n", input->ofile);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);
    snprintf(str_buf, sizeof(str_buf), "fileuse = %s\n", input->fileuse);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "sday = %d\n", input->sday);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);
    snprintf(str_buf, sizeof(str_buf), "eday = %d\n", input->eday);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "latnorth = %f\n", input->latnorth);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);
    snprintf(str_buf, sizeof(str_buf), "latsouth = %f\n", input->latsouth);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);
    snprintf(str_buf, sizeof(str_buf), "loneast = %f\n", input->loneast);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);
    snprintf(str_buf, sizeof(str_buf), "lonwest = %f\n", input->lonwest);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "resolve = %s\n", input->resolve);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "rowgroup = %d\n", input->rowgroup);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "flaguse = %s\n", input->flaguse);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "l3bprod = %s\n", input->l3bprod);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "output_wavelengths = %s\n", input->output_wavelengths);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "prodtype = %s\n", input->prodtype);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "pversion = %s\n", input->pversion);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "suite = %s\n", input->suite);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "night = %d\n", input->night);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "verbose = %d\n", input->verbose);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "minobs = %d\n", input->minobs);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "delta_crossing_time = %f\n", input->deltaeqcross);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "deflate = %d\n", input->deflate);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "qual_prod = %s\n", input->qual_prod);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "composite_prod = %s\n", input->composite_prod);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "composite_scheme = %s\n", input->composite_scheme);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "qual_max = %d\n", input->qual_max);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "area_weighting = %d\n", input->area_weighting);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "doi = %s\n", input->doi);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);

    snprintf(str_buf, sizeof(str_buf), "oprodname = %s\n", input->output_product_names);
    strncat(input->parms, str_buf, sizeof(input->parms) - strlen(input->parms) - 1);
    
    clo_deleteList(list);

    return 0;
}
