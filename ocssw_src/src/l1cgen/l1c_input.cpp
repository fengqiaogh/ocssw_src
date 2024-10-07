//******************************************
// l1c_input.cpp
//  Created by Martin Montes on 8/15/2022
//  lasr version on 12/14/2022
//****************************************

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <netcdf>
#include <unistd.h>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <clo.h>
#include <filetype.h>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

#include "l1c_input.h"
#include "genutils.h"
#include "passthebuck.h"
#include "sensorInfo.h"

#include "l1c_filehandle.h"
#include <iostream>

// store the name of program we are running.
static char mainProgramName[50];

namespace l1c {

static char *l1cgen_optionKeys[] = {
    "-help", "-version", "-dump_options", "-dump_options_paramfile", "-dump_options_xmlfile", "par",
    "pversion","doi","verbose","ifile","ofile","outlist","l1c_grid", "l1c_anc", "mode", "south", "north", "west",
    "east", "history","l2prod","ix_l1cprod",
    "selgran",          // selected granules up to 10 files, they are ids not indexes!!
    "selyear",          // selected year
    "selmon",           // selected month
    "selday",           // selected day
    "grid_resolution",  // grid resolution in km
    "sensor",           // SPEX 1, OCI 2 and HARP 3
    "gransize", "grantype", "bintype", "start_timeflag",
    "start_time",  // initial time as ISO for selecting granules
    "end_time",    // final time as ISO for selecting granules
    "projection",       // projection type,"swath_grid":0 (default-Fred) or "socea=1"
    "sort_method",      // binning sorting type: 0 orbital fred search, 1: SBS
    "cloud_height",     // cloud top height in km
    "demfile",
    "terrain_correct",  // terrain distortion correction , 1: yes
    "cloud_correct",    // 0: no parallax, L1C at cth=0, 1: L1C at cth = k, 2: L1B at cht=variable
    "cloud_type",
    "demcloud_flag",
    // multi  attributes (view, pol, bands)
    "overlap_vflag",  // tells if we want merged views
    "overlap_pflag",  // tells if we want merged polarizations
    "overlap_bflag",  // tells if we want merged spectral bands
    // uncertainty params l1c merged products
    "unc_meth",     // uncertainity calculation method
    "unc_thres_v",  // uncertainity threshold of angular merged products as %
    "unc_thres_p",  // same but for polarization
    "unc_thres_b",  // sam
    NULL};

L1C_input::L1C_input(){};  // constructor
L1C_input::~L1C_input(){};

int32_t L1C_input::l1c_usage(const char *prog, const char *ver) {
    clo_optionList_t *list;

    list = clo_createList();
    l1c_init_options(list, prog, ver);
    clo_printUsage(list);

    return 0;
}

int32_t L1C_input::l1c_init_options(clo_optionList_t *list, const char *prog, const char *version) {
    char tmpStr[2048];
    clo_option_t *option;

    if (strcmp(prog, "l1cgen"))
        clo_setSelectOptionKeys(l1cgen_optionKeys);

    // set the min program name
    strcpy(mainProgramName, prog);

    sprintf(tmpStr, "Usage: %s argument-list\n\n", prog);
    strcat(tmpStr, "  The argument-list is a set of keyword=value pairs. The arguments can\n");
    strcat(tmpStr, "  be specified on the commandline, or put into a parameter file, or the\n");
    strcat(tmpStr, "  two methods can be used together, with commandline over-riding.\n\n");
    strcat(tmpStr, "  return value: 0=ALL GOOD, 1=ERROR, 110=NO PIXELS BINNED\n");
    strcat(tmpStr, "  file data.\n\n");
    strcat(tmpStr, "The list of valid keywords follows:\n");
    clo_setHelpStr(tmpStr);

    // add the parfile alias for backward compatibility
    clo_addAlias(list, "par", "parfile");
    strcpy(tmpStr, "input L1b or L2 file name");
    option = clo_addOption(list, "ifile", CLO_TYPE_IFILE, NULL, tmpStr);
    clo_addOptionAlias(option, "infile");
    strcpy(tmpStr, "input L1C file name");
    option = clo_addOption(list, "l1c_grid", CLO_TYPE_IFILE, NULL, tmpStr);
    option = clo_addOption(list, "l1c_anc", CLO_TYPE_IFILE, NULL, tmpStr);
    clo_addOption(list, "ofile", CLO_TYPE_OFILE, "output", "output for 1 file name");
    clo_addOption(list, "outlist", CLO_TYPE_OFILE, "output", "output list file name");
    clo_addOption(list, "verbose", CLO_TYPE_BOOL, "off", "Allow more verbose screen messages");
    clo_addOption(list, "pversion", CLO_TYPE_STRING, NULL, "processing version string");
    clo_addOption(list, "doi", CLO_TYPE_STRING, NULL, "Digital Object Identifier (DOI) string");

    strcpy(tmpStr, "flags masked [see /SENSOR/l2bin_defaults.par]");
    strcpy(tmpStr, "l2prod = bin products [default=all products]\n");
    strcat(tmpStr, "    Set to \"ALL\" or \"all\" for all L2 products in 1st input file.\n");
    strcat(tmpStr, "    Use ',' as delimiters.\n");
    clo_addOption(list, "l2prod", CLO_TYPE_STRING, "ALL", tmpStr);
    // l1c options--------------------------
    //**************L1C options *****************************************
    strcpy(tmpStr, "L1C processing flag\n");
    strcat(tmpStr, "        5: L1C grid creation from HKT telemetry\n");
    strcat(tmpStr,
           "        7: CTH-corrected L1B (cloud height parallax) and L1C grid at cloud height from L1C "
           "granules with CTH=0--\n");
    strcat(tmpStr,
           "        8: L1C FULL file creation from L1B granules-and input L1C grid SOCEA-L1 readers\n");

    clo_addOption(list, "mode", CLO_TYPE_INT, "0", tmpStr);

    // it has a different in l2gen
    strcpy(tmpStr, "L1C grid min binning latitude");
    clo_addOption(list, "south", CLO_TYPE_FLOAT, "-90", tmpStr);
    strcpy(tmpStr, "L1C grid max binning latitude");
    clo_addOption(list, "north", CLO_TYPE_FLOAT, "+90", tmpStr);
    strcpy(tmpStr, "L1C grid min binning longitude");
    clo_addOption(list, "west", CLO_TYPE_FLOAT, "-180", tmpStr);
    strcpy(tmpStr, "L1C grid max binning longitude");
    clo_addOption(list, "east", CLO_TYPE_FLOAT, "+180", tmpStr);

    strcpy(tmpStr, "L1C processing of granules");
    strcpy(tmpStr, "granule id (1 to 10) note: not indexes!");
    clo_addOption(list, "selgran", CLO_TYPE_INT, "[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]", tmpStr);

    // not implemented---
    strcpy(tmpStr, "Index L1C product [0,0,0]");
    strcat(tmpStr, "        0: pc");
    strcat(tmpStr, "        1: vsf");
    strcat(tmpStr, "       2: dpr");
    clo_addOption(list, "ix_l1cprod", CLO_TYPE_INT, "[0,0,0]", tmpStr);

    // fixed bearing projection---
    strcpy(tmpStr, "Day of the year for processing L1C swath");
    strcat(tmpStr, "           units in day number (1-365/366)");
    clo_addOption(list, "selday", CLO_TYPE_INT, "-1", tmpStr);

    strcpy(tmpStr, "Month of the year for processing L1C swath");
    strcat(tmpStr, "           units in month number (1-12)");
    clo_addOption(list, "selmon", CLO_TYPE_INT, "-1", tmpStr);

    strcpy(tmpStr, "Year for processing L1C swath");
    strcat(tmpStr, "           units in year");
    clo_addOption(list, "selyear", CLO_TYPE_INT, "-1", tmpStr);

    //******************************************************************
    strcpy(tmpStr, "Common grid resolution");
    strcat(tmpStr, "           units in km");
    clo_addOption(list, "grid_resolution", CLO_TYPE_FLOAT, "5.2", tmpStr);

    strcpy(tmpStr, "PACE sensor to be gridded");
    strcat(tmpStr, "     e.g. SPEXONE, OCI, HARP2, MISR");
    clo_addOption(list, "sensor", CLO_TYPE_STRING, "OCI", tmpStr);  // default
    strcpy(tmpStr, "granule size for telemetry-derived L1C files");
    strcat(tmpStr, "    in minutes--5' by default");
    clo_addOption(list, "gransize", CLO_TYPE_INT, "5", tmpStr);

    strcpy(tmpStr, "granule type processing for telemetry-derived L1C files");
    strcat(tmpStr, "  0: granules, 1: swath, 0  by default");
    clo_addOption(list, "grantype", CLO_TYPE_INT, "0", tmpStr);

    strcpy(tmpStr, "binning type for binned L1C products");
    strcat(tmpStr, "  0: discrete, 1: area-weighting, 0  by default");
    clo_addOption(list, "bintype", CLO_TYPE_INT, "0", tmpStr);

    strcpy(tmpStr, "Initial time as ISO for processing HKT granules ->L1C files");
    strcat(tmpStr, "  units as ISO");
    clo_addOption(list, "start_time", CLO_TYPE_STRING, "2022-03-21T00:00:00", tmpStr);

    strcpy(tmpStr, "End time as ISO for processing HKT granules ->L1C files");
    strcat(tmpStr, "  units as ISO ");
    clo_addOption(list, "end_time", CLO_TYPE_STRING, "2022-03-21T00:00:00", tmpStr);

    strcpy(tmpStr, "log command history");
    clo_addOption(list, "history", CLO_TYPE_STRING, "l1cgen ifile=xx mode=5 ofile=xx", tmpStr);

    strcpy(tmpStr, "initial time flag for L1C granule ");
    strcat(tmpStr,
           "    0: time zero seconds of the day 00:00:00, 1: starting coverage time for the HKT file by "
           "default, 0 by default");
    clo_addOption(list, "start_timeflag", CLO_TYPE_INT, "0", tmpStr);

    strcpy(tmpStr, "swath # for 1 full orbit processing for telemetry-derived L1C files");
    strcat(tmpStr, "  1 or 2, ascending or descending");
    clo_addOption(list, "swath_num", CLO_TYPE_INT, "1", tmpStr);

    strcpy(tmpStr, "Projection type");
    strcat(tmpStr, "      0: SOCEA(default)");
    strcat(tmpStr, "      1: SOCEA-2");
    clo_addOption(list, "projection", CLO_TYPE_INT, "0", tmpStr);

    strcpy(tmpStr, "Binning sorting type");
    strcat(tmpStr, "      0: Orbital-vectorsdefault)");
    strcat(tmpStr, "      1: SADDLEBACK SEARCH");
    clo_addOption(list, "sort_method", CLO_TYPE_INT, "0", tmpStr);

    strcpy(tmpStr, "height type flag");
    strcat(tmpStr, "        0: geoid or L1C height");
    strcat(tmpStr, "        1: orthometric or dem height");
    clo_addOption(list, "demcloud_flag", CLO_TYPE_INT, "0", tmpStr);

    strcpy(tmpStr, "Cloud top height for L1C corrections");
    strcat(tmpStr, "      (km)");
    clo_addOption(list, "cloud_height", CLO_TYPE_INT, "0", tmpStr);

    strcpy(tmpStr, "Terrain correct flag");
    strcat(tmpStr, "        DEM correction");
    strcat(tmpStr, "        0: off");
    strcat(tmpStr, "        1: on");
    clo_addOption(list, "terrain_correct", CLO_TYPE_BOOL, "0", tmpStr);

    strcpy(tmpStr, "Cloud correct flag for L1C");
    strcat(tmpStr, "        0: CTH=0 no parallax correction (default)");
    strcat(tmpStr, "        1: CTH=k, constant CTH from ANC");
    strcat(tmpStr, "        2: CTH=k, constant CTH from L1C grid height");
    clo_addOption(list, "cloud_correct", CLO_TYPE_INT, "0", tmpStr);

    strcpy(tmpStr, "Cloud type for parallax correction");
    strcat(tmpStr, "        0: water (default)");
    strcat(tmpStr, "        1: ice");
    clo_addOption(list, "cloud_type", CLO_TYPE_INT, "0", tmpStr);

    strcpy(tmpStr, "overlap views flag");
    strcat(tmpStr, "        0: off");
    strcat(tmpStr, "        1: on");
    clo_addOption(list, "overlap_vflag", CLO_TYPE_BOOL, "0", tmpStr);

    strcpy(tmpStr, "overlap polarizations flag");
    strcat(tmpStr, "        0: off");
    strcat(tmpStr, "        1: on");
    clo_addOption(list, "overlap_pflag", CLO_TYPE_BOOL, "0", tmpStr);

    strcpy(tmpStr, "overlap spectral bands flag");
    strcat(tmpStr, "        0: off");
    strcat(tmpStr, "        1: on");
    clo_addOption(list, "overlap_bflag", CLO_TYPE_BOOL, "0", tmpStr);

    strcpy(tmpStr, "Uncertainty calculation method");
    strcat(tmpStr, "        0: error propagation");
    strcat(tmpStr, "        1: Monte Carlo");
    clo_addOption(list, "unc_meth", CLO_TYPE_INT, "0", tmpStr);

    strcpy(tmpStr, "Uncertainty threshold for angular merged product");
    strcat(tmpStr, "        as percentage");
    clo_addOption(list, "unc_thres_v", CLO_TYPE_FLOAT, "10", tmpStr);

    strcpy(tmpStr, "Uncertainty threshold for polarization  merged product");
    strcat(tmpStr, "        as percentage");
    clo_addOption(list, "unc_thres_p", CLO_TYPE_FLOAT, "10", tmpStr);

    strcpy(tmpStr, "Uncertainty threshold for spectral bands  merged product");
    strcat(tmpStr, "        as percentage");
    clo_addOption(list, "unc_thres_b", CLO_TYPE_FLOAT, "10", tmpStr);

    strcpy(tmpStr, "Digital Elevation Model file");
    strcat(tmpStr, " *.nc file");
    clo_addOption(list, "demfile", CLO_TYPE_STRING, "$OCDATAROOT/common/gebco_ocssw_v2020.nc", tmpStr);

    //*************************************************************************-

    clo_setVersion2(prog, version);
    return 0;
}

// copy input info from list into instr structure
int32_t L1C_input::l1c_load_input(clo_optionList_t *list, L1C_input *l1ccli) {
    char *tmp_str;
    char keyword[50];
    char *parm_str;
    char tmp_file[FILENAME_MAX];
    int numOptions, optionId;
    clo_option_t *option;

    int *iArray;
    int count = -1;
    char **cArray;

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
            if (isupper(*tmp_str))
                *tmp_str = tolower(*tmp_str);
            tmp_str++;
        }

        if (strcmp(keyword, "help") == 0) {
        } else if (strcmp(keyword, "version") == 0) {
        } else if (strncmp(keyword, "dump_options", 12) == 0) {
        } else if (strncmp(keyword, "par", 3) == 0) {
        } else if (strcmp(keyword, "ifile") == 0) {
            if (clo_isOptionSet(option)) {
                parm_str = clo_getOptionString(option);
                parse_file_name(parm_str, tmp_file);
                strcpy(l1ccli->infile, tmp_file);
            }
        } else if (strcmp(keyword, "ofile") == 0) {
            if (clo_isOptionSet(option)) {
                parm_str = clo_getOptionString(option);
                parse_file_name(parm_str, tmp_file);
                strcpy(l1ccli->ofile, tmp_file);
            }
        } else if (strcmp(keyword, "outlist") == 0) {
            if (clo_isOptionSet(option)) {
                parm_str = clo_getOptionString(option);
                parse_file_name(parm_str, tmp_file);
                strcpy(l1ccli->outlist, tmp_file);
            }
        } else if (strcmp(keyword, "l1c_grid") == 0) {
            if (clo_isOptionSet(option)) {
                parm_str = clo_getOptionString(option);
                parse_file_name(parm_str, tmp_file);
                strcpy(l1ccli->l1c_grid, tmp_file);
            }

        }
         else if (strcmp(keyword, "l1c_anc") == 0) {
            if (clo_isOptionSet(option)) {
                parm_str = clo_getOptionString(option);
                parse_file_name(parm_str, tmp_file);
                strcpy(l1ccli->l1c_anc, tmp_file);
            }

        } else if (strcmp(keyword, "l2prod") == 0) {
            parm_str = clo_getOptionRawString(option);
            strcpy(l1ccli->l2prod, parm_str);

        } else if (strcmp(keyword, "verbose") == 0) {
            l1ccli->verbose = clo_getOptionBool(option);

        } else if (strcmp(keyword, "pversion") == 0) {
            if (clo_isOptionSet(option)) {
                strcpy(l1ccli->pversion, clo_getOptionString(option));
            }
        } else if (strcmp(keyword, "doi") == 0) {
            if (clo_isOptionSet(option)) {
                strcpy(l1ccli->doi, clo_getOptionString(option));
            }
        }

        // L1C new options-----------------------------------------
        else if (strcmp(keyword, "mode") == 0) {
            l1ccli->l1c_pflag = clo_getOptionInt(option);
        } else if (strcmp(keyword, "south") == 0) {
            l1ccli->south = clo_getOptionFloat(option);
        } else if (strcmp(keyword, "north") == 0) {
            l1ccli->north = clo_getOptionFloat(option);
        } else if (strcmp(keyword, "west") == 0) {
            l1ccli->west = clo_getOptionFloat(option);
        } else if (strcmp(keyword, "east") == 0) {
            l1ccli->east = clo_getOptionFloat(option);
        } else if (strcmp(keyword, "selgran") == 0) {
            iArray = clo_getOptionInts(option, &count);
            if (count <= 10) {
                for (int i = 0; i < count; i++) {
                    if (iArray[i] <= 10)
                        l1ccli->selgran[i] = iArray[i];
                    else {
                        printf("-E- %s: Granule Id cant be larger than 10.\n", __FILE__);
                        exit(1);
                    }
                }
            } else {
                printf("-E- %s: Max number of granules to be processed is 10.\n", __FILE__);
                exit(1);
            }
        } else if (strcmp(keyword, "ix_l1cprod") == 0) {
            iArray = clo_getOptionInts(option, &count);
            for (int i = 0; i < count; i++)
                l1ccli->ix_l1cprod[i] = iArray[i];
        }

        else if (strcmp(keyword, "grid_resolution") == 0) {
            l1ccli->grid_resolution = clo_getOptionFloat(option);
        } else if (strcmp(keyword, "sensor") == 0) {
            if (clo_isOptionSet(option)) {
                const char *sname = clo_getOptionString(option);
                l1ccli->sensor = sensorName2SensorId(sname);              
            }
        } else if (strcmp(keyword, "gransize") == 0) {
            l1ccli->gransize = clo_getOptionInt(option);
        } else if (strcmp(keyword, "grantype") == 0) {
            l1ccli->grantype = clo_getOptionInt(option);
        } else if (strcmp(keyword, "bintype") == 0) {
            l1ccli->bintype = clo_getOptionInt(option);
        } else if (strcmp(keyword, "history") == 0) {
            if (clo_isOptionSet(option)) {
                cArray = clo_getOptionStrings(option, &count);
                string history;
                for (int i = 0; i < count; i++) {
                    string s1(cArray[i]);
                    history += s1;
                }
                strcpy(l1ccli->history, history.c_str());
            }
        } else if (strcmp(keyword, "start_time") == 0) {
            if (clo_isOptionSet(option)) {
                cArray = clo_getOptionStrings(option, &count);
                string s1(cArray[0]), s2(cArray[1]), s3(cArray[2]);
                string start_time = s1 + ":" + s2 + ":" + s3.substr(0, 2) + "Z";
                strcpy(l1ccli->start_time, start_time.c_str());
            }
        } else if (strcmp(keyword, "end_time") == 0) {
            if (clo_isOptionSet(option)) {
                cArray = clo_getOptionStrings(option, &count);
                string s1(cArray[0]), s2(cArray[1]), s3(cArray[2]);
                string gran_end_time = s1 + ":" + s2 + ":" + s3.substr(0, 2) + "Z";
                strcpy(l1ccli->end_time, gran_end_time.c_str());
            }
        } 
          else if (strcmp(keyword, "demfile") == 0) {
            if (clo_isOptionSet(option)) {
                cArray = clo_getOptionStrings(option, &count);
                string s1(cArray[0]);
                string gran_demfile = s1;
                strcpy(l1ccli->demfile, gran_demfile.c_str());
        }
  /*          else if (strcmp(keyword, "demfile") == 0) {
            parm_str = clo_getOptionString(option);   
            strcpy(l1ccli->demfile,parm_str);
        } */
          }    
            else if (strcmp(keyword, "start_timeflag") == 0) {
            l1ccli->start_timeflag = clo_getOptionInt(option);         
        } else if (strcmp(keyword, "selday") == 0) {
            l1ccli->selday = clo_getOptionInt(option);
        } else if (strcmp(keyword, "selmon") == 0) {
            l1ccli->selmon = clo_getOptionInt(option);
        } else if (strcmp(keyword, "selyear") == 0) {
            l1ccli->selyear = clo_getOptionInt(option);
        } else if (strcmp(keyword, "swath_num") == 0) {
            l1ccli->swath_num = clo_getOptionInt(option);
        } else if (strcmp(keyword, "projection") == 0) {
            l1ccli->projection = clo_getOptionInt(option);
        } else if (strcmp(keyword, "sort_method") == 0) {
            l1ccli->sort_method = clo_getOptionInt(option);
        } else if (strcmp(keyword, "demcloud_flag") == 0) {
            l1ccli->demcloud_flag = clo_getOptionInt(option);
        } else if (strcmp(keyword, "cloud_height") == 0) {
            l1ccli->cloud_height = clo_getOptionInt(option);
        } else if (strcmp(keyword, "terrain_correct") == 0) {
            l1ccli->terrain_correct = clo_getOptionBool(option);
        } else if (strcmp(keyword, "cloud_correct") == 0) {
            l1ccli->cloud_correct = clo_getOptionInt(option);         
        } else if (strcmp(keyword, "cloud_type") == 0) {
            l1ccli->cloud_type = clo_getOptionInt(option);
        } else if (strcmp(keyword, "overlap_vflag") == 0) {
            l1ccli->overlap_vflag = clo_getOptionBool(option);
        } else if (strcmp(keyword, "overlap_pflag") == 0) {
            l1ccli->overlap_pflag = clo_getOptionBool(option);
        } else if (strcmp(keyword, "overlap_bflag") == 0) {
            l1ccli->overlap_bflag = clo_getOptionBool(option);
        } else if (strcmp(keyword, "unc_meth") == 0) {
            l1ccli->unc_meth = clo_getOptionInt(option);
        } else if (strcmp(keyword, "unc_thres_v") == 0) {
            l1ccli->unc_thres_v = clo_getOptionFloat(option);
        } else if (strcmp(keyword, "unc_thres_p") == 0) {
            l1ccli->unc_thres_p = clo_getOptionFloat(option);
        } else if (strcmp(keyword, "unc_thres_b") == 0) {
            l1ccli->unc_thres_b = clo_getOptionFloat(option);
        } else {
            printf("-E- Invalid argument \"%s\"\n", keyword);
            clo_dumpOption(option);
            exit(1);
        }

    }  // end for

    return 0;
}

// it is called by l1c_input method
// init defaults, makes sense?
int32_t L1C_input::l1c_input_init(L1C_input *l1ccli) {
    l1ccli->infile[0] = '\0';
    l1ccli->ofile[0] = '\0';
    l1ccli->outlist[0] = '\0';
    l1ccli->l1c_grid[0] = '\0';
    l1ccli->l1c_anc[0] = '\0';
    l1ccli->verbose = 0;
    l1ccli->pversion[0] = '\0';
    l1ccli->doi[0] = '\0';
    l1ccli->l2prod[0] = '\0';
    for (int i = 0; i < 3; i++) {
        l1ccli->ix_l1cprod[i] = i+1;
    }  // 3x1 array with selected l1c products, 1: selected
    l1ccli->l1c_pflag = 0;
    l1ccli->south = -90;  // latitude in degrees
    l1ccli->north = 90;
    l1ccli->west = -180;
    l1ccli->east = 180;

    for (int i = 0; i < 10; i++) {
        l1ccli->selgran[i] = -1;  // first file of the list
    }
    l1ccli->swath_num = 1;
    l1ccli->grid_resolution = 5.2;  // grid resolution in km
    l1ccli->sensor = 30;            // 30 is OCI, 31 is OCIS
    l1ccli->gransize = 5;           // in minutes
    l1ccli->grantype = 0;
    l1ccli->bintype = 0;
    l1ccli->start_time[0] = '\0';
    l1ccli->end_time[0] = '\0';
    l1ccli->demfile[0]='\0';  // as ISO
    l1ccli->history[0] = '\0';
    l1ccli->start_timeflag = 0;  // flag=0 or to=0 in seconds of the day
    l1ccli->swath_num = 1;       // in minutes
    l1ccli->selyear = -1;
    l1ccli->selmon = -1;
    l1ccli->selday = -1;
    l1ccli->projection = 0;  // projection type,"swath_grid":0 or "socea=1"
    l1ccli->sort_method = 0;
    l1ccli->demcloud_flag = 0;
    l1ccli->cloud_height = 0;
    l1ccli->terrain_correct = 0;  // terrain distortion correction , 1: yes
    l1ccli->cloud_correct = 0;    // cloud distortion correction , 1: yes
    l1ccli->cloud_type = 0;
    // multi  attributes (view, pol, bands)
    l1ccli->overlap_vflag = 0;  // tells if we want merged views
    l1ccli->overlap_pflag = 0;  // tells if we want merged polarizations
    l1ccli->overlap_bflag = 0;  // tells if we want merged spectral bands
    // uncertainty params l1c merged products
    l1ccli->unc_meth = 0;          // uncertainity calculation method
    l1ccli->unc_thres_v = -999.0;  // uncertainity threshold of angular merged products as %
    l1ccli->unc_thres_p = -999.0;  // same but for polarization
    l1ccli->unc_thres_b = -999.0;  // same but for multispectral products, same view and polarization

    return 0;
}

/*-----------------------------------------------------------------------------
    Function:  l1c_input

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

int32_t L1C_input::l1c_inputmain(int argc, char **argv, L1C_input *l1cinput, l1c_filehandle *l1cfile,
                                 const char *prog, const char *version) {
    char *dataRoot;

    char localIfile[FILENAME_MAX], localIfile_l1c[FILENAME_MAX];
    char *ifile;
    string ifile_str;

    /* hold all of the command line options */
    clo_optionList_t *list;

    list = clo_createList();

    /* initialize the option list with descriptions and default values */
    l1c_init_options(list, prog, version);

    if (argc == 1) {
        clo_printUsage(list);
        exit(1);
    }

    // disable the dump option until we have read all of the files
    clo_setEnableDumpOptions(0);
    clo_readArgs(list, argc, argv);

    if (l1c_input_init(l1cinput) != 0) {
        printf("-E- %s: Error initializing l1c input structure.\n", __FILE__);
        return (-1);
    }

    // get list of input files
    strcpy(localIfile, clo_getString(list, "ifile"));
    l1cinput->files = readFileList(localIfile);  // files stored in a vector container in input struc

    if (l1cinput->files.size() == 0) {
        printf("No NetCDF input files found in %s.\n", localIfile);
        exit(EXIT_FAILURE);
    }

    ifile_str = l1cinput->files[0];
    ifile = (char*)ifile_str.c_str();

    l1cfile->l1b_name = ifile_str;

    if(l1cinput->verbose) cout<<"ifile.."<<ifile<<endl;
    file_format format = getFormat(ifile);  // reading from a nc file ---
    l1cfile->format = format.type;

    if(l1cinput->verbose)    
    {
        cout<<"l1cinput->files[0].."<<l1cfile->l1b_name<<endl;
        printf("format.type....%d..",format.type);
        printf("sensor id.....%d..",format.sensor_id);
    }   

    l1cfile->format = format.type;

    if(l1cinput->verbose) cout<<"sensor id.."<<format.sensor_id<<"sensor name.."<<sensorId2SensorName(format.sensor_id)<<endl;

    if (format.type == FT_INVALID) {
        printf("-E- %s Line %d: Could not find type for file %s.\n", __FILE__, __LINE__, ifile);
        return (-1);
    }

    if (format.sensor_id == -1) {
        printf("-E- Can not look up sensor ID for PLEASE PROVIDE PLATFORM FOR OCIS--!! %s.\n ", localIfile);
        cout << "forcing to be HARP2 when HARP2 L1B is beta from Meng/Richard" << endl;
        return (1);
    }

    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        printf("OCDATAROOT environment variable is not defined.\n");
        return (1);
    }

    // re-load the command line and par file

    clo_setEnableDumpOptions(1);
    clo_readArgs(list, argc, argv);

    // re-load the command line and par file
    if (want_verbose)
        printf("Loading command line parameters for L1C processing\n\n");
    // load input struct with command line arguments
    if (l1c_load_input(list, l1cinput) != 0) {
        printf("-E- %s: Error loading options into input structure.\n", __FILE__);
        return (-1);
    }

    if (l1cinput->l1c_pflag == 8 || l1cinput->l1c_pflag == 7 || l1cinput->l1c_pflag == 3) {
        strcpy(localIfile_l1c, clo_getString(list, "l1c_grid"));
        l1cinput->files_l1c = readFileList(localIfile_l1c);  // files stored in a vector container in input struc
        if(l1cinput->verbose)
        {
            cout << localIfile_l1c << endl;
            cout << l1cinput->files_l1c.size() << endl;
        }

        if (l1cinput->files_l1c.size() == 0) {
            printf("No NetCDF L1C input files found in %s.\n", localIfile_l1c);
            exit(EXIT_FAILURE);
        }
    }

     if(l1cinput->verbose) cout << "processing mode...." << l1cinput->l1c_pflag << endl;

    clo_deleteList(list);

    return 0;
}

}  // namespace l1c
