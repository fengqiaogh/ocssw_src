#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <genutils.h>
#include <clo.h>
#include <sensorInfo.h>
#include "netcdf.h"

#include "l2qc_viirs.h"

/** add all of the accepted command line options to list */
int l2qcviirs_init_options(clo_optionList_t* list, const char* softwareVersion) {
    char tmpStr[2048];
    char *dataRoot;

    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        fprintf(stderr, "-E- OCDATAROOT environment variable is not defined.\n");
        return (-1);
    }

    clo_setVersion2("l2qc_viirs", softwareVersion);

    sprintf(tmpStr, "Usage: l2qc_viirs argument-list\n\n");

    strcat(tmpStr, "  This program checks metadata of a VIIRS L2 data.\n");

    clo_setHelpStr(tmpStr);

    clo_addOption(list, "ifile", CLO_TYPE_IFILE, NULL, "input L2 file name");    
    
    clo_addOption(list, "cfile", CLO_TYPE_IFILE, "default", "input L2qc configuration file name");

    clo_addOption(list, "ofile", CLO_TYPE_OFILE, NULL, "output filename");

    return 0;
}

/* 
   Read the command line option and all of the default parameter files.

   This is the order for loading the options:
    - load the main program defaults file
    - load the command line (including specified par files)
    - re-load the command line disabling file descending so command
       line arguments will over ride

 */
int l2qcviirs_read_options(clo_optionList_t* list, int argc, char* argv[]) {
    char *dataRoot;
    char tmpStr[FILENAME_MAX];

    assert(list);

    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        fprintf(stderr, "-E- OCDATAROOT environment variable is not defined.\n");
        return (-1);
    }

    clo_readArgs(list, argc, argv);
    
    
    // find default config file (cfile) if not set in options
    if (strcmp(clo_getString(list, "cfile"),"default")==0) {
    	
    	// read platform and instrument directories from L2 file metadata
		int ncid; // netCDF ID 
		char *str_L2file = NULL; // VIIRS L2 file
		char platformStr[100];
		char instrumentStr[100];
		str_L2file = clo_getString(list, "ifile");
		nc_open(str_L2file, NC_NOWRITE, &ncid);
		nc_get_att_text(ncid, NC_GLOBAL, "platform", platformStr);
		nc_get_att_text(ncid, NC_GLOBAL, "instrument", instrumentStr);
		nc_close(ncid);
		int sensorId = instrumentPlatform2SensorId(instrumentStr, platformStr);
		int subsensorId = sensorId2SubsensorId(sensorId);
		
		// load the sensor specific default config file
		sprintf(tmpStr, "%s/%s/%s/l2qc_viirs.conf", dataRoot, sensorId2SensorDir(sensorId),subsensorId2SubsensorDir(subsensorId));				
		clo_setString(list, "cfile", tmpStr, "internal");   
	}

    return 0;
}
