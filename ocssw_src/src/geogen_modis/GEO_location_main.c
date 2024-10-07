#include "imsl_wrap.h"
#include "GEO_DEM.h"
#include "GEO_main.h"
#include "GEO_main_func.h"
#include "GEO_earth.h"
#include "version.h"

int main (void)
/*****************************************************************************
!C

!Description:
	Routine in Main group of the Level-1A geolocation software to call
	routines to read the parameter file, retrieve the file/path name for
	each input/output file pair, and sequence through each granule.

!Input Parameters:
	None

!Output Parameters:
	None

Return Values:
	EXIT_SUCCESS if all granules processed and output successfully,
	EXIT_FAIL otherwise

Externally Defined:
	__DATE__		defined by compiler
	__TIME__		defined by compiler
	DEM_resolutions		"GEO_DEM.h"
	EXIT_FAIL               "GEO_basic.h"
	FAIL                    "hdf.h"
	geoID                   "GEO_main.h"
	l1aID                   "GEO_main.h"
	MAKEFILE_REVISION	MOD_PR03.mk
	MODIS_E_GEO             "PGS_MODIS_35251.h"
	MODIS_U_GEO_BEGIN       "PGS_MODIS_35251.h"
	MODIS_U_GEO_END         "PGS_MODIS_35251.h"
	PGS_S_SUCCESS		"PGS_SMF.h"
	PGSd_DEM_WATER_LAND     "PGS_DEM.h"
	PGSd_DEM_30ARC          "PGS_DEM.h"
	PGSd_PC_LINE_LENGTH_MAX "PGS_PC.h"
	PGSd_PC_FILE_PATH_MAX   "PGS_PC.h"
	RESOLUTIONS		"GEO_DEM.h"
	TERRAIN_CORRECT_LUN     "GEO_geo.h"
	TERRAIN_CORRECT         "GEO_geo.h"

Called by:
	None

Routines called:
	imsl_error_options	"imsl_wrap.h"
	modsmf			"smfio.h"
	PGS_DEM_Close		"PGS_DEM.h"
	PGS_DEM_Open		"PGS_DEM.h"
	PGS_PC_GetConfigData	"PGS_PC.h"
	PGS_PC_GetNumberOfFiles	"PGS_PC.h"
	PGS_PC_GetReference	"PGS_PC.h"
	GEO_close_DEM		"GEO_earth.h"
	GEO_initialize_DEM	"GEO_earth.h"
	GEO_locate_one_granule	"GEO_main_func.h"
	GEO_read_param_file	"GEO_main_func.h"

!Revision History:
GEO_location_main.c,v
Revision 6.2  2010/08/03 20:12:34  kuyper
Changed to use 15 arc second DEM files, falling back to 30 arc second if
  necessary.
Changed to use GEO_DEM.h
Changed to expect a status code from GEO_initialize_DEM().

Revision 6.1  2010/06/18 19:54:27  kuyper
Helped resolve Bug 2470 by removing an obsolete structure and the function
  that filled it in.
Changed point of definition for many objects, to restrict their scope.

Revision 5.2  2006/01/20 16:08:55  kuyper
Reinstated Revision keyword.

Revision 5.1  2004/10/03 22:18:32  vlin
 Created a maneuver_list, filled it in with a call to GEO_read_maneuver_file,
and passed it on to GEO_locate_one_granule.
vlin@saicmodis.com

Revision 4.3  2003/12/10 22:30:56  kuyper
Corrected typo.

Revision 4.2  2003/08/13 18:12:03  kuyper
Added more information to startup message.
Changed to use imsl_wrap.h.

Revision 4.1  2002/06/11 13:33:30  vlin
Updated according to v4.1 GEO_location_main.PDL

Requirements:
		Call's PR03-CSC-2 through 34 and 46, either directly or
		indirectly, to meet requirements.

!Team-unique Header:

                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END**************************************************************************/
{
    GEO_param_struct geo_params;	/* Geolocation parameters.	*/
    int	returnStatus = EXIT_SUCCESS; /* Value to be returned by this program. */
    char msgbuf[PGSd_PC_FILE_PATH_MAX];
    char filefunc[] = __FILE__ ", main";

    modsmf(MODIS_U_GEO_BEGIN, "6.2 compiled on " __DATE__ " at "
	__TIME__ " using makefile " MAKEFILE_REVISION , filefunc);

    imsl_error_options(
	IMSL_SET_PRINT,	IMSL_WARNING,	0,
	IMSL_SET_PRINT,	IMSL_FATAL,	0,
	IMSL_SET_PRINT,	IMSL_TERMINAL,	0,
	IMSL_SET_STOP,	IMSL_FATAL,	0,
	IMSL_SET_STOP,	IMSL_TERMINAL,	0,
	0);
    
    setlinebuf (stdout);
    printf("MODIS GEO version %s, built %s %s\n",
           PROCESSVERSION, __DATE__, __TIME__);

    /* Execute GEO_read_param_file */
    if (GEO_read_param_file(&geo_params) != PGS_S_SUCCESS) {
	modsmf(MODIS_E_GEO, "GEO_read_param_file()", filefunc);
	returnStatus = EXIT_FAIL;
    }
    else
    {
	/* Resolution of Digital Elevation Model */
	PGSt_integer land_sea_layer=PGSd_DEM_WATER_LAND; /* land/sea mask tag */
	PGSt_integer numFiles;	/* Number of L1A input files to be processed.*/
	char correct_terrain[PGSd_PC_LINE_LENGTH_MAX]=TERRAIN_CORRECT;

	if (PGS_DEM_Open(DEM_resolutions, RESOLUTIONS, &land_sea_layer, 1)
	    != PGS_S_SUCCESS)
	    modsmf(MODIS_E_GEO, "PGS_DEM_Open() for land-sea mask", filefunc);

	if (PGS_PC_GetConfigData(TERRAIN_CORRECT_LUN, correct_terrain)
	    !=PGS_S_SUCCESS)
	{
	    sprintf(msgbuf,"PGS_PC_GetConfigData(%ld)",
		(long)TERRAIN_CORRECT_LUN);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);
	}
       
	if ( /* Check if terrain correction is to be performed */
	  strncmp(correct_terrain, TERRAIN_CORRECT, sizeof(TERRAIN_CORRECT))==0
	  && GEO_initialize_DEM() != PGS_S_SUCCESS)
	{
	    modsmf(MODIS_E_GEO, "GEO_initialize_DEM()", filefunc);
	    returnStatus = EXIT_FAIL;
	}
	/* Get number of level 1a input files */
	else if (PGS_PC_GetNumberOfFiles(l1aID, &numFiles)!=PGS_S_SUCCESS)
	{
	    sprintf(msgbuf,"PGS_PC_GetNumberOfFiles(%ld)", (long)l1aID);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);
	    returnStatus = EXIT_FAIL;
	}
	else
	{
	    PGSt_integer version;	/* L1A input file version.	*/

	    for (version = 1; version <= numFiles; version++)
	    {
		/* PGS_PC_GetReference changes copy_vers */
		PGSt_integer copy_vers = version; /* copy of version.	*/
		/* Input/Output file name. */
		char l1a_file_name[PGSd_PC_FILE_PATH_MAX];

		/* Get level 1A file names (Product Input) */
		if (PGS_PC_GetReference(l1aID, &copy_vers, l1a_file_name)
		    != PGS_S_SUCCESS)
		{
		    sprintf(msgbuf,"PGS_PC_GetReference(%ld,%ld)",(long)l1aID,
			   (long)version);
		    modsmf(MODIS_E_GEO, msgbuf, filefunc);
		    returnStatus = EXIT_FAIL;
		}
		else
		{
		    /* Output file name. */
		    char geo_file_name[PGSd_PC_FILE_PATH_MAX];

		    copy_vers = version;
		    /* PGS_PC_GetReference() changes copy_vers */
	       
		    /* Get geolocation file names (Product Output) */
		    if (PGS_PC_GetReference(geoID, &copy_vers, geo_file_name)
			!= PGS_S_SUCCESS)
		    {
			 sprintf(msgbuf, "PGS_PC_GetReference(%ld,%ld)",
			     (long)geoID, (long)version);
			 modsmf(MODIS_E_GEO, msgbuf, filefunc);
			 returnStatus = EXIT_FAIL;
		    }

		    /* Execute GEO_locate_one_granule */
		    else if(GEO_locate_one_granule(l1a_file_name, geo_file_name,
			&geo_params, version) != PGS_S_SUCCESS)
		    {
			sprintf(msgbuf, "GEO_locate_one_granule(%s,%s)", 
			    l1a_file_name, geo_file_name);
			modsmf(MODIS_E_GEO, msgbuf, filefunc);
			returnStatus = EXIT_FAIL;
		    }

		    modsmf(MODIS_U_GEO_GRANULE_ID, l1a_file_name, filefunc);
		} /* Else get l1a_file_name succeeded.	*/
	    }	/* for each version */
	} /* Else GetNumFiles worked. */

	if (strncmp(correct_terrain, TERRAIN_CORRECT, sizeof(TERRAIN_CORRECT))
	    ==0 && GEO_close_DEM() != PGS_S_SUCCESS)
	    modsmf(MODIS_E_GEO, "GEO_close_DEM()", filefunc);

	if (PGS_DEM_Close(DEM_resolutions, RESOLUTIONS, &land_sea_layer, 1)
	    !=PGS_S_SUCCESS)
	    modsmf(MODIS_E_GEO, "PGS_DEM_Close(land_sea_layer)", filefunc);

    } /* Else parameter file read in OK. */
    
    sprintf(msgbuf, " = %d", returnStatus);
    modsmf(MODIS_U_GEO_END, msgbuf, filefunc);
    return returnStatus;
}

