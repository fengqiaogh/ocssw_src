#include "PGS_SMF.h"
#include "version.h"
#include "GEO_output.h"
#include "PGS_MODIS_35251.h"

PGSt_SMF_status GEO_get_version_metadata(
                version_metadata_struct *version_metadata)

/***********************************************************************
!C

!Description:  Routine for loading most of the elements of version_metadata.

!Input Parameters:
	None

!Output Parameters:
	version_metadata	A structure to be filled in with the 
                                version metadata.

Return Parameters:
	MODIS_E_BAD_INPUT_ARG	If version_metadata is NULL
	MODIS_E_GEO		If either function call fails.
	PGS_S_SUCCESS		Otherwise

Externally Defined:
	MODIS_E_BAD_INPUT_ARG	"PGS_MODIS_35251.h"
	MODIS_E_GEO		"PGS_MODIS_35251.h"
	PGS_S_SUCCESS		"PGS_SMF.h"
	PROCESSVERSION		"version.h"

Called by:
	GEO_write_granule_metadata

Routines Called:
	modsmf			To log status messages.
	PGS_PC_GetConfigData	To retrieve runtime configuration parameters.

!Revision History:
$Log: GEO_get_version_metadata.c,v $
Revision 4.3  2003/08/12 19:03:30  kuyper
Corrected initialization of retval.

Revision 4.2  2003/08/07 14:50:26  kuyper
Changed handling of errors to simplify testing.

Revision 4.1  2003/08/05 18:52:54  kuyper
Added processing metadata fields.

Revision 3.1  2002/06/13 23:02:32  kuyper
Removed unnecessary NCSA acknowledgement.

Revision 1.1  2001/01/18 13:56:37  vlin
Initial revision


Requirements
	PR03-F-1-3
	PR03-F-1.3.1
	PR03-F-4.1-1
	PR03-I-1
	PR03-I-2
	PR03-S-1

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits

!END*********************************************************************/

#define LOCALVERSIONID_LUN		600001
#define PGEVERSION_LUN			800500 
#define PROCESSINGENVIRONMENT_LUN	800550
#define REPROCESSINGACTUAL_LUN		800600
#define REPROCESSINGPLANNED_LUN		800605

{
    const char *env_list[] = {"OSTYPE", "HOST", "REVISION", "MACHINE"};
#define NUM_ENVS (sizeof env_list/sizeof env_list[0])
    const char *val;
    int env;
    char msgbuf[128] = "";
    PGSt_SMF_status retval=PGS_S_SUCCESS;
    char filefunc[] = __FILE__ ", GEO_get_version_metadata";

    if (version_metadata == NULL){
	modsmf(MODIS_E_BAD_INPUT_ARG, ".", filefunc);
	return MODIS_E_BAD_INPUT_ARG;
    }

    (void)strncpy(version_metadata->processversion, PROCESSVERSION,
	sizeof(PROCESSVERSION));

    if(PGS_PC_GetConfigData(LOCALVERSIONID_LUN,version_metadata->localversionid)
	!= PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "PGS_PC_GetConfigData(%d)", LOCALVERSIONID_LUN);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
	version_metadata->localversionid[0] = '\0';
	retval = MODIS_E_GEO;
    }

    if(PGS_PC_GetConfigData(PGEVERSION_LUN,version_metadata->pgeversion)
	!= PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "PGS_PC_GetConfigData(%d)", PGEVERSION_LUN);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
	version_metadata->pgeversion[0] = '\0';
	retval = MODIS_E_GEO;
    }

    if(PGS_PC_GetConfigData(PROCESSINGENVIRONMENT_LUN,
	version_metadata->processingenvironment) != PGS_S_SUCCESS)
    {
	version_metadata->processingenvironment[0] = '\0';
	for(env=0; env<NUM_ENVS; env++)
	{
	    val = getenv(env_list[env]);
	    if(val)
	    {
		if(version_metadata->processingenvironment[0])
		    strncat(version_metadata->processingenvironment, " ",
			sizeof(version_metadata->processingenvironment));
		strncat(version_metadata->processingenvironment, val,
		    sizeof(version_metadata->processingenvironment));
	    }
	}
	version_metadata->processingenvironment
	    [sizeof(version_metadata->processingenvironment)-1] = '\0';
    }

    if(PGS_PC_GetConfigData(REPROCESSINGACTUAL_LUN,
	version_metadata->reprocessingactual) != PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "PGS_PC_GetConfigData(%d)", REPROCESSINGACTUAL_LUN);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
	version_metadata->reprocessingactual[0] = '\0';
	retval = MODIS_E_GEO;
    }

    if(PGS_PC_GetConfigData(REPROCESSINGPLANNED_LUN,
	version_metadata->reprocessingplanned) != PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "PGS_PC_GetConfigData(%d)", REPROCESSINGPLANNED_LUN);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
	version_metadata->reprocessingplanned[0] = '\0';
	retval = MODIS_E_GEO;
    }

    return retval;
}
