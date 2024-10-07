#include "SDST_TK.h"
#include "GEO_main.h"
#include "GEO_output.h"
#include "GEO_product.h"
#include "GEO_geo.h"
#include "L1a_data.h"
#include "GEO_parameters.h"
#include "PGS_MODIS_35251.h"

PGSt_SMF_status GEO_write_ECS_metadata(
        MODFILE                 * const geo_file,
	ECS_metadata_struct	* const ECS_metadata,
        EPH_metadata_struct     * const EPH_metadata,
	GEO_bcoord_struct       * const bounding_coords,
	GEO_GRing_struct        * const GRing_points,
	pointer_metadata_struct * const pointer_metadata,
	qa_metadata_struct      * const qa_metadata,
        char                    * const sci_state,
        char                    * const sci_abnorm
	)
/*
!C******************************************************************************
!Description:   
	Routine in Output group of the Level-1A geolocation software to set ECS
	Inventory and Archive metadata in MCF memory.  The metadata are written
	to the geolocation product just before the end of the routine.

	Warning messages will be recorded (to the SMF) if metadata inputs are
	missing, but the routine will not abort, which may result in
	geolocation products with incomplete or missing ECS metadata.
		
!Input Parameters:
	geo_file		the MAPI structure for the product
	ECS_metadata		ECS Core Metadata from l1a file
	EPH_metadata		Spacecraft ephemeris metadata
	bounding_coords		granule's bounding coordinates
	GRing_points		granule's GRing vertice's data
	pointer_metadata	universal references of geolocation data
				    sources and sinks
	qa_metadata		Granule-level quality assurance metadata.
	sci_abnorm		spacecraft in unusual state during granule
				    (i.e. maneuver) flag
				    "0" = maneuver
				    "1" = normal
	sci_state		instrument test mode during granule flag
				    "0" = MODIS test mode
				    "1" = normal

!Output Parameters:
	none

Return parameter:
	MODIS_E_BAD_INPUT_ARG	If geo_file is null (the other pointer
				    arguments are all allowed to be null)
	MODIS_E_GEO		If any subroutine failed.
	PGS_S_SUCCESS		Otherwise

Externally Defined:
	ADDITIONALATTRIBUTENAME		"GEO_product.h"
	ARCHIVEMETADATA			"GEO_product.h"
	ARCHIVED_METADATA		"mapi.h"
	AUTOMATICQUALITYFLAG		"GEO_product.h"
	AUTOMATICQUALIFYFLAGEXPLANATION	"GEO_product.h"
	COREMETADATA			"GEO_product.h"
	DAYNIGHTFLAG			"GEO_product.h"
	DOIAUTHORITY			"GEO_product.h"
	EASTBOUNDINGCOORDINATE		"GEO_product.h"
	EQUATORCROSSINGLONGITUDE	"GEO_product.h"
	EQUATORCROSSINGDATE		"GEO_product.h"
	EQUATORCROSSINGTIME		"GEO_product.h"
	EXCLUSIONGRINGFLAG		"GEO_product.h"
	geo_MCFID			"GEO_main.h"
	GRINGPOINTLATITUDE		"GEO_product.h"
	GRINGPOINTLONGITUDE		"GEO_product.h"
	GRINGPOINTSEQUENCENO		"GEO_product.h"
	IDENTIFIERPRODDOI		"GEO_product.h"
	IDENTIFIERPRODDOIAUTH		"GEO_product.h"
	INPUTPOINTER			"GEO_product.h"
	INVENTORY_METADATA		"mapi.h"
	LOCALGRANULEID			"GEO_product.h"
	LOCALINPUTGRANULEID		"GEO_product.h"
	LOCALVERSIONID			"GEO_product.h"
	MODIS_E_BAD_INPUT_ARG		"PGS_MODIS_35251.h"
	MODIS_E_GEO			"PGS_MODIS_35251.h"
	NORTHBOUNDINGCOORDINATE		"GEO_product.h"
	ORBITNUM			"GEO_product.h"
	PARAMETERNAME			"GEO_product.h"
	PARAMETERVALUE			"GEO_product.h"
	PGS_S_SUCCESS			"PGS_SMF.h"
	PGEVERSION			"GEO_product.h"
	PROCESSINGENVIRONMENT		"GEO_product.h"
	PRODUCTIONHISTORY		"GEO_product.h"
	QAPERCENTMISSINGDATA		"GEO_product.h"
	QAPERCENTOUTOFBOUNDSDATA	"GEO_product.h"
	RANGEBEGINNINGDATE		"GEO_product.h"
	RANGEBEGINNINGTIME		"GEO_product.h"
	RANGEENDINGDATE			"GEO_product.h"
	RANGEENDINGTIME			"GEO_product.h"
	REPROCESSINGACTUAL		"GEO_product.h"
	REPROCESSINGPLANNED		"GEO_product.h"
	SOUTHBOUNDINGCOORDINATE		"GEO_product.h"
	WESTBOUNDINGCOORDINATE		"GEO_product.h"

Called by:
	GEO_write_granule_metadata

Routines Called:
	modsmf				"smfio.h"
	PGS_MET_GetSetAttr		"PGS_MET.h"
	PGS_MET_Init			"PGS_MET.h"
	PGS_MET_Remove			"PGS_MET.h"
	PGS_MET_SetAttr			"PGS_MET.h"
	PGS_MET_Write			"PGS_MET.h"
	putMODISfileinfo		"mapi.h"

!Revision History:
$Log: GEO_write_ECS_metadata.c,v $
Revision 6.3  2014/10/22 21:59:23  jkuyper
Changed to insert ".NRT" in DOI metadata if running on the Near Real Time
  system.

Revision 6.2  2012/06/28 21:42:15  kuyper
Corrected format specifier for versionid.

Revision 6.1  2012/06/27 19:16:18  kuyper
Changed to write DOI metadata to Geolocation product.

Revision 5.1  2004/06/23 20:44:42  vlin
QAPercent items set to 100 when number of pixels is zero.

Revision 4.2  2003/08/20 18:02:14  kuyper
Corrected to handle the increased number of input pointers.

Revision 4.1  2003/08/12 20:46:02  kuyper
Changed to write ReprocessingActual from version_metadata.
Added ReprocessingPlanned and ProcessingEnvironment.
Change to write ProductionHistory based upon the PGE and PGE version,
  rather than upon the process and the LocalVersionID.
Changed to return status codes.

Revision 3.2  2002/07/08 21:32:52  kuyper
Change to use SDST toolkit to write LocaleGranuleID.

Revision 3.1  2001/04/07 22:56:43  kuyper
Changed to check L1A productionhistory for "NOT SET".
Made productionhistory length checks more flexible.

James Kuyper - kuyper@ltpmail.gsfc.nasa.gov

Requirements:
	PR03-F-4.1-1
	PR03-F-4.1-2

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.
!END****************************************************************************
*/

{
#define NUM_GRING_PTS 4
#define NUM_IVALS (NUM_GRING_PTS + 3)
#define NUM_DVALS (2 * NUM_GRING_PTS + 5)
#define NUM_STRVALS (36 + INPUTPOINTERS)

    PGSt_integer ivals[NUM_IVALS];
    PGSt_double dvals[NUM_DVALS];
    char *strings[NUM_STRVALS];
    int obj=0, str=0, ival=0, dval=0;
    int i;
    PGSt_SMF_status retval = PGS_S_SUCCESS;
    PGSt_MET_all_handles  mdHandles;
    ECSattr_names_for_all_handles HDFattrNames;
    char productionhistory[64] = "";
    char shortname[9], *pshort=shortname;
    PGSt_integer versionid;
    char doi[64];
    char msgbuf[128] = "";
    char *valid_flag = msgbuf;
    /* value of valid_flag is not important, so long as not NULL */
    char filefunc[] = __FILE__ ", GEO_write_ECS_metadata";

    struct {
	char *attrNameStr;
	/* handle number of the metadata group to write the parameter to */
	int   metadata_group;
	void *attrValue;
	/* Note: only the following are permitted: PGSt_integer *attrValue,
	 * PGSt_double *attrValue, PGSt_real *attrValue, and char **attrValue
	 */
    } MCF_object[NUM_IVALS + NUM_DVALS + NUM_STRVALS];

    /* handle numbers to the metadata groups accessed through mdHandles */
#define INVENTORY_METADATA 1
#define ARCHIVED_METADATA  2
#define min(x,y) ((x) < (y) ? (x) : (y))

#define STROBJ(name, handle, base, value) { MCF_object[obj].attrNameStr=name; \
    MCF_object[obj].metadata_group = handle;  \
    if(base!=NULL){	\
	strings[str] = value; MCF_object[obj++].attrValue = &strings[str++]; \
    } else MCF_object[obj++].attrValue = NULL;}
      
#define INTOBJ(name, handle, base, value) { MCF_object[obj].attrNameStr=name; \
    MCF_object[obj].metadata_group = handle;  \
    if(base!=NULL){	\
	ivals[ival] = (PGSt_integer)(value);	\
	MCF_object[obj++].attrValue = &ivals[ival++];	\
    } else MCF_object[obj++].attrValue = NULL;}
      
#define DBLOBJ(name, handle, base, value) { MCF_object[obj].attrNameStr=name; \
    MCF_object[obj].metadata_group = handle;  \
    if(base!=NULL){	\
	dvals[dval] = (PGSt_double)(value);	\
	MCF_object[obj++].attrValue = &dvals[dval++];	\
    } else MCF_object[obj++].attrValue = NULL;}

    if( geo_file == NULL)
    {
	modsmf(MODIS_E_BAD_INPUT_ARG,".", filefunc);
	return MODIS_E_BAD_INPUT_ARG;
    }

    sprintf(productionhistory, "PGE01:%.*s",
	(int)(sizeof(productionhistory)-sizeof("PGE01:")),
	ECS_metadata->version_metadata.pgeversion);

    /* WARNING: do not add any *OBJ() macros below, without first checking that
      strings[], ivals[], or dvals[] is big enough.	*/

    STROBJ(PGEVERSION, INVENTORY_METADATA, valid_flag,
	ECS_metadata->version_metadata.pgeversion);
    STROBJ(LOCALVERSIONID, INVENTORY_METADATA, valid_flag,
	ECS_metadata->version_metadata.localversionid);
    STROBJ(RANGEBEGINNINGDATE, INVENTORY_METADATA, ECS_metadata, 
	ECS_metadata->rangebeginningdate);
    STROBJ(RANGEBEGINNINGTIME, INVENTORY_METADATA, ECS_metadata, 
	ECS_metadata->rangebeginningtime);
    STROBJ(RANGEENDINGDATE, INVENTORY_METADATA, ECS_metadata, 
	ECS_metadata->rangeendingdate);
    STROBJ(RANGEENDINGTIME, INVENTORY_METADATA, ECS_metadata, 
	ECS_metadata->rangeendingtime);
    STROBJ(DAYNIGHTFLAG, INVENTORY_METADATA, ECS_metadata, 
	ECS_metadata->operationmode);
    STROBJ(EQUATORCROSSINGDATE ".1", INVENTORY_METADATA, EPH_metadata, 
	EPH_metadata->equatorcrossingdate);
    STROBJ(EQUATORCROSSINGTIME ".1", INVENTORY_METADATA, EPH_metadata, 
	EPH_metadata->equatorcrossingtime);
    STROBJ(EXCLUSIONGRINGFLAG, INVENTORY_METADATA, GRing_points, 
	GRing_points->exclusionflag);
    STROBJ(INPUTPOINTER, INVENTORY_METADATA, pointer_metadata, 
	pointer_metadata->inputpointer[0]);
    for(i=1; i<INPUTPOINTERS && pointer_metadata->inputpointer[i][0]; i++)
	strings[str++] =  pointer_metadata->inputpointer[i];
    if(i<INPUTPOINTERS)
	strings[str++] = NULL;
    STROBJ(REPROCESSINGACTUAL, INVENTORY_METADATA, ECS_metadata,
	ECS_metadata->version_metadata.reprocessingactual);
    STROBJ(REPROCESSINGPLANNED, INVENTORY_METADATA, ECS_metadata,
	ECS_metadata->version_metadata.reprocessingplanned);

    STROBJ(LOCALINPUTGRANULEID, ARCHIVED_METADATA, ECS_metadata, 
	ECS_metadata->localinputgranuleid);
    STROBJ(MECS_PRODHISTORY, ARCHIVED_METADATA, ECS_metadata,
	productionhistory);
    STROBJ(PROCESSINGENVIRONMENT, ARCHIVED_METADATA, ECS_metadata,
	ECS_metadata->version_metadata.processingenvironment);
    DBLOBJ(EASTBOUNDINGCOORDINATE, ARCHIVED_METADATA, bounding_coords, 
	bounding_coords->eastcoord*RAD2DEG);
    DBLOBJ(NORTHBOUNDINGCOORDINATE, ARCHIVED_METADATA, bounding_coords,
	bounding_coords->northcoord*RAD2DEG);
    DBLOBJ(SOUTHBOUNDINGCOORDINATE, ARCHIVED_METADATA, bounding_coords,
	bounding_coords->southcoord * RAD2DEG);
    DBLOBJ(WESTBOUNDINGCOORDINATE, ARCHIVED_METADATA, bounding_coords,
	bounding_coords->westcoord * RAD2DEG);
    DBLOBJ(EQUATORCROSSINGLONGITUDE ".1", INVENTORY_METADATA, EPH_metadata,
	EPH_metadata->equatorcrossinglongitude * RAD2DEG);

  /* The following are arrays, so will require different treatment. */
    if(GRing_points==NULL)
    {
	MCF_object[obj].attrNameStr = GRINGPOINTLATITUDE;
	MCF_object[obj].metadata_group = INVENTORY_METADATA;
	MCF_object[obj++].attrValue = NULL;

	MCF_object[obj].attrNameStr = GRINGPOINTLONGITUDE;
	MCF_object[obj].metadata_group = INVENTORY_METADATA;
	MCF_object[obj++].attrValue = NULL;

	MCF_object[obj].attrNameStr = GRINGPOINTSEQUENCENO;
	MCF_object[obj].metadata_group = INVENTORY_METADATA;
	MCF_object[obj++].attrValue = NULL;
    }
    else
    {
	MCF_object[obj].attrNameStr = GRINGPOINTLATITUDE;
	MCF_object[obj].metadata_group = INVENTORY_METADATA;
	MCF_object[obj++].attrValue = &dvals[dval];
	for(i=0; i<NUM_GRING_PTS ; i++)
	    dvals[dval++] = (PGSt_double)GRing_points->latitude[i];

	MCF_object[obj].attrNameStr = GRINGPOINTLONGITUDE;
	MCF_object[obj].metadata_group = INVENTORY_METADATA;
	MCF_object[obj++].attrValue = &dvals[dval];
	for(i=0; i<NUM_GRING_PTS ; i++)
	    dvals[dval++] = (PGSt_double)GRing_points->longitude[i];

	MCF_object[obj].attrNameStr = GRINGPOINTSEQUENCENO;
	MCF_object[obj].metadata_group = INVENTORY_METADATA;
	MCF_object[obj++].attrValue = &ivals[ival];
	for(i=0; i<NUM_GRING_PTS ; i++)
	    ivals[ival++] = (PGSt_integer)GRing_points->sequenceno[i];
    }

    INTOBJ(ORBITNUMBER".1", INVENTORY_METADATA, EPH_metadata, 
	EPH_metadata->orbitnumber);
    STROBJ(ADDITIONALATTRIBUTENAME".1", INVENTORY_METADATA, valid_flag,
	"GRANULENUMBER");
    STROBJ(PARAMETERVALUE".1", INVENTORY_METADATA, ECS_metadata, 
	ECS_metadata->granulenumber);
    STROBJ(ADDITIONALATTRIBUTENAME".2", INVENTORY_METADATA, valid_flag,
	"SCI_STATE");
    STROBJ(PARAMETERVALUE".2", INVENTORY_METADATA, sci_state,
	sci_state);
    STROBJ(ADDITIONALATTRIBUTENAME".3", INVENTORY_METADATA, valid_flag,
	"SCI_ABNORM");
    STROBJ(PARAMETERVALUE".3", INVENTORY_METADATA, sci_abnorm,
	sci_abnorm);
    STROBJ(PARAMETERVALUE".5", INVENTORY_METADATA, ECS_metadata,
	ECS_metadata->version_metadata.processversion);

    /*  added elements.  */
    
    if(qa_metadata==NULL){
	STROBJ(AUTOMATICQUALITYFLAG".1", INVENTORY_METADATA, valid_flag,
	    "Failed");
	INTOBJ(QAPERCENTMISSINGDATA".1", INVENTORY_METADATA, valid_flag, 0);
	INTOBJ(QAPERCENTOUTOFBOUNDSDATA".1", INVENTORY_METADATA, valid_flag, 0);
	STROBJ(PARAMETERVALUE".4", INVENTORY_METADATA, valid_flag, "");
    }
    else
    {
	if(qa_metadata->retval == SUCCESS)
	    STROBJ(AUTOMATICQUALITYFLAG".1", INVENTORY_METADATA, valid_flag,
		"Passed")
	else
	    STROBJ(AUTOMATICQUALITYFLAG".1", INVENTORY_METADATA, valid_flag,
		"Failed")

	if(qa_metadata->no_of_pixels != 0)
	{
	    INTOBJ(QAPERCENTMISSINGDATA".1", INVENTORY_METADATA, valid_flag, 
		(int)floor(0.5+100.0*(double)qa_metadata->missingdata/
		(double)qa_metadata->no_of_pixels));
	    INTOBJ(QAPERCENTOUTOFBOUNDSDATA".1", INVENTORY_METADATA, valid_flag,
		(int)floor(0.5+100.0*(double)qa_metadata->outofboundsdata/
		(double)qa_metadata->no_of_pixels));
	}else {
	    INTOBJ(QAPERCENTMISSINGDATA".1", INVENTORY_METADATA, valid_flag, 
                   100);
	    INTOBJ(QAPERCENTOUTOFBOUNDSDATA".1", INVENTORY_METADATA, valid_flag,
		   100);
	}
	STROBJ(PARAMETERVALUE".4", INVENTORY_METADATA, valid_flag,
	    qa_metadata->rms_error);
    }
    
    STROBJ(ADDITIONALATTRIBUTENAME".4", INVENTORY_METADATA, valid_flag,
	"GEO_EST_RMS_ERROR");
    STROBJ(ADDITIONALATTRIBUTENAME".5", INVENTORY_METADATA, valid_flag,
	"PROCESSVERSION");

    STROBJ(ADDITIONALATTRIBUTENAME".6", INVENTORY_METADATA, valid_flag,
	IDENTIFIERPRODDOI);
    STROBJ(PARAMETERVALUE".6", INVENTORY_METADATA, valid_flag, doi);

    STROBJ(ADDITIONALATTRIBUTENAME".7", INVENTORY_METADATA, valid_flag,
	IDENTIFIERPRODDOIAUTH);
    STROBJ(PARAMETERVALUE".7", INVENTORY_METADATA, valid_flag, DOIAUTHORITY);
	
    STROBJ(PARAMETERNAME".1", INVENTORY_METADATA, valid_flag, "Geolocation");
    STROBJ(AUTOMATICQUALITYFLAGEXPLANATION".1", INVENTORY_METADATA, valid_flag,
    "Set to 'Failed' if processing error occurred, set to 'Passed' otherwise");

    if(PGS_MET_Init(geo_MCFID, mdHandles) != PGS_S_SUCCESS)
    {
      sprintf(msgbuf, "PGS_MET_Init(%ld)", (long)geo_MCFID);
      modsmf(MODIS_E_GEO, msgbuf, filefunc);
      return MODIS_E_GEO;
    }

    if(PGS_MET_GetSetAttr(mdHandles[INVENTORY_METADATA], SHORTNAME, &pshort)
	!= PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "PGS_MET_GetSetAttr(\"%s\",\"" SHORTNAME "\")",
	    mdHandles[INVENTORY_METADATA]);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
	return MODIS_E_GEO;
    }

    if(PGS_MET_GetSetAttr(mdHandles[INVENTORY_METADATA], VERSIONID, &versionid)
	!= PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "PGS_MET_GetSetAttr(\"%s\",\"" VERSIONID "\")",
	    mdHandles[INVENTORY_METADATA]);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
	return MODIS_E_GEO;
    }

    sprintf(doi, "10.5067/MODIS/%s%s.%03ld", shortname,
	strcmp(ECS_metadata->version_metadata.reprocessingactual,
	"Near Real Time") ? "" : ".NRT", (long)versionid);
    if(putMODISfileinfo(geo_file, IDENTIFIERPRODDOI, TXT, strlen(doi), doi)
	!= MAPIOK)
    {
	sprintf(msgbuf, "putMODISfileinfo(\"%s\",\"" IDENTIFIERPRODDOI "\", \""
	    TXT ", %d, \"%s\")", geo_file->filename, (int)strlen(doi), doi);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
	return MODIS_E_GEO;
    }

    if(putMODISfileinfo(geo_file, IDENTIFIERPRODDOIAUTH, TXT,
	sizeof(DOIAUTHORITY)-1, DOIAUTHORITY)
	!= MAPIOK)
    {
	sprintf(msgbuf, "putMODISfileinfo(\"%s\",\"" IDENTIFIERPRODDOIAUTH
	    "\", \"" TXT ", %d, \"" DOIAUTHORITY "\")",
	    geo_file->filename, (int)sizeof(DOIAUTHORITY)-1);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
	return MODIS_E_GEO;
    }
					  /* Set each parameter with its value 
					     in the designated metadata group */
    for(i=0; i<obj; i++)
    {
      if(MCF_object[i].attrValue == NULL)
	modsmf(MODIS_W_MISSING_OUTPUT, MCF_object[i].attrNameStr, filefunc);
      else if(PGS_MET_SetAttr(mdHandles[MCF_object[i].metadata_group],
	MCF_object[i].attrNameStr, MCF_object[i].attrValue) != PGS_S_SUCCESS)
      {
	sprintf(msgbuf, "PGS_MET_SetAttr(%s)", MCF_object[i].attrNameStr);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
	retval = MODIS_E_GEO;
      }
    } /* END FOR each MCF_object */

    if(SDST_SetLocalGranId(GEO, mdHandles)!=MODIS_S_SDST_SUCCESS)
    {
	modsmf(MODIS_E_GEO, "SDST_SetLocalGranId", filefunc);
	retval =  MODIS_E_GEO;
    }

    (void)strncpy(HDFattrNames[INVENTORY_METADATA], COREMETADATA,
	  sizeof(HDFattrNames[0]));
    (void)strncpy(HDFattrNames[ARCHIVED_METADATA], ARCHIVEMETADATA,
	  sizeof(HDFattrNames[0]));

    if(PGS_MET_Write(mdHandles[INVENTORY_METADATA], 
       HDFattrNames[INVENTORY_METADATA], (PGSt_integer)geo_file->sd_id)
       !=  PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "PGS_MET_Write(%s)",COREMETADATA);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
	retval = MODIS_E_GEO;
    }

    if(PGS_MET_Write(mdHandles[ARCHIVED_METADATA],
	HDFattrNames[ARCHIVED_METADATA], (PGSt_integer)geo_file->sd_id)
	!= PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "PGS_MET_Write(%s)",ARCHIVEMETADATA);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
	retval = MODIS_E_GEO;
    }

    PGS_MET_Remove();

    return retval;
}
