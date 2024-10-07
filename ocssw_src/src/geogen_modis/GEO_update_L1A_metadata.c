#include "PGS_MET.h"
#include "GEO_product.h"
#include "GEO_main.h"
#include "GEO_output.h"
#include "GEO_geo.h"
#include "L1a_data.h"

PGSt_SMF_status GEO_update_L1A_metadata(
			   MODFILE 	      * const l1a_file,
			   GEO_GRing_struct   * const GRing_points,
                           EPH_metadata_struct * const EPH_metadata,
                           GEO_bcoord_struct   * const bounding_coords
			   )

/*
!C******************************************************************************
!Description:   
	Routine in Output group of the Level-1A geolocation software to write
	the geolocated metadata to the L1A product file.  Writing GRing equator
	crossing, and orbitnumber metadata to the ECS Core Inventory Metadata
	requires re-writing the entire Core Metadata content back to the file.

	Writes ECS Inventory Metadata:
	    ECSDataGranule
		LOCALGRANULEID
		DAYNIGHTFLAG
		REPROCESSINGACTUAL
		LOCALVERSIONID
		REPROCESSINGPLANNED
	    OrbitCalculatedSpatialDomain
		OrbitCalculatedSpatialDomainContainer
		    EQUATORCROSSINGDATE
		    EQUATORCROSSINGTIME
		    ORBITNUMBER
		    EQUATORCROSSINGLONGITUDE
	    InputGranule
		INPUTPOINTER
	    SpatialDomainContainer
		HorizontalSpatialDomainContainer
		    GPolygon
			GPolygonContainer
			    GRingPoint
				GRINGPOINTLATITUDE
				GRINGPOINTLONGITUDE
				GRINGPOINTSEQUENCENO
			    EXCLUSIONGRINGFLAG
	    RangeDateTime
		RANGEENDINGDATE
		RANGEENDINGTIME
		RANGEBEGINNINGDATE
		RANGEBEGINNINGTIME
	    AdditionalAttributes
		AdditionalAttributesContainer
		    ADDITIONALATTRIBUTENAME.1
		    ADDITIONALATTRIBUTENAME.2
		    ADDITIONALATTRIBUTENAME.3
		    ADDITIONALATTRIBUTENAME.4
		    PARAMETERVALUE.1
		    PARAMETERVALUE.2
		    PARAMETERVALUE.3
		    PARAMETERVALUE.4
	    PGEVersionClass
		PGEVERSION

        Writing the East, North, South, and West BoundingCoordinate fields to
        the Archived Metadata requires copying over all of the other
        Source = "PGE" archived metadata:
            ProductionHistory
            ProcessingEnvironment



!Input Parameters:
	l1a_file	MAPI structure for Level 1A file
	GRing_points	G-Ring polygon latitudes and longitudes
	EPH_metadata	Spacecraft ephemeris metadata.
        bounding_coords Bounding values for granule in latitude and longitude

!Output Parameters:
	None

Return values:
	MODIS_E_BAD_INPUT_ARG		If any input pointer is null
	MODIS_E_GEO			If certain subroutines fail
	PGS_S_SUCCESS			Otherwise

Externally Defined:
	ADDITIONALATTRIBUTENAME		"GEO_product.h"
	DATATYPELENMAX			"mapi.h"
	DAYNIGHTFLAG			"GEO_product.h"
	DOIAUTHORITY			"GEO_product.h"
	EQUATORCROSSINGDATE		"GEO_product.h"
	EQUATORCROSSINGLONGITUDE	"GEO_product.h"
	EQUATORCROSSINGTIME		"GEO_product.h"
	EXCLUSIONGRINGFLAG		"GEO_product.h"
	GRINGPOINTLATITUDEG		"GEO_product.h"
	GRINGPOINTLONGITUDEG		"GEO_product.h"
	GRINGPOINTSEQUENCENO 		"GEO_product.h"
	IDENTIFIERPRODDOI		"GEO_product.h"
	IDENTIFIERPRODDOIAUTH		"GEO_product.h"
	INPUTPOINTER			"GEO_product.h"
	LOCALGRANULEID			"GEO_product.h"
	LOCALVERSIONID			"GEO_product.h"
	MAPIOK				"mapi.h"
        MECS_ARCHIVE                    "mapi.h"
	MECS_CORE			"mapi.h"
        MECS_PRODHISTORY                "L1a_data.h"
	MODIS_E_BAD_INPUT_ARG		"PGS_MODIS_35251.h"
	MODIS_E_GEO			"PGS_MODIS_35251.h"
	MODIS_W_MISSING_OUTPUT		"PGS_MODIS_35251.h"
	ORBITNUMBER			"GEO_product.h"
	PARAMETERVALUE			"GEO_product.h"
	PGS_S_SUCCESS			"PGS_SMF.h"
	PGEVERSION			"GEO_product.h"
        PROCESSINGENVIRONMENT           "GEO_product.h"
	PRODUCTIONDATETIME		"GEO_product.h"
	RAD2DEG				"GEO_geo.h"
	RANGEBEGINNINGDATE		"GEO_product.h"
	RANGEBEGINNINGTIME		"GEO_product.h"
	RANGEENDINGDATE			"GEO_product.h"
	RANGEENDINGTIME			"GEO_product.h"
	REPROCESSINGPLANNED		"GEO_product.h"
	REPROCESSINGACTUAL		"GEO_product.h"
	TXT				"mapi.h"

Routines called:
	getMODISECSinfo			"mapi.h"
	MODISsizeof			"mapi.h"
	modsmf				"smfio.h"
	PGS_MET_GetSetAttr		"PGS_MET.h"
	PGS_MET_Init			"PGS_MET.h"
	PGS_MET_Remove			"PGS_MET.h"
	PGS_MET_SetAttr			"PGS_MET.h"
	PGS_MET_Write			"PGS_MET.h"
	putMODISfileinfo		"mapi.h"
	substrMODISECSinfo		"mapi.h"

Called by:
	GEO_write_granule_metadata

!Revision History:
$Log: GEO_update_L1A_metadata.c,v $
Revision 6.4  2014/12/12 21:45:13  jkuyper
Corrected check for NRT processing.

Revision 6.3  2014/10/23 20:21:51  jkuyper
Corrected assignment operator to equality operator.

Revision 6.2  2014/10/22 21:57:53  jkuyper
Changed to insert ".NRT" in DOI metadata of L1A file if it was produced on
  the Near Real Time system.

Revision 6.1  2012/06/28 21:40:25  kuyper
Changed to write DOI metadata to the L1A product.

Revision 5.4  2004/11/18 22:26:24  kuyper
Corrected use of NULL attribute arguments for getMODISECSinfo(), and handling
  of error returns.

Revision 5.3  2004/11/05 01:05:16  kuyper
Corrected to write bounding coordinate values in degrees, not radians.

Revision 5.2  2004/10/20 19:14:52  seaton
Made walkthrough changes
,

Revision 5.1  2004/10/07 18:55:20  seaton
Added metadata for Collection 5 processing.
This includes Bounding Coordinates, ProductionHistory and ProcessingEnvironment.
seaton@saicmodis.com

Revision 4.5  2003/11/12 18:31:15  kuyper
Updated MAX_STRINGS to accomodate excess value of NUM_VALS for InputPointer
in MYD01.mcf.

Revision 4.4  2003/10/27 01:24:47  vlin
message buffer size expanded.

Revision 4.3  2003/10/01 18:19:53  kuyper
Changed to update missing metadata with "Missing value".

Revision 4.2  2003/08/01 21:37:56  kuyper
Corrected several typos.

Revision 4.1  2003/07/28 22:10:02  kuyper
Added new metadata fields.
Changed to return a status code.
Cleaned up messaging.
Simplified by recognising that all fields currently used have a data type
  of TXT.


Requirements:
		PR03-F-4.1-1
		PR03-F-4.1-2

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END****************************************************************************
*/

#define INVENTORY_METADATA 1

/* MAX_STRINGS is set to one more than the maximum number of INPUTPOINTERs
  documented in the L1A product spec. to ensure that the INPUTPOINTER
  metadata are terminated with a NULL pointer. */

#define MAX_STRINGS 6
#define META_BUF_SIZE ((MAX_STRINGS-1) * PGSd_PC_UREF_LENGTH_MAX)

{
    char shortname[9], *pshort=shortname;
    PGSt_integer versionid;
    char *pident = IDENTIFIERPRODDOI;
    char doi[64], *pdoi = doi;
    char *pident_auth = IDENTIFIERPRODDOIAUTH;
    char *pdoi_auth = DOIAUTHORITY;
    const char *NRT_infix = ".NRT";

    struct {
        char *attribute;
        char *object;

    } echoed_metadata[] = {
	{MECS_CORE, ADDITIONALATTRIBUTENAME".1"},
	{MECS_CORE, ADDITIONALATTRIBUTENAME".2"},
	{MECS_CORE, DAYNIGHTFLAG},
	{MECS_CORE, INPUTPOINTER},
	{MECS_CORE, LOCALGRANULEID},
	{MECS_CORE, LOCALVERSIONID},
	{MECS_CORE, PARAMETERVALUE".1"},
	{MECS_CORE, PARAMETERVALUE".2"},
	{MECS_CORE, PGEVERSION},
	{MECS_CORE, PRODUCTIONDATETIME},
	{MECS_CORE, RANGEBEGINNINGDATE},
	{MECS_CORE, RANGEBEGINNINGTIME},
	{MECS_CORE, RANGEENDINGDATE},
	{MECS_CORE, RANGEENDINGTIME},
#define ACTUAL_META 14
	{MECS_CORE, REPROCESSINGACTUAL},
	{MECS_CORE, REPROCESSINGPLANNED},
#define ARCHIVE_META 16
        {MECS_ARCHIVE, MECS_PRODHISTORY},
        {MECS_ARCHIVE, PROCESSINGENVIRONMENT}
    };

    struct {
        char *attribute;
	char *object;
	void *data;
    }geo_metadata[] = {
#define EQUAT_META 0	/* Starting index for equator metadata.	*/
	{MECS_CORE,	EQUATORCROSSINGLONGITUDE".1",	NULL},
	{MECS_CORE,	ORBITNUMBER".1",		NULL},
	{MECS_CORE,	EQUATORCROSSINGDATE".1",	NULL},
	{MECS_CORE,	EQUATORCROSSINGTIME".1",	NULL},
#define GRING_META 4	/* Starting index for GRing metadata.	*/
	{MECS_CORE,	EXCLUSIONGRINGFLAG,		NULL},
	{MECS_CORE,	GRINGPOINTLATITUDE,		NULL},
	{MECS_CORE,	GRINGPOINTLONGITUDE,		NULL},
	{MECS_CORE,	GRINGPOINTSEQUENCENO,		NULL},
#define BOUND_META 8    /* starting index for bounding coords metadata */ 
        {MECS_ARCHIVE,	EASTBOUNDINGCOORDINATE,		NULL},
        {MECS_ARCHIVE,	WESTBOUNDINGCOORDINATE,		NULL},
        {MECS_ARCHIVE,	NORTHBOUNDINGCOORDINATE,	NULL},
        {MECS_ARCHIVE,	SOUTHBOUNDINGCOORDINATE,	NULL},
#define DOI_META 12
	{MECS_CORE,	ADDITIONALATTRIBUTENAME".3",	&pident},
	{MECS_CORE,	PARAMETERVALUE".3",		&pdoi},
	{MECS_CORE,	ADDITIONALATTRIBUTENAME".4",	&pident_auth},
	{MECS_CORE,	PARAMETERVALUE".4",		&pdoi_auth}
    };

    PGSt_MET_all_handles mdHandles;
    const char* lastattribute="";
    /* to transfer metadata from the L1A to the MCF file */
    /* needed for proper alignment              */
    int retval = PGS_S_SUCCESS; /* routine state return value */
    int i, handle;
    char *metadata_strings[MAX_STRINGS]={NULL};
    char msgbuf[256];
    char filefunc[] = __FILE__ ", GEO_update_L1A_metadata";

    if(l1a_file == NULL) {
	modsmf(MODIS_E_BAD_INPUT_ARG, "", filefunc);
	return MODIS_E_BAD_INPUT_ARG;
    }
   
    if(PGS_MET_Init(l1a_MCFID, mdHandles) != PGS_S_SUCCESS)
    {
	modsmf(MODIS_E_GEO, "PGS_MET_Init()", filefunc);
	return MODIS_E_GEO;
    }

    for(i=0; i<(int)(sizeof(echoed_metadata)/sizeof(echoed_metadata[0])); i++)
    {
	char datatype[DATATYPELENMAX];
	float64 metadata_buffer[(META_BUF_SIZE/sizeof(float64)) + 1];
	long n_strings;
	long n_elements = META_BUF_SIZE / MODISsizeof(TXT);

	/* Pass NULL for i!=0 to avoid rebuilding of PVL table */
	strncpy(datatype, TXT, sizeof(datatype));
	if(getMODISECSinfo(l1a_file, strncmp(echoed_metadata[i].attribute,
	    lastattribute, sizeof echoed_metadata[i].attribute) ?
	    echoed_metadata[i].attribute : NULL,
	    echoed_metadata[i].object, datatype, &n_elements, metadata_buffer)
	    != MAPIOK)
	{
	    sprintf(msgbuf, "getMODISECSinfo(\"%s\",\"%s\",\"%s\")",
	        l1a_file->filename, echoed_metadata[i].attribute, echoed_metadata[i].object);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);

	    /* It is not an error for non-mandatory metadata to be missing, and
	     * PGS_MET_Write() will take care of complaining about the absence
	     * of mandatory metadata. Therefore, we fail at this point only if
	     * there's some other kind of problem.
	     */
	    if(n_elements!=0 && n_elements != META_BUF_SIZE / MODISsizeof(TXT))
		retval = MODIS_E_GEO;
		/* getMODISECSinfo has failed for some reason other than
				not finding echoed_metadata.attrNameStr */

	    /* Setting an invalid value at this point allows PGE01B 4.1.0+ to
	     * run on input MOD01 files created by earlier versions of PGE01.
	     * The output MOD01 file would fail to ingest due to not having
	     * a valid value for the metadata, but the MOD01 file isn't
	     * ingested from PGE01B.
	     */
	    strncpy((char*)metadata_buffer, "Missing value",
		sizeof(metadata_buffer));
	    n_elements = (long)strlen((char*)metadata_buffer)+1L;
	}
	else
	    lastattribute = echoed_metadata[i].attribute;
	/* Break up the text metadata into its component strings before
	  re-writing it to the MCF: */
	n_strings = MAX_STRINGS;

	if(substrMODISECSinfo((char *)metadata_buffer, n_elements,
	    &n_strings, metadata_strings)!= MAPIOK)
	{
	    sprintf(msgbuf, "substrMODISECSinfo(\"%s\")",
		echoed_metadata[i].object);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);
	    retval = MODIS_E_GEO;

	}
	else
	{
	    int str;
	    for(str=(int)n_strings; str<MAX_STRINGS; str++)
		metadata_strings[str] = NULL;

	    if (i == ACTUAL_META && strcmp(metadata_strings[0], "Near Real Time"))
		NRT_infix = "";
            if (i < ARCHIVE_META)
               handle = INVENTORY_METADATA;
            else
               handle = ARCHIVED_METADATA;

	    if(PGS_MET_SetAttr(mdHandles[handle],
		echoed_metadata[i].object, metadata_strings) != PGS_S_SUCCESS)
	    {
		sprintf(msgbuf, "PGS_MET_SetAttr(\"%s\")",
		    echoed_metadata[i].object);
		modsmf(MODIS_E_GEO, msgbuf, filefunc);
		retval = MODIS_E_GEO;
	    }
	}
    }	/* for each Inventory metadata */

    /* Release memory ocupied by PVL tree stucture.	*/
    if(getMODISECSinfo(l1a_file, NULL, NULL, NULL, NULL, NULL) != MAPIOK)
    {
	sprintf(msgbuf, "getMODISECSinfo(\"%s\",NULL)", l1a_file->filename);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
    }

    if ( retval == PGS_S_SUCCESS)
    {
	double equat_crossing, longitude[4], latitude[4], bound_coord[4];
	int  orb_num;
	char *p_flag, *time, *date;

	if( EPH_metadata != NULL)
	{
	    i = EQUAT_META;
	    equat_crossing =
		(double)(EPH_metadata->equatorcrossinglongitude * RAD2DEG);
	    geo_metadata[i++].data = &equat_crossing;
	    orb_num = (int)EPH_metadata->orbitnumber;
	    geo_metadata[i++].data = &orb_num;
	    date = EPH_metadata->equatorcrossingdate;
	    geo_metadata[i++].data = &date;
	    time = EPH_metadata->equatorcrossingtime;
	    geo_metadata[i++].data = &time;
	}else
	    modsmf(MODIS_W_MISSING_OUTPUT, "EPH_metadata", filefunc);

	if( GRing_points != NULL)
	{
	    int j;

	    i = GRING_META;
	    p_flag = GRing_points->exclusionflag;
	    geo_metadata[i++].data = &p_flag;
	    for (j = 0; j < 4; j++){
	       latitude[j]  = (double)GRing_points->latitude[j];
	       longitude[j] = (double)GRing_points->longitude[j];
	    }
	    geo_metadata[i++].data = latitude;
	    geo_metadata[i++].data = longitude;
	    geo_metadata[i++].data = GRing_points->sequenceno;

	}else
	    modsmf(MODIS_W_MISSING_OUTPUT, "GRing_points", filefunc);

        if ( bounding_coords != NULL )
        {
            i = BOUND_META;
            bound_coord[0] = (double)bounding_coords->eastcoord * RAD2DEG;
            geo_metadata[i++].data = &bound_coord[0];
            bound_coord[1] = (double)bounding_coords->westcoord * RAD2DEG;
            geo_metadata[i++].data = &bound_coord[1];
            bound_coord[2] = (double)bounding_coords->northcoord * RAD2DEG;
            geo_metadata[i++].data = &bound_coord[2];
            bound_coord[3] = (double)bounding_coords->southcoord * RAD2DEG;
            geo_metadata[i++].data = &bound_coord[3];
        }else
            modsmf(MODIS_W_MISSING_OUTPUT, "bounding_coords", filefunc);

	if(PGS_MET_GetSetAttr(mdHandles[INVENTORY_METADATA], SHORTNAME,
	    &pshort) != PGS_S_SUCCESS)
	{
	    sprintf(msgbuf, "PGS_MET_GetSetAttr(\"%s\",\"" SHORTNAME "\")",
		mdHandles[INVENTORY_METADATA]);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);
	    return MODIS_E_GEO;
	}

	if(PGS_MET_GetSetAttr(mdHandles[INVENTORY_METADATA], VERSIONID,
	    &versionid) != PGS_S_SUCCESS)
	{
	    sprintf(msgbuf, "PGS_MET_GetSetAttr(\"%s\",\"" VERSIONID "\")",
		mdHandles[INVENTORY_METADATA]);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);
	    return MODIS_E_GEO;
	}

	sprintf(doi, "10.5067/MODIS/%s%s.%03ld",
	    shortname, NRT_infix, (long)versionid);
	if(putMODISfileinfo(l1a_file, IDENTIFIERPRODDOI, TXT, strlen(doi), doi)
	    != MAPIOK)
	{
	    sprintf(msgbuf, "putMODISfileinfo(\"%s\",\"" IDENTIFIERPRODDOI
		"\", \"" TXT ", %d, \"%s\")",
		l1a_file->filename, (int)strlen(doi), doi);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);
	    return MODIS_E_GEO;
	}

	if(putMODISfileinfo(l1a_file, IDENTIFIERPRODDOIAUTH, TXT,
	    sizeof(DOIAUTHORITY)-1, DOIAUTHORITY)
	    != MAPIOK)
	{
	    sprintf(msgbuf, "putMODISfileinfo(\"%s\",\"" IDENTIFIERPRODDOIAUTH
		"\", \"" TXT ", %d, \"" DOIAUTHORITY "\")",
		l1a_file->filename, (int)sizeof(DOIAUTHORITY)-1);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);
	    return MODIS_E_GEO;
	}

	for(i=0; i<(int)(sizeof(geo_metadata)/sizeof(geo_metadata[0])); i++)
	{
            if(strncmp(geo_metadata[i].attribute, MECS_CORE, sizeof(MECS_CORE))
		== 0)
               handle = INVENTORY_METADATA;
            else
               handle = ARCHIVED_METADATA;

	    if(geo_metadata[i].data != NULL &&
	      PGS_MET_SetAttr(mdHandles[handle],
	      geo_metadata[i].object, geo_metadata[i].data)
	      != PGS_S_SUCCESS)
	    {
		sprintf(msgbuf, "PGS_MET_SetAttr(\"%s\")",
		    geo_metadata[i].object);
		modsmf(MODIS_E_GEO, msgbuf, filefunc);
		retval = MODIS_E_GEO;
		break;
	    }
	} /* for each geo_metadata */

	if(PGS_MET_Write(mdHandles[INVENTORY_METADATA], MECS_CORE,
	    (PGSt_integer)l1a_file->sd_id) !=  PGS_S_SUCCESS)
	{
	    modsmf(MODIS_E_GEO, "PGS_MET_Write(\"" MECS_CORE "\")", filefunc);
	    retval = MODIS_E_GEO;
	}

        if(PGS_MET_Write(mdHandles[ARCHIVED_METADATA], MECS_ARCHIVE,
            (PGSt_integer)l1a_file->sd_id) !=  PGS_S_SUCCESS)
        {
            modsmf(MODIS_E_GEO, "PGS_MET_Write(\"" MECS_ARCHIVE "\")", filefunc);
            retval = MODIS_E_GEO;
        }

    }

    PGS_MET_Remove();

    return retval;
}

