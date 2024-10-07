#include "GEO_output.h"
#include "GEO_product.h"
#include "PGS_MODIS_35251.h"
#include "granule_extract_metadata.h"

PGSt_SMF_status  GEO_write_geospecific_metadata(
	MODFILE				* const geo_file, 
        l1a_metadata_struct const	* const granule_metadata,
	int				const number_of_scans,
        GEO_param_struct const		* const geo_parameter,
	qa_metadata_struct const	* const qa_metadata,
	utcpole_metadata_struct const	* const utcpole_metadata
)
/*!C****************************************************************************
!Description:   
	Routine in Output group of the Level-1A geolocation software to write
	geolocation specific granule metadata to the output product.  These
	metadata are stored as individual global HDF attributes using MAPI
	calls.

	Warning messages will be recorded (to the SMF) if metadata inputs are
	missing, but the routine will not abort, which may result in
	geolocation products with incomplete granule metadata.

!Input Parameters:
	geo_file                the MAPI structure for the product
	granule_metadata        l1a specific granule metadata
	number_of_scans         the number of scans in the granule
	geo_parameter           program parameters read from (or derived
	directly from) the geolocation parameter file
	qa_metadata             Quality metadata
	utcpole_metadata        metadata identifying the utcpole.dat file

!Output Parameters:
	None

Return Values:
	MODIS_E_GEO		If putMODISfileinfo() fails.
	MODIS_E_BAD_INPUT_ARG	If geo_file is null.
	MODIS_E_MISSING_OUTPUT	If an optional data source was unavailable.
	PGS_S_SUCCESS		Otherwise

Externally Defined:
	BAD_PCKTS			"GEO_product.h"
	CUM_GFLAGS			"GEO_product.h"
	DATATYPELENMAX			"mapi.h"
	DISCRD_PKTS			"GEO_product.h"
	EA_SOURCE			"GEO_product.h"
	EA_SOURCE_SELECT_LUN		"GEO_geo.h"
	FROFF_FRAME_2                   "GEO_product.h"
	FROFF_FRAME_4                   "GEO_product.h"
	FROFF_SCAN_20                   "GEO_product.h"
	FROFF_SCAN_40                   "GEO_product.h"
	GEO_EST_RMS_ERROR               "GEO_product.h"
	PARVERS				"GEO_product.h"
	INCOMP_SCANS                    "GEO_product.h"
	I32				"mapi.h"
	MAPIOK				"mapi.h"
	MAX_EFRM	                "GEO_product.h"
	MAX_SFRM                        "GEO_product.h"
	MAX_SV_FRM                      "GEO_product.h"
	MISS_PCKTS			"GEO_product.h"
	MODIS_E_BAD_INPUT_ARG		"PGS_MODIS_35251.h"
	MODIS_E_GEO			"PGS_MODIS_35251.h"
	MODIS_E_MISSING_OUTPUT		"PGS_MODIS_35251.h"
	NUMSCN                          "GEO_product.h"
	PGSd_PC_VALUE_LENGTH_MAX	"PGS_PC.h"
	POLAR_MOTION                    "GEO_product.h"
	R32				"mapi.h"
	TERRAIN_CORRECT_LUN		"GEO_geo.h"
	TERRAIN_CORRECTION		"GEO_product.h"
	TXT				"mapi.h"
	UI32				"mapi.h"
	UTCPOLE_FH                      "GEO_product.h"

Called by:
	GEO_write_granule_metadata

Routines Called:
	putMODISfileinfo                "mapi.h"
	modsmf                          "smfio.h"
	PGS_PC_GetConfigData            "PGS_PC.h"


!Revision History:
$Log: GEO_write_geospecific_metadata.c,v $
Revision 6.4  2011/02/14 21:59:46  kuyper
Removed cast that is no longer needed with current version of M-API toolkit.

Revision 6.3  2010/05/28 21:34:38  kuyper
Resolved Bug 2249 by writing value of terrain correction flag as a file
  attribute.
Helped resolve Bug 2470 by dropping a parameter.
Changed to return status code.
Improved const-safety.

Revision 6.3  2010/05/28 21:20:17  kuyper
Resolved Bug 2249 by writing value of terrain correction flag as a file
  attribute.
Helped resolve Bug 2470 by dropping a parameter.
Improved const-safety.

Revision 6.2  2009/05/30 22:53:03  kuyper
Corrected metadata name, case number, and message contents.

Revision 6.1  2009/05/28 22:59:01  kuyper
Added FractionalOffset file attributes.
Corrected prolog.
Improved error messages.

James Kuyper Jr.	James.R.Kuyper@nasa.gov

Revision 5.1  2004/09/13 18:06:50  vlin
1. input parameter maneuver_list added.
2. maneuver_list->revision written as a file attribute to output file.

Revision 4.2  2003/10/27 18:49:23  vlin
Size of message buffer was changed.

Revision 4.1  2003/02/21 22:54:51  kuyper
Corrected to use void* pointers with %p format code.

Requirements:
		PR03-F-4.1-1
		PR03-F-4.1-1.6
		PR03-F-4.1-2
		PR03-I-1
		PR03-I-2
		PR03-S-1

!Team-unique Header:

    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

References and Credits
	None

Design Notes

!END****************************************************************************
*/

{
    static float32 zero = 0.0, half = 0.5;
    struct {
       char		*name;                      /* the attribute's name. */
       char		data_type[DATATYPELENMAX];  /* M-API data type. */
       long		n_elements;			/* number of elements */
       void const	*data;                      /* pointer to data */
    }attribute_object[] = {
	{ NUMSCN,		I32,	1, NULL},
	{ MAX_EFRM,		I32,	1, NULL},
	{ MAX_SFRM,	        I32,	1, NULL},
	{ MAX_SV_FRM,		I32,	1, NULL},
	{ INCOMP_SCANS,		I32,	1, NULL},
	{ MISS_PCKTS,		I32,	1, NULL},
	{ BAD_PCKTS,		I32,	1, NULL},
	{ DISCRD_PCKTS,		I32,	1, NULL},
	{ EA_SOURCE,		TXT,	0, NULL},
	{ TERRAIN_CORRECTION,	TXT,	0, NULL},
	{ PARVERS,		TXT,	0, NULL},
	{ GEO_EST_RMS_ERROR,	R32,	1, NULL},
	{ CUM_GFLAGS,		UI32,	8, NULL},
	{ UTCPOLE_FH,		TXT,	0, NULL},
	{ POLAR_MOTION,		R64,	3, NULL},
	{ NULL,			R32,	1, &half},
	{ NULL,			R32,	1, &zero}
    };

    static char ea_source[PGSd_PC_VALUE_LENGTH_MAX] = "";
    static char terrain_correction[PGSd_PC_VALUE_LENGTH_MAX] = "";
    int32 tmp_int32;
    int num_objects;
    int obj = 0;
    char filefunc[] = __FILE__ ", GEO_write_geospecific_metadata";
    char msgbuf[PGS_SMF_MAX_MSGBUF_SIZE] = "";
    PGSt_SMF_status ret_val = PGS_S_SUCCESS;

    if (geo_file  == NULL )
    {
	sprintf(msgbuf, "geo_file = %p ", (void*)geo_file);
	modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
	return MODIS_E_BAD_INPUT_ARG;
    }

    /*Identify the source of spacecraft kinematic data */
    if( ea_source[0] == '\0' && 
	(PGS_PC_GetConfigData(EA_SOURCE_SELECT_LUN, ea_source)!= PGS_S_SUCCESS))
    {
       sprintf(msgbuf, "PGS_PC_GetConfigData(%d)", EA_SOURCE_SELECT_LUN);
       modsmf(MODIS_E_GEO, msgbuf, filefunc);
       ea_source[0] = '\0';
    }

    /*Identify the source of spacecraft kinematic data */
    if( terrain_correction[0] == '\0' && 
	(PGS_PC_GetConfigData(TERRAIN_CORRECT_LUN, terrain_correction)
	!= PGS_S_SUCCESS))
    {
       sprintf(msgbuf, "PGS_PC_GetConfigData(%d)", TERRAIN_CORRECT_LUN);
       modsmf(MODIS_E_GEO, msgbuf, filefunc);
       terrain_correction[0] = '\0';
    }

    /*
    Finish initializing the attribute_object array with the contents 
	  defined above.  IF any of the input arguments are NULL or 
	  ea_source was not successfully retrieved, set
	  the dependent .data('s) to NULL.
    */

    /*  int32 */
    tmp_int32 = (int32)number_of_scans;
    attribute_object[obj].data = &tmp_int32;
    obj++;

    if (granule_metadata != NULL){

       /* element 1 int32 */
       attribute_object[obj].data = &granule_metadata->max_earth_frames;
       obj++;

       /* element 2 int32 */
       attribute_object[obj].data = &granule_metadata->max_sd_frames;
       obj++;

       /*element 3  int32 */
       attribute_object[obj].data = &granule_metadata->max_sv_frames;
       obj++;

       /*element 4  int32 */
       attribute_object[obj].data = &granule_metadata->incomplete_scans;
       obj++;

       /*element 5  int32 */
       attribute_object[obj].data = &granule_metadata->missing_packets;
       obj++;

       /*element 6  int32 */
       attribute_object[obj].data = &granule_metadata->bad_CRC_packets;
       obj++;

       /* element 7  int32 */
       attribute_object[obj].data = &granule_metadata->discarded_packets;
       obj++;
    }else
       obj += 7;

    /* TXT */
    attribute_object[obj].n_elements = (long)strlen(ea_source);
    attribute_object[obj].data =
	attribute_object[obj].n_elements ? ea_source : NULL;
    obj++;

    /* TXT */
    attribute_object[obj].n_elements = (long)strlen(terrain_correction);
    attribute_object[obj].data =
	attribute_object[obj].n_elements ? terrain_correction : NULL;
    obj++;

    if ( geo_parameter != NULL){
       /* TXT */
       attribute_object[obj].n_elements = (long)strlen(geo_parameter->revision);
       attribute_object[obj].data =  geo_parameter->revision;
       obj++;

       /* float32 */
       attribute_object[obj].data = &geo_parameter->RMS_error;
       obj++;
    }else
       obj +=2;

    if (qa_metadata != NULL) 
       /* uint32 */
       attribute_object[obj].data = qa_metadata->cumulated_gflags;
    obj++;

    if (utcpole_metadata != NULL)
    {
	/* TXT */
	attribute_object[obj].n_elements =
	   (long)strlen(utcpole_metadata->header);
	attribute_object[obj++].data = utcpole_metadata->header;

	/* R32 */
	attribute_object[obj++].data = utcpole_metadata->polar_motion;
    } else
	obj += 2;

    switch(geo_parameter->geometry_params.N_samp
	[geo_parameter->geometry_params.band_number])
    {
	default:
	    break;
	case 2:
	    attribute_object[obj++].name = FROFF_SCAN_20;
	    attribute_object[obj++].name = FROFF_FRAME_2;
	    break;
	case 4:
	    attribute_object[obj++].name = FROFF_SCAN_40;
	    attribute_object[obj++].name = FROFF_FRAME_4;
	    break;
    }
    num_objects = obj;

    for (obj =0; obj < num_objects; obj++)
    {
	if (attribute_object[obj].data == NULL)
	{
	    /* "Output value not available for" */
	    modsmf(MODIS_E_MISSING_OUTPUT, attribute_object[obj].name,
		filefunc);

	    ret_val = MODIS_E_MISSING_OUTPUT;
	}else if(putMODISfileinfo(geo_file, attribute_object[obj].name,
	    attribute_object[obj].data_type, attribute_object[obj].n_elements,
	    attribute_object[obj].data) != MAPIOK)
	{
	    sprintf(msgbuf, "putMODISfileinfo(\"%s\", \"%s\", \"%s\", %ld)",
		geo_file->filename, attribute_object[obj].name,
		attribute_object[obj].data_type,
		attribute_object[obj].n_elements);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);

	    ret_val = MODIS_E_GEO;
	}
    }

    /* Write Extract Metadata */
    write_extract_metadata(geo_file->sd_id,
                           granule_metadata->Extract_Pixel_Offset,
                           granule_metadata->Extract_Pixel_Count,
                           granule_metadata->Extract_Line_Offset,
                           granule_metadata->Extract_Line_Count);

    return ret_val;
}
