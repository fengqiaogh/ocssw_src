#include "PGS_TD.h"
#include "GEO_earth.h"
#include "GEO_output.h"
#include "GEO_product.h"
#include "PGS_MODIS_35251.h"
#include "L1a_data.h"

PGSt_SMF_status GEO_write_granule_metadata(
	MODFILE			* const geo_file,
	MODFILE			* const l1a_file,
	GEO_param_struct const	* const geo_parameter,
	GEO_bcoord_struct const	* const bounding_coords,
	int			const version,
	l1a_data_struct		* const l1a_data,
	qa_metadata_struct	* const qa_metadata
)

/**********************************************************************
!C 

!Description:   
	Routine in Output group of the Level-1A geolocation software to write
	the granule-level metadata to the geolocation product file.

!Input Parameters:
	geo_file	MAPI structure for geolocation file
	l1a_file	MAPI structure for Level 1A file
	geo_parameter	program parameters read from (or derived directly from)
			the geolocation parameter file
	bounding_coords	the granule's bounding coordinates.
	version		version number of input/output files.

!Output Parameters:
	None

!Input/Output Parameters
	l1a_data	data read from L1A file
	qa_metadata	Granule-level quality assurance metadata. Optional,
			can be null.

Return Values:
	MODIS_E_GEO		If any subroutine fails.
	MODIS_E_GEO_EPH		If the ephemeris file contains an orbit number
				    that is inconsistent with the orbit time.
	MODIS_E_BAD_INPUT_ARG	If any pointer argument other than qa_metadata
				    is null.
	PGS_S_SUCCESS		Otherwise.

Externally Defined:
	AVERAGE_TEMPERATURES		"GEO_product.h"
	FILL_INT8			"GEO_output.h"
	L1A_ENGINEERING_GRP		"L1a_data.h"
	MAPIOK				"mapi.h"
	MAX_SCAN_NUMBER			"GEO_geo.h"
	MODIS_E_BAD_INPUT_ARG		"PGS_MODIS_35251.h"
	MODIS_E_GEO			"PGS_MODIS_35251.h"
	MODIS_E_GEO_EPH                 "PGS_MODIS_35251.h"
	MODIS_W_GEO_EPHEMERIS		"PGS_MODIS_35251.h"
	MODIS_W_NO_GEO			"PGS_MODIS_35251.h"
	PGSd_EOS_AM			"PGS_TD.h"
	PGSd_EOS_PM			"PGS_TD.h"
	PGSd_PC_FILE_PATH_MAX           "PGS_PC.h"
	PGSd_PC_LINE_LENGTH_MAX		"PGS_PC.h"
	PGSEPH_E_NO_SC_EPHEM_FILE	"PGS_EPH.h"
	PGS_S_SUCCESS			"PGS_SMF.h"
	SUCCESS				"GEO_basic.h"
	TIMECODEASIZE			"smfio.h"

Called by:
	GEO_locate_one_granule


Routines Called:
	GEO_get_GRing_points		"GEO_earth.h"
	GEO_get_utcpole_metadata	"GEO_output.h"
	GEO_get_version_metadata	"GEO_output.h"
	GEO_update_L1A_metadata		"GEO_output.h"
	GEO_write_ECS_metadata		"GEO_output.h"
	GEO_write_geospecific_metadata	"GEO_output.h"
	GEO_write_input_metadata	"GEO_output.h"
	GEO_write_parameters		"GEO_output.h"
	modsmf				"smfio.h"
	PGS_EPH_GetEphMet		"PGS_EPH.h"
        PGS_TD_UTCtoTAI                 "PGS_TD_Prototypes.h"
	putMODIStable			"mapi.h"
	PGS_SMF_GenerateStatusReport    "PGS_SMF.h"


!Revision History:
$Log: GEO_write_granule_metadata.c,v $
Revision 6.3  2010/06/30 21:28:26  kuyper
Added casts of cumulated_gflags to a printable type.

Revision 6.2  2010/06/18 20:57:58  kuyper
Corrected format of LogReport message.

Revision 6.1  2010/05/28 21:50:04  kuyper
Helped resolve Bug 2470 by removing a parameter.
Changed to return an SMF status code.
Improved several error messages.
Corrected sprintf() format codes for C90 compatibility.
Improved const-safety.

James Kuyper Jr.		James.R.Kuyper@NASA.gov

Revision 5.6  2008/12/19 19:07:29  kuyper
Corrected size of orbit_validation array.

Revision 5.5  2008/12/15 21:12:59  kuyper
Changed to make orbit number validation optional, controlled by a PCF
  parameter.

Revision 5.4  2008/09/21 19:06:38  kuyper
Increased size of message buffer to accomodate long path names.
Changed to pass sc_tag to GEO_write_input_metadata().

Revision 5.3  2005/10/05 16:51:09  vlin
Code added to check for invalid Orbit Number from Ephemeris file.
MOD_PR03 will fail if the orbit number is incorrect.
vlin@saicmodis.com

Revision 5.2  2004/10/15 22:03:17  kuyper
Corrected position of bounding_coords parameter.

Revision 5.1  2004/09/22 19:44:25  vlin
Changed to pass bounding_coords to GEO_update_L1A_metadata().
Added maneuver_list parameter, passed on to GEO_write_geospecific_metadata.
Added a status report.    vlin@saicmodis.com

Revision 4.5  2004/04/09 22:17:45  kuyper
Changed FILL_BYTE to FILL_INT8.

Revision 4.4  2004/01/23 21:40:35  vlin
PGS_PC_GetUniversalReference not called and removed from prologue.

Revision 4.3  2004/01/20 19:44:27  kuyper
Corrected error messages.
Added validity checks on file handles.

Revision 4.2  2003/10/07 13:53:29  kuyper
Corrected setting of equatorcrossingtime.

Revision 4.1  2003/08/12 23:05:33  kuyper
Corrected to identify l1a_data and qa_metadata as input/output.
Added code to write the AverageTemperatures vdata record.
Moved input pointer handling into GEO_write_input_metadata().
Dropped code setting modisproductfilename.
Simplified handling of orbitdescendtime.
Changed to expect status codes from GEO_update_l1a_metadata() and
  GEO_write_ECS_metadata.
Changed format for RMS error.

Requirements:
      PR03-F-1-3
      PR03-F-1-3.1
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

!END
*********************************************************************************/
{
    static pointer_metadata_struct pointer_metadata;

    EPH_metadata_struct EPH_metadata={-1, -999.0*DEG2RAD, "", ""};
    GEO_GRing_struct GRing_points = {"N", {1, 2, 3, 4},
      {0.01, 0.01, 0.0,  0.0},	/* Latitude */
      {0.0,  0.01, 0.01, 0.0}	/* Longitude. */};

    utcpole_metadata_struct utcpole_metadata = {"", {0}};
    PGSt_double offset = 0.0;	/* offset from rangedatetime, seconds.*/
    PGSt_integer numOrbits=1;
    PGSt_tag sc_tag=PGSd_EOS_AM;
    int i;
    PGSt_SMF_status	retval = PGS_S_SUCCESS;
    double pred_dtime; 
    char msgbuf[2*PGSd_PC_FILE_PATH_MAX];
    char orbitAscendTime[TIMECODEASIZE];
    char orbitdescendtime[TIMECODEASIZE];
    char rangedatetime[TIMECODEASIZE];
    char sci_abnorm[2]="1", sci_state[2]="1";
    PGSt_integer param_version=1;
    static char filefunc[] = __FILE__ ", GEO_write_granule_metadata";
    int rtcode;

    if (l1a_file==NULL || geo_file==NULL || l1a_data == NULL ||
	geo_parameter == NULL || bounding_coords == NULL)
    {
	sprintf(msgbuf, "l1a_file = %p, geo_file=%p, l1a_data = %p,\n"
            "geo_parameter = %p, bounding_coords=%p",
	    (void*)l1a_file, (void*)geo_file, (void*)l1a_data,
	    (void*)geo_parameter, (void*)bounding_coords);
	modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);

	return MODIS_E_BAD_INPUT_ARG;
    }

    /* Read l1a_data.ECS_metadata and fill in uptcpole_metadata */
    if (GEO_get_utcpole_metadata(&(l1a_data->ECS_metadata), &utcpole_metadata)
	!= PGS_S_SUCCESS) {
	modsmf(MODIS_E_GEO, "GEO_get_utcpole_metadata()", filefunc);

	retval = MODIS_E_GEO;
    }

    /* Write the geolocation geo_parameter's to the geolocation file*/
    if( SUCCESS != GEO_write_parameters(geo_file,
	(GEO_param_struct *)geo_parameter) )
    {
      sprintf(msgbuf, "GEO_write_parameters(\"%s\")", geo_file->filename);
      modsmf(MODIS_E_GEO, msgbuf, filefunc);

      retval = MODIS_E_GEO;
    }

    /* Calculate spatial metadata */
    rtcode = GEO_get_GRing_points(&GRing_points);

    /* Since the V2.1 GEO_get_GRing_points will return WARNING for several
     * 'sparse data' conditions (which are more likely to occur than the
     * functional failures) retval is set to FAIL only for a FAIL return.
     */
    if(rtcode != PGS_S_SUCCESS && rtcode != MODIS_W_NO_GEO)
    {
	modsmf(MODIS_E_GEO, "GEO_get_GRing_points()", filefunc);

	retval = MODIS_E_GEO;
    }

    sprintf(rangedatetime, "%.10sT%.15sZ",
	l1a_data->ECS_metadata.rangebeginningdate,
	l1a_data->ECS_metadata.rangebeginningtime);
    
    if(strncmp(geo_parameter->spacecraft_ID, "Terra", sizeof("Terra")) )
    {
	sc_tag = PGSd_EOS_PM;	/* Already validated; must be "Aqua" */
	param_version = 2;
    }

    rtcode = PGS_EPH_GetEphMet(sc_tag, 1, rangedatetime, &offset, &numOrbits,
	&EPH_metadata.orbitnumber, &orbitAscendTime, &orbitdescendtime,
	&EPH_metadata.equatorcrossinglongitude);

    if(rtcode != PGS_S_SUCCESS)
    {
	if(rtcode == PGSEPH_E_NO_SC_EPHEM_FILE)
	    modsmf(MODIS_W_GEO_EPHEMERIS, rangedatetime, filefunc);
	else
	{
	    sprintf(msgbuf, "PGS_EPH_GetEphMet(\"%s\")", rangedatetime);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);

	    retval = MODIS_E_GEO;
	}
    }
    else
    {
      static char orbit_validation[PGSd_PC_LINE_LENGTH_MAX];
      PGSt_double TAI_orbitdescendtime;

      if(orbit_validation[0] == '\0' && 
	(PGS_PC_GetConfigData(VALIDATE_ORBIT_NO_LUN, orbit_validation)
	!= PGS_S_SUCCESS ||
	(strcmp(orbit_validation, "TRUE") && 
	strcmp(orbit_validation, "FALSE"))) )
      {
	sprintf(msgbuf, "PGS_PC_GetConfigData():\"%s\"", orbit_validation);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
	orbit_validation[0] = '\0';
      }


	/* Assumes that orbitdescendtime has the format of a UTC A time code:
	 * "YYYY-MM-DDTHH:MM:SS.SSSSSSZ"
	 */
      if (PGS_TD_UTCtoTAI(orbitdescendtime, &TAI_orbitdescendtime) !=
	 PGS_S_SUCCESS)
      {
	 sprintf(msgbuf, "PGS_TD_UTCtoTAI(\"%s\")", orbitdescendtime);
	 modsmf(MODIS_E_GEO, msgbuf, filefunc);

	 retval = MODIS_E_GEO;
      }
      else if(strcmp(orbit_validation, "TRUE")==0)
      {
	int step = (EPH_metadata.orbitnumber < 
	       geo_parameter->orbit_valid_params.transition_orbit) ? 0 : 1;
	pred_dtime = geo_parameter->orbit_valid_params.descend_time_0[step] +
	  geo_parameter->orbit_valid_params.period[step]*
	  EPH_metadata.orbitnumber;
	if (abs(TAI_orbitdescendtime - pred_dtime) > 
	       geo_parameter->orbit_valid_params.orbit_tolerance[step])
	{
	  sprintf(msgbuf,". Orbit number = %d, orbit descend time = %s",
		  EPH_metadata.orbitnumber, orbitdescendtime);
	  modsmf(MODIS_E_GEO_EPH, msgbuf, filefunc);

	  retval = MODIS_E_GEO_EPH;
	}
      }

      memcpy(EPH_metadata.equatorcrossingdate, orbitdescendtime, 10);
      EPH_metadata.equatorcrossingdate[10] = '\0';
      memcpy(EPH_metadata.equatorcrossingtime, orbitdescendtime+11, 15);
      EPH_metadata.equatorcrossingtime[15] = '\0';
    }

    if( GEO_update_L1A_metadata(l1a_file, &GRing_points, &EPH_metadata,
	(GEO_bcoord_struct *)bounding_coords) != PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "GEO_update_L1A_metadata(%s)", l1a_file->filename);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);

	retval = MODIS_E_GEO;
    }  

    if(GEO_write_input_metadata(geo_file, version, param_version,
	&pointer_metadata, sc_tag) != PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "GEO_write_input_metadata(\"%s\", %ld, %ld )",
	    geo_file->filename, (long)version, (long)param_version);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);

	retval = MODIS_E_GEO;
    }

    /* Granule-level SCI_STATE and SCI_ABNORM */
    if (l1a_data->num_scans <= MAX_SCAN_NUMBER) /* array bounds protection */
	for(i=0; i<l1a_data->num_scans; i++)
	{
	    if((l1a_data->frame_data[i].SCI_STATE!=1)
		&& (l1a_data->frame_data[i].SCI_STATE != FILL_INT8) )
		sci_state[0] = '0';
	
	    if((l1a_data->frame_data[i].SCI_ABNORM!=1)
		&& (l1a_data->frame_data[i].SCI_ABNORM!= FILL_INT8) )
		sci_abnorm[0] = '0';
	}
    else
	sci_state[0] = sci_abnorm[0] = '0';

    if (qa_metadata)
    {
	sprintf(qa_metadata->rms_error,"%-8.3g", geo_parameter->RMS_error);

	if (retval != SUCCESS)
	    qa_metadata->retval = retval;
    }

    if (GEO_get_version_metadata(&(l1a_data->ECS_metadata.version_metadata))
	!= PGS_S_SUCCESS)
    {
	modsmf(MODIS_E_GEO, "GEO_get_version_metadata()", filefunc);

	retval = MODIS_E_GEO;
    }

    if(GEO_write_ECS_metadata(geo_file, &l1a_data->ECS_metadata, &EPH_metadata,
	(GEO_bcoord_struct *)bounding_coords, &GRing_points,
	&pointer_metadata, qa_metadata, sci_state, sci_abnorm) != PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "GEO_write_ECS_metadata(\"%s\")", geo_file->filename);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);

	retval = MODIS_E_GEO;
    }

    /* Write geolocation specific metadata */
    if(GEO_write_geospecific_metadata( geo_file, &l1a_data->granule_metadata,
	l1a_data->num_scans, geo_parameter, qa_metadata, &utcpole_metadata)
	!= PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "GEO_write_geospecific_metadata(\"%s\")",
	    geo_file->filename);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);

	retval = MODIS_E_GEO;
    }

    if(putMODIStable(geo_file, AVERAGE_TEMPERATURES, L1A_ENGINEERING_GRP,
	0, 1, (unsigned char*)l1a_data->temperatures) != MAPIOK)
    {
	modsmf(MODIS_E_GEO, "putMODIStable(\"" AVERAGE_TEMPERATURES "\")",
	    filefunc);

	retval = MODIS_E_GEO;
    }

    /* Write to status report */
    sprintf(msgbuf,
	"Input file name: %s\n"
	"Output file name: %s\n"
	"Local Granule ID: %s\n"
	"scans: %d\n"
	"pixels: %d\n"
	"\tmissing: %d\n"
	"\tout of bounds: %d\n"
	"Cumulated gflags: %lu, %lu, %lu, %lu, %lu, %lu, %lu, %lu\n"
	"Bounding coordinates\n"
	"\teast: %f\n"
	"\twest: %f\n"
	"\tnorth: %f\n"
	"\tsouth: %f",
	l1a_file->filename, geo_file->filename,
	l1a_data->ECS_metadata.localinputgranuleid, l1a_data->num_scans,
	qa_metadata->no_of_pixels, qa_metadata->missingdata,
	qa_metadata->outofboundsdata,
	(unsigned long)qa_metadata->cumulated_gflags[0],
	(unsigned long)qa_metadata->cumulated_gflags[1],
	(unsigned long)qa_metadata->cumulated_gflags[2],
	(unsigned long)qa_metadata->cumulated_gflags[3],
	(unsigned long)qa_metadata->cumulated_gflags[4],
	(unsigned long)qa_metadata->cumulated_gflags[5],
	(unsigned long)qa_metadata->cumulated_gflags[6],
	(unsigned long)qa_metadata->cumulated_gflags[7],
	bounding_coords->eastcoord,
	bounding_coords->westcoord, bounding_coords->northcoord,
	bounding_coords->southcoord);
    PGS_SMF_GenerateStatusReport(msgbuf);

    return retval;
}
