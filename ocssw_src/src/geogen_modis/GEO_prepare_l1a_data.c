#include "smfio.h"
#include "GEO_global_arrays.h"
#include "GEO_input.h"
#include "PGS_MODIS_35251.h"

int num_impulse[MAX_SCAN_NUMBER];

PGSt_SMF_status GEO_prepare_l1a_data(
	MODFILE			* const l1a_file,
	GEO_param_struct const	* const geo_params,
	l1a_data_struct		* const l1a_data
)
/*
!C*****************************************************************************
!Description:   
		subroutine in input group of the Level-1A geolocation
                software to read and process geolocation input data from
		the Level 1A file.  First it reads relevant Level 1A 
                metadata which are required for geolocation.  It then reads 
                scan start times, mirror sides and number of frames for the 
                scans, followed by the the mirror encoder times, sector 
                start times and spacecraft ancillary data as packed arrays. 
                The packed data are unpacked, converted and validated for 
                use in the geolocation processing.

!Input Parameters:
                l1a_file - the MAPI structure for the Level 1A file 
                geo_params - structure containing geolocation process
                        parameters from the parameter file.

!Output Parameters:
		l1a_data - structure for data read from L1A file

Return parameter:
		MODIS_E_BAD_INPUT_ARG		If any pointer argument is null
		MODIS_E_GEO			If any subroutine fails.
		MODIS_E_GEO_SCANNO_INPUT	Invalid number of scans in file
		MODIS_E_GEO_WRONG_PLATFORM	L1A file is for wrong platform
                PGS_S_SUCCESS			Otherwise

Externally Defined:
                CHAN_A                          "GEO_parameters.h"
                CHAN_B                          "GEO_parameters.h"
                ELEC_SIDES                      "GEO_parameters.h" 
                ENCODER_LENGTH                  "GEO_geo.h" 
                MAX_SCAN_NUMBER                 "GEO_geo.h" 
                MODIS_E_BAD_INPUT_ARG           "PGS_MODIS_325251.h" 
                MODIS_E_GEO                     "PGS_MODIS_325251.h"
                MODIS_E_GEO_FORMATTER           "PGS_MODIS_325251.h" 
                MODIS_E_GEO_SCANNO_INPUT        "PGS_MODIS_325251.h"
                MODIS_E_GEO_WRONG_PLATFORM      "PGS_MODIS_325251.h" 
                PGS_S_SUCCESS                   "PGS_SMF.h"
                SECTOR_LENGTH                   "GEO_geo.h" 
                SUCCESS                         "GEO_basic.h" 

Called By:
                GEO_locate_one_granule


Routines called:
                GEO_prepare_mirr_data           "GEO_input.h"
                GEO_prepare_ancil_data          "GEO_input.h"
                GEO_read_L1AECS_metadata        "GEO_input.h"
                GEO_read_L1Apacket_data         "GEO_input.h"
                GEO_read_L1Ascan_metadata       "GEO_input.h"
                GEO_read_L1Aspecific_metadata   "GEO_input.h"
                GEO_read_L1Atemp_data           "GEO_input.h"
                modsmf                          "smfio.h"

!Revision History:
$Log: GEO_prepare_l1a_data.c,v $
Revision 6.2  2010/06/29 20:18:56  kuyper
Corrected test for error return from GEO_prepare_ancil_data().

Revision 6.1  2010/06/01 23:03:49  kuyper
Helped resolve Bug 2472 by defining ss_cp_mode, passing it to
  GEO_read_L1Apacket_data() to be filled in, and to GEO_prepare_ancil_data() to
  be used.
Corrected to test for != PGS_SUCCESS, rather than == FAIL.
Changed to return status codes.
Changed to expect status codes to be returned by GEO_read_L1Apacket_data(),
  GEO_prepare_ancil_data(), and GEO_prepare_mirr_data().
Made error messages more informative.

James Kuyper Jr.		James.R.Kuyper@NASA.gov

Revision 4.3  2003/08/12 14:16:38  kuyper
Corrected to move mirr_side from globals to locals before calling
   GEO_prepare_mirr_data().

Revision 4.2  2003/06/03 13:57:56  vlin
Moved GEO_prepare_ancil_data() ahead of GEO_prepare_mirr_data()
call GEO_read_L1Atemp_data() to fill in l1a_data.temperatures array.
vlin@saicmodis.com

Revision 4.1  2003/03/12 17:42:05  kuyper
Corrected pointer types passed to %p.

Revision 3.4  2001/06/05 23:54:48  kuyper
Changed GEO_prepare_ancil_data failures back to being non-fatal, per PDL.

10/10/95
Tracey W. Holmes
Added debug option. 

Requirements:	
                PR03-F-2.1-1
                PR03-F-2.1-2
                PR03-F-2.1-3
                PR03-F-2.2-1
                PR03-F-2.3-1
                PR03-F-2.4-1
                PR03-F-2.4-2
                PR03-F-2.5-1
                PR03-F-2.5-2
                PR03-F-2.5-3
                PR03-F-2.5-5
                PR03-I-1

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.
                
!END
*****************************************************************************/

{
    int    scan;  			    

    /* array of View Sector Start encoder count and Vernier Count words */
    int16  view_sector_start[MAX_SCAN_NUMBER][SECTOR_LENGTH] = {0};

    /* array of Earth Encoder Time Data (the time between every 100th mirror
     *encoder pulse over the Earth View)
     */
    uint16 earth_encoder_times[MAX_SCAN_NUMBER][ENCODER_LENGTH] = {0};

    /* data identifying the selected mirror assembly channel */
    int16  mirr_side[MAX_SCAN_NUMBER];
    uint16 FRside[MAX_SCAN_NUMBER][ELEC_SIDES];
    uint16 SAside[MAX_SCAN_NUMBER][ELEC_SIDES];
    uint16 ss_cp_mode[MAX_SCAN_NUMBER];
    sc_ancil_struct sc_ancillary_data[2];
    char msgbuf[512];    	
    char filefunc[] = __FILE__ ", GEO_prepare_l1a_data";
    
    /* Begin program logic */

    if ((l1a_file == NULL) || (geo_params == NULL) || (l1a_data == NULL)) {
	sprintf(msgbuf, " l1a_file = %p, geo_params = %p, l1a_data = %p",
	  (void*)l1a_file, (void*)geo_params, (void*)l1a_data);
	modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);

	return MODIS_E_BAD_INPUT_ARG;
    } 

    if (GEO_read_L1Aspecific_metadata(l1a_file, &l1a_data->granule_metadata, 
	&l1a_data->num_scans)!= SUCCESS)
    {
	sprintf(msgbuf, "GEO_read_L1Aspecific_metadata(\"%s\")",
	    l1a_file->filename);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);

	return MODIS_E_GEO;
    } 

    if (GEO_read_L1AECS_metadata(l1a_file, &l1a_data->ECS_metadata) != SUCCESS)
    {
	sprintf(msgbuf, "GEO_read_L1AECS_metadata(\"%s\")",
	    l1a_file->filename);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);

	return MODIS_E_GEO;
    }

    if (strncmp((l1a_data->ECS_metadata).platformshortname,
	geo_params->spacecraft_ID,sizeof(geo_params->spacecraft_ID))!=0)
    {
	sprintf(msgbuf, "platformshortname:%s\tdoesn't match spacecraft_ID: %s",
	    (l1a_data->ECS_metadata).platformshortname,
	    geo_params->spacecraft_ID);
	modsmf(MODIS_E_GEO_WRONG_PLATFORM, msgbuf, filefunc);

	return MODIS_E_GEO_WRONG_PLATFORM;
    }

    if ((l1a_data->num_scans < 0) || (l1a_data->num_scans > MAX_SCAN_NUMBER))
    {
	sprintf(msgbuf, "%d",l1a_data->num_scans);
	modsmf(MODIS_E_GEO_SCANNO_INPUT, msgbuf, filefunc);

	return MODIS_E_GEO_SCANNO_INPUT;
    } 
    else if (l1a_data->num_scans == 0) 
	return SUCCESS;

    if (GEO_read_L1Ascan_metadata(l1a_file, l1a_data->num_scans, 
	&l1a_data->granule_metadata, l1a_data->frame_data, mirr_side)
	!= SUCCESS)
    {
	sprintf(msgbuf, "GEO_read_L1Ascan_metadata(\"%s\", %ld)",
	    l1a_file->filename, (long)l1a_data->num_scans);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);

	return MODIS_E_GEO;
    }

    if (GEO_read_L1Apacket_data(l1a_file, l1a_data->num_scans, 
	earth_encoder_times, sc_ancillary_data, view_sector_start, FRside,
	SAside, ss_cp_mode) != PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "GEO_read_L1Apacket_data(\"%s\", %ld)",
	    l1a_file->filename, (long)l1a_data->num_scans);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);

	return MODIS_E_GEO;
    }

    if (GEO_read_L1Atemp_data(geo_params, l1a_file, l1a_data->num_scans, 
	l1a_data->temperatures) != PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "GEO_read_L1Atemp_data(\"%s\", %ld)",
	    l1a_file->filename, (long)l1a_data->num_scans);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);

	return MODIS_E_GEO;
    }

    for (scan = 0; scan < l1a_data->num_scans; scan++) {
      if (SAside[scan][CHAN_A]) {
	if (SAside[scan][CHAN_B]) {
	  if (FRside[scan][CHAN_A]) {
	      l1a_data->mirr_data[scan].mirr_chan = CHAN_A;
	      if (FRside[scan][CHAN_B]) {
		  sprintf(msgbuf, "%d", scan);
		  modsmf(MODIS_E_GEO_FORMATTER, msgbuf, filefunc);
	      }
	  } else l1a_data->mirr_data[scan].mirr_chan = CHAN_B;
	} else l1a_data->mirr_data[scan].mirr_chan = CHAN_A;
      } else l1a_data->mirr_data[scan].mirr_chan = CHAN_B;
    }
	    
  /* Note: the ancillary data must be initialized before the mirror data can
   * be prepared, because part of the preparation depends upon correcting
   * encoder times that happen at a fixed time after the reciept of ancillary
   * data packets.
   */ 

    if(GEO_prepare_ancil_data(l1a_data->num_scans, geo_params,
	sc_ancillary_data, ss_cp_mode) != PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "GEO_prepare_ancil_data(%ld)",
	    (long)l1a_data->num_scans);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
    }

    /* Needed by GEO_prepare_mirr_data. */
    for (scan = 0; scan < l1a_data->num_scans; scan++) 
	l1a_data->mirr_data[scan].mirr_side = (uint16)mirr_side[scan];

    if (GEO_prepare_mirr_data(earth_encoder_times, view_sector_start,
	&geo_params->mirror_prep_params, geo_params->geometry_params.t_frame,
	l1a_data) != SUCCESS)
    {
	sprintf(msgbuf, "GEO_prepare_mirr_data(%g)",
	    geo_params->geometry_params.t_frame);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);

	return MODIS_E_GEO;
    }

    /* Set by GEO_prepare_mirr_data. */
    for (scan = 0; scan < l1a_data->num_scans; scan++)
	num_impulse[scan] = (int)l1a_data->mirr_data[scan].num_impulse;

    return PGS_S_SUCCESS;
}

