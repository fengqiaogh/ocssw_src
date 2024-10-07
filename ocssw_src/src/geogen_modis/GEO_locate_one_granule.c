#include <time.h>
#include "hdf.h"
#include "HdfEosDef.h"
#include "mapi.h"
#include "GEO_earth.h"
#include "GEO_input.h"
#include "GEO_main_func.h"
#include "GEO_output.h"
#include "PGS_MODIS_35251.h"

PGSt_SMF_status GEO_locate_one_granule(
	char			l1a_file_name[PGSd_PC_FILE_PATH_MAX],
	char			geo_file_name[PGSd_PC_FILE_PATH_MAX],
	GEO_param_struct const	* geo_params,
	int			const version
	)
/*
!C*****************************************************************************
!Description:   
		subroutine in input group of the Level-1A geolocation
                software to geolocate one granule.  It opens the Level 1A 
		file and creates the geocation output file using the
		MAPI routine.  It calls routines to read and prepare the 
		Level 1A data for geolocation, geolocate each scan and
		write to the file, write the metadata and close the files.

!Input Parameters:
                l1a_file_name	path/file name string of the l1a input file
				for the granule
                geo_file_name	path/file name string of the geolocation output
                        	file for the granule
                geo_params	geolocation program parameters
                version		version of L1A file accessed through the L1A LUN

!Output Parameters:
                None

Return values:
		MODIS_E_BAD_INPUT_ARG	If any pointer argument is NULL
		MODIS_E_GEO		If any subroutine fails
		PGS_S_SUCCESS		Otherwise

Externally defined:
		MAPIOK			mapi.h
		MODIS_E_BAD_INPUT_ARG	PGS_MODIS_35251.h
		MODIS_E_GEO		PGS_MODIS_35251.h
		NUM_TEMPS               GEO_parameters.h
		PGS_S_SUCCESS		PGS_SMF.h
		SUCCEED			hdf.h
		TEMP_FVALUE             GEO_geo.h

Called by:      main

Routines called:
	 openMODISfile                  "mapi.h"
	 GEO_prepare_l1a_data           "GEO_input.h"
	 GEO_locate_one_scan            "GEO_output.h"
	 GEO_initialize_product         "GEO_output.h"
	 GEO_set_T_inst2sc              "GEO_earth.h"
	 GEO_write_granule_metadata     "GEO_output.h"
	 closeMODISfile                 "mapi.h"
         modsmf                         "smfio.h"
	 SWopen                         "HdfEosDef.h"
	 SWclose                        "HdfEosDef.h"
	 createMAPIfilehandle           "mapi.h"
	 releaseMAPIfilehandle          "mapi.h"

!Revision History:
 * $Log: GEO_locate_one_granule.c,v $
 * Revision 6.4  2013/06/13 21:21:20  jkuyper
 * Expanded buffer to avoid buffer overflow.
 *
 * Revision 6.3  2011/12/12 23:27:17  kuyper
 * Changed to not open the Geolocation product file unless and until metadata
 * is successfully extracted from the L1A file.
 * Changed to fail if a fatal error occurs in any scan.
 *
 * Revision 6.2  2011/02/14 21:07:28  kuyper
 * Corrected const-qualification of *geo_params.
 *
 * Revision 6.1  2010/06/07 17:59:25  kuyper
 * Helped resolve Bug 17 by initializing fill value for raw_mir_enc.
 * Helped resolve Bug 2470 by dropping maneuver_list.
 * Changed to expect status codes from all subroutines.
 * Changed order of parameters passed to GEO_write_granule_metadata().
 *
 * Revision 5.1  2004/09/10 14:41:18  vlin
 * 1. input parameter maneuver_list added
 * 2. initialize l1a_data.fill_values with fill values
 * 3. initialize l1a_data.temperatures with fill values
 * 4. Update calls to functions GEO_initialize_product,
 *    GEO_locate_one_scan, and GEO_write_granule_metadata.
 *
 * Revision 4.2  2004/06/01 16:35:43  kuyper
 * Changed to treat single-scan failures as non-fatal.
 *
 * Revision 4.1  2003/02/21 21:14:17  kuyper
 * Corrected to use void* pointers with %p format code.

Requirements:
                PR03-F-2.1-1
                PR03-F-2.1-2
                PR03-F-4.1-1
                PR03-F-4.1-2
                PR03-F-4.2-1
                PR03-F-4.2-2
                PR03-F-4.3-1
                PR03-F-4.3-2

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/
{
    GEO_bcoord_struct bounding_coords;
    l1a_data_struct l1a_data;
    MODFILE *l1a_file; /* MAPI structure for the Level 1A file */
    qa_metadata_struct qa_metadata;
    int         i;
    int32	swfid;		/* HDF-EOS file ID. */
    int retval = PGS_S_SUCCESS;
    char msgbuf[1024];
    const char *function = NULL;
    const char *argument = "";
    char filefunc[] = __FILE__ ", GEO_locate_one_granule";
    time_t tnow;

    /* Begin program logic */
    if (l1a_file_name == NULL || geo_file_name==NULL || geo_params==NULL)
    {
      sprintf(msgbuf,
	"l1a_file_name = %p, geo_file_name = %p, geo_params = %p",
	(void *)l1a_file_name, (void *)geo_file_name, (void *)geo_params);
      modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);

      return MODIS_E_BAD_INPUT_ARG;
    }

    /* Open the input L1A granule file */
    if ((l1a_file = openMODISfile(l1a_file_name, "a")) == NULL)
    {
      sprintf(msgbuf, "openMODISfile(\"%s\")", l1a_file_name);
      modsmf(MODIS_E_GEO, msgbuf, filefunc);

      return MODIS_E_GEO;
    }
    
    /* Reset QA counters.	*/
    qa_metadata.no_of_pixels = 0;
    qa_metadata.missingdata = 0;
    qa_metadata.outofboundsdata = 0;
    memset(qa_metadata.cumulated_gflags, 0,
	sizeof(qa_metadata.cumulated_gflags));

    /* Initialize boundaries.	*/
    bounding_coords.northcoord = -0.5*PGS_PI;
    bounding_coords.southcoord = 0.5*PGS_PI;
    bounding_coords.eastcoord = PGS_PI;
    bounding_coords.westcoord = -PGS_PI;
    bounding_coords.easthemi_ebc = 0.0;
    bounding_coords.easthemi_wbc = PGS_PI;
    bounding_coords.westhemi_ebc = -PGS_PI;
    bounding_coords.westhemi_wbc = 0.0;

    l1a_data.fill_values.scan_number = (int16)0;
    l1a_data.fill_values.EV_start_time = (float64)-2E9;
    l1a_data.fill_values.SD_start_time = (float64)-2E9;
    l1a_data.fill_values.SV_start_time = (float64)-2E9;
    l1a_data.fill_values.SD_sun_zenith = (float32)GEO_DOUBLE_FILLVALUE;
    l1a_data.fill_values.SD_sun_azimuth = (float32)GEO_DOUBLE_FILLVALUE;
    l1a_data.fill_values.moon_vector = (float32)GEO_DOUBLE_FILLVALUE;
    l1a_data.fill_values.sun_ref = (float32)GEO_DOUBLE_FILLVALUE;
    l1a_data.fill_values.mirr_side = (uint16)-1;
    l1a_data.fill_values.raw_mir_enc = MAX_UINT16_VAL;
    l1a_data.fill_values.impulse_enc = (float64)GEO_DOUBLE_FILLVALUE;
    l1a_data.fill_values.impulse_time = (float64)GEO_DOUBLE_FILLVALUE;
    l1a_data.fill_values.L1_scan_quality = (int32)-1;
    l1a_data.fill_values.geo_scan_quality = (int8)-127;
    l1a_data.fill_values.EV_center_time = (float64)-2E9;
    l1a_data.fill_values.orb_pos = (float64)GEO_DOUBLE_FILLVALUE;
    l1a_data.fill_values.orb_vel = (float64)GEO_DOUBLE_FILLVALUE;
    l1a_data.fill_values.T_inst2ECR = (float64)GEO_DOUBLE_FILLVALUE;
    l1a_data.fill_values.attitude_angels =(float64)GEO_DOUBLE_FILLVALUE;
    l1a_data.fill_values.EV_frames = (int32)0;
    l1a_data.fill_values.SD_frames = (int32)0;
    l1a_data.fill_values.SV_frames = (int32)0;
    l1a_data.fill_values.num_impulse = (uint8)0;
    strcpy((char*)l1a_data.fill_values.Scan_type, "");
    for (i=0; i<NUM_TEMPS; i++)
	l1a_data.temperatures[i] = (float32)TEMP_FVALUE;

    /* Read and validate geolocation data from the L1A file */
    if(GEO_prepare_l1a_data(l1a_file, geo_params, &l1a_data)
	!= PGS_S_SUCCESS)
    {
	function = "GEO_prepare_l1a_data";
	argument = l1a_file_name;
    }
    else if ( FAIL == (swfid = SWopen(geo_file_name, DFACC_CREATE)) )
    {
	function = "SWopen";
	argument = geo_file_name;
    }
    else
    {	/* Open the output geolocation product file */
	MODFILE *geo_file = createMAPIfilehandle(swfid);

	if(NULL == geo_file)
	{
	    function = "createMAPIfilehandle";
	    argument = geo_file_name;
	}
	else
	{
	    int scan_number;

	    if(GEO_initialize_product(l1a_data.num_scans,
		&l1a_data.fill_values, geo_file, swfid, geo_params)
		!= PGS_S_SUCCESS)
	    {
		function = "GEO_initialize_product";
		argument = l1a_file_name;
	    }
	    else if(GEO_set_T_inst2sc(&geo_params->coord_trans,
		&l1a_data.ECS_metadata) != PGS_S_SUCCESS)
		function = "GEO_set_T_inst2sc";
	    else for (scan_number = 0; scan_number < l1a_data.num_scans;
		    scan_number++)
	    {
              if (scan_number % 10 == 0) {
                  time(&tnow);
                  printf("scan: %d out of %d %s",
                         scan_number, l1a_data.num_scans, asctime(localtime(&tnow)));
                }

		/* Geolocate one scan */
		PGSt_SMF_status status = GEO_locate_one_scan(geo_params,
		    &l1a_data, scan_number, &qa_metadata, &bounding_coords,
		    geo_file);
		if(status != PGS_S_SUCCESS)
		{
		    sprintf(msgbuf, "GEO_locate_one_scan(%d)", scan_number);
		    modsmf(MODIS_E_GEO, msgbuf, filefunc);
		    /* No fatal errors are currently returned, but we plan to
		     * add some in the future. */
		    if(PGS_SMF_TestFatalLevel(status) != PGS_FALSE)
			retval = MODIS_E_GEO;
		} 
	    }

	    qa_metadata.retval = retval;

	    /* Generate and write output metadata to geolocation and L1A files*/
	    if(GEO_write_granule_metadata(geo_file, l1a_file, geo_params, 
		&bounding_coords, version, &l1a_data, &qa_metadata)
		!= PGS_S_SUCCESS)
	    {
		sprintf(msgbuf, "GEO_write_granule_metadata(\"%s\",\"%s\")",
		    geo_file_name, l1a_file_name);
		modsmf(MODIS_E_GEO, msgbuf, filefunc);

		retval = MODIS_E_GEO;
	    }

	    if(MAPIOK != releaseMAPIfilehandle(&geo_file))
	    {
		sprintf(msgbuf, "releaseMAPIfilehandle(\"%s\")", geo_file_name);
		modsmf(MODIS_E_GEO, msgbuf, filefunc);
		
		retval = MODIS_E_GEO;
	    }
	}

	if(SUCCEED != SWclose(swfid))
	{
	    sprintf(msgbuf, "SWclose(\"%s\")", geo_file_name);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);
	  
	    retval = MODIS_E_GEO;
	}
    }

    if(function)
    {
	sprintf(msgbuf, "%s(\"%s\")", function, argument);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
      
	retval = MODIS_E_GEO;
    }

    if(closeMODISfile(&l1a_file)!=MAPIOK)
    {
	sprintf(msgbuf, "closeMODISfile(\"%s\")", l1a_file_name);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
      
	retval = MODIS_E_GEO;
    }

    return retval;
}
