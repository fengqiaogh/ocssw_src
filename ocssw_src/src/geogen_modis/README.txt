1. README for MOD_PR03 (V6.1.0)



2. POINTS OF CONTACT:

   Developers:
   ----------
   James Kuyper Jr.
   MODIS Science Data Support Office
   10210 Greenbelt Road, Room 563
   Lanham, MD 20706
   tel: 301-867-2064, fax: 301-867-6246
   e-mail: James.R.Kuyper@nasa.gov

3. BUILD:
  a) This software optionally uses the IMSL library, if it is available. If it
      is used, you need to be careful to protect the LM_LICENSE_FILE
      environment variable against changes by the IMSL setup script. Type:

        setenv OLD_LM $LM_LICENSE_FILE
        source /usr/local/lib/vni/ipt/bin/iptsetup.csh
        setenv LM_LICENSE_FILE $OLD_LM\:$LM_LICENSE_FILE
        unsetenv OLD_LM
        source $VNI_DIR/ipt/bin/IRIX64.csh

     The *.csh scripts  will set the AC_INCLUDE and LINK_CNL_STATIC environment
     variables.


  b) Set all of the remaining environment variables required by the makefile.
     At the TLCF, this may be done by typing:

        source /data/common/util/bin/setup-5.2.19.csh -f90

     NOTE: The IMSL setup scripts interact with these scripts through the CC,
     CFLAGS, and PATH environment variables. As a result, performing this step
     before step a) would cause the build to fail.


  c) Be sure to use matching up-to-date SDPTK message and include files
     PGS_MODIS_35251.h and PGS_35251. The PGS_MODIS_35251.h file must be on the
     include file search path when compiling the code. The PGS_35251 file must
     be located in the $PGSMSG directory when running the code.

     If those files doe not match, build PGS_MODIS_35251.h and PGS_35251 using
     the SDPTK smfcompile facility and the Geolocation MODIS_35251.t seed file
     included with this delivery.  At the UNIX command line type

	$PGSBIN/smfcompile -f MODIS_35251.t -r

     The -r option redirects the message file to the $PGSMSG directory. The
     PGS_MODIS_35251.h file will be created in the current directory, which is
     acceptable if that's the source code directory.

     NOTE:
     The PGS_MODIS_39501.h and PGS_39501 message file also need to be
     accessible.  They should be installed as part of the SDST Toolkit, with
     PGS_MODIS_39501.h being found in $SDST_INC, and PGS_39501 being found in
     $PGSMSG.
     The PGS_MODIS_35000.h and PGS_35000 message file also need to be
     accessible.  They should come with the SDST in-house library, under the
     shared_source directory.


  d) To build the geolocation code, type

        make -f MOD_PR03.mk MOD_PR03.exe |& fold -s > make.run

     The 'fold -s' command is used solely to render the make.run file more
     readable when using text editors that don't automatically wrap long lines,
     or which allow wrapping to occur in the middle of a word.


  e) Set the SDPTK environment variable PGS_PC_INFO_FILE to the full path name
     of the PCF, either MOD_PR03.pcf or MYD_PR03.pcf, depending upon which
     satellite is the source of the input data. Edit the PCF to reflect the
     names and directories of the locations of input and output files on the
     host system. There must be entries for both the MOD03LUT.coeff and the
     MYD03LUT.coeff in the PCF, regardless of which satellite is being used.
     For each of the other LUNS listed below, the correct version for each
     satellite, and only that version, must be listed in the PCF.

     Files with "*" represent file naming conventions for dynamic files,
     rather than specific file names:

     LUN     Description		Terra		Aqua
     ======  ========================	============	============
     10501   Real Ephemeris files	AM1EPH*		PM1EPH*
     10501   Simulated Ephemeris files	EOSAM1_*.eph	EOSPM1_*eph
     10502   Real Attitude files	AM1ATT*		PM1ATT*
     10502   Simulated Attitude files	EOSAM1_*.att	EOSPM1_*.att
     500100  L1A input files		MOD01.*		MYD01.*
     500500  L1A MCF			MOD01.mcf	MYD01.mcf
     600000  Geolocation output files	MOD03.*		MYD03.*
     600020  Geolocation parameters	MOD03LUT.coef*	MYD03LUT.coef*
     600111  GEO MCF			MOD03.mcf	MYD03.mcf

     The LocalVersionID configuration parameter(LUN 600001) should be set to
     match the the version of the Geolocation filespec that is implemented by
     this version of the code.

     The maneuver list configuration parameter (LUN 600002) must be set to
     contain a space delimited list of the start and stop times for any
     maneuvers that overlap the processing time, in CCSDS format. Example:
     2010-06-01T15:30:00.000000Z 2010-06-01T18:03:09.0000000Z

     The parameter file revision configuration parameter (LUN 600021) must be
     set to match the numeric value of the RCS Revision keyword from the
     appropriate parameter file.

     The EA Source configuration parameter (LUN 600280) must be "SDP Toolkit"
     for Aqua input files; it can be either "SDP Toolkit" or "MODIS Packet" for
     Terra input files.

     The orbit validation configuration parameter (LUN 600281) must be either
     "TRUE" or "FALSE". It needs to be "FALSE" when running in Near Real Time
     mode, because the validation criteria it uses are too strict for the 
     low-precision estimation methods we use in that mode.

     The PGEVersion configuration parameter (LUN 800500) should be set to
     the correct PGE version.

     The SatelliteInstrument parameter (LUN 800510) must be set to AM1M for
     Terra data, or PM1M for Aqua data.

     The ProcessingEnvironment configuration parameter (LUN 800550) should be
     set to the output of a UNIX 'uname -a' command, as executed in the
     processing environment.

     The ReprocessingActual (LUN 800600) and ReprocessingPlanned (LUN 800605)
     should be set to the values that should appear in the corresponding ECS
     metadata items of the same name. Typically, that would be 
     "reprocessed once" and "further update anticipated", respectively. In
     LANCE-MODIS processing, ReprocessingActual is set to "Near Real Time".

     NOTE: The SDP Toolkit time conversion utilities will produce several
     warnings per pixel if there is any use of predicted UT1-UTC values. Since
     the UT1-UTC value for a given date remains 'predicted' for several weeks
     after that date, our forward production runs will always be using predicted
     UT1-UTC values, which would result in Gigabytes of status log files. To
     avoid this, our PCF file differs from the template version, by replacing
     '0' with 'PGSCSC_W_PREDICTED_UT1' for the 10119 entry.
     PGS_DEM_M_FILLVALUE_INCLUDED is also on the message suppression list, since
     some versions of the DEM files have all open ocean pixels marked with
     fill values.

  f) MOD_PR03 updates the ECS geolocation metadata in it's input L1A product
     file.  That is why the "input" file is listed in the output section of the
     PCF, and is also why that file must be writable. Make sure that you copy
     it into that location from a safe original before each run of MOD_PR03.
     Never use the original input file directly, and never re-use the output
     file.

  g) Type:

        MOD_PR03.exe
    
  h) Normal status log warnings:
  The following messages, while apparently indicating a problem, are completely
  normal for a MOD_PR03 run.

  If no maneuvers occurred during the time period covered by the run, it's
  perfectly acceptable for LUN 600002 to be empty. However, if it is, the
  following warning messages will be generated:

      PGS_PC_GetPCSDataRetrieveData():PGSPC_W_NO_CONFIG_VALUE
      PGS_PC_GetPCSData():PGSPC_W_NO_CONFIG_VALUE
      PGS_PC_GetConfigData():PGSPC_W_NO_CONFIG_FOR_ID

  Those messages can be avoided, if you wish, by providing start and end times
  for a fake maneuver that doesn't overlap the time period of the run.

  As of PGE01 v6.0.5, we are currently using a new 15 arc second Digital
  Elevation Model (DEM) that does not cover Greenland or areas south of 65S.
  Inside those areas, MOD_PR03 has set up the SDP Toolkit's DEM routines to
  automatically fall back to using the 30 arcsec DEM. When fallback occurs, the
  following messages will appear in in the status log:

      PGS_DEM_RecursiveSearchDeg():PGSDEM_E_CANNOT_ACCESS_DATA
      PGS_DEM_GetRegion():PGSDEM_M_MULTIPLE_RESOLUTIONS

  We're also using a new 30 arc second DEM that does cover those areas, and
  should be an improvement over the old one.

  There should also be four copies of the following pair of messages for every
  successfully geolocated output granule:

	PGS_CSC_ECRtoGEO():PGSCSC_W_INVALID_ALTITUDE
	PGS_CSC_GrazingRay():PGSCSC_W_HIT_EARTH 

  The position being converted is one meter from the center of the earth, which
  is considered an invalid altitude, but is processed correctly. The ray being
  extrapolated is supposed to hit the Earth. If either of these messages is NOT
  observed, something has gone wrong, and there should be other messages
  telling you what that is.


4. ENVIRONMENT VARIABLES:
   used by MOD_PR03.mk
   See "!Env Variables" in MOD_PR03.mk prologue for list of system environment
   variables.



5. MAPI: version 6.0.2c7
 


6. HDF: version 4.2r10
 


7. SDP TOOLKIT: version 5.2.19v1.01



8. This code was developed on an Intel Xeon processor using CentOS Linux.




9. ANCILLARY DATA:

  a) AM1EPHN0#001010120011400000000000, AM1EPHN0#001010120011600000000000,
     AM1ATTN0#001010120011400000000000, AM1ATTN0#001010120011600000000000,
     PM1EPHND.P2002193.1200.001.2002194114016,
     PM1EPHND.P2002194.1200.001.2002195150308,
     PM1ATTNR.P2002194.1000.001.2002194155029,
     PM1ATTNR.P2002194.1200.001.2002195150512

    These are Spacecraft Ephemeris and Attitude files created by EDOS.
  b) leapsec.dat
    This file must have actual leap second dates listed up to and beyond the
    date covered by the input L1A product file.
  c) utcpole.dat 
    Ideally, to accurately simulate the production environment, this file
    should have only predicted UTC-UT1 corrections for the time period covered
    by the L1A input file. In production, actual values for those corrections
    will not be available until three weeks after the processing date. Messages
    complaining about the use of predicted values are supposed to be suppressed
    by the 10119 entry in the PCF, but you won't be able to verify that if the
    utcpole.dat file doesn't have them.
  d) DEM Elevation data files in SDP Toolkit format:
    30 arc second resolution:
    dem30ARC_E60N90.hdf  dem30ARC_W60N90.hdf  dem30ARC_W180N90.hdf
    dem30ARC_E60N0.hdf   dem30ARC_W60N0.hdf   dem30ARC_W180N0.hdf.

    15 arc second resolution:
    dem15ARC_W180N90.hdf dem15ARC_W120N90.hdf dem15ARC_W60N90.hdf
    dem15ARC_E0N90.hdf   dem15ARC_E60N90.hdf  dem15ARC_E120N90.hdf
    dem15ARC_W180N45.hdf dem15ARC_W120N45.hdf dem15ARC_W60N45.hdf
    dem15ARC_E0N45.hdf   dem15ARC_E60N45.hdf  dem15ARC_E120N45.hdf
    dem15ARC_W180N0.hdf  dem15ARC_W120N0.hdf  dem15ARC_W60N0.hdf
    dem15ARC_E0N0.hdf    dem15ARC_E60N0.hdf   dem15ARC_E120N0.hdf
    dem15ARC_W180S45.hdf dem15ARC_W120S45.hdf dem15ARC_W60S45.hdf
    dem15ARC_E0S45.hdf   dem15ARC_E60S45.hdf  dem15ARC_E120S45.hdf

  e) All SDP Toolkit message files.
  f) M-API message files: PGS_39602, PGS_39603
  g) SDST message files: PGS_39501, PGS_39604


10. MODIS INPUT PRODUCTS:
    L1A files: MOD01 or MYD01


11. OTHER inputs:
    MOD03LUT.coeff and MYD03LUT.coeff
    maneuver_terra.coeff or maneuver_aqua.coeff


12. MODIS OUTPUT PRODUCTS:
    Geolocation files: MOD03 or MYD03


13. PROBLEMS:
  a) A granule which has a long stretch of bad data while passing near either
    pole, and crossing either the Greenwich meridian or the international
    dateline, can be incorrectly identified as a polar granule. This produces
    an unnecessarily large (but still valid) bounding rectangle.

  b) A granule containing no geolocatable pixels locations will have bounding
    coordinates with a north boundary of -90 degrees, and a south boundary of
    +90 degrees. The archiving system will reject such such a bounding box if
    it appears in inventory metadata.
    The data model does not allow for the use of any bounding box for a
    non-geolocatable granule, that can't be mistaken for the bounding box of a
    correctly geolocated granule.

  c) If for any reason the correct GRing cannot be calculated, a GRing
    consisting of 4 points, 0.01 degrees long on each side, with the southwest
    corner at 0N 0E, will be produced instead. This looks like a valid GRing,
    but does not correctly represent the geographic extent of the granule.

  d) The occurrence of IMSL errors is noted in the status log, but details of
    IMSL error messages are ignored.

  e) Some memory leaks have been found in the code package, but none of them
    can be traced directly to the science code. All the ones we've successfully
    traced were due to memory allocation occurring inside SDP Toolkit functions.

  f) If the output from MOD_PR03 is to be compared with simulated Geolocation
    product file, it is important to use the exact same leap seconds (TAI-UTC)
    file, and the exact same polar motion and UTC-UT1 file for all three of the
    relevant programs: the orbsim program used to create the simulated
    ephemeris/attitude files, the simulator used to create the simulated L1A
    and Geolocation product files, and MOD_PR03.exe itself. If any one of these
    programs is run using a different leap seconds (TAI-UTC) file or polar
    motion and UTC-UT1 file from the other two, discrepancies may occur between
    the simulated Geolocation product files, and the actual Geolocation product
    files created from the simulated L1A product files.

  g) The reported Terra spacecraft orbit does not pass our original validation
    tests.  We've been unable to get an answer as to whether the actual orbit,
    the reported orbit, or our validation criteria are incorrect. In the
    meantime, we've relaxed our validation criteria to accept the reported
    orbit.

  h) The T_inst2SD matrix is believed to be incorrect. We're still waiting for
    details on the corrected values for this matrix. In the meantime, as
    requested by MCST, the algorithm for calculating solar diffuser angles has
    been changed to ignore our estimates of post-launch shifts in the MODIS
    orientation relative to the Aqua spacecraft.

  j) Due to the hurried delivery of this code, there are likely to be an
    unknown number of additional problems that should have been caught during 
    testing. See section 16.


14. DELIVERY SIZE:
  The total size of this delivery is 881 MB. A total of ~2GB will be needed to
  compile and run the test.



15. OUTPUT FILE SIZE:
  The expected HDF output Geolocation product file is 89.3 MB.
  The L1A product file's size is expected to increase by 31 KB for these test
  files.
  The actual size of the increase varies from file to file.


16. TESTS PERFORMED:
  This code has been regression tested against the previous version.

17. EXIT CODES:
  MOD_PR03 will return an exit code of 0 if it succeeds, and 1 if it fails for
  any granule. In general, failure within a given granule does NOT terminate
  processing of that granule unless necessary, nor does it prevent processing
  of any later ones.


18. ERROR LIST:

These are the error messages that the MOD_PR03 Geolocation process may
generate directly.  Error messages that may be generated by supporting
utility routines (PGS Toolkit and M-API) are not included in this listing.

Each message is identified by its Message mnemonic and the module name(s)
that may generate it.  Both pieces of information are included as part of
the message in the LogStatus listing.  A brief (hopefully not too inscrutable)
explanation of why the process would generate the message is given.  Sometimes
there is a different explanation for each module that might generate the
message.  In other cases, the explanation for the message is the same
irrespective of the routine generating it.

The recommended action to take should a particular message occur is described
in MOD_PR03's Fault Isolation Procedures (FIPs).  For each message, there
is one or more recommended FIP referenced that will describe what action to
take.  In some cases, the FIP referenced will depend upon which module
generated the message.  The FIP selected should normally be based upon the
first error message MOD_PR03 generates.


Message mnemonic	module	Explanation	Fault Isolation Procedure

MODIS_E_BAD_DEM        GEO_read_DEM()
	Height retrieved from Digital Elevation Model is
	outside the expected Range				MOD_PR03_FIP-1

MODIS_E_BAD_INPUT_ARG Function arguments failed input checks	MOD_PR03_FIP-1
		GEO_create_swath()
		GEO_del_limit_check()
		GEO_ellip_position()
		GEO_cumulate_GRing()
		GEO_get_GRing_points()
		GEO_get_T_inst2ecr()
		GEO_get_bounding_coords()
		GEO_get_ephatt_inputs()
		GEO_get_utcpole_metadata()
		GEO_get_version_metadata()
		GEO_in_maneuver()
		GEO_initialize_product()
		GEO_interp_mirr_enc()
		GEO_interp_ephemeris_attitude()
		GEO_landsea_mask()
		GEO_locate_one_granule()
		GEO_locate_one_scan()
		GEO_prepare_ancil_data()
		GEO_prepare_l1a_data()
		GEO_read_L1AECS_metadata()
		GEO_read_L1Apacket_data()
		GEO_read_L1Ascan_metadata()
		GEO_read_L1Aspecific_metadata()
		GEO_read_maneuver_file()
		GEO_read_param_file()
		GEO_set_T_inst2sc()
		GEO_solar_and_lunar_vectors()
		GEO_terrain_correct()
		GEO_update_L1A_metadata()
		GEO_validate_derived_products()
		GEO_validate_earth_location()
		GEO_write_ECS_metadata()
		GEO_write_geospecific_metadata()
		GEO_write_granule_metadata()
		GEO_write_one_scan()
		GEO_write_parameters()
		GEO_write_scan_data()
		GEO_write_scan_metadata()

MODIS_E_BAD_SCAN        GEO_read_param()
	Unix error reading Geolocation Parameter text file	MOD_PR03_FIP-1

MODIS_E_BAD_VIEW_VEC	GEO_ellip_position()			MOD_PR03_FIP-1

MODIS_E_CREATE_MODIS_ARRAY
	GEO_write_parameters()
		Unable to create output SDS's			MOD_PR03_FIP-1

MODIS_E_DATA_SCAN       GEO_prepare_mirr_data()
	No valid mirror encoder data for scan			MOD_PR03_FIP-11

MODIS_E_DEM_DATA_SIZE   GEO_read_DEM()
	DEM data of unexpected data size			MOD_PR03_FIP-6

MODIS_E_DEM_IN_MEM      GEO_latlon2height()
	Requested DEM data not in memory			MOD_PR03_FIP-1

MODIS_E_DEM_METADATA    GEO_initialize_DEM()
	Unexpected DEM position or data units			MOD_PR03_FIP-6

MODIS_E_GEO	Error indication returned by subroutine

	GEO_create_swath()     
		Unable to create Geolocation swath.		MOD_PR03_FIP-5
		Can't detach from Geolocation swath.		MOD_PR03_FIP-5
	GEO_cumulate_GRing()
		PGS_CSC_quatRotate failed			MOD_PR03_FIP-1
		PGS_CSC_Norm failed				MOD_PR03_FIP-1
		imsl_d_lin_sol_posdef failed			MOD_PR03_FIP-1
	GEO_derived_products()
		PGS_CBP_Earth_CB_Vector				MOD_PR03_FIP-15
		PGS_CSC_ECItoECR				MOD_PR03_FIP-17
	GEO_ellip_position()
		PGS_CSC_ECRtoGEO failed				MOD_PR03_FIP-11
		PGS_CSC_GetEarthFigure				MOD_PR03_FIP-16
	GEO_get_ephatt_inputs()
		PGS_IO_Gen_Open()				MOD_PR03_FIP-9
		fread()						MOD_PR03_FIP-11
		PGS_IO_Gen_Close()				MOD_PR03_FIP-11
	GEO_get_geoid()
		PGS_DEM_GetSize failed				MOD_PR03_FIP-6
		PGS_DEM_GetQualityFlags failed			MOD_PR03_FIP-6
	GEO_get_GRing_points()
		PGS_CSC_GrazingRay failed			MOD_PR03_FIP-1
		PGS_CSC_ECRtoGEO failed				MOD_PR03_FIP-11
	GEO_get_T_inst2ecr()
		PGS_CBP_Earth_CB_Vector				MOD_PR03_FIP-15
		PGS_CSC_ECItoECR				MOD_PR03_FIP-17
		PGS_CSC_EulerToQuat				MOD_PR03_FIP-1
		PGS_CSC_quatRotate				MOD_PR03_FIP-1
		PGS_TD_TAItoUTC()				MOD_PR03_FIP-18
	GEO_get_utcpole_metadata()
		PGS_IO_Gen_Open failed				MOD_PR03_FIP-19
		PGS_TD_UTCtoUTCjd failed			MOD_PR03_FIP-19.2
		PGS_CSC_UTC_UT1Pole failed			MOD_PR03_FIP-19.2
	GEO_get_version_metadata()
		PGS_PC_GetConfigData failed			MOD_PR03_FIP-27
	GEO_get_view_vec()
		GEO_get_inst_mirr_normal failed			MOD_PR03_FIP-3
	GEO_in_maneuver()
		PGS_PC_GetConfigData()				MOD_PR03_FIP-32
		PGS_TD_UTCtoTAI()				MOD_PR03_FIP-18
	GEO_initialize_product()
		createMODISarray failed				MOD_PR03_FIP-5
		putMODISarinfo failed				MOD_PR03_FIP-5
		putMODISdimname failed				MOD_PR03_FIP-5
		GEO_create_swath failed				MOD_PR03_FIP-3
		VSattach failed					MOD_PR03_FIP-5
		VSdetach failed					MOD_PR03_FIP-5
		VSfdefine failed				MOD_PR03_FIP-5
		VSsetattr failed				MOD_PR03_FIP-5
		VSsetfields failed				MOD_PR03_FIP-5
		VSsetname failed				MOD_PR03_FIP-5
	GEO_interp_ephemeris_attitude()
		GEO_poly_coef1 failed				MOD_PR03_FIP-1
		GEO_poly_fit failed				MOD_PR03_FIP-3
		PGS_CSC_EulerToQuat()				MOD_PR03_FIP-1
		PGS_EPH_EphemAttit()				MOD_PR03_FIP-20
		PGS_PC_GetConfigData()				MOD_PR03_FIP-12
		PGS_TD_TAItoUTC()				MOD_PR03_FIP-18
	GEO_interp_mirr_enc()
		imsl_d_spline_interp()				MOD_PR03_FIP-1
		imsl_d_spline_value()				MOD_PR03_FIP-1
	GEO_locate_one_granule()
		Unable to open L1A file				MOD_PR03_FIP-5
		Unable to open Geolocation file			MOD_PR03_FIP-5
		Unable to create M-API handle			MOD_PR03_FIP-1
		Unable to prepare L1A data			MOD_PR03_FIP-3
		Unable to initialize Geolocation product	MOD_PR03_FIP-3
		Unable to set T_inst2sc				MOD_PR03_FIP-3
		Unable to locate one scan			MOD_PR03_FIP-3
		Unable to write granule metadata		MOD_PR03_FIP-3
		Unable to release M-API handle			MOD_PR03_FIP-1
		Unable to close Geolocation file		MOD_PR03_FIP-1
		Unable to close L1A file			MOD_PR03_FIP-1
	GEO_locate_one_scan()
		GEO_get_bounding_coords failed			MOD_PR03_FIP-1
		GEO_earth_location failed			MOD_PR03_FIP-3
		GEO_validate_earth_location failed		MOD_PR03_FIP-1
		GEO_derived_products failed			MOD_PR03_FIP-3
		GEO_interp_ephemeris failed			MOD_PR03_FIP-3
		GEO_get_T_inst2ecr failed			MOD_PR03_FIP-3
		GEO_validate_derived_products failed		MOD_PR03_FIP-1
		GEO_write_one_scan failed			MOD_PR03_FIP-3
		GEO_cumulate_GRing failed			MOD_PR03_FIP-3
		PGS_TD_TAItoUTC()				MOD_PR03_FIP-18
	GEO_prepare_ancil_data()
                PGS_TD_EOSAMtoTAI failed			MOD_PR03_FIP-2
	GEO_prepare_l1a_data()
		GEO_read_L1Aspecific_metadata failed		MOD_PR03_FIP-3
		GEO_read_L1AECS_metadata failed			MOD_PR03_FIP-3
		GEO_read_L1A_scan_metadata failed		MOD_PR03_FIP-3
		GEO_read_L1A_packet_metadata failed		MOD_PR03_FIP-3
		GEO_prepare_mirr_data failed			MOD_PR03_FIP-1
		GEO_prepare_ancil_data failed			MOD_PR03_FIP-12
	GEO_read_DEM()
		GEO_DEMalloc failed				MOD_PR03_FIP_19
		GEO_get_geoid failed				MOD_PR03_FIP_6
		PGS_DEM_GetSize failed				MOD_PR03_FIP_10
	GEO_read_L1AECS_metadata()
		getMODISECSinfo failed				MOD_PR03_FIP-11
	GEO_read_L1Apacket_data()
		getMODISarray failed				MOD_PR03_FIP-1
		getMODIStable failed				MOD_PR03_FIP-8
	GEO_read_L1Ascan_metadata()
		getMODISarray failed				MOD_PR03_FIP-1
	GEO_read_L1Aspecific_metadata()
		getMODISfileinfo failed				MOD_PR03_FIP-8
	GEO_read_param_file()
		PGS_PC_GetConfigData failed			MOD_PR03_FIP-28
		PGS_IO_GenOpen failed				MOD_PR03_FIP-13
		GEO_read_param failed				MOD_PR03_FIP-3
		malloc failed					MOD_PR03_FIP-29
		PGS_TD_UTCtoTAI failed				MOD_PR03_FIP-24
		PGS_IO_GenClose failed				MOD_PR03_FIP-11
	GEO_read_maneuver_file()
		PGS_IO_Gen_Open	failed				MOD_PR03_FIP-13
		PGS_PC_GetConfigData failed			MOD_PR03_FIP-31
		PGS_TD_UTCtoTAI failed				MOD_PR03_FIP-24
		sscanf failed					MOD_PR03_FIP-24
	GEO_set_T_inst2sc()
		PGS_TD_UTCtoTAI failed				MOD_PR03_FIP-18
	GEO_solar_and_lunar_vectors()
		GEO_interp_ephemeris_attitude failed		MOD_PR03_FIP-3
		GEO_get_T_inst2ecr failed			MOD_PR03_FIP-3
		GEO_vec_unit3 failed				MOD_PR03_FIP-1
		PGS_CBP_Earth_CB_Vector				MOD_PR03_FIP-15
		PGS_CSC_ECItoECR				MOD_PR03_FIP-17
		PGS_TD_TAItoUTC()				MOD_PR03_FIP-18
        GEO_update_L1A_metadata() 
                getMODISECSinfo failed                          MOD_PR03_FIP-11
                PGS_MET_init failed                             MOD_PR03_FIP-22
                PGS_MET_SetAttr failed                          MOD_PR03_FIP-22 
                PGS_MET_Write failed to write metadata          MOD_PR03_FIP-5
                substrMODISECSinfo failed                       MOD_PR03_FIP-1 
	GEO_write_ECS_metadata()
		PGS_MET_Init failed				MOD_PR03_FIP-22
		PGS_MET_SetAttr failed				MOD_PR03_FIP-22
		PGS_MET_Write failed				MOD_PR03_FIP-5
	GEO_write_geospecific_metadata()
		PGS_PC_GetConfigData failed			MOD_PR03_FIP-12
		putMODISfileinfo failed				MOD_PR03_FIP-5
	GEO_write_granule_metadata()
		GEO_get_GRing_points failed			MOD_PR03_FIP-11
		GEO_get_utcpole_metadata failed			MOD_PR03_FIP-3
		GEO_get_version_metadata failed			MOD_PR03_FIP-3
		GEO_update_L1A_metadata failed			MOD_PR03_FIP-3
		GEO_write_ECS_metadata failed			MOD_PR03_FIP-3
		GEO_write_geospecific_metadata failed		MOD_PR03_FIP-3
		GEO_write_parameters failed			MOD_PR03_FIP-3
		PGS_EPH_GetEphMet failed			MOD_PR03_FIP-9
		PGS_PC_GetUniversalRef failed			MOD_PR03_FIP-7
		PGS_PC_GetConfigData failed			MOD_PR03_FIP-12.1
	GEO_write_input_data()
		GEO_get_ephatt_inputs()				MOD_PR03_FIP-3
	GEO_write_one_scan()
		GEO_landsea_mask failed				MOD_PR03_FIP-10
	main()
		GEO_initialize_product failed			MOD_PR03_FIP-10
		PGE_DEM_Open failed				MOD_PR03_FIP-10
		GEO_close_DEM failed				MOD_PR03_FIP-11
		Failed to close DEM's Land/Sea mask		MOD_PR03_FIP-11
		Failed to geo-locate a MODIS granule		MOD_PR03_FIP-3
		Geolocation Parameter file ingest failed	MOD_PR03_FIP-3
		GEO_read_maneuver_file failed			MOD_PR03_FIP-3
		Call to PGS_PC_GetNumberOfFiles fails		MOD_PR03_FIP-7
		Call to PGS_PC_GetReference for L1A files fails	MOD_PR03_FIP-7
		Call to PGS_PC_GetReference for GEO files fails	MOD_PR03_FIP-7

MODIS_E_GEO_ANCIL_DATA
	GEO_prepare_ancil_data()				MOD_PR03_FIP-12
	Almost all spacecraft ancillary data from the L1A product is invalid.

MODIS_E_GEO_BAD_SIZE	GEO_get_geoid()
	The geoid data in the DEM file has the wrong size	MOD_PR03_FIP-6

MODIS_E_GEO_CB_VECTORS  GEO_write_one_scan()
	Error computing solar or lunar vectors for scan		MOD_PR03_FIP-3

MODIS_E_GEO_ECR2_LAT_LON GEO_terrain_correct()  
	Call to PGS_CSC_ECRtoGEO did not return PGS_S_SUCCESS	MOD_PR03_FIP-11

MODIS_E_GEO_ELLIP_POS   GEO_earth_location()
	Error in computing view vector's Earth ellipse position or terrain
	correction						MOD_PR03_FIP-3

MODIS_E_GEO_FORMATTER	GEO_prepare_l1a_data()
	Both sides of the Formatter electronics are turned on	MOD_PR03_FIP-1

MODIS_E_GEO_GET_STIME   GEO_get_view_vec()
	Failed to compute frame's sample time			MOD_PR03_FIP-1

MODIS_E_GEO_GET_VIEW_VEC	GEO_earth_location()
	Failed to compute a frame's view vectors		MOD_PR03_FIP-3

MODIS_E_GEO_INT_MIRR_ANG	GEO_get_inst_mirr_normal()
	Failed to interpolate the instrument's mirror angle	MOD_PR03_FIP-3

MODIS_E_GEO_INT_MIRR_ENC	GEO_get_inst_mirr_normal()
	Failed to interpolate the instrument's mirror encoder	MOD_PR03_FIP-3

MODIS_E_GEO_LATLON2HEIGHT	GEO_terrain_correct()
	Failed to retrieve terrain height for location		MOD_PR03_FIP-3

MODIS_E_GEO_MIRR_MOTION         GEO_interp_mirr_enc()
        Implausible mirror motion data                          MOD_PR03_FIP-1

MODIS_E_GEO_MISSING_INPUT	GEO_check_ea_headers()
	Needed Ephemeris/Attitude files were not staged		MOD_PR03_FIP-4

MODIS_E_GEO_NO_ZERO_ENCODER	GEO_interp_mirr_enc()
	Encoder times are probably corrupt                      MOD_PR03_FIP-1

MODIS_E_GEO_READ_DEM            GEO_terrain_correct()
	DEM data not retrieved					MOD_PR03_FIP-6

MODIS_E_GEO_SCANNO_INPUT	GEO_prepare_l1a_data()
	Invalid number of scans					MOD_PR03_FIP-1

MODIS_E_GEO_SWDEFINEDIM         GEO_create_swath()
	SWdefdim() returned FAIL defining dimension name	MOD_PR03_FIP-5

MODIS_E_GEO_SWFLDS      GEO_create_swath()
	Failed to set Swath metadata field names		MOD_PRO3_FIP-5

MODIS_E_GEO_VEC_UNIT3 Unable to generate unit vector from 0-vector
	GEO_get_T_inst2ecr()					MOD_PR03_FIP-1
	GEO_get_inst_mirr_normal()				MOD_PR03_FIP-1
	GEO_terrain_correct()					MOD_PR03_FIP-1

MODIS_E_GEO_WRITE_SCAN_DATA	GEO_write_one_scan()
	Error writing scan data					MOD_PR03_FIP-3

MODIS_E_GEO_WRITE_SCAN_MDATA	GEO_write_one_scan()
	Error writing scan metadata				MOD_PR03_FIP-3

MODIS_E_GEO_WRONG_PLATFORM
	GEO_read_param_file: File staged as Geolocation parameter file is
		for the wrong platform				MOD_PR03_FIP-25

	GEO_prepare_l1a_data: The L1A file staged is for the wrong platform
								MOD_PR03_FIP-26

MODIS_E_INSUFFICIENT_PKTS	GEO_interp_ephemeris_attitude()
	No valid spacecraft ancillary data available		MOD_PR03_FIP-12

MODIS_E_MISSING_OUTPUT  GEO_write_geospecific_metadata()
	Missing output metadata					MOD_PR03_FIP-1

MODIS_E_ONE_UNFLAG Insufficient data to perform measurement delta validation.
	GEO_del_limit_check()					MOD_PR03_FIP-11

MODIS_E_PGS_PC_GETCONFIGDATA
	GEO_terrain_correct() PGS_PC_GetConfigData failed to
		retrieve spacecraft Terrain correction
		performance selection runtime parameter		MOD_PR03_FIP-23

MODIS_E_PREMATURE_EOF
	GEO_read_param() Corrupted Geolocation parameters	MOD_PR03_FIP-24.1
	GEO_read_param_file() Corrupted Geolocation parameters	MOD_PR03_FIP-24.1
	GEO_get_utcpole_metadata() Corrupted utcpole.dat	MOD_PR03_FIP-19.2

MODIS_E_PUT_MODIS_ARINFO					MOD_PR03_FIP-5
	Failed to write local attribute to Geolocation product file.
	GEO_write_parameters()

MODIS_E_PUT_MODIS_ARRAY
	Failed to write data to Geolocation product file	MOD_PR03_FIP-5
	GEO_initialize_product()
	GEO_write_parameters()
	GEO_write_scan_data()

MODIS_E_PUT_MODIS_DIMNAME
	Failed to name SDS dimension in Geolocation product	MOD_PR03_FIP-5
	GEO_write_parameters()

MODIS_E_UNKNOWN_PARAMETER	GEO_interp_ephemeris_attitude()
	Ephemeris and Attitude source identifier runtime
	parameter's value is invalid				MOD_PR03_FIP-12

MODIS_E_UTIL_MAX_POLY   GEO_poly_fit()
	The degree of the polynomial GEO_poly_fit is to
	evaluate is too large					MOD_PR03_FIP-1

MODIS_E_UTIL_POLY
	GEO_poly_fit()	The degree of the polynomial
		GEO_poly_fit is to evaluate is invalid		MOD_PR03_FIP-1
	GEO_interp_mirr_ang()	Interpolation of the MODIS
		mirror's state failed				MOD_PR03_FIP-3

MODIS_E_UTIL_VEC        GEO_vec_unit3()
	Input vector has zero magnitude				MOD_PR03_FIP-11

MODIS_E_WRONG_FIELD     GEO_read_param()
	A geolocation parameter name in Geolocation parameter file was not
	valid							MOD_PR03_FIP-24

MODIS_E_WRONG_FILE
	GEO_read_param_file()
	File staged as Geolocation parameter file not recognized as
	correct parameter file					MOD_PR03_FIP-24

	GEO_read_maneuver_file()
	File staged as maneuver list not recognized as
	correct file						MOD_PR03_FIP-24

MODIS_E_WRONG_LUT
	GEO_read_param_file()
	The Geolocation parameter file has the wrong version, incompatible with
	the current PCF.					MOD_PR03_FIP-24.1
PGS_E_UNIX
	GEO_get_ephatt_inputs()
	Can't read ephemeris or attitude file.			MOD_PR03_FIP-3

	GEO_get_utcpole_metadata()
	Can't read the utcpole file header line			MOD_PR03_FIP-19

	GEO_read_param()					MOD_PR03_FIP-13
	Can't read line of parameter file

	GEO_read_parameter_file()				MOD_PR03_FIP-13
	Can't read header of parameter file

                        Fault Isolation Procedures

The MOD_PR03 Fault Isolation Procedures (FIPs) document specific courses
of actions to take to respond to various MOD_PR03 casualty states.
The casualty state will typically be identified by error messages MOD_PR03
generates into the LogStatus listing during processing.

Each FIP describes an action to take in response to the casualty state.
In some FIPs, a query is made as to the outcome of the action.  The
response to the query will normally lead to further actions.  Often
one FIP will reference another for these follow-up actions.


MOD_PR03_FIP-1:         Notify SDST of fault.  Retain LogStatus file and
                        L1A and Geolocation products generated when error
			message occurred for SDST investigation.

MOD_PR03_FIP-2:         Check PGE's return code:

                        Is return code 0?
                                Y:      Go to MOD_PR03_FIP-2.1
                                N:      Go to MOD_PR03_FIP-1.

MOD_PR03_FIP-2.1:       Notify SDST of event.  Error is generally non-fatal to
                        the process and the associated L1A and Geolocation
			products may be archived pending SDST review.  Retain
			granule identification and LogStatus file generated
			when error message occurred for SDST investigation.

MOD_PR03_FIP-3:         Inspect preceding LogStatus message:

                        Is it a MODIS_E_* message?
                                Y:      Perform procedure associated with that
                                        message.
                                N:      Go to MOD_PR03_FIP-1.

MOD_PR03_FIP-4:         Correct production environment and re-run PGE.

MOD_PR03_FIP-5:         Inspect file and directory permissions of the output
                        products.

                        Does process have write permissions to the output
                        products?
                                Y:      Go to MOD_PR03_FIP-7.
                                N:      Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-6:         Inspect Toolkit DEM file dates.

                        Have DEM files just been replaced?
                                Y:      Go to MOD_PR03_FIP-6.1.
                                N:      Go to MOD_PR03_FIP-1.

MOD_PR03_FIP-6.1:       Replace DEM files with previous version for PGE01
                        processing.  Report problem with DEM files to SDST.
                        Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-7:         Inspect PGE's Process control.

                        Does PGE's Process control correctly reference output
                                files?
                                Y:      Go to MOD_PR03_FIP-1.
                                N:      Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-8:         Are MOD_PR03 processing messages preceded by Error
                        messages from MOD_PR01?

                                Y:      Perform procedure associated with
                                        MOD_PR01 Error message.
                                N:      Go to MOD_PR03_FIP-8.1

MOD_PR03_FIP-8.1:       MOD_PR01 did not create L1A product.  Inspect size
                        of L1A product (identified in the Geolocation product's
                        INPUTPOINTER metadata):

                        Is the L1A product larger than 1MB?
                                Y:      Go to MOD_PR03_FIP-5.
                                N:      Go to MOD_PR03_FIP-7.

MOD_PR03_FIP-9:         Inspect Toolkit Ephemeris and attitude files.

                        Are Ephemeris and Attitude files covering the granule's
                        observation period staged?
                                Y:      Go to MOD_PR03_FIP-9.1.
                                N:      Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-9.1:       Inspect PGE's Process control.

                        Does PGE's Process control correctly reference
                                Ephemeris and Attitude files?
                                Y:      Go to MOD_PR03_FIP-1.
                                N:      Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-10:        Inspect Toolkit DEM files.

                        Are DEM files properly installed?
                                Y:      Go to MOD_PR03_FIP-10.1.
                                N:      Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-10.1:      Inspect PGE's Process control.

                        Does PGE's Process control correctly reference DEM
                                files?
                                Y:      Go to MOD_PR03_FIP-6.
                                N:      Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-11:        Inspect remaining LogStatus messages:

                        Is this message followed by an error message from the
                                PGE run?
                                Y:      Perform procedure associated with that
                                        message.
                                N:      Go to MOD_PR03_FIP-2.

MOD_PR03_FIP-12:        Inspect LUN 600280 in the PGE's Process control.

                        Is the input parameter LUN set to "SDP Toolkit"?
                                Y:      Go to MOD_PR03_FIP-11.
                                N:      Set LUN 600280 to "SDP Toolkit" and
                                        rerun PGE.  Report action to SDST,
                                        identifying granule action was taken on.

MOD_PR03_FIP-12.1:      Inspect LUN 600280 in the PGE's Process control.

                        Is the input parameter LUN set to either "TRUE" or
			"FALSE"?
                                Y:      Go to MOD_PR03_FIP-11.
                                N:      Correct LUN value.

MOD_PR03_FIP-13:        Inspect file and directory permissions of the input
                        (Geolocation LUT) file.

                        Does process have read permissions to the input file?
                                Y:      Go to MOD_PR03_FIP-4.
                                N:      Go to MOD_PR03_FIP-14.

MOD_PR03_FIP-14:        Inspect PGE's Process control.

                        Does PGE's Process control correctly reference the input
                                (Geolocation LUT) file?
                                Y:      Go to MOD_PR03_FIP-1.
                                N:      Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-15:        Inspect Toolkit planetary ephemeris file.

                        Is planetary ephemeris file properly installed?
                                Y:      Go to MOD_PR03_FIP-15.1.
                                N:      Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-15.1:      Inspect PGE's Process control.

                        Does PGE's Process control correctly reference the
                                planetary ephemeris file?
                                Y:      Go to MOD_PR03_FIP-15.2.
                                N:      Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-15.2:      Inspect Toolkit planetary ephemeris file's date.

                        Has the planetary ephemeris file just been replaced?
                                Y:      Go to MOD_PR03_FIP-15.3.
                                N:      Go to MOD_PR03_FIP-1.

MOD_PR03_FIP-15.3:      Replace the planetary ephemeris file with the previous
                        version for PGE01 processing.  Report problem with the
                        planetary ephemeris file to SDST.
                        Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-16:        Inspect Toolkit Earth Figure file.

                        Is the Earth Figure file properly installed?
                                Y:      Go to MOD_PR03_FIP-16.1.
                                N:      Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-16.1:      Inspect PGE's Process control.

                        Does PGE's Process control correctly reference the
                                Earth Figure file?
                                Y:      Go to MOD_PR03_FIP-16.2.
                                N:      Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-16.2:      Inspect Toolkit Earth Figure file's date.

                        Has the Earth Figure file just been replaced?
                                Y:      Go to MOD_PR03_FIP-16.3.
                                N:      Go to MOD_PR03_FIP-1.

MOD_PR03_FIP-16.3:      Replace the Earth Figure file with the previous
                        version for PGE01 processing.  Report problem with the
                        Earth Figure file to SDST.
                        Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-17:        Inspect previous LogStatus message from PGS_CSC_ECItoECR
                        routine:

                        Is the message from PGS_CSC_ECItoECR the following?
                                PGSTD_E_NO_LEAP_SECS
                                Y:      Go to MOD_PR03_FIP-18
                                N:      Go to MOD_PR03_FIP-17.1

MOD_PR03_FIP-17.1:      Inspect previous LogStatus message from PGS_CSC_ECItoECR
                        routine:

                        Is the message from PGS_CSC_ECItoECR the following?
                                PGSCSC_W_PREDICTED_UT1
                                PGSTD_E_NO_UT1_VALUE
                                Y:      Go to MOD_PR03_FIP-19
                                N:      Go to MOD_PR03_FIP-1.

MOD_PR03_FIP-18:        Inspect Toolkit leap second file.

                        Is the leap second file properly installed?
                                Y:      Go to MOD_PR03_FIP-18.1.
                                N:      Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-18.1:      Inspect PGE's Process control.

                        Does PGE's Process control correctly reference the
                                leap second file?
                                Y:      Go to MOD_PR03_FIP-18.2.
                                N:      Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-18.2:      Inspect Toolkit leap second file's date.

                        Has the leap second file just been replaced?
                                Y:      Go to MOD_PR03_FIP-18.3.
                                N:      Go to MOD_PR03_FIP-18.4.

MOD_PR03_FIP-18.3:      Replace the leap second file with the previous
                        version for PGE01 processing.  Report problem with the
                        leap second file to SDST.
                        Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-18.4:      Inspect Toolkit leap second file's date.

                        Is the leap second file out of date according to
                        Toolkit leap second file maintenance directives?
                                Y:      Go to MOD_PR03_FIP-18.5.
                                N:      Go to MOD_PR03_FIP-1.

MOD_PR03_FIP-18.5:      Replace the stale leap second file with an updated
                        one in accordance with Toolkit leap second file
                        maintenance directives.
                        Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-19:        Inspect Toolkit 'utcpole' file.

                        Is the 'utcpole' file properly installed?
                                Y:      Go to MOD_PR03_FIP-19.1.
                                N:      Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-19.1:      Inspect PGE's Process control.

                        Does PGE's Process control correctly reference the
                                utcpole file?
                                Y:      Go to MOD_PR03_FIP-19.2.
                                N:      Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-19.2:      Inspect Toolkit 'utcpole' file's date.

                        Has the 'utcpole' file just been replaced?
                                Y:      Go to MOD_PR03_FIP-19.3.
                                N:      Go to MOD_PR03_FIP-19.4.

MOD_PR03_FIP-19.3:      Replace the 'utcpole' file with the previous
                        version for PGE01 processing.  Report problem with the
                        leap second file to SDST.
                        Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-19.4:      Inspect Toolkit 'utcpole' file's date.

                        Is the 'utcpole' file out of date according to
                        Toolkit leap second file maintenance directives?
                                Y:      Go to MOD_PR03_FIP-19.5.
                                N:      Go to MOD_PR03_FIP-1.

MOD_PR03_FIP-19.5:      Replace the stale 'utcpole' file with an updated
                        one in accordance with Toolkit leap second file
                        maintenance directives.
                        Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-20:        Inspect previous LogStatus message from
                        PGS_EPH_EphemAttit routine:

                        Is the message from PGS_EPH_EphemAttit the following?
                                PGSTD_E_NO_LEAP_SECS
                                Y:      Go to MOD_PR03_FIP-18
                                N:      Go to MOD_PR03_FIP-20.1

MOD_PR03_FIP-20.1:      Inspect previous LogStatus message from
                        PGS_EPH_EphemAttit routine:

                        Is the message from PGS_EPH_EphemAttit the following?
                                PGSEPH_E_NO_SC_EPHEM_FILE
                                Y:      Go to MOD_PR03_FIP-9
                                N:      Go to MOD_PR03_FIP-1.

MOD_PR03_FIP-21:        Inspect Toolkit land/sea mask files.

                        Are land/sea mask files properly installed?
                                Y:      Go to MOD_PR03_FIP-21.1.
                                N:      Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-21.1:      Inspect PGE's Process control.

                        Does PGE's Process control correctly reference the
                                land/sea mask files?
                                Y:      Go to MOD_PR03_FIP-1.
                                N:      Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-22:        Inspect L1A and Geolocation Metadata Control Files(MCF):

                        Are MCF files properly staged? (Are they in the
                        correct location?  Does the PGE have read access
                        to them?)
                                Y:      Go to MOD_PR03_FIP-22.1.
                                N:      Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-22.1:      Inspect PGE's Process control.

                        Does PGE's Process control correctly reference the
                                L1A and Geolocation MCF files?
                                Y:      Go to MOD_PR03_FIP-1.
                                N:      Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-23:        Inspect LUN 600310 in the PGE's Process control.

                        Is the input parameter LUN set to "TRUE"?
                                Y:      Go to MOD_PR03_FIP-1.
                                N:      Set LUN 600310 to "TRUE" and
                                        rerun PGE.  Report action to SDST,
                                        identifying granule action was taken on.

MOD_PR03_FIP-24:        Inspect Geolocation LUT file date.

                        Has Geolocation LUT file just been replaced?
                                Y:      Go to MOD_PR03_FIP-24.1.
                                N:      Go to MOD_PR03_FIP-1.

MOD_PR03_FIP-24.1:      Make sure that Geolocation LUT file's Revision
			keyword matches the corresponding entry in the PCF.
                        Go to MOD_PR03_FIP-4.

MOD_PR03_FIP-25		Look at the Satellite Instrument Configuration parameter
			in the PCF.
				"AM1M":	Go to MOD_PR03_FIP-25.1
				"PM1M": Go to MOD_PR03_FIP-25.2
				Other: Replace the PCF, which is corrupt.

MOD_PR03_FIP-25.1	Examine the spacecraft_ID parameter in the first
			Geolocation LUT listed in the PCF. Does it say "Terra"?
				Yes: Go to MOD_PR03_FIP-1
				No: Replace with the Terra version of the LUT.

MOD_PR03_FIP-25.2	Examine the spacecraft_ID parameter in the second
			Geolocation LUT listed in the PCF. Does it say "Aqua"?
				Yes: Go to MOD_PR03_FIP-1
				No: Replace with the Aqua version of the LUT.

MOD_PR03_FIP-26		Look at the Satellite Instrument Configuration parameter
			in the PCF.
				"AM1M":	Go to MOD_PR03_FIP-26.1
				"PM1M": Go to MOD_PR03_FIP-26.2
				Other: Replace the PCF, which is corrupt.

MOD_PR03_FIP-26.1	Look at AssociatedPlatformShortname.1 in the
			corresponding MOD01 file. Does it say "Terra"?
				Yes: Go to MOD_PR03_FIP-1
				No: Change either the PCF to match the input

MOD_PR03_FIP-26.2	Look at AssociatedPlatformShortname.1 in the
			corresponding MOD01 file. Does it say "Aqua"?
				Yes: Go to MOD_PR03_FIP-1
				No: Change either the PCF to match the input
				    file, or vice versa.

MOD_PR03_FIP-27:        Inspect the PGE's Process Control File, looking for the
			specified LUN.

                        Is there a correctly formatted entry for that LUN?
                                Y:      Go to MOD_PR03_FIP-1.
                                N:      Set LUN 600310 to "TRUE" and
                                        rerun PGE.  Report action to SDST,
                                        identifying granule action was taken on.

MOD_PR03_FIP-28:        Inspect the 800500 LUN in the PGE's Process control file

                        Is there a correctly formatted entry for that LUN?
                                Y:      Go to MOD_PR03_FIP-1.
                                N:      Set LUN 600310 to "TRUE" and
                                        rerun PGE.  Report action to SDST,
                                        identifying granule action was taken on.

MOD_PR03_FIP-29:	Increase the amount of memory available to the process,
			and re-run it.
MOD_PR03_FIP-30		Look at the Satellite Instrument Configuration parameter
			in the PCF.
				"AM1M":	Go to MOD_PR03_FIP-30
				"PM1M": Go to MOD_PR03_FIP-30
				Other: Replace the PCF, which is corrupt.

MOD_PR03_FIP-30.1:	Examine the spacecraft parameter in the maneuver list
			Does it say "Terra"?
				Yes: Go to MOD_PR03_FIP-1
				No: Replace with the Terra version of the LUT.

MOD_PR03_FIP-30.2:	Examine the spacecraft parameter in the maneuver list
			Does it say "Aqua"?
				Yes: Go to MOD_PR03_FIP-1
				No: Replace with the Aqua version of the LUT.

MOD_PR03_FIP-31:        Inspect the 600003 LUN in the PGE's Process control

                        Is there a correctly formatted entry for that LUN?
                                Y:      Go to MOD_PR03_FIP-1.
                                N:      Set LUN 600003 to the correct value and
                                        rerun PGE.  Report action to SDST,
                                        identifying granule action was taken on.

MOD_PR03_FIP-32:        Inspect the 600002 LUN in the PGE's Process control

                        Is there a correctly formatted entry for that LUN?
                                Y:      Go to MOD_PR03_FIP-1.
                                N:      Set LUN 600002 to the correct value and
                                        rerun PGE.  Report action to SDST,
                                        identifying granule action was taken on.


19. PLATFORM NOTES.
See 3f for platform dependencies.

