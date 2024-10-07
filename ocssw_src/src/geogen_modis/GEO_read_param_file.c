/*
 * $Log: GEO_read_param_file.c,v $
 * Revision 6.4  2010/04/23 14:48:52  kuyper
 * Bug 2472: Converted new numbers controlling setting of quality flags into
 *   parameters retrieved from the parameter file.
 *
 * Revision 6.3  2009/09/15 20:11:30  kuyper
 * Corrected #include list.
 *
 * James Kuyper		James.R.Kuyper@nasa.gov
 *
 * Revision 6.2  2009/06/12 20:34:50  kuyper
 * Changed to allow modification to the characteristics of band 0, without
 *   inappropriate consequences for the calculation of u_tel.
 *
 * Revision 6.1  2009/05/26 18:35:21  ltan
 * Added hires_scale.
 *
 * Revision 5.3  2005/10/11 14:55:15  vlin
 * Code updated to read in orbit descend time, orbit tolerance, orbit period, and
 * transition orbit number.  vlin@saicmodis.com
 *
 * Revision 5.2  2005/03/24 17:31:52  kuyper
 * Corrected handling of macros with empty arguments.
 *
 * Revision 5.1  2004/08/03 15:40:17  vlin
 * Add new spacecraft orientation correction matrix dependent
 * in a piece-wise linear manner on solar elevation angle.
 *
 * Revision 4.8  2003/09/09 19:34:29  kuyper
 * Corrected handling of rpy_count==0 case.
 *
 * Revision 4.7  2003/08/22 22:19:22  kuyper
 * Changed T_inst2SD to T_sc2SD.
 *
 * Revision 4.6  2003/04/02 21:41:37  vlin
 * Updated after code walk through
 *
 * Revision 4.5  2003/01/29 21:14:45  vlin
 * Updated to get the expected LUT RCS revision value from PCF.
 * vlin@saicmodis.com
 *
 * Revision 4.4  2002/06/06 14:27:12  vlin
 * Read the value of "rpy_count" from the parameter file.
 *
 * Revision 4.3  2002/06/04 14:04:39  vlin
 * Updated after code walkthrough.
 *
 * Revision 4.2  2002/06/03 18:32:43  vlin
 * Allow return value "PGSTD_E_NO_LEAP_SECS" from function PGS_TD_UTCtoTAI().
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>
#include <stddef.h>
#include "PGS_TD.h"
#include "GEO_main.h"
#include "GEO_main_func.h"


typedef enum{integral, floating, revise, string} datatype_enum;

static PGSt_SMF_status GEO_read_param(
       FILE          * const param_file,
       const char      name[],
       datatype_enum   type,
       size_t	       size,
       int             elements,
       char            data[])

{
/***************************************************************************
!C 
  
!Description: 
        Routine for reading in a single parameter from the parameter file.
  
!Input Parameters:
  	param_file		Points at the file handle for parameter file.
  	name			Name of the parameter.
  	type			The category of the data type.
  	size			The size of a single element, in bytes.
        elements                The number of elements in the parameter
  
!Output Parameters:
  	data			char pointer to the place where the parameter
  				is to be stored.
  
Return Values:
     MODIS_E_BAD_SCAN	        If a call to fscanf() failed.
     MODIS_E_PREMATURE_EOF	If end of file reached unexpectedly.
     MODIS_E_WRONG_FIELD	If a field doesn't have the expected name.
     PGS_E_UNIX		        If a standard library function failed.
     PGS_S_SUCCESS		Otherwise
  
Externally Defined:
     errno                      errno.h
     MODIS_E_BAD_SCAN           PGS_MODIS_35251.h
     MODIS_E_PREMATURE_EOF      PGS_MODIS_35251.h
     MODIS_E_WRONG_FIELD        PGS_MODIS_35251.h 
     PGS_E_UNIX                 PGS_SMF.h
     PGS_S_SUCCESS              PGS_SMF.h 
  
Called by:
     GEO_read_param_file()
  
Routines Called:
     PGS_SMF_SetUnixMsg      Reports message based upon 'errno'.
     modsmf                  Logs status messages.
  
!Revision History:
     Please see the top of the file.

Requirements:
     PR03-F-1-1 
     PR03-F-1-2 
     PR03-F-1-3 
     PR03-F-1-4 
     PR03-F-1-5 
     PR03-F-1-6

!Team-unique Header:

     This software is developed by the MODIS Science Data Support
     Team for the National Aeronautics and Space Administration,
     Goddard Space Flight Center, under contract NAS5-32373.
  
Design Notes

!END
*****************************************************************************/

/* No need to validate arguments, since caller, which sets their value, is in
 * the same file.
 */

  int i, next_char=0;
  int ret_val;	      		/* value returned to this function	*/
  const char *scanstr="";
  char buffer[128] = "";
  char msgbuf[128] = "";
  char filefunc[] = __FILE__ ", GEO_read_param";

  while ((next_char=fgetc(param_file))=='#') {  /* Skip comment lines. */
     if (fgets(buffer, sizeof(buffer), param_file)!=buffer) {
	if (ferror(param_file)) {
	   PGS_SMF_SetUNIXMsg(errno, "fgets()", filefunc);
           return PGS_E_UNIX;
        }
	else {
           sprintf(msgbuf, "parameter file, for parameter %s \n", name);
	   modsmf(MODIS_E_PREMATURE_EOF, msgbuf, filefunc);
           return MODIS_E_PREMATURE_EOF;
        }
     }
  }

  if (next_char==EOF) { /* Failed to read characters. */
     if (ferror(param_file)) {
	PGS_SMF_SetUNIXMsg(errno, "fgetc()", filefunc);
        return PGS_E_UNIX;
     }
     else {
	modsmf(MODIS_E_PREMATURE_EOF, "parameter file.", filefunc);
        return MODIS_E_PREMATURE_EOF;
     }
  }
  ungetc(next_char, param_file);   /* Push character read back into input buffer */

  ret_val = fscanf(param_file,"%32s = " , buffer);   /* Read in parameter name. */
  if (ret_val!=1) {
     if (ferror(param_file))
	 PGS_SMF_SetUNIXMsg(errno, "fscanf()", filefunc);
     sprintf(msgbuf, "(\"%.20s\") returned %d", name, ret_val);
     modsmf(MODIS_E_BAD_SCAN, msgbuf, filefunc);
     return MODIS_E_BAD_SCAN;
  }
      
  if (strcmp(name, buffer)!=0) {/* Compare with expected name. */
     sprintf(msgbuf, "%.20s, actual:%.20s", name, buffer);
     modsmf(MODIS_E_WRONG_FIELD, msgbuf, filefunc);
     return MODIS_E_WRONG_FIELD;
  }

    /* NOTE: use of sizeof() below allows correct handling of typedefed
      parameters, such as int32, or PGSt_double. */;

  if (type == revise) {
      /* the 31 below must be < sizeof(revision). */
      scanstr = "%*[$]Revision: %31s $ ";
      /* The %*[] serves solely to prevent interpretation of this scanstr
	value as an RCS keyword. The scanstr is supposed to match an RCS
	Revision keyword string; it's not supposed to be one itself. */
  }
  else if (type == string) {
      scanstr = "%s\n";
  }
  else if (type == integral) {
      if(size == sizeof(int))
	scanstr = "%d ";
      else if(size == sizeof(long))
	scanstr = "%ld ";
      else if(size == sizeof(short))
	scanstr = "%hd ";
      else
	scanstr = "invalid integral size";
  }
  else if (type == floating) {
      if(size == sizeof(double))
	scanstr = "%lf ";
      else if(size == sizeof(float))
	scanstr = "%f ";
      else
	scanstr = "invalid floating point size";
  }
  else
      scanstr = "invalid type";

  /* Note: the 'invalid' scanstr values above will cause fscanf() to fail
     below, and the resulting error message will contain 'scanstr'. This
     should only occur if initialization of members[] is invalid; never
     during production runs. */
    
  /* Scan in each element. */
  for (i=0; i<elements; i++) {
      ret_val = fscanf(param_file, scanstr, data);
      if (ret_val != 1) {
	if (ferror(param_file))
	   PGS_SMF_SetUNIXMsg(errno, "fscanf()", filefunc);
	else if (feof(param_file))
	   modsmf(MODIS_E_PREMATURE_EOF, "parameter file. while scanning element.", 
                  filefunc);
	sprintf(msgbuf, "(param_file, \"%.16s\",&%.20s[%d]) returned %d",
	        scanstr, name, i, ret_val);
	modsmf(MODIS_E_BAD_SCAN, msgbuf, filefunc);
	return MODIS_E_BAD_SCAN;
      }
      data += size;
  }	/* For each element	*/

  return PGS_S_SUCCESS;
}

/*============================================================================*/

PGSt_SMF_status GEO_read_param_file(GEO_param_struct *param)
{

/******************************************************************************
!C

!Description:   Subroutine in Main group of the Level-1A geolocation software
    to open, read, and close the parameter file, and to compute the telescope
    view vector from the other parameters.

!Input Parameters: None

!Output Parameters: GEO_param_struct *param

Return parameters:
    MODIS_E_BAD_INPUT_ARG       If any argument was detected as invalid.
    MODIS_E_GEO                 If any subroutine (other than PGS_IO_GenClose)
                                fails.
    MODIS_E_GEO_WRONG_PLATFORM  If the platform from the PCF doesn't match the
                                parameters file.
    MODIS_E_PREMATURE_EOF       Reached the end of the parameter file without
                                finding all expected contents.
    MODIS_E_WRONG_FILE          The PCF points to a file which isn't a valid
                                Geolocation parameters file.
    MODIS_E_WRONG_LUT           The incorrect version of the parameters file
                                was used.
    PGS_E_UNIX                  An I/O error occurred, for which the value of
                                errno may be a useful diagnostic.
    PGS_S_SUCCESS               Otherwise

Externally Defined:
    errno                       errno.h
    BASE_SAMPLES                GEO_geo.h
    MAX_BAND_NUMBER             GEO_geo.h
    MAX_SAMPLES                 GEO_geo.h
    MODIS_E_BAD_INPUT_ARG       PGS_MODIS_35251.h
    MODIS_E_GEO                 PGS_MODIS_35251.h
    MODIS_E_GEO_WRONG_PLATFORM  PGS_MODIS_35251.h
    MODIS_E_PREMATURE_EOF       PGS_MODIS_35251.h
    MODIS_E_WRONG_FILE          PGS_MODIS_35251.h
    MODIS_E_WRONG_LUT           PGS_MODIS_35251.h
    PARAM                       GEO_main.h
    PARAM_TAG                   GEO_parameters.h
    PGS_E_UNIX                  PGS_SMF.h
    PGS_S_SUCCESS               PGS_SMF.h
    PGSd_IO_Gen_Read            PGS_IO_Gen.h
    PGSd_PC_LINE_LENGTH_MAX     PGS_PC.h
    TIMECODEASIZE               smfio.h

Called by:
    main()
         
Routines Called:
    GEO_read_param - reads in one parameter (which may be an array).
    modsmf - writes error status messages to log
    PGS_IO_Gen_Close - closes an ascii file
    PGS_IO_Gen_Open - opens a general file
    PGS_SMF_SetUNIXMsg - writes UNIX error status messages to log.
    PGS_SMF_GetConfigData - reads configuration values from PCF

!Revision History:
    Please see the top of the file.
    
Requirements:
    PR03-F-1-1 
    PR03-F-1-2 
    PR03-F-1-3 
    PR03-F-1-4 
    PR03-F-1-5 
    PR03-F-1-6

!Team-unique Header:

    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

!END
****************************************************************************/

#define LUT_REVISION_LUN 600021

#define SCALAR_PARM(sub,d,name,type)  \
 GEO_read_param(fp, #name, type, sizeof(param->sub d name), 1, \
 (char *)&param->sub d name)

#define VECTOR_PARM(sub, d, name,type)  \
 GEO_read_param(fp, #name, type, sizeof(param->sub d name[0]), \
 (int)(sizeof(param->sub d name)/sizeof(param->sub d name[0])), \
 (char *)&param->sub d name[0])

#define ARRAY_PARM(sub, d, name,type)   \
 GEO_read_param(fp, #name, type, sizeof(param->sub d name[0][0]), \
 (int)(sizeof(param->sub d name)/sizeof(param->sub d name[0][0])), \
 (char *)&param->sub d name[0][0])

/* Local variable declarations */

  char satellite_instrument[PGSd_PC_LINE_LENGTH_MAX]; 
                      /* satellite instrument: AM1M for terra, PM1M for aqua */
  char expected_revision[PGSd_PC_VALUE_LENGTH_MAX];
  PGSt_IO_Gen_FileHandle *fp;
  PGSt_PC_Logical SATELLITE_INSTRUMENT_LUN = 800510; 
                      /* Logical Unit Number for satellite instrument */
  PGSt_integer version;         /* parameter file version
				   version = 1 for Terra satellite
				           = 2 for Aqua  satellite */
  double center_detector;
  int i = 0;    		/* iteration parameter */ 
  int det;			/* detector number	*/
  int band;
  PGSt_SMF_status status, ret_code = PGS_S_SUCCESS;
  char (*asciiUTC)[TIMECODEASIZE] = NULL;
  char msgbuf[64] = "";
  char buffer[128] = "";
  char filefunc[] = __FILE__ ", GEO_read_param_file";

/* Take configuration parameter value (AM1M/Terra or PM1M/Aqua) from pcf file */
   status = PGS_PC_GetConfigData(SATELLITE_INSTRUMENT_LUN,
       satellite_instrument);
   if (status != PGS_S_SUCCESS) {
     sprintf(msgbuf, "PGS_PC_GetConfigData()");
     modsmf(MODIS_E_GEO, msgbuf, filefunc);
     return MODIS_E_GEO;
   }

   if (strcmp(satellite_instrument, "AM1M") == 0)
     version = 1;
   else if (strcmp(satellite_instrument, "PM1M") == 0)
     version = 2;
   else {
     sprintf(msgbuf, " SatelliteInstrument is %s", satellite_instrument);
     modsmf(MODIS_E_GEO_WRONG_PLATFORM, msgbuf, filefunc);
     return MODIS_E_GEO_WRONG_PLATFORM;
   }
   
/* Open parameter file for reading */
   status = PGS_IO_Gen_Open(PARAM, PGSd_IO_Gen_Read, &fp, version);
   if (status != PGS_S_SUCCESS) {
     sprintf(msgbuf, "PGS_IO_Gen_Open(%ld) ", (long)PARAM);
     modsmf(MODIS_E_GEO, msgbuf, filefunc);
     return MODIS_E_GEO;
   }

/* Read the first line of parameter file */ 
  if (fgets(buffer, sizeof(buffer), fp)!=buffer) {
    if (ferror(fp)) {
        PGS_SMF_SetUNIXMsg(errno, "fgets()", filefunc);
        ret_code = PGS_E_UNIX;
    }
    else {
        modsmf(MODIS_E_PREMATURE_EOF, "parameter file.", filefunc);
        ret_code = MODIS_E_PREMATURE_EOF;
    }
  }
  else {   /* Skip leading white space */
    for(i=0; buffer[i]!='\0' && isspace((int)buffer[i]); i++);

    if(strncmp(PARAM_TAG, buffer+i, sizeof(PARAM_TAG)-1)!=0)
    {	/* First line of parameter file doesn't contain required PARAM_TAG */
      modsmf(MODIS_E_WRONG_FILE, buffer+i, filefunc);
      ret_code = MODIS_E_WRONG_FILE;
    }
  }

/* Read each geolocation parameter */
  if ((SCALAR_PARM(, , revision, revise) != PGS_S_SUCCESS) || 
      (SCALAR_PARM(, , spacecraft_ID, string) != PGS_S_SUCCESS) ||
      (SCALAR_PARM(geometry_params, ., band_number, integral)!=PGS_S_SUCCESS) ||
      (VECTOR_PARM(geometry_params, ., N_samp, integral) != PGS_S_SUCCESS) ||
      (VECTOR_PARM(geometry_params,.,band_position, floating)!=PGS_S_SUCCESS) ||
      (ARRAY_PARM(geometry_params,.,det_position, floating) != PGS_S_SUCCESS) ||
      (VECTOR_PARM(geometry_params,.,det_space, floating) != PGS_S_SUCCESS) ||
      (VECTOR_PARM(geometry_params,.,focal_length, floating)!=PGS_S_SUCCESS) ||
      (SCALAR_PARM(mirror_prep_params,., N_reset, floating) != PGS_S_SUCCESS) ||
      (SCALAR_PARM(geometry_params,., t_reset, floating) != PGS_S_SUCCESS) ||
      (SCALAR_PARM(geometry_params,., t_frame, floating) != PGS_S_SUCCESS) ||
      (SCALAR_PARM(mirror_prep_params,.,t_vernier, floating)!=PGS_S_SUCCESS) ||
      (SCALAR_PARM(mirror_prep_params,.,t_encoder, floating)!=PGS_S_SUCCESS) ||
      (SCALAR_PARM(mirror_prep_params,.,sample_impulse,integral)!=PGS_S_SUCCESS) ||
      (VECTOR_PARM(geometry_params,.,F_offset, floating) != PGS_S_SUCCESS) ||
      (VECTOR_PARM(mirror_model,.,mirr_side1_range,floating)!=PGS_S_SUCCESS) ||
      (SCALAR_PARM(, , poly_degree, integral) != PGS_S_SUCCESS) ||
      (ARRAY_PARM(, , poly_coef, floating) != PGS_S_SUCCESS) ||
      (SCALAR_PARM(mirror_model, ., alpha, floating) != PGS_S_SUCCESS) ||
      (SCALAR_PARM(mirror_model, ., beta, floating) != PGS_S_SUCCESS) ||
      (SCALAR_PARM(mirror_model, ., gammaa, floating) != PGS_S_SUCCESS) ||
      (ARRAY_PARM(coord_trans, ., T_tel2inst, floating) != PGS_S_SUCCESS) ||
      (ARRAY_PARM(coord_trans, ., T_mirr2inst, floating) != PGS_S_SUCCESS) ||
      (ARRAY_PARM(coord_trans, ., T_inst2sc, floating) != PGS_S_SUCCESS) ||
      (ARRAY_PARM(coord_trans, ., T_sc2SD, floating) != PGS_S_SUCCESS) ||
      (SCALAR_PARM(, , max_extrap, floating) != PGS_S_SUCCESS) ||
      (SCALAR_PARM(ancil_params.ancil_scale_factors, ., S_angvel, floating)
	  !=PGS_S_SUCCESS)||
      (SCALAR_PARM(ancil_params.ancil_scale_factors, ., S_attitude, floating)
	  !=PGS_S_SUCCESS)||
      (SCALAR_PARM(ancil_params.ancil_scale_factors, .,S_position, floating)
	  !=PGS_S_SUCCESS) ||
      (SCALAR_PARM(ancil_params.ancil_scale_factors, .,S_velocity, floating)
	  !=PGS_S_SUCCESS) ||
      (SCALAR_PARM(, , angle_scale, floating) != PGS_S_SUCCESS) ||
      (SCALAR_PARM(, , hires_scale, floating) != PGS_S_SUCCESS) ||
      (SCALAR_PARM(, , range_scale, floating) != PGS_S_SUCCESS) ||
      (VECTOR_PARM(ancil_params.ancil_words, ., attit_words, integral)
	  != PGS_S_SUCCESS) ||
      (VECTOR_PARM(ancil_params.ancil_words, ., orbit_words, integral)
	  != PGS_S_SUCCESS) ||
      (VECTOR_PARM(ancil_params.ancil_words, ., time_words, integral)
	  != PGS_S_SUCCESS) ||
      (VECTOR_PARM(mirror_prep_params, ., sector_word, integral)
	  != PGS_S_SUCCESS) ||
      (VECTOR_PARM(mirror_prep_params, ., vernier_word, integral)
	  != PGS_S_SUCCESS) ||
      (SCALAR_PARM(, , max_non_gap, floating) != SUCCESS) ||
      (VECTOR_PARM(orbit_valid_params, ., ang_mom_limit, floating)
	  != PGS_S_SUCCESS) ||
      (VECTOR_PARM(orbit_valid_params, ., ang_mom_z_limit, floating)
	  != PGS_S_SUCCESS) ||
      (SCALAR_PARM(orbit_valid_params, ., orbit_consistency, floating)
	  != PGS_S_SUCCESS) ||
      (VECTOR_PARM(orbit_valid_params, ., descend_time_0, floating)
	  !=PGS_S_SUCCESS) ||
      (VECTOR_PARM(orbit_valid_params, ., orbit_tolerance, floating)
	  !=PGS_S_SUCCESS) ||
      (VECTOR_PARM(orbit_valid_params, ., period, floating) !=PGS_S_SUCCESS) ||
      (SCALAR_PARM(orbit_valid_params, ., transition_orbit, integral)
          !=PGS_S_SUCCESS) ||
      (VECTOR_PARM(orbit_valid_params, ., position_abs_limit, floating)
	  !=PGS_S_SUCCESS) ||
      (VECTOR_PARM(orbit_valid_params, ., position_mag_limit, floating)
	  !=PGS_S_SUCCESS) ||
      (VECTOR_PARM(orbit_valid_params, ., velocity_abs_limit, floating)
	  !=PGS_S_SUCCESS) ||
      (VECTOR_PARM(orbit_valid_params, ., velocity_mag_limit, floating)
	  !=PGS_S_SUCCESS) ||
      (SCALAR_PARM(orbit_valid_params, ., eph_max_short_gap, floating)
	  !=PGS_S_SUCCESS) ||
      (VECTOR_PARM(attit_valid_params, ., angvel_abs_limit, floating)
	  != PGS_S_SUCCESS) ||
      (SCALAR_PARM(attit_valid_params, ., angvel_del_limit, floating)
	  != PGS_S_SUCCESS) ||
      (SCALAR_PARM(attit_valid_params, ., attit_consistency, floating)
	  !=PGS_S_SUCCESS) ||
      (VECTOR_PARM(attit_valid_params, ., attitude_abs_limit, floating)
	  !=PGS_S_SUCCESS) ||
      (SCALAR_PARM(attit_valid_params, ., attitude_del_limit, floating)
	  !=PGS_S_SUCCESS) ||
      (VECTOR_PARM(attit_valid_params, ., att_valid_range, integral)
	  !=PGS_S_SUCCESS) ||
      (SCALAR_PARM(attit_valid_params, ., att_max_short_gap, floating)
	  !=PGS_S_SUCCESS) ||
      (VECTOR_PARM(mirror_prep_params, ., mirr_abs_limit, floating)
	  != PGS_S_SUCCESS) ||
      (SCALAR_PARM(mirror_prep_params, ., mirr_del_limit, floating)
	  != PGS_S_SUCCESS) ||
      (VECTOR_PARM(mirror_prep_params, ., encoder_gap, floating)
	  != PGS_S_SUCCESS) ||
      (SCALAR_PARM(mirror_prep_params, ., encoder_adjustment, floating)
	  != PGS_S_SUCCESS) ||
      (SCALAR_PARM(mirror_prep_params, ., packet_interval, floating)
	  != PGS_S_SUCCESS) ||
      (SCALAR_PARM(, , RMS_error, floating) != PGS_S_SUCCESS) ||
      (ARRAY_PARM(, , temp_coeff, floating) != PGS_S_SUCCESS) || 
      (ARRAY_PARM(, , sol_elev_cor, floating) != PGS_S_SUCCESS) || 
      (SCALAR_PARM(coord_trans, ., rpy_count, integral) != PGS_S_SUCCESS) )
  {
      modsmf(MODIS_E_GEO, "GEO_read_param()", filefunc);
      ret_code = MODIS_E_GEO;
  }
  else {
      if (PGS_PC_GetConfigData(LUT_REVISION_LUN, expected_revision) != PGS_S_SUCCESS) {
          sprintf(msgbuf, "PGS_PC_GetConfig_data(%d)", LUT_REVISION_LUN);
          modsmf(MODIS_E_GEO, msgbuf, filefunc);  
          ret_code = MODIS_E_GEO;
      }
      else if (strncmp(param->revision, expected_revision, sizeof(expected_revision))){
          sprintf(msgbuf, "%s, expected: %s", param->revision, expected_revision);
          modsmf(MODIS_E_WRONG_LUT, msgbuf, filefunc);
          ret_code = MODIS_E_WRONG_LUT; 
      }

      if(param->coord_trans.rpy_count)
      {
	  asciiUTC = (char(*)[TIMECODEASIZE])
	      malloc(param->coord_trans.rpy_count*TIMECODEASIZE);
	  param->coord_trans.rpy_times = (double *)
	      malloc(param->coord_trans.rpy_count*sizeof(double));
	  param->coord_trans.rpy_inst2sc = (double(*)[3])
	      malloc(param->coord_trans.rpy_count*sizeof(double)*3);
	  if (asciiUTC == NULL || param->coord_trans.rpy_times == NULL ||
	      param->coord_trans.rpy_inst2sc == NULL)
	  {
	      modsmf(MODIS_E_GEO, "malloc()", filefunc);
	      ret_code = MODIS_E_GEO;
	  }
	  else {
	      if (GEO_read_param(fp, "rpy_times", string, TIMECODEASIZE,
		  param->coord_trans.rpy_count, (char *)asciiUTC)
		  != PGS_S_SUCCESS ||
		  GEO_read_param(fp, "rpy_inst2sc",floating, sizeof(double),
		  param->coord_trans.rpy_count*3,
		  (char *)param->coord_trans.rpy_inst2sc) != PGS_S_SUCCESS ||
		  GEO_read_param(fp, "end_of_parameters", integral, (size_t)0, 0,
		  0) != PGS_S_SUCCESS)
	      {
		  modsmf(MODIS_E_GEO, "GEO_read_param(rpy params)", filefunc);
		  ret_code = MODIS_E_GEO;
	      }

	      for (i = 0; i < param->coord_trans.rpy_count; i++)
	      {
		  status = PGS_TD_UTCtoTAI(asciiUTC[i],
		      &param->coord_trans.rpy_times[i]);
		  if (status != PGS_S_SUCCESS && status != PGSTD_E_NO_LEAP_SECS)
		  {
		      sprintf(msgbuf, "PGS_TD_UTCtoTAI(%s)", asciiUTC[i]);
		      modsmf(MODIS_E_GEO, msgbuf, filefunc);
		      ret_code = MODIS_E_GEO;
		  }
	     }

	     free(asciiUTC);
	  }
      }
  }

  if (PGS_IO_Gen_Close(fp)!=PGS_S_SUCCESS)
      modsmf(MODIS_E_GEO, "PGS_IO_GenClose()", filefunc);

  /* Delayed exit for errors detected while file was open. */
  if (ret_code != PGS_S_SUCCESS)
      return ret_code;
  
  /* tested the matches between satellite_instruments and spacecraft_IDs.
     the right couplings are "AM1M/Terra" and "PM1M/Aqua"       */
  if ( (strcmp(satellite_instrument, "AM1M") == 0 &&
        strcmp(param->spacecraft_ID, "Terra") != 0) ||
       (strcmp(satellite_instrument, "PM1M") == 0 &&
	strcmp(param->spacecraft_ID, "Aqua") != 0))
  {
     sprintf(msgbuf, ": parameter file is for %s", param->spacecraft_ID);
     modsmf(MODIS_E_GEO_WRONG_PLATFORM, msgbuf, filefunc);
     return MODIS_E_GEO_WRONG_PLATFORM;
  }

  /* Validity checks of parameters needed for calculation of telescope view
    vector.	*/
  band = param->geometry_params.band_number;
  if (band <0 || band > MAX_BAND_NUMBER) {
      sprintf(msgbuf, "band_number = %d", band);
      modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
      return MODIS_E_BAD_INPUT_ARG;
  }

  if (param->geometry_params.N_samp[band] < 1 ||
      param->geometry_params.N_samp[band] > MAX_SAMPLES) {
      sprintf(msgbuf, "N_samp[%d] = %hd", band,
      (short)param->geometry_params.N_samp[band]);
      modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
      return MODIS_E_BAD_INPUT_ARG;
  }

  if (param->geometry_params.focal_length[0]<=0.0) {
      sprintf(msgbuf, "focal_length[0]=%f",
      param->geometry_params.focal_length[0]);
      modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
      return MODIS_E_BAD_INPUT_ARG;
  }

  param->num_detectors = BASE_SAMPLES * param->geometry_params.N_samp[band];
  center_detector = 0.5*(double)(param->num_detectors-1);

  /* Calculate telescope view vector	*/
  for (det=0; det<param->num_detectors; det++)  
  { 
    param->u_tel[det][0] = param->geometry_params.det_position[band][0]
      -(param->geometry_params.band_position[band]
      - 0.5/(double)param->geometry_params.N_samp[band])*
      param->geometry_params.focal_length[band] * 0.00054/0.380859;
    /* The band_position is offset by 0.5, because it gives the trailing edge
      of the detector, rather than the center.  */
 
    param->u_tel[det][1] = param->geometry_params.det_position[band][1]
      + (center_detector-(double)det) * param->geometry_params.det_space[band];
    /* We add (center_detector-det), instead of subtracting it, to implement the
      detector renumbering required by CCR-184. The purpose is to have
      increasing detector numbers in the positive track direction.      */
 
    param->u_tel[det][2] = param->geometry_params.focal_length[band];
  }
 
  return PGS_S_SUCCESS;
}
