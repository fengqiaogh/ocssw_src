# Start of 35251.t
#*************************************************************************
# 35251.t 
#
#!Description: MODIS Geolocation SMF file
#
#!Input Parameters: N/A  
#
#!Output Parameters: N/A
#
#!Revision History:
#MODIS_35251.t,v
#Revision 6.2  2010/03/31 19:54:16  kuyper
#Modified EXCESS_FILES message to allow expected number to vary.
#
#Revision 6.1  2009/05/30 18:47:39  kuyper
#Removed obsolete error messages.
#
#Revision 5.5  2006/10/04 17:40:30  kuyper
#Corrected zero message.
#
#Revision 5.4  2006/01/27 19:08:00  kuyper
#Reinstated RCS Log keyword.
#
# MODIS_35251.t,v
# Revision 5.1  2004/07/16 22:55:17  kuyper
# Added MODIS_N_GEO_MANEUVER.
#
# Revision 4.12  2004/04/09 15:59:37  vlin
# MODIS_E_GEO_GET_POLY_COEF & MODIS_E_GEO_POLY_FIT removed.
#
# Revision 4.11  2004/03/16 00:00:42  kuyper
# Added MODIS_E_GEO_MIRR_MOTION.
#
# Revision 4.10  2003/11/14 17:45:21  vlin
# corrected a typo for MODIS_W_GEO_MISSING_MET
#
# Revision 4.9  2003/08/20 16:04:17  kuyper
# Corrected start-up message.
#
# Revision 4.8  2003/05/14 18:15:38  kuyper
# Changed EXCESS_FILES from _E_ to _N_.
# Added MISSING_MET message, dropped obsolete ones.
#
# Revision 4.7  2003/03/03 18:44:20  kuyper
# Corrected typos in messages.
#
# Revision 4.6  2003/02/20 18:09:23  vlin
# MODIS_E_GEO_MISSING_INPUTS added
#
#!Team-unique Header:
#
#  This software is developed by the MODIS Science Data Support
#  Team for the National Aeronautics and Space Administration,
#  Goddard Space Flight Center, under contract NAS5-32373.
#
#!END
#**************************************************************************
#
%INSTR = MODIS
%LABEL = MODIS
%SEED = 35251

#Dynamic timestamp message
MODIS_E_BAD_DEM %sInvalid input argument: %s
MODIS_E_BAD_INPUT_ARG %sInvalid input argument %s
MODIS_E_BAD_SCAN %sfscanf %s
MODIS_E_BAD_VIEW_VEC %sUnable to compute sample view vector %s
MODIS_E_DATA_SCAN  %sNo valid encoder data in scan %s
MODIS_E_DEM_DATA_SIZE %sUnexpected data size: %s
MODIS_E_DEM_IN_MEM %sDEM in memory does not include lat/lon %s 
MODIS_E_DEM_METADATA %sUnexpected DEM Metadata - units: %s
MODIS_E_GEO %sError returned by function %s
MODIS_E_GEO_ANCIL_DATA  %sFewer than two unflagged samples: %s
MODIS_E_GEO_BAD_SIZE %sUnexpected size for geoid data:%s
MODIS_E_GEO_EPH %sThe ephemeris file is defective and should be replaced%s
MODIS_E_GEO_FORMATTER %sBoth sides of the formatter electronics were turned on for scan %s
MODIS_E_GEO_INT_MIRR_ANG  %sError in GEO_interp_mirr_ang.
MODIS_E_GEO_INT_MIRR_ENC  %sError in GEO_interp_mirr_enc.
MODIS_E_GEO_MIRR_MOTION %sImplausible mirror motion data for scan %s
MODIS_E_GEO_MISSING_INPUTS %sInput value not available for %s
MODIS_E_GEO_NO_ZERO_ENCODER %sCan't extrapolate to encoder count of 0.0 for scan %s
MODIS_E_GEO_SCANNO_INPUT %sUnable to support processing of number of scans:%s
MODIS_E_GEO_VEC_UNIT3 %sError in GEO_vec_unit3 %s
MODIS_E_GEO_WRONG_PLATFORM %sIncorrect platform parameter %s
MODIS_E_INSUFFICIENT_PKTS %sInsufficient Ancillary Packet Data for interpolation. %s
MODIS_E_MISSING_OUTPUT %sOutput value not available for %s
MODIS_E_NO_VAL_ATT_SAMP  %sNo valid attitude samples.
MODIS_E_NO_VAL_AV_SAMP  %sNo valid angular velocity samples.
MODIS_E_NO_VAL_ORB_SAMP  %sNo valid orbit samples.
MODIS_E_NO_VAL_POS_SAMP  %sNo valid position samples.
MODIS_E_NO_VAL_VEL_SAMP  %sNo valid velocity samples.
MODIS_E_ONE_UNFLAG  %sOnly one unflagged value in array.
MODIS_E_PREMATURE_EOF %sPremature end of file on %s
MODIS_E_UNKNOWN_PARAMETER %sUnknown configuration parameter: %s
MODIS_E_UTIL_MAX_POLY  %sPolynomial degree greater than MAX_POLY_DEGREE.
MODIS_E_UTIL_POLY  %sPolynomial degree less than one.
MODIS_E_UTIL_VEC  %sVector has zero magnitude.
MODIS_E_WRONG_FIELD %sexpected:%s
MODIS_E_WRONG_FILE %sFirst line:%s
MODIS_E_WRONG_LUT %sLUT has incorrect version:%s
MODIS_N_GEO_EXCESS_FILES %sMore than %s files were staged.
MODIS_N_GEO_MANEUVER Return value (only): spacecraft was manuevering.
MODIS_U_GEO_BEGIN %s Seed file: 6.2 Running MOD_PR03 MODIS Geolocation main() %s
MODIS_U_GEO_END %sMOD_PR03 controlled exit, with exit code %s
MODIS_U_GEO_GRANULE_ID %sGRANULE_ID:%s
MODIS_W_GEO_EPHEMERIS %sNo s/c ephemeris/attitude files could be found for input time %s
MODIS_W_GEO_FRAMENO_INPUT %sInvalid %s
MODIS_W_GEO_MIRR_CHAN %sNo valid mirror channel data available.
MODIS_W_GEO_MISSING_MET %sMetadata could not be retrieved for LUN %s
MODIS_W_MISSING_OUTPUT %sOutput value not available for %s
MODIS_W_NO_GEO %sGranule contains no geolocatable pixels.%s
MODIS_W_SAMP_TIME_OOR  %sSample time %s
MODIS_W_SCAN_OFF  %sScan off edge of Earth. %s
# End of MODIS_35251.t
