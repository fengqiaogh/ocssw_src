/*
!C-INC*************************************************************************
!Description:   Defines the information needed to write the geolocation
                product.

!Input Parameters: N/A

!Output Parameters: N/A

!Revision History:
 * $Log: GEO_product.h,v $
 * Revision 6.12  2012/06/26 21:02:00  kuyper
 * Added DOI authority metadata values.
 *
 * Revision 6.11  2012/06/21 20:31:51  kuyper
 * Added Digital Object Identifier.
 *
 * Revision 6.10  2011/02/09 19:04:43  kuyper
 * Added WaterPresent SDS.
 *
 * Revision 6.9  2010/09/09 18:57:59  kuyper
 * Updated to use macros from mapiL1Bgeo.h that were added in M-API 5.0.0 or
 *   later.
 *
 * Revision 6.8  2010/06/30 20:48:45  kuyper
 * Backed out changes dependent upon M-API 6.0.0, which isn't quite ready yet.
 *
 * Revision 6.7  2010/05/28 18:02:32  kuyper
 * Moved duplicate and necessarily matching definitions of MAX_EA_FILES and
 *   MAX_EA_INPUTS from GEO_check_ea_headers() and GEO_write_input_metadata() into
 *   this header
 * Reverted more macros to form compatible with M-API 2.3.4.
 *
 * Revision 6.6  2010/05/19 20:17:34  kuyper
 * Removed references to macros that were added in M-API 6.0.0; this delivery will
 *   have to use M-API 2.3.4.
 * Dropped macros for quality flag values that will never actually appear in the
 *   product files.
 *
 * Revision 6.5  2010/04/01 20:28:24  james
 * Corrected to make quality flag values explicitly unsigned.
 *
 * Revision 6.4  2010/03/31 17:02:45  kuyper
 * Helped resolve Bug 1969 by defining THERMCORR.
 * Helped resolve Bug 2249 by defining TERRAIN_CORRECTION
 * Helped resolve Bug 2472 by defining names for the ephemeris and attitutde
 *   quality SDSs.
 *
 * Revision 6.3  2009/05/31 20:08:32  kuyper
 * Corrected capitalization of offsets.
 *
 * Revision 6.2  2009/04/06 20:18:18  kuyper
 * Added NUMBANDS
 *
 * Revision 6.1  2009/03/20 21:46:19  kuyper
 * Added dimension names and SDS names for high resolution offsets.
 * Added an alias for the swath name.
 * Added names for band number and fractional offsets file attributes.
 *
 * James Kuyper Jr.	James.R.Kuyper@nasa.gov
 *
 * Revision 5.2  2005/04/19 21:47:39  kuyper
 * Changed header guard macro name to avoid reserved name space.
 *
 * Revision 5.1  2004/08/26 16:18:21  vlin
 * Macro "MANEUVERS" defined.
 *
 * Revision 4.5  2003/08/01 21:20:33  kuyper
 * Corrected typo.
 *
 * Revision 4.4  2003/07/28 19:49:43  kuyper
 * Added REPROCESSINGPLANNED.
 * Reorganized, and removed duplicate entry.
 *
 * Revision 4.3  2003/04/24 19:57:44  kuyper
 * Removed GEO_TERRA, GEO_AQUA. Reorganised other #defines.
 *
 * Revision 4.2  2002/12/16 17:56:06  kuyper
 * Corrected definition of NUMENC.
 *
 * Revision 4.1  2002/11/22 21:31:47  kuyper
 * Added names for temperatures vdata.
 * Added roll, pitch, and yaw SDS attributes.
 * Added dimension names.

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END**************************************************************************
*/

#ifndef GEO_PRODUCT_H
#define GEO_PRODUCT_H
#include "mapi.h"
#include "mapiL1Bgeo.h"

/* Define ECS metadata Inventory and Archive */

#define ADDITIONALATTRIBUTENAME		MCORE_ADDATTRIBUTENAME
#define ARCHIVEMETADATA			MECS_ARCHIVE
#define ASSOCIATEDPLATFORMSHORTNAME	MCORE_APSHORTNAME
#define AUTOMATICQUALITYFLAG		MCORE_AUTO_QUALITY
#define AUTOMATICQUALITYFLAGEXPLANATION	MCORE_AUTO_QUALITY_FLG
#define COREMETADATA			MECS_CORE
#define DAYNIGHTFLAG			MCORE_DAYNIGHTFLAG
#define EASTBOUNDINGCOORDINATE		MCORE_EAST_BOUND
#define EQUATORCROSSINGDATE		MCORE_EQUATCROSSINGDATE
#define EQUATORCROSSINGLONGITUDE	MCORE_EQUATCROSSINGLONG
#define EQUATORCROSSINGTIME		MCORE_EQUATCROSSINGTIME
#define EXCLUSIONGRINGFLAG		MCORE_EXCLUS_GRING_FLG
#define GRANULENUMBER			MPROD_GRANULE_NUM
#define GRINGPOINTLATITUDE		MCORE_GRING_POINT_LAT
#define GRINGPOINTLONGITUDE		MCORE_GRING_POINT_LON
#define GRINGPOINTSEQUENCENO		MCORE_GRING_POINT_NUM
#define INPUTPOINTER			MCORE_INPUT_POINTER
#define L1A_COREMETADATA		MECS_CORE
#define LOCALGRANULEID			MCORE_LOCALGRANULEID
#define LOCALINPUTGRANULEID		MCORE_LOCALINPUTGRANULEID
#define LOCALVERSIONID			MCORE_LOCALVERSIONID
#define NORTHBOUNDINGCOORDINATE		MCORE_NORTH_BOUND
#define ORBITNUMBER			MCORE_ORBIT_NUM
#define PARAMETERNAME			MCORE_PARAMETERNAME
#define PARAMETERVALUE			MCORE_PARAMETERVALUE
#define PGEVERSION			MCORE_PGEVERSION
#define PROCESSINGENVIRONMENT		MCORE_PROCESSING_ENV
#define PRODUCTIONDATETIME		MCORE_PRODUCTIONDATETIME
#define QAPERCENTMISSINGDATA		MCORE_PERCENT_MISSING
#define QAPERCENTOUTOFBOUNDSDATA	MCORE_PERCENT_OUT
#define RANGEBEGINNINGDATE		MCORE_RANGE_BEG_DATE
#define RANGEBEGINNINGTIME		MCORE_RANGE_BEG_TIME
#define RANGEENDINGDATE			MCORE_RANGE_ENDING_DATE
#define RANGEENDINGTIME			MCORE_RANGE_ENDING_TIME
#define REPROCESSINGACTUAL		MCORE_ACTUALLY_REDONE
#define REPROCESSINGPLANNED		MCORE_TO_BE_REDONE
#define SHORTNAME			MCORE_SHORT_NAME
#define SOUTHBOUNDINGCOORDINATE		MCORE_SOUTH_BOUND
#define VERSIONID			MCORE_VERSIONID
#define WESTBOUNDINGCOORDINATE		MCORE_WEST_BOUND
#define IDENTIFIERPRODDOI		"identifier_product_doi"
#define IDENTIFIERPRODDOIAUTH		"identifier_product_doi_authority"
#define DOIAUTHORITY			"http://dx.doi.org"

/* Define SDS names for geolocation HDF format */
#define SWATH_NAME	M03_SWATH_NAME

  /* Group names (none defined at present */
#define SCAN_GRP	""
#define SCAN_META_GRP	""
#define PARM_GRP	""

  /* Processing and geometric parameters */
#define FOCAL_LENGTH		M03FOCAL_LENGTH
#define BAND_NUMBER		M03BAND_NUMBER
#define BAND_POSITION		M03BAND_POSITION
#define DETECTOR_SPACE		M03DETECTOR_SPACE
#define DETECTOR_OFFSETS	M03DETECTOR_OFFSETS
#define T_OFFSET		M03T_OFFSET
#define NUM_SAMPLES		M03NUM_SAMPLES

  /* Scan line metadata */
#define S_NUM	     M03SCAN_NO
#define S_TYPE       M03SCAN_TYPE
#define EV_FRAMES    M03EV_FRAMES
#define SD_FRAMES    M03SD_FRAMES
#define SV_FRAMES    M03SV_FRAMES
#define EVTIME       M03EV_START_TIME
#define SDTIME       M03SD_START_TIME
#define SVTIME       M03SV_START_TIME
#define SCTIME       M03EV_CENTER_TIME
#define MSIDE	     M03MIR_SIDE
#define SUN_ZENITH   M03SD_SUN_ZENITH
#define SUN_AZIMUTH  M03SD_SUN_AZIMUTH
#define MOON_VECTOR  M03MOON_VECTOR
#define L1_QUALITY   M03L1_SCAN_QUALITY
#define GEO_QUALITY  M03GEO_SCAN_QUALITY
#define ORB_POS      M03ORB_POS
#define ORB_VEL      M03ORB_VEL
#define T_INST2ECR   M03T_INST2ECR
#define ATTIT_ANG    M03ATTITUDE_ANGLES
#define SUN_REF      M03SUN_REF
#define NUM_IMPULSE  M03NUM_IMPULSE
#define IMPULSE_ENC  M03IMPULSE_ENC
#define IMPULSE_TIME M03IMPULSE_TIME
#define NSCANS_10    M03NSCANS_10
#define NSCANS_20    M03NSCANS_20
#define NSCANS_40    M03NSCANS_40
#define MFRAMES      M03MFRAMES
#define MFRAMES_2    M03MFRAMES_2
#define MFRAMES_4    M03MFRAMES_4
#define THERMCORR    M03THERM_CORR
#define ATT_QUALITY  M03ATT_QUALITY
#define EPH_QUALITY  M03EPH_QUALITY

  /* Quality Flag Values. */
#define EPH_OVERALL		0x00000001U
#define EPH_DATA		0x00000002U
#define EPH_REDLO		0x00000004U
#define EPH_YELLOWLO		0x00000008U
#define EPH_YELLOWHI		0x00000010U
#define EPH_REDHI		0x00000020U
#define EPH_LONG_FOLLOW		0x00000040U
#define	EPH_SHORT_FOLLOW	0x00000080U
#define EPH_SHORT_PRECEED	0x00000100U
#define EPH_LONG_PRECEED	0x00000200U
#define EPH_REPAIRED		0x00000400U
#define EPH_QUALITY_UNAVAIL	0x00000800U
/* PGSd_NO_DATA:		0x00001000U	Used as fill value. */	
/* PGSd_INTERPOLATED_POINT	0x00004000U	*/
/* PGSd_PLATFORM_FATAL		0x00010000U	*/

  /* Quality Flag SDS Attribute Names. */
#define QFL_OVERALL	M03QFL_OVERALL
#define QFL_DATASUM	M03QFL_DATASUM
#define QFL_RED_LO	M03QFL_RED_LO
#define QFL_YEL_LO	M03QFL_YEL_LO
#define QFL_YEL_HI	M03QFL_YEL_HI
#define QFL_RED_HI	M03QFL_RED_HI
#define QFL_LONG_FOLL	M03QFL_LONG_FOLL
#define QFL_SHRT_FOLL	M03QFL_SHRT_FOLL
#define QFL_SHRT_PREC	M03QFL_SHRT_PREC
#define QFL_LONG_PREC	M03QFL_LONG_PREC
#define QFL_REPAIRED	M03QFL_REPAIRED
#define QFL_QFL_PROB	M03QFL_QFL_PROB
#define QFL_INTERP	M03QFL_INTERP

  /* Pixel-level Geolocation Data */
#define LONGITUDE	M03LONGITUDE
#define LATITUDE	M03LATITUDE
#define HEIGHT		M03HEIGHT
#define SEN_ZENITH	M03SENSOR_ZEN
#define SEN_AZIMUTH	M03SENSOR_AZ
#define RANGE		M03RANGE
#define SOL_ZENITH	M03SOLAR_ZENITH
#define SOL_AZIMUTH	M03SOLAR_AZIMUTH
#define LAND_SEAMASK	M03LAND_SEAMASK
#define WATER_PRESENT   "WaterPresent"
#define GFLAGS		M03GFLAGS

  /* High-resolution Geolocation data */
#define SCAN_OFFSET	M03SCAN_OFFSET
#define TRACK_OFFSET	M03TRACK_OFFSET
#define HEIGHT_OFFSET	M03HEIGHT_OFFSET

  /* Geolocation Specific Metadata */
#define ATTIT_INPUT	M03ATTIT_INPUT
#define BAD_PCKTS	M03BAD_PACKETS
#define CUM_GFLAGS	M03CUM_GFLAGS
#define DISCRD_PCKTS	M03DISCARD_PACKETS
#define EPHEM_INPUT	M03EPHEM_INPUT
#define INCOMP_SCANS	M03INCOMP_SCANS
#define MAX_EA_FILES	3
#define MAX_EA_INPUTS	6
#define MAX_EFRM        M03MAX_EARTH_FRAMES
#define MAX_SFRM        M03MAX_SD_FRAMES
#define MAX_SV_FRM      M03MAX_SV_FRAMES
#define MISS_PCKTS	M03MISSING_PACKETS
#define NUMSCN          M03NUMBER_OF_SCANS
#define PARVERS		M03PARAM_VERS
#define POLAR_MOTION	M03POLAR_MOTION
#define UTCPOLE_FH	M03UTCPOLE
#define FROFF		"HDFEOS_FractionalOffset_"
#define FROFF_SCAN_20	FROFF NSCANS_20 "_" SWATH_NAME
#define FROFF_FRAME_2	FROFF MFRAMES_2 "_" SWATH_NAME
#define FROFF_SCAN_40	FROFF NSCANS_40 "_" SWATH_NAME
#define FROFF_FRAME_4	FROFF MFRAMES_4 "_" SWATH_NAME
#define TERRAIN_CORRECTION	M03TERR_CORR

/* Vdata names. */
#define AVERAGE_TEMPERATURES	M03AVG_TEMP
#define TA_RC_SMIR_CFPA		M03TEMP_TRS
#define TP_AO_SMIR_OBJ		M03TEMP_TAS
#define TP_MF_CALBKHD_SR	M03TEMP_TMC
#define TP_MF_Z_BKHD_BB		M03TEMP_TMZ
#define TP_SA_RCT1_MIR		M03TEMP_TSR
#define TP_SR_SNOUT		M03TEMP_TSS

/* Input source for spacecraft kinematic, and possible values. */
#define EA_SOURCE		M03EPHEMERIS_ATT_SRCE
#define TOOLKIT_EA		"SDP Toolkit"
#define L1A_PKT_EA		"MODIS Packet"
#define GEO_EST_RMS_ERROR	M03GEO_EST_ERROR

/* SDS attributes */
#define SCANDIM		M03SCAN_D
#define TRACKDIM	M03TRACK_D
#define ROLL_ELEM	M03ROLL_ELEM
#define PITCH_ELEM	M03PITCH_ELEM
#define YAW_ELEM	M03YAW_ELEM

/* Dimension names */
#define NSCANS		M03NSCANS
#define VECDIM		M03VECDIM
#define NUMQUAL 	M03NUMQUAL
#define NUMENC		M03NUMENC
#define NUMBANDS	M03NUMBANDS
#endif
