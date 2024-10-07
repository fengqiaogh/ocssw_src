/* file: L1a_data.h */

/*
!C-INC*************************************************************************
!Description:   Defines the information needed by the geolocation software 
                to read the Level 1A product.

!Input Parameters: N/A

!Output Parameters: N/A

!Revision History:
		$Log: L1a_data.h,v $
		Revision 6.2  2010/06/30 20:18:00  kuyper
		Backed out a change that relied upon M-API 6.0.0, which is not quite ready.

		Revision 6.1  2010/05/04 19:29:49  kuyper
		Added MAJCYCALL1 and SS_CP_MODE.

		Revision 5.1  2005/03/16 21:36:28  kuyper
		Changed header guard macro name to avoid reserved name space.

		Revision 4.2  2003/12/17 21:06:37  kuyper
		Removed obsolete macro.

		Revision 4.1  2003/04/24 19:28:40  kuyper
		Removed obsolete ESDT macros.

		Revision 2.12  2001/04/02 16:43:25  seaton
		Added macros for GEO_prepare_l1a_data.c

 * Revision 2.11  2001/03/31  15:40:34  seaton
 * Entered macros used by GEO_read_L1Apacket_data.c.
 *
 * Revision 2.10  2001/01/17  13:45:29  vlin
 * Added "PRODUCTIONHISTORY" macro
 *
 * Revision 2.9  2000/08/18  00:45:54  kuyper
 * Changed to support Aqua as well at Terra data.
 *
 * Revision 2.8  2000/08/14  20:18:24  fhliang
 * Added macro CORE_ASSOCIATEDPLATFORMSHORTNAME.
 *
 * Revision 2.7  2000/06/13  18:01:05  lma
 * changed definaition of macro SCAN_TYPES to M01SCAN_TYPE
 *
 * Revision 2.6  2000/06/09  18:52:58  lma
 * changed definaition of macro SCAN_TYPES to M02SW_SCAN_TYPE.
 *
 * Revision 2.5  2000/05/04  15:38:00  lma
 * added macro SCAN_TYPES
 *
 * Revision 2.4  1999/03/12  17:48:37  kuyper
 * Capitalized Prolog Sections
 *
 * Revision 2.3  1999/02/05  17:59:04  seaton
 * Added M-API #defines for the Spacecraft Ancillary Data fields
 * used in GEO_read_L1Apacket_data.c.
 *
 * Revision 2.2  1998/03/04  03:34:17  jjb
 * Added L1A ESDT macro.
 *
 * Revision 2.1  1997/10/21  18:15:47  kuyper
 * Returned from ClearCase
 *
 * Revision 1.6  1997/07/18  21:58:00  kuyper
 * Baselined Version 1
 *
 * Revision 1.6  1997/03/26  19:12:34  fhliang
 * Initial revision of SDST delivery of L1a_data.h.
 *
		Revision 1.5  1997/01/14 21:17:19  kuyper
		Added macros neede by GEO_read_L1A*data() functions.

		James Kuyper (kuyper@ltpmail.gsfc.nasa.gov)

                6/20/95
                Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
                Finished coding

                9/20/95
                Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
                Modified to update Level 1A field names
		
                10/31/95
                Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
                Modified to update Level 1A field names
		

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END**************************************************************************
*/

#ifndef L1A_DATA_H
#define L1A_DATA_H
#include "mapiL1A.h"

/* Define the group names for the Level 1A data */

#define L1A_SCAN_META_GRP ""
#define L1A_ENGINEERING_GRP ""

/* Define the SDS names for the Level 1A data */

# define DISCARD_PACKETS	M01DISCARD_PACKETS
# define SCAN_START_TIME	M01EV_START_TIME
# define SCAN_TYPES		M01SCAN_TYPE
# define EARTH_SECTOR_FRAMES	M01FRAME_COUNT_ARRAY
# define INCOMPL_SCANS		M01INCOMPL_SCANS
# define MAX_EARTH_FRAMES	M01MAX_EARTH_FRAMES
# define MAX_SD_FRAMES		M01MAX_SD_FRAMES
# define MAX_SV_FRAMES		M01MAX_SV_FRAMES
# define MIRROR_SIDE		M01MIRROR_SIDE
# define MISSING_PACKETS	M01MISSING_PACKETS
# define NUMBER_OF_SCANS	M01NUMBER_OF_SCANS
# define PACKTS_BAD_CRC		M01PACKTS_BAD_CRC
# define EARTH_ENCODER_TIMES	M01RAW_MIR_ENC
# define SC_ANCILLARY_DATA	M01RAW_SC_ANCIL
# define VIEW_SECTOR_START	M01RAW_VS_START
# define SCAN_NUMBER		M01SCAN_NUMBER
# define SCAN_QUALITY_ARRAY	M01SCAN_QUALITY_ARRAY
# define SD_START_TIME		M01SD_START_TIME
# define SV_START_TIME		M01SV_START_TIME
# define SCIENCE_STATE          M01SCIENCE_STATE
# define SCIENCE_ABNORM         M01SCIENCE_ABNORM
# define CORE_RANGE_BEG_DATE    MCORE_RANGE_BEG_DATE
# define CORE_RANGE_BEG_TIME    MCORE_RANGE_BEG_TIME
# define CORE_RANGE_ENDING_DATE MCORE_RANGE_ENDING_DATE
# define CORE_RANGE_ENDING_TIME MCORE_RANGE_ENDING_TIME
# define CORE_DAYNIGHTFLAG      MCORE_DAYNIGHTFLAG
# define CORE_LOCALGRANULEID    MCORE_LOCALGRANULEID
# define CORE_PARAMETERVALUE    MCORE_PARAMETERVALUE
# define CORE_ASSOCIATEDPLATFORMSHORTNAME MCORE_APSHORTNAME
# define MECS_PRODHISTORY	"PRODUCTIONHISTORY"

/* Define the S/C Ancillary Data names for level 1A Data */
#define TIME_STAMP              M01TIME_STAMP
#define SC_POSITION_X          M01SC_POSITION_X
#define SC_POSITION_Y          M01SC_POSITION_Y
#define SC_POSITION_Z          M01SC_POSITION_Z
#define SC_VELOCITY_X          M01SC_VELOCITY_X
#define SC_VELOCITY_Y          M01SC_VELOCITY_Y
#define SC_VELOCITY_Z          M01SC_VELOCITY_Z
#define ATTITUDE_ANGLE_ROLL    M01ATTITUDE_ANGLE_ROLL
#define ATTITUDE_ANGLE_PITCH   M01ATTITUDE_ANGLE_PITCH
#define ATTITUDE_ANGLE_YAW     M01ATTITUDE_ANGLE_YAW
#define ATTITUDE_RATE_ROLL     M01ATTITUDE_RATE_ROLL
#define ATTITUDE_RATE_PITCH    M01ATTITUDE_RATE_PITCH
#define ATTITUDE_RATE_YAW      M01ATTITUDE_RATE_YAW
#define PRIOR_SC_ANCIL_DATA    M01PRIOR_SC_ANCIL_DATA
#define CURR_SC_ANCIL_DATA     M01CURR_SC_ANCIL_DATA
#define CR_FR_A_ON             M01CR_FR_A_ON
#define CR_FR_B_ON             M01CR_FR_B_ON
#define CR_SA_A_SCAN_ON        M01CR_SA_A_SCAN_ON
#define CR_SA_B_SCAN_ON        M01CR_SA_B_SCAN_ON
#define LAST_VALID_SCAN        M01LAST_VALID_SCAN
#define MAJCYCALL1		M01MAJCYCALL1
#define MAJCYC3COF7            M01MAJCYC3COF7
#define MAJCYC5BOF7            M01MAJCYC5BOF7
#define	SS_CP_MODE		M01SS_CP_MODE
#endif

