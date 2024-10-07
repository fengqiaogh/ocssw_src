#include <stdio.h>
#include "smfio.h"
#include "GEO_geo.h"
#include "GEO_earth.h"
#include "GEO_product.h"
#include "PGS_EPH.h"
#include "PGS_SMF.h"
#include "PGS_TD.h"
#include "PGS_MODIS_35251.h"

#define TIMEPOS sizeof("2000-01-01")

PGSt_SMF_status GEO_check_ea_headers(
	PGSt_double	in_time,
	PGSt_scTagInfo	*scTagInfo
)
/******************************************************************************
!C

!Description:   
      Determines whether or not the appropriate ephemeris and attitude files 
      have been staged to cover a specified time. Note: for the purpose of 
      this routine, a file is treated as covering the entire two hour (or 
      for PM1EPH, 24 hour) time period it was intended to cover, regardless 
      of any data gaps that might reduce the actual coverage.

!Input Parameters:
      in_time	The time being queried (TAI seconds)
      scTagInfo	Information about the spacecraft whose headers are needed.

!Output Parameters:
      scTagInfo	A structure containing information about the spacecraft.

Return Values:
      MODIS_E_BAD_INPUT_ARG	  If scTagInfo is null
      MODIS_E_GEO                 If any subroutine failed.
      MODIS_N_GEO_EXCESS_FILES    If too many files were staged.
      MODIS_E_GEO_MISSING_INPUTS  ephemeris/attitude files weren't staged
      PGS_S_SUCCESS               Otherwise

Externally Defined:
      MAX_EA_FILES		  "GEO_product.h"
      MODIS_E_BAD_INPUT_ARG	  "PGS_MODIS_35251.h"
      MODIS_E_GEO                 "PGS_MODIS_35251.h"
      MODIS_N_GEO_EXCESS_FILES    "PGS_MODIS_35251.h"
      MODIS_E_GEO_MISSING_INPUTS  "PGS_MODIS_35251.h"
      PGS_S_SUCCESS               "PGS_SMF.h"
      PGSd_EOS_AM                 "PGS_TD.h"
      PGSd_EOS_AM1                "PGS_TD.h"
      PGSd_EOS_PM                 "PGS_TD.h"
      PGSe_TAG_SEARCH             "PGS_TD.h"
      TIMECODEASIZE               "smfio.h"

Called by:
      GEO_interp_ephemeris_attitude()

Routines Called:
      PGS_EPH_getAttitHeaders      "PGS_EPH.h"
      PGS_EPH_getEphemHeaders      "PGS_EPH.h"
      PGS_TD_TAItoUTC              "PGS_TD.h"
      PGS_TD_UTCtoTAI              "PGS_TD.h"

!Revision History:
      $Log: GEO_check_ea_headers.c,v $
      Revision 6.3  2010/05/27 21:02:59  kuyper
      Changed MAX_EA_FILES to come from a single header file.

      Revision 6.2  2010/04/09 20:37:32  kuyper
      Removed parameter, corrected places where it was used.

      Revision 6.1  2010/03/31 19:58:09  kuyper
      Helped resolve Bug 2473 by changing scTagInfo into an output argument.

      James Kuyper	james.kuyper@sigmaspace.com

      Revision 5.1  2008/12/16 18:06:35  kuyper
      Increased number of ephemeris and attitude files allowed, to allow use of
        the 5-minute files used in Near Real Time production.

      Revision 4.4  2003/11/05 16:24:58  kuyper
      Initialize header_array to NULL before retrieving headers.

      Revision 4.3  2003/06/24 18:36:37  vlin
      Changed to match new interface for PGS_EPH_getEphemHeaders and
      PGS_EPH_getAttitHeaders in SDPTK V5.2.9

      Revision 4.2  2003/05/15 16:30:08  vlin
      change EXCESS_FILES from _E_ to _N_

      Revision 4.1  2003/04/10 16:46:33  vlin
      updated after code walkthrough

      Revision 4.0  2003/03/05 21:12:31  vlin
      initial revision
      vlin@saicmodis.com

Requirements:

!Team-unique Header:
      This software is developed by the MODIS Science Data Support
      Team for the National Aeronautics and Space Administration,
      Goddard Space Flight Center, under contract NAS5-32373.

!END
******************************************************************************/

{
  static PGSt_integer    ephem_files, attit_files;
  static PGSt_hdrSummary ephem_summary[MAX_EA_FILES];
  static PGSt_hdrSummary attit_summary[MAX_EA_FILES];
  static int             initialized = 0;
  PGSt_integer           file;
  PGSt_hdrSummary        *header_array=(PGSt_hdrSummary*)NULL;
  PGSt_integer           lendcheck;
  char                   asciiUTC[TIMECODEASIZE];
  PGSt_SMF_status        status = PGS_S_SUCCESS;
  int                    hour, noon_time, zero_time;
  PGSt_double            offset;
  char                   msg[128];
  char                   filefunc[] = __FILE__ ", GEO_check_ea_headers";

  if(scTagInfo == NULL)
  {
      modsmf(MODIS_E_BAD_INPUT_ARG, "scTagInfo is null", filefunc);
      return MODIS_E_BAD_INPUT_ARG;
  }

  if (initialized == 0)
  {
      if (PGS_EPH_getEphemHeaders(scTagInfo, &header_array, &ephem_files,
	  &lendcheck) != PGS_S_SUCCESS) {
	  modsmf(MODIS_E_GEO, "PGS_EPH_getEphemHeaders()", filefunc);
	  status = MODIS_E_GEO_MISSING_INPUTS;
      }

      if (ephem_files > MAX_EA_FILES) {
	  sprintf(msg, "%d ephemeris", ephem_files);
	  modsmf(MODIS_N_GEO_EXCESS_FILES, msg, filefunc);
	  return MODIS_N_GEO_EXCESS_FILES;
      }

      for (file = 0; file < ephem_files; file++) {
      /* Round start times downward */
	   if (PGS_TD_TAItoUTC(header_array[file].startTAI93, asciiUTC) !=
	       PGS_S_SUCCESS) {
	       sprintf(msg, "PGS_TD_TAItoUTC(%.6f) at ephemeris start",
		       header_array[file].startTAI93);
	       modsmf(MODIS_E_GEO, msg, filefunc);
	       status = MODIS_E_GEO;
	   }

	   hour = atoi(asciiUTC+TIMEPOS);
	   if (scTagInfo->spacecraftTag == PGSd_EOS_AM)
	   {	/* Round down to nearest even hour */
	       sprintf(asciiUTC+TIMEPOS,"%02d:00:00.000000Z", hour);
	       if (hour%2)     /* hour is odd */
		   offset = -3600.0;
	       else
		   offset = 0.0;
	   }
	   else {  /* Aqua ephemeris files - round down to the nearest noon */
	       sprintf(asciiUTC+TIMEPOS,"12:00:00.000000Z");
	       if (hour < 12)
		   offset = -86400.0;
	       else
		   offset = 0.0;
	   }

	   if (PGS_TD_UTCtoTAI(asciiUTC, &ephem_summary[file].startTAI93) !=
	       PGS_S_SUCCESS) {
	       sprintf(msg, "PGS_TD_UTCtoTAI(%s) at ephemeris start", asciiUTC);
	       modsmf(MODIS_E_GEO, msg, filefunc);
	       status = MODIS_E_GEO;
	   }

	   ephem_summary[file].startTAI93 += offset;

      /* Round stop times upward. */
	   if (PGS_TD_TAItoUTC(header_array[file].stopTAI93, asciiUTC) !=
	       PGS_S_SUCCESS) {
	       sprintf(msg, "PGS_TD_TAItoUTC(%.6f) at ephemeris stop", 
		       header_array[file].stopTAI93);
	       modsmf(MODIS_E_GEO, msg, filefunc);
	       status = MODIS_E_GEO;
	   }

	   hour = atoi(asciiUTC+TIMEPOS);
	   noon_time = strcmp(asciiUTC+TIMEPOS, "12:00:00.000000Z");
	   zero_time = strcmp(asciiUTC+TIMEPOS+2, ":00:00.000000Z");

	   if (scTagInfo->spacecraftTag ==  PGSd_EOS_AM1)
	   {	/* Round upward to nearest even hour */
	      if (hour%2)
		  offset = 3600.0;
	      else if (zero_time == 0)
		  offset = 0.0;
	      else
		  offset = 7200.0;
	      sprintf(asciiUTC+TIMEPOS,"%02d:00:00.000000Z", hour);
	   }
	   else {  /* Aqua ephemeris files - round upward to the nearest noon */
	       if (hour < 12 || noon_time == 0) 
		   offset = 0.0;
	       else 
		   offset = 86400.0;
	       sprintf(asciiUTC+TIMEPOS,"12:00:00.000000Z");
	   }

	   if (PGS_TD_UTCtoTAI(asciiUTC, &ephem_summary[file].stopTAI93) != 
	       PGS_S_SUCCESS) {
	       sprintf(msg, "PGS_TD_UTCtoTAI(%s) at ephemeris stop", asciiUTC);
	       modsmf(MODIS_E_GEO, msg, filefunc);
	       status = MODIS_E_GEO;
	   }

	   ephem_summary[file].stopTAI93 += offset;

      } /* End for */

      header_array = (PGSt_hdrSummary*)NULL;
      if (PGS_EPH_getAttitHeaders(scTagInfo, &lendcheck, &header_array, 
	  &attit_files) != PGS_S_SUCCESS) {
	  modsmf(MODIS_E_GEO, "PGS_EPH_getAttitHeaders()", filefunc);
	  status = MODIS_E_GEO_MISSING_INPUTS;
      }

      if (attit_files > MAX_EA_FILES) {
	  sprintf(msg, "%d attitude", attit_files);
	  modsmf(MODIS_N_GEO_EXCESS_FILES, "attitude", filefunc);
	  return MODIS_N_GEO_EXCESS_FILES;
      }

      for (file = 0; file < attit_files; file++) { 
      /* Round start times downward */
	   if (PGS_TD_TAItoUTC(header_array[file].startTAI93, asciiUTC) !=
	       PGS_S_SUCCESS) {
	       sprintf(msg, "PGS_TD_TAItoUTC(%.6f) at attitude start",
		       header_array[file].startTAI93);
	       modsmf(MODIS_E_GEO, msg, filefunc);
	       status = MODIS_E_GEO;
	   }

	   hour = atoi(asciiUTC+TIMEPOS);
	   sprintf(asciiUTC+TIMEPOS,"%02d:00:00.000000Z", hour);

	   if (hour%2)     /* hour is odd */
	       offset = -3600.0;
	   else
	       offset = 0.0;

	   if (PGS_TD_UTCtoTAI(asciiUTC, &attit_summary[file].startTAI93) !=
	       PGS_S_SUCCESS) {
	       sprintf(msg, "PGS_TD_UTCtoTAI(%s) at attitude start", asciiUTC);
	       modsmf(MODIS_E_GEO, msg, filefunc);
	       status = MODIS_E_GEO;
	   }

	   attit_summary[file].startTAI93 += offset;

      /* Round stop times upward. */
	   if (PGS_TD_TAItoUTC(header_array[file].stopTAI93, asciiUTC) !=
	       PGS_S_SUCCESS) {
	       sprintf(msg, "PGS_TD_TAItoUTC(%.6f) at attitude stop", 
		       header_array[file].stopTAI93);
	       modsmf(MODIS_E_GEO, msg, filefunc);
	       status = MODIS_E_GEO;
	   }

	   hour = atoi(asciiUTC+TIMEPOS);
	   noon_time = strcmp(asciiUTC+TIMEPOS, "12:00:00.000000Z");
	   zero_time = strcmp(asciiUTC+TIMEPOS+2, ":00:00.000000Z");

	   /* Round upward to nearest even hour. */
	   if (hour%2)
	       offset = 3600.0;
	   else if (zero_time == 0)
	       offset = 0.0;
	   else
	       offset = 7200.0;
	   sprintf(asciiUTC+TIMEPOS, "%02d:00:00.000000Z", hour);

	   if (PGS_TD_UTCtoTAI(asciiUTC, &attit_summary[file].stopTAI93) != 
	       PGS_S_SUCCESS) {
	       sprintf(msg, "PGS_TD_UTCtoTAI(%s) at attitude stop", asciiUTC);
	       modsmf(MODIS_E_GEO, msg, filefunc);
	       status = MODIS_E_GEO;
	   }

	   attit_summary[file].stopTAI93 += offset;

      } /* End for */

      initialized = 1;
  }

  if (ephem_files < 1 || attit_files < 1) 
      return MODIS_E_GEO_MISSING_INPUTS;

  for (file = 0; file < ephem_files; file++) {
      if (ephem_summary[file].startTAI93 <= in_time &&
          in_time <= ephem_summary[file].stopTAI93)
          break;
  }

  if (file > ephem_files - 1) { /* A needed ephemeris file was not staged. */
      sprintf(msg, "%.6f from ephemeris files.", in_time);
      modsmf(MODIS_E_GEO_MISSING_INPUTS, msg, filefunc);
      status = MODIS_E_GEO_MISSING_INPUTS;
  }

  for (file = 0; file < attit_files; file++) {
       if (attit_summary[file].startTAI93 <= in_time &&
           in_time <= attit_summary[file].stopTAI93)
           break;
  }

  if (file > attit_files - 1) { /* A needed attitude file was not staged. */
      sprintf(msg, "%.6f from attitude files.", in_time);
      modsmf(MODIS_E_GEO_MISSING_INPUTS, msg, filefunc);
      status = MODIS_E_GEO_MISSING_INPUTS;
  }

  return status;
}
