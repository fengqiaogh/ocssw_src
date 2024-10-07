#include <errno.h>
#include "PGS_MODIS_35251.h"
#include "PGS_CSC.h"
#include "PGS_SMF.h"
#include "PGS_IO.h"
#include "smfio.h"
#include "GEO_output.h"
#define MAXMSGLEN 120

PGSt_SMF_status GEO_get_utcpole_metadata(
	const ECS_metadata_struct	*ECS_metadata,
	utcpole_metadata_struct		*utcpole_metadata)
/*****************************************************************************
!C
 
!Description:	
    Fills in utcpole_metadata with values read/interpolated from
    the utpole.data file.
 
!Input Parameters:
 	ECS_metadata		Used to obtain RangeBeginningDateTime
 
!Output Parameters:
 	utcpole_metadata	Contains metadata identifying the utcpole file.
 
Return Parameters:
 	MODIS_E_BAD_INPUT_ARG	If either argument is NULL
 	PGS_E_UNIX		If the call to fgets() fails.
 	MODIS_E_GEO		If any other subroutine fails.
	MODIS_E_PREMATURE_EOF   if the utcpole.dat file is empty.
 	PGS_S_SUCCESS		Otherwise


Externally Defined:
 	errno			<errno.h>
 	MODIS_E_BAD_INPUT_ARG	"PGS_MODIS_35251.h"
 	MODIS_E_GEO		"PGS_MODIS_35251.h"
	MODIS_E_PREMATURE_EOF   "PGS_MODIS_35251.h"
	PGSCSC_W_PREDICTED_UT1	"PGS_CSC.h"
	PGSd_IO_Gen_Read        "PGS_IO.h"
 	PGSd_UTCPOLE		"PGS_CSC.h"
 	PGS_E_UNIX		"PGS_SMF.h"
 	PGS_S_SUCCESS		"PGS_SMF.h"
	TIMECODEASIZE           "smfio.h" 
 
Called by:
 	GEO_write_granule_metadata
 
Routines Called:
 	modsmf			To log status messages.
 	PGS_IO_Gen_Open		To open a file, identified by a LUN
 	PGS_IO_Gen_Close	To close afile opened with PGS_IO_GenOpen
 	PGS_SMF_SetUnixMessage	To log Unix error messages as status messages.
 	PGS_TD_UTCtoUTCjd	To calculate the Julian Day for a UTC time.
 	PGS_CSC_UTC_UT1Pole	To calculate utcpole value for a given time.
 
!Revision History:
$Log: GEO_get_utcpole_metadata.c,v $
Revision 4.2  2003/08/28 16:24:04  kuyper
Corrected prolog.

Revision 4.1  2003/02/21 22:40:11  kuyper
Corrected to use void* pointers with %p format code.

Revision 3.4  2001/06/12 11:12:34  kuyper
Corrected to treat PGSCSC_W_PREDICTED_UT1 as a normal return from
  PGS_CSC_UTC_UT1Pole.

Revision 3.3  2001/03/28  19:34:20  pliu
Corrected handling of error returns from fgets().

Revision 3.2  2001/03/15  21:43:08  pliu
Corrected based on walkthrough. Change direct use of 79 to sizeof.

Revision 3.1  2001/03/01  19:43:57  pliu
Initial revision.


Ping Liu  (pliu@ltpmail.gsfc.nasa.gov)


Requirements:
 	PR03-F-4.1-1.6
 	PR03-I-1
 	PR03-I-2
 	PR03-S-1
 
!Team-unique Header:
 	This software is developed by the MODIS Science Data Support
 	Team for the National Aeronautics and Space Administration,
 	Goddard Space Flight Center, under contract NAS5-32373.
 
 References and Credits

!END**************************************************************************
 */
{
  /* define local variables  */

  PGSt_SMF_status retval = PGS_S_SUCCESS;
  PGSt_SMF_status status = PGS_S_SUCCESS;
  char asciiUTC[TIMECODEASIZE] = "";
  PGSt_double jdUTC[2];

  char msgbuf[MAXMSGLEN];
  static char filefunc[] = __FILE__ ", GEO_get_utcpole_metadata";
  PGSt_integer version = 1;
  PGSt_IO_Gen_FileHandle *filehandle;
  PGSt_double jdtable;

  /* Checking for NULL input  */

  if (ECS_metadata == NULL || utcpole_metadata == NULL) {
    sprintf(msgbuf, "ECS_metadata: %p utcpole_metadata: %p", 
	    (void*)ECS_metadata, (void*)utcpole_metadata);
    modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
    return MODIS_E_BAD_INPUT_ARG;
  }


  /* Call PGS_IO_Gen_Open() to open the file utcpole.dat */

  if (PGS_IO_Gen_Open(PGSd_UTCPOLE, PGSd_IO_Gen_Read, &filehandle, version)
      != PGS_S_SUCCESS) {
    sprintf(msgbuf, "PGS_IO_Gen_Open(%d)", (int)PGSd_UTCPOLE);
    modsmf(MODIS_E_GEO, msgbuf, filefunc);
    return MODIS_E_GEO;
  }


  /* Call fgets to read in the first line of the utcpole file. */

  if ( fgets(utcpole_metadata->header, sizeof(utcpole_metadata->header), filehandle) 
       != utcpole_metadata->header) {
    if (ferror(filehandle)) {
      PGS_SMF_SetUNIXMsg(errno, "fgets(utcpole.dat)", filefunc);
      retval = PGS_E_UNIX;
    }
    else {
      modsmf(MODIS_E_PREMATURE_EOF, "utcpole.dat", filefunc);
      retval = MODIS_E_PREMATURE_EOF;
    }
  }

  /* Call PGS_IO_Gen_Close() to close the utcpole file. */
  status = PGS_IO_Gen_Close(filehandle);
    
  /* Fill in asciiUTC */
  sprintf(asciiUTC, "%.10sT%.15s", ECS_metadata->rangebeginningdate,
	  ECS_metadata->rangebeginningtime);
    
  /* Convert asciiUTC into jdUTC. If success, load the elements of the 
     utcpole.polar_motion array with the x and y values for the polar wander, 
     and the UT1-UTC difference, in that order, corresponding to Julian date 
     jdUTC[0]+jdUTC[1]. */
  status = PGS_TD_UTCtoUTCjd(asciiUTC, jdUTC);
  if (status != PGS_S_SUCCESS) {
    sprintf(msgbuf, "PGS_TD_UTCtoUTCjd(%s)", asciiUTC);
    modsmf(MODIS_E_GEO, msgbuf, filefunc);
    return MODIS_E_GEO;
  }
  else {
    status = PGS_CSC_UTC_UT1Pole(jdUTC, &(utcpole_metadata->polar_motion[0]),
				 &(utcpole_metadata->polar_motion[1]), 
				 &(utcpole_metadata->polar_motion[2]), 
				 &jdtable);
    if (status != PGSCSC_W_PREDICTED_UT1 && status != PGS_S_SUCCESS) {
      sprintf(msgbuf, "PGS_CSC_UTC_UT1Pole(%s)", asciiUTC);
      modsmf(MODIS_E_GEO, msgbuf, filefunc);
      return MODIS_E_GEO;
    }
  }

  return retval;
}

