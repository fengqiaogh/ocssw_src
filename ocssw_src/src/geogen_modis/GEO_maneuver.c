#include "PGS_TD.h"
#include "GEO_main_func.h"
#include "PGS_MODIS_35251.h"

#define MANEUVER_LUN            600002

PGSt_SMF_status	GEO_in_maneuver(
	PGSt_double		tai1993
)
/*!C****************************************************************************

!Description:
	GEO_in_maneuver checks to determine whether the spacecraft was in a
	listed maneuver at the specified time.

!Input Parameters:
	tai1993		Specified time, in TAI seconds since 1993-01-01T00:00:00

!Output Parameters:
	None

Return Values:
	MODIS_E_GEO		A subroutine failed
	MODIS_E_BAD_INPUT_ARG	The maneuver list is incorrectly formatted.
	MODIS_N_GEO_MANEUVER	If the spacecraft was in a listed maneuver
	PGS_S_SUCCESS		Otherwise.

Externally Defined:
	MODIS_E_GEO		"PGS_MODIS_35251.h"
	MODIS_E_BAD_INPUT	"PGS_MODIS_35251.h"
	MODIS_N_GEO_MANEUVER	"PGS_MODIS_35251.h"
	PGS_S_SUCCESS		"PGS_SMF.h"

Called by:
	GEO_write_scan_metadata	"GEO_output.h"

Routines Called:
	modsmf			"smfio.h"
	PGS_PC_GetConfigData	"PGS_PC.h"
	PGS_TD_UTCtoTAI		"PGS_TD.h"

!Revision History:
$Log: GEO_maneuver.c,v $
Revision 6.3  2010/12/16 22:34:46  kuyper
Corrected parsing of maneuver list configuration parameter.

Revision 6.2  2010/06/24 19:53:25  kuyper
Corrected to treat a PGSPC_W_NO_CONFIG_FOR_ID return from
  PGS_PC_GetConfigData() as a successful return that requires no handling.

Revision 6.1  2010/05/26 18:32:26  kuyper
Changed MANEUVER_LUN from a PCF file reference to a PCF configuration
  parameter.
Changed to read the start and stop times of the maneuvers from that
  parameter rather than from that file.

James Kuyper Jr.	James.R.Kuyper@NASA.gov

Revision 5.8  2006/12/26 22:29:44  kuyper
Corrected to check against last maneuver in list.

Revision 5.7  2005/01/07 21:51:29  vlin
LogStatus message corrected.

Revision 5.6  2004/11/17 20:50:14  kuyper
Corrected handling of maneuver times that are too far in future for the
  leapsec.dat file.

Revision 5.5  2004/11/15 20:42:11  kuyper
Added check for null pointer arguments.

Revision 5.4  2004/11/05 02:13:29  kuyper
Corrected to initialize retval.

Revision 5.3  2004/11/05 00:49:59  kuyper
Corrected range of MANEUVER_TAG comparison.

Revision 5.2  2004/10/27 21:59:55  kuyper
Corrected to list maneuver_list as an output of GEO_read_maneuver_file.
Added validity checks on null pointer parameters.
Corrected use of maneuver length.
Corrected to list maneuver_list as an input/output parameter of
  GEO_in_maneuver.

Revision 5.1  2004/10/15 21:38:54  kuyper
Initial revision.

kuyper@saicmodis.com

James Kuyper Jr. kuyper@saicmodis.com

Requirements:	N/A

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits

Design Notes
!END**************************************************************************/
{
#define MAX_MANEUVERS \
    (PGSd_PC_VALUE_LENGTH_MAX/(2* sizeof "2010-05-25T12:57:04.000"))

    static PGSt_double maneuver_list[MAX_MANEUVERS][2];
    static int	count = -1;
    int man;
    char msgbuf[PGSd_PC_VALUE_LENGTH_MAX + 32];
    char filefunc[] = __FILE__ ", GEO_in_maneuver";

    if(count == -1)
    {
	char maneuver_string[PGSd_PC_VALUE_LENGTH_MAX];
	char const *pstart = maneuver_string;

	count = 0;
	switch(PGS_PC_GetConfigData(MANEUVER_LUN, maneuver_string))
	{
	case PGSPC_W_NO_CONFIG_FOR_ID:
	    return PGS_S_SUCCESS;	/* Empty maneuver list. */

	case PGS_S_SUCCESS:
	    break;

	default:
	    /* Unexpected return value, probably PGSPC_E_DATA_ACCESS_ERROR */
	    sprintf(msgbuf, "PGS_PC_GetConfigData(%ld)", (long)MANEUVER_LUN);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);

	    return MODIS_E_GEO;
	}

	while(*pstart && count < MAX_MANEUVERS)
	{
	    int val;
	    for(val=0; val<2; val++)
	    {
		char *pblank = strchr(pstart, ' ');

		if(!pblank)
		{
		    sprintf(msgbuf, "time string = \"%s\"", pstart);
		    modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);

		    return MODIS_E_BAD_INPUT_ARG;
		}

		*pblank = 'T';
		pblank = strchr(pblank+1, ' ');
		if(pblank)
		    *pblank = '\0';

		if(PGS_TD_UTCtoTAI((char*)pstart, maneuver_list[count] + val)
		    != PGS_S_SUCCESS)
		{
		    sprintf(msgbuf, "PGS_TD_UTCtoTAI(\"%s\")", pstart);
		    modsmf(MODIS_E_GEO, msgbuf, filefunc);

		    return MODIS_E_GEO;
		}

		if(pblank)
		    pstart = pblank+1;
		else
		{
		    if(val == 0)
		    {	/* Second time value is not correctly recorded. */
			sprintf(msgbuf, "maneuver_string = \"%s\"",
			    maneuver_string);
			modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);

			return MODIS_E_BAD_INPUT_ARG;
		    }
		    while(*pstart)
			pstart++;
		}
	    }
	    count++;
	}
    } 

    /* The count will almost always be 0, and when it is not, it will almost
     * always be 1, so the simplest search algorithm is also the best:
     */
    for(man = 0; man < count; man++)
	if(maneuver_list[man][0] < tai1993 && tai1993 < maneuver_list[man][1])
	    return MODIS_N_GEO_MANEUVER;
  
    return PGS_S_SUCCESS;
}

