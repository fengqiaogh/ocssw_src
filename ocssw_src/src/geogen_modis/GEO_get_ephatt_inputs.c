#include <errno.h>
#include <stdio.h>
#include "PGS_IO.h"
#include "smfio.h"
#include "GEO_earth.h"
#include "GEO_product.h"
#include "PGS_MODIS_35251.h"

/* Used internally by the EPH portion of the SDP Toolkit library. */
extern int byteswap(char *, int);

PGSt_SMF_status GEO_get_ephatt_inputs(
	PGSt_PC_Logical	file_logical,
	PGSt_integer	file_version,
	char		universal_references[][PGSd_UR_FIELD_SIZE]
){
/*!C****************************************************************************

!Description:   
	A routine for extracting from files in the SDP Toolkit's ephemeris or
	attitude file formats the names of the input files used to create them.

!Input Parameters:
	file_logical		Identifies the file type to be read.
	file_version		Identifies which one of those files to read.

!Output Parameters:
	universal_references	An array for storing the universal reference
				strings extracted from that file.

Return Values:
	MODIS_E_BAD_INPUT_ARG	If universal_references is null.
	MODIS_E_GEO		If an SDP Toolkit subroutine failed.
	PGS_E_UNIX		If a POSIX subroutine failed.
	PSG_S_SUCCESS		Otherwise

Externally Defined:
	errno			<errno.h>
	MODIS_E_BAD_INPUT_ARG	"PGS_MODIS_35251.h"
	MODIS_E_GEO		"PGS_MODIS_35251.h"
	PGS_E_UNIX		"PGS_SMF.h"
	PSG_S_SUCCESS	`	"PGS_SMF.h"
	PGSd_IO_Gen_Read	"PGS_IO_Gen.h"
	PGSd_UR_FIELD_SIZE	"PGS_EPH.h"

Called by:
	GEO_write_input_metadata()	"GEO_output.h"

Routines Called:
	PGS_IO_Gen_Close()		"PGS_IO.h"
	PGS_IO_Gen_Open()		"PGS_IO.h"
	byteswap()			"PGS_EPH_getEphemRecords.c"

!Revision History:
$Log: GEO_get_ephatt_inputs.c,v $
Revision 6.2  2010/06/29 19:47:10  kuyper
Corrected order of arguments to fread().
Added an endianness check.

Revision 6.1  2010/05/27 21:15:34  kuyper
Initial revision, MOD_PR03 6.0.3.

James Kuyper Jr.		James.R.Kuyper@nasa.gov

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits

Design Notes
	The ephemeris and attitude files were designed to have the same layout,
	at least with regard to all of the fields used by this program.
	The layout is defined by a C struct type. As such, it is inherently
	non-portable, because the C standard gives implementations a great
	deal of latitude for deciding how to lay out fields in a structure.
	However, the SDP Toolkit's own routines use this same approach, so if
	this rountine can't access those files correctly, neither will the
	Toolkit's routines.
!END**************************************************************************
*/
    PGSt_IO_Gen_FileHandle	*handle;
    PGSt_ephemHeader		fileHeader;
    PGSt_SMF_status		retval = PGS_S_SUCCESS;
    char			filefunc[] = __FILE__ ", GEO_get_ephatt_inputs";
    char			msgbuf[256];

    if(universal_references == NULL)
    {
	modsmf(MODIS_E_BAD_INPUT_ARG, "universal_references is null",
	    filefunc);
	
	return MODIS_E_BAD_INPUT_ARG;
    }

    if(PGS_IO_Gen_Open(file_logical, PGSd_IO_Gen_Read, &handle, file_version)
	!= PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "PGS_IO_Gen_Open(%ld, %ld)",
	    (long)file_logical, (long)file_version);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
	
	return MODIS_E_GEO;
    }

    errno = 0;
    if(fread(&fileHeader, sizeof fileHeader, 1, handle) != 1)
    {
	if(errno)
	    PGS_SMF_SetUNIXMsg(errno, "during fread(fileHeader)", filefunc);
	else
	    modsmf(MODIS_E_GEO, "fread(fileHeader)", filefunc);

	retval = PGS_E_UNIX;
    }
    else
    {
	/* Simplistic byte-order check, borrowed from PGS_EPH_GetEphMet() */
	if(fileHeader.nURs > 100)
	    byteswap((char*)&fileHeader.nURs, sizeof fileHeader.nURs);

	if(fileHeader.nURs > MAX_EA_INPUTS)
	    fileHeader.nURs = MAX_EA_INPUTS;

	if(fread(universal_references, PGSd_UR_FIELD_SIZE, fileHeader.nURs,
	    handle) != fileHeader.nURs)
	{
	    if(errno)
		PGS_SMF_SetUNIXMsg(errno, "during fread(universal_references",
		    filefunc);
	    else
		modsmf(MODIS_E_GEO, "fread(universal_references)", filefunc);

	    retval = PGS_E_UNIX;
	}
	else if(fileHeader.nURs < MAX_EA_INPUTS)
	    /*Terminate the array with an empty string. */
	    universal_references[fileHeader.nURs][0] = '\0';
    }

    if(PGS_IO_Gen_Close(handle) != PGS_S_SUCCESS)
	modsmf(MODIS_E_GEO, "PGE_IO_Gen_Close()", filefunc);

    return retval;
}

