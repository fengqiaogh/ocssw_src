#include "PGS_EPH.h"
#include "GEO_earth.h"
#include "GEO_main.h"
#include "GEO_output.h"
#include "GEO_product.h"
#include "PGS_MODIS_35251.h"

PGSt_SMF_status GEO_write_input_metadata(
	MODFILE			* const geo_file,
	int			const l1a_version,
	PGSt_integer		const param_version,
	pointer_metadata_struct	* const pointer_metadata,
	PGSt_tag		sc_tag
)

/*******************************************************************************
!C
!Description:   
	Routine in Output group of the Level-1A geolocation software to write
	the input pointer metadata.

!Input Parameters:
	geo_file		MAPI structure for geolocation file
	l1a_version		Version number of input/output files.
	param_version		Version number of the parameter file.
	sc_taq			Spacecraft tag

!Output Parameters:
	None.

!Input/Output Parameters:
	pointer_metadata	the universal references of geolocation inputs

Return Parameters:
	MODIS_E_BAD_INPUT_ARG	If pointer_metadata is null.
	MODIS_E_GEO		If any subroutine failed
	PGS_S_SUCCESS		Otherwise

Externally Defined:
	ATTIT_INPUT			"GEO_product.h"
	EPHEM_INPUT			"GEO_product.h"
	l1aID				"L1a_data.h"
	MAPIOK				"mapi.h"
	MAX_EA_FILES			"GEO_product.h"
	MAX_EA_INPUTS			"GEO_product.h"
	MAX_NC_NAME			"netcdf.h"
	MODIS_E_BAD_INPUT_ARG		"PGS_MODIS_35251.h"
	MODIS_E_GEO			"PGS_MODIS_35251.h"
	MODIS_N_GEO_EXCESS_FILES	"PGS_MOIDS_35251.h"
	MODIS_W_GEO_MISSING_MET		"PGS_MOIDS_35251.h"
	PARAM				"GEO_main.h"
	PGSd_EOS_AM			"PGS_TD.h"
	PGSd_SC_ATTIT_DATA		"PGS_EPH.h"
	PGSd_SC_EPHEM_DATA		"PGS_EPH.h"
	PGSd_UR_FIELD_SIZE		"PGS_EPH.h"
	PGS_S_SUCCESS			"PGS_SMF.h"

Called by:
	GEO_write_granule_metadata()

Routines Called:
	GEO_get_ephatt_inputs		"GEO_earth.h"
	modsmf				"smfio.h"
	PGS_PC_GetNumberOfFiles		"PGS_PC.h"
	PGS_PC_GetUniversalRef		"PGS_PC.h"
	putMODIStable			"mapi.h"

!Revision History:
$Log: GEO_write_input_metadata.c,v $
Revision 6.2  2012/06/25 15:45:02  kuyper
Corrected to insert null characters only where needed.

Revision 6.1  2010/05/28 17:27:59  kuyper
Helped resolve Bug 2471 by changing to get epheris and atttitude input files
  from a PCF entry rather than from the corresponding *.met files.
Moved MAX_EA_FILES, MAX_EA_INPUTS to GEO_product.h.
Changed to no longer assume that input_strings are nul-terminated, and to
  ensure that output strings are nul-delimited.

James Kuyper Jr.		James.R.Kuyper@NASA.gov

Revision 5.1  2008/09/21 19:00:19  kuyper
Added sc_tag argument, in order to enable supression of the error message
  for missing attitude files for Terra runs.

Revision 4.8  2004/04/09 21:09:36  kuyper
Corrected typo.

Revision 4.7  2003/12/30 22:03:02  kuyper
Changed to initialize input_strings.

Revision 4.6  2003/11/12 20:35:39  vlin
Function returns PGSt_SMF_status, check for null inputs.

Revision 4.5  2003/10/23 21:40:05  kuyper
Increased value of MAX_EA_INPUTS. (sufficiently?)

Revision 4.4  2003/08/28 16:02:50  kuyper
Corrected prolog.

Revision 4.3  2003/08/26 21:18:10  kuyper
Minor re-write, to silence warnings from GNU compiler.

Revision 4.2  2003/08/20 17:54:12  kuyper
Corrected direction of copies when merging input pointer strings.

Revision 4.1  2003/08/12 22:11:22  kuyper
Initial revision.

James Kuyper Jr. (kuyper@saicmodis.com)

Requirements:
	None

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits
	None

Design Notes
!END**************************************************************************
*/
{
    typedef char ephem_inputs_t[MAX_EA_INPUTS*PGSd_UR_FIELD_SIZE];
    static  ephem_inputs_t *ephem_inputs;
    typedef char attit_inputs_t[MAX_EA_INPUTS*PGSd_UR_FIELD_SIZE];
    static attit_inputs_t *attit_inputs;
    if(!ephem_inputs) {
        ephem_inputs = (ephem_inputs_t*) calloc(MAX_EA_FILES, sizeof(ephem_inputs_t));
        attit_inputs = (attit_inputs_t*) calloc(MAX_EA_FILES, sizeof(attit_inputs_t));
    }
    static PGSt_integer num_ephem, num_attit, num_ptrs;
    static ptrdiff_t ephem_length[MAX_EA_FILES], attit_length[MAX_EA_FILES];

    char input_strings[MAX_EA_INPUTS][PGSd_UR_FIELD_SIZE];

    PGSt_SMF_status retval = PGS_S_SUCCESS;
    PGSt_integer ver;
    /* temp_ver is needed, because PGS_PC_GetReference() updates the number
     * pointed at by the pointer passed to it.
     */
    PGSt_integer temp_ver;
    char attribname[MAX_NC_NAME];
    int str;
    char *p, *q;
    char msgbuf[PGSd_UR_FIELD_SIZE*2];
    char filefunc[] = __FILE__ ", GEO_write_input_metadata";

    if (geo_file == NULL || pointer_metadata == NULL)
    {
        sprintf(msgbuf,"geo_file: %p, pointer_metadata: %p", (void*)geo_file,
                (void*)pointer_metadata);
	modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
	return MODIS_E_BAD_INPUT_ARG;
    }

    temp_ver = (PGSt_integer)l1a_version;
    if(PGS_PC_GetUniversalRef(l1aID, &temp_ver,
	pointer_metadata->inputpointer[0]) != PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "PGS_PC_GetUniversalRef(%ld,%d)",
	    (long)l1aID, l1a_version);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);

	retval = MODIS_E_GEO;
    }

    if(num_ptrs == 0)
    {	/* Haven't yet initialized the constant inputs. */

	temp_ver = param_version;
	if(PGS_PC_GetUniversalRef(PARAM, &temp_ver,
	    pointer_metadata->inputpointer[1]) != PGS_S_SUCCESS)
	{
	    sprintf(msgbuf, "PGS_PC_GetUniversalRef(%ld,%d)",
		(long)PARAM, param_version);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);

	    retval = MODIS_E_GEO;
	}

	num_ptrs = 2;


	/* Check ephemeris files. */
	if(PGS_PC_GetNumberOfFiles(PGSd_SC_EPHEM_DATA, &num_ephem)
	    != PGS_S_SUCCESS)
	{
	    sprintf(msgbuf, "PGS_PC_GetNumberOfFiles(%ld)",
		(long)PGSd_SC_EPHEM_DATA);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);

	    retval = MODIS_E_GEO;
	}

	if(num_ephem > MAX_EA_FILES)
	{
	    /* It is perfectly legal to stage more files than needed. This
	     * shouldn't even be a warning.
	     */
	    sprintf(msgbuf, "%d ephemeris", MAX_EA_FILES);
	    modsmf(MODIS_N_GEO_EXCESS_FILES, msgbuf, filefunc);

	    num_ephem = MAX_EA_FILES;
	}

	for(ver=0; ver<num_ephem; ver++)
	{
	    temp_ver = ver+1;
	    if(PGS_PC_GetUniversalRef(PGSd_SC_EPHEM_DATA, &temp_ver,
		pointer_metadata->inputpointer[num_ptrs]) != PGS_S_SUCCESS)
	    {
		sprintf(msgbuf, "PGS_PC_GetUniversalRef(%ld,%d)",
		    (long)PGSd_SC_EPHEM_DATA, ver+1);
		modsmf(MODIS_E_GEO, msgbuf, filefunc);

		retval = MODIS_E_GEO;
	    }
	    else
		++num_ptrs;


	    if(GEO_get_ephatt_inputs(PGSd_SC_EPHEM_DATA, (PGSt_integer) (ver+1),
		input_strings) != PGS_S_SUCCESS)
	    {
		sprintf(msgbuf, "GEO_get_ephatt_inputs(%ld, %d)",
		  (long)PGSd_SC_EPHEM_DATA, ver+1);
		modsmf(MODIS_E_GEO, msgbuf, filefunc);

		return MODIS_E_GEO;
	    }

	    /* Catenate input strings into a single null-delimited string */
	    p = ephem_inputs[ver];
	    for(str=0; str<MAX_EA_INPUTS && input_strings[str][0]; str++)
	    {
		/* Append input_strings[str] to ephem_inputs[ver]	*/
		q = input_strings[str];

		while(q < (char*)(input_strings + str + 1) &&
		    (*p++ = *q++));
		if(p[-1])
		    *p++ = '\0';	/* nul terminate */
	    }

	    ephem_length[ver] = p - ephem_inputs[ver];
	}

	/* Check Attitude Files */
	if(PGS_PC_GetNumberOfFiles(PGSd_SC_ATTIT_DATA, &num_attit)
	    != PGS_S_SUCCESS && sc_tag!=PGSd_EOS_AM)
	{
	    /* Attitude files are not always necessary for Terra runs; LUN
	     * 600280 might be set to "MODIS Packet".
	     */ 
	    sprintf(msgbuf, "PGS_PC_GetNumberOfFiles(%ld)",
		(long)PGSd_SC_ATTIT_DATA);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);

	    retval = MODIS_E_GEO;
	}

	if(num_attit > MAX_EA_FILES)
	{
	    /* It is perfectly legal to stage more files than needed. This
	     * shouldn't even be a warning.
	     */
	    sprintf(msgbuf, "%d attitude", MAX_EA_FILES);
	    modsmf(MODIS_N_GEO_EXCESS_FILES, msgbuf, filefunc);

	    num_attit = MAX_EA_FILES;
	}

	for(ver=0; ver<num_attit; ver++)
	{
	    temp_ver = ver+1;
	    if(PGS_PC_GetUniversalRef(PGSd_SC_ATTIT_DATA, &temp_ver,
		pointer_metadata->inputpointer[num_ptrs]) != PGS_S_SUCCESS)
	    {
		sprintf(msgbuf, "PGS_PC_GetUniversalRef(%ld,%d)",
		    (long)PGSd_SC_ATTIT_DATA, ver+1);
		modsmf(MODIS_E_GEO, msgbuf, filefunc);

		retval = MODIS_E_GEO;
	    }
	    else
		++num_ptrs;

	    if(GEO_get_ephatt_inputs(PGSd_SC_ATTIT_DATA, (PGSt_integer)(ver+1),
		input_strings) != PGS_S_SUCCESS)
	    {
		sprintf(msgbuf, "GEO_get_ephatt_inputs(%ld, %d)",
		  (long)PGSd_SC_EPHEM_DATA, ver+1);
		modsmf(MODIS_E_GEO, msgbuf, filefunc);

		return MODIS_E_GEO;
	    }

	    /* Catenate input strings into a single null-delimited string */
	    p = attit_inputs[ver];
	    for(str=0; str<MAX_EA_INPUTS; str++)
	    {
		/* Append input_strings[str] to attit_inputs[ver]	*/
		q = input_strings[str];

		while(q < (char*)(input_strings + str + 1) &&
		    (*p++ = *q++));
		if(p[-1])
		    *p++ = '\0';	/* nul terminate */
	    }

	    attit_length[ver] = p - attit_inputs[ver];
	}
    }

    /* Write out the ephemeris/attitude input strings */
    for(ver=0; ver<num_ephem; ver++)
    {
	sprintf(attribname, "%s.%d", EPHEM_INPUT, ver+1);

	if(putMODISfileinfo(geo_file, attribname, TXT, ephem_length[ver],
	    ephem_inputs[ver]) != MAPIOK)
	{
	    sprintf(msgbuf, "putMODISfileinfo(\"%s\", \"%s\")",
		geo_file->filename, attribname);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);

	    retval = MODIS_E_GEO;
	}
    }

    for(ver=0; ver<num_attit; ver++)
    {
	sprintf(attribname, "%s.%d", ATTIT_INPUT, ver+1);

	if(putMODISfileinfo(geo_file, attribname, TXT, attit_length[ver],
	    attit_inputs[ver]) != MAPIOK)
	{
	    sprintf(msgbuf, "putMODISfileinfo(\"%s\", \"%s\")",
		geo_file->filename, attribname);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);

	    retval = MODIS_E_GEO;
	}
    }

    return retval;
}

