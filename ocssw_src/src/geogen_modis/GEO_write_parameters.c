#include <stddef.h>
#include "PGS_MODIS_35251.h"
#include "GEO_product.h"
#include "GEO_output.h"

int GEO_write_parameters(
	MODFILE			* const geo_file,
	GEO_param_struct const	* parameter
	)
/*!C****************************************************************************
!Description:   
	Routine in Output group of the Level-1A geolocation software to write
	geolocation processing parameters to the output product. These
	parameters, which are read from the geolocation parameter file, are
	output using MAPI calls.

		
!Input Parameters:
	geo_file	the MAPI structure for the product
	parameter	geolocation processing parameters structure

!Output Parameters:
                None

Return values:
	SUCCESS if all geolocation processing parameters written to geo_file
	FAIL otherwise.

Externally Defined:
	BAND_NUMBER             "GEO_product.h"
	BAND_POSITION           "GEO_product.h"
	DATATYPELENMAX          "mapi.h"
	DETECTOR_SPACE          "GEO_product.h"
	DETECTOR_OFFSETS        "GEO_product.h"
	FAIL                    "GEO_basic.h"
	FOCAL_LENGTH            "GEO_product.h"
	I16                     "mapi.h"
	MAPIOK                  "mapi.h"
	MAX_BAND_NUMBER         "GEO_geo.h"
	MODIS_E_BAD_INPUT_ARG   "PGS_MODIS_35251.h"
	MODIS_E_GEO             "PGS_MODIS_35251.h"
	MUNITS                  "mapi.h"
	NUMBANDS                "GEO_product.h"
	NUM_SAMPLES             "GEO_product.h"
	PARM_GRP                "GEO_product.h"
	R64                     "mapi.h"
	SCANDIM                 "GEO_product.h"
	SUCCESS                 "GEO_basic.h"
	T_OFFSET                "GEO_product.h"
	TRACKDIM                "GEO_product.h"
	TXT                     "mapi.h"
	UI16                    "mapi.h"
	UI32                    "mapi.h"

Called by:
	GEO_write_granule_metadata()
	
Routines Called:
	createMODISarray        "mapi.h"
	modsmf                  "smfio.h"
	putMODISarinfo          "mapi.h"
	putMODISarray           "mapi.h"
	putMODISdimname         "mapi.h"
	putMODISfileinfo        "mapi.h"

!Revision History:
 * $Log: GEO_write_parameters.c,v $
 * Revision 6.2  2011/02/14 21:11:40  kuyper
 * Corrected const-qualification of *parameter.
 *
 * Revision 6.1  2009/05/28 23:09:04  kuyper
 * Added BAND_NUMBER file attribute.
 * Corrected prolog.
 * Made error messages more informative.
 *
 * James Kuyper Jr.	James.R.Kuyper@nasa.gov
 *
 * Revision 4.1  2003/02/21 23:08:11  kuyper
 * Corrected to use void* pointers with %p format code.
 *
 * Revision 3.1  2002/06/13 22:57:47  kuyper
 * Removed unnecessary NCSA acknowledgement.
 *
 * Revision 2.5  1998/11/16 19:37:46  kuyper
 * Removed scaling of band_position to meters.
 *
 * Revision 2.4  1998/02/08  22:30:28  jjb
 * Merged from V2.0 DAAC Delivery.
 *
   1997/12/30   Liqun Ma
   Added NCSA acknowledgement

 * Revision 2.3  1997/11/03  20:20:08  jjb
 * Corrected 'units' strings for band_position and T_offset parameters.
 *
 * Revision 2.2  1997/10/27  16:08:25  ding@ltpmail.gsfc.nasa.gov
 * Added setting of parameter SDS dimension names.
 *
 * Revision 2.1  1997/10/21  18:16:22  kuyper
 * Returned from ClearCase
 *
 * Revision /main/GEO_V2_DEV/1 1997/10/03 kuyper
 * Removed mirror model parameters, poly_coeff, and coordinate transformation
 *   matrices. Added num_samples parameter.
 * Expanded INVALID_INPUT_ARG message.
 * Added scandim and trackdim attributes for detector_offsets SDS.
 *
 * Revision 1.10  1997/07/21  16:24:34  kuyper
 * Baselined Version 1
 *
 * Revision 1.11  1997/04/10  19:36:40  fhliang
 * Initial revision of SDST re-delivery of GEO_write_parameters.c.
 *
 * Revision 1.10  1997/04/09  14:23:25  kuyper
 * Added T_tel2inst, per filespec.
 *
 * Revision 1.9  1997/03/26  18:18:48  fhliang
 * Initial revision of SDST delivery of GEO_write_parameters.c
 *
		Revision 1.8  1997/02/13 20:03:42  kuyper
		Changed POLY_M1,POLY_M2 to fixed length of 3,
		  to match filespec.

		Revision 1.7  1996/12/26 23:07:24  kuyper
		Used new local array, save_band[], instead of global band_position,
		to save temporary band_positions.

		Revision 1.6  1996/12/06 22:29:01  kuyper
		Changed to exit on first failure in loop.

		Revision 1.5  1996/12/06 22:09:43  kuyper
		Corrected T_mirr2inst line.

		Revision 1.4  1996/12/06 21:19:57  kuyper
		Get parameters from structure passed by reference, rather than globals.
		Attach 'units' attribute to some of the arrays.
		Changed list of parameters.
		Collapse repeated code to a single loop, using an array of structures to
		control the loop.
		Validate input arguments.

		Revision 1.3  1996/07/24 21:59:19  kuyper
		Standardized order of #include files.
		Declared arguments const.
		Inserted required '!'s in comments.
		Removed ret_val.
		Converted constants to double, to skip implied conversion.
		Made implicit casts explicit.

		Revision 1.2  1996/07/18 21:33:50  kuyper
		Included self-checking header file.
		Replaced extern declarations with corresponding header files.
		Added needed header file.
		James Kuyper Jr. (kuyper@ltpmail.gsfc.nasa.gov)

		10/10/95
		Tracey W. Holmes
		Added debug option.

		9/27/95
		Frederick S. Patt
		Corrected declaration of N_samp; corrected indexing of
		  band offset arrays

		7/3/95
		Tracey W. Holmes (holmes@modis-xl.gsfc.nasa.gov)
		Added SDP error messages.

		6/21/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Finished coding.

!Team-unique Header:
	This software is developed by the MODIS Science Data Support Team for
	the National Aeronautics and Space Administration, Goddard Space Flight
	Center, under contract NAS5-32373.

!END*************************************************************************
*/
{

  /* Local variable declarations */

  #define BANDS (long)(MAX_BAND_NUMBER+1)
    struct {
      char *name;
      long rank;
      long dimensions[2];
      char data_type[DATATYPELENMAX];
      size_t data;
      char *units;
    } SDS_object[] = {
	{FOCAL_LENGTH,     1L, {BANDS},	R64,
	  offsetof(GEO_param_struct, geometry_params.focal_length), "meters"},
	{BAND_POSITION,    1L, {BANDS},	R64,
	  offsetof(GEO_param_struct, geometry_params.band_position),
	  "frame sample time"},
	{DETECTOR_SPACE,   1L, {BANDS},	R64,
	  offsetof(GEO_param_struct, geometry_params.det_space), "meters"},
	{DETECTOR_OFFSETS, 2L, {BANDS, 2L},	R64,
	  offsetof(GEO_param_struct, geometry_params.det_position), "meters"},
	{T_OFFSET,	       1L, {BANDS},	R64,
	  offsetof(GEO_param_struct, geometry_params.F_offset), 
	  "frame sample time"},
	{NUM_SAMPLES,      1L, {BANDS},     UI16,
	  offsetof(GEO_param_struct, geometry_params.N_samp), NULL}
    };
    long start[] = {0L, 0L, 0L};
    int obj;
    int retval = SUCCESS;
    uint32 dim_index[] = {0,1};
    char msgbuf[PGS_SMF_MAX_MSGBUF_SIZE];
    char filefunc[] = __FILE__ ", GEO_write_parameters";

   
    if(parameter==NULL || geo_file==NULL)
    {
	sprintf(msgbuf, "parameter=%p, geo_file=%p",
	    (void*)parameter, (void*)geo_file);
	modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);

	return FAIL;
    }

    /* Begin output of parameters */
    for(obj=0; obj<(int)(sizeof(SDS_object)/sizeof(SDS_object[0]))
      && retval==SUCCESS; obj++)
    {
      if(createMODISarray(geo_file, SDS_object[obj].name, PARM_GRP,
	  SDS_object[obj].data_type, SDS_object[obj].rank,
	  SDS_object[obj].dimensions) !=MAPIOK)
      {
	retval = FAIL;
	sprintf(msgbuf, "createMODISarray(\"%s\", \"" PARM_GRP
	    "\", \"%s\", %ld)", SDS_object[obj].name, SDS_object[obj].data_type,
	    SDS_object[obj].rank);
	    
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
      }
      else if( putMODISdimname(geo_file, SDS_object[obj].name, PARM_GRP, 0,
	  NUMBANDS) != MAPIOK)
      {
	  retval = FAIL;
	  sprintf(msgbuf, "putMODISdimname(\"%s\", \"%s\", \"" PARM_GRP
	      "\", 0, \"" NUMBANDS "\")",
	      geo_file->filename, SDS_object[obj].name);
	  modsmf(MODIS_E_GEO, msgbuf, filefunc);
      }
      else if( putMODISarray(geo_file, SDS_object[obj].name, PARM_GRP, start,
	  SDS_object[obj].dimensions, (char *)parameter+SDS_object[obj].data)
	  != MAPIOK)
      {
	  retval = FAIL;
	  sprintf(msgbuf, "putMODISarray(\"%s\", \"%s\", \"" PARM_GRP "\")",
	      geo_file->filename, SDS_object[obj].name);
	  modsmf(MODIS_E_GEO, msgbuf, filefunc);
      }
      else if(SDS_object[obj].units!=NULL)
      {
	  if(putMODISarinfo(geo_file, SDS_object[obj].name, PARM_GRP, MUNITS,
	      TXT, (long)strlen(SDS_object[obj].units), SDS_object[obj].units)
	      != MAPIOK)
	  {
	      retval = FAIL;
	      sprintf(msgbuf, "putMODISarinfo(\"%s\", \"%s\", \"" PARM_GRP
		  "\", \"" MUNITS "\", \"" TXT "\", %ld, \"%s\")",
		  geo_file->filename, SDS_object[obj].name,
		  (long)strlen(SDS_object[obj].units), SDS_object[obj].units);
	      modsmf(MODIS_E_GEO, msgbuf, filefunc);
	  }
      }
    } /* for each SDS_object	*/

    if(putMODISfileinfo(geo_file, BAND_NUMBER, I16, 1,
	&parameter->geometry_params.band_number) != MAPIOK)
    {
	retval = FAIL;
	sprintf(msgbuf,
	    "putMODISfileinfo(\"%s\", \"" BAND_NUMBER "\", \"" I16 "\", 1, %d)",
	    geo_file->filename, (int)parameter->geometry_params.band_number);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
    }

    if(putMODISarinfo(geo_file, DETECTOR_OFFSETS, PARM_GRP, SCANDIM, UI32,
	1L, &dim_index[0]) != MAPIOK)
    {
	retval = FAIL;
	sprintf(msgbuf, "putMODISarinfo(\"%s\", \"" PARM_GRP "\", \"" SCANDIM
	    "\", \"" UI32 "\", 1, %lu)", geo_file->filename,
	    (unsigned long)dim_index[0]);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
    }

    if(putMODISarinfo(geo_file, DETECTOR_OFFSETS, PARM_GRP, TRACKDIM, UI32,
	1L, &dim_index[1]) != MAPIOK)
    {
	retval = FAIL;
	sprintf(msgbuf, "putMODISarinfo(\"%s\", \"" PARM_GRP "\", \"" TRACKDIM
	    "\", \"" UI32 "\", 1, %lu)",
	    geo_file->filename, (unsigned long)dim_index[1]);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);
    }

    return retval;
}
