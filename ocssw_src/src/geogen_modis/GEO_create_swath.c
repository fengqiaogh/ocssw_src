#include "PGS_MODIS_35251.h"
#include "GEO_output.h"
#include "GEO_product.h"
#include "HdfEosDef.h"
#include "GEO_geo.h"
#include "smfio.h"

PGSt_SMF_status GEO_create_swath(
        int     const number_of_scans,
        int     const num_detectors,
        int32   const swfid
)

/*
***********************************************************************
!C
!Description:
        Routine in Output group of the Level-1A geolocation software to
        initialize the HDF-EOS objects in the output product for subsequent
        scan-by-scan output.  It calls HDF-EOS routines to initialize the
        geolocation metadata.

!Input Parameters:
        number_of_scans         the number of scans.
        num_detectors           the number of detectors per sample.
        swfid                   swath file id.

!Output Parameters:
        None

Return Values:
        MODIS_E_BAD_INPUT_ARG   If any input parameter is invalid
        MODIS_E_GEO             If any subroutine fails
        PGS_S_SUCCESS           Otherwise


Externally Defined:
        DETECTORS_1KM           "GEO_geo.h"
        DFNT_FLOAT32            "hdf.h"
        DFNT_INT16              "hdf.h"
        DFNT_UINT16             "hdf.h"
        DFNT_UINT8              "hdf.h"
        FAIL                    "hdf.h"
        GFLAGS                  "GEO_product.h"
        HDFE_NOMERGE            "HdfEosDef.h"
        HEIGHT                  "GEO_product.h"
        HEIGHT_OFFSET           "GEO_product.h"
	LAND_SEAMASK		"GEO_product.h"
        LATITUDE                "GEO_product.h"
        LONGITUDE               "GEO_product.h"
        MAX_FRAMES              "GEO_geo.h"
        MFRAMES                 "GEO_product.h"
        MFRAMES_2               "GEO_product.h"
        MFRAMES_4               "GEO_product.h"
        MODIS_E_BAD_INPUT_ARG   "PGS_MODIS_35251.h"
        MODIS_E_GEO             "PGS_MODIS_35251.h"
        NSCANS			"GEO_product.h"
        NSCANS_10               "GEO_product.h"
        NSCANS_20               "GEO_product.h"
        NSCANS_40               "GEO_product.h"
        PGS_S_SUCCESS           "PGS_SMF.h"
        RANGE                   "GEO_product.h"
        SCAN_OFFSET             "GEO_product.h"
        SEN_ZENITH              "GEO_product.h"
        SEN_AZIMUTH             "GEO_product.h"
        SOL_ZENITH              "GEO_product.h"
        SOL_AZIMUTH             "GEO_product.h"
        SUCCEED                 "hdf.h"
        SWATH_NAME              "GEO_product.h"
        TRACK_OFFSET            "GEO_product.h"
	WATER_PRESENT		"GEO_product.h"

Called by:
        GEO_initialize_product  "GEO_output.h"

Routines Called:
        modsmf                  "smfio.h"
        SWcreate                "HdfEosDef.h"
        SWdefdatafield          "HdfEosDef.h"
        SWdefdim                "HdfEosDef.h"
        SWdefdimmap             "HdfEosDef.h"
        SWdefgeofield           "HdfEosDef.h"
        SWdetach                "HdfEosDef.h"

!Revision History:
 * $Log: GEO_create_swath.c,v $
 * Revision 6.9  2011/02/11 21:20:34  kuyper
 * Added WaterPresent SDS.
 *
 * Revision 6.8  2010/06/23 17:55:37  kuyper
 * Removed attitude and ephemeris quality flags from swath.
 *
 * Revision 6.7  2010/06/18 20:46:40  kuyper
 * Corrected formation of message about SWdefdimmap() failure.
 *
 * Revision 6.6  2010/05/17 18:01:07  kuyper
 * Helped resolve Bug 2472 by adding attitude and ephemeris quality as
 *   datafields.
 *
 * James Kuyper Jr.		James.R.Kuyper@NASA.gov
 *
 * Revision 6.5  2009/06/12 18:51:34  kuyper
 * Corrected the comparison used to determine whether dimension maps need to
 *   be created.
 *
 * Revision 6.4  2009/05/31 21:06:42  kuyper
 * Corrected handling of the dimension list for the offsets Datafields.
 * Corrected setting of the 1km scan dimension.
 *
 * Revision 6.3  2009/05/31 01:47:12  ltan
 * Minor corrections.
 *
 * Revision 6.2  2009/05/31 01:28:36  ltan
 * Minor corrections.
 *
 * Revision 6.1  2009/05/29 15:44:41  ltan
 * Changed to match code by specifying use of macros for swath, dimension, and field names.
 * Changed length of MFRAMES dimension to MAX_FRAMES, and NSCANS_10 dimension to DETECTORS_1KM.
 * Added high resolution offsets, corresponding dimensions, and corresponding dimension maps.
 * Clarified where default zero-initialization is being relied upon.
 *
 * Revision 4.1  2003/02/21 20:25:16  kuyper
 * Simplified confusion notation used for dimlist fields.
 *
 * Revision 3.1  2002/06/13 22:47:26  kuyper
 * Removed unnecessary NCSA acknowledgment.
 *
 * Revision 2.4  1999/02/02 18:03:23  kuyper
 * Make an implicit conversion explicit.
 *
 * Revision 2.3  1998/02/08  22:20:07  jjb
 * Merged from V2.0 DAAC delivery.
 *
  1997/12/20   Liqun Ma
  Added NCSA acknowledgement

 * Revision 2.2  1997/10/24  13:13:59  kuyper
 * Included header file which declares GEO_create_swath().
 *
 * Revision 2.1  1997/10/21  18:16:22  kuyper
 * Returned from ClearCase
 *
 * Revision GEO_create_swath.c@@/main/GEO_V2_DEV/1
 * "version 2 development"
 *
# Revision 2.1  1997/07/07  18:36:25  fshaw
# made changes reported by 7/3/97 W/T
#
# Revision 2.0  1997/07/07  17:06:09  fshaw
# Initial version
#
!Requirements:
		PR03-I-4

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.
!END
*************************************************************************
*/
{
 /*  Declare local variables                         */

   int32 swathID;
   PGSt_SMF_status return_status = PGS_S_SUCCESS;
   int32 dim, fld;
   uint16 N_samp = num_detectors/DETECTORS_1KM;
   int32 dimensions =4;
   int32 startfield = 0;
   char  msgbuf[256];

   static struct{
      char *fieldname;
      int32 dim;
      int32 offset;
      int32 increment;
   }Swath_Dimname[] = 
   { 
     {NSCANS_10, 0, 0, 0},
     {MFRAMES, MAX_FRAMES, 0, 0},
     {NULL,0,0,0},
     {NULL,0,0,0},
   };

   static struct{
      char *fieldname;
      char *dimlist;
      int32 datatype;
   }Swath_Geofields[] =
   { {LATITUDE,		NSCANS_10","MFRAMES,	DFNT_FLOAT32},
     {LONGITUDE,	NSCANS_10","MFRAMES,	DFNT_FLOAT32},
   };

   static char offsets_dimlist[] = NSCANS_20","MFRAMES_2;
   static struct{
      char *fieldname;
      char *dimlist;
      int32 datatype;
   }Swath_Datafields[] =
   { {SCAN_OFFSET,      offsets_dimlist,	DFNT_INT8},
     {TRACK_OFFSET,     offsets_dimlist,	DFNT_INT8},
     {HEIGHT_OFFSET,    offsets_dimlist,	DFNT_INT8},
     {HEIGHT,		NSCANS_10 "," MFRAMES,	DFNT_INT16},
     {SEN_ZENITH,	NSCANS_10 "," MFRAMES,	DFNT_INT16},
     {SEN_AZIMUTH,	NSCANS_10 "," MFRAMES,	DFNT_INT16},
     {RANGE,		NSCANS_10 "," MFRAMES,	DFNT_UINT16},
     {SOL_ZENITH,	NSCANS_10 "," MFRAMES,	DFNT_INT16},
     {SOL_AZIMUTH,	NSCANS_10 "," MFRAMES,	DFNT_INT16},
     {LAND_SEAMASK,	NSCANS_10 "," MFRAMES,	DFNT_UINT8},
     {WATER_PRESENT,	NSCANS_10 "," MFRAMES,  DFNT_UINT8},
     {GFLAGS,		NSCANS_10 "," MFRAMES,	DFNT_UINT8},
   };

   char filefunc[] = __FILE__ ",GEO_create_swath";

  if( number_of_scans == 0)
    return PGS_S_SUCCESS;

  if ( (number_of_scans < 0) || (num_detectors <= 0 )
   || (swfid == FAIL ))
     {
        sprintf(msgbuf, "\n number_of_scans = %d,\n"
           " num_detectors = %d,\n swfid = %ld\n",
            number_of_scans, num_detectors, (long)swfid);
        modsmf(MODIS_E_BAD_INPUT_ARG,  msgbuf,
         "GEO_create_swath, GEO_create_swath.c");
        return MODIS_E_BAD_INPUT_ARG;
      }

  if (N_samp == 1)
  {
    dimensions = 2;
    startfield = 3;
  }
  else if(N_samp == 2)
  {
    Swath_Dimname[2].fieldname = NSCANS_20;
    Swath_Dimname[2].dim = number_of_scans*20;
    Swath_Dimname[2].increment = 2;
    Swath_Dimname[3].fieldname = MFRAMES_2;
    Swath_Dimname[3].dim = MAX_FRAMES*2;
    Swath_Dimname[3].increment = 2;
  }
  else if(N_samp == 4)
  {
    Swath_Dimname[2].fieldname = NSCANS_40;
    Swath_Dimname[2].dim = number_of_scans*40;
    Swath_Dimname[2].offset = 1;
    Swath_Dimname[2].increment = 4;
    Swath_Dimname[3].fieldname = MFRAMES_4;
    Swath_Dimname[3].dim = MAX_FRAMES*4;
    Swath_Dimname[3].increment = 4;

    strncpy(offsets_dimlist, NSCANS_40","MFRAMES_4, sizeof offsets_dimlist);
  }
 
   /* Create the Swath    */

    if ((swathID = SWcreate(swfid, SWATH_NAME)) == FAIL){
       modsmf(MODIS_E_GEO,  
              "SWcreate(\""SWATH_NAME"\")", 
              filefunc);
        return MODIS_E_GEO;
    }

   /* Name the Swath dimensions   */
   
   Swath_Dimname[0].dim = number_of_scans * DETECTORS_1KM;

   for(dim=0; dim<dimensions; dim++ ){
      if (SWdefdim(swathID, 
                 Swath_Dimname[dim].fieldname, 
                 Swath_Dimname[dim].dim) != SUCCEED){
         /* call SDP SMF function to report error */
	 sprintf(msgbuf, "SWdefdim(\"%s\")", Swath_Dimname[dim].fieldname);
         modsmf(MODIS_E_GEO, msgbuf, filefunc);
         return_status = MODIS_E_GEO;
      } 
   }

   if (N_samp >1) {
     for(dim=0; dim < 2; dim++){
        if (SWdefdimmap(swathID, 
                      Swath_Dimname[dim].fieldname, 
                      Swath_Dimname[dim+2].fieldname, 
                      Swath_Dimname[dim+2].offset, 
                      Swath_Dimname[dim+2].increment) != SUCCEED){
         /* call SDP SMF function to report error */
	 sprintf(msgbuf, "SWdefdimmap(\"%s\",\"%s\",\"%ld\",\"%ld\")", 
                 Swath_Dimname[dim].fieldname, 
                 Swath_Dimname[dim+2].fieldname, 
                 (long)Swath_Dimname[dim+2].offset, 
                 (long)Swath_Dimname[dim+2].increment);
         modsmf(MODIS_E_GEO, msgbuf, filefunc);
         return_status = MODIS_E_GEO;
        }
     } 
   }

  /* Set up Geolocation Fields  */

   for(fld=0; fld<(int)(sizeof Swath_Geofields/sizeof Swath_Geofields[0]); fld++){
      if (SWdefgeofield(swathID,
                      Swath_Geofields[fld].fieldname,
                      Swath_Geofields[fld].dimlist,
                      Swath_Geofields[fld].datatype,
                      HDFE_NOMERGE) != SUCCEED){
         /* call SDP SMF function to report error */
	 sprintf(msgbuf, "SWdefgeofield(\"%s\")", Swath_Geofields[fld].fieldname);
         modsmf(MODIS_E_GEO, msgbuf, filefunc);
         return_status = MODIS_E_GEO;
      }
   }

   for(fld=startfield; fld<(int)(sizeof Swath_Datafields/sizeof Swath_Datafields[0]); fld++){
      if (SWdefdatafield(swathID,
                      Swath_Datafields[fld].fieldname,
                      Swath_Datafields[fld].dimlist,
                      Swath_Datafields[fld].datatype,
                      HDFE_NOMERGE) != SUCCEED){
         /* call SDP SMF function to report error */
	 sprintf(msgbuf, "SWdefdatafield(\"%s\")", Swath_Datafields[fld].fieldname);
         modsmf(MODIS_E_GEO, msgbuf,filefunc);
         return_status = MODIS_E_GEO;
      }
   }

   /* Detach swath   */

   if( SWdetach(swathID) != SUCCEED){
      /* call SDP SMF function to report error */
       modsmf(MODIS_E_GEO, 
              "SWdetach()",
              filefunc);
       return_status =  MODIS_E_GEO;
    }

   return(return_status);

}
