#include "L1A_prototype.h"
#include "SC_scan.h"
#include "PH_pkt_hdr.h"


void   initialize_pixel_qual_data (SC_PIXEL_QUALITY_DATA_t    *scan_pixel)

/*
!C************************************************************************

!Description:  This function fills the Pixel Quality Data structure (section
               3 of the MODIS Level 1A Data Product Format) with the 
               appropriate fill data (missing packet (1)).

!Input Parameters:
               None

!Output Parameters:
               SC_PIXEL_QUALITY_DATA_t  *scan_pixel  ** Pixel Quality Data 
                                                        structure        **

Return Values: 
               None

Externally Defined:  
               SC_PIXEL_QUALITY_DATA_t                    (SC_scan.h)
               SC_MISSING_PACKET                          (SC_scan.h)
               PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX         (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT     (PH_pkt_hdr.h)

Called By:
               initialize_scan

Routines Called:
               None

!Revision History:      
               Revision 2.1  1997/08/16  12:04 EDT 
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated CODE            

               Revision 2.0  1997/00/00  00:00
               Tom Johnson/GSC (johnson@ltpmail.gsfc.nasa.gov)
               Original design

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes: 
               The CODE below was developed in C language.

!END**********************************************************************
*/

 {
  /**************************************************************************/
  /*                                                                        */
  /*                Declare and Initialize Local Variables                  */
  /*                                                                        */
  /**************************************************************************/

  int    frame;     /*  Loop variable to count frames                       */


  /**************************************************************************/
  /*                                                                        */
  /* Set all elements of SC_PIXEL_QUALITY_DATA_t.SD_pix_qual to             */
  /*    SC_MISSING_PACKET                                                   */
  /*                                                                        */
  /* Set all elements of SC_PIXEL_QUALITY_DATA_t.SRCA_pix_qual to           */
  /*    SC_MISSING_PACKET                                                   */
  /*                                                                        */
  /* Set all elements of SC_PIXEL_QUALITY_DATA_t.BB_pix_qual to             */
  /*    SC_MISSING_PACKET                                                   */
  /*                                                                        */
  /* Set all elements of SC_PIXEL_QUALITY_DATA_t.SV_pix_qual to             */
  /*    SC_MISSING_PACKET                                                   */
  /*                                                                        */
  /**************************************************************************/

  for (frame = 0; frame < PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX; frame++)
    {
     scan_pixel->SD_pix_qual[frame][0] = SC_MISSING_PACKET;
     scan_pixel->SD_pix_qual[frame][1] = SC_MISSING_PACKET;

     scan_pixel->SRCA_pix_qual[frame][0] = SC_MISSING_PACKET;
     scan_pixel->SRCA_pix_qual[frame][1] = SC_MISSING_PACKET;

     scan_pixel->BB_pix_qual[frame][0] = SC_MISSING_PACKET;
     scan_pixel->BB_pix_qual[frame][1] = SC_MISSING_PACKET;

     scan_pixel->SV_pix_qual[frame][0] = SC_MISSING_PACKET;
     scan_pixel->SV_pix_qual[frame][1] = SC_MISSING_PACKET;
    }


  /**************************************************************************/
  /*                                                                        */
  /* Set all elements of SC_PIXEL_QUALITY_DATA_t.EV_pix_qual to             */
  /*    SC_MISSING_PACKET                                                   */
  /*                                                                        */
  /**************************************************************************/

  for (frame = 0; frame < PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT; frame++)
    {
     scan_pixel->EV_pix_qual[frame][0] = SC_MISSING_PACKET;
     scan_pixel->EV_pix_qual[frame][1] = SC_MISSING_PACKET;
    }


 } /* End of routine initialize_pixel_qual_data */
