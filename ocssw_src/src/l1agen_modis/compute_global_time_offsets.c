#include "L1A_prototype.h"
#include "PGS_TYPES.h"
#include "SC_scan.h"

void   compute_global_time_offsets  (PGSt_double  scan_rate)

/*
!C************************************************************************

!Description:  This function will compute the time offsets from the SD time
               for each sector within a scan based on the scan rate.

!Input Parameters:
               PGSt_double   scan_rate   ** The length of each scan in    **
                                         ** seconds                       **

!Output Parameters:
               None

Return Values: 
               None

Externally Defined:  
               PGSt_double                        (PGS_types.h)
               SC_SD_TIME_OFFSET_PERCENTAGE       (SC_scan.h)
               SC_SRCA_TIME_OFFSET_PERCENTAGE     (SC_scan.h)
               SC_BB_TIME_OFFSET_PERCENTAGE       (SC_scan.h)
               SC_SV_TIME_OFFSET_PERCENTAGE       (SC_scan.h)
               SC_EV_TIME_OFFSET_PERCENTAGE       (SC_scan.h)

Called By:     
               initialize_level1a

Routines Called:
               None

!Revision History:
               Revision 2.0  1997/08/26 15:03 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated Code. 
            
               Revision 1.0  1997/07/22 15:16 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Original design  

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes: 
               The CODE below is written in C language.
    
!END**********************************************************************
*/

 {
  /***************************************************************************/
  /*                                                                         */
  /*                      Define Global Variables                            */
  /*                                                                         */
  /***************************************************************************/

  extern PGSt_double  global_time_offset_array[5]; /* array to hold the time */
                                                  /* offsets for each sector */


  /***************************************************************************/
  /* Set global_time_offset_array[0] equal to                                */
  /*    SC_SD_TIME_OFFSET_PERCENTAGE * scan_rate                             */
  /*      (to calculate the SD time offset from the SD sector)               */
  /*                                                                         */
  /* Set global_time_offset_array[1] equal to                                */
  /*    SC_SRCA_TIME_OFFSET_PERCENTAGE * scan_rate                           */
  /*      (to calculate the SRCA time offset from the SD sector)             */
  /*                                                                         */
  /* Set global_time_offset_array[2] equal to                                */
  /*    SC_BB_TIME_OFFSET_PERCENTAGE * scan_rate                             */
  /*      (to calculate the BB time offset from the SD sector)               */
  /*                                                                         */
  /* Set global_time_offset_array[3] equal to                                */
  /*    SC_SV_TIME_OFFSET_PERCENTAGE * scan_rate                             */
  /*      (to calculate the SV time offset from the SD sector)               */
  /*                                                                         */
  /* Set global_time_offset_array[4] equal to                                */
  /*    SC_EV_TIME_OFFSET_PERCENTAGE * scan_rate                             */
  /*      (to calculate the EV time offset from the SD sector)               */
  /*                                                                         */
  /***************************************************************************/

  global_time_offset_array[0] = SC_SD_TIME_OFFSET_PERCENTAGE * scan_rate;

  global_time_offset_array[1] = SC_SRCA_TIME_OFFSET_PERCENTAGE * scan_rate;

  global_time_offset_array[2] = SC_BB_TIME_OFFSET_PERCENTAGE * scan_rate;

  global_time_offset_array[3] = SC_SV_TIME_OFFSET_PERCENTAGE * scan_rate; 

  global_time_offset_array[4] = SC_EV_TIME_OFFSET_PERCENTAGE * scan_rate;

 } /* End of routine compute_global_time_offsets */
