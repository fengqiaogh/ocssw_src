#include "L1A_prototype.h"
#include "PGS_IO.h"
#include "PD_pkt_data.h"
#include "SC_scan.h"


void   output_eng2_pkt2_to_scan (PGSt_IO_L0_Packet  *pkt,
                                 SC_SCAN_DATA_t     *L1A_scan)

/*
!C************************************************************************

!Description:  This routine extracts the engineering2, packet2 data contents
               and puts it in the scan structure.

!Input Parameters:
               PGSt_IO_L0_Packet   *pkt          ** The structure containing **
                                                 ** the current packet's     **
                                                 ** unpacked contents        **

!Output Parameters:
               None

!Input/Output Parameters:
               SC_SCAN_DATA_t      *L1A_scan     ** The MODIS scan structure **
                                                 ** currently being built    **

Return Values: 
               None

Externally Defined:
               PGSt_IO_L0_Packet               (PGS_IO.h)
               PD_E2P2_NUM_PV_GAINS            (PD_pkt_data.h)
               PD_E2P2_PV_GAINS_BYTE_OFFSET    (PD_pkt_data.h)
               SC_SCAN_DATA_t                  (SC_scan.h)

Called By:
               output_eng_data_to_scan

Routines Called:
               None

!Revision History:
               Revision 2.0  1997/09/09  13:19 EDT
               Timi Adelekan/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated Code.

               Revision 1.1  1997/08/27  10:00
               Tom Johnson  (johnson@ltpmail.gsfc.nasa.gov)
               Incorporate PDL walkthru comments

               Revision 1.0  1997/08/14  16:30 EDT
               Timi Adelekan/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Original design

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes: The CODE below was developed in C language.     

               This routine was designed totally under the assumption that
               the packet header has been previously validated and there
               is no need to test for error conditions.

!END************************************************************************
*/

 {
  /**************************************************************************/
  /*                                                                        */
  /*              Define and Initialize Local Variables                     */
  /*                                                                        */
  /**************************************************************************/

  int  i;      /* loop variable */


  /**************************************************************************/
  /*                                                                        */
  /* Use loop to move engineering 2, packet 2 contents to scan structure    */
  /* Suggested loop: FOR loop                                               */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /* FOR index equals 0 up to PD_E2P2_NUM_PV_GAINS                          */
  /*    set SC_SCAN_DATA_t.raw_pv_gains[index] equal to                     */
  /*                   pkt[PD_E2P2_PV_GAINS_BYTE_OFFSET + index]            */
  /* ENDFOR                                                                 */
  /*                                                                        */
  /**************************************************************************/


     for (i = 0; i < PD_E2P2_NUM_PV_GAINS; i++)
        L1A_scan->raw_pv_gains[i] = pkt[PD_E2P2_PV_GAINS_BYTE_OFFSET + i];

 } /* End of routine output_engr2_pkt2_to_scan */
