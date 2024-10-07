#include "L1A_prototype.h"
#include "PGS_IO.h"
#include "PD_pkt_data.h"
#include "SC_scan.h"
#include "hdfi.h"


void   output_eng2_pkt1_to_scan (PGSt_IO_L0_Packet  *pkt,
                                 SC_SCAN_DATA_t     *L1A_scan)

/*
!C************************************************************************

!Description:  This routine extracts the engineering2, packet1 data contents
               and puts it in the scan structure.

!Input Parameters:
               PGSt_IO_L0_Packet    *pkt         ** The structure containing **
                                                 ** the current packet's     **
                                                 ** unpacked contents        **

!Output Parameters:
               None

!Input/Output Parameters:
               SC_SCAN_DATA_t       *L1A_scan    ** The MODIS scan structure **
                                                 ** currently being built    **

Return Values: 
               None

Externally Defined:
               PGSt_IO_L0_Packet             (PGS_IO.h)
               PD_E2P1_NUM_HK_TELEM_BYTES    (PD_pkt_data.h)
               PD_E2P1_CURR_HK_BYTE_OFFSET   (PD_pkt_data.h)
               PD_E2P1_PRIOR_HK_BYTE_OFFSET  (PD_pkt_data.h)
               PD_E2P1_NUM_SC_ANCIL_WORDS    (PD_pkt_data.h)
               PD_E2P1_SC_ANCIL_BYTE_OFFSET  (PD_pkt_data.h)
               PD_E2P1_NUM_PARAM_BYTES       (PD_pkt_data.h)
               PD_E2P1_PARAM_BYTE_OFFSET     (PD_pkt_data.h)
               SC_SCAN_DATA_t                (SC_scan.h)

Called By:
               output_eng_data_to_scan

Routines Called:
               None

!Revision History:
               Revision 2.0  1997/09/09  13:16 EDT
               Timi Adelekan/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated Code.

               Revision 1.1  1997/08/27  10:00
               Tom Johnson  (johnson@ltpmail.gsfc.nasa.gov)
               Incorporate PDL walkthru comments

               Revision 1.0  1997/08/14  16:10 EDT
               Timi Adelekan/GSC (adelekan@ltpmail.gsfc.nasa.gov)
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

               This routine was designed totally under the assumption that
               the packet header has been previously validated and there
               is no need to test for error conditions.

!END************************************************************************
*/

 {
  /**************************************************************************/
  /*              Define and Initialize Local Variables                     */
  /*                                                                        */
  /**************************************************************************/

  int  start_pos;               /* Offset to structure location within pkt  */
  int  i;                       /* loop variable                            */


  /**************************************************************************/
  /*                                                                        */
  /* Use loop to move engineering 2, packet 1 contents to scan structure    */
  /* Suggested loop: FOR loop                                               */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /* FOR index equals 0 up to PD_E2P1_NUM_HK_TELEM_BYTES                    */
  /*    set SC_SCAN_DATA_t.raw_hk_telem[index] equal to                     */
  /*                   pkt[PD_E2P1_CURR_HK_BYTE_OFFSET + index]             */
  /* ENDFOR                                                                 */
  /*                                                                        */
  /**************************************************************************/

     for (i = 0; i < PD_E2P1_NUM_HK_TELEM_BYTES; i++)
        L1A_scan->raw_hk_telem[i] = pkt[PD_E2P1_CURR_HK_BYTE_OFFSET + i];


  /**************************************************************************/
  /*                                                                        */
  /* set start_pos equal to PD_E2P1_SC_ANCIL_BYTE_OFFSET                    */
  /*                                                                        */
  /* FOR index equal 0 upto PD_E2P1_NUM_SC_ANCIL_WORDS                      */
  /*    SC_SCAN_DATA_t.raw_sc_ancil[index] = pkt[start_pos] shifted left 8  */
  /*       bits OR pkt[start_pos+1]                                         */
  /*                                                                        */
  /*    increment start_pos by PD_NUM_BYTES_IN_WORD                         */
  /* ENDFOR                                                                 */
  /*                                                                        */
  /**************************************************************************/

     start_pos = PD_E2P1_SC_ANCIL_BYTE_OFFSET;

     for (i = 0; i < PD_E2P1_NUM_SC_ANCIL_WORDS; i++) {
        L1A_scan->raw_sc_ancil[i] = 
            ((int16)(pkt[start_pos]))<<8 | ((int16)(pkt[start_pos+1]));
        start_pos += PD_NUM_BYTES_IN_WORD;
     }

  /**************************************************************************/
  /*                                                                        */
  /* FOR index equals 0 up to PD_E2P1_NUM_PARAM_BYTES                       */
  /*    set SC_SCAN_DATA_t.raw_param[index] equal to                        */
  /*                   pkt[PD_E2P1_PARAM_BYTE_OFFSET + index]               */
  /* ENDFOR                                                                 */
  /*                                                                        */
  /**************************************************************************/

     for (i = 0; i < PD_E2P1_NUM_PARAM_BYTES; i++)
        L1A_scan->raw_param[i] = pkt[PD_E2P1_PARAM_BYTE_OFFSET + i];

 } /* End of routine output_engr2_pkt1_to_scan */
