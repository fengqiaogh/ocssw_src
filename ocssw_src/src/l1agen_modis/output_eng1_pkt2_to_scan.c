#include "L1A_prototype.h"
#include "PGS_IO.h"
#include "PD_pkt_data.h"
#include "SC_scan.h"
#include "hdfi.h"

void   output_eng1_pkt2_to_scan (PGSt_IO_L0_Packet  *pkt,
                                 SC_SCAN_DATA_t     *L1A_scan)

/*
!C************************************************************************

!Description:  This routine extracts the engineering1, packet2 data contents
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
               PGSt_IO_L0_Packet                           (PGS_IO.h)
               PD_FIRST_BIT_IN_BYTE                        (PD_pkt_data.h)
               PD_NUM_BYTES_IN_WORD                        (PD_pkt_data.h)
               PD_E1P2_NUM_EARTH_ENCODER_TIMES             (PD_pkt_data.h)
               PD_E1P2_NUM_VIEW_SECTOR_DEFINITIONS         (PD_pkt_data.h)
               PD_E1P2_NUM_VIEW_SECTOR_ACTUALS             (PD_pkt_data.h)
               PD_E1P2_NUM_SCI_ENG_BYTES                   (PD_pkt_data.h)
               PD_E1P2_EARTH_ENCODER_TIMES_BYTE_OFFSET     (PD_pkt_data.h)
               PD_E1P2_VIEW_SECTOR_DEFINITIONS_BYTE_OFFSET (PD_pkt_data.h)
               PD_E1P2_VIEW_SECTOR_ACTUALS_BYTE_OFFSET     (PD_pkt_data.h)
               PD_E1P2_SCI_ENG_BYTE_OFFSET                 (PD_pkt_data.h)
               SC_SCAN_DATA_t                              (SC_scan.h)

Called By:
               output_eng_data_to_scan

Routines Called:
               None

!Revision History:
               Revision 3.0  2001/04/13
               John Seaton (seaton@ltpmail.gsfc.nasa.gov)
               Changed raw_mir_enc to uint16

               Revision 2.0  1997/09/09  13:10 EDT
               Timi Adelekan/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated Code.

               Revision 1.1  1997/08/27   9:50
               Tom Johnson  (johnson@ltpmail.gsfc.nasa.gov)
               Incorporate PDL walkthru comments

               Revision 1.0  1997/08/14  15:45 EDT
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
  /*************************************************************************/
  /*                                                                       */
  /*              Define and Initialize Local Variables                    */
  /*                                                                       */
  /*************************************************************************/

  int  start_pos;               /* Offset to structure location within pkt */
  int  i;                       /* loop variable                           */



  /*************************************************************************/
  /*                                                                       */
  /* Use loop to move engineering 1, packet 2 contents to scan structure   */
  /* Suggested loop: FOR loop                                              */
  /*                                                                       */
  /*************************************************************************/
  /*                                                                       */
  /* set start_pos equal to PD_E1P2_EARTH_ENCODER_TIMES_BYTE_OFFSET        */
  /*                                                                       */
  /* FOR index equal 0 upto PD_E1P2_NUM_EARTH_ENCODER_TIMES                */
  /*    SC_SCAN_DATA_t.raw_mir_enc[index] = pkt[start_pos] shifted left 8  */
  /*       bits OR pkt[start_pos+1]                                        */
  /*                                                                       */
  /*    increment start_pos by PD_NUM_BYTES_IN_WORD                        */
  /* ENDFOR                                                                */
  /*                                                                       */
  /*************************************************************************/

     start_pos = PD_E1P2_EARTH_ENCODER_TIMES_BYTE_OFFSET;

     for (i = 0; i < PD_E1P2_NUM_EARTH_ENCODER_TIMES; i++) {
        L1A_scan->raw_mir_enc[i] = 
            ((int16)(pkt[start_pos]))<<8 | ((int16)(pkt[start_pos+1]));

        start_pos += PD_NUM_BYTES_IN_WORD;
     }



  /*************************************************************************/
  /*                                                                       */
  /* set start_pos equal to PD_E1P2_VIEW_SECTOR_DEFINITIONS_BYTE_OFFSET    */
  /*                                                                       */
  /* FOR index equal 0 upto PD_E1P2_NUM_VIEW_SECTOR_DEFINITIONS            */
  /*    SC_SCAN_DATA_t.raw_vs_def[index] = pkt[start_pos] shifted left 8   */
  /*       bits OR pkt[start_pos+1]                                        */
  /*                                                                       */
  /*    increment start_pos by PD_NUM_BYTES_IN_WORD                        */
  /* ENDFOR                                                                */
  /*                                                                       */
  /*************************************************************************/

     start_pos = PD_E1P2_VIEW_SECTOR_DEFINITIONS_BYTE_OFFSET;

     for (i = 0; i < PD_E1P2_NUM_VIEW_SECTOR_DEFINITIONS; i++) {
        L1A_scan->raw_vs_def[i] = 
            ((int16)(pkt[start_pos]))<<8 | ((int16)(pkt[start_pos+1]));
        start_pos += PD_NUM_BYTES_IN_WORD;
     }



  /*************************************************************************/
  /*                                                                       */
  /* set start_pos equal to PD_E1P2_VIEW_SECTOR_ACTUALS_BYTE_OFFSET        */
  /*                                                                       */
  /* FOR index equal 0 upto PD_E1P2_NUM_VIEW_SECTOR_ACTUALS                */
  /*    SC_SCAN_DATA_t.raw_vs_act[index] = pkt[start_pos] shifted left 8   */
  /*       bits OR pkt[start_pos+1]                                        */
  /*                                                                       */
  /*    increment start_pos by PD_NUM_BYTES_IN_WORD                        */
  /* ENDFOR                                                                */
  /*                                                                       */
  /*************************************************************************/

     start_pos = PD_E1P2_VIEW_SECTOR_ACTUALS_BYTE_OFFSET;

     for (i = 0; i < PD_E1P2_NUM_VIEW_SECTOR_ACTUALS; i++) {
        L1A_scan->raw_vs_act[i] = 
            ((int16)(pkt[start_pos]))<<8 | ((int16)(pkt[start_pos+1]));

        start_pos += PD_NUM_BYTES_IN_WORD;
     }



  /*************************************************************************/
  /*                                                                       */
  /* FOR index equals 0 up to PD_E1P2_NUM_SCI_ENG_BYTES                    */
  /*   set SC_SCAN_DATA_t.raw_sci_eng[index] equal to                      */
  /*                              pkt[PD_E1P2_SCI_ENG_BYTE_OFFSET + index] */
  /* ENDFOR                                                                */
  /*                                                                       */
  /*************************************************************************/

     for (i = 0; i < PD_E1P2_NUM_SCI_ENG_BYTES; i++)
        L1A_scan->raw_sci_eng[i] = pkt[PD_E1P2_SCI_ENG_BYTE_OFFSET + i];

 } /* End of routine output_eng1_pkt2_to_scan */
