#include "L1A_prototype.h"
#include "PGS_IO.h"
#include "PH_pkt_hdr.h"
#include "SC_scan.h"


void   output_eng_data_to_scan (PH_PACKET_HEADER_t  *pkt_header,
                                PGSt_IO_L0_Packet   *pkt,
                                SC_SCAN_DATA_t      *L1A_scan)

/*
!C************************************************************************

!Description:  This routine determines what type of engineering packet the 
               current packet is, and calls the appropriate routines to output 
               the packet's data to the scan's structure

!Input Parameters:
               PH_PACKET_HEADER_t  *pkt_header   ** Buffer that contains pkt **
                                                 ** header info of the pckt  **

               PGSt_IO_L0_Packet   *pkt          ** The structure containing **
                                                 ** the current packet's     **
                                                 ** unpacked contents        **

!Output Parameters:
               None

!Input/Output Parameters:
               SC_SCAN_DATA_t      *L1A_scan     ** The MODIS scan structure **
                                                 ** currently being built    **

Return Values: None

Externally Defined:
               PGSt_IO_L0_Packet                      (PGS_IO.h)
               PH_PACKET_HEADER_t                     (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_ENG1_GROUP             (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_ENG2_GROUP             (PH_pkt_hdr.h)
               PH_PRI_SEQUENCE_FIRST_PKT_IN_GROUP     (PH_pkt_hdr.h)
               PH_PRI_SEQUENCE_SECOND_PKT_IN_GROUP    (PH_pkt_hdr.h)
               SC_SCAN_DATA_t                         (SC_scan.h)

Called By:
               put_pkt_cont_in_scan

Routines Called:
               output_eng1_pkt1_to_scan          
               output_eng1_pkt2_to_scan          
               output_eng2_pkt1_to_scan          
               output_eng2_pkt2_to_scan          

!Revision History:
               Revision 2.0  1997/08/15  07:30 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated Code for version 2.

               Revision 1.1  1997/08/27   9:45
               Tom Johnson  (johnson@ltpmail.gsfc.nasa.gov)
               Incorporate PDL walkthru comments

               Revision 1.0  1997/06/18  16:40 EDT  
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Original design.          

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

!END*************************************************************************
*/

 {
  /**************************************************************************/
  /*                                                                        */
  /* Determine what type of engineering packet the current packet is then   */
  /* call the appropiate routine to output the data to the scan structure.  */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /* IF PH_PACKET_HEADER_t.pkt_type is equal to PH_SEC_PKT_TYPE_ENG1_GROUP  */
  /* THEN                                                                   */
  /*    IF PH_PACKET_HEADER_t.sequence_flag is equal to                     */
  /*                                    PH_PRI_SEQUENCE_FIRST_PKT_IN_GROUP  */
  /*    THEN                                                                */
  /*       CALL output_eng1_pkt1_to_scan to output the Engineering group 1  */
  /*            packet 1 data to the scan structure                         */
  /*         INPUTS:  pkt, SC_SCAN_DATA_t                                   */
  /*         OUTPUTS: SC_SCAN_DATA_t                                        */
  /*         RETURN:  NONE                                                  */
  /*    ELSE                                                                */
  /*       CALL output_eng1_pkt2_to_scan to output the Engineering group 1  */
  /*            packet 2 data to the scan structure                         */
  /*         INPUTS:  pkt, SC_SCAN_DATA_t                                   */
  /*         OUTPUTS: SC_SCAN_DATA_t                                        */
  /*         RETURN:  NONE                                                  */
  /*    ENDIF                                                               */
  /* ELSE                                                                   */
  /*    IF PH_PACKET_HEADER_t.sequence_flag is equal to                     */
  /*                                    PH_PRI_SEQUENCE_FIRST_PKT_IN_GROUP  */
  /*    THEN                                                                */
  /*       CALL output_eng2_pkt1_to_scan to output the Engineering group 2  */
  /*            packet 1 data to the scan structure                         */
  /*         INPUTS:  pkt, SC_SCAN_DATA_t                                   */
  /*         OUTPUTS: SC_SCAN_DATA_t                                        */
  /*         RETURN:  NONE                                                  */
  /*    ELSE                                                                */
  /*       CALL output_eng2_pkt2_to_scan to output the Engineering group 2  */
  /*            packet 2 data to the scan structure                         */
  /*         INPUTS:  pkt, SC_SCAN_DATA_t                                   */
  /*         OUTPUTS: SC_SCAN_DATA_t                                        */
  /*         RETURN:  NONE                                                  */
  /*    ENDIF                                                               */
  /* ENDIF                                                                  */
  /*                                                                        */
  /**************************************************************************/

     if (pkt_header->pkt_type == PH_SEC_PKT_TYPE_ENG1_GROUP) {
        if (pkt_header->sequence_flag == PH_PRI_SEQUENCE_FIRST_PKT_IN_GROUP)
           output_eng1_pkt1_to_scan (pkt, L1A_scan);
        else
           output_eng1_pkt2_to_scan (pkt, L1A_scan);
     }

     else {
        if (pkt_header->sequence_flag == PH_PRI_SEQUENCE_FIRST_PKT_IN_GROUP)
           output_eng2_pkt1_to_scan (pkt, L1A_scan);
        else
           output_eng2_pkt2_to_scan (pkt, L1A_scan);
     }

 } /* End of routine output_eng_data_to_scan */
