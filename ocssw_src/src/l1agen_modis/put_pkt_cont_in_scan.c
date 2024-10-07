#include "L1A_prototype.h"
#include "PGS_IO.h"
#include "PH_pkt_hdr.h"
#include "SC_scan.h"


void   put_pkt_cont_in_scan (PH_PACKET_HEADER_t  pkt_header,
                             PGSt_IO_L0_Packet   *pkt,
                             uint16              *pkt_contents,
                             SC_SCAN_DATA_t      *L1A_scan)

/*
!C************************************************************************

!Description:  This routine will fill the appropriate part of the current scan 
               with the data from the current packet.  Depending on what type 
               the current packet is, the routine will call the appropriate 
               output routine to output the packet data to the scan.

!Input Parameters:
               PGSt_IO_L0_Packet   *pkt          ** The structure containing **
                                                 ** the current packet's     **
                                                 ** packed data              **

               uint16              *pkt_contents ** The structure containing **
                                                 ** the current packet's     **
                                                 ** unpacked contents        **

               PH_PACKET_HEADER_t  pkt_header    ** Buffer that contains pkt **
                                                 ** header info of the pkt   **

!Output Parameters:
               None

!Input/Output Parameters:
               SC_SCAN_DATA_t      *L1A_scan     ** The MODIS scan that is   **
                                                 ** currently being built    **

Return Values: 
               None

Externally Defined:
               PGSt_IO_L0_Packet                      (PGS_IO.h)
               SC_SCAN_DATA_t                         (SC_scan.h)
               PH_PACKET_HEADER_t                     (PH_pkt_hdr.h)
               PH_MOD_FPA_AEM_CONFIG_NUM_ELEMENTS     (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_DAY_GROUP              (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_NIGHT_GROUP            (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_ENG1_GROUP             (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_ENG2_GROUP             (PH_pkt_hdr.h)

Called By:
               process_a_scan

Routines Called:
               output_daymode_data_to_scan    
               output_nightmode_data_to_scan  
               output_eng_data_to_scan       

!Revision History:
               Revision 2.0  1997/09/09  10:40 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Oringinated code for version 2.

               Revision 1.1  1997/08/27   9:30
               Tom Johnson  (johnson@ltpmail.gsfc.nasa.gov)
               Incorporate PDL walkthru comments

               Revision 1.0  1997/06/18  16:40 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Baseline from Version 1.

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes: 
               The CODE below was developed in C language.

               This routine assumes that the data contained in the
               input parameters have been validated before this
               routine is executed.  No validation of data is performed
               in this routine.

!END*************************************************************************
*/

 {
  /**************************************************************************/
  /*                                                                        */
  /*              Declare and Initialize Local Variables                    */
  /*                                                                        */
  /**************************************************************************/

  int  i;


  /**************************************************************************/
  /*                                                                        */
  /* Set SC_SCAN_DATA_t.science_state to PH_PACKET_HEADER_t.sci_state       */
  /*                                                                        */
  /* Set SC_SCAN_DATA_t.science_abnormal to PH_PACKET_HEADER_t.sci_abnorm   */
  /*                                                                        */
  /* Set all elements of SC_SCAN_DATA_t.fpa_aem_config to                   */
  /*    PH_PACKET_HEADER_t.fpa_aem_config                                   */
  /*                                                                        */
  /**************************************************************************/

  L1A_scan->science_state = pkt_header.sci_state;

  L1A_scan->science_abnormal = pkt_header.sci_abnorm;

  for (i=0; i<PH_MOD_FPA_AEM_CONFIG_NUM_ELEMENTS; i++)
     L1A_scan->fpa_aem_config[i] = pkt_header.fpa_aem_config[i];


  /**************************************************************************/
  /*                                                                        */
  /* Determine which group the current packet falls under and then call the */
  /* appropriate routine to output the packet data to the scan structure.   */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /* IF PH_PACKET_HEADER_t.pkt_type is equal to PH_SEC_PKT_TYPE_DAY_GROUP   */
  /* THEN                                                                   */
  /*    CALL output_daymode_data_to_scan to output the packet data to the   */
  /*        scan structure if the packet is a daymode packet                */
  /*      INPUTS:  PH_PACKET_HEADER_t, pkt_contents                         */
  /*      OUTPUTS: SC_SCAN_DATA_t                                           */
  /*      RETURN:  NONE                                                     */
  /* ELSE                                                                   */
  /*    IF PH_PACKET_HEADER_t.pkt_type is equal to                          */
  /*                                  PH_SEC_PKT_TYPE_NIGHT_GROUP           */
  /*    THEN                                                                */
  /*       CALL output_nightmode_data_to_scan to output the packet data to  */
  /*            the scan structure if the pkt is a nightmode packet         */
  /*         INPUTS:  PH_PACKET_HEADER_t, pkt_contents                      */
  /*         OUTPUTS: SC_SCAN_DATA_t.EV_1km_night                           */
  /*         RETURN:  NONE                                                  */
  /*    ENDIF                                                               */
  /* ELSE                                                                   */
  /*    IF ( (PH_PACKET_HEADER_t.pkt_type is equal to                       */
  /*                                    PH_SEC_PKT_TYPE_ENG1_GROUP)  OR     */
  /*         (PH_PACKET_HEADER_t.pkt_type is equal to                       */
  /*                                    PH_SEC_PKT_TYPE_ENG2_GROUP) )       */
  /*    THEN                                                                */
  /*       CALL output_eng_data_to_scan to output the packet data to the    */
  /*          scan structure if the packet is an engineering packet         */
  /*         INPUTS:  PH_PACKET_HEADER_t, pkt                               */
  /*         OUTPUTS: SC_SCAN_DATA_t                                        */
  /*         RETURN:  NONE                                                  */
  /*    ENDIF                                                               */
  /* ENDIF                                                                  */
  /*                                                                        */
  /**************************************************************************/

     if (pkt_header.pkt_type == PH_SEC_PKT_TYPE_DAY_GROUP)
        output_daymode_data_to_scan (&pkt_header, pkt_contents, L1A_scan);

     else if (pkt_header.pkt_type == PH_SEC_PKT_TYPE_NIGHT_GROUP)
        output_nightmode_data_to_scan (&pkt_header, 
                                       pkt_contents, 
                                        L1A_scan->EV_1km_night);

     else if ((pkt_header.pkt_type == PH_SEC_PKT_TYPE_ENG1_GROUP) ||
           (pkt_header.pkt_type == PH_SEC_PKT_TYPE_ENG2_GROUP))
        output_eng_data_to_scan (&pkt_header, pkt, L1A_scan); 

 } /* End of routine put_pkt_cont_in_scan */
