#include "L1A_prototype.h"
#include "hdf.h"
#include "PGS_IO_L0.h"
#include "EN_eng_data.h"
#include "PGS_SMF.h"
#include "PH_pkt_hdr.h"
#include "PGS_MODIS_35005.h"


void   process_eng_packet (EN_VDATA_TYPE_t       *eng_data, 
                           int                   scan_number, 
                           PH_PACKET_HEADER_t    *pkt_header,
                           PGSt_IO_L0_Packet     *pkt)

/*
!C*****************************************************************************

!Description:  This function is a router that processes eng packets.

!Input Parameters:
               uint16             scan_number    ** The scan number        **
	                                         ** (counting from 1) of   **
                                                 ** the current scan within**
                                                 ** the current granule    **

               PH_PACKET_HEADER_t *pkt_header    ** The header of the eng  **
	                                         ** packet                 **

               PGSt_IO_LO_Packet  *pkt           ** The current eng packet **

!Output Parameters:
               None

!Input/Output Parameters:
               EN_VDATA_TYPE_t     eng_data      ** The Vdata array        **
               					 ** structure              **

Return Values:
               None

Externally Defined: 
               EN_VDATA_TYPE_t              (EN_eng_data.h)
               PGSt_SMF_status              (PGS_SMF.h)
               PGSt_IO_LO_Packet            (PGS_IO_L0.h)
               PH_PACKET_HEADER_t           (PH_pkt_hdr.h)
               MODIS_E_NULL_POINTER         (PGS_MODIS_35005.h)

Called By:
               load_eng_data
               process_a_scan

Routines Called:
               process_cp_hk_tlmy
               process_sci_eng_data
               log_fmt_msg

!Revision History:
               Revision 2.1  2001/01/04
               John Seaton  (seaton@ltpmail.gsfc.nasa.gov)
               Simplified code. Removed 4 subroutine calls.

               Revision 2.0  1998/10/26  10:22 EST
               John Seaton/SAIC/GSC  (seaton@ltpmail.gsfc.nasa.gov)
               Added calls to process Current/Prior S/C Ancillary
               Data Vdatas, and Command Parameters Vdata.

               revision 1.0 1997/09/11  17:30:00
               Qi Huang/RDC    (qhuang@ltpmail.gsfc.nasa.gov)
               Original development

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes:
               None

!END******************************************************************************
*/
{

  /*****************************************************************************/
  /* Set char pointer routine to "process_eng_packet"                          */
  /*****************************************************************************/

  char *routine="process_eng_packet";

  /*****************************************************************************/
  /* IF pkt_header equals NULL                                                 */
  /* THEN                                                                      */
  /*    CALL log_fmt_msg to report the error in LogStatus file                 */
  /*      INPUTS: MODIS_E_NULL_POINTER, routine, "Packet Header is NULL"       */
  /*      OUTPUTS: None                                                        */
  /*      RETURN:                                                              */
  /*    return                                                                 */
  /* ENDIF                                                                     */
  /*****************************************************************************/

  if (pkt_header == NULL) {
    log_fmt_msg(MODIS_E_NULL_POINTER, routine, "Packet Header is NULL");  
    return;
  }

  /*****************************************************************************/
  /* ELSE                                                                      */
  /*   IF ( pkt_header->sequence_flag == PH_PRI_SEQUENCE_FIRST_PKT_IN_GROUP)   */
  /*       AND (pkt_header->pkt_type == PH_SEC_PKT_TYPE_ENG2_GROUP)            */
  /*   THEN                                                                    */
  /*      CALL process_cp_hk_tlmy to extract the current and prior HK telemetry*/
  /*         information from this packet and place it into the appropriate    */
  /*         part of the eng_data structure                                    */
  /*        INPUTS:  eng_data, pkt, scan_number                                */
  /*        OUTPUTS: eng_data                                                  */
  /*        RETURN:  None                                                      */
  /*                                                                           */
  /*      CALL process_group1_packet2_vdata to process the Prior/Current       */
  /*           S/C Ancillary data along with the Command Procedures SDS.       */
  /*           These data will be stored as vdatas.                            */
  /*       INPUTS: pkt, scan_number, eng_data                                  */
  /*       OUTPUTS: eng_data                                                   */
  /*       RETURN: None                                                        */
  /*                                                                           */
  /*      CALL update_eng_data to write the Command Parameters data to the     */
  /*           Vdata array (eng_data).                                         */
  /*        INPUTS: EN_COMM_PROC_VDATA_NUMBER (49), pkt, scan_number, eng_data,*/
  /*                FALSE                                                      */
  /*        OUTPUTS: eng_data                                                  */
  /*        RETURNS: None                                                      */
  /*                                                                           */
  /*   ELSE                                                                    */
  /*                                                                           */
  /*      IF (pkt_header->pkt_type == PH_SEC_PKT_TYPE_ENG1_GROUP) AND          */
  /*         (pkt_header->sequence_flag == PH_PRI_SEQUENCE_SECOND_PKT_IN_GROUP)*/
  /*      THEN                                                                 */
  /*                                                                           */
  /*         CALL process_sci_eng_data to process this packet                  */
  /*           INPUTS:  eng_data, pkt, scan_number                             */
  /*           OUTPUTS: eng_data                                               */
  /*           RETURN:  None                                                   */
  /*                                                                           */
  /*      ENDIF                                                                */
  /*   ENDIF                                                                   */
  /* ENDELSE                                                                   */
  /*****************************************************************************/

  else { 
    if ((pkt_header->sequence_flag == PH_PRI_SEQUENCE_FIRST_PKT_IN_GROUP) &&
        (pkt_header->pkt_type == PH_SEC_PKT_TYPE_ENG2_GROUP)) {
        process_cp_hk_tlmy(eng_data,pkt,scan_number);
        process_group2_packet1_vdata(pkt,eng_data); 
        update_eng_data(EN_COMM_PROC_VDATA_NUMBER,pkt,scan_number,eng_data,FALSE);
    }
    else
        if ((pkt_header->pkt_type == PH_SEC_PKT_TYPE_ENG1_GROUP) &&
            (pkt_header->sequence_flag == PH_PRI_SEQUENCE_SECOND_PKT_IN_GROUP))
          process_sci_eng_data(eng_data,pkt,scan_number);
  }

  return;
}
