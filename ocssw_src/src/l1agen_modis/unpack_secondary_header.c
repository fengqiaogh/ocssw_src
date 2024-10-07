#include "PGS_MODIS_35005.h"
#include "PGS_SMF.h"
#include "PGS_TD.h"
#include "PGS_IO.h"
#include "PGS_IO_L0.h"
#include "PGS_TYPES.h"
#include "PH_pkt_hdr.h"
#include "PD_pkt_data.h"
#include "L1A_prototype.h"

PGSt_SMF_status   unpack_secondary_header (PGSt_IO_L0_Packet   *pkt,
                                           PH_PACKET_HEADER_t  *packet_header)

/*
!C**********************************************************************

!Description:  This function extracts the secondary header information 
               contained in the MODIS packet and places each item into the 
               appropriate variable in the PH_PACKET_HEADER_t structure.

!Input Parameters:
               PGSt_IO_L0_Packet   *pkt            ** The MODIS packet         **

!Output Parameters:
               PH_PACKET_HEADER_t  *packet_header  ** Pointer to structure that 
                                                   contains unpacked contents 
                                                   of the primary header,
                                                   secondary header, and MODIS 
                                                   header                      **

Return Values: 
               MODIS_S_SUCCESS                    (PGS_MODIS_35005.h)
               MODIS_E_FAILED_TIMECODE_CONV       (PGS_MODIS_35005.h)
               MODIS_E_INV_QL_FLAG                (PGS_MODIS_35005.h)
               MODIS_E_INV_PKT_TYPE               (PGS_MODIS_35005.h)

Externally Defined:
               PGSt_IO_L0_Packet                  (PGS_IO.h)
               PGSt_SMF_status                    (PGS_SMF.h)
               PH_PACKET_HEADER_t                 (PH_pkt_hdr.h)
               PH_SEC_TIME_TAG_NUM_BYTES          (PH_pkt_hdr.h)
               PH_SEC_QUICK_LOOK_FLAG_BYTE_OFFSET (PH_pkt_hdr.h)
               PH_SEC_QUICK_LOOK_FLAG_BIT_OFFSET  (PH_pkt_hdr.h)
               PH_SEC_QUICK_LOOK_FLAG_NUM_BITS    (PH_pkt_hdr.h)
               PH_SEC_QUICK_LOOK_FLAG_SET         (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_BYTE_OFFSET        (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_BIT_OFFSET         (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_NUM_BITS           (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_ENG2_GROUP         (PH_pkt_hdr.h)
               PH_SEC_SCAN_CNT_BYTE_OFFSET        (PH_pkt_hdr.h)
               PH_SEC_SCAN_CNT_BIT_OFFSET         (PH_pkt_hdr.h)
               PH_SEC_SCAN_CNT_NUM_BITS           (PH_pkt_hdr.h)
               PH_SEC_MIRROR_SIDE_BYTE_OFFSET     (PH_pkt_hdr.h)
               PH_SEC_MIRROR_SIDE_BIT_OFFSET      (PH_pkt_hdr.h)
               PH_SEC_MIRROR_SIDE_NUM_BITS        (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_SPARE              (PH_pkt_hdr.h)

Called By:
               unpack_packet_header

Routines Called:
               extr_bits                       
               log_fmt_msg
               PGS_TD_EOSAMtoTAI
               PGS_SMF_TestStatusLevel 

!Revision History:
               Revision 2.2  2000/07/14
               John Seaton  (seaton@ltpmail.gsfc.nasa.gov)
               Implemented Aqua changes.

               Revision 2.1  1999/08/20  10:25
               John Seaton/SAIC/GSC (seaton@ltpmail.gsfc.nasa.gov)
               Updated the verification of the packet type in the
               secondary header to check for all spare values.

               Revision 2.0  1997/07/14  14:25
               Tom Johnson/SAIC/GSC (johnson@ltpmail.gsfc.nasa.gov)
               Created from unpack_packet_header and validate_packet

               Revision 1.0  1997/06/18  16:40 EDT
               Timi Adelekan/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Baseline from Version 1.

!Team-unique Header:
               This software is developed by the MODIS Science 
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes: 
               The details for the packet header data locations were taken
               from Hughes Santa Barbara Remote Sensing (SBRS) Contract Data
               Requirements List (CDRL) 305, MODIS Engineering Telemetry 
               Description, Figures 3-7 and 30-8.


!END************************************************************************
*/

{
  /***************************************************************************/
  /*              Declare and Initialize Local Variables                     */
  /***************************************************************************/

  PGSt_SMF_status  returnStatus;   /* SMF-style message returned by function */
  PGSt_SMF_status  PGS_status;     /* SMF-style message returned by function */
  char             msg[300];       /* amplifying message                     */
  int              i;              /* loop counter                           */

  char             *routine = "unpack_secondary_header";

  /***************************************************************************/


  /***************************************************************************/
  /*                                                                         */
  /*  Set returnStatus equal to MODIS_S_SUCCESS                              */
  /*                                                                         */
  /*  set routine to "unpack_secondary_header" (done during declaration)     */
  /*                                                                         */
  /***************************************************************************/

  returnStatus = MODIS_S_SUCCESS;

  
  /***************************************************************************/
  /*                                                                         */
  /*  FOR i equals all 8 bytes (PH_SEC_TIME_TAG_NUM_BYTES) in secondary      */
  /*     header time tag                                                     */
  /*     set packet_header->pkt_time_code[i] to pkt secondary header         */
  /*        time tag[i]                                                      */
  /*  ENDFOR                                                                 */
  /*                                                                         */
  /***************************************************************************/

  for (i = 0; i < PH_SEC_TIME_TAG_NUM_BYTES; i++)
     packet_header->pkt_time_code[i] = 
                       (PGSt_scTime) pkt[i+PH_SEC_TIME_TAG_BYTE_OFFSET];

  /***************************************************************************/
  /*                                                                         */
  /* The MODIS instrument on both the Terra and Aqua spacecraft use segmented*/
  /* time codes. This allows the toolkit function PGS_TD_EOSAMtoTAI to be    */
  /* used for both terra and aqua spacecrafts.                               */
  /*                                                                         */
  /* CALL PGS_TD_EOSAMtoTAI to convert the time tag from EOS AM format to TAI*/
  /*     INPUT:  packet_header->pkt_time_code                                */
  /*     OUTPUT: packet_header->pkt_TAI_time                                 */
  /*     RETURN: returnStatus                                                */
  /*                                                                         */
  /***************************************************************************/

    PGS_status = PGS_TD_EOSAMtoTAI(packet_header->pkt_time_code,
       &packet_header->pkt_TAI_time); 

  /***************************************************************************/
  /*                                                                         */
  /*  CALL PGS_SMF_TestStatusLevel to determine if the time conversion was   */
  /*     successful                                                          */
  /*    INPUT:  returnStatus                                                 */
  /*    OUTPUT: None                                                         */
  /*    RETURN: PGS_status                                                   */
  /*                                                                         */
  /*  IF PGS_status greater than or equal to PGS_SMF_MASK_LEV_E              */
  /*  THEN                                                                   */
  /*     set returnStatus to MODIS_E_FAILED_TIMECODE_CONV                    */
  /*     set msg to "unable to convert secondary header time tag from        */
  /*        spacecraft time to TAI"                                          */
  /*     CALL log_fmt_msg to report that the secondary header time tag was   */
  /*        not converted from Spacecraft time to TAI                        */
  /*       INPUTS: returnStatus, routine, msg                                */
  /*       OUTPUT: None                                                      */
  /*       RETURN: None                                                      */
  /*  ENDIF                                                                  */
  /*                                                                         */
  /***************************************************************************/

  if (PGS_SMF_TestStatusLevel(PGS_status) >= PGS_SMF_MASK_LEV_E) 
    {
     returnStatus = MODIS_E_FAILED_TIMECODE_CONV;
     strcpy (msg, "unable to convert secondary header time tag from Spacecraft time to TAI");
     log_fmt_msg (MODIS_E_FAILED_TIMECODE_CONV, routine, msg);
    }


  /***************************************************************************/
  /*                                                                         */
  /*  CALL extr_bits to extract the quick look flag from the packet secondary*/
  /*     header                                                              */
  /*    INPUTS: pkt, PH_SEC_QUICK_LOOK_FLAG_BIT_OFFSET,                      */
  /*            PH_SEC_QUICK_LOOK_FLAG_BYTE_OFFSET,                          */
  /*            PH_SEC_QUICK_LOOK_FLAG_NUM_BITS                              */
  /*    OUTPUT: None                                                         */
  /*    RETURN: packet_header->QL_flag                                       */
  /*                                                                         */
  /*  IF packet_header->QL_flag is equal to PH_SEC_QUICK_LOOK_FLAG_SET       */
  /*  THEN                                                                   */
  /*     set returnStatus to MODIS_E_INV_QL_FLAG                             */
  /*     set msg to "QL flag set to %d in the packet secondary header"       */
  /*     CALL log_fmt_msg to report that an invalid QL flag set has been     */
  /*        detected in the packet secondary header"                         */
  /*       INPUTS: returnStatus, routine, msg                                */
  /*       OUTPUT: None                                                      */
  /*       RETURN: None                                                      */
  /*  ENDIF                                                                  */
  /*                                                                         */
  /***************************************************************************/

  packet_header->QL_flag = extr_bits (pkt, PH_SEC_QUICK_LOOK_FLAG_BIT_OFFSET,
     PH_SEC_QUICK_LOOK_FLAG_BYTE_OFFSET, PH_SEC_QUICK_LOOK_FLAG_NUM_BITS);

  if (packet_header->QL_flag == PH_SEC_QUICK_LOOK_FLAG_SET)
    {
     returnStatus = MODIS_E_INV_QL_FLAG;
     sprintf(msg, "QL flag set to %d in the packet secondary header", 
             packet_header->QL_flag);
     log_fmt_msg (MODIS_E_INV_QL_FLAG, routine, msg);
    }


  /***************************************************************************/
  /*                                                                         */
  /*  CALL extr_bits to extract the packet type from the packet secondary    */
  /*     header                                                              */
  /*    INPUTS: pkt, PH_SEC_PKT_TYPE_BIT_OFFSET, PH_SEC_PKT_TYPE_BYTE_OFFSET,*/ 
  /*            PH_SEC_PKT_TYPE_NUM_BITS                                     */
  /*    OUTPUT: None                                                         */
  /*    RETURN: packet_header->pkt_type                                      */
  /*                                                                         */
  /***************************************************************************/

  packet_header->pkt_type = extr_bits (pkt, PH_SEC_PKT_TYPE_BIT_OFFSET, 
     PH_SEC_PKT_TYPE_BYTE_OFFSET, PH_SEC_PKT_TYPE_NUM_BITS);

  /***************************************************************************/
  /*                                                                         */
  /* The pkt_type is a 3 bit field in the secondary header. The valid values */
  /* for this field are 0, 1, 2, and 4. The other possible values 3, 5, 6,   */
  /* and 7 are spare values.                                                 */
  /***************************************************************************/
  /*                                                                         */
  /* IF the value of the pkt_type is a spare value                           */
  /* THEN                                                                    */
  /*   set returnStatus to MODIS_E_INV_PKT_TYPE                              */
  /*   set msg to "Packet Type: %d in the packet secondary header"           */
  /*   CALL log_fmt_msg to print an error message in LogStatus stating that  */
  /*      an invalid packet type has been detected in the packet secondary   */
  /*       header.                                                           */
  /*     INPUTS: MODIS_E_INV_PKT_TYPE, routine, msg                          */
  /*     OUTPUTS: None                                                       */
  /*     RETURNS: None                                                       */
  /* ENDIF                                                                   */
  /*                                                                         */
  /***************************************************************************/

  if (packet_header->pkt_type > PH_SEC_PKT_TYPE_ENG2_GROUP ||
      packet_header->pkt_type == PH_SEC_PKT_TYPE_SPARE)
    {
     returnStatus = MODIS_E_INV_PKT_TYPE;
     sprintf(msg, "Packet Type: %d in the packet secondary header", 
             packet_header->pkt_type);
     log_fmt_msg (MODIS_E_INV_PKT_TYPE, routine, msg);
    }


  /***************************************************************************/
  /*                                                                         */
  /*  CALL extr_bits to extract the scan count from the packet secondary     */
  /*     header                                                              */
  /*     INPUTS: pkt, PH_SEC_SCAN_CNT_BIT_OFFSET,                            */
  /*             PH_SEC_SCAN_CNT_BYTE_OFFSET, PH_SEC_SCAN_CNT_NUM_BITS       */
  /*     OUTPUT: None                                                        */
  /*     RETURN: packet_header->scan_cnt                                     */
  /*                                                                         */
  /***************************************************************************/

  packet_header->scan_cnt = extr_bits (pkt, PH_SEC_SCAN_CNT_BIT_OFFSET,
     PH_SEC_SCAN_CNT_BYTE_OFFSET, PH_SEC_SCAN_CNT_NUM_BITS);


  /***************************************************************************/
  /*                                                                         */
  /*  CALL extr_bits to extract the mirror side from the packet secondary    */
  /*     header                                                              */
  /*    INPUTS: pkt, PH_SEC_MIRROR_SIDE_BIT_OFFSET,                          */
  /*            PH_SEC_MIRROR_SIDE_BYTE_OFFSET, PH_SEC_MIRROR_SIDE_NUM_BITS  */
  /*    OUTPUT: None                                                         */
  /*    RETURN: packet_header->mirror_side                                   */
  /*                                                                         */
  /***************************************************************************/

  packet_header->mirror_side = extr_bits (pkt, PH_SEC_MIRROR_SIDE_BIT_OFFSET,
     PH_SEC_MIRROR_SIDE_BYTE_OFFSET, PH_SEC_MIRROR_SIDE_NUM_BITS);


  /***************************************************************************/
  /*                                                                         */
  /*  RETURN  returnStatus                                                   */
  /*                                                                         */
  /***************************************************************************/

  return (returnStatus);

}
