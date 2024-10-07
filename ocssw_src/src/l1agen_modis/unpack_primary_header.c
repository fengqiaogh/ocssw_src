#include "PGS_IO.h"
#include "PGS_SMF.h"
#include "PGS_MODIS_35005.h"
#include "PH_pkt_hdr.h"
#include "L1A_prototype.h"

PGSt_SMF_status   unpack_primary_header (PGSt_IO_L0_Packet   *pkt,
                                         PH_PACKET_HEADER_t  *packet_header)

/*
!C************************************************************************

!Description:  This function extracts the primary header information contained
               in the MODIS packet and places each item into the appropriate
               variable in the PH_PACKET_HEADER_t structure. It also validates each 
               item.

!Input Parameters:
               PGSt_IO_L0_Packet   *pkt             ** The MODIS packet       **

!Output Parameters:
               PH_PACKET_HEADER_t  *packet_header   ** Pointer to structure that 
                                                    contains unpacked contents 
                                                    of the primary header, secondary
                                                    header, and MODIS header  **

Return Values: 
               MODIS_S_SUCCESS                    (PGS_MODIS_35005.h)
               MODIS_E_INV_VERSION                (PGS_MODIS_35005.h)
               MODIS_E_INV_TYPE                   (PGS_MODIS_35005.h)
               MODIS_E_INV_SEC_HDR_FLAG           (PGS_MODIS_35005.h)
               MODIS_E_INV_APID                   (PGS_MODIS_35005.h)
               MODIS_E_INV_APID_TEST              (PGS_MODIS_35005.h)
               MODIS_E_INV_PKT_SEQ_FLAG           (PGS_MODIS_35005.h)

Externally Defined:
               PGSt_IO_L0_Packet                  (PGS_IO.h)
               PGSt_SMF_status                    (PGS_SMF.h)
               PH_PACKET_HEADER_t                 (PH_pkt_hdr.h)
               PH_PRI_VERSION_BYTE_OFFSET         (PH_pkt_hdr.h)
               PH_PRI_VERSION_BIT_OFFSET          (PH_pkt_hdr.h)
               PH_PRI_VERSION_NUM_BITS            (PH_pkt_hdr.h)
               PH_PRI_VERSION_VALUE               (PH_pkt_hdr.h)
               PH_PRI_TYPE_BYTE_OFFSET            (PH_pkt_hdr.h)
               PH_PRI_TYPE_BIT_OFFSET             (PH_pkt_hdr.h)
               PH_PRI_TYPE_NUM_BITS               (PH_pkt_hdr.h)
               PH_PRI_TYPE_VALUE                  (PH_pkt_hdr.h)
               PH_PRI_SEC_HDR_FLAG_BYTE_OFFSET    (PH_pkt_hdr.h)
               PH_PRI_SEC_HDR_FLAG_BIT_OFFSET     (PH_pkt_hdr.h)
               PH_PRI_SEC_HDR_FLAG_NUM_BITS       (PH_pkt_hdr.h)
               PH_PRI_SEC_HDR_PRESENT             (PH_pkt_hdr.h)
               PH_PRI_APID_BYTE_OFFSET            (PH_pkt_hdr.h)
               PH_PRI_APID_BIT_OFFSET             (PH_pkt_hdr.h)
               PH_PRI_APID_NUM_BITS               (PH_pkt_hdr.h)
               PH_PRI_MIN_MODIS_APID_AM1          (PH_pkt_hdr.h) 
               PH_PRI_MAX_MODIS_APID_AM1          (PH_pkt_hdr.h)
               PH_PRI_APID_TEST_PACKET            (PH_pkt_hdr.h)
               PH_PRI_SEQUENCE_FLAG_BYTE_OFFSET   (PH_pkt_hdr.h)
               PH_PRI_SEQUENCE_FLAG_BIT_OFFSET    (PH_pkt_hdr.h)
               PH_PRI_SEQUENCE_FLAG_NUM_BITS      (PH_pkt_hdr.h)
               PH_PRI_SEQUENCE_NOT_USED           (PH_pkt_hdr.h)
               PH_PRI_SOURCE_SEQ_CNT_BYTE_OFFSET  (PH_pkt_hdr.h)
               PH_PRI_SOURCE_SEQ_CNT_BIT_OFFSET   (PH_pkt_hdr.h)
               PH_PRI_SOURCE_SEQ_CNT_NUM_BITS     (PH_pkt_hdr.h)
               PH_PRI_PKT_LENGTH_BYTE_OFFSET      (PH_pkt_hdr.h)
               PH_PRI_PKT_LENGTH_BIT_OFFSET       (PH_pkt_hdr.h)
               PH_PRI_PKT_LENGTH_NUM_BITS         (PH_pkt_hdr.h)

Called By:
               unpack_packet_header

Routines Called:
               extr_bits   
               log_fmt_msg

!Revision History:
               Revision 2.0  1997/07/14  12:30
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
  /*                                                                         */
  /*              Declare and Initialize Local Variables                     */
  /*                                                                         */
  /***************************************************************************/

  PGSt_SMF_status  returnStatus;   /* SMF-style message returned by function */
  char             *routine = "unpack_primary_header";
  char             msg[300];


  /***************************************************************************/
  /*                                                                         */
  /*  Set returnStatus equal to MODIS_S_SUCCESS                              */
  /*                                                                         */
  /*  set routine to "unpack_primary_header" (done during declaration)       */
  /*                                                                         */
  /***************************************************************************/

  returnStatus = MODIS_S_SUCCESS;

  
  /***************************************************************************/
  /*                                                                         */
  /*  CALL extr_bits to extract the version from the packet primary header   */
  /*    INPUTS:  PH_PRI_VERSION_BYTE_OFFSET, PH_PRI_VERSION_BIT_OFFSET,      */ 
  /*             PH_PRI_VERSION_NUM_BITS, pkt                                */
  /*    OUTPUTS: None                                                        */
  /*    RETURN:  PH_PACKET_HEADER_t.version                                  */
  /*                                                                         */
  /*  IF PH_PACKET_HEADER_t.version is not equal to PH_PRI_VERSION_VALUE     */
  /*  THEN                                                                   */
  /*     set returnStatus equal to MODIS_E_INV_VERSION                       */
  /*     CALL log_fmt_msg to report that an invalid version value has been   */
  /*        detected in the packet primary header                            */
  /*       INPUTS: returnStatus, routine, msg                                */
  /*  ENDIF                                                                  */
  /*                                                                         */
  /***************************************************************************/

  packet_header->version = extr_bits(pkt, PH_PRI_VERSION_BIT_OFFSET,
     PH_PRI_VERSION_BYTE_OFFSET, PH_PRI_VERSION_NUM_BITS);

  if (packet_header->version != PH_PRI_VERSION_VALUE)
    {
     returnStatus = MODIS_E_INV_VERSION;
     sprintf(msg, "version number from packet primary header: %d valid version number %d", 
             packet_header->version, PH_PRI_VERSION_VALUE);
     log_fmt_msg(MODIS_E_INV_VERSION, routine, msg);
    }


  /***************************************************************************/
  /*                                                                         */
  /*  CALL extr_bits to extract the type from the packet primary header      */
  /*    INPUTS:  PH_PRI_TYPE_BYTE_OFFSET, PH_PRI_TYPE_BIT_OFFSET,            */
  /*             PH_PRI_TYPE_NUM_BITS, pkt                                   */
  /*    OUTPUTS: None                                                        */
  /*    RETURN:  PH_PACKET_HEADER_t.type                                     */
  /*                                                                         */
  /*  IF PH_PACKET_HEADER_t.type is not equal to PH_PRI_TYPE_VALUE           */
  /*  THEN                                                                   */
  /*     set returnStatus equal to MODIS_E_INV_TYPE                          */
  /*     CALL log_fmt_msg to report that a test packet has been detected from*/
  /*        information in the type field of the packet primary header       */
  /*       INPUTS: returnStatus, routine, msg                                */
  /*  ENDIF                                                                  */
  /*                                                                         */
  /***************************************************************************/

  packet_header->type = extr_bits(pkt, PH_PRI_TYPE_BIT_OFFSET,
     PH_PRI_TYPE_BYTE_OFFSET, PH_PRI_TYPE_NUM_BITS);

  if (packet_header->type != PH_PRI_TYPE_VALUE)
    {
     returnStatus = MODIS_E_INV_TYPE;
     sprintf(msg, "Type from primary header: %d Normal Type: 0 Test Packet: 1",
             packet_header->type);
     log_fmt_msg(MODIS_E_INV_TYPE, routine, msg);
    }


  /***************************************************************************/
  /*                                                                         */
  /*  CALL extr_bits to extract the secondary header flag from the packet    */
  /*     primary header                                                      */
  /*    INPUTS:  PH_PRI_SEC_HDR_FLAG_BYTE_OFFSET,                            */
  /*             PH_PRI_SEC_HDR_FLAG_BIT_OFFSET,                             */
  /*             PH_PRI_SEC_HDR_FLAG_NUM_BITS, pkt                           */
  /*    OUTPUTS: None                                                        */
  /*    RETURN:  PH_PACKET_HEADER_t.sec_hdr_flag                             */
  /*                                                                         */
  /*  IF PH_PACKET_HEADER_t.sec_hdr_flag is not equal to                     */
  /*     PH_PRI_SEC_HDR_PRESENT                                              */
  /*  THEN                                                                   */
  /*     set returnStatus equal to MODIS_E_INV_SEC_HDR_FLAG                  */
  /*     CALL log_fmt_msg to report that an invalid secondary header flag    */
  /*        value has been detected in the packet primary header             */
  /*       INPUTS: returnStatus, routine, msg                                */
  /*  ENDIF                                                                  */
  /*                                                                         */
  /***************************************************************************/

  packet_header->sec_hdr_flag = extr_bits(pkt, PH_PRI_SEC_HDR_FLAG_BIT_OFFSET,
     PH_PRI_SEC_HDR_FLAG_BYTE_OFFSET, PH_PRI_SEC_HDR_FLAG_NUM_BITS);

  if (packet_header->sec_hdr_flag != PH_PRI_SEC_HDR_PRESENT)
    {
     returnStatus = MODIS_E_INV_SEC_HDR_FLAG;
     sprintf(msg, "Secondary Header Flag value: %d Normal Value: %d", 
                  packet_header->sec_hdr_flag, PH_PRI_SEC_HDR_PRESENT);
     log_fmt_msg(MODIS_E_INV_SEC_HDR_FLAG, routine, msg);
    }



  /***************************************************************************/
  /*                                                                         */
  /*   CALL extr_bits to extract the application process ID (APID) from the  */ 
  /*      packet primary header                                              */
  /*     INPUTS:  pkt, PH_PRI_APID_BIT_OFFSET, PH_PRI_APID_BYTE_OFFSET,      */ 
  /*              PH_PRI_APID_NUM_BITS                                       */
  /*     OUTPUTS: None                                                       */
  /*     RETURN:  PH_PACKET_HEADER_t.apid                                    */
  /*                                                                         */
  /*   IF PH_PACKET_HEADER_t.apid is less than PH_PRI_MIN_MODIS_APID_AM1 OR  */
  /*      greater than PH_PRI_MAX_MODIS_APID_AM1                             */
  /*   THEN                                                                  */
  /*      set returnStatus equal to MODIS_E_INV_APID                         */
  /*      CALL log_fmt_msg to report that an invalid APID value has been     */
  /*         detected in the packet primary header                           */
  /*        INPUTS: returnStatus, routine, msg                               */
  /*   ELSE                                                                  */
  /*      IF PH_PACKET_HEADER_t.apid is equal to PH_PRI_APID_TEST_PACKET     */
  /*      THEN                                                               */
  /*         set returnStatus equal to MODIS_E_INV_APID_TEST                 */
  /*         set msg to "a test packet has been detected from information in */
  /*            the APID field of the packet primary header"                 */
  /*         CALL log_fmt_msg to report that a test packet has been detected */
  /*            from information in the APID field of the packet primary     */
  /*            header                                                       */
  /*           INPUTS: returnStatus, routine, msg                            */
  /*      ENDIF                                                              */
  /*   ENDIF                                                                 */
  /*                                                                         */
  /***************************************************************************/

  packet_header->apid = extr_bits(pkt, PH_PRI_APID_BIT_OFFSET,
     PH_PRI_APID_BYTE_OFFSET, PH_PRI_APID_NUM_BITS);

  if ((packet_header->apid < PH_PRI_MIN_MODIS_APID_AM1) ||
     (packet_header->apid > PH_PRI_MAX_MODIS_APID_AM1))
    {
     returnStatus = MODIS_E_INV_APID;
     sprintf(msg, "APID value: %d MODIS APID range 64-127", packet_header->apid);
     log_fmt_msg(MODIS_E_INV_APID, routine, msg);
    }
  else
    {
     if (packet_header->apid == PH_PRI_APID_TEST_PACKET)
       {
        returnStatus = MODIS_E_INV_APID_TEST;
        sprintf(msg, "APID Value: %d", packet_header->apid);
        log_fmt_msg(MODIS_E_INV_APID_TEST, routine, msg);
       }
    }


  /***************************************************************************/
  /*                                                                         */
  /*   CALL extr_bits to extract the sequence flag from the packet primary   */
  /*      header                                                             */
  /*     INPUTS:  pkt, PH_PRI_SEQUENCE_FLAG_BIT_OFFSET,                      */
  /*              PH_PRI_SEQUENCE_FLAG_BYTE_OFFSET,                          */
  /*              PH_PRI_SEQUENCE_FLAG_NUM_BITS                              */
  /*     OUTPUTS: None                                                       */
  /*     RETURN:  PH_PACKET_HEADER_t.sequence_flag                           */
  /*                                                                         */
  /*   IF PH_PACKET_HEADER_t.sequence_flag is equal to                       */
  /*      PH_PRI_SEQUENCE_NOT_USED                                           */
  /*   THEN                                                                  */
  /*      set returnStatus equal to MODIS_E_INV_SEQ_FLAG                     */
  /*      CALL log_fmt_msg to report that an invalid sequence flag value has */
  /*         been detected in the packet primary header                      */
  /*        INPUTS: returnStatus, routine, msg                               */
  /*   ENDIF                                                                 */
  /*                                                                         */
  /***************************************************************************/

  packet_header->sequence_flag = extr_bits(pkt, PH_PRI_SEQUENCE_FLAG_BIT_OFFSET,
     PH_PRI_SEQUENCE_FLAG_BYTE_OFFSET, PH_PRI_SEQUENCE_FLAG_NUM_BITS);

  if (packet_header->sequence_flag == PH_PRI_SEQUENCE_NOT_USED)
    {
     returnStatus = MODIS_E_INV_PKT_SEQ_FLAG;
     sprintf(msg, "sequence flag value: %d", packet_header->sequence_flag);
     log_fmt_msg(MODIS_E_INV_PKT_SEQ_FLAG, routine, msg);
    }
     

  /***************************************************************************/
  /*                                                                         */
  /*   CALL extr_bits to extract the source sequence count from the packet   */
  /*      primary header                                                     */
  /*     INPUTS:  pkt, PH_PRI_SOURCE_SEQ_CNT_BIT_OFFSET,                     */
  /*              PH_PRI_SOURCE_SEQ_CNT_BYTE_OFFSET,                         */
  /*              PH_PRI_SOURCE_SEQ_CNT_NUM_BITS,                            */ 
  /*     OUTPUTS: None                                                       */
  /*     RETURN:  PH_PACKET_HEADER_t.pkt_seq_cnt                             */
  /*                                                                         */
  /***************************************************************************/

  packet_header->pkt_seq_cnt = extr_bits(pkt, PH_PRI_SOURCE_SEQ_CNT_BIT_OFFSET,
     PH_PRI_SOURCE_SEQ_CNT_BYTE_OFFSET, PH_PRI_SOURCE_SEQ_CNT_NUM_BITS);



  /***************************************************************************/
  /*                                                                         */
  /*  CALL extr_bits to extract the packet length from the packet primary    */
  /*     header                                                              */
  /*    INPUTS:  pkt, PH_PRI_PKT_LENGTH_BIT_OFFSET,                          */
  /*             PH_PRI_PKT_LENGTH_BYTE_OFFSET, PH_PRI_PKT_LENGTH_NUM_BITS   */
  /*    OUTPUTS: None                                                        */
  /*    RETURN:  PH_PACKET_HEADER_t.pkt_length                               */
  /*                                                                         */
  /***************************************************************************/

  packet_header->pkt_length = extr_bits(pkt, PH_PRI_PKT_LENGTH_BIT_OFFSET,
     PH_PRI_PKT_LENGTH_BYTE_OFFSET, PH_PRI_PKT_LENGTH_NUM_BITS);


  /***************************************************************************/
  /*                                                                         */
  /*  RETURN returnStatus                                                    */
  /*                                                                         */
  /***************************************************************************/

  return (returnStatus);

}
