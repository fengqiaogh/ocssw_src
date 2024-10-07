#include "L1A_prototype.h"
#include "PGS_MODIS_35005.h"
#include "PGS_IO_L0.h"
#include "PGS_SMF.h"
#include "MS_misc.h"
#include "PC_pcf_info.h"


unsigned int bytetoint (void *bitstream, unsigned int length);

PGSt_SMF_status validate_L0_header (PGSt_IO_L0_VirtualDataSet L0_file)

/*
!C******************************************************************************

!Description:  Function validate_L0_header gets the L0 header and validate 
               the contents of the header.

!Input Parameters: 
               PGSt_IO_L0_VirtualDataSet  L0_file    **  L0 file descriptor  **

!Output Parameters: 
               None

Return Values: 
               MODIS_S_SUCCESS                          (PGS_MODIS_35005.h)
               MODIS_F_L0_HEADER_VAL_FAILED             (PGS_MODIS_35005.h)

Externally Defined: 
               PGSt_IO_L0_VirtualDataSet                (PGS_IO_L0.h)
               MS_HEADER_BUFF_SIZE                      (MS_misc.h)
               MS_FOOTER_BUFF_SIZE                      (MS_misc.h)
               PC_INSTRUMENT                            (PC_pcf_info.h)
               AQUA_SCID                                (PC_pcf_info.h)
               TERRA_SCID                               (PC_pcf_info.h)
               INSTRUMENT_AQUA                          (PC_pcf_info.h)
               NSTRUMENT_TERRA                          (PC_pcf_info.h)

Called By:
               get_valid_L0_file
               read_a_packet
               set_start_position
  
Routines Called:
               PGS_IO_L0_GetHeader
               PGS_SMF_TestSuccessLevel
               log_fmt_msg
               PGS_PC_GetConfigData

!Revision History:    
$Log: validate_L0_header.c,v $
Revision 5.2  2004/09/29 19:32:53  seaton
Fixed bug found while unit testing.

Revision 5.1  2004/09/23 19:06:00  seaton
Added code to check the L0 construction record SCID to make sure we are processing the
correct spacecraft L0 data. This is a collection 5 modification.
seaton@saicmodis.com

               revision 3.0
               John Seaton
               Added logic to not fail if the CR header buffer overflows.

               revision 1.0 1997/09/19  17:30:00
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

!END***************************************************************************
*/

{
  /***************************************************************************/
  /*                                                                         */
  /*               Declare and Initialize Local Variables                    */
  /*                                                                         */
  /***************************************************************************/

  char                  *routine = "validate_L0_header";
  PGSt_SMF_status       returnStatus=MODIS_S_SUCCESS;
  PGSt_IO_L0_Header     head_buff[MS_HEADER_BUFF_SIZE];
  PGSt_IO_L0_Footer     foot_buff[MS_FOOTER_BUFF_SIZE];
  PGSt_SMF_status       PGS_status;
  char                  buffer[PGSd_PC_VALUE_LENGTH_MAX];
  int                   scid_val;
  char                  msg[300];
  int                   offset=0;
  int                   value,i;

  /***************************************************************************/
  /*                                                                         */
  /*  CALL PGS_IO_L0_GetHeader to get the L0 header                          */
  /*    INPUT:  L0_file, header_buff_size, foot_buff_size                    */
  /*    OUTPUT: head_buff, foot_buff                                         */
  /*    RETURN: PGS_status                                                   */
  /*                                                                         */
  /***************************************************************************/

  PGS_status = PGS_IO_L0_GetHeader(L0_file,MS_HEADER_BUFF_SIZE,head_buff,
                                   MS_FOOTER_BUFF_SIZE,foot_buff);


  /***************************************************************************/
  /*                                                                         */
  /*  CALL PGS_SMF_TestSuccessLevel to determine if getting L0 header was    */
  /*     successful                                                          */
  /*    INPUT:  PGS_status                                                   */
  /*    ONTPUT: None                                                         */
  /*    RETURN: TestStatus                                                   */
  /*                                                                         */
  /*  IF TestStatus is not equal to PGS_TRUE                                 */
  /*  THEN                                                                   */
  /*    IF PGS_status not warning that header buffer was truncated           */
  /*    THEN                                                                 */
  /*      Set returnStatus to MODIS_F_L0_HEADER_VAL_FAILED                   */
  /*      CALL log_fmt_msg to report that getting L0 header was not successful*/
  /*        INPUT:  PGS_status, routine, msg                                 */
  /*        OUTPUT: None                                                     */
  /*        RETURN: None                                                     */
  /*    ENDIF                                                                */
  /*  ENDIF                                                                  */
  /*                                                                         */
  /***************************************************************************/

  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      if (PGS_status != PGSIO_W_L0_HDR_BUF_TRUNCATE) { 
        returnStatus = MODIS_F_L0_HEADER_VAL_FAILED;
        log_fmt_msg(MODIS_F_L0_HEADER_VAL_FAILED, routine, 
                  "L0 header was not successfully retrieved");
      }
  } else {
    PGS_status = PGS_PC_GetConfigData(PC_INSTRUMENT, buffer);
                                                                                                                                      
    if (!PGS_SMF_TestSuccessLevel(PGS_status)) {
      returnStatus = MODIS_E_GETCONFIG_FAILED;
      log_fmt_msg(MODIS_E_GETCONFIG_FAILED, routine, 
                  "The satellite instrument string could not be retrieved from the pcf file");
    }
    else {
      /* parse cr to get scid here */
      /* See CDRL B301 for Layout of the L0 Construction Record file */

       offset += 50;
                                                                                                                                      
       value = bytetoint(&head_buff[offset], 2);
       offset+=2;
                                                                                                                                      
       for(i=0; i<value; i++)
         offset += 16;
                                                                                                                                      
       offset += 81;
       scid_val = bytetoint(&head_buff[offset], 1);

       if ((!((scid_val == TERRA_SCID) && (strcmp(buffer,INSTRUMENT_TERRA) == 0))) &&
           (!((scid_val == AQUA_SCID) && (strcmp(buffer,INSTRUMENT_AQUA) == 0)))) {
         returnStatus = MODIS_F_L0_HEADER_VAL_FAILED;
         sprintf(msg, "SCID from Construction Record %d does not match PC_INSTRUMENT tag %s in the pcf file", scid_val, buffer); 
         log_fmt_msg(MODIS_F_L0_HEADER_VAL_FAILED, routine, msg);
       }
    }
  }


  /***************************************************************************/
  /*                                                                         */
  /*  RETURN returnStatus                                                    */
  /*                                                                         */
  /***************************************************************************/

  return (returnStatus);

} /* End of routine validate_L0_header */

unsigned int bytetoint (void *bitstream, unsigned int length) {
                                                                                                                                      
  unsigned char  *bitPtr;
  unsigned int number=0;
                                                                                                                                      
  bitPtr = ((unsigned char*) bitstream);
                                                                                                                                      
  while (length--)
        number = number*256 + *(bitPtr++);
                                                                                                                                      
  return number;
                                                                                                                                      
}

