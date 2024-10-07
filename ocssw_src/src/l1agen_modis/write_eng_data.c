#include "L1A_prototype.h"
#include "EN_eng_data.h"
#include "PGS_MODIS_35005.h"
#include "PGS_SMF.h"
#include "hdfi.h"
#include "hdf.h"
#include "hntdefs.h"


PGSt_SMF_status  write_eng_data (EN_VDATA_TYPE_t  *eng_data)

/*
!C****************************************************************************

!Description:  This function writes the eng data values for the current scan 
               to the L1A file .

!Input Parameters:
               EN_VDATA_TYPE_t  *eng_data   ** The eng_data array structure **

!Output Parameters:
               None

Return Values:
               MODIS_S_SUCCESS              (PGS_MODIS_35005.h)
               MODIS_F_WRITE_ENG_DATA_FAIL  (PGS_MODIS_35005.h)

Externally Defined:  
               EN_VDATA_TYPE_t               (EN_eng_data.h)
               EN_NUM_VDATAS                 (EN_eng_data.h)
               EN_MAX_FIELDS_PER_VDATA       (EN_eng_data.h)
               PGSt_SMF_status               (PGS_SMF.h)
               MODIS_E_WRITE_VDATA           (PGS_MODIS_35005.h)
               EN_SC_ANCILLARY_VDATA_START   (EN_eng_data.h)
               EN_SC_ANCILLARY_VDATA_STOP    (EN_eng_data.h)
               EN_ANCIL_VDATA_BUFFER_SIZE    (EN_eng_data.h)
               MODIS_E_NULL_POINTER          (PGS_MODIS_35005.h)
               MODIS_E_VDATA_BUFFER_OVERFLOW (PGS_MODIS_35005.h)

Called By:
               write_scan
               handle_missing_scans

Routines Called:
               write_Vdata
               log_fmt_msg

!Revision History:
               Revision 2.0  1998/10/26  11:05  EST
               John Seaton/SAIC/GSC  (seaton@ltpmail.gsfc.nasa.gov)
               Code added to write the Current/Prior S/C Ancillary
               Vdatas, which contain multiple size and order data.

               Revision 1.1  1997/09/03  10:55
               Tom Johnson/GSC     (johnson@ltpmail.gsfc.nasa.gov)
               Incorporate walkthrough comments

               Revision 1.0  1997/07/14  15:58 EDT
               David Catozzi/SAIC/GSC (cato@ltpmail.gsfc.nasa.gov)
               Original design.

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
  /*                      Declare Local Variables                            */
  /*                                                                         */
  /*     ---------------------------------------------------------------     */
  /*                                                                         */
  /*     Allocate record_values[EN_MAX_FIELDS_PER_VDATA] as temporary        */
  /*     space to hold the field values to be written for the current        */
  /*     Vdata (eng_data[i])                                                 */
  /*                                                                         */
  /***************************************************************************/

  int i, j, k; 

  char             *routine = "write_eng_data";

  char             msg[300];

  PGSt_SMF_status  returnStatus; 

  PGSt_SMF_status  tempStatus; 

  uint16 record_values[EN_MAX_FIELDS_PER_VDATA];

  unsigned char *record = (unsigned char *) record_values;  

  unsigned char buffer[EN_ANCIL_VDATA_BUFFER_SIZE]="";

  unsigned char *p;

  int cur_pos_in_buf=0;



  /***************************************************************************/
  /*                                                                         */
  /*  Set routine to "write_eng_data"  (done during declaration)             */
  /*                                                                         */
  /*  Set returnStatus to MODIS_S_SUCCESS                                    */
  /*                                                                         */
  /***************************************************************************/
 
  returnStatus = MODIS_S_SUCCESS;

  /***************************************************************************/
  /* if input parameters equal NULL                                          */
  /* THEN                                                                    */
  /*   CALL log_fmt_msg to record error in LogStatus file                    */
  /*     INPUTS: MODIS_E_NULL_POINTER, routine, " "                          */
  /*     OUTPUTS: None                                                       */
  /*     RETURNS: None                                                       */
  /*                                                                         */
  /*     return MODIS_F_WRITE_ENG_DATA_FAIL                                  */
  /* ENDIF                                                                   */
  /***************************************************************************/

  if (eng_data == NULL) {
    log_fmt_msg(MODIS_E_NULL_POINTER, routine, " ");
    return MODIS_F_WRITE_ENG_DATA_FAIL;
  }


  /***************************************************************************/
  /*                                                                         */
  /*  DO FOR ( i = 0 to EN_NUM_VDATAS-1 )                                    */
  /*                                                                         */
  /*    IF i = EN_SC_ANCILLARY_VDATA_START or i = EN_SC_ANCILLARY_VDATA_END  */
  /*    THEN                                                                 */
  /*      initialize unsigned char buffer and a pointer into that buffer     */
  /*                                                                         */
  /*      initialize cur_pos_in_buf equal to zero                            */
  /*                                                                         */
  /*      FOR each field in the Vdata (j)                                    */
  /*                                                                         */
  /*        FOR each order within the field (k)                              */
  /*                                                                         */
  /*          IF cur_pos_in_buf > EN_ANCIL_DATA_BUFFER_SIZE                  */
  /*          THEN                                                           */
  /*            CALL log_fmt_msg to record error inthe LogStatus file        */
  /*              INPUTS: MODIS_E_VDATA_BUFFER_OVERFLOW, routine, " "        */
  /*              OUTPUTS: None                                              */
  /*              RETURNS: None                                              */
  /*                                                                         */
  /*            return MODIS_F_WRITE_ENG_DATA_FAIL                           */
  /*          ENDIF                                                          */
  /*                                                                         */
  /*          SWITCH on eng_data[i].field[j].type                            */
  /*                                                                         */
  /*             CASE DFNT_INT8 : copy data from                             */
  /*                              eng_data[i].fields[j].union_value.i8type   */
  /*                              into the buffer                            */
  /*                                                                         */
  /*                     increment the pointer into the buffer by the size of*/ 
  /*                      int8 type                                          */
  /*                                                                         */
  /*                     increment cur_pos_in_buf by size of int8 type       */
  /*                                                                         */
  /*                     break CASE                                          */
  /*                                                                         */
  /*            CASE DFNT_UINT8 : copy data from                             */
  /*                              eng_data[i].fields[j].union_value.ui8type  */
  /*                              into the buffer                            */
  /*                                                                         */
  /*                     increment the pointer into the buffer by the size of*/ 
  /*                      uint8 type                                         */
  /*                                                                         */
  /*                     increment cur_pos_in_buf by size of uint8 type      */
  /*                                                                         */
  /*                     break CASE                                          */
  /*                                                                         */
  /*            CASE DFNT_INT16 : copy data from                             */
  /*                              eng_data[i].fields[j].union_value.i16type  */
  /*                              into the buffer                            */
  /*                                                                         */
  /*                     increment cur_pos_in_buf by size of int16 type      */
  /*                                                                         */
  /*                     increment the pointer into the buffer by the size of*/ 
  /*                      int16 type                                         */
  /*                                                                         */
  /*                     break CASE                                          */
  /*                                                                         */
  /*            CASE DFNT_UINT16 : copy data from                            */
  /*                               eng_data[i].fields[j].union_value.ui16type*/ 
  /*                               into the buffer                           */
  /*                                                                         */
  /*                     increment cur_pos_in_buf by size of uint16 type     */
  /*                                                                         */
  /*                     increment the pointer into the buffer by the size of*/ 
  /*                      uint16 type                                        */
  /*                                                                         */
  /*                     break CASE                                          */
  /*                                                                         */
  /*            CASE DFNT_INT32 : copy data from                             */
  /*                              eng_data[i].fields[j].union_value.i32type  */
  /*                              into the buffer                            */
  /*                                                                         */
  /*                     increment cur_pos_in_buf by size of int32 type      */
  /*                                                                         */
  /*                     increment the pointer into the buffer by the size of*/ 
  /*                      int32 type                                         */
  /*                                                                         */
  /*                     break CASE                                          */
  /*                                                                         */
  /*            CASE DFNT_UINT32 : copy data from                            */
  /*                               eng_data[i].fields[j].union_value.ui32type*/ 
  /*                               into the buffer                           */
  /*                                                                         */
  /*                     increment cur_pos_in_buf by size of uint32 type     */
  /*                                                                         */
  /*                     increment the pointer into the buffer by the size of*/ 
  /*                      uint32 type                                        */
  /*                                                                         */
  /*                     break CASE                                          */
  /*                                                                         */
  /*             CASE default : CALL log_fmt_msg to record error in LogStatus*/
  /*                             INPUTS: MODIS_E_INVALID_VDATA_TYPE, routine,*/
  /*                                     " "                                 */
  /*                             OUTPUTS: None                               */
  /*                             RETURNS: None                               */
  /*                                                                         */
  /*                      break CASE                                         */
  /*                                                                         */
  /*          END SWITCH                                                     */
  /*                                                                         */
  /*        END FOR (k)                                                      */
  /*                                                                         */
  /*      END FOR (j)                                                        */
  /***************************************************************************/
 
  for (i = 0; i < EN_NUM_VDATAS; i++)
    {
     if ((i == EN_SC_ANCILLARY_VDATA_START) || (i == EN_SC_ANCILLARY_VDATA_END)) {
       memset(buffer, 0, sizeof(buffer));
       p = buffer;
       cur_pos_in_buf = 0;
       for (j = 0; j < eng_data[i].num_fields; j++) {
         for (k = 0; k < eng_data[i].field[j].order; k++ ) {
           if (cur_pos_in_buf > EN_ANCIL_VDATA_BUFFER_SIZE) {
             sprintf(msg, "Current Position: %d Max Buffer Size: %d", 
                     cur_pos_in_buf, EN_ANCIL_VDATA_BUFFER_SIZE);
             log_fmt_msg(MODIS_E_VDATA_BUFFER_OVERFLOW, routine, msg);
             return MODIS_F_WRITE_ENG_DATA_FAIL;
           }

           switch (eng_data[i].field[j].type) {
             case DFNT_INT8   :
                       memcpy(p, (const void *)&eng_data[i].field[j].union_value[k].i8type,
                                  sizeof(eng_data[i].field[j].union_value[k].i8type));
                       p += sizeof(eng_data[i].field[j].union_value[k].i8type);
                       cur_pos_in_buf += sizeof(eng_data[i].field[j].union_value[k].i8type);
                       break;
             case DFNT_UINT8   :
                       memcpy(p, (const void *)&eng_data[i].field[j].union_value[k].ui8type,
                                  sizeof(eng_data[i].field[j].union_value[k].ui8type));
                       p += sizeof(eng_data[i].field[j].union_value[k].ui8type);
                       cur_pos_in_buf += sizeof(eng_data[i].field[j].union_value[k].ui8type);
                       break;
             case DFNT_INT16   :
                       memcpy(p, (const void *)&eng_data[i].field[j].union_value[k].i16type,
                                  sizeof(eng_data[i].field[j].union_value[k].i16type));
                       p += sizeof(eng_data[i].field[j].union_value[k].i16type);
                       cur_pos_in_buf += sizeof(eng_data[i].field[j].union_value[k].i16type);
                       break;
             case DFNT_UINT16    :
                       memcpy(p, (const void *)&eng_data[i].field[j].union_value[k].ui16type,
                                  sizeof(eng_data[i].field[j].union_value[k].ui16type));
                       p += sizeof(eng_data[i].field[j].union_value[k].ui16type);
                       cur_pos_in_buf += sizeof(eng_data[i].field[j].union_value[k].ui16type);
                       break;
             case DFNT_INT32    :
                       memcpy(p, (const void *)&eng_data[i].field[j].union_value[k].i32type,
                                  sizeof(eng_data[i].field[j].union_value[k].i32type));
                       p += sizeof(eng_data[i].field[j].union_value[k].i32type);
                       cur_pos_in_buf += sizeof(eng_data[i].field[j].union_value[k].i32type);
                       break;
             case DFNT_UINT32    :
                       memcpy(p, (const void *)&eng_data[i].field[j].union_value[k].ui32type,
                                  sizeof(eng_data[i].field[j].union_value[k].ui32type));
                       p += sizeof(eng_data[i].field[j].union_value[k].ui32type);
                       cur_pos_in_buf += sizeof(eng_data[i].field[j].union_value[k].ui32type);
                       break;
             default             :
                       sprintf(msg, "Engineering Data Field Type: %d Valid Range: %d to %d", 
                         eng_data[i].field[j].type, EN_MIN_VDATA_TYPE, EN_MAX_VDATA_TYPE);
                       log_fmt_msg(MODIS_E_INVALID_VDATA_TYPE, routine, msg);
                       break;
           }  /* switch */

         } /* end for k */
        
       } /* for j */

  /***************************************************************************/
  /*      CALL write_Vdata to update all fields (write a single consolidated */
  /*          record of different data types for the current Vdata to the    */   
  /*          file) for this scan                                            */
  /*        INPUTS:  eng_data[i].vdata_name, buffer, 1 (Number of Records)   */
  /*        OUTPUTS: None                                                    */
  /*        RETURN:  tempStatus                                              */
  /***************************************************************************/

       tempStatus = write_Vdata(eng_data[i].vdata_name, (unsigned char *)buffer, 1);

  /***************************************************************************/
  /*      IF ( tempStatus is FAIL )                                          */  
  /*      THEN                                                               */
  /*        set Status to MODIS_E_WRITE_VDATA                                */
  /*        CALL log_fmt_msg to report that the fields for a Vdata could not */
  /*           be written to a scan                                          */ 
  /*          INPUTS:  Status, routine, msg                                  */
  /*          OUTPUTS: None                                                  */
  /*          RETURN:  None                                                  */
  /*        set returnStatus to FAIL                                         */
  /*      ENDIF                                                              */  
  /***************************************************************************/

       if (tempStatus == FAIL)
         {
          sprintf(msg, "Vdata Name = %s", eng_data[i].vdata_name);
          log_fmt_msg(MODIS_E_WRITE_VDATA, routine, msg);
          returnStatus = MODIS_F_WRITE_ENG_DATA_FAIL;
         }

     }

  /***************************************************************************/
  /*   ELSE                                                                  */
  /*                                                                         */ 
  /*     set all elements of record_values to 0                              */
  /*                                                                         */
  /***************************************************************************/

     else {
       memset(record_values, 0, sizeof(record_values));


  /***************************************************************************/
  /*                                                                         */
  /*     DO FOR ( j = 0 to eng_data[i].num_fields-1 )                        */
  /*        set record_values[j] to eng_data[i].field[j].value               */
  /*     END DO                                                              */
  /*     (To consolidate a record (a set of field values) for writing.)      */
  /*                                                                         */
  /***************************************************************************/

       for (j = 0; j < eng_data[i].num_fields; j++)
          record_values[j] = eng_data[i].field[j].value;


  /***************************************************************************/
  /*                                                                         */
  /*     CALL write_Vdata to update all fields (write a single consolidated  */
  /*        record for the current Vdata to the file) for this scan          */
  /*       INPUTS:  eng_data[i].vdata_name, record_values, 1                 */
  /*       OUTPUTS: None                                                     */
  /*       RETURN:  tempStatus                                               */
  /*                                                                         */
  /*     IF (tempStatus is FAIL)                                             */
  /*     THEN                                                                */
  /*        set Status to MODIS_E_WRITE_VDATA                                */
  /*        CALL log_fmt_msg to report that the fields for a Vdata could     */
  /*           not be written to a scan                                      */
  /*          INPUTS:  Status, routine, msg                                  */
  /*          OUTPUTS: None                                                  */
  /*          RETURN:  None                                                  */
  /*        set returnStatus to MODIS_F_WRITE_ENG_DATA_FAIL                  */
  /*     ENDIF                                                               */
  /*                                                                         */
  /*  END DO                                                                 */
  /*                                                                         */
  /***************************************************************************/

       tempStatus = write_Vdata(eng_data[i].vdata_name, record, 1);
       if (tempStatus == FAIL)
         {
          sprintf(msg, "Vdata Name = %s", eng_data[i].vdata_name);
          log_fmt_msg(MODIS_E_WRITE_VDATA, routine, msg);
          returnStatus = MODIS_F_WRITE_ENG_DATA_FAIL;
         }
     }
 
    }


  /***************************************************************************/
  /*                                                                         */
  /*  RETURN returnStatus                                                    */
  /*                                                                         */
  /***************************************************************************/

  return (returnStatus);

}

