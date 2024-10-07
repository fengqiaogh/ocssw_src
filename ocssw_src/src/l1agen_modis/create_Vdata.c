#include "L1A_prototype.h"
#include "PGS_SMF.h"
#include "PGS_MODIS_35005.h"
#include "VU_vdata_utility.h"
#include "hdfi.h"

#define	 WRITE	"w"

PGSt_SMF_status  create_Vdata(char  *Vdata_name, 
                              char  field_names[][VU_MAX_NAME_LENGTH], 
                              char  data_types[][VU_MAX_DATA_TYPE_STRING_LENGTH],
                              int16 num_fields,
                              uint16 order[])
/*
!C**************************************************************************

!Description:  Creates the specified Vdata in the currently open hdf file.

!Input Parameters:
               char  *Vdata_name              ** name of the Vdata          **
               char  field_names[][VU_MAX_NAME_LENGTH]
                                              ** an array of field names    **
                                              ** (each string's length is   **
                                              ** <= MAX_NAME_LENGTH)        **

               char  data_types[][VU_MAX_DATA_TYPE_STRING_LENGTH] 
               				      ** an array of data type      **
                                              ** strings which correspond to**
                                              ** to the field names (i.e.   **
                                              ** field_names[i]'s data type **
                                              ** is data_types[i]).         **

               int16  num_fields              ** the rank of the two arrays.** 

               uint16 order                   ** the number of data type    **
                                              ** size entries per record    **

!Output Parameters:
               None

!Input/Output Parameters:
               None

Return Values:
               MODIS_S_SUCCESS                 (PGS_MODIS_35005.h)
               FAIL                            (HDF)

Externally Defined:  
               PGSt_SMF_status                 (PGS_SMF.h)
               VU_MAX_NAME_LENGTH              (VU_vdata_utility.h)
               VU_MAX_DATA_TYPE_STRING_LENGTH  (VU_vdata_utility.h)
               VU_MAX_FIELD_NAME_LIST_SIZE     (VU_vdata_utility.h)
               VU_NEW_VDATA                    (VU_vdata_utility.h)
               int16                           (hdfi.h)
               int32                           (hdfi.h)
               MODIS_E_NULL_POINTER            (PGS_MODIS_35005.h)
               global_H_ID                     (level1a)
               MODIS_E_VSSETFIELDS             (PGS_MODIS_35005.h)

Called By:   
               init_L1A_HDF_vdatas

Routines Called:
               VSattach
               VSsetname
               VSsetfields
               remember
               create_Vdata_field
               log_fmt_msg

!Revision History:
               Revision 2.0  1998/10/26  09:58 EST
               John Seaton/SAIC/GSC (seaton@ltpmail.gsfc.nasa.gov)
               Allows for an array of orders to be passed in to create
               Vdata fields that have different orders.

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
	       HDF portions developed at the National Center for 
               Supercomputing Applications at the University of Illinois 
               at Urbana-Champaign.

!Design Notes: 
               This code was developed in C language.

!END************************************************************************
*/

 {
  /**************************************************************************/
  /*                                                                        */
  /*                  Declare the global variables.                         */
  /*                                                                        */
  /**************************************************************************/

  extern   int32   global_H_ID;  


  /**************************************************************************/
  /*                                                                        */
  /*          Declare the local variables and initialize them.              */
  /*                                                                        */
  /**************************************************************************/
  /* declare field_name_list to be an array of characters                   */
  /*         VU_MAX_FIELD_NAME_LIST_SIZE long                               */
  /*                                                                        */
  /* declare Vdata_id to be a variable of type int32                        */
  /*                                                                        */
  /* Set routine to "create_Vdata"                                          */
  /*                                                                        */
  /* Set access_mode to "w" for write access                                */
  /*                                                                        */
  /* Set returnStatus to MODIS_S_SUCCESS                                    */
  /*                                                                        */
  /**************************************************************************/

  PGSt_SMF_status  returnStatus;  /* SMF-style message returned by function */

  PGSt_SMF_status  status;        /* SMF-style message returned by function */

  char  field_name_list[VU_MAX_FIELD_NAME_LIST_SIZE];
                                  /* List to hold the VSsetfields           */
  char  msg[300];                 /* Array to hold error messages           */

  int32 Vdata_id;                 /* Vdata access indentifier returned from */
                                  /* routine VSattach                       */

  char  *routine;                 /* Variable to hold routine name          */

  char  *access_mode;             /* Mode in which the file is accessed     */

  int  i;                         /* Loop Variable                          */
  unsigned long cur_field_name_list_size=0; /* stores size of field array   */

  returnStatus    = MODIS_S_SUCCESS;
  routine         = "create_Vdata";
  access_mode     = WRITE;

  /**************************************************************************/
  /*                                                                        */
  /* IF Vdata_name == NULL                                                  */
  /* THEN                                                                   */
  /*   CALL log_fmt_msg to report that the new Vdata name passed in is NULL.*/
  /*     INPUTS: MODIS_E_NULL_POINTER, routine, " "                         */
  /*     OUTPUTS: None                                                      */
  /*     RETURNS: None                                                      */
  /*                                                                        */
  /*   return FAIL                                                          */
  /* ENDIF                                                                  */
  /*                                                                        */
  /**************************************************************************/

  if (Vdata_name == NULL) {
    log_fmt_msg(MODIS_E_NULL_POINTER, routine, " ");
    return FAIL;
  }

  /**************************************************************************/
  /*                                                                        */
  /* CALL VSattach to attach the new Vdata to the file                      */
  /*   INPUTS:  global_H_ID, VU_NEW_VDATA, access_mode                      */
  /*   OUTPUTS: None                                                        */
  /*   RETURNS: Vdata_id                                                    */
  /*                                                                        */
  /**************************************************************************/

     Vdata_id = VSattach(global_H_ID, VU_NEW_VDATA, access_mode);


  /**************************************************************************/
  /*                                                                        */
  /* IF (Vdata_id equals FAIL)                                              */
  /* THEN                                                                   */
  /*    CALL log_fmt_msg to report that the new Vdata could not be attached */
  /*         to the file                                                    */
  /*      INPUTS:  MODIS_E_ATTACHED_VDATAS, routine, "unable to attach the  */
  /*               new Vdata to the file"                                   */
  /*      OUTPUTS: None                                                     */
  /*      RETURNS: None                                                     */
  /*                                                                        */
  /*    set returnStatus to FAIL                                            */
  /*                                                                        */
  /**************************************************************************/

     if (Vdata_id == FAIL) {
        sprintf(msg, "Vdata Name: %s", Vdata_name);
        log_fmt_msg(MODIS_E_CREATE_VDATA, routine, msg);
        returnStatus = FAIL;
     }


  /**************************************************************************/
  /*                                                                        */
  /* ELSE                                                                   */
  /*    CALL VSsetname to set the Vdata name                                */
  /*      INPUTS:  Vdata_id, Vdata_name                                     */
  /*      OUTPUTS: None                                                     */
  /*      RETURNS: None                                                     */
  /*                                                                        */
  /*    CALL remember to store the Vdata's name and id as an associated pair*/
  /*        INPUTS:  Vdata_name, Vdata_id                                   */
  /*        OUTPUTS: None                                                   */
  /*        RETURNS: None                                                   */
  /*                                                                        */
  /**************************************************************************/

     else {
        VSsetname(Vdata_id, Vdata_name);
        remember(Vdata_name, Vdata_id);


  /**************************************************************************/
  /*    initialize array to hold field names to zeros                       */
  /*    initialize the counter for field name array overrun to 0            */
  /*    FOR (i = 0;i is less than num_fields) AND (returnStatus is equal to */
  /*     MODIS_S_SUCCESS); i++)                                             */
  /*       CALL create_Vdata_field to create a field within this Vdata      */
  /*         INPUTS:  Vdata_name, Vdata_id, field_names[i], data_types[i],  */
  /*                  order[i]                                              */
  /*         OUTPUTS: None                                                  */
  /*         RETURNS: returnStatus                                          */
  /*                                                                        */
  /*       IF field name array overflows it bounds                          */
  /*       THEN                                                             */
  /*         set returnStatus = FAIL                                        */
  /*         set msg to hold the Vdata name that caused failure             */
  /*                                                                        */
  /*         CALL log_fmt_msg to record the error in LogStatus file         */
  /*          INPUTS: MODIS_E_FIELD_NAME_LIST_OVERRUN, routine, msg         */
  /*          OUTPUTS: None                                                 */
  /*          RETURNS: None                                                 */
  /*                                                                        */
  /*       ELSE                                                             */
  /*         Append a field name (field_names[i]) to the list of field names*/
  /*           (field_name_list)                                            */
  /*                                                                        */
  /*         IF (i is less than (num_fields-1))                             */
  /*         THEN                                                           */
  /*            Append a comma to the list of field names (field_name_list) */
  /*         ENDIF                                                          */
  /*       ENDIF                                                            */
  /*    END FOR                                                             */
  /*                                                                        */
  /**************************************************************************/

        memset(field_name_list, '\0', sizeof(field_name_list));
        cur_field_name_list_size = 0;
        for (i=0; ((i < num_fields) && (returnStatus == MODIS_S_SUCCESS)); i++) {
           returnStatus = create_Vdata_field(Vdata_name,
                                             Vdata_id,
                                             field_names[i],
                                             data_types[i],
                                             (int32) order[i]);

           cur_field_name_list_size += strlen(field_names[i]);
           if (cur_field_name_list_size > VU_MAX_FIELD_NAME_LIST_SIZE-1) {
             returnStatus = MODIS_E_FIELD_NAME_LIST_OVERRUN;
             sprintf(msg, "Vdata name = %s", Vdata_name);
             log_fmt_msg(MODIS_E_FIELD_NAME_LIST_OVERRUN, routine, msg);

           }
           else {
             strcat(field_name_list, field_names[i]);
             if (i < num_fields-1)
                strcat(field_name_list, ",");
           }
        }


  /**************************************************************************/
  /*                                                                        */
  /*    IF (returnStatus == FAIL)                                           */
  /*    THEN                                                                */
  /*      CALL log_fmt_msg to report that create_Vdata_field failed.        */
  /*        INPUTS: MODIS_E_CREATE_VDATA_FIELD, routine, " "                */
  /*        OUTPUTS: None                                                   */
  /*        RETURNS: None                                                   */
  /*    ENDIF                                                               */
  /*                                                                        */
  /*    IF (returnStatus == MODIS_E_FIELD_NAME_LIST_OVERRUN)                */
  /*    THEN                                                                */
  /*      set returnStatus equal to FAIL                                    */
  /*    ENDIF                                                               */
  /*                                                                        */
  /*    field_name_list is a comma separated list of field names            */
  /*       (e.g. "X,Y,Z")                                                   */
  /*                                                                        */
  /*    CALL VSsetfields to set the fields of the Vdata                     */
  /*         INPUTS:  Vdata_id, field_name_list                             */
  /*         OUTPUTS: None                                                  */
  /*         RETURNS: Status                                                */
  /*                                                                        */
  /*    IF ( Status equals FAIL )                                           */
  /*    THEN                                                                */
  /*       CALL log_fmt_msg to report that the fields of the Vdata could    */
  /*             not set                                                    */
  /*         INPUTS:  MODIS_E_VSSETFIELDS, routine, "unable to set the      */
  /*                  fields of the Vdata"                                  */
  /*         OUTPUTS: None                                                  */
  /*         RETURNS: None                                                  */
  /*       set returnStatus to FAIL                                         */
  /*    ENDIF                                                               */
  /* ENDIF                                                                  */
  /*                                                                        */
  /**************************************************************************/

        if (returnStatus == FAIL) {
          sprintf(msg, "Vdata Name: %s", Vdata_name);
          log_fmt_msg(MODIS_E_CREATE_VDATA_FIELD, routine, msg);
        }

        if (returnStatus == MODIS_E_FIELD_NAME_LIST_OVERRUN)
           returnStatus = FAIL;

        status = VSsetfields(Vdata_id, field_name_list);

        if (status == FAIL) {
           sprintf(msg, "Vdata Name: %s", Vdata_name);
           log_fmt_msg(MODIS_E_VSSETFIELDS, routine, msg);
           returnStatus = FAIL;
        }
     }

  /**************************************************************************/
  /*                                                                        */
  /* RETURN returnStatus                                                    */
  /*                                                                        */
  /**************************************************************************/

     return returnStatus;

 }   /* End of create_Vdata */
