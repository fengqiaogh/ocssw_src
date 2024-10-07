#include "L1A_prototype.h"
#include "PGS_SMF.h"
#include "PGS_MODIS_35005.h"
#include "hdf.h"
#include "hdfi.h"
#include "VU_vdata_utility.h"

void  remember (char    *Vdata_name, 
                int32    Vdata_id )

/*
!C************************************************************************

!Description:  Stores the Vdata name/id pair in the id_table for future 
               reference.

!Input Parameters:
               char    *Vdata_name            ** The Vdata name to store  **
               int32   Vdata_id               ** The Vdata id to store    **

!Output Parameters:
               None

Return Values:
               None

Externally Defined:  
               PGS_SMF_status                    (PGS_SMF.h)
               int32                             (hdfi.h)
               VU_INCREMENT                      (VU_vdata_utility.h)
               MODIS_E_MAX_VDATAS_EXCEEDED       (MODIS_35003.t)
               MODIS_E_ID_TABLE_OVERFLOW         (MODIS_35003.t)
               FAIL                              (hdf.h)
               global_VU_VDATA_ID                (level1a)
               global_VU_ID_TABLE_READY          (level1a)

Called By:
               create_Vdata

Routines Called:
               initialize_id_table
               attached_Vdata_counter
               get_empty_slot
               log_fmt_msg

!Revision History:
               Revision 2.0  1997/10/01  20:45 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated Code.

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

!END************************************************************************
*/

 {
  /**************************************************************************/
  /*                                                                        */
  /*                  Declare the global variables.                         */
  /*                                                                        */
  /**************************************************************************/

  extern  VU_ID_TABLE  global_VU_VDATA_ID[VU_MAX_NUMBER_OF_VDATAS];

  extern   int   global_VU_ID_TABLE_READY;


  /**************************************************************************/
  /*                                                                        */
  /*          Declare the local variables and initialize them.              */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /* Set routine to "remember"                                              */
  /*                                                                        */
  /**************************************************************************/

  char  *routine = "remember";            /* Variable to hold routine name  */

  char  msg[300];                         /* Variable to hold error message */

  int16 number_attached_Vdatas;
  
  int16 empty_slot_index;


  /**************************************************************************/
  /*                                                                        */
  /* {if this is the first entry then initialize the id table}              */
  /* IF ( global_VU_ID_TABLE_READY equals False )                           */
  /* THEN                                                                   */
  /*    CALL initialize_id_table to prepare the id table for use            */
  /*      INPUTS:  None                                                     */
  /*      OUTPUTS: None                                                     */
  /*      RETURNS: None                                                     */
  /* ENDIF                                                                  */
  /*                                                                        */
  /* CALL attached_Vdata_counter to increment the number of attached Vdatas */
  /*   INPUTS:  VU_INCREMENT                                                */
  /*   OUTPUTS: None                                                        */
  /*   RETURNS: number_attached_Vdatas                                      */
  /*                                                                        */
  /* IF ( number_attached_Vdatas equals FAIL )                              */
  /* THEN                                                                   */
  /*    Set Status to MODIS_E_MAX_VDATAS_EXCEEDED                           */
  /*    set msg to "Program Error: Attempting to attach more than           */
  /*                VU_MAX_NUMBER_OF_VDATAS Vdatas to the file"             */
  /*    CALL log_fmt_msg to report the attempt to attach more that the      */
  /*         maximum number of Vdatas to the file                           */
  /*      INPUTS:  Status, routine, msg                                     */
  /*      OUTPUTS: None                                                     */
  /*      RETURNS: None                                                     */
  /* ENDIF                                                                  */
  /*                                                                        */
  /**************************************************************************/

     if (!global_VU_ID_TABLE_READY)
        initialize_id_table();

     number_attached_Vdatas = attached_Vdata_counter(VU_INCREMENT);
     if (number_attached_Vdatas == FAIL) {
        sprintf(msg, "Vdata Name: %s", Vdata_name);
        log_fmt_msg(MODIS_E_MAX_VDATAS_EXCEEDED, routine, msg);
     }



  /**************************************************************************/
  /*                                                                        */
  /* CALL get_empty_slot to get an empty id_table slot to use               */
  /*   INPUTS:  None                                                        */
  /*   OUTPUTS: None                                                        */
  /*   RETURNS: empty_slot_index                                            */
  /*                                                                        */
  /* IF ( empty_slot_index equals FAIL )                                    */
  /* THEN                                                                   */
  /*    set Status to MODIS_E_ID_TABLE_OVERFLOW                             */
  /*    set msg to "Program Error: Attempting to insert into a full id      */
  /*                table."                                                 */
  /*    CALL log_fmt_msg to report the attempt to insert into a full id     */
  /*         table                                                          */
  /*      INPUTS:  Status, routine, msg                                     */
  /*      OUTPUTS: None                                                     */
  /*      RETURNS: None                                                     */
  /* ELSE                                                                   */
  /*    set global_VU_VDATA_ID[empty_slot_index].id to Vdata_id             */
  /*    set global_VU_VDATA_ID[empty_slot_index].name to Vdata_name         */
  /* ENDIF                                                                  */
  /*                                                                        */
  /**************************************************************************/


     if ((empty_slot_index = get_empty_slot()) == FAIL){
        sprintf(msg, "No empty slots for Vdata %s exist", Vdata_name);
        log_fmt_msg(MODIS_E_ID_TABLE_OVERFLOW, routine, msg);
     }
     else {
        global_VU_VDATA_ID[empty_slot_index].id = Vdata_id;
        strcpy(global_VU_VDATA_ID[empty_slot_index].name, Vdata_name);
     }

 }  /* End of routine remember */
