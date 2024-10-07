#include "L1A_prototype.h"
#include "PGS_SMF.h"
#include "PGS_MODIS_35005.h"
#include "hdf.h"
#include "hdfi.h"
#include "VU_vdata_utility.h"

void  forget (char  *Vdata_name)

/*
!C************************************************************************

!Description:  Removes the Vdata name/id pair from the id_table.

!Input Parameters:
               char  *Vdata_name            ** The Vdata name to remove **

!Output Parameters:
               None

Return Values:
               None

Externally Defined:  
               int16                                  (hdfi.h)
               VU_EMPTY                               (VU_vdata_utility.h)
               VU_DECREMENT                           (VU_vdata_utility.h)
               MODIS_W_ID_TABLE_UNINITIALIZED         (PGS_MODIS_35005.t)
               MODIS_W_VDATA_NOT_FOUND_IN_TABLE       (PGS_MODIS_35005.t)
               MODIS_W_NEGATIVE_VDATA_COUNTER         (PGS_MODIS_35005.t)
               FAIL                                   (hdf.h)
               global_VU_VDATA_ID                     (level1a)
               global_VU_ID_TABLE_READY               (level1a)

Called By:
               close_Vdata

Routines Called:
               log_fmt_msg
               get_index
               attached_Vdata_counter

!Revision History:
               Revision 2.0  1997/10/01  21:23 EDT
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

  extern VU_ID_TABLE  global_VU_VDATA_ID[VU_MAX_NUMBER_OF_VDATAS];

  extern   int   global_VU_ID_TABLE_READY;


  /**************************************************************************/
  /*                                                                        */
  /*          Declare the local variables and initialize them.              */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /* Set routine to "forget"                                                */
  /*                                                                        */
  /**************************************************************************/

  char  *routine = "forget";              /* Variable to hold routine name  */

  char  msg[300];                         /* Variable to hold error message */

  int16 this_Vdata_index;

  int16 number_attached_Vdatas;



  /**************************************************************************/
  /*                                                                        */
  /* {first ensure that the id table has been initialized}                  */
  /* IF ( global_VU_ID_TABLE_READY equals False )                           */
  /* THEN                                                                   */
  /*    CALL log_fmt_msg to report that the id table was not initialized    */
  /*      INPUTS:  Status, routine, msg                                     */
  /*      OUTPUTS: None                                                     */
  /*      RETURNS: None                                                     */
  /*                                                                        */
  /**************************************************************************/

     if (!global_VU_ID_TABLE_READY)
        log_fmt_msg(MODIS_W_ID_TABLE_NOT_INIT, routine, 
                  "Warning: forget called before the id_table was initialized");


  /**************************************************************************/
  /*                                                                        */
  /* ELSE                                                                   */
  /*    CALL  get_index to get the Vdata index                              */
  /*      INPUTS:  Vdata_name                                               */
  /*      OUTPUTS: None                                                     */
  /*      RETURNS: this_Vdata_index                                         */
  /*                                                                        */
  /*    IF ( this_Vdata_index equals FAIL )                                 */
  /*    THEN                                                                */
  /*       set msg to "Vdata Name: <name>                                   */
  /*       CALL log_fmt_msg to report that the Vdata was not found in the   */
  /*            id table                                                    */
  /*         INPUTS:  Status, routine, msg                                  */
  /*         OUTPUTS: None                                                  */
  /*         RETURNS: None                                                  */
  /*    ENDIF                                                               */
  /*                                                                        */
  /**************************************************************************/

     else {
        if ((this_Vdata_index = get_index(Vdata_name)) == FAIL){
           sprintf(msg, "Vdata Name: %s", Vdata_name);
           log_fmt_msg(MODIS_W_VDATA_NOT_IN_TABLE, routine, msg);
        }

  /**************************************************************************/
  /*                                                                        */
  /*    CALL attached_Vdata_counter to decrement the number of attached     */
  /*         Vdatas                                                         */
  /*      INPUTS:  VU_DECREMENT                                             */
  /*      OUTPUTS: None                                                     */
  /*      RETURNS: number_attached_Vdatas                                   */
  /*                                                                        */
  /*    IF ( number_attached_Vdatas equals FAIL )                           */
  /*    THEN                                                                */
  /*       set Status to MODIS_W_NEGATIVE_VDATA_CTR                         */
  /*       CALL log_fmt_msg to report that the Vdata counter is negative    */
  /*         INPUTS:  Status, routine, Attempting to decrement              */
  /*                         attached_Vdata counter to a negative value     */
  /*         OUTPUTS: None                                                  */
  /*         RETURNS: None                                                  */
  /*    ENDIF                                                               */
  /*                                                                        */
  /*    {remove the item from the id table}                                 */
  /*    set global_VU_VDATA_ID[this_Vdata_index].id to VU_EMPTY             */
  /*    set global_VU_VDATA_ID[this_Vdata_index].name to "EMPTY_SLOT"       */
  /* ENDIF                                                                  */
  /*                                                                        */
  /**************************************************************************/

        number_attached_Vdatas = attached_Vdata_counter(VU_DECREMENT);
        if (number_attached_Vdatas == FAIL) {
           log_fmt_msg(MODIS_W_NEGATIVE_VDATA_CTR, routine, 
              "Attempting to decrement attached_Vdata counter to a negative value");
        }

        global_VU_VDATA_ID[this_Vdata_index].id = VU_EMPTY;
        strcpy(global_VU_VDATA_ID[this_Vdata_index].name, EN_EMPTY_SLOT);
     }

 } /* End of routine forget */
