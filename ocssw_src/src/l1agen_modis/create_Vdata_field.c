#include "L1A_prototype.h"
#include "PGS_SMF.h"
#include "PGS_MODIS_35005.h"
#include "hdf.h"
#include "hdfi.h"


PGSt_SMF_status  create_Vdata_field (char    *Vdata_name, 
                                     int32    Vdata_id, 
                                     char    *field_name,
                                     char    *data_type,
                                     int32    order)
/*
!C************************************************************************

!Description:  This routine creates the specified Vdata field inside the 
               specified Vdata within the currently open hdf file.

!Input Parameters:
               char    *Vdata_name      ** the name of the Vdata in which  **
                                        ** to create the new Vdata field.  **
 
               int32   Vdata_id         ** the Vdata's id                  **

               char    *field_name      ** the name of the field to create **

               char    *data_type       ** the field's data type           **

               int32   order            ** the number of data type size    **
                                        ** entries per record              **
!Output Parameters:
               None

Return Values:
               MODIS_S_SUCCESS                 (PGS_MODIS_35005.t)
               MODIS_E_VSFDEFINE               (PGS_MODIS_35005.t)
               MODIS_E_DATATYPE_FAIL           (PGS_MODIS_35005.t)
               FAIL                            (hdf.h)

Externally Defined:
               PGSt_SMF_status                 (PGS_SMF.h)
               int32                           (hdfi.h)   

Called By:
               create_Vdata

Routines Called:
               L1A_datatype_to_DFNT
               VSfdefine
               log_fmt_msg
               
!Revision History:
               Revision 2.0  1997/10/01  11:25 EDT
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
	       HDF portions developed at the National Center for 
               Supercomputing Applications at the University of Illinois 
               at Urbana-Champaign.               

!Design Notes:
               None

!END************************************************************************
*/

 {
  /**************************************************************************/
  /*                                                                        */
  /*          Declare the local variables and initialize them.              */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /* Set routine to "create_Vdata_field"                                    */
  /* declare number_type to be a variable of type int32                     */
  /* Set returnStatus to MODIS_S_SUCCESS                                    */
  /*                                                                        */
  /**************************************************************************/

  char  *routine = "create_Vdata_field";
  char  msg[300];

  PGSt_SMF_status  Status;        /* SMF-style message returned by function */
  PGSt_SMF_status  returnStatus;  /* SMF-style message returned by function */

  int32  number_type;

  returnStatus = MODIS_S_SUCCESS;   

  /**************************************************************************/
  /*                                                                        */
  /* CALL L1A_datatype_to_DFNT to assign the data type                      */
  /*   INPUTS:  data_type                                                   */
  /*   OUTPUTS: None                                                        */
  /*   RETURNS: number_type                                                 */
  /*                                                                        */
  /* IF number_type is not equal to MFAIL                                   */
  /* THEN                                                                   */
  /*    CALL VSfdefine to define the new field                              */
  /*      INPUTS:  Vdata_id, field_name, number_type, order                 */
  /*      OUTPUTS: None                                                     */
  /*      RESULT:  Status                                                   */
  /*                                                                        */
  /*    IF (Status equals FAIL)                                             */
  /*    THEN                                                                */
  /*       set Status to MODIS_E_VSFDEFINE                                  */
  /*       set msg to "unable to define the new field"                      */
  /*       CALL log_fmt_msg to report that the new field could not be       */
  /*            defined                                                     */
  /*         INPUTS:  Status, routine msg                                   */
  /*         OUTPUTS: None                                                  */
  /*         RETURNS: None                                                  */
  /*                                                                        */
  /*       set returnStatus to FAIL                                         */
  /*    ENDIF                                                               */
  /*                                                                        */
  /**************************************************************************/

     number_type = L1A_datatype_to_DFNT(data_type);

     if (number_type != MFAIL) {
        Status = VSfdefine(Vdata_id, field_name, number_type, order);

        if (Status == FAIL) {
           sprintf(msg, "unable to define the new field, %s for the Vdata %s",
                   field_name, Vdata_name);
           log_fmt_msg(MODIS_E_VSFDEFINE, routine, msg);
           returnStatus = FAIL;
        }
     }


  /**************************************************************************/
  /*                                                                        */
  /* ELSE                                                                   */
  /*    set Status to MODIS_E_DATATYPE_FAIL                                 */
  /*    set msg to "unable to associate the datatype for Vdata_name"        */
  /*    CALL log_fmt_msg to report that the datatype could not be associated*/
  /*      INPUTS:  Status, routine msg                                      */
  /*      OUTPUTS: None                                                     */
  /*      RETURNS: None                                                     */
  /*                                                                        */
  /*    set returnStatus to FAIL                                            */
  /* ENDIF                                                                  */
  /*                                                                        */
  /* return returnStatus                                                    */
  /*                                                                        */
  /**************************************************************************/

     else {
        sprintf(msg, "unable to associate the datatype for Vdata Name: %s",
                Vdata_name);
        log_fmt_msg(MODIS_E_DATATYPE_FAIL, routine, msg);
        returnStatus = FAIL;
     }

     return returnStatus;

 } /* End of routine create_Vdata_field */
