#include "L1A_prototype.h"
#include "PGS_SMF.h"
#include "PGS_MODIS_35005.h"
#include "hdf.h"
#include "hdfi.h"


PGSt_SMF_status  write_Vdata ( char           *Vdata_name, 
                               unsigned char  *data, 
                               int32           num_records )

/*
!C************************************************************************

!Description:  This function writes to the specified Vdata in the currently 
               open hdf file, the Vdata values specified in data.

!Input Parameters:
               char          *Vdata_name    ** The name of the Vdata to write **
               unsigned char *data          ** The data to write (stored in a **
                                            **   byte array)                  **
               int32         num_records    ** The number of records to write **

!Output Parameters:
               None

Return Values: 
               MODIS_S_SUCCESS              (PGS_MODIS_35005.h)
               MODIS_E_RECALL_ID            (PGS_MODIS_35005.h)
               MODIS_E_WRITE_VDATA          (PGS_MODIS_35005.h)
               FAIL                         (hdf.h)

Externally Defined:  
               PGSt_SMF_status              (PGS_SMF.h)
               int32                        (hdfi.h)
               FULL_INTERLACE               (hdf.h)

Called By:
               write_eng_data
               write_failed_packets

Routines Called:
               recall_id
               VSwrite
               log_fmt_msg

!Revision History:
               Revision 2.0  1997/10/01  13:35 EDT
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
  /* Set routine to "write_Vdata"                                           */
  /* Set returnStatus to MODIS_S_SUCCESS                                    */
  /* declare num_rcds_written and Vdata_id to be variables of type int32    */
  /*                                                                        */
  /**************************************************************************/

  PGSt_SMF_status  returnStatus;  /* SMF-style message returned by function */

  char  *routine = "write_Vdata";

  char  msg[300];

  int32  num_rcds_written;

  int32  Vdata_id;


  returnStatus = MODIS_S_SUCCESS;



  /**************************************************************************/
  /*                                                                        */
  /* CALL recall_id to get the Vdata id for the Vdata                       */
  /*   INPUTS:  Vdata_name                                                  */
  /*   OUTPUTS: None                                                        */
  /*   RETURNS: Vdata_id                                                    */
  /*                                                                        */
  /* IF Vdata_id is not equal to FAIL                                       */
  /* THEN                                                                   */
  /*  {FULL_INTERLACE writes record by record (i.e all fields in a record}  */
  /*  {are written before moving on to the next record)}                    */
  /*    CALL VSwrite to write the data to file                              */
  /*      INPUTS:  Vdata_id, data, num_records, FULL_INTERLACE              */
  /*      OUTPUTS: None                                                     */
  /*      RETURNS: returnStatus                                             */
  /*                                                                        */
  /*    set num_rcds_written to returnStatus                                */
  /*                                                                        */
  /**************************************************************************/

     Vdata_id = recall_id(Vdata_name);

     if (Vdata_id != FAIL) {
        returnStatus = VSwrite(Vdata_id, data, num_records, FULL_INTERLACE);

        num_rcds_written = returnStatus;


  /**************************************************************************/
  /*                                                                        */
  /*    IF ( returnStatus equals FAIL )                                     */
  /*    THEN                                                                */
  /*       CALL log_fmt_msg to report that the Vdata could not be written   */
  /*            to the L1A granule                                          */
  /*         INPUTS:  Status, routine, msg                                  */
  /*         OUTPUTS: None                                                  */
  /*         RETURN:  None                                                  */
  /*    ELSE                                                                */
  /*       IF ( num_rcds_written is NOT equal to num_records )              */
  /*       THEN                                                             */
  /*          CALL log_fmt_msg to report that not all of the records        */
  /*               requested were written to the Vdata in the L1A granule   */
  /*            INPUTS:  Status, routine, msg                               */
  /*            OUTPUTS: None                                               */
  /*            RETURN:  None                                               */
  /*                                                                        */
  /*          set returnStatus to FAIL                                      */
  /*       ENDIF                                                            */
  /*    ENDIF                                                               */
  /*                                                                        */
  /**************************************************************************/

        if (returnStatus == FAIL) {
           sprintf(msg, "Vdata Name: %s", Vdata_name);
           log_fmt_msg (MODIS_E_WRITE_VDATA, routine, msg);
        }
        else if (num_rcds_written != num_records) {
           sprintf(msg, "VSwrite wrote fewer records than requested: "
	      "Num Written %ld  Num Requested %ld",
              (long)num_rcds_written, (long)num_records);
           log_fmt_msg(MODIS_E_WRITE_VDATA, routine, msg);

           returnStatus = FAIL;
        }
     }


  /**************************************************************************/
  /*                                                                        */
  /* ELSE                                                                   */
  /*    CALL log_fmt_msg to report that the Vdata id could not be retrieved */
  /*      INPUTS:  Status, routine, msg                                     */
  /*      OUTPUTS: None                                                     */
  /*      RETURN:  None                                                     */
  /*                                                                        */
  /*    set returnStatus to FAIL                                            */
  /* ENDIF                                                                  */
  /*                                                                        */
  /* return returnStatus                                                    */
  /*                                                                        */
  /**************************************************************************/

     else {
        sprintf(msg, "Vdata Name: %s", Vdata_name);
        log_fmt_msg(MODIS_E_RECALL_ID, routine, msg);
       
        returnStatus = FAIL;
     }

     return returnStatus;

 }  /* End of routine write_Vdata */
