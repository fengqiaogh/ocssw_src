#include "L1A_prototype.h"
#include "hdf.h"
#include "MD_metadata.h"
#include "mapi.h"
#include "PGS_MODIS_35005.h"


PGSt_SMF_status write_global_metadata (MODFILE                *mfile,
                                       PGSt_MET_all_handles   md_handles,
                                       MD_ECS_GRA_INV_MET_t   *ecs_gra_inv_met,
                                       MD_L1A_SPECIFIC_MET_t  *l1a_specific_met)
/*
!C***************************************************************************

!Description:   Function write_global_metadata assigns values of the ECS 
                Granule Inventory Metadata (MD_ECS_INV_MET_t)to metadata 
                attributes in memory and stores L1A Specific Metadata 
                (MD_L1A_SPECIFIC_MET_t) in the L1A file.
                 
!Input Parameters:
                MODIFILE               *mfile             ** MODFILE structure      **
                PGSt_MET_all_handles   md_handles         ** metadata group in MCF  **
                MD_ECS_INV_MET_t       *ecs_gra_inv_met   ** ECS Granule Metadata   **
                MD_L1A_SPECIFIC_MET_t  *l1a_specific_met  ** L1A Specific Metadata  **

!Output Parameters:
                None

Return Values:                
                MODIS_S_SUCCESS                 (PGS_MODIS_35005.h)
                MODIS_E_WRITE_GLOBAL_METADATA   (PGS_MODIS_35005.h)

Externally Defined:
                MODFILE                         (mapi.h)
                PGSt_MET_all_handles            (PGS_MET.h)
                MD_ECS_INV_MET_t                (MD_metadata.h)
                MD_L1A_SPECIFIC_MET_t           (MD_metadata.h)

Called By:
                process_a_granule

Routines Called:
                write_ECS_metadata
                write_specific_granule_metadata
                log_fmt_msg

!Revision History:

                revision 1.0 1997/08/21  17:30:00
                Qi Huang/RDC    (qhuang@ltpmail.gsfc.nasa.gov)
                Original development

!Team-unique Header:
                This software is developed by the MODIS Science Data Support 
                Team (SDST) for the National Aeronautics and Space Administration 
                (NASA), Goddard Space Flight Center (GSFC), under contract 
                NAS5-32373.

!References and Credits:
                None

!Design Notes:
                None


!END***************************************************************************
*/
{
  /****************************************************************************/
  /*                                                                          */
  /*              Define and Initialize Local Variables                       */
  /*                                                                          */
  /****************************************************************************/

  char              *routine = "write_global_metadata";
  PGSt_SMF_status   returnStatus;
  PGSt_SMF_status   L1A_status;


  /****************************************************************************/
  /*                                                                          */
  /*  Set routine to "write_global_metadata"  (done during declaration)       */
  /*                                                                          */
  /*  Set returnStatus to MODIS_S_SUCCESS                                     */
  /*                                                                          */
  /****************************************************************************/

  returnStatus = MODIS_S_SUCCESS;


  /****************************************************************************/
  /*                                                                          */
  /*  CALL write_ECS_metadata to set the ECS Granule Metadata attributes      */
  /*    INPUT:  md_handles, ecs_gra_inv_met                                   */
  /*    OUTPUT: None                                                          */
  /*    RETURN: L1A_status                                                    */
  /*                                                                          */
  /*  IF L1A_status is not equal to MODIS_S_SUCCESS                           */
  /*  THEN                                                                    */
  /*    Set returnStatus to MODIS_E_WRITE_GLOBAL_METADATA                     */
  /*    Set msg to "The ECS Granule Metadata attributes could not be set"     */
  /*    CALL log_fmt_msg to report that not all of the ECS Granule Metadata   */
  /*      attributes could be set                                             */
  /*      INPUT:  L1A_status, routine, msg                                    */
  /*      OUTPUT: None                                                        */
  /*      RETURN: None                                                        */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/
  
  L1A_status = write_ECS_metadata(md_handles,ecs_gra_inv_met);
  if ( L1A_status != MODIS_S_SUCCESS)
    {
      returnStatus = MODIS_E_WRITE_GLOBAL_METADATA;
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, 
                  "The ECS Granule Metadata attributes could not be set");
    }
    
    
  /****************************************************************************/
  /*                                                                          */
  /*  CALL write_specific_granule_metadata to store the L1A Specific metadata */
  /*    in the L1A granule                                                    */
  /*    INPUT:  mfile, l1a_specific_met                                       */
  /*    OUTPUT: None                                                          */
  /*    RETURN: L1A_status                                                    */
  /*                                                                          */
  /*  IF L1A_status is not equal to MODIS_S_SUCCESS                           */
  /*  THEN                                                                    */
  /*    Set returnStatus to MODIS_E_WRITE_GLOBAL_METADATA                     */
  /*    Set msg to "The L1A Specific Metadata could not be written to the     */
  /*       L1A granule."                                                      */
  /*    CALL log_fmt_msg to report that not all of the L1A Specific Metadata  */
  /*       could be written to the L1A granule                                */
  /*      INPUT:  L1A_status, routine, msg                                    */
  /*      OUTPUT: None                                                        */
  /*      RETURN: None                                                        */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  L1A_status = write_specific_granule_metadata(mfile,l1a_specific_met);
  if ( L1A_status != MODIS_S_SUCCESS)
    {
      returnStatus = MODIS_E_WRITE_GLOBAL_METADATA;
      log_fmt_msg(MODIS_E_GATTRIB_FAILED, routine, 
            "The L1A Specific Metadata could not be written to the L1A granule");
    }


  /****************************************************************************/
  /*                                                                          */
  /* RETURN returnStatus                                                      */
  /*                                                                          */
  /****************************************************************************/
  return(returnStatus);
}
