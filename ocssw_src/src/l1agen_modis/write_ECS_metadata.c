#include "L1A_prototype.h"
#include "hdf.h"
#include "MD_metadata.h"
#include "mapi.h"
#include "PGS_SMF.h"
#include "PGS_MET.h"
#include "PGS_MODIS_35005.h"
#include "SDST_SetLGranId.h"
#include "PGS_MODIS_39604.h"

PGSt_SMF_status write_ECS_metadata (PGSt_MET_all_handles md_handles,
                                    MD_ECS_GRA_INV_MET_t *ecs_gra_inv_met)

/****************************************************************************
!C

!Description:  Function write_ECS_metadata assigns values of the ECS Granule 
               Inventory Metadata to metadata attributes after an MCF file 
               is initialized into memory.
                
!Input Parameters:
       PGSt_MET_all_handles  md_handles        metadata group in MCF
       MD_ECS_GRA_INV_MET_t  ecs_gra_inv_met   ECS Inventory Metadata structure

!Output Parameters:         None

Return Values:                
       MODIS_E_GETCONFIG_FAILED         If the universal reference number of
                                        engineering data list file was not found.
       MODIS_E_MET_SETATTRA_FAILED      If PGS_MET_SetAttr() failed.
       MODIS_E_NULL_POINTER             If either input argument is NULL
       MODIS_S_SUCCESS                  Otherwise

Externally Defined:
                ARCHIVED_METADATA               (mapi.h)
                global_input_pointer            (L1A_prototype.h)
                INVENTORY_METADATA              (mapi.h)
                MCORE_ACTUALLY_REDONE           (mapi.h)
                MCORE_ADDATTRIBUTENAME          (mapi.h)
                MCORE_DAYNIGHTFLAG              (mapi.h)
                MCORE_EQUATCROSSINGDATE         (mapi.h)
                MCORE_EQUATCROSSINGLONG         (mapi.h)
                MCORE_EQUATCROSSINGTIME         (mapi.h)
                MCORE_EXCLUS_GRING_FLG          (mapi.h)
                MCORE_GRING_POINT_LAT           (mapi.h)
                MCORE_GRING_POINT_LON           (mapi.h)
                MCORE_GRING_POINT_NUM           (mapi.h)
                MCORE_INPUT_POINTER             (mapi.h)
                MCORE_LOCALGRANULEID            (mapi.h)
                MCORE_LOCALVERSIONID            (mapi.h)
                MCORE_ORBIT_NUM                 (mapi.h)
                MCORE_PARAMETERVALUE            (mapi.h)
                MCORE_PGEVERSION                (mapi.h)
                MCORE_PRODUCTIONDATETIME        (mapi.h)
                MCORE_RANGE_BEG_DATE            (mapi.h)
                MCORE_RANGE_BEG_TIME            (mapi.h)
                MCORE_RANGE_ENDING_DATE         (mapi.h)
                MCORE_RANGE_ENDING_TIME         (mapi.h)
                MCORE_TO_BE_REDONE              (mapi.h)
                MD_ECS_GRA_INV_MET_t            (MD_metadata.h)
                MECS_PRODHISTORY                (MD_metadata.h)
                EASTBOUNDINGCOORDINATE          (MD_metadata.h)
                WESTBOUNDINGCOORDINATE          (MD_metadata.h)
                SOUTHBOUNDINGCOORDINATE         (MD_metadata.h)
                NORTHBOUNDINGCOORDINATE         (MD_metadata.h)
                EASTBOUNDVALUE                  (MD_metadata.h)
                WESTBOUNDVALUE                  (MD_metadata.h)
                SOUTHBOUNDVALUE                 (MD_metadata.h)
                NORTHBOUNDVALUE                 (MD_metadata.h)
                MODIS_E_GETCONFIG_FAILED        (PGS_MODIS_35005.h)
                MODIS_E_MET_SETATTRA_FAILED     (PGS_MODIS_35005.h)
                MODIS_E_NULL_POINTER            (PGS_MODIS_35005.h)
                MODIS_S_SUCCESS                 (PGS_MODIS_35005.h)
                PC_L1A_ENG_DATA_LIST_FILE       (PC_pcf_info.h)
                PC_L1A_PROCESSINGENVIRONMENT    (PC_pcf_info.h)
                PGS_SMF_MAX_MSGBUF_SIZE         (PGS_SMF.h)
                PGSt_SMF_status                 (PGS_SMF.h)
                PGS_TRUE                        (PGS_SMF.h)
                PGSd_PC_VALUE_LENGTH_MAX        (PGS_PC.h)

Called By:
               write_global_metadata

Routines Called:
               PGS_MET_SetAttr
               PGS_PC_GetUniversalRef
               PGS_SMF_TestSuccessLevel
               log_fmt_msg

!Revision History:
  $Log: write_ECS_metadata.c,v $
  Revision 5.1  2004/09/23 19:17:54  seaton
  Added Bounding Coordinates to the Archive Metadata set in the L1A product file.
  seaton@saicmodis.com

  Revision 4.1  2003/11/12 21:08:12  kuyper
  Corrected to always provide a terminating NULL pointer in input_pointers.

  Revision 4.0  2003/01/07 21:37:01  vlin
  Updated according to write_ECS_metadata.pdl v4.1
  vlin@saicmodis.com

               revision 1.0 1997/08/20  17:30:00
               Qi Huang/RDC    (qhuang@ltpmail.gsfc.nasa.gov)
               Original development

!Team-unique Header:

         This software is developed by the MODIS Science Data Support Team
         for the National Aeronautics and Space Administration,
         Goddard Space Flight Center, under contract NAS5-32373.

References and Credits:

Design Notes:

!END
******************************************************************************/

{
  char              *routine = "write_ECS_metadata";
  char              msg[PGS_SMF_MAX_MSGBUF_SIZE];
  char              *pchar, *input_ptr[MAX_INPUTS+1];
  PGSt_integer      version_num;
  PGSt_SMF_status   returnStatus = MODIS_S_SUCCESS;
  PGSt_SMF_status   PGS_status;
  double fval[1];

  if (md_handles == NULL || ecs_gra_inv_met == NULL) {
      log_fmt_msg(MODIS_E_NULL_POINTER, routine, "");
      return MODIS_E_NULL_POINTER;
  }

  version_num = 1;
  PGS_status = PGS_PC_GetUniversalRef(PC_L1A_ENG_DATA_LIST_FILE, &version_num, 
               global_input_pointer[0]);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      sprintf(msg, "unable to retrieve eng data list file name from pcf file");
      log_fmt_msg(MODIS_E_GETCONFIG_FAILED, routine, msg);
      return MODIS_E_GETCONFIG_FAILED;
  }

  input_ptr[0] = global_input_pointer[0];
  input_ptr[1] = global_input_pointer[1];
  input_ptr[2] = global_input_pointer[2];
  input_ptr[MAX_INPUTS] = (char*)NULL;

  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], MCORE_PGEVERSION,
                               &ecs_gra_inv_met->pgeversion);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s attribute could not be set", MCORE_PGEVERSION);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }
  
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
               MCORE_LOCALVERSIONID, &ecs_gra_inv_met->localversionid);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s attribute could not be set", MCORE_LOCALVERSIONID);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }

  pchar = ecs_gra_inv_met->productionhistory;
  PGS_status = PGS_MET_SetAttr(md_handles[ARCHIVED_METADATA], 
                               MECS_PRODHISTORY, &pchar);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s attribute could not be set", MECS_PRODHISTORY);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }

  pchar = ecs_gra_inv_met->rangebeginningdate;
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                               MCORE_RANGE_BEG_DATE, &pchar);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s attribute could not be set", MCORE_RANGE_BEG_DATE);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }
  
  pchar = ecs_gra_inv_met->rangebeginningtime;
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                               MCORE_RANGE_BEG_TIME, &pchar);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s attribute could not be set", MCORE_RANGE_BEG_TIME);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }
    
  pchar = ecs_gra_inv_met->rangeendingdate;
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                               MCORE_RANGE_ENDING_DATE, &pchar);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s attribute could not be set", MCORE_RANGE_ENDING_DATE);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }
    
  pchar = ecs_gra_inv_met->rangeendingtime;
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                               MCORE_RANGE_ENDING_TIME, &pchar);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s attribute could not be set", MCORE_RANGE_ENDING_TIME);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
    }
    
  pchar = ecs_gra_inv_met->day_night_flag;
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                               MCORE_DAYNIGHTFLAG, &pchar);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s attribute could not be set", MCORE_DAYNIGHTFLAG);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }
    
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                   MCORE_ORBIT_NUM".1", &ecs_gra_inv_met->orbit_num_1);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s.1 attribute could not be set", MCORE_ORBIT_NUM);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }
    
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                               MCORE_EQUATCROSSINGLONG".1", 
                       &ecs_gra_inv_met->equatorcrossinglongitude_1);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s.1 attribute could not be set", 
                  MCORE_EQUATCROSSINGLONG);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }
    
  pchar = ecs_gra_inv_met->equatorcrossingdate_1;
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                   MCORE_EQUATCROSSINGDATE".1", &pchar);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s.1 attribute could not be set", 
                   MCORE_EQUATCROSSINGDATE);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }
    
  pchar = ecs_gra_inv_met->equatorcrossingtime_1;
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                   MCORE_EQUATCROSSINGTIME".1", &pchar);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s.1 attribute could not be set", 
                   MCORE_EQUATCROSSINGTIME);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }
  
  pchar = ecs_gra_inv_met->exclusiongringflag_1;
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                   MCORE_EXCLUS_GRING_FLG, &pchar);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s.1 attribute could not be set", MCORE_EXCLUS_GRING_FLG);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }
  
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                   MCORE_GRING_POINT_LAT,
                   ecs_gra_inv_met->gringpointlatitude_1);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s.1 attribute could not be set", MCORE_GRING_POINT_LAT);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }
    
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                   MCORE_GRING_POINT_LON,
                   ecs_gra_inv_met->gringpointlongitude_1);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s.1 attribute could not be set", MCORE_GRING_POINT_LON);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }
    
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                   MCORE_GRING_POINT_NUM,
                   ecs_gra_inv_met->gringpointsequenceno_1);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s.1 attribute could not be set", MCORE_GRING_POINT_NUM);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }
    
  pchar = ecs_gra_inv_met->additionalattributename_1;
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                   MCORE_ADDATTRIBUTENAME".1", &pchar);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s.1 attribute could not be set", 
                   MCORE_ADDATTRIBUTENAME);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }
    
  pchar = ecs_gra_inv_met->parametervalue_1;
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                   MCORE_PARAMETERVALUE".1", &pchar);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s.1 attribute could not be set", MCORE_PARAMETERVALUE);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }
  
  pchar = ecs_gra_inv_met->additionalattributename_2;
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                   MCORE_ADDATTRIBUTENAME".2", &pchar);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s.2 attribute could not be set", 
                   MCORE_ADDATTRIBUTENAME);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }

  pchar = ecs_gra_inv_met->parametervalue_2;
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                   MCORE_PARAMETERVALUE".2", &pchar);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s.2 attribute could not be set", MCORE_PARAMETERVALUE);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }
  
  pchar = ecs_gra_inv_met->reprocessingactual;
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                   MCORE_ACTUALLY_REDONE, &pchar);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s attribute could not be set", MCORE_ACTUALLY_REDONE);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }

  pchar = ecs_gra_inv_met->reprocessingplanned;
  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                   MCORE_TO_BE_REDONE, &pchar);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s attribute could not be set", MCORE_TO_BE_REDONE);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }

  pchar = ecs_gra_inv_met->processingenvironment;
  PGS_status = PGS_MET_SetAttr(md_handles[ARCHIVED_METADATA],
                   PC_L1A_PROCESSINGENVIRONMENT, &pchar);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s attribute could not be set", 
                   PC_L1A_PROCESSINGENVIRONMENT);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  } 

  PGS_status = PGS_MET_SetAttr(md_handles[INVENTORY_METADATA], 
                   MCORE_INPUT_POINTER, &input_ptr);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s attribute could not be set", MCORE_INPUT_POINTER);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }

  PGS_status = SDST_SetLocalGranId(L1A, md_handles);
  if (PGS_status != MODIS_S_SDST_SUCCESS) {
      sprintf(msg,"The %s attribute could not be set", MCORE_LOCALGRANULEID);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }

  *fval = EASTBOUNDVALUE;
  PGS_status = PGS_MET_SetAttr(md_handles[ARCHIVED_METADATA],
                    EASTBOUNDINGCOORDNIATE, &fval);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s attribute could not be set", EASTBOUNDINGCOORDNIATE);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }

  *fval = WESTBOUNDVALUE;
  PGS_status = PGS_MET_SetAttr(md_handles[ARCHIVED_METADATA],
                    WESTBOUNDINGCOORDNIATE, &fval);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s attribute could not be set", WESTBOUNDINGCOORDNIATE);       
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }

  *fval = SOUTHBOUNDVALUE;
  PGS_status = PGS_MET_SetAttr(md_handles[ARCHIVED_METADATA],
                     SOUTHBOUNDINGCOORDNIATE, &fval);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s attribute could not be set", SOUTHBOUNDINGCOORDNIATE);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  } 

  *fval = NORTHBOUNDVALUE;
  PGS_status = PGS_MET_SetAttr(md_handles[ARCHIVED_METADATA],
                     NORTHBOUNDINGCOORDNIATE, &fval);
  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_E_MET_SETATTR_FAILED;
      sprintf(msg,"The %s attribute could not be set", NORTHBOUNDINGCOORDNIATE);
      log_fmt_msg(MODIS_E_MET_SETATTR_FAILED, routine, msg);
  }

  return returnStatus;
}
