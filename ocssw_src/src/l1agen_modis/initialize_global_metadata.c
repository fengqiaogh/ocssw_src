#include <stdlib.h>
#include "L1A_prototype.h"
#include "hdf.h"
#include "PGS_MODIS_35005.h"
#include "PGS_TYPES.h"
#include "PGS_SMF.h"
#include "MD_metadata.h"
#include "PC_pcf_info.h"
#include "version.h"
#include "PGS_PC.h"
#include "smfio.h"

PGSt_SMF_status initialize_global_metadata(
                           PGSt_double gran_start_time,
                           PGSt_double gran_end_time, 
                           int nscans,
                           PCF_CONFIG_t *pcf_config,
                           MD_ECS_GRA_INV_MET_t *ecs_gra_inv_met,
                           MD_L1A_SPECIFIC_MET_t *l1a_specific_met)

/******************************************************************************
!C

!Description:  Function initialize_global_metadata initializes the ECS Granule
               Inventory Metadata (MD_ECS_GRA_INV_MET_t) and the MODIS Level1A
               Specific Granule Metadata (MD_L1A_SPECIFIFC_MET_t).
               
!Input Parameters:
       PGSt_double    gran_start_time   Start time of the L1A granule
       PGSt_double    gran_end_time     Stop time of the L1A granule.
       int            nscans            Number of scans             
       PCF_CONFIG_t   pcf_config        PCF configuration parameters

!Output Parameters:
       MD_ECS_GRA_INV_MET_t  ecs_gra_inv_met  ECS granule inventory metadata. 
       MD_L1A_SPECIFIC_MET_t l1a_specific_met Level 1A specific granule metadata

Return Values:                
       MODIS_E_NULL_POINTER             If any pointer argument is NULL
       MODIS_S_SUCCESS                  If successful
       PGS_E_TOOLKIT                    From subroutines
       PGS_S_SUCCESS                    From subroutines
       PGSTD_E_NO_LEAP_SECS             From subroutines

Externally Defined:
       MD_ECS_GRA_INV_MET_t             (MD_metadata.h)
       MD_GRANULENUMBER                 (MD_metadata.h)
       MD_L1A_SPECIFIC_MET_t            (MD_metadata.h)
       MD_MAX_MISSING_PKTS_IN_SCAN      (MD_metadata.h)
       MD_MIDNIGHT                      (MD_metadata.h)
       MD_NA                            (MD_metadata.h)
       MD_OTHER_STRING                  (MD_metadata.h)
       MD_PROCESSVERSION                (MD_metadata.h)
       MODIS_E_NULL_POINTER             (PGS_MODIS_35005.h)
       MODIS_E_TAI_TO_UTC_FAILED        (PGS_MODIS_35005.h)
       MODIS_E_UTC_TO_TAI_FAILED        (PGS_MODIS_35005.h)
       MODIS_S_SUCCESS                  (PGS_MODIS_35005.h)
       PCF_CONFIG_t                     (PC_pcf_info.h)
       PGS_S_SUCCESS                    (PGS_SMF.h)
       PGS_SMF_MASK_LEV_E               (PGS_SMF.h)
       PGS_E_TOOLKIT                    (PGS_SMF.h)
       PGSt_SMF_status                  (PGS_SMF.h)
       PGSt_double                      (PGS_types.h)
       PGSTD_E_NO_LEAP_SECS             (PGS_TD_3.h)
       PROCESSVERSION                   (version.h)
       TIMECODEASIZE                    (smfio.h)

Called By:
       process_a_granule

Routines Called:
       PGS_TD_TAItoUTC
       PGS_TD_UTCtoTAI
       PGS_SMF_TestStatusLevel
       log_fmt_msg

!Revision History:
       $Log: initialize_global_metadata.c,v $
       Revision 6.1  2010/08/24 16:08:49  kuyper
       Corrected "%lf" sprintf() format specifiers to "%f".

       Revision 5.1  2006/01/02 19:38:32  kuyper
       Changed calculation of granule number to work properly even for granules
        covering leap seconds.

       Revision 4.2  2003/01/08 19:52:43  vlin
       Initialization of field name "operation_mode" removed.

       Revision 4.1  2002/12/04 20:11:54  vlin
       Use TIMECODEASIZE for buffer size

       Revision 4.0  2002/10/30 14:58:38  vlin
       Updated to match  initialize_global_metadata.pdl  v4.2
       vlin@saicmodis.com

               revision 1.0 1997/08/20  17:30:00
               Qi Huang/RDC    (qhuang@ltpmail.gsfc.nasa.gov)
               Original development

!Team-unique Header:

       This software is developed by the MODIS Science Data Support Team 
       for the National Aeronautics and Space Administration, Goddard 
       Space Flight Center, under contract NAS5-32373.

References and Credits:
       None

Design Notes:
       None

!END
****************************************************************************/

{
  char             *routine = "initialize_global_metadata";
  PGSt_SMF_status  returnStatus = MODIS_S_SUCCESS;
  PGSt_SMF_status  tempStatus, level;
  char             temp_string[TIMECODEASIZE], midnight_UTC[TIMECODEASIZE];
  char             msg[128];
  PGSt_double      midnight_TAI;
  int              i, temp_granule_number;

  if (pcf_config==NULL || ecs_gra_inv_met==NULL || l1a_specific_met==NULL) {
      log_fmt_msg(MODIS_E_NULL_POINTER, routine, "");
      return MODIS_E_NULL_POINTER; 
  }

  tempStatus = PGS_TD_TAItoUTC(gran_start_time,temp_string);

  level = PGS_SMF_TestStatusLevel(tempStatus);

  if (level >= PGS_SMF_MASK_LEV_E) {
      sprintf(msg,"The granule start time %f could not be converted from TAI "
              "to UTC", gran_start_time);
      returnStatus = tempStatus;
      log_fmt_msg(MODIS_E_TAI_TO_UTC_FAILED, routine, msg);
  }
  else {
      strcpy(ecs_gra_inv_met->rangebeginningdate,strtok(temp_string,"T"));
      strcpy(ecs_gra_inv_met->rangebeginningtime,strtok(NULL,"Z"));
  }

  ecs_gra_inv_met->additionalattributename_1 = MD_GRANULENUMBER;
  ecs_gra_inv_met->additionalattributename_2 = MD_PROCESSVERSION;
  ecs_gra_inv_met->parametervalue_2 = PROCESSVERSION;
  tempStatus = PGS_TD_TAItoUTC(gran_end_time,temp_string);
  level = PGS_SMF_TestStatusLevel(tempStatus);

  if (level >= PGS_SMF_MASK_LEV_E) {
      sprintf(msg, "The granule end time %f could not be converted from "
              "TAI to UTC", gran_end_time);
      returnStatus = tempStatus;
      log_fmt_msg(MODIS_E_TAI_TO_UTC_FAILED, routine, msg);
      strcpy(ecs_gra_inv_met->parametervalue_1, "0");
  }
  else {
      strcpy(ecs_gra_inv_met->rangeendingdate,strtok(temp_string,"T"));
      strcpy(ecs_gra_inv_met->rangeendingtime,strtok(NULL,"Z"));
 
      sprintf(midnight_UTC,ecs_gra_inv_met->rangeendingdate, MD_MIDNIGHT);
      tempStatus = PGS_TD_UTCtoTAI(midnight_UTC,&midnight_TAI);
      level = PGS_SMF_TestStatusLevel(tempStatus);
      if (level >= PGS_SMF_MASK_LEV_E) {
         sprintf(msg, "Midnight %s could not be converted from UTC to TAI "
                 "therefore no granule number can be computed", midnight_UTC);
         returnStatus = tempStatus;
         log_fmt_msg(MODIS_E_UTC_TO_TAI_FAILED, routine, msg);
         strcpy(ecs_gra_inv_met->parametervalue_1, "0");
      }
      else {
	 temp_granule_number = (int)((gran_end_time-midnight_TAI)/300) + 1;
         sprintf(ecs_gra_inv_met->parametervalue_1, "%d", temp_granule_number);
      }
  } 

  strcpy(ecs_gra_inv_met->day_night_flag, MD_NA);
  ecs_gra_inv_met->orbit_num_1 = -1;
  ecs_gra_inv_met->equatorcrossinglongitude_1 = (double)-1.0;
  ecs_gra_inv_met->equatorcrossingdate_1 = ecs_gra_inv_met->rangebeginningdate;
  ecs_gra_inv_met->equatorcrossingtime_1 = ecs_gra_inv_met->rangebeginningtime;
  strcpy(ecs_gra_inv_met->exclusiongringflag_1, "");

  for (i=0; i<4; i++) {
      ecs_gra_inv_met->gringpointlatitude_1[i] = (double)-1.0;
      ecs_gra_inv_met->gringpointlongitude_1[i] = (double)-1.0;
      ecs_gra_inv_met->gringpointsequenceno_1[i] = -1;
  }

  ecs_gra_inv_met->processingenvironment = pcf_config->processingenvironment;
  ecs_gra_inv_met->localversionid = pcf_config->localversionid;
  strcpy(ecs_gra_inv_met->productionhistory, "PGE01:");
  strcat(ecs_gra_inv_met->productionhistory, pcf_config->pgeversion);
  ecs_gra_inv_met->pgeversion = pcf_config->pgeversion;
  ecs_gra_inv_met->reprocessingactual = pcf_config->reprocessingactual;
  ecs_gra_inv_met->reprocessingplanned = pcf_config->reprocessingplanned;

/* Set Level 1 A specific granule metadata */
  l1a_specific_met->num_scans = 0;
  l1a_specific_met->num_day_scans = 0;
  l1a_specific_met->num_night_scans = 0;
  l1a_specific_met->max_total_frames = 0;
  l1a_specific_met->max_earth_frames = 0;
  l1a_specific_met->max_sd_frames = 0;
  l1a_specific_met->max_srca_frames = 0;
  l1a_specific_met->max_bb_frames = 0;
  l1a_specific_met->max_sv_frames = 0;
  memset(l1a_specific_met->scan_types_product, '\0', 
         sizeof(l1a_specific_met->scan_types_product));
  strcpy(l1a_specific_met->scan_types_product,MD_OTHER_STRING);
  l1a_specific_met->incomplete_scans = 0;
  l1a_specific_met->missing_packets = nscans * MD_MAX_MISSING_PKTS_IN_SCAN;
  l1a_specific_met->packets_bad_crc = 0;
  l1a_specific_met->discarded_packets = 0;

  return returnStatus;
}
