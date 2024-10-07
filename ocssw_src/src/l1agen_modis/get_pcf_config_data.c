#include <stdio.h>
#include <string.h>
#include "L1A_prototype.h"
#include "PC_pcf_info.h"
#include "PGS_MODIS_35005.h"
#include "PGS_PC.h"
#include "PGS_SMF.h"
#include "PGS_TD.h"
#include "PGS_TYPES.h"

PGSt_SMF_status  get_pcf_config_data(PCF_CONFIG_t *pcf_config)

/***************************************************************************
!C

!Description:  This function will read the parameters from the PCF and convert
               them from strings to the appropriate PGSt types.

!Input Parameters: 
               None.

!Output Parameters:
               PCF_CONFIG_t  *pcf_config ** Structure to be filled in with  **
                                         ** PCF configuration parameters    **

Return Values: MODIS_S_SUCCESS                 (PGS_MODIS_35005.h)
               MODIS_E_GETCONFIG_FAILED        (PGS_MODIS_35005.h)
               MODIS_E_NULL_POINTER            (PGS_MODIS_35005.h)
               MODIS_E_UTC_TO_TAI_FAILED       (PGS_MODIS_35005.h)

Externally Defined:  
               INSTRUMENT_AQUA                  (PC_pcf_info.h)
               INSTRUMENT_TERRA                 (PC_pcf_info.h)
               MODIS_E_GETCONFIG_FAILED         (PGS_MODIS_35005.h)
               MODIS_E_NULL_POINTER             (PGS_MODIS_35005.h)
               MODIS_E_UTC_TO_TAI_FAILED        (PGS_MODIS_35005.h)
               MODIS_S_SUCCESS                  (PGS_MODIS_39501.h)
               PC_L1A_FIRST_GRAN_START_TIME     (PC_pcf_info.h)
               PC_L1A_GRAN_TIME_LENGTH          (PC_pcf_info.h)
               PC_L1A_LAST_GRAN_STOP_TIME       (PC_pcf_info.h)
               PC_L1A_PROCESSING_ENVIRONMENT    (PC_pcf_info.h)
               PC_L1A_REPROCESS_ACTUAL          (PC_pcf_info.h)
               PC_L1A_REPROCESS_PLANNED         (PC_pcf_info.h)
               PC_L1A_SCAN_RATE                 (PC_pcf_info.h)
               PC_L1A_PGE_VERSION               (PC_pcf_info.h)
               PC_L1A_VERSION                   (PC_pcf_info.h)
               PGS_TRUE                         (PGS_SMF.h)
               PGS_SMF_MASK_LEV_E               (PGS_SMF.h)
               PGSPC_W_NO_CONFIG_FOR_ID         (PGS_PC.h)
               PGSt_SMF_status                  (PGS_SMF.h)
               PGSd_PC_VALUE_LENGTH_MAX         (PGS_PC.h)
               PGSd_EOS_AM1                     (PGS_TD.h)
               PGSd_EOS_PM1                     (PGS_TD.h)
               SYS_NMLN                         (sys/utsname.h)

Called By:     initialize_level1a.c

Routines Called:
               PGS_TD_UTCtoTAI         
               PGS_PC_GetConfigData    
               PGS_SMF_TestSuccessLevel
               PGS_SMF_TestStatusLevel
               log_fmt_msg               

!Revision History:
 $Log: get_pcf_config_data.c,v $
 Revision 4.1  2002/10/03 16:14:30  kuyper
 Changed to work correctly, even if getenv() reuses the same buffer on every
   call.

 Revision 3.3  2002/10/03 14:27:07  vlin
 Updated to use getenv() instead of uname().  uname() is a prohibited function.

 Revision 3.2  2002/10/01 14:53:31  vlin
 Updated after code walkthrough, PDL v3.5
 vlin@saicmodis.com

 Revision 3.1  2002/09/19 17:52:41  vlin
 Updated according to get_pcf_config_data.PDL revision 3.4

               Revision 2.2  2001/01/04
               John Seaton  (seaton@ltpmail.gsfc.nasa.gov)
               Added call to get PGE version from the pcf file.

!Team-unique Header:

         This software is developed by the MODIS Science Data Support Team 
         for the National Aeronautics and Space Administration, Goddard 
         Space Flight Center, under contract NAS5-32373.

!END
****************************************************************************/

{
  char  *routine = "get_pcf_config_data";
  char  buffer[PGSd_PC_VALUE_LENGTH_MAX] = "";
  const char  *env_var[] = { "OSTYPE", "HOST", "REVISION", "MACHINE"};
#define ENVSTRS (sizeof env_var/sizeof env_var[0])
  char	*p, *pend;
  const char  *env_str;
  int	i;
  PGSt_SMF_status  returnStatus = MODIS_S_SUCCESS;  
  PGSt_SMF_status  tempStatus;   

  if (pcf_config == NULL) {
      log_fmt_msg(MODIS_E_NULL_POINTER, routine, "");
      return MODIS_E_NULL_POINTER;
  }

  tempStatus = PGS_PC_GetConfigData(PC_L1A_FIRST_GRAN_START_TIME, buffer);
 
  if (PGS_SMF_TestSuccessLevel(tempStatus)) {
     tempStatus = PGS_TD_UTCtoTAI(buffer, &pcf_config->first_gran_start_time);
     if (PGS_SMF_TestStatusLevel(tempStatus) >= PGS_SMF_MASK_LEV_E) {
        returnStatus = MODIS_E_UTC_TO_TAI_FAILED;
        strcat(buffer, "could not be converted");
        log_fmt_msg(MODIS_E_UTC_TO_TAI_FAILED, routine, buffer);
     }
  }
  else {
     returnStatus = MODIS_E_GETCONFIG_FAILED;
     log_fmt_msg(MODIS_E_GETCONFIG_FAILED, routine,
                 "The first granule start time could not be retrieved");
  }

  tempStatus = PGS_PC_GetConfigData(PC_L1A_LAST_GRAN_STOP_TIME, buffer);
 
  if (PGS_SMF_TestSuccessLevel(tempStatus)) {
     tempStatus = PGS_TD_UTCtoTAI(buffer, &pcf_config->last_gran_stop_time);
     if (PGS_SMF_TestStatusLevel(tempStatus) >= PGS_SMF_MASK_LEV_E) {
        returnStatus = MODIS_E_UTC_TO_TAI_FAILED;
        strcat(buffer, "could not be converted"); 
        log_fmt_msg(MODIS_E_UTC_TO_TAI_FAILED, routine, buffer);
     }
  }
  else {
     returnStatus = MODIS_E_GETCONFIG_FAILED;
     log_fmt_msg(MODIS_E_GETCONFIG_FAILED, routine,
                 "The last granule stop time could not be retrieved");
  }

  tempStatus = PGS_PC_GetConfigData(PC_L1A_GRAN_TIME_LENGTH, buffer);
  
  if (PGS_SMF_TestSuccessLevel(tempStatus)) 
     pcf_config->gran_time_length = (PGSt_double) atof(buffer);
  else {
     returnStatus = MODIS_E_GETCONFIG_FAILED;
     log_fmt_msg(MODIS_E_GETCONFIG_FAILED, routine,
                 "The granule time length could not be retrieved");
  }

  tempStatus = PGS_PC_GetConfigData(PC_L1A_SCAN_RATE, buffer);
 
  if (PGS_SMF_TestSuccessLevel(tempStatus)) 
     pcf_config->scan_rate = (PGSt_double) atof(buffer);
  else {
     returnStatus = MODIS_E_GETCONFIG_FAILED;
     log_fmt_msg(MODIS_E_GETCONFIG_FAILED, routine,
                 "The scan rate could not be retrieved");
  }

  tempStatus = PGS_PC_GetConfigData(PC_L1A_VERSION, pcf_config->localversionid);
  
  if (!PGS_SMF_TestSuccessLevel(tempStatus)) {
     returnStatus = MODIS_E_GETCONFIG_FAILED;
     log_fmt_msg(MODIS_E_GETCONFIG_FAILED, routine,
       "The L1A product version could not be retrieved");
  }

  tempStatus = PGS_PC_GetConfigData(PC_L1A_LUT_REVISION, pcf_config->lutrevision);
  
  if (!PGS_SMF_TestSuccessLevel(tempStatus)) {
     returnStatus = MODIS_E_GETCONFIG_FAILED;
     log_fmt_msg(MODIS_E_GETCONFIG_FAILED, routine,
       "The MOD_PR01 LUT RCS Revision number could not be retrieved");
  }

  tempStatus = PGS_PC_GetConfigData(PC_L1A_PGE_VERSION, pcf_config->pgeversion);
  
  if (!PGS_SMF_TestSuccessLevel(tempStatus)) {
     returnStatus = MODIS_E_GETCONFIG_FAILED;
     log_fmt_msg(MODIS_E_GETCONFIG_FAILED, routine,
       "The PGE version could not be retrieved");
  }

  tempStatus = PGS_PC_GetConfigData(PC_INSTRUMENT, buffer);
  
  if (!PGS_SMF_TestSuccessLevel(tempStatus)) {
     pcf_config->instrument = 0;
     returnStatus = MODIS_E_GETCONFIG_FAILED;
     log_fmt_msg(MODIS_E_GETCONFIG_FAILED, routine,
       "The satellite instrument string could not be retrieved");
  }
  else if (strcmp(buffer,INSTRUMENT_AQUA) == 0)
     pcf_config->instrument = PGSd_EOS_PM1;
  else if (strcmp(buffer,INSTRUMENT_TERRA) == 0)
     pcf_config->instrument = PGSd_EOS_AM1;
  else 
     pcf_config->instrument = 0;

  tempStatus = PGS_PC_GetConfigData(PC_L1A_PROCESSING_ENVIRONMENT, 
               pcf_config->processingenvironment);
  
  if (tempStatus == PGSPC_W_NO_CONFIG_FOR_ID){
     /* Would like to call uname() here, but that's a prohibited function.
      * Instead, we'll rely upon certain environment variables being set up.
      * It's not really very important what this string contains, except in the
      * MODAPS environment, where the tempStatus should be PGS_S_SUCCESS.
      */
     p = pcf_config->processingenvironment;
     pend = p+sizeof(pcf_config->processingenvironment);

     for(i=0; i<ENVSTRS; i++)
     {
	 env_str = getenv(env_var[i]);
	 if(env_str)
	 {	/* The string does have a value */
	    if(p!=pcf_config->processingenvironment)
	    {	/* We've already found one previous string */
		*p++=' ';	/* So we need a space between them. */
		*p='\0';
	    }

	    strncat(p, env_str, (int)(pend-p) );
	    p += strlen(p);
	 }
      }
  }
  else if (tempStatus != PGS_S_SUCCESS) {
     returnStatus = MODIS_E_GETCONFIG_FAILED;
     log_fmt_msg(MODIS_E_GETCONFIG_FAILED, routine,
        "The ProcessingEnvironment string could not be retrieved");
  }

  tempStatus = PGS_PC_GetConfigData(PC_L1A_REPROCESS_ACTUAL, 
               pcf_config->reprocessingactual);
  
  if (!PGS_SMF_TestSuccessLevel(tempStatus)) {
     returnStatus = MODIS_E_GETCONFIG_FAILED;
     log_fmt_msg(MODIS_E_GETCONFIG_FAILED, routine,
       "The ReprocessingActual string could not be retrieved");
  }

  tempStatus = PGS_PC_GetConfigData(PC_L1A_REPROCESS_PLANNED, 
               pcf_config->reprocessingplanned);
  
  if (!PGS_SMF_TestSuccessLevel(tempStatus)) {
     returnStatus = MODIS_E_GETCONFIG_FAILED;
     log_fmt_msg(MODIS_E_GETCONFIG_FAILED, routine,
       "The ReprocessingPlanned string could not be retrieved");
  }

  return returnStatus;

 } /* End of routine get_pcf_config_data */
