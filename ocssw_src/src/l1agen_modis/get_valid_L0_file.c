#include "L1A_prototype.h"
#include "PGS_TYPES.h"
#include "PGS_MODIS_35005.h"
#include "PGS_IO_L0.h"
#include "PGS_SMF.h"
#include "PC_pcf_info.h"


PGSt_SMF_status  get_valid_L0_file (PGSt_tag                   spacecraft_tag,
                                    PGSt_IO_L0_VirtualDataSet  *L0_file,
                                    PGSt_double                *start_time, 
                                    PGSt_double                *stop_time)
/*                       
!C*****************************************************************************

!Description:  Function get_valid_L0_file opens the L0 file and validates the
               L0 header.

!Input Parameters: 
               PGSt_tag                spacecraft_tag    **  Spacecraft tag  **
               
!Output Parameters: 
               PGSt_double             start_time        **  Start time of 
                                                             L0 file         **
               PGSt_double             stop_time         **  Stop time of
                                                             L0 file         **
               PGEt_IO_L0_virtualDataSet  L0_file        opened L0 file structure

Return Values: 
               MODIS_S_SUCCESS                          (PGS_MODIS_35005.h)
               MODIS_E_GET_VALID_L0_FILE                (PGS_MODIS_35005.h)

Externally Defined:
               PGSt_PC_Logical                          (PGS_TYPES.h)
               PGSt_tag                                 (PGS_TYPES.h)
               PGSt_IO_L0_VirtualDataSet                (PGS_IO_L0.h)
               PGSt_double                              (PGS_TYPES.h)
               PC_CURRENT_L0_PCF_ID                     (PC_pcf_info.h)
               global_last_gran_stop_time               (level1a)
               global_first_gran_start_time             (level1a)
               PC_L1A_SCAN_RATE                         (PC_pcf_info.h)

Called By:
               initialize_level1a
               read_a_packet
               set_start_position

Routines Called:
               PGS_IO_L0_Open                          
               PGS_SMF_TestSuccessLevel
               validate_L0_header
               log_fmt_msg
               PGS_PC_GetConfigData
               PGS_PC_GetUniversalRef

!Revision History:    
  $Log: get_valid_L0_file.c,v $
  Revision 6.1  2010/08/24 13:54:43  kuyper
  Corrected inappropriate use of NULL in an arithmetic context.

  Revision 4.2  2003/11/12 21:05:37  kuyper
  Removed declarations that duplicate L1A_prototype.h entries.

  Revision 4.1  2002/10/21 19:01:59  vlin
  Updated buffer format for variable "global_L0_logical".

               revision 1.1 1999/11/18  10:32:00
               John Seaton/GSC  (seaton@ltpmail.gsfc.nasa.gov)
               Added scan rate on to last_gran_stop_time to pick
               up data that crosses L) file boundaries from processing
               599001 to 599002.

               revision 1.0 1997/09/18  17:30:00
               Qi Huang/RDC    (qhuang@ltpmail.gsfc.nasa.gov)
               Original development

!Team-unique Header:  
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

References and Credits:
               None

Design Notes: 
               None

!END
**************************************************************************/

{
  char                      *routine = "get_valid_L0_file";
  PGSt_SMF_status           returnStatus;
  PGSt_SMF_status           PGS_status;
  char                      msg[300];
  PGSt_integer              version_num;
  char                      file_name[PGSd_PC_VALUE_LENGTH_MAX];
  PGSt_double               scan_rate;
  char                      scan_rate_string[PGSd_PC_VALUE_LENGTH_MAX];


  returnStatus = MODIS_S_SUCCESS;

  version_num = 2;
  PGS_status = PGS_PC_GetUniversalRef(global_L0_logical, 
                                      &version_num,
                                      file_name);

  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
     sprintf(msg, "Unable to retrieve L0 file name. LUN: %u Version: %d",
             (unsigned)global_L0_logical, version_num);
     log_fmt_msg(MODIS_E_GETCONFIG_FAILED, routine, msg);
  }

  else {
     if (global_input_pointer[1][0] == '\0')
        strcpy(global_input_pointer[1], file_name);
     else
        strcpy(global_input_pointer[2], file_name);
  }

  PGS_status = PGS_IO_L0_Open(global_L0_logical, spacecraft_tag, L0_file,
                              start_time, stop_time);

  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
     returnStatus = MODIS_E_GET_VALID_L0_FILE;
     sprintf(msg, "PGS_IO_L0_Open could not open the L0 file successfully LUN: %u",
             (unsigned)global_L0_logical);
     log_fmt_msg(MODIS_E_GET_VALID_L0_FILE, routine, msg);
  }

  else {
     PGS_status = PGS_PC_GetConfigData(PC_L1A_SCAN_RATE, scan_rate_string);
     if (PGS_SMF_TestSuccessLevel(PGS_status) == PGS_TRUE)
        scan_rate = (PGSt_double) atof(scan_rate_string);
     else {
        scan_rate = 0.0;
        log_fmt_msg(MODIS_E_GETCONFIG_FAILED, routine,
                    "The scan rate could not be retrieved from pcf file");
     }

     if ((*start_time > global_last_gran_stop_time + scan_rate) ||
         ((*stop_time < global_first_gran_start_time) && 
          (global_L0_logical == PC_CURRENT_L0_PCF_ID))) {
        returnStatus = MODIS_E_GET_VALID_L0_FILE;
        sprintf(msg, "L0 data set does not correspond with granule times\n"
        "L0 Start/Stop: %f %f\nPCF Start/Stop: %f %f", *start_time, *stop_time, 
        global_first_gran_start_time, global_last_gran_stop_time);
        log_fmt_msg(MODIS_E_GET_VALID_L0_FILE, routine, msg);
     }
     /* Comment out for now; use "inhouse" tools another time.
     else
        if ((PGS_status = validate_L0_header(*L0_file)) != MODIS_S_SUCCESS) {
           returnStatus = MODIS_E_GET_VALID_L0_FILE;
           sprintf(msg, "L0 LUN: %u", (unsigned)global_L0_logical);
           log_fmt_msg(MODIS_F_L0_HEADER_VAL_FAILED, routine, msg);
        }
     */
  }
      
  return (returnStatus);

}
