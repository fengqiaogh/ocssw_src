#include "PGS_TD.h"
#include "PGS_IO.h"
#include "PGS_SMF.h"
#include "PGS_PC_9.h"
#include "PGS_IO_L0.h"
#include "PGS_TYPES.h"
#include "PH_pkt_hdr.h"
#include "EN_eng_data.h"
#include "PC_pcf_info.h"
#include "PGS_MODIS_35005.h"
#include "L1A_prototype.h"


PGSt_SMF_status  initialize_level1a (PCF_CONFIG_t        *pcf_config,
                                     EN_VDATA_TYPE_t     *eng_data,
                                     PGSt_IO_L0_Packet   *pkt,
                                     PH_PACKET_HEADER_t  *pkt_header,
                                     PGSt_IO_L0_VirtualDataSet *L0_file)

/*************************************************************************
!C

!Description:  This function will call various routines to get the PCF
               parameters, open the L0 file, read the Engineering Data
               structures from a file. It will also set the L0 file 
               pointer to a point before the start of the first L1A granule 
               to be produced, and also preload the eng_data structure. 

!Input Parameters: None

!Output Parameters:
        PCF_CONFIG_t    *pcf_config
                                  ** Configuration parameters        **
                                  ** extracted from the PCF          **
        EN_VDATA_TYPE_t  *eng_data
                                  ** The structure containing the    **
                                  ** de-commutated engineering data  **
        PGSt_IO_L0_Packet   *pkt  ** Packet buffer which contains    **
                                  ** packed packet                   **
        PH_PACKET_HEADER_t  *pkt_header
                                  ** Buffer that contains the pkt    **
                                  ** header info of the packet       **
        PGSt_IO_L0_VirtualDataSet  *L0_file
                                  ** The L0 file that contains       **
                                  ** packet information              **

Return Values: 
        MODIS_E_NULL_POINTER
        MODIS_F_GETREF_FAILED
        MODIS_F_INVAL_L0_FILE_RET
        MODIS_F_L0_SETSTART_FAILED
        MODIS_F_NO_PCF_CONFIG_DATA
        MODIS_F_PKT_READ_FAILED
        MODIS_S_SUCCESS
        MODIS_W_NO_MORE_PACKETS
        MODIS_W_PRELOAD_ENG_DATA 

Externally Defined:  
        EN_VDATA_TYPE_t              (EN_eng_data.h)
        global_first_gran_start_time (L1A_prototype.h)
        global_last_gran_stop_time   (L1A_prototype.h)
        global_L0_logical            (PC_pcf_info.h)
        MODIS_E_NULL_POINTER         (PGS_MODIS_35005.h)
        MODIS_F_GETREF_FAILED        (PGS_MODIS_35005.h)
        MODIS_F_INVAL_L0_FILE_RET    (PGS_MODIS_35005.h)
        MOIDS_F_INVALID_ENG_DATA_LIST(PGS_MODIS_35005.h)
        MODIS_F_L0_SETSTART_FAILED   (PGS_MODIS_35005.h)
        MODIS_F_NO_PCF_CONFIG_DATA   (PGS_MODIS_35005.h)
        MODIS_F_PKT_READ_FAILED      (PGS_MODIS_35005.h)
        MODIS_S_SUCCESS              (PGS_MODIS_35005.h)
        MODIS_W_NO_MORE_PACKETS      (PGS_MODIS_35005.h)
        MODIS_W_PRELOAD_ENG_DATA     (PGS_MODIS_35005.h)
        PC_CURRENT_L0_PCF_ID         (PC_pcf_info.h)
        PC_PRIOR_L0_PCF_ID           (PC_pcf_info.h)
        PCF_CONFIG_t                 (PC_pcf_info.h)
        PGSPC_E_DATA_ACCESS_ERROR    (PGS_PC_9.h)
        PGSPC_W_NO_REFERENCE_FOUND   (PGS_PC_9.h)
        PGSd_EOS_AM1                 (PGS_TD.h)
        PGSd_FILE_PATH_MAX           (PGS_PC.h)
        PGSt_double                  (PGS_types.h)
        PGSt_integer                 (PGS_types.h)
        PGSt_IO_L0_Packet            (PGS_IO.h)
        PGSt_IO_L0_VirtualDataSet    (PGS_IO_LO.h)
        PGSt_Logical                 (PGS_types.h)
        PGSt_SMF_status              (PGS_SMF.h)
        PH_PACKET_HEADER_t           (PH_pkt_hdr.h) 

Called By:     level1a 

Routines Called:
        PGS_PC_GetReference
        get_pcf_config_data
        compute_global_time_offsets
        get_valid_L0_file
        load_eng_data
        parse_eng_data_list
        set_start_position
        log_fmt_msg

!Revision History:
  $Log: initialize_level1a.c,v $
  Revision 4.2  2003/01/10 15:43:14  vlin
  Check for NULL pointers in the beginning.

  Revision 4.1  2002/10/09 17:49:24  vlin
  Updated after code walkthrough.

  Revision 3.1  2002/09/16 21:26:02  vlin
  Updated according to initialize_level1a.pdl revision 3.2

               Revision 2.2  2000/07/14
               John Seaton/SAIC/GSC  (seaton@ltpmail.gsfc.nasa.gov)
               Added changes for Aqua instrument.

               Revision 2.1  1997/08/25  14:20 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Original Code Development.

               Revision 2.0  1997/07/15  14:40 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Recreated Module per Version2 development.

               Revision 1.0  1997/06/18  16:40 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Baseline from Version 1.

!Team-unique Header:

     This software is developed by the MODIS Science Data Support 
     Team for the National Aeronautics and Space Administration, 
     Goddard Space Flight Center, under contract NAS5-32373.

!END
***********************************************************************/

{
  char    routine[] = "initialize_level1a";
  char    msg[300];
  PGSt_SMF_status  returnStatus = MODIS_S_SUCCESS;
  PGSt_SMF_status  tempStatus;    /* SMF-style message returned by function */
  PGSt_double    L0_start_time, L0_stop_time; 
                                  /* L0 start time and L0 stop time (TAI)   */ 
                                  /* returned by PGS_IO_L0_Open function    */
  PGSt_double    preload_start_time = 0.0;
                                  /* Time to start loading engineering data */
                                  /* prior to first product                 */
  PGSt_integer   version = 1;     /* Version number needed by               */
                                  /* PGS_PC_GetReference                    */
  char           refID[PGSd_PC_FILE_PATH_MAX];
                                  /* Char buffer needed by                  */
                                  /* PGS_PC_GetReference                    */

  if (pcf_config == NULL || eng_data == NULL || pkt == NULL || 
      pkt_header == NULL || L0_file == NULL) {
      log_fmt_msg(MODIS_E_NULL_POINTER, routine, "");
      return MODIS_E_NULL_POINTER;
  }

  returnStatus = get_pcf_config_data(pcf_config);
  global_first_gran_start_time = pcf_config->first_gran_start_time;
  global_last_gran_stop_time = pcf_config->last_gran_stop_time;
  if (returnStatus == MODIS_S_SUCCESS) {
     compute_global_time_offsets(pcf_config->scan_rate);
     tempStatus = PGS_PC_GetReference(PC_PRIOR_L0_PCF_ID, &version, refID);
     if (tempStatus == PGSPC_E_DATA_ACCESS_ERROR) {
        returnStatus = MODIS_F_GETREF_FAILED;
        sprintf(msg, "LUN Number: %d", PC_PRIOR_L0_PCF_ID);
        log_fmt_msg (MODIS_F_GETREF_FAILED, routine, msg);
     }
     else {
        if (tempStatus == PGS_S_SUCCESS)
            global_L0_logical = PC_PRIOR_L0_PCF_ID;
        else if (tempStatus == PGSPC_W_NO_REFERENCE_FOUND) 
           global_L0_logical = PC_CURRENT_L0_PCF_ID;

        returnStatus = parse_eng_data_list(eng_data);
        if (returnStatus != MODIS_S_SUCCESS) 
           log_fmt_msg(MODIS_F_INVALID_ENG_DATA_LIST, routine,
           "error occurred trying to parse the engineering data list file");
        else if (pcf_config->instrument != eng_data->instrument) {
           returnStatus = MODIS_F_INVALID_ENG_DATA_LIST;
           sprintf(msg, "PCF instrument = %d, LUT instrument = %d",
                   pcf_config->instrument, eng_data->instrument);
           log_fmt_msg(MODIS_F_INVALID_ENG_DATA_LIST, routine, msg);
        }
        else if (strcmp(pcf_config->lutrevision, eng_data->revision)) {
           returnStatus = MODIS_F_INVALID_ENG_DATA_LIST; 
           sprintf(msg, "Expected LUT revision = %s, actual revision = %s",
                  pcf_config->lutrevision, eng_data->revision);
           log_fmt_msg(MODIS_F_INVALID_ENG_DATA_LIST, routine, msg);
        }
        else
           returnStatus = get_valid_L0_file(PGSd_EOS_AM1, L0_file,
                          &L0_start_time, &L0_stop_time);

        if (returnStatus == MODIS_S_SUCCESS) {
           preload_start_time = pcf_config->first_gran_start_time - (65 * 1.024);
           returnStatus = set_start_position (*L0_file, 
                          &preload_start_time, L0_start_time, L0_stop_time);
           if (returnStatus == MODIS_S_SUCCESS) {
              returnStatus = load_eng_data(pcf_config->scan_rate, eng_data,
                             pkt, pkt_header, *L0_file);
              if (returnStatus == MODIS_W_PRELOAD_ENG_DATA)
                 log_fmt_msg (MODIS_W_PRELOAD_ENG_DATA, routine,
                     "Warning: Could not preload any Eng Data");
              else if (returnStatus == MODIS_W_NO_MORE_PACKETS)
                 log_fmt_msg(MODIS_W_NO_MORE_PACKETS, routine, "Reached "
                 "the end of L0 data while trying to preload eng data");
              else if (returnStatus == MODIS_F_PKT_READ_FAILED)
                 log_fmt_msg(MODIS_F_PKT_READ_FAILED, routine,
                   "error reading packet data during engineering preload");
           }
           else {
              returnStatus = MODIS_F_L0_SETSTART_FAILED;
              log_fmt_msg(MODIS_F_L0_SETSTART_FAILED, routine, "");
           }
        }
        else {
           returnStatus = MODIS_F_INVAL_L0_FILE_RET;
           sprintf(msg, "L0 LUN: %u", (unsigned)global_L0_logical);
           log_fmt_msg(MODIS_F_INVAL_L0_FILE_RET, routine, msg);
        } /* else get_valid_L0_file fails */
     }
  }  
  else {
     returnStatus = MODIS_F_NO_PCF_CONFIG_DATA;
     log_fmt_msg(MODIS_F_NO_PCF_CONFIG_DATA, routine, "error occured in "
     "trying to get PCF config data values Granule Length or Scan Rate");
  }

  return returnStatus;

}
