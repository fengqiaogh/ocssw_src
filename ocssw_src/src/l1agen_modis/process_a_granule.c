#include <time.h>
#include "PGS_IO.h"
#include "PGS_SMF.h"
#include "mapi.h"
#include "L1A_prototype.h"
#include "packet_stats.h"
#include "PGS_MODIS_35005.h"

PGSt_SMF_status  process_a_granule (PGSt_IO_L0_VirtualDataSet   L0_file,
                                    PGSt_double                 gran_start_time,
                                    PGSt_double                 gran_end_time,
                                    PCF_CONFIG_t                *pcf_config,
                                    EN_VDATA_TYPE_t             *eng_data,
                                    PH_PACKET_HEADER_t          *pkt_header,
                                    PGSt_IO_L0_Packet           *pkt,
                                    FP_QUEUE_t                  *failed_pkts)

/******************************************************************************
!C

!Description:   Function process_a_granule processes scans within a granule,
                collects global metadata about the scans for the granule, and 
                writes the scans into a L1A granule.
                  
!Input Parameters:
       L0_file                L0 file descriptor    
       gran_start_time        Start time of the L1A granule to be made    
       gran_end_time          End time of the L1A granule to be made    
       pcf_config             PCF configuration parameters

!Output Parameters:           None

!Input/Output Parameters:
       *eng_data              Engineering data    
       *pkt_header            Packet header information      
       *pkt                   Packet buffer       
       *failed_pkts           Discarded packets

Return Values:                
                MODIS_E_NULL_POINTER
                MODIS_F_PROCESSAGRANULE
                MODIS_F_UNABLE_TO_INIT_MCF
                MODIS_F_PKT_READ_FAILED
                MODIS_F_WRITE_ENG_DATA_FAIL
                MODIS_S_SUCCESS
                MODIS_W_NO_MORE_PACKETS

Externally Defined:
                ARCHIVED_METADATA               (mapi.h)
                ECSattr_names_for_all_handles   (mapi.h)
                EN_VDATA_TYPE_t                 (EN_eng_data.h)
                FALSE                           (hdf.h)
                FP_QUEUE_t                      (FP_failed_pkt_queue.h)
                INVENTORY_METADATA              (mapi.h)
                MAPIOK                          (mapi.h)
                MD_ECS_GRA_INV_MET_t            (MD_metadata.h)
                MD_L1A_SPECIFIC_MET_t           (MD_metadata.h)
                MD_SCAN_MET_t                   (MD_metadata.h)
                MECS_ARCHIVE                    (mapi.h)
                MECS_CORE                       (mapi.h)
                MODFILE                         (mapi.h)
		MODIS_E_CHECKSUM_NOT_VALID	(PGS_MODIS_35005.h)
                MODIS_E_CREATE_L1A_GRANULE      (PGS_MODIS_35005.h)
		MODIS_E_INV_PKT_SEQ_FLAG	(PGS_MODIS_35005.h)
		MODIS_E_INV_PKT_TIME		(PGS_MODIS_35005.h)
		MODIS_E_L1A			(PGS_MODIS_35005.h)
                MODIS_E_NULL_POINTER            (PGS_MODIS_35005.h)
		MODIS_E_SCANCNT_NOT_VALID	(PGS_MODIS_35005.h)
                MODIS_E_WRITE_GLOBAL_METADATA   (PGS_MODIS_35005.h)
                MODIS_F_PROCESSAGRANULE         (PGS_MODIS_35005.h)
                MODIS_F_UNABLE_TO_INIT_MCF      (PGS_MODIS_35005.h)
                MODIS_F_PKT_READ_FAILED         (PGS_MODIS_35005.h)
                MODIS_F_WRITE_ENG_DATA_FAIL     (PGS_MODIS_35005.h)
		MODIS_M_PKT_NOT_IN_SCAN		(PGS_MODIS_35005.h)
                MODIS_S_SUCCESS                 (PGS_MODIS_35005.h)
		MODIS_W_INV_PKTLEN		(PGS_MODIS_35005.h)
                MODIS_W_NO_MORE_PACKETS         (PGS_MODIS_35005.h)
                PC_MCF_PCF_ID                   (PC_pcf_info.h)
                PCF_CONFIG_t                    (PC_pcf_info.h)
                PGS_SMF_MASK_LEV_F              (PGS_SMF.h)
                PGS_TRUE                        (PGS_SMF.h)
                PGSt_double                     (PGS_TYPES.h)
                PGSt_IO_L0_Packet               (PGS_IO.h)
                PGSt_IO_L0_VirtualDataSet       (PGS_IO.h)
                PGSt_MET_all_handles            (PGS_MET.h)
                PGSt_PC_Logical                 (PGS_IO.h)
                PGSt_SMF_status                 (PGS_SMF.h)
                PH_PACKET_HEADER_t              (PH_pkt_hdr.h)
                SC_SCAN_DATA_t                  (SC_scan.h)
                SC_PIXEL_QUALITY_DATA_t         (SC_scan.h)
		stats				(packet_stats.h)
                TRUE                            (hdf.h)

Called By:
                level1a

Routines Called:
                completeMODISfile
                compute_SD_start_time
                create_L1A_granule
                end_eng_data_access_to_file
                handle_missing_scans
                initialize_global_metadata
                log_fmt_msg
                PGS_MET_Init
                PGS_SMF_TestSuccessLevel
                reset_last_valid_scan
                update_global_metadata
                process_a_scan
                write_global_metadata
                write_scan

!Revision History:
   $Log: process_a_granule.c,v $
   Revision 6.5  2014/09/10 21:44:06  jkuyper
   Resolve Bug 5555 by handling records with invalid sector values correctly.

   Revision 6.4  2012/06/22 20:02:50  kuyper
   Corrected to initialize gran_start_time_used.

   Revision 6.3  2012/02/02 23:19:12  kuyper
   Merged changes made between 5.2 and 6.1 with changes made between 5.2 and 5.2.1.2. Did it manually
     this time, hopefully with better results.

   Revision 6.2  2011/12/12 23:24:45  kuyper
   Merged changes made in 5.2.1.2 back to main branch.

   Revision 6.1  2010/08/25 16:41:35  kuyper
   Corrected to handle all valid returns from process_a_scan(), and to treat
     invalid return values as errors.
   Added initialization and printout of packet filtering statistics.

   Revision 5.2  2007/01/26 21:15:11  kuyper
   Changed to pass gran_start_time to handle_missing_scans(), allowing a full
     resolution to Bug 199.

   Revision 5.1  2006/11/09 18:00:39  kuyper
   Changed to prevent processing one more than the maximum number of scans.

   Revision 4.2  2006/01/05 16:38:56  kuyper
   Corrected to make md_handles static.

   Revision 4.1  2003/04/15 18:23:55  vlin
   calculate scan number after null pointers were checked.

   Revision 4.0  2002/12/26 20:22:20  vlin
   Updated to match  process_a_granule.pdl V4.0
   vlin@saicmodis.com

		Revision 1.3 2000/07/11
                John Seaton/GSC  (seaton@ltpmail.gsfc.nasa.gov)
                Removed revision 1.2 changes. Code has been developed to
                fix the split scan problem

                revision 1.0 1997/08/27  17:30:00
                Qi Huang/RDC    (qhuang@ltpmail.gsfc.nasa.gov)
                Original development

!Team-unique Header:
                This software is developed by the MODIS Science Data Support 
                Team for the National Aeronautics and Space Administration, 
                Goddard Space Flight Center, under contract NAS5-32373.

References and Credits:      None

Design Notes:                None

!END
******************************************************************************/

{
  const char			filefunc[] = "process_a_granule";
  PGSt_SMF_status               returnStatus = MODIS_S_SUCCESS;
  PGSt_SMF_status               PGS_status = MODIS_S_SUCCESS, L1A_status;
  int                           nscans;
  ECSattr_names_for_all_handles ECS_attr_names;

  static PGSt_MET_all_handles   md_handles;
  static int                    MET_Init_flag = TRUE;

  memset(&stats, 0, sizeof(stats));
  strcpy(ECS_attr_names[INVENTORY_METADATA], MECS_CORE);
  strcpy(ECS_attr_names[ARCHIVED_METADATA], MECS_ARCHIVE);

  if (pcf_config == NULL || eng_data == NULL || pkt_header == NULL ||
      pkt == NULL || failed_pkts == NULL) {
      log_fmt_msg(MODIS_E_NULL_POINTER, filefunc, "");
      return MODIS_E_NULL_POINTER;
  }

  nscans = ((gran_end_time - gran_start_time) / pcf_config->scan_rate) + 5;

  if (MET_Init_flag == TRUE) {
     PGS_status = PGS_MET_Init(PC_MCF_PCF_ID, md_handles);
     MET_Init_flag = FALSE;
  }

  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE) {
      returnStatus = MODIS_F_UNABLE_TO_INIT_MCF;
      log_fmt_msg(MODIS_F_UNABLE_TO_INIT_MCF, filefunc, "MCF LUN: %d",
	  PC_MCF_PCF_ID);
  }
   
  else  /*  MCF initialization is OK   */
  {
      time_t tnow;
      PGSt_double               SD_start_time;
      int			prev_scan_number;
      MODFILE			*L1A_file_ptr = NULL;
      MD_ECS_GRA_INV_MET_t	ecs_gra_inv_met; 
      MD_L1A_SPECIFIC_MET_t	l1a_specific_met;
      static int		last_granule_was_empty = TRUE;
      static PGSt_double	pre_SD_time = -1.0;
      int			Mapi_R_Status = MAPIOK;
      int			gran_start_time_used = FALSE;

     L1A_status = initialize_global_metadata(gran_start_time, gran_end_time,
	 nscans, pcf_config, &ecs_gra_inv_met, &l1a_specific_met);

     if (L1A_status != MODIS_S_SUCCESS)
        log_fmt_msg(MODIS_E_INITIALIZE_GLOBAL_MD, filefunc, "");

     reset_last_valid_scan(eng_data);


    if ((pre_SD_time < 0.0) || (last_granule_was_empty == TRUE))
      {
       pre_SD_time = gran_start_time;
       gran_start_time_used = TRUE;
      }

    compute_SD_start_time(pkt_header, &SD_start_time);
    prev_scan_number = 0;

    if ((SD_start_time >= gran_start_time) &&
	(SD_start_time < gran_end_time))
    {
       int continue_processing = TRUE;

       last_granule_was_empty = FALSE;
       if(create_L1A_granule(eng_data, nscans, &L1A_file_ptr)
	   != MODIS_S_SUCCESS)
       {
	   returnStatus = MODIS_F_PROCESSAGRANULE;
	   log_fmt_msg(MODIS_E_CREATE_L1A_GRANULE, filefunc, " ");
       }
       else while ((SD_start_time >= gran_start_time) &&
	      (SD_start_time < gran_end_time) &&
	      (l1a_specific_met.num_scans < nscans) &&
	      (continue_processing == TRUE))
       {
         if (prev_scan_number % 10 == 0) {
             time(&tnow);
             printf("scan: %d out of %d %s",
                    prev_scan_number, nscans, asctime(localtime(&tnow)));
           }

	  MD_SCAN_MET_t	scan_metadata;   

	  L1A_status = handle_missing_scans(L1A_file_ptr, SD_start_time,
		       &pre_SD_time, pcf_config->scan_rate, 
		       &prev_scan_number, &ecs_gra_inv_met, 
		       &l1a_specific_met, &gran_start_time_used,
		       gran_start_time, eng_data);

	  if (L1A_status == MODIS_F_WRITE_ENG_DATA_FAIL) {
	     returnStatus = L1A_status;
	     log_fmt_msg(MODIS_F_WRITE_ENG_DATA_FAIL, filefunc,
		 "Could not write data for the missing scans at either the\n"
		 "beginning or the middle of the granule");
	     continue_processing = FALSE;
	  }
	  else
	  {
	      SC_SCAN_DATA_t		scan_data;        
	      SC_PIXEL_QUALITY_DATA_t	pixel_quality_data;

	     if (L1A_status != MODIS_S_SUCCESS)
		log_fmt_msg(MODIS_E_HANDLE_MISSING_SCANS, filefunc, "The scan "
		    "metadata could not be written for the missing scans at\n"
		    "either the beginning or the middle of the granule");

	     L1A_status = process_a_scan(&prev_scan_number, pkt, 
			  &pcf_config->scan_rate, &SD_start_time,
			  &scan_data, &scan_metadata, eng_data, failed_pkts,
			  pkt_header, &pixel_quality_data, &L0_file);

	     switch (L1A_status) {
		case MODIS_E_CHECKSUM_NOT_VALID:
		case MODIS_E_INV_PKT_SEQ_FLAG:
		case MODIS_E_INV_PKT_TIME:
		case MODIS_E_L1A:
		case MODIS_E_SCANCNT_NOT_VALID:
		case MODIS_M_PKT_NOT_IN_SCAN:
		case MODIS_S_SUCCESS:
		case MODIS_W_INV_PKTLEN:
		case MODIS_E_INV_SECTOR:
		   break;

		case MODIS_F_PKT_READ_FAILED:
		   log_fmt_msg(MODIS_F_PKT_READ_FAILED, filefunc,
		       "A fatal error occurred while creating scan number %d ",
		       scan_metadata.scan_num);
		   continue_processing = FALSE;
		   returnStatus = L1A_status;
		   break;

		case MODIS_W_NO_MORE_PACKETS:
		   log_fmt_msg(MODIS_W_NO_MORE_PACKETS, filefunc,
		       "Unable to finish creating L1A scan number %d",
		       scan_metadata.scan_num);
		   continue_processing = FALSE;
		   returnStatus = MODIS_W_NO_MORE_PACKETS;
		   break;

		default:
		   log_fmt_msg(MODIS_E_SCAN_PROCESS, filefunc, "Scan number %d",
		       scan_metadata.scan_num);
		   returnStatus = MODIS_E_L1A;
		   continue_processing = FALSE;
	     }

	     L1A_status = write_scan(L1A_file_ptr, &scan_data,
				     &scan_metadata, &pixel_quality_data,
				     *failed_pkts, eng_data);

	     if (L1A_status == MODIS_F_WRITE_ENG_DATA_FAIL) {
		returnStatus = L1A_status;
		log_fmt_msg(MODIS_F_WRITE_ENG_DATA_FAIL, filefunc,
		    "The engineering data could not be written for\n"
		    "scan number %d", scan_metadata.scan_num);
		continue_processing = FALSE;
	     }
	     else  if (L1A_status != MODIS_S_SUCCESS)
		 log_fmt_msg(MODIS_E_WRITE_SCAN_FAIL, filefunc,
		     "Not all of the data could be written for scan number %d",
		     scan_metadata.scan_num);
	  }
	    
	  update_global_metadata(&scan_metadata, &ecs_gra_inv_met,
				 &l1a_specific_met);

       }  /*  End while */
       if(l1a_specific_met.num_scans >= nscans && SD_start_time < gran_end_time)
	   log_fmt_msg(MODIS_W_TOO_MANY_SCANS, filefunc, "");

       L1A_status = PGS_SMF_TestStatusLevel(returnStatus);

       if (L1A_status < PGS_SMF_MASK_LEV_F && L1A_file_ptr)
       {
	  L1A_status = handle_missing_scans(L1A_file_ptr, gran_end_time,
		       &pre_SD_time, pcf_config->scan_rate, 
		       &prev_scan_number, &ecs_gra_inv_met, 
		       &l1a_specific_met, &gran_start_time_used,
		       gran_start_time, eng_data);

	  if (L1A_status == MODIS_F_WRITE_ENG_DATA_FAIL) {
	      returnStatus = L1A_status;
	      log_fmt_msg(MODIS_F_WRITE_ENG_DATA_FAIL, filefunc,
		  "The engineering data could not be written successfully for\n"
		  "the missing scans at the end of the granule");
	  }
	  else if (L1A_status != MODIS_S_SUCCESS) 
	      log_fmt_msg(MODIS_E_HANDLE_MISSING_SCANS, filefunc,
		  "The scan metadata could not be written successfully for\n"
		  "the missing scans at the end of the granule");
       }
    }  /*  gran_start_time < SD_start_time < gran_end_time  */
    else
       last_granule_was_empty = TRUE;

    if(L1A_file_ptr)
    {
        const long int                      num_handles = 2L;

	L1A_status = write_global_metadata(L1A_file_ptr, md_handles,
			  &ecs_gra_inv_met, &l1a_specific_met);
	if (L1A_status != MODIS_S_SUCCESS)
	    log_fmt_msg(MODIS_E_WRITE_GLOBAL_METADATA, filefunc,
		"The global metadata could not be written successfully for\n"
		"granule number %s", ecs_gra_inv_met.parametervalue_1);
	L1A_status = end_eng_data_access_to_file (L1A_file_ptr, eng_data);
	if (L1A_status == FAIL) {
	    returnStatus = MODIS_F_PROCESSAGRANULE;
	    log_fmt_msg(MODIS_F_PROCESSAGRANULE, filefunc, "Either the "
		"Discarded Packets vdata or the Engineering Data vdata could\n"
		"not end its access to the L1A granule for granule number %s",
		ecs_gra_inv_met.parametervalue_1);
	}
	   
	Mapi_R_Status = completeMODISfile(&L1A_file_ptr, md_handles, 
					  ECS_attr_names, num_handles);
    }
    else
	log_fmt_msg(MODIS_W_CREATE_L1A_GRANULE, filefunc, "%s %s",
	    ecs_gra_inv_met.rangebeginningdate,
	    ecs_gra_inv_met.rangebeginningtime);

    if (Mapi_R_Status != MAPIOK) {
	returnStatus = MODIS_F_PROCESSAGRANULE;
	log_fmt_msg(MODIS_F_PROCESSAGRANULE, filefunc,
	    "The L1A granule could not be closed for granule number %s",
	    ecs_gra_inv_met.parametervalue_1);
    }
  }      /*  MCF Initialization check    */ 
  print_stats();

  return returnStatus;
}
