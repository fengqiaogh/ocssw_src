#include <math.h>
#include "EN_eng_data.h"
#include "FP_failed_pkt_queue.h"
#include "hdfi.h"
#include "L1A_prototype.h"
#include "MS_misc.h" 
#include "PC_pcf_info.h"
#include "PGS_IO_L0.h"
#include "PGS_MODIS_35005.h"
#include "PGS_PC.h"
#include "PGS_SMF.h"
#include "PGS_TD.h"
#include "PGS_TYPES.h"
#include "PH_pkt_hdr.h"
#include "VU_vdata_utility.h"
#include "version.h"
#include "smfio.h"
 
/* 
 * Define Global Variables 
 */

char            global_input_pointer[MAX_INPUTS][PGSd_PC_VALUE_LENGTH_MAX];   
                             /* Input file names: Eng Data List file name,
                                prior L0 file name, current L0 file name  */
int             global_VU_ID_TABLE_READY = FALSE;
int32           global_H_ID;
                   /* This is shared between the functions in this file so the 
                      user of this package doesn't have to deal with it.  */
PGSt_double     global_first_gran_start_time; /* first granule start time */
PGSt_double     global_last_gran_stop_time;   /* last granule stop time   */
PGSt_double     global_time_offset_array[5];  
                             /* The time offsets for each sector from 
                                the beginning of the scan (ie. SD sector) */
PGSt_PC_Logical global_L0_logical;                 /* L0 file unit number */
VU_ID_TABLE     global_VU_VDATA_ID[VU_MAX_NUMBER_OF_VDATAS];
                            /* used to store the Vdata IDs so that functions 
                                        need only refer to Vdatas by name */

int  main () 
/******************************************************************************
!C

!Description:  This is the main function of the MODIS Level 1A processing
               system. 

!Input Parameters:   None

!Output Parameters:  None

Return Values: 
                MS_MT_EXIT_FATAL_FAILURE 
                MS_MT_EXIT_SUCCESS      

Externally Defined:  
                EN_NUM_VDATAS                   (EN_eng_data.h)
                EN_VDATA_TYPE_t                 (EN_eng_data.h)
                FALSE                           (hdfi.h)
                FP_QUEUE_t                      (FP_failed_pkt_queue.h)
                int32                           (hdfi.h)
                MODIS_E_INITIALIZATION_FAILED   (PGS_MODIS_35005.h)
		MODIS_E_L1A			(PGS_MODIS_35005.h)
                MODIS_E_MALLOC_FAILED           (PGS_MODIS_35005.h)
                MODIS_E_NO_SCANS_IN_PRODUCT     (PGS_MODIS_35005.h)
                MODIS_F_PROCESSAGRANULE         (PGS_MODIS_35005.h)
                MODIS_U_L1A_BEGIN               (PGS_MODIS_35005.h)
                MODIS_U_L1A_END                 (PGS_MODIS_35005.h)
                MS_MT_EXIT_FATAL_FAILURE        (MS_misc.h)
                MS_MT_EXIT_SUCCESS              (MS_misc.h)
                PCF_CONFIG_t                    (PC_pcf_info.h)
                PD_PKT_BUF_MAX                  (PD_pkt_data.h)
                PGSd_PC_VALUE_LENGTH_MAX        (PGS_PC.h)
                PGS_SMF_MAX_MSGBUF_SIZE         (PGS_SMF.h)
                PGS_SMF_MASK_LEV_F              (PGS_SMF.h)
                PGS_TRUE                        (PGS_SMF.h)
                PGSt_double                     (PGS_TYPES.h)
                PGSt_IO_L0_VirtualDataSet       (PGS_IO_L0.h)
                PGSt_PC_Logical                 (PGS_TYPES.h)
                PGSt_SMF_status                 (PGS_SMF.h)
                PH_PACKET_HEADER_t              (PH_pkt_hdr.h)
                PROCESSVERSION                  (version.h)
                VU_ID_TABLE                     (VU_vdata_utility.h)
                VU_MAX_NUMBER_OF_VDATAS         (VU_vdata_utility.h)

Called By:     ECS

Routines Called:
               PGS_SMF_TestSuccessLevel
               initialize_level1a           
               process_a_granule
               make_queue
               free_queue            
               log_fmt_msg               
               close_processing_run         

!Revision History:
               $Log: level1a.c,v $
               Revision 5.2  2006/01/04 17:12:19  kuyper
               Changed to avoid daylight savings time problems by using a method of
                 skipping leap seconds that works directly on the asciiUTC representation.

               Revision 5.1  2006/01/04 00:02:44  kuyper
               Change to use minutes rather than seconds for granule sizes, for better
                 handling of leap seconds, resolving MODur00122 and GSFcd02888.

               Revision 4.3  2003/11/12 21:06:50  kuyper
               Changed to use MAX_INPUTS macro.

               Revision 4.2  2003/04/24 18:25:02  vlin
               updated after code walkthrough

               Revision 4.1  2003/02/24 20:21:04  kuyper
               Removed inappropriate uses of NULL.

               Revision 4.0  2003/01/08 19:35:56  vlin
               updated according to level1a.pdl revision 4.2
               vlin@saicmodis.com

               Revision 2.2  2001/01/04
               John Seaton  (seaton@ltpmail.gsfc.nasa.gov)
               Added new mnemonic for L1A version, not PGE version.

               Revision 2.1  2000/07/17
               John Seaton  (seaton@ltpmail.gsfc.nasa.gov)
               Added changes for Aqua spacecraft

               Revision 2.0  1997/08/08  15:15 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Remodeled PDL per Version 2 development.

!Team-unique Header:

          This software is developed by the MODIS Science Data Support Team 
          for the National Aeronautics and Space Administration, 
          Goddard Space Flight Center, under contract NAS5-32373.

References and Credits:

Design Notes: 

!END
******************************************************************************/

{
  PGSt_SMF_status  Status;       
  int              exitStatus = MS_MT_EXIT_SUCCESS;  
  char             routine[] = "level1a";
  PCF_CONFIG_t     pcf_config;
  PGSt_double      gran_start_time;      
  PGSt_double      gran_stop_time;      
  PGSt_IO_L0_VirtualDataSet   L0_file;     /* File containing L0 data        */
  PGSt_IO_L0_Packet  pkt[PD_PKT_BUF_MAX];  /* Packet of L0 data              */
  PH_PACKET_HEADER_t pkt_header;           /* Unpacked header of L0 packet   */
  EN_VDATA_TYPE_t  eng_data[EN_NUM_VDATAS];/* Engineering vdata structure    */
  FP_QUEUE_t         failed_pkts;          /* Discarded packets queue        */

  setlinebuf (stdout);
  printf("MODIS L1A version %s, built %s %s\n",
         PROCESSVERSION, __DATE__, __TIME__);

  log_fmt_msg(MODIS_U_L1A_BEGIN, routine, "L1A version %s built on %s, at %s",
      PROCESSVERSION, __DATE__, __TIME__);
  Status = initialize_level1a(&pcf_config, eng_data, pkt, 
                              &pkt_header, &L0_file);

  if (PGS_SMF_TestStatusLevel(Status) < PGS_SMF_MASK_LEV_F) {
     failed_pkts = make_queue();
     if (failed_pkts != NULL)
     {
        gran_start_time = pcf_config.first_gran_start_time;
        gran_stop_time  = pcf_config.first_gran_start_time +
                          pcf_config.gran_time_length*60.0;

        while (fabs(gran_start_time - pcf_config.last_gran_stop_time) > 0.1)
	{
	    char asciiUTC[TIMECODEASIZE];

	    //gran_stop_time = gran_start_time + pcf_config.gran_time_length*60.0
		//+ 1.0;	/* to cover leapseconds, if present. */

	    if(PGS_SMF_TestErrorLevel(PGS_TD_TAItoUTC(gran_stop_time, asciiUTC))
		!= PGS_FALSE)
	    {
		log_fmt_msg(MODIS_E_L1A, routine, "PGS_TD_TAItoUTC(%.1f)",
		    gran_start_time);
		exitStatus = MS_MT_EXIT_FATAL_FAILURE;
	    }
	    else
	    {
	       /* Cancel out the +1.0, if no leap second was present, by
		* setting the last seconds digit to 0.
		*/
              // JMG: This assumes that before adding the leap second above,
              //      the stop_time had a 0 in the first digit to the left
              //      of the decimal point.  This assumption does not
              //      always hold.

              //      I have modified the source code by commenting out
              //      the addition of the leapsec above and removing the
              //      call to UTCtoTAI just below.

	       asciiUTC[sizeof "2001-01-01T00:00:"] = '0';

               /*
	       if(PGS_SMF_TestErrorLevel(PGS_TD_UTCtoTAI(asciiUTC,
		   &gran_stop_time)) != PGS_FALSE)
               */
               if(0)
	       {
		  log_fmt_msg(MODIS_E_L1A, routine, "PGS_TD_UTCtoTAI(\"%s\")",
		      asciiUTC);
		  exitStatus = MS_MT_EXIT_FATAL_FAILURE;
	       }
	       else if(PGS_SMF_TestStatusLevel(process_a_granule(L0_file,
		   gran_start_time, gran_stop_time, &pcf_config, eng_data,
		   &pkt_header, pkt, &failed_pkts)) != PGS_SMF_MASK_LEV_F)
	       {
		   /** Advance the start/stop times for the next granule and **/
		   /** check to see if the last granule was just created     **/
		  gran_start_time = gran_stop_time;
		  if (global_input_pointer[2][0] != '\0') {
		     memset(global_input_pointer[1], '\0', 
			    sizeof(global_input_pointer[1]));
		     strcpy(global_input_pointer[1], global_input_pointer[2]);
		     memset(global_input_pointer[2], '\0', 
			    sizeof(global_input_pointer[2]));
		  }
	       }
	       else {
		  log_fmt_msg(MODIS_F_PROCESSAGRANULE, routine,
		      "Fatal Error: Could not process a granule");
		  gran_start_time = pcf_config.last_gran_stop_time;
		  exitStatus = MS_MT_EXIT_FATAL_FAILURE;
	       }
	    }
	}  /*  end while  */
        free_queue (failed_pkts);
     }
     else {
        log_fmt_msg(MODIS_E_MALLOC_FAILED, routine,
                   "unable to allocate memory for the discarded packets queue");
        exitStatus = MS_MT_EXIT_FATAL_FAILURE;
     }
  }
  else {
     log_fmt_msg(MODIS_E_INITIALIZATION_FAILED, routine, 
                 "Routine initialize_level1a failed");
     exitStatus = MS_MT_EXIT_FATAL_FAILURE;
  }

  close_processing_run(L0_file);
  log_fmt_msg(MODIS_U_L1A_END, routine, "L1A return code = %d", exitStatus);
  return exitStatus;

}
