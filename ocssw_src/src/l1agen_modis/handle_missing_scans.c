#include "PGS_SMF.h"
#include "mapi.h"
#include "EN_eng_data.h"
#include "L1A_prototype.h"
#include "MD_metadata.h"
#include "PGS_MODIS_35005.h"
#include "SC_scan.h"

PGSt_SMF_status  handle_missing_scans (MODFILE                 *L1A_file_ptr,
                                       PGSt_double             SD_start_time, 
                                       PGSt_double             *pre_SD_time,
                                       PGSt_double             scan_rate,
                                       int                     *scan_number,
                                       MD_ECS_GRA_INV_MET_t    *ecs_gran_meta,
                                       MD_L1A_SPECIFIC_MET_t   *L1A_specific_meta,
                                       PGSt_SMF_boolean        *gran_start_time_used,
				       PGSt_double	       gran_start_time,
                                       EN_VDATA_TYPE_t         *eng_data )
/*
!C******************************************************************************

!Description:   Function handle_missing_scans determines if there are any missing
                scans in this granule and if so, creates those scans.
                  
!Input Parameters:
                PGSt_double           SD_start_time       ** SD start time for the
                                                              current packet       **
                PGSt_double           scan_rate           ** Scan rate              **
                EN_VDATA_TYPE_t       *eng_data           ** Engineering data      **
 
!Output Parameters:
                None

!Input/Output Parameters:
                PGSt_double           *pre_SD_time        ** SD start time for the
                                                              last packet          *
                int16                 *scan_number        ** Previous scan number  **
                MD_ECS_GRA_INV_MET_t, *ecs_gran_meta      ** ECS granule metadata  **
                MD_L1A_SPECIFIC_MET_t *L1A_specific_meta  ** L1A specific metadata **
                PGSt_SMF_boolean      *gran_start_time_used ** Flag indicating if
                                                                the SD start time for
                                                                the last packet is an
                                                                estimate           **
                MODFILE               *L1A_file_ptr       ** L1A file pointer      **

Return Values:                
                MODIS_S_SUCCESS
		MODIS_E_ARRAY_OUTPUT_ERR
                MODIS_E_HANDLE_MISSING_SCANS
                MODIS_F_WRITE_ENG_DATA_FAIL

Externally Defined:
                EN_VDATA_TYPE_t                   (EN_eng_data.h)
                MD_ECS_GRA_INV_MET_t              (MD_metadata.h)
                MD_L1A_SPECIFIC_MET_t             (MD_metadata.h)
                MD_SCAN_MET_t                     (MD_metadata.h)
		PGS_TRUE			  (PGS_SMF.h)
		PGS_FALSE			  (PGS_SMF.h)
                PGSt_double			  (PGS_TYPES.h)
		PGSt_SMF_boolean		  (PGS_SMF.h)
                SC_SCAN_RATE_TOLERANCE            (SC_scan.h)

Called By:
                process_a_granule

Routines Called:
                create_missing_scans
                write_scan_metadata
                write_eng_data
                update_global_metadata
                log_fmt_msg

!Revision History:
		$Log: handle_missing_scans.c,v $
		Revision 5.1  2007/01/26 20:56:31  kuyper
		Fixed to resolve Bug 199 by preventing the insertion of too many missing
		  scans.

		Revision 1.3  2001/04/13 18:09:20  seaton
		Changed raw_mir_enc to type uint16 and fixed numerous prologs.

		Revision 1.2  2000/02/08 17:19:49  seaton
		Bug fixes and error message corrections version 2.1.4

		Revision 1.1  1999/08/10 14:05:25  seaton
		Initial revision

                revision 1.0 1997/08/28  17:30:00
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

!END****************************************************************************
*/
{
  char                               *routine = "handle_missing_scans";
  char                               msg[300];
  PGSt_SMF_status                    returnStatus;
  PGSt_SMF_status                    L1A_status;

  int                                num_missing_scans;
  int                                i;
  MD_SCAN_MET_t                      scan_meta;
  SC_PIXEL_QUALITY_DATA_t            pix_qual_data;

  returnStatus = MODIS_S_SUCCESS;

  if ((SD_start_time - *pre_SD_time) / (scan_rate + SC_SCAN_RATE_TOLERANCE) > 0)
    {
      num_missing_scans = ((SD_start_time - gran_start_time) /
	  (scan_rate+SC_SCAN_RATE_TOLERANCE)) - 0.5
	  - L1A_specific_meta->num_scans;
      if(num_missing_scans < 0)
	  num_missing_scans = 0;
      if (*gran_start_time_used == PGS_TRUE)
        *pre_SD_time = SD_start_time - (scan_rate * (num_missing_scans + 1));
       
      i = 0;
      while ((i < num_missing_scans) && 
             (returnStatus != MODIS_F_WRITE_ENG_DATA_FAIL))
      {
        create_missing_scans(*scan_number, 
                             scan_rate,
                             pre_SD_time,
                             &scan_meta, 
                             &pix_qual_data);
        
        *scan_number = scan_meta.scan_num;

        L1A_status = write_scan_metadata(L1A_file_ptr,&scan_meta);
        if (L1A_status != MODIS_S_SUCCESS)
          {
            returnStatus = MODIS_E_HANDLE_MISSING_SCANS;
            sprintf(msg,"The scan metadata could not be written succesfully for missing scan number %d",
                    scan_meta.scan_num);
            log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
          }

        L1A_status = write_pix_qual(L1A_file_ptr, &pix_qual_data, scan_meta.scan_num);

        if (L1A_status != MODIS_S_SUCCESS)
          {
            returnStatus = MODIS_E_HANDLE_MISSING_SCANS;
            sprintf(msg,"The pixel quality data could not be written succesfully for missing scan number %d",
                    scan_meta.scan_num);
            log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
          }

        L1A_status = write_eng_data(eng_data);
        
        if (L1A_status != MODIS_S_SUCCESS)
          {
            returnStatus = L1A_status;
            sprintf(msg,"The engineering data could not be written succesfully for missing scan number %d",
                    scan_meta.scan_num);
            log_fmt_msg(MODIS_F_WRITE_ENG_DATA_FAIL, routine, msg);
          }

        update_global_metadata(&scan_meta,ecs_gran_meta,L1A_specific_meta);
        i++;
    }
  }
 

  *pre_SD_time = SD_start_time;
  *gran_start_time_used = PGS_FALSE;

  return (returnStatus);

}
