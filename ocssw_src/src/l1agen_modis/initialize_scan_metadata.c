#include "L1A_prototype.h"
#include "SC_scan.h"
#include "MD_metadata.h"

void   initialize_scan_metadata (MD_SCAN_MET_t    *scan_meta)

/****************************************************************************
!C

!Description:  This function fills the scan level metadata data structures 
               with the appropriate data.

!Input Parameters: None

!Output Parameters:
               MD_SCAN_MET_t  *scan_meta    ** Scan level metadata structure **

Return Values: None

Externally Defined:  
               MD_BAD_CHECKSUM_PACKET          (MD_metadata.h)
               MD_DISCARDED_PACKET             (MD_metadata.h)
               MD_MAX_MISSING_PKTS_IN_SCAN     (MD_metadata.h)
               MD_MISSING_PACKET               (MD_metadata.h)
               MD_NO_VALID_DATA_IN_SCAN        (MD_metadata.h)
               MD_OTHER_STRING                 (MD_metadata.h)
               MD_SCAN_DATA_PRESENCE           (MD_metadata.h)
               MD_SCAN_MET_t                   (MD_metadata.h)
               SC_FILL_VALUE                   (SC_scan.h)
               TIME_FILL_VALUE                 (SC_scan.h)

Called By:
               initialize_scan
               create_missing_scans

Routines Called: None

!Revision History:
  $Log: initialize_scan_metadata.c,v $
  Revision 5.1  2004/09/23 18:53:23  seaton
  Modified fill value of time SDSs to be -2E9.
  seaton@saicmodis.com

  Revision 4.2  2003/04/30 19:18:54  vlin
  corrected the size for function strncpy.

  Revision 4.1  2003/03/07 20:52:41  vlin
  Updated after code walkthrough

  Revision 4.0  2002/12/02 20:55:53  vlin
  Called memcmp before calling memset
  vlin@saicmodis.com

               Revision 2.1 1997/09/08  19:49
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               initial revision

!Team-unique Header:

      This software is developed by the MODIS Science Data Support Team 
      for the National Aeronautics and Space Administration, 
      Goddard Space Flight Center, under contract NAS5-32373.

References and Credits: None

Design Notes: 
               The CODE below was developed in C language.

!END
**************************************************************************/

{
  int   i;

  scan_meta->scan_num = 0;
  strncpy(scan_meta->scan_type, MD_OTHER_STRING, sizeof(scan_meta->scan_type));
  scan_meta->sd_start_time = TIME_FILL_VALUE;
  scan_meta->srca_start_time = TIME_FILL_VALUE;
  scan_meta->bb_start_time = TIME_FILL_VALUE;
  scan_meta->sv_start_time = TIME_FILL_VALUE;
  scan_meta->ev_start_time = TIME_FILL_VALUE;
  scan_meta->srca_cal_mode = SC_FILL_VALUE;
  scan_meta->packet_scan_count = SC_FILL_VALUE;
  scan_meta->packet_expedited_data_flag = SC_FILL_VALUE;
  scan_meta->mirror_side = SC_FILL_VALUE;
  scan_meta->scan_qual_array[MD_SCAN_DATA_PRESENCE] = MD_NO_VALID_DATA_IN_SCAN;
  scan_meta->scan_qual_array[MD_MISSING_PACKET] = MD_MAX_MISSING_PKTS_IN_SCAN;
  scan_meta->scan_qual_array[MD_BAD_CHECKSUM_PACKET] = 0;
  scan_meta->scan_qual_array[MD_DISCARDED_PACKET] = 0;

  scan_meta->ccsds_apids[0] = SC_FILL_VALUE;
  memset(scan_meta->ccsds_apids+1, SC_FILL_VALUE, sizeof(scan_meta->ccsds_apids[0]));

  if (memcmp(scan_meta->ccsds_apids, scan_meta->ccsds_apids+1,
     sizeof(scan_meta->ccsds_apids[0]))) {         /* memset didn't work. */
     for (i=1; i<3; i++)
          scan_meta->ccsds_apids[i] = SC_FILL_VALUE;
     for (i=0; i<6; i++)
          scan_meta->frame_count_array[i] = 0; 
  }
  else {
     memset(scan_meta->frame_count_array, 0, sizeof(scan_meta->frame_count_array));
     memset(scan_meta->ccsds_apids+2, SC_FILL_VALUE, 
            sizeof(scan_meta->ccsds_apids)-2*sizeof(scan_meta->ccsds_apids[0]));
  }
}
