#include "PGS_MODIS_35251.h"
#include "GEO_input.h"
#include "L1a_data.h"
#include "smfio.h"

int GEO_read_L1Ascan_metadata(
        MODFILE  * const l1a_file,
	int const number_of_scans,
	l1a_metadata_struct const * const granule_metadata,
	frame_data_struct frame_data[MAX_SCAN_NUMBER],
	int16 mirr_side[MAX_SCAN_NUMBER]
	)

/*
!C***********************************************************************
 
!Description:   
		Subroutine in the input group of the Level-1A geolocation
                software to read L1A scan-level metadata from the L1A
		product file using the MAPI.

!Input Parameters:
                l1a_file - MAPI structure for l1a file
                number_of_scans - number of scans in the granule
                granule_metadata - l1a specific granule metadata

!Output Parameters:
		frame_data_struct frame_data - L1A scan start times and
			scan frame sizes
		mirr_side - mirror side for each scan

Return parameter:
		SUCCESS if all specified data
		are successfully retrieved,
		FAIL otherwise.

Global variables:
		None.

Routins called:
               modsmf - writes error status messages to log
	       int getMODISarray - MAPI routine to read from an array           
Called by:
               GEO_prepare_L1Ascan_metadata()

!Revision History:
 * $Log: GEO_read_L1Ascan_metadata.c,v $
 * Revision 4.1  2003/02/21 22:46:47  kuyper
 * Corrected to use void* pointers with %p format code.
 *
 * Revision 3.1  2002/06/13 22:53:22  kuyper
 * Removed unnecessary NCSA acknowledgement.
 *
 * Revision 2.6  2000/07/26 19:31:31  vlin
 * Modified the content of the error message for calls to
 * function "modsmf"
 *
 * Revision 2.5  2000/05/05  15:30:20  lma
 * updated after walkthrough
 *
 * Revision 2.4  2000/05/02  20:48:48  lma
 * added one SDS(Scan Type) per CCR507
 *
 * Revision 2.3  1999/02/02  18:07:08  kuyper
 * Make sprintf() format string consistent with the argument list.
 *
 * Revision 2.2  1998/02/08  23:00:57  jjb
 * Corrected defect described in DDTS MODxl00596:
 * Check number of Space View frames in each scan
 * against max_sv_frames.
 *
 * Revision 2.1  1997/10/21  18:16:22  kuyper
 * Returned from ClearCase
 *
 * Revision /main/GEO_V2_DEV/3 1997/10/03 kuyper
 * Corrected MODIS_W_GEO_FRAMENO_INPUT messages.
 *
 * Revision /main/GEO_V2_DEV/2  1997/08/29  Ding
 * Added in SDS sci_state and sci_abnorm.
 * ding@ltpmail.gsfc.nasa.gov
 *
 * Revision 1.4.1.1  1997/07/21  22:20:17  kuyper
 * Merged in out-of-sequence changes.
 *
 * Revision 1.4  1997/07/21  16:24:34  kuyper
 * Baselined Version 1
 *
 *Parallel development:
 * Revision 1.3  1997/06/02  18:01:11  kuyper
 * Merged seed files.
 *
 * Revision 1.3  1997/03/26  18:13:34  fhliang
 * Initial revision of SDST delivery of GEO_read_L1Ascan_metadata.c.
 *
		Revision 1.2  1997/01/15 17:15:49  kuyper
		Changed to use M01 M-API macros indirectly.
		Corrected scan_id, frame_count to int16.

		James Kuyper kuyper@ltpmail.gsfc.nasa.gov

		Revision 1.1  1997/01/08 18:29:35  mikej
		Initial revision

		
!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.
                
!END****************************************************************************
*/

{

/* Note: Earth_Frames, SD_Frames, and SV_Frames are not all capitals to avoid
conflict with string literals in GEO_product.h */
#define Earth_Frames 1
#define SD_Frames 2 
#define SV_Frames 5 
/* Magic number, the index of the "Mirror side" entry in the Scan_metadata
array. If additional entries are added this must be changed */
#define MIRR_SIDE_INDEX 8
#define NUM_METADATA (sizeof(Scan_metadata) / sizeof(Scan_metadata[0]))

  static float64 *EV_start;
  static float64 *SD_start;
  static float64 *SV_start;
  typedef int32 scan_quality_t[NUM_L1A_QUALITY_FLAGS];
  static scan_quality_t *scan_quality;
  static int16 *scan_id;
  typedef int16 frame_count_t[NUM_L1A_SECTOR_VIEWS];
  static frame_count_t *frame_count;
  static int8  *sci_state;
  static int8  *sci_abnorm;
  typedef char8 scan_type_t[SCAN_TYPE_LEN];
  static scan_type_t *scan_type;

  // allocate static memory
  static int firstRun = 1;
  if(firstRun) {
    firstRun = 0;
    EV_start = (float64*) calloc(MAX_SCAN_NUMBER, sizeof(float64));
    SD_start = (float64*) calloc(MAX_SCAN_NUMBER, sizeof(float64));
    SV_start = (float64*) calloc(MAX_SCAN_NUMBER, sizeof(float64));
    scan_quality = (scan_quality_t*) calloc(MAX_SCAN_NUMBER, sizeof(scan_quality_t));
    scan_id = (int16*) calloc(MAX_SCAN_NUMBER, sizeof(int16));
    frame_count = (frame_count_t*) calloc(MAX_SCAN_NUMBER, sizeof(frame_count_t));
    sci_state = (int8*) calloc(MAX_SCAN_NUMBER, sizeof(int8));
    sci_abnorm = (int8*) calloc(MAX_SCAN_NUMBER, sizeof(int8));
    scan_type = (scan_type_t*) calloc(MAX_SCAN_NUMBER, sizeof(scan_type_t));
  }

  struct scan_meta{
	  char *name; /* name of SDS object in the L1A file */
	  long int dims[2]; /* extent of locations in SDS */
	  void * data; /* pointer to the data */
  } Scan_metadata[] = {
	  {SCAN_NUMBER, {MAX_SCAN_NUMBER, 0L}, scan_id},
	  {EARTH_SECTOR_FRAMES, {MAX_SCAN_NUMBER, NUM_L1A_SECTOR_VIEWS}
	      , frame_count},
	  {SCAN_START_TIME, {MAX_SCAN_NUMBER, 0L}, EV_start},
	  {SD_START_TIME, {MAX_SCAN_NUMBER, 0L}, SD_start},
	  {SV_START_TIME, {MAX_SCAN_NUMBER, 0L}, SV_start},
	  {SCAN_QUALITY_ARRAY, {MAX_SCAN_NUMBER, NUM_L1A_QUALITY_FLAGS},
	       scan_quality},
	  {SCIENCE_STATE, {MAX_SCAN_NUMBER, 0L}, sci_state},
	  {SCIENCE_ABNORM, {MAX_SCAN_NUMBER, 0L}, sci_abnorm},
	  {MIRROR_SIDE, {MAX_SCAN_NUMBER, 0L}, NULL},
	  {SCAN_TYPES, {MAX_SCAN_NUMBER, SCAN_TYPE_LEN}, scan_type}
  };

  long start[2] = {0L, 0L};	/* offset for getMODISarray */
  int retval=SUCCESS;
  int i, scan;	/* loop indices */
  char msgbuf[PGS_SMF_MAX_MSGBUF_SIZE];	/* scratch string buffer */

/*  Check input values */

  if((l1a_file == NULL) || (granule_metadata == NULL) || (number_of_scans < 0)
    || (number_of_scans > MAX_SCAN_NUMBER) || (frame_data == NULL) || 
    (mirr_side == NULL))
  {
    sprintf(msgbuf, "Invalid input argument l1a_file = %p, " 
        "granule_metadata = %p, number_of_scans = %d, frame_data = %p, " 
    	"mirror_side = %p", (void*)l1a_file, (void*)granule_metadata,
	number_of_scans, (void*)frame_data, (void*)mirr_side );
    modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf,
      "GEO_read_L1Ascan_metadata.c, GEO_read_L1Ascan_metadata"); 
    return FAIL;
  }
  else if(number_of_scans == 0)
    return SUCCESS; 


/*  
  Finish initializing the scan_metadata array with the contents
	defined above.  Data is retrieved into data types  
	conforming with the SDS's data type.
	
	See magic number warning in define section.
*/
	Scan_metadata[MIRR_SIDE_INDEX].data = (void *) mirr_side;
	for(i = 0; i < (int) NUM_METADATA; i++) {
		Scan_metadata[i].dims[0] = (long int) number_of_scans;
	}
	
/* read a scan_metadata[i].dims array into packet_data[i].data */

	for (i = 0; i < (int) NUM_METADATA; i++) {
	    if (getMODISarray(l1a_file, Scan_metadata[i].name, L1A_SCAN_META_GRP, 
	       start, Scan_metadata[i].dims, Scan_metadata[i].data) != MAPIOK) {
		sprintf(msgbuf, " getMODISarray reading %s", Scan_metadata[i].name);
		modsmf(MODIS_E_GEO, msgbuf,
		"GEO_read_L1Ascan_metadata.c, GEO_read_L1Ascan_metadata");
		retval=FAIL;
		}
	}

/* 
  load the output structure 
  FOR each of the number_of_scans scan's
*/
	for (scan = 0; scan < number_of_scans; scan++) {
		frame_data[scan].L1A_scan_id      = (int32)scan_id[scan];		
		frame_data[scan].EV_start         = EV_start[scan];
                for(i = 0; i < NUM_L1A_QUALITY_FLAGS; i++) {
                        frame_data[scan].L1A_scan_quality[i] = scan_quality[scan][i];
                }
              if((frame_count[scan][Earth_Frames]>=0) &&
                 (frame_count[scan][Earth_Frames] <=
		   granule_metadata->max_earth_frames))
	      {
		 frame_data[scan].EV_frames =
		   (int32)frame_count[scan][Earth_Frames];
              }
	      else
	      {
                 sprintf(msgbuf, "%d Earth frames in scan %d",
		   (int)frame_count[scan][Earth_Frames], scan);
                 modsmf(MODIS_W_GEO_FRAMENO_INPUT, msgbuf,
		   "GEO_read_L1Ascan_metadata.c, GEO_read_L1Ascan_metadata");

                 frame_data[scan].EV_frames = 0;
              }

	      frame_data[scan].SD_start = SD_start[scan];
              if((frame_count[scan][SD_Frames] >=0 ) &&
                 (frame_count[scan][SD_Frames] <=
		 granule_metadata->max_sd_frames))
	      {
		frame_data[scan].SD_frames =
		  (int32)frame_count[scan][SD_Frames];
	      }
              else
	      {
                sprintf(msgbuf, "%d SD Frames in scan %d",
		  (int)frame_count[scan][SD_Frames], scan);
                modsmf(MODIS_W_GEO_FRAMENO_INPUT, msgbuf,
		  "GEO_read_L1Ascan_metadata.c, GEO_read_L1Ascan_metadata");

                frame_data[scan].SD_frames = 0;
              }

	      frame_data[scan].SV_start	= SV_start[scan];
              if((frame_count[scan][SV_Frames] >= 0) &&
                 (frame_count[scan][SV_Frames] <=
		   granule_metadata->max_sv_frames))
	      {
                frame_data[scan].SV_frames =
		  (int32)frame_count[scan][SV_Frames];
	      }
              else
	      {
                sprintf(msgbuf, "%d SV frames in scan %d",
		  (int)frame_count[scan][SV_Frames], scan);
                modsmf(MODIS_W_GEO_FRAMENO_INPUT, msgbuf,
                  "GEO_read_L1Ascan_metadata.c, GEO_read_L1Ascan_metadata");

		frame_data[scan].SV_frames = 0;
              }

            frame_data[scan].SCI_STATE = sci_state[scan];
            frame_data[scan].SCI_ABNORM = sci_abnorm[scan];
            memcpy(frame_data[scan].Scan_type, scan_type[scan],sizeof(scan_type[scan]));
	}	
 
	return retval;
}
