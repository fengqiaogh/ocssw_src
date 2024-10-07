#include "PGS_MODIS_35251.h"
#include "GEO_global_arrays.h"
#include "GEO_output.h"
#include "GEO_product.h"
#include "GEO_validation.h"

#define DIMSDIM 3
#define MAXMSGBUFFLEN 200 
#define CB_VECTOR_NUMVALUES 9

PGSt_SMF_status GEO_write_scan_metadata(
	int				const scan_number,
        l1a_data_struct const		* const l1a_data,
        celestial_bodies_struct const	* const cb_vectors,
	MODFILE 			* const geo_file,
	PGSt_double const		rpy[],
	uint32				frame_quality[][MAX_FRAMES]
)
/*****************************************************************************
!C
!Description:   
  		Subroutine in output group of the Level-1A geolocation
                software to write one scan of scan-level metadata
  		to the geolocation product file using the MAPI.
 
!Input Parameters:
 		scan_number	the scan (array index) number
 		l1a_data	data read from L1A file
 		cb_vectors	solar and lunar position vectors.
 		geo_file	MAPI structure for geolocation file
		rpy		Roll, pitch, and yaw thermal corrections
		frame_quality	Ephemeris and attitude quality flags
  
!Output Parameters:
 		None
  
Return parameter:
 		SUCCESS if all metadata for the scan written to the 
  			geolocation file,
  		FAIL otherwise.
  
Externally defined:
		ATT_IDX			"GEO_global_arrays.h"
		ATTIT_ANG		"GEO_product.h"
		EPH_IDX			"GEO_global_arrays.h"
		EV_FRAMES		"GEO_product.h"
		EVTIME			"GEO_product.h"
		FILL_INT8		"GEO_output.h"
  		GEO_DOUBLE_FILLVALUE	"GEO_output.h"
		GEO_QUALITY		"GEO_product.h"
		IMPULSE_ENC		"GEO_product.h"
		IMPULSE_TIME		"GEO_product.h"
		L1_QUALITY		"GEO_product.h"
		MAX_FRAMES		"GEO_geo.h"
  		MODIS_E_BAD_INPUT_ARG	"PGS_MODIS_35251.h"
  		MODIS_E_GEO		"PGS_MODIS_35251.h"
                MODIS_N_GEO_MANEUVER    "PGS_MODIS_35251.h"
		MOON_VECTOR		"GEO_product.h"
		MSIDE			"GEO_product.h"
		NUM_IMPULSE		"GEO_product.h"
		NUM_L1A_QUALITY_FLAGS	"GEO_geo.h"
		ORB_POS			"GEO_product.h"
		ORB_VEL			"GEO_product.h"
                PGS_S_SUCCESS           "PGS_SMF.h"
		S_NUM			"GEO_product.h"
		S_TYPE			"GEO_product.h"
  		SCAN_TYPE_LEN		"GEO_parameters.h"	
  		SCAN_META_GRP		"L1A_data.h"
		SCTIME			"GEO_product.h"
		SD_FRAMES		"GEO_product.h"
		SDTIME			"GEO_product.h"
		SUN_AZIMUTH		"GEO_product.h"
		SUN_REF			"GEO_product.h"
		SUN_ZENITH		"GEO_product.h"
		SV_FRAMES		"GEO_product.h"
		SVTIME			"GEO_product.h"
		T_INST2ECR		"GEO_product.h"
		THERMCORR		"GEO_product.h"
  
Called by:
 		GEO_write_one_scan	"GEO_output.h"

Routines called:
		modsmf			"smfio.h"
		putMODISarray		"mapi.h"

!Revision History:
 * $Log: GEO_write_scan_metadata.c,v $
 * Revision 6.3  2011/02/14 22:05:45  kuyper
 * Removed casts that are no longer necessary with current version of M-API.
 *
 * Revision 6.2  2010/05/27 15:22:09  kuyper
 * Corrected dimension of frame_quality.
 *
 * Revision 6.1  2010/05/26 22:47:47  kuyper
 * Changed to return a status code.
 * Helped resolve Bug 1969 by adding rpy as input, and THERMCORR as output
 * Helped resolve Bug 2470 by dropping maneuver_list.
 * Helped resolve Bug 2472 by adding frame_quality as input, and ATT_QUALITY and
 *   EPH_QUALITY as output.
 * Added SDS name macros to the Externally Defined section of the Prolog.
 * Improved naming of loop counters.
 * Changed 'name' and 'data' members to be pointers to const,
 * Corrected error message to bring code and PDL in sync.
 *
 * Revision 5.3  2004/10/25 18:36:38  vlin
 * typo for GEO_in_maneuver fixed.  vlin@saicmodis.com
 *
 * Revision 5.2  2004/08/24 15:22:39  vlin
 * Input parameter maneuver_list added, call to GEO_in_maneuver updated.
 *
 * Revision 5.1  2004/08/13 15:16:55  vlin
 * Added C5 feature: write fourth quality flag as the manuever flag
 *                   based upon GEO_in_maneuver().
 *
 * Revision 4.2  2004/04/09 22:18:39  kuyper
 * Replaced FILL_BYTE with FILL_INT8.
 *
 * Revision 4.1  2003/02/21 23:10:32  kuyper
 * Corrected to use void* pointers with %p format code.

!Team-unique Header:

        This software is developed by the MODIS Science Data Support
        Team for the National Aeronautics and Space Administration,
        Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/
{
  /* Local variable declarations */
    float32 moon_vector[3];
    float32 rpy_angles[3];
    float32 sun_vector[3];
    float32 sun_azimuth, sun_zenith;
    float64 imp_enc[MAX_IMPULSE_NUMBER], imp_time[MAX_IMPULSE_NUMBER];
    int obj, dim; /* Loop index */
    PGSt_SMF_status retval = PGS_S_SUCCESS; 
    uint16 scan_id, EV_frames, SD_frames, SV_frames;
    char msgbuf[MAXMSGBUFFLEN];
    int8 geoqual[4];
    char filefunc[] = __FILE__ ", GEO_write_scan_metadata";

    struct {
	char const	*name; /* Name of the SDS object in the product file. */
	long		start[3];	/* Starting location for write. */
	long		dims[DIMSDIM];	/* dimensions of array to write.  */
	void const	*data;		/* data to write to the SDS.         */
    }scan_metadata[] = 
    { {S_NUM,        {0,0},	{1},				NULL},
      {EV_FRAMES,    {0,0},	{1},				NULL},
      {SD_FRAMES,    {0,0},	{1},				NULL},
      {SV_FRAMES,    {0,0},	{1},				NULL},
      {EVTIME,       {0,0},	{1},				NULL},
      {SDTIME,       {0,0},	{1},				NULL},
      {SVTIME,       {0,0},	{1},				NULL},
      {S_TYPE,       {0,0},	{1,SCAN_TYPE_LEN},		NULL},
      {MSIDE,        {0,0},	{1},				NULL},
      {L1_QUALITY,   {0,0},	{1,NUM_L1A_QUALITY_FLAGS},	NULL},
      {GEO_QUALITY,  {0,0},	{1,4},				NULL},
      {SCTIME,       {0,0},	{1},				NULL},
      {SUN_ZENITH,   {0,0},	{1},				NULL},
      {SUN_AZIMUTH,  {0,0},	{1},				NULL},
      {MOON_VECTOR,  {0,0},	{1,3},				NULL},
      {ORB_POS,      {0,0},	{1,3},				NULL},
      {ORB_VEL,      {0,0},	{1,3},				NULL},
      {T_INST2ECR,   {0,0,0},	{1,3,3},			NULL},
      {ATTIT_ANG,    {0,0},	{1,3},				NULL},
      {SUN_REF,      {0,0},	{1,3},  			NULL},
      {NUM_IMPULSE,  {0,0},	{1},				NULL},
      {IMPULSE_ENC,  {0,0},	{1,MAX_IMPULSE_NUMBER},		NULL},
      {IMPULSE_TIME, {0,0},	{1,MAX_IMPULSE_NUMBER},		NULL},
      {THERMCORR,   {0,0},	{1,3},				NULL},
      {ATT_QUALITY,  {0,0},	{1,MAX_FRAMES},			NULL},
      {EPH_QUALITY,  {0,0},	{1,MAX_FRAMES},			NULL},
    };

    if ((l1a_data == NULL) || scan_number < 0 ||
	scan_number >= l1a_data->num_scans || cb_vectors == NULL ||
	geo_file == NULL || rpy==NULL || frame_quality == NULL )
    {
       sprintf(msgbuf, "l1a_data = %p, scan_number = %i, cb_vectors = %p, "
	  "geo_file = %p, rpy = %p, frame_quality = %p",
	  (void*)l1a_data, scan_number, (void*)cb_vectors, (void*)geo_file,
	  (void*)rpy, (void*)frame_quality);
       modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);

       return MODIS_E_BAD_INPUT_ARG;
    }

    /* Finish initializing the scan_metadata array with the contents
       defined above.  Copy data into local variables, where
       necessary, to ensure conformity with the SDS's data type.         */

    obj = 0;
    /* uint16 */
    scan_metadata[obj].start[0] = (long)scan_number;
    scan_id = (uint16)l1a_data->frame_data[scan_number].L1A_scan_id;
    scan_metadata[obj].data = &scan_id;
    obj++;

    /* uint16 */
    scan_metadata[obj].start[0] = (long)scan_number;
    EV_frames = (uint16)l1a_data->frame_data[scan_number].EV_frames;
    scan_metadata[obj].data = &EV_frames;
    obj++;

    /* uint16 */
    scan_metadata[obj].start[0] = (long)scan_number;
    SD_frames = (uint16)l1a_data->frame_data[scan_number].SD_frames;
    scan_metadata[obj].data = &SD_frames;
    obj++;

    /* uint16 */
    scan_metadata[obj].start[0] = (long)scan_number;
    SV_frames = (uint16)l1a_data->frame_data[scan_number].SV_frames;
    scan_metadata[obj].data = &SV_frames;
    obj++;

    /* float64 */
    scan_metadata[obj].start[0] = (long)scan_number;
    scan_metadata[obj].data = &l1a_data->frame_data[scan_number].EV_start;
    obj++;

    /* float64 */
    scan_metadata[obj].start[0] = (long)scan_number;
    scan_metadata[obj].data = &l1a_data->frame_data[scan_number].SD_start;
    obj++;

    /* float64 */
    scan_metadata[obj].start[0] = (long)scan_number;
    scan_metadata[obj].data = &l1a_data->frame_data[scan_number].SV_start;
    obj++;

    /* char8 */
    scan_metadata[obj].start[0] = (long)scan_number;
    scan_metadata[obj].data = &l1a_data->frame_data[scan_number].Scan_type;
    obj++;

    /* uint16 */
    scan_metadata[obj].start[0] = (long)scan_number;
    scan_metadata[obj].data = &l1a_data->mirr_data[scan_number].mirr_side;
    obj++;

    /* int32 */
    scan_metadata[obj].start[0] = (long)scan_number;
    scan_metadata[obj].data =
	&l1a_data->frame_data[scan_number].L1A_scan_quality;
    obj++;

    /* int8 */
    geoqual[0] =  (int8)l1a_data->mirr_data[scan_number].impulse_flag;

    /* int8 */
    geoqual[1] = l1a_data->frame_data[scan_number].SCI_ABNORM;

    /* int8 */
    geoqual[2] = l1a_data->frame_data[scan_number].SCI_STATE;

    /* C5 feature: call subroutine GEO_in_maneuver() to check if the          */
    /*             spacecraft was in a listed maneuver at the specified time. */
    
    switch(GEO_in_maneuver(l1a_data->frame_data[scan_number].EV_start))
    {
    case PGS_S_SUCCESS:
	geoqual[3] = (int8)0;
	break;

    case MODIS_N_GEO_MANEUVER:
	geoqual[3] = (int8)1;
	break;

    default:
	geoqual[3] = FILL_INT8;
	break;
    }

    scan_metadata[obj].start[0] = (long)scan_number;
    scan_metadata[obj].start[1] =  0L;
    scan_metadata[obj].data = geoqual;
    obj++;

    if ( cb_vectors->sc.qa_flag == GOOD_DATA ) {

       /* float64  EV center time */
       scan_metadata[obj].start[0] = (long)scan_number;
       scan_metadata[obj].data = &cb_vectors->sc.time;
       obj++;

       /* float32  SD sun zenith */
       scan_metadata[obj].start[0] = (long)scan_number;
       sun_zenith = (float32)cb_vectors->SD_sun_zenith;
       scan_metadata[obj].data = &sun_zenith;
       obj++;

       /* float32  SD sun azimuth */
       scan_metadata[obj].start[0] = (long)scan_number;
       sun_azimuth = (float32)cb_vectors->SD_sun_azimuth;
       scan_metadata[obj].data = &sun_azimuth;
       obj++;

       /* float32  moon unit vector */
       scan_metadata[obj].start[0] = (long)scan_number;
       for(dim=0; dim<3; dim++)
	  moon_vector[dim] = (float32)cb_vectors->moon_unit_vector[dim];
       scan_metadata[obj].data = &moon_vector;
       obj++;

       /* float64 array[3] orb_pos  */
       scan_metadata[obj].start[0] = (long)scan_number;
       scan_metadata[obj].start[1] =  0L;
       scan_metadata[obj].data = cb_vectors->sc.positionECR;
       obj++;

       /* float64 array[3]  orb_vec */
       scan_metadata[obj].start[0] = (long)scan_number;
       scan_metadata[obj].start[1] =  0L;
       scan_metadata[obj].data = cb_vectors->sc.velocityECR;
       obj++;

       /* float64 array[3][3] T_inst2ecr */
       scan_metadata[obj].start[0] = (long)scan_number;
       scan_metadata[obj].start[1] =  0L;
       scan_metadata[obj].start[2] =  0L;
       scan_metadata[obj].data = cb_vectors->sc.T_inst2ecr;
       obj++;

       /* float64 array[3] attitude angles */
       scan_metadata[obj].start[0] = (long)scan_number;
       scan_metadata[obj].start[1] =  0L;
       scan_metadata[obj].data = cb_vectors->sc.eulerAngles;
       obj++;

       /* float32 array[3]  sun unit vector */
       scan_metadata[obj].start[0] = (long)scan_number;
       for(dim=0; dim<3; dim++)
	  sun_vector[dim] = (float32)cb_vectors->sun_unit_vector[dim];
       scan_metadata[obj].data = &sun_vector;
       obj++;

    } else  
	  /* obj must be incremented by number of members above. 
	   * This value will need to be changed if the members 
	   * of the scan_metadata struct changes.              */
	  obj += CB_VECTOR_NUMVALUES;

    /* uint8 */
    scan_metadata[obj].start[0] = (long)scan_number;
    scan_metadata[obj].data = &l1a_data->mirr_data[scan_number].num_impulse;
    obj++;

    if(l1a_data->mirr_data[scan_number].num_impulse != 0)
    {
	int imp;

	scan_metadata[obj].start[0] = (long)scan_number;
	scan_metadata[obj].dims[1] =
	  (long)l1a_data->mirr_data[scan_number].num_impulse;

	for (imp=0; imp < MAX_IMPULSE_NUMBER; imp++)
	    imp_enc[imp] =
		(float64)l1a_data->mirr_data[scan_number].mirr_impulse_enc[imp];
	scan_metadata[obj].data = &imp_enc;
	obj++;

	scan_metadata[obj].start[0] = (long)scan_number;
	scan_metadata[obj].dims[1] =
	  (long)l1a_data->mirr_data[scan_number].num_impulse;

	for (imp=0; imp<MAX_IMPULSE_NUMBER; imp++)
	    imp_time[imp] = (float64)l1a_data->mirr_data[scan_number].
		mirr_impulse_time[imp];
	scan_metadata[obj].data = &imp_time;
	obj++;
    }
    else
	obj += 2;
    
    scan_metadata[obj].start[0] = (long)scan_number;
    for(dim = 0; dim < 3; dim++)
	rpy_angles[dim] = (float32)(rpy[dim] *
	    (rpy[dim] <= THERMCORR_FVALUE ? 1.0 : RAD2DEG));
    scan_metadata[obj].data = &rpy_angles;
    obj++;

    scan_metadata[obj].start[0] = (long)scan_number;
    scan_metadata[obj].data = frame_quality + ATT_IDX;
    obj++;

    scan_metadata[obj].start[0] = (long)scan_number;
    scan_metadata[obj].data = frame_quality + EPH_IDX;
    obj++;

    for(obj = 0; obj < (int)(sizeof(scan_metadata)/sizeof(scan_metadata[0]));
	obj++)
    {
	if (scan_metadata[obj].data)
	{
	    if (putMODISarray(geo_file, scan_metadata[obj].name,
		SCAN_META_GRP, scan_metadata[obj].start,
		scan_metadata[obj].dims, scan_metadata[obj].data)
		!= MAPIOK)
	    {
		 sprintf(msgbuf, "putMODISarray(%.*s, %.*s)",
		     (MAXMSGBUFFLEN-12)/2, geo_file->filename,
		     (MAXMSGBUFFLEN-12)/2, scan_metadata[obj].name);
		 modsmf(MODIS_E_GEO, msgbuf, filefunc);

		 retval = MODIS_E_GEO;
	    }
	}
    }

    return retval;
}
