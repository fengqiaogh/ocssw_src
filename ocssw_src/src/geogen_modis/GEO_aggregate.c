#include "GEO_geo.h"
#include "GEO_earth.h"
#include "GEO_output.h"
#include "GEO_product.h"
#include "GEO_global_arrays.h"
#include "PGS_MODIS_35251.h"
#include "PGS_SMF.h"
#include "PGS_CSC.h"
#include "smfio.h"

/* Revision History:
$Log: GEO_aggregate.c,v $
Revision 6.8  2011/02/21 19:57:16  kuyper
In order to resolve feature-request Bug 3446, changed to aggregate landsea
  mask, and to create new waterpresent array.

Revision 6.7  2010/12/14 22:01:14  kuyper
Moved responsibility for conversion from height above ellipsoid to height
  above geoid into this module.

Revision 6.6  2010/06/28 19:17:05  kuyper
Corrected QFL_IDXS loop conditions.

Revision 6.5  2010/05/27 17:50:26  kuyper
Corrected dimensions of frame_quality.
Removed references to quality flag values that could not be present.

Revision 6.4  2010/04/02 21:35:41  kuyper
Helped resolve Bug 2472 by aggregating ephemeris/attitude quality flags.
Dropped unused parameter.

James Kuyper	james.kuyper@sigmaspace.com

Revision 6.3  2009/05/31 23:08:54  kuyper
Corrected to use DETECTORS_1KM as loop limit.
Corrected order of arguments to PGS_CSC_GEOtoECR().

Revision 6.2  2009/05/22 17:01:25  kuyper
Corrected some misuses of N_samp.
Corrected error messages.

Revision 6.1  2009/05/19 16:59:16  kuyper
Initial revision.

James Kuyper	James.R.Kuyper@nasa.gov
*/

static int GEO_aggregate_lw(
	int const	lwmask_count[],
	uint8		*water
)
/*!C****************************************************************************

!Description:   
	Given a count of the number of pixels by Land/Water Mask value,
	computes an aggregate Land/Water Mask value to describe the entire set
	of pixels. The total of the water pixels is written to the location
	pointed at by the "water" paramet

!Input Parameters:
	lwmask_count		The counts.

!Output Parameters:
	water			A pointer to the place where the water total
				should be written.

Return Value:
	The Land/Water mask value describing the set of pixels.
	CONTINENTAL	Ocean > 5km from coast AND > 50m deep AND < 500m deep.
	DEEP_INLAND	Inland water > 5km from shoreline AND > 50m deep
	DEEP_OCEAN	Ocean > 500m deep
	DRYLAND		not anything else
	EPHEMERAL	Ephemeral (intermittent) Water.
	SHALLOW_INLAND	Inland Water < 5km from shore OR < 50m deep
	SHALLOW_OCEAN	Ocean <5k from coast OR <50m deep

Externally Defined:
	COAST		"GEO_geo.h"
	CONTINENTAL	"GEO_geo.h"
	DEEP_INLAND	"GEO_geo.h"
	DEEP_OCEAN	"GEO_geo.h"
	DRYLAND		"GEO_geo.h"
	EPHEMERAL	"GEO_geo.h"
	NUM_LWMASK	"GEO_geo.h"
	SHALLOW_INLAND	"GEO_geo.h"
	SHALLOW_OCEAN	"GEO_geo.h"

Called by:
	GEO_aggregate()

Routines Called:
	None

!Revision History:
	See top of file.

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits

Design Notes
!END**************************************************************************
*/
{
    /* Treat Coastal pixels as land */
    int land = lwmask_count[DRYLAND] + lwmask_count[COAST];
    int inland = lwmask_count[SHALLOW_INLAND] + lwmask_count[DEEP_INLAND];
    int ocean = lwmask_count[SHALLOW_OCEAN] + lwmask_count[CONTINENTAL] +
	lwmask_count[DEEP_OCEAN];
    int moderate = lwmask_count[DEEP_INLAND] + lwmask_count[CONTINENTAL];
    int lwmask_sum;
    *water = inland + ocean + lwmask_count[EPHEMERAL];
    lwmask_sum = land + *water;

    if(lwmask_count[EPHEMERAL] > lwmask_sum / 2)
	 return EPHEMERAL;

    if(land + lwmask_count[EPHEMERAL] > lwmask_sum / 2)
	return DRYLAND;	/* GEO_aggregate will convert to COAST if *water>0. */
    else
    {
	int shallow = lwmask_count[EPHEMERAL] + lwmask_count[SHALLOW_INLAND]
	    + lwmask_count[SHALLOW_OCEAN];

	if(shallow > *water/2)
	{
	    if(inland > *water/2)
		return SHALLOW_INLAND;
	    else
		return SHALLOW_OCEAN;
	}
	if(inland > *water/2)
	    return DEEP_INLAND;

	/* Redefine "shallow" downwards: */
	shallow += moderate;
	if(shallow > *water/2)
	    return CONTINENTAL;
	else
	    return DEEP_OCEAN;
    }
}

PGSt_SMF_status GEO_aggregate(
	int32		EV_frames,
	uint16		N_samp,
	double		hires_scale,
	unsigned char	sample_flags[][MAX_PADDED],
	double  	ecr_sample_position[][MAX_PADDED][3],
	double		ecr_sc_sample_position[][3],
	double		terrain_sample_position[][MAX_PADDED][3],
	uint32		sample_quality[][QFL_IDXS],
	uint8		sample_landsea[][MAX_PADDED],
	double		ecr_frame_position[][MAX_FRAMES][3],
	double		terrain_frame_position[][MAX_FRAMES][3],
	uint8		frame_flags[][MAX_FRAMES],
	double		ecr_sc_frame_position[MAX_FRAMES][3],
	int8		hires_offsets[][DETECTORS_QKM][SAMPLES_QKM],
	uint32		frame_quality[][MAX_FRAMES],
	uint8		frame_landsea[][MAX_FRAMES],
	uint8		frame_waterpresent[][MAX_FRAMES]
)
/*!C****************************************************************************

!Description:   

!Input Parameters:
	EV_frames		Number of frames of data to process.
	N_samp			Samples per frame for high-resolution band.
	hires_scale		Scale factor for high resolution offsets (km)
	sample_flags		High resolution quality/error flags
	ecr_sample_position	High resolution ground position, ECR coordinates
	ecr_sc_sample_position	High resolution spacecraft position, ECR coords
	terrain_sample_position	High resolution ground position, geodetic ""
	sample_quality		High resolution ephemeris/attitude quality flags

!Output Parameters:
	ecr_frame_position	Low resolution ground position, ECR coordinates
	terrain_frame_position	Low resolution ground position, geodetic ""
	frame_flags		Low resolution quality/error flags
	ecr_sc_frame_position	Low resolution spacecraft position, ECR coords
	hires_offsets		High resolution offsets.
	frame_quality		Low resolution ephemeris/attitude quality flags

Return Values:
	MODIS_E_BAD_INPUT_ARG	If any pointer parameter is NULL
	MODIS_E_GEO		If GEO_hires() fails.
	PGS_S_SUCCESS		Otherwise

Externally Defined:
	ATT_IDX			"GEO_global_arrays.h"
	COAST			"GEO_geo.h"
	DETECTORS_QKM		"GEO_geo.h"
	DRYLAND			"GEO_geo.h"
	EARTH_MODEL		"GEO_earth.h"
	EPH_IDX			"GEO_global_arrays.h"
	EPH_LONG_FOLLOW		"GEO_product.h"
	EPH_LONG_PRECEED	"GEO_product.h"
	EPH_SHORT_FOLLOW	"GEO_product.h"
	PEH_SHORT_PRECEED	"GEO_product.h"
	GEO_DOUBLE_FILLVALUE	"GEO_output.h"
	HT_IDX			"GEO_global_arrays.h"
	MAX_FRAMES		"GEO_geo.h"
	MODIS_E_BAD_INPUT_ARG	"PGS_MODIS_35251.h"
	MODIS_E_GEO		"GEO_geo.h"
	NO_ELLIPSE_INTERSECT	"GEO_geo.h"
	MAX_PADDED		"GEO_geo.h"
	PGS_S_SUCCESS		"PGS_SMF.h"
	SAMPLES_QKM		"GEO_geo.h"

Called by:
	GEO_locate_one_scan		"GEO_output.h"

Routines Called:
	GEO_aggregate_lw()		"GEO_aggregate.c"
	GEO_get_geoid()			"GEO_earth.h"
	GEO_hires			"GEO_earth.h"
	modsmf				"smfio.h"
	PGS_CSC_ECRtoGEO		"PGS_CSC.h"
	PGS_CSC_GEOtoECR		"PGS_CSC.h"

!Revision History:
	See top of file.

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits
	N/A

Design Notes
!END**************************************************************************
*/
{
    int frame, efat; 	/* Loop counters */
    char filefunc[] = __FILE__ ", GEO_aggregate";
    char msgbuf[PGS_SMF_MAX_MSGBUF_SIZE];

    if(!sample_flags || !ecr_sample_position || !ecr_sc_sample_position ||
	!terrain_sample_position || !ecr_frame_position ||
	!terrain_frame_position || !frame_flags || !ecr_sc_frame_position ||
	!hires_offsets)
    {

	sprintf(msgbuf, "sample_flags:%p ecr_sample_position:%p\n"
	    "ecr_sc_sample_position:%p terrain_sample_position:%p\n"
	    "sample_quality:%p ecr_frame_position:%p\n"
	    "terrain_frame_position:%p frame_flags:%p\n"
	    "ecr_sc_frame_position:%p hires_offsets:%p frame_quality:%p",
	    (void*)sample_flags, (void*)ecr_sample_position,
	    (void*)ecr_sc_sample_position, (void*)terrain_sample_position,
	    (void*)sample_quality, (void*)ecr_frame_position,
	    (void*)terrain_frame_position, (void*)frame_flags,
	    (void*)ecr_sc_frame_position, (void*)hires_offsets,
	    (void*)frame_quality);
	modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);

	return MODIS_E_BAD_INPUT_ARG;
    }

    for(frame = 0; frame < EV_frames; frame++)
    {
	int det, off;
	/* Subsample the spacecraft positions. */
	memcpy(ecr_sc_frame_position[frame],
	    ecr_sc_sample_position[N_samp * (frame + 1) - 1],
	    sizeof *ecr_sc_frame_position);

	for(efat=EPH_IDX; efat<QFL_IDXS; efat++)
	{
	    frame_quality[efat][frame] = sample_quality[frame*N_samp][efat];

	    for(off=1; off < 2*N_samp-1; off++)
		frame_quality[efat][frame] |=
		    sample_quality[frame*N_samp+off][efat];
	}
	for(det = 0; det < DETECTORS_1KM; det++)
	{
	    int idet, ind;
	    double latitude, longitude, altitude, height=0.0;
	    PGSt_double posecr[3];	/* position in ECR coordinates. */
	    int lwmask_count[NUM_LWMASK]={0};

	    for(ind = 0; ind < 3; ind++)
		ecr_frame_position[det][frame][ind] = 0.0;
	    frame_flags[det][frame] = 0;

	    if(N_samp == 1)
	    { /* No aggregation is needed */
		frame_flags[det][frame] = sample_flags[det][frame];
		memcpy( ecr_frame_position[det][frame],
		    ecr_sample_position[det][frame],
		    sizeof ecr_frame_position[det][frame]);
		memcpy( terrain_frame_position[det][frame],
		    terrain_sample_position[det][frame],
		    sizeof terrain_frame_position[det][frame]);
		continue;
	    }

	    for(idet = N_samp*det; idet < N_samp*(det+1); idet++)
	    {
		for(off = 0; off < 2*N_samp-1; off++)
		{
		    /* Weighting factors along the scan direction for N_samp
		     * from 1 to 4. The factors that would apply for N_samp==3
		     * are included, even though they won't be used, to
		     * simplify use of this table.
		     */
		    static const double weight[4][7]= {
			{1.0},
			{1.0/8.0,  2.0/8.0,  1.0/8.0},
			{1.0/27.0, 2.0/27.0, 3.0/27.0, 2.0/27.0, 1.0/27.0},
			{1.0/64.0, 2.0/64.0, 3.0/64.0, 4.0/64.0, 3.0/64.0,
			    2.0/64.0, 1.0/64.0}
		    };

		    /* Aggregate the frame flags */
		    frame_flags[det][frame] |=
			sample_flags[idet][frame*N_samp+off];

		    /* Count the landsea_mask values. */
		    if(sample_landsea[idet][frame*N_samp+off] < NUM_LWMASK)
			lwmask_count[sample_landsea[idet][frame*N_samp+off]]
			    += (int)(N_samp*N_samp*N_samp*weight[N_samp-1][off]);

		    /* Calculate the weighted average position */
		    for(ind = 0; ind < 3; ind++)
			ecr_frame_position[det][frame][ind] +=
			    weight[N_samp-1][off] *
			    ecr_sample_position[idet][frame*N_samp+off][ind];

		    /* Calculate the weighted average height above the
		     * ellipsoid */
		    height += weight[N_samp-1][off]*
			terrain_sample_position[idet][frame*N_samp+off][HT_IDX];
		}
	    } 

	    /* Determine aggregated landsea_mask. */
	    frame_landsea[det][frame] = GEO_aggregate_lw(lwmask_count,
		frame_waterpresent[det]+frame);

	    if(frame_flags[det][frame] >= NO_ELLIPSE_INTERSECT)
	    {
		for(ind = 0; ind < 3; ind++)
		    ecr_frame_position[det][frame][ind] = GEO_DOUBLE_FILLVALUE;
		continue;
	    }

	    for(ind = 0; ind < 3; ind++)
		posecr[ind] = ecr_frame_position[det][frame][ind];

	    if(PGS_CSC_ECRtoGEO(posecr, EARTH_MODEL, &longitude, &latitude,
		&altitude) != PGS_S_SUCCESS)
	    {
		frame_flags[det][frame] = INVALID_INPUT_DATA;
		continue;
	    }

	    terrain_frame_position[det][frame][LAT_IDX] = latitude;
	    terrain_frame_position[det][frame][LON_IDX] = longitude;

	    /* Note that the altitude of the weighted average ECR position will
	     * be systematically biased on the order of 0.1 m lower than the
	     * average height, due to the curvature of the earth. Therefore we
	     * use the average height.
	     */
	    terrain_frame_position[det][frame][HT_IDX] = height;

	    /* Which means we need to correct the ECR coordinates: */
	    if(PGS_CSC_GEOtoECR(longitude, latitude, height, EARTH_MODEL,
		posecr) != PGS_S_SUCCESS)
	    {
		frame_flags[det][frame] = INVALID_INPUT_DATA;
		continue;
	    }

	    for(ind = 0; ind < 3; ind++)
		ecr_frame_position[det][frame][ind] = posecr[ind];
 
	    /* Convert to height above geoid. */
	    terrain_frame_position[det][frame][HT_IDX]
		-= GEO_get_geoid(latitude, longitude);
	    /* Note: can only get here if GEO_read_DEM() has already
	     * successfully called GEO_get_geoid(). Therefore, this call
	     * is guaranteed to not return BAD_GEOID.
	     */
	}
    }

    /* Set coastlines for frame_landsea. */
    for(frame = 0; frame < EV_frames; frame++)
    {
	int det;

	for(det = 0; det < DETECTORS_1KM; det++)
	{
	    if(frame_landsea[det][frame] == DRYLAND)
	    {
		int min_det = det - 1;
		int max_det = det + 2;
		int min_frame = frame - 1;
		int max_frame = frame + 2;
		int d;

		if(min_det < 0)
		    min_det = 0;
		else if(max_det > DETECTORS_1KM)
		    max_det = DETECTORS_1KM;
		if(min_frame < 0)
		    min_frame = 0;
		else if(max_frame > EV_frames)
		    max_frame = EV_frames;

		for(d=min_det; d<max_det; d++)
		{
		    int f;

		    for(f=min_frame; f<max_frame; f++)
		    {
			if(frame_waterpresent[d][f])
			{
			    frame_landsea[det][frame] = COAST;
			    break;
			}
		    }
		    if(f<max_frame)
			break;	/* the pixel is coastal. */
		}
	    }
	}
    }


    /* Fix up frame flags for internal consistency */
    for(frame=0; frame<EV_frames; frame++)
    {
	for(efat=EPH_IDX; efat<QFL_IDXS; efat++)
	{
	    if(frame_quality[efat][frame] & (uint32)PGSd_NO_DATA)
		frame_quality[efat][frame] = (uint32)PGSd_NO_DATA;

	    if(frame_quality[efat][frame] & EPH_LONG_PRECEED)
		frame_quality[efat][frame] &= ~EPH_SHORT_PRECEED;

	    if(frame_quality[efat][frame] & EPH_LONG_FOLLOW)
		frame_quality[efat][frame] &= ~EPH_SHORT_FOLLOW;
	}
    }

    if(N_samp > 1 && GEO_hires(N_samp, N_samp * EV_frames + N_samp - 1,
	hires_scale, terrain_sample_position, ecr_sample_position,
	ecr_frame_position, sample_flags, frame_flags, hires_offsets)
	!= PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "GEO_hires(%d, %d, %f)", (int)N_samp,
	    (int)(N_samp * EV_frames + N_samp - 1), hires_scale);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);

	return MODIS_E_GEO;
    }

    return PGS_S_SUCCESS;
}
