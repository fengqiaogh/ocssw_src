/* 
 *Revision History:
 * $Log: GEO_hires.c,v $
 * Revision 6.5  2010/12/14 21:12:18  kuyper
 * Split samp into osamp and psamp, thereby simplifying the logic, enabling
 *   identification and correction of some off-by-one errors.
 *
 * Revision 6.4  2010/06/18 20:38:31  kuyper
 * Corrected to initialize 'triplets'.
 *
 * Revision 6.3  2009/06/12 16:01:46  kuyper
 * Corrected sign error in formula for the interpolation factor in the scan
 *   direction.
 * Corrected handling of the case where the first call to GEO_calculate_offsets
 *   fails.
 *
 * Revision 6.2  2009/05/22 18:15:20  kuyper
 * Removed design notes; the ones in the PDL are sufficient.
 * Added sample_flags to invalid input arguments message.
 * Corrected several misuses of N_samp.
 *
 * Revision 6.1  2009/04/28 18:05:37  kuyper
 * Initial revision.
 *
 * James Kuyper	James.R.Kuyper@nasa.gov
 */

#include "PGS_CSC.h"
#include "smfio.h"
#include "GEO_earth.h"
#include "GEO_geo.h"
#include "PGS_MODIS_35251.h"

#define KILO 1000.0

static double GEO_triple_product(
	double x[3],
	double y[3],
	double z[3]
)
/*!C****************************************************************************

!Description:   
	Computes the scalar triple product of three three-dimensional vectors.
	Geometrically, the scalar product gives the oriented volume of the
	parrallelopiped whose edges are made up of those vectors.

!Input Parameters:
	x, y, z:		Vectors to be multiplied.

!Output Parameters:
	None
	
Return Value:
	The triple product.

Externally Defined:
	None.

Called by:
	GEO_hires		"GEO_earth.h"

Routines Called:
	None

!Revision History:
See top of file.

James Kuyper	James.R.Kuyper@nasa.gov

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
    return x[0] * (y[1] * z[2] - y[2] * z[1]) +
	   x[1] * (y[2] * z[0] - y[0] * z[2]) +
	   x[2] * (y[0] * z[1] - y[1] * z[0]);
}

static int GEO_calculate_offsets(
	double		hires_scale,
	double		num,
	double		den,
	double		nst,
	double		nds,
	double		ndt,
	double		nqs,
	double		nqt,
	const double	delta[],
	const double	quad[],
	const double	normal[],
	const double	scan[],
	const double	track[],
	int8		sth[]
)

/*!C****************************************************************************

!Description:   
	This routine calculates the value of the scan, track, and height offsets
	for a particular value of num and den, such that the scan offset should
	be equal to num/den. The calculations are protected against possible
	over flows, and will be halted if the result at any step is out of
	range.

!Input Parameters:
	hires_scale	Scale factor for high resolution offsets.
	num		The numerator in the calculation of the scan offset
	den		The denominator in the calculation of the scan offset
	nst		Triple product of normal, scan, and track
	nds		Triple product of normal, delta, and scan
	ndt		Triple product of normal, delta, and track
	nqs		Triple product of normal, quad, and scan
	nqt		Triple product of normal, quad, and track
	delta		Difference between estimated and actual position
	quad		Vector describing deviation from parallelogram
	normal		Vector normal to the ellipsoid.
	scan		Vector in the scan direction.
	track		Vector in the track direction.

!Output Parameters:
	sth		The scan, track, and height offsets.

Return value:		Used to keep a count of successful calls to this routine
	0		If unsuccessful
	1		If successful

Externally Defined:
	HIRES_FVALUE	"GEO_geo.h"

Called by:
	GEO_hires	"GEO_earth.h"

Routines Called:

!Revision History:
See top of file.

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits
	N/A

Design Notes
	See prolog for GEO_hires() for derivation of the relevant equations.
!END**************************************************************************
*/
{
    int ind;		/* Loop counter */
    double doff[3]; 	/* temporary offsets */
    double den1, den2;	/* possible denominators */
    const double max_offset = (-0.5 - HIRES_FVALUE) * hires_scale;

    if(fabs(den) <= fabs(num)/DBL_MAX)
	return 0;

    doff[0] = num/den;
    if(fabs(doff[0]) >= max_offset)
	return 0;

    den1 = doff[0]*nqs - nst;
    den2 = doff[0]*nqt;
    if(fabs(den1) > fabs(den2))
    {
	num = nds;
	den = den1;
    }
    else
    {
	num = ndt - doff[0]*nst;
	den = den2;
    }

    if(fabs(den) <= fabs(num)/DBL_MAX)
	return 0;

    doff[1] = num/den;

    doff[2] = 0.0;
    for(ind=0; ind<3; ind++)
	doff[2] += normal[ind] * (delta[ind] - doff[0]*scan[ind] -
	    doff[1]*track[ind] - doff[0]*doff[1]*quad[ind]) / (KILO * KILO);

    if(fabs(doff[2]) >= max_offset)
	return 0;

    for(ind=0; ind<3; ind++)
	sth[ind] = (int)floor(doff[ind]/hires_scale + 0.5);

    return 1;
}


PGSt_SMF_status GEO_hires(
        uint16	N_samp,
	int	padded_samples,
	double	hires_scale,
        double	terrain_sample_position[][MAX_PADDED][3],
        double	ecr_sample_position[][MAX_PADDED][3],
        double	ecr_frame_position[][MAX_FRAMES][3],
	uint8   sample_flags[][MAX_PADDED],
	uint8   frame_flags[][MAX_FRAMES],
	int8	hires_offsets[][DETECTORS_QKM][SAMPLES_QKM]
)
/*!C****************************************************************************

!Description:   
	Routine for calculating high resolution offsets that can be applied to
	bilinear interpolation of the low resolution positions, to give the
	high resolution positions.

!Input Parameters:
        N_samp				Number of samples per frame for band
	padded_samples			Total number of samples calculated.
	hires_scale			Scale factor used for high resolution
					offsets.
        terrain_sample_position		High resolution positions in geodetic
					coordinates.
        ecr_sample_position		High resolution positions in ECR
					coordinates.
        ecr_frame_position		Low resolution positions in ECR
					coordinates.
	sample_flags			High resolution pixel flags
	frame_flags			Low resolution pixel flags

!Output Parameters:
	hires_offsets			High resolution scan, track, and height
					offsets.

Return values:
	MODIS_E_GEO_BAD_INPUT_ARG	If any pointer argument is null.
	PGS_S_SUCCESS			Otherwise

Externally Defined:
	DETECTORS_QKM			"GEO_geo.h"
	DETECTORS_1KM			"GEO_geo.h"
	INVALID_INPUT_DATA		"GEO_geo.h"
	HIRES_FVALUE			"GEO_geo.h"
	MAX_FRAMES			"GEO_geo.h"
	MODIS_E_GEO_BAD_INPUT_ARG	"PGS_MODIS_35251.h"
	MAX_PADDED			"GEO_geo.h"
	PGS_S_SUCCESS			"PGS_SMF.h"

Called by:
	GEO_aggregate			"GEO_earth.h"

Routines Called:
	GEO_calculate_offsets		"GEO_hires.c"
	GEO_triple_product		"GEO_hires.c"

!Revision History:
See top of file.

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits
	N/A

Design Notes
	See PDL.
!END**************************************************************************
*/
{
    /* loop counters to identify which high resolution pixel is being
     * processed.
     */
    int hdet;
    int samples = padded_samples - N_samp + 1;
    int frames = (samples+1)/N_samp;

    if(!terrain_sample_position || !ecr_sample_position || !ecr_frame_position
	|| !sample_flags || !frame_flags || !hires_offsets)
    {
	char filefunc[] = __FILE__ ", GEO_hires()";
	char msgbuf[PGS_SMF_MAX_MSGBUF_SIZE];

	sprintf(msgbuf, "terrain_sample_position=%p, "
	    "ecr_sample_position=%p, ecr_frame_position=%p, sample_flags=%p, "
	    "frame_flags=%p, hires_offsets=%p", (void*)terrain_sample_position,
	    (void*)ecr_sample_position, (void*) ecr_frame_position,
	    (void*)sample_flags, (void*)frame_flags, (void*)hires_offsets);
	modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);

	return MODIS_E_BAD_INPUT_ARG;
    }

    for(hdet = 0; hdet < DETECTORS_1KM * N_samp; hdet++)
    {
	int ind;  /* loop counter to index across coordinates. */
	/* The addition of N_samp cancels final -1, mathematically , but in C
	 * it is needed to ensure the same behavior. */
	int det = (hdet + N_samp - N_samp/2)/N_samp - 1;
	double dt;	/* Fractional offset in track direction. */
	int osamp;	/* output sample number	*/

	if(det < 0)
	    det = 0;
	else if(det > DETECTORS_1KM - 2)
	    det = DETECTORS_1KM - 2;

	dt = (hdet + 0.5)/N_samp - (det + 0.5);

	/* The bracketing detectors are det and det+1 */
	for(osamp =0; osamp < samples; osamp++)
	{	/* Determine which 1km frames to use for the 	*/
		/* estimate, and the corresponding values of ds.*/
	    double  estimated[3], delta[3], scan[3], track[3];
	    double normal[3], quad[3];
	    double nds, nqs, nst, ndt, nqt;  /* triple products */
	    double a, b, c, q;	/* quadratic coefficients. */
	    int8 sth[2][3];	/* tentative s/t/h triplets. */
	    int triplets=0;	/* number of triplets selected. */
	    double ds;		/* Fractional offset in scan direction */
	    /* sample number in padded hi-resolution arrays. */
	    int psamp = osamp + N_samp -1;
	    int frame = osamp/N_samp;

	    if(frame > frames - 2)
		frame = frames - 2;
	    ds = osamp/(double)N_samp - frame;
	    /* The bracketing frames are frame and frame+1 */

	    if(sample_flags[hdet ][psamp  ] >= INVALID_INPUT_DATA || 
		frame_flags[det  ][frame  ] >= INVALID_INPUT_DATA ||
		frame_flags[det  ][frame+1] >= INVALID_INPUT_DATA ||
		frame_flags[det+1][frame  ] >= INVALID_INPUT_DATA ||
		frame_flags[det+1][frame+1] >= INVALID_INPUT_DATA)
	    {
		for(ind = 0; ind<3; ind++) /* hires_offsets is NOT padded */
		    hires_offsets[ind][hdet][osamp] = HIRES_FVALUE;
	    }

	    /*Perform bilinear interpolation/extrapolation */

	    for(ind =0; ind<3; ind++)
	    {
		estimated[ind] =
		    (1 - dt)*((1 - ds)*ecr_frame_position[det  ][frame  ][ind] +
				   ds *ecr_frame_position[det  ][frame+1][ind])+
			 dt *((1 - ds)*ecr_frame_position[det+1][frame  ][ind] +
				   ds *ecr_frame_position[det+1][frame+1][ind]);
		delta[ind] =
		    ecr_sample_position[hdet][psamp][ind] - estimated[ind];
		scan[ind] = 
		    (1 - dt)*(ecr_frame_position[det  ][frame+1][ind]  -
			      ecr_frame_position[det  ][frame  ][ind]) +
			 dt *(ecr_frame_position[det+1][frame+1][ind]  -
			      ecr_frame_position[det+1][frame  ][ind]);
		track[ind] = 
		    (1 - ds)*(ecr_frame_position[det+1][frame  ][ind]  -
			      ecr_frame_position[det  ][frame  ][ind]) +
			 ds *(ecr_frame_position[det+1][frame+1][ind]  -
			      ecr_frame_position[det  ][frame+1][ind]);
		quad[ind] = ecr_frame_position[det  ][frame  ][ind]-
			    ecr_frame_position[det  ][frame+1][ind]-
			    ecr_frame_position[det+1][frame  ][ind]+
			    ecr_frame_position[det+1][frame+1][ind];
	    }

	    normal[0] = KILO *
		cos(terrain_sample_position[hdet][psamp][LAT_IDX]) *
		cos(terrain_sample_position[hdet][psamp][LON_IDX]);
	    normal[1] = KILO *
		cos(terrain_sample_position[hdet][psamp][LAT_IDX]) *
		sin(terrain_sample_position[hdet][psamp][LON_IDX]);
	    normal[2] = KILO *
		sin(terrain_sample_position[hdet][psamp][LAT_IDX]);

	    nds = GEO_triple_product(normal, delta, scan);
	    ndt = GEO_triple_product(normal, delta, track);
	    nst = GEO_triple_product(normal, scan, track);
	    nqs = GEO_triple_product(normal, quad, scan);
	    nqt = GEO_triple_product(normal, quad, track);

	    /* Calculate quadratic coefficients */
	    a = nst * nqs;
	    b = nds * nqt - nst * nst - ndt * nqs;
	    c = ndt * nst;

	    if(b*b >= 4.0*a*c)
	    {
		q = -0.5*(b + (b<0 ? -1 : 1)*sqrt(b*b - 4.0*a*c));

		triplets = GEO_calculate_offsets(hires_scale, q, a, nst, nds,
		    ndt, nqs, nqt, delta, quad, normal, scan, track, sth[0]);
		triplets += GEO_calculate_offsets(hires_scale, c, q, nst, nds,
		    ndt, nqs, nqt, delta, quad, normal, scan, track,
		    sth[triplets]);
	    }

	    if(triplets > 1 && 
		sth[1][0]*sth[1][0]+sth[1][1]*sth[1][1]+sth[1][2]*sth[1][2] <
		sth[0][0]*sth[0][0]+sth[0][1]*sth[0][1]+sth[0][2]*sth[0][2]) 
		    memcpy(sth[0], sth[1], sizeof sth[0]);

	    if(triplets > 0)
	    {
		for(ind=0; ind < 3; ind++) /* hires_offsets is NOT padded */
		    hires_offsets[ind][hdet][osamp] = sth[0][ind];
	    }
	    else
	    {
		for(ind=0; ind < 3; ind++) /* hires_offsets is NOT padded */
		    hires_offsets[ind][hdet][osamp] = HIRES_FVALUE;
	    }
	}
    }

    return PGS_S_SUCCESS;
}

