#include "imsl_wrap.h"
#include "PGS_MODIS_35251.h"
#include "GEO_global_arrays.h"
#include "GEO_inst.h"
#include "GEO_util.h"
#include "GEO_output.h"
#include "smfio.h"

int GEO_interp_mirr_enc(
	int const scan_number,
	int const sample_number,
	double * const sample_enc
	)
/*
!C*****************************************************************************
!Description:   
	Routine in the instrument model of the Level-1A geolocation software
	interpolate the mirror encoder number. Once per scan, it uses spline
	interpolation to estimate the mirror position for all frames in that
	scan. Each time it's called, it retrieves the interpolated value for
	the specified frame.

!Input Parameters:
        scan_number     The scan number
        sample_number   The sample number for which the encoder values

!Output Parameters:
        sample_enc      The pointer to sample encoder number

Return value:
        FAIL            If sample_enc is NULL, or the interpolation failed.
        SUCCESS         Otherwise.

Externally defined:
        Global Variables (input):
                num_impulse             The number of impulses in this scan.
                mirr_impulse_enc        The encoder values for each impulse
                mirr_impulse_time       The time values for each impulse
                sample_time             The sample time for each sample
                scan_start_time         The start time for each scan
        DBL_MAX                         <float.h>
        FAIL                            "GEO_basic.h"
        GEO_DOUBLE_FILLVALUE            "GEO_output.h"
        MAX_IMPULSE_NUMBER              "GEO_geo.h"
        MAX_SCAN_NUMBER                 "GEO_geo.h"
        MAX_PADDED                      "GEO_geo.h"
        MODIS_E_GEO                     "PGS_MODIS_35251.h"
        MODIS_E_BAD_INPUT_ARG           "PGS_MODIS_35251.h"
	MODIS_E_GEO_MIRR_MOTION		"PGS_MODIS_35251.h"
        MODIS_E_GEO_NO_ZERO_ENCODER     "PGS_MODIS_35251.h"
        SUCCESS                         "GEO_basic.h"

Called by:
        GEO_get_inst_mirr_normal

Routines Called:
        imsl_d_spline_interp            "imslwrap.h"
        imsl_d_spline_value             "imslwrap.h"
        imsl_d_machine                  "imslwrap.h"
        modsmf                          "smfio.h"

!Revision History:
   $Log: GEO_interp_mirr_enc.c,v $
   Revision 6.1  2009/05/11 17:49:26  xgeng
   Changed MAX_SCAN_SAMPLE to MAX_PADDED.
   Changed to determine padded_samples by examination of sample_time values,
   and to interpolate only for that many values.
   The above updates match with PDL revision 6.1 and 6.2.

   Revision 5.1  2006/10/05 16:10:31  kuyper
   Changed to call imsl_d_spline_value only once per scan, for all of the
     frames in the scan, rather than calling it seperately for each frame.

   Revision 4.6  2004/04/08 19:39:45  kuyper
   Corrected knot sequence to satisfy constraints.

   Revision 4.5  2004/03/17 18:35:20  kuyper
   Changed to fit encoder data to a linear trend plus a quadratic B-spline, with
     the knots of the B-spline chosen to force a smooth return to 0 at the ends
     of the knot sequence.

   Revision 4.4  2003/08/22 17:12:48  kuyper
   Corrected handling of starting point.

   Revision 4.3  2003/08/19 19:20:32  vlin
   Updated to avoid division by zero.

   Revision 4.2  2003/08/11 19:53:17  kuyper
   Corrections discovered during testing.

   Revision 4.1  2003/08/07 15:59:28  vlin
   Changed to interpolate intermediate values using splines, rather than fitting
    the data to chebyshev polynomials.
   Changed to use the interpolated spline to find the time where the encoder
    count is 0.0, rather than relying on GEO_prepare_mirr_data to extrapolate
    the first (unknown) encoder time from the later ones.

		4/5/95
		Ruiming Chen
		Finished coding.

Requirements:
                PR03-F-3.1.2-1

!Team-unique Header:

        This software is developed by the MODIS Science Data Support
        Team for the National Aeronautics and Space Administration,
        Goddard Space Flight Center, under contract NAS5-32373.

Design notes:
    For a linear least-squares fit to a set of N points (x,y), use the following
    formulae:

		    sum(x^2)*sum(y) - sum(x)*sum(x*y)
	intercept = ---------------------------------
			N*sum(x^2)-(sum(x))^2

		N*sum(x*y)-sum(x)*sum(y)
	slope = ------------------------
		N*sum(x^2)-(sum(x))^2

!END
*****************************************************************************/
{
#define QUADRATIC 3
#define DERIVS 3
    static double	last_scan_time = -DBL_MAX;
    static int scan_success=0;
    static double *enc_array;
    if(!enc_array) {
        enc_array = calloc(MAX_PADDED, sizeof(double));
    }
    static double intercept, slope;
    const  double aNaN = imsl_d_machine(6);
    char	 msgbuf[128]="";
    char   filefunc[] = __FILE__ ", GEO_interp_mirr_enc"; 


    if (sample_enc == NULL || scan_number < 0 || scan_number >= MAX_SCAN_NUMBER 
      || sample_number < 0 || sample_number >= MAX_PADDED ||
      num_impulse[scan_number] < 0 ||
      num_impulse[scan_number] > MAX_IMPULSE_NUMBER)
    {
	sprintf(msgbuf, "sample_enc = %p, scan = %d, frame = %d, impulse = %d",
		 (void *)sample_enc, scan_number, sample_number, 
		 num_impulse[scan_number]);
	modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
	return FAIL;
    }

    /* The following test uses the scan start time, rather than scan number, to
     * protect against the possibility that the last valid scan in one granule
     * might have the same scan number as the first one in the next granule.
     */
    if (fabs(scan_start_time[scan_number]-last_scan_time)>=DBL_MIN)
    {
	int		i, d, ndata, order;
	double		q, tzero, sqrt_arg;
	double		knots[MAX_IMPULSE_NUMBER+QUADRATIC], derivs[DERIVS];
	double		diff[MAX_IMPULSE_NUMBER-1];
	double		sumt=0.0, sumt2=0.0, sume=0.0, sumet=0.0;
	double		denom, delta_t;
        int             padded_samples = 0;
	Imsl_d_spline	*encoder_spline=NULL;

	last_scan_time = scan_start_time[scan_number];
	scan_success = 0;

	/* Don't use mirr_impulse_time[scan_number][0], it's bogus */
	ndata = num_impulse[scan_number] - 1;
	if (ndata < QUADRATIC)
	     order = ndata;
	else
	     order = QUADRATIC;

	/* Calculate coefficients for a least-squares linear fit. */
	for(i=1; i<num_impulse[scan_number]; i++)
	{
	    sumt += mirr_impulse_time[scan_number][i];
	    sumt2 += mirr_impulse_time[scan_number][i] *
		mirr_impulse_time[scan_number][i];
	    sume += mirr_impulse_enc[scan_number][i];
	    sumet += mirr_impulse_time[scan_number][i] *
		mirr_impulse_enc[scan_number][i];
	}

	denom = ndata*sumt2-sumt*sumt;
	intercept = sumt2*sume - sumt*sumet;
	slope = ndata*sumet-sumt*sume;
	if( fabs(denom) > fabs(intercept)/DBL_MAX &&
	     fabs(denom) > fabs(slope)/DBL_MAX) 
	{
	    intercept /= denom;
	    slope /= denom;
	}
	else
	{
	    /* Shouldn't happen for num_impulse > 1 */
	    sprintf(msgbuf, "%d", scan_number);
	    modsmf(MODIS_E_GEO_MIRR_MOTION, msgbuf, filefunc);
	    return FAIL;
	}

	/* Calculate differences from the least-squares linear fit. */
	for(i=1; i<num_impulse[scan_number]; i++)
	    diff[i-1] = mirr_impulse_enc[scan_number][i] -
		(intercept+slope*mirr_impulse_time[scan_number][i]);

	/* Set up knots that never double up. This means that the first order-2
	 * derivatives of the B-spline are forced to go continuously to 0 at the
	 * ends. That ensures that the combination of linear trend+B-spline
	 * returns as smoothly as possible to the linear trend, outside the
	 * range of the available data.
	 */

	delta_t =
	    mirr_impulse_time[scan_number][1]-mirr_impulse_time[scan_number][0];

	/* Knots evenly spaced from the second data point. The first data point
	 * is not used, because it's an arbitrary insert, rather than real data.
	 */
	for (i = 0; i < order; i++)
	    knots[i] = mirr_impulse_time[scan_number][1] + delta_t*(i+1-order);

	if(order&1) /* Knots are halfway between data points. */
	    for(i=order; i<ndata; i++)
		knots[i] = 0.5*(mirr_impulse_time[scan_number][i-order/2] +
		    mirr_impulse_time[scan_number][i-order/2+1]);
	else /* Knots coincide with data points. */
	    for(i=order; i<num_impulse[scan_number]; i++)
		knots[i] = mirr_impulse_time[scan_number][i-order/2];
	
	/* Knots evenly spaced from the last data point. */
	for(i=ndata; i<ndata+order; i++)
	    knots[i] = mirr_impulse_time[scan_number][ndata] +
		delta_t*(i-ndata);

	/* Calculate a spline interpolant of diff. */
	/* Don't use mirr_impulse_time[0], it's bogus */
	encoder_spline = imsl_d_spline_interp(ndata, 
	    mirr_impulse_time[scan_number]+1, diff,
	    IMSL_ORDER, order, IMSL_KNOTS, knots, 0);
	if (imsl_error_code()) {
	    sprintf(msgbuf, "imsl_d_spline_interp():%ld",
		(long)imsl_error_code());
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);
	    free(encoder_spline);
	    return FAIL;
	}
	 
	for (d = 0; d < DERIVS; d++)
	{
	    derivs[d] = imsl_d_spline_value(mirr_impulse_time[scan_number][0],
		encoder_spline, IMSL_DERIV, d, 0);
	    if ((aNaN != aNaN && derivs[d] != derivs[d]) ||
		(aNaN == aNaN && derivs[d] == aNaN) )
	    {
		sprintf(msgbuf, "imsl_d_spline_value() for derivs [%d],\n"
		       "imsl_error_code() = %ld", d, (long)imsl_error_code());
		modsmf(MODIS_E_GEO, msgbuf, filefunc);
		free(encoder_spline);
		return FAIL;
	    }
	}

	derivs[0] += intercept;
	derivs[1] += slope;

	/* Estimate tzero, the offset from mirr_impulse_time[scan_number][0] for
	 * which the interpolated curve is at the point * where scan start should
	 * occur:
	 *          derivs[0]+derivse[1]*tzero+0.5*derivs[2]*tzero*tzero
	 *                  = mirr_impulse_enc[scan_number][0]
	 */
	derivs[0] -= mirr_impulse_enc[scan_number][0];
	sqrt_arg = derivs[1]*derivs[1] - 2*derivs[0]*derivs[2];
	if (sqrt_arg < 0 || derivs[1] < 0 ) {
	    sprintf(msgbuf, "%d", scan_number);
	    modsmf(MODIS_E_GEO_NO_ZERO_ENCODER, msgbuf, filefunc); 
	    free(encoder_spline);
	    return FAIL;
	}

	/* See "Numerical Recipes in C" by Preuss et. al., section 5.5 */
	q = -(derivs[1]+sqrt(sqrt_arg))*0.5;

	/* Take the smaller of the two roots; the curve will normally be a very
	 * straight line, and 0.0 will be a very good estimate of the proper
	 * value, so the smaller root will normally be very much smaller than
	 * the larger one. Therefore don't calculate both values; the larger one
	 * might overflow!
	 * Don't test q*q vs fabs(derivs[0]*derivs[2]*0.5), because that might
	 * overflow (though this is pretty unlikely).
	 */

	/* In a typical implementation of C, DBL_MAX*DBL_EPSILON*DBL_EPSILON is
	 * much greater than 1.0. Under those circumstances, if the program
	 * manages to reach this point in the code, it is not possible for the
	 * tests involving DBL_MAX below to actually fail. However, proving
	 * that fact is complicated, so it's safer to leave these tests in.
	 */
	if (fabs(q) < sqrt(fabs(derivs[0]))*sqrt(fabs(derivs[2]*0.5)))
	{
	    if(fabs(derivs[2]) <= 2.0*fabs(q)/DBL_MAX)
	    {
		sprintf(msgbuf, "%d: derivs[2] too small", scan_number);
		modsmf(MODIS_E_GEO_NO_ZERO_ENCODER, msgbuf, filefunc); 
		free(encoder_spline);
		return FAIL;
	    }
	    else
		tzero = 2.0*q/derivs[2];
	}
	else
	{
	    if(fabs(q) <= fabs(derivs[0])/DBL_MAX)
	    {
		sprintf(msgbuf, "%d, q too small", scan_number);
		modsmf(MODIS_E_GEO_NO_ZERO_ENCODER, msgbuf, filefunc); 
		free(encoder_spline);
		return FAIL;
	    }
	    else
		tzero = derivs[0]/q;
	}

	/* Shift encoder times so imsl_d_spline_value(
	 * mirr_impulse_time[scan_number][0], encoder_spline) will be almost
	 * exactly mirr_impulse_enc[scan_number][0]. If DERIVS<=order-1,
	 * it is exact.
	 */
	for (i = 0; i < ndata+order; i++) 
	     encoder_spline->knots[0][i] -= tzero;
	for (i = 1; i<num_impulse[scan_number]; i++)
	    mirr_impulse_time[scan_number][i] -= tzero; 

        while (padded_samples < MAX_PADDED && sample_time[padded_samples] < GEO_DOUBLE_FILLVALUE) {
          padded_samples++;
        }

	(void)imsl_d_spline_value(0.0, encoder_spline,
	    IMSL_GRID_USER, padded_samples, sample_time, enc_array, 0);
	/* When using IMSL_GRID_USER, return value does not indicate error. */
	if(imsl_error_code())
	    modsmf(MODIS_E_GEO, "imsl_d_spline_value()", filefunc);
	else
	    scan_success = 1;

	free(encoder_spline);
    } /* end if scan_start_time != last_scan_time */

    if(!scan_success)
	return FAIL;
	/* First call for this scan failed; no point in further messages. */

    *sample_enc = enc_array[sample_number] +
	intercept + slope*sample_time[sample_number];
    return SUCCESS;
}

