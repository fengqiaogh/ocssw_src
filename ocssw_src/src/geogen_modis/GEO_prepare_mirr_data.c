#include <stdio.h>
#include "PGS_MODIS_35251.h"
#include "GEO_global_arrays.h"
#include "GEO_input.h"
#include "GEO_validation.h"
#include "GEO_earth.h"
#include "smfio.h"

static int compare_double(
	const void *first,
	const void *second
){	/* Designed for use with qsort() or bsearch() */
    const double *left = (double*)first;
    const double *right = (double*)second;
    if(*left<*right)
	return -1;
    if(*left>*right)
	return 1;
    return 0;
}

PGSt_SMF_status GEO_prepare_mirr_data(
        uint16		earth_encoder_times[MAX_SCAN_NUMBER][ENCODER_LENGTH],
        int16		view_sector_start[MAX_SCAN_NUMBER][SECTOR_LENGTH],
        const mirror_preparation_struct *mirr_params,
        double const                    t_frame,
        l1a_data_struct			*l1a_data
)
/*****************************************************************************
!C

!Description:
		Routine in Input group of the Level-1A geolocation
                software to unpack the mirror encoder data from the
		encoder and sector start segments.  For each scan, it
		unpacks mirror encoder times from the encoder data,
		computes encoder values from the sector start and mirror
		side, and validates the data

!Input Parameters:
        earth_encoder_times     Array of Earth View Encoder Time Data (time
                                samples at every 100th mirror encoder pulse
                                over the earth view)
        view_sector_start       array of View Sector Start Encoder Count and
                                Vernier Count words
        mirr_params             parameters describing the mirror operation.
        t_frame                 The time interval between frames.

!Input/Output Parameters:
	l1a_data		Data from the L1A file.

!Output Parameters:
	None

Return Values:
        MODIS_E_BAD_INPUT_ARG   if input parameter check fails
        PGS_S_SUCCESS           otherwise
 
Externally Defined:
        BAD_DATA                "GEO_validation.h"
        CHAN_A                  "GEO_parameters.h"
        ENCODER_LENGTH          "GEO_geo.h"
        GOOD_DATA               "GEO_validation.h"
        MAX_IMPULSE_NUMBER      "GEO_geo.h"
        MAX_SCAN_NUMBER         "GEO_geo.h"
        MODIS_E_BAD_INPUT_ARG   "PGS_MODIS_35251.h"
        MODIS_E_DATA_SCAN       "PGS_MODIS_35251.h"
        SECTOR_LENGTH           "GEO_geo.h"
        SUCCESS                 "GEO_basic.h"
        PGS_S_SUCCESS           "PGS_SMF.h"

Called by: 	GEO_prepare_l1a_data

Routines called:
        GEO_del_limit_check             "GEO_validation.h"
        GEO_abs_limit_check             "GEO_validation.h"
        GEO_get_ancil_packet_time       "GEO_earth.h"
        modsmf                          "smfio.h"

!Revision History:
   $Log: GEO_prepare_mirr_data.c,v $
   Revision 6.2  2010/06/18 20:41:04  kuyper
   Corrected parameters that were pointers to array by removing their leading
     'const' qualifiers.

   Revision 6.1  2010/04/08 18:49:59  kuyper
   Corrected to have encoder_adjustment default to 0.
   Replaced arguments derived from an l1a_data_struct with a pointer to the
     l1a_data_struct that they came from.
   Partially addressed Bug 17 by validating earth_encoder_times and EV_start
     times using l1a_data->fill_values.
   Corrected buffer underrun error.

   Revision 4.15  2004/06/14 22:21:41  kuyper
   Corrected to not adjust bogus first time interval.

   Revision 4.14  2004/06/08 20:09:24  vlin
   added if condition around line 263

   Revision 4.13  2004/06/04 23:35:46  kuyper
   Corrected to not include bogus sample in statistics.
   Corrected to re-initialize time_count for each scan during second pass.

   Revision 4.12  2004/06/01 22:17:28  kuyper
   Corrected way that gap correction was applied.

   Revision 4.11  2004/05/21 22:59:26  kuyper
   Corrected handling of missing mirror encoder times.

   Revision 4.9  2004/04/09 21:59:28  kuyper
   Corrected handling of low-count cases.

   Revision 4.8  2004/04/09 17:50:43  vlin
   replaced PACKET_PERIOD with mirr_params-> packet_interval.

   Revision 4.7  2004/03/30 22:56:33  kuyper
   Corrected gap in flow of encoder values.

   Revision 4.6  2004/03/19 22:48:59  vlin
   Added fill values check, changed to calibrate jump detection parameters
   using current granule.  Validate twice instead of once:  once before calibration,
   and second time after applying jump correction.

   Revision 4.5  2003/08/11 19:10:16  kuyper
   Changed handling of first impulse time and encoder; needed by
     GEO_interp_mirr_enc(), even if following time count is bad.

   Revision 4.4  2003/08/07 19:08:19  kuyper
   Added GEO_earth.h to the #include list.

   Revision 4.3  2003/06/03 18:49:28  vlin
   Corrected limits on final impulse loop according to PDL change

   Revision 4.2  2003/05/28 15:59:41  vlin
   Changed to partially use and otherwise skip over bad data, rather than
   replacing it with interpolated or extrapolated values.
   Adjust encoder times for the difference in clock rates of the
   encoder timer and the spacecraft clock.
   Changed to return SMF status codes.
   vlin@saicmodis.com

   Revision 4.1  2003/02/21 22:23:47  kuyper
   Corrected pointer types used with %p.


Requirements:
                PR03-F-2.5-1
                PR03-F-2.5-2
                PR03-F-2.5-3
                PR03-F-2.5-4

!Team-unique Header:

                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

Design Notes:
		Assumes that sector_word[i] and vernier_word[i] < SECTOR_LENGTH 
                for i = 0 and i = 1.

!END
**************************************************************************/

{
    const double TIME_ROLL_CORRECTION = 64000.0; 
				/* encoder time rollover correction */
    const int SIDEB_OFF = 2;	/* the offset in words of the encoder data */
    int vernier = 0;		/* vernier offset count */
    int offset, scan, sample;
    int num_impulse;		/* count of good encoder values */
    size_t best_low = 0, best_high = 0;
    size_t num_times = 0, nnorm = 0, njump = 0;
    size_t low, high;
    double best = DBL_MAX;
    double cum_sum[MAX_SCAN_NUMBER*MAX_IMPULSE_NUMBER];
    double cum_sum2[MAX_SCAN_NUMBER*MAX_IMPULSE_NUMBER];
    double sorted_times[MAX_SCAN_NUMBER*MAX_IMPULSE_NUMBER][2];
    double norm_sum, jump_sum, norm_sum2;
    double jump_sum2, criterion, phase;
    int encoder_flags[MAX_IMPULSE_NUMBER] =  {0};
			/* validation flag for each encoder time delta	*/
    double encoder_adjustment=0; /* The adjustment appropriate for mirr_chan. */
    double factor=1.0;	/* correction needed for the difference in speed */
			/* of the encoder timer and the spacecraft clock */
    double gap_start = 0.0, gap_end = 0.0;
			/* range of impulse times that need adjustment	*/
    double impulse_enc[MAX_IMPULSE_NUMBER]={0.0}; 
			/* the Earth sector start encoder count state	*/
    double impulse_time[MAX_IMPULSE_NUMBER]={0.0};
			/* the current impulse time			*/
    double time_count[MAX_IMPULSE_NUMBER] = {0.0}; /* encoder time deltas */
			/* between every  100 encoder counts		*/
    double mirr_del_limit, mirr_abs_limit[2];
    char filefunc[] = __FILE__ ", GEO_prepare_mirr_data";
    char msgbuf[256] = "";   


    if ((earth_encoder_times == NULL) ||  (view_sector_start == NULL) || 
	(mirr_params == NULL) || (l1a_data == NULL))
    {
       sprintf(msgbuf, "earth_encoder_times:%p, view_sector_start:%p,\n"
	   "mirr_params:%p, l1a_data:%p", (void*)earth_encoder_times,
	   (void*)view_sector_start, (void*)mirr_params, (void*)l1a_data);
       modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
       return MODIS_E_BAD_INPUT_ARG;
    }
    if (l1a_data->num_scans > MAX_SCAN_NUMBER ||  l1a_data->num_scans < 0)
    {
	sprintf(msgbuf, "num_scans:%d", l1a_data->num_scans);
	modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
	return MODIS_E_BAD_INPUT_ARG;
    }

    /* The first encoder time reported after the creation of a new engineering
     * data packet usually shows a jump by a very consistent amount of time.
     * the interval of time between packets is extremely reliable, so the
     * affected encoder times can be reliably identified by taking the
     * remainder of the absolute time at which the encoder time was recorded
     * when divided by the packet interval. The encoder jumps are small enough
     * that, over the time period covered by one scan, they can be ignored for
     * purposes of calculating those absolute times:
     */
    for (scan = 0; scan < l1a_data->num_scans; scan++)
    {
	time_count[0] = (mirr_params->mirr_abs_limit[0] +
                	 mirr_params->mirr_abs_limit[1]) / 2.0;
	for (sample = 1; sample < MAX_IMPULSE_NUMBER-1; sample++)
	{
	    if (earth_encoder_times[scan][sample]==MAX_UINT16_VAL ||
		earth_encoder_times[scan][sample-1]==MAX_UINT16_VAL)
	    {
		encoder_flags[sample] = BAD_DATA;
		time_count[sample] = time_count[0];	/* Simplistic. */
	    }
	    else
	    {
		encoder_flags[sample] = GOOD_DATA;
		/* cast unsigned integers into larger, signed integers */
		/* to ensure computation performed correctly.          */
		time_count[sample] =
		    (double)((int32)earth_encoder_times[scan][sample] -
		    (int32)earth_encoder_times[scan][sample-1]);
		if(time_count[sample] < 0.0)
		    /* earth_encoder_time has rolled over  */
		    time_count[sample] += TIME_ROLL_CORRECTION;
	    }
	}

	/* Validate mirror side and encoder times */
	if ((l1a_data->mirr_data[scan].mirr_side != 0 &&
		l1a_data->mirr_data[scan].mirr_side != 1) ||
	    l1a_data->frame_data[scan].EV_start ==
		l1a_data->fill_values.EV_start_time)
	{
	    sprintf(msgbuf, "in scan %d", scan);
	    modsmf(MODIS_E_DATA_SCAN, msgbuf, filefunc);
	}
	else
	{   /* Unpack Earth sector start and vernier count */
	    offset = (l1a_data->mirr_data[scan].mirr_chan == CHAN_A
		? 0 : SIDEB_OFF);
	    impulse_enc[0] = (double)view_sector_start[scan]
		[mirr_params->sector_word[l1a_data->mirr_data[scan].mirr_side]
		    + offset];
	    vernier = (int)view_sector_start[scan]
		[mirr_params->vernier_word[l1a_data->mirr_data[scan].mirr_side]
		    + offset];
	    impulse_time[0] = - (double)(mirr_params->N_reset*t_frame 
		+ vernier*mirr_params->t_vernier);

	    /* For our initial estimates, we want to include those encoder
	     * times that contain jumps. Therefore, we relax the validation
	     * criteria:
	     */
	    mirr_abs_limit[0] = mirr_params->mirr_abs_limit[0]-32.0;
	    mirr_abs_limit[1] = mirr_params->mirr_abs_limit[1]+32.0;
	    mirr_del_limit = mirr_params->mirr_del_limit+32.0;

	    /* No need to validate magnitude of first time. However, the delta
	     * between the first and the second might still need checking.
	     */
	    if(GEO_abs_limit_check(time_count+1, MAX_IMPULSE_NUMBER-2,
		mirr_abs_limit, encoder_flags+1) != SUCCESS
	    || GEO_del_limit_check(time_count, MAX_IMPULSE_NUMBER-1,
		mirr_del_limit, encoder_flags) != SUCCESS)
	    {
		sprintf(msgbuf, "in scan %d", scan);
		modsmf(MODIS_E_DATA_SCAN, msgbuf, filefunc);
	    }
	    else for (sample = 0; sample < MAX_IMPULSE_NUMBER; sample++)
	    {
		if(sample)
		{
		    impulse_enc[sample] = impulse_enc[sample-1] +
			(double)mirr_params->sample_impulse;
		    impulse_time[sample] = impulse_time[sample-1] +
			time_count[sample-1]*mirr_params->t_encoder;
		    /* Phase of encoder sample */
                    if (sample < (MAX_IMPULSE_NUMBER - 1))
                    {
		    sorted_times[num_times][0] =
			fmod(l1a_data->frame_data[scan].EV_start 
			+ impulse_time[sample], mirr_params->packet_interval);
		    sorted_times[num_times++][1] =
			time_count[sample]*mirr_params->t_encoder;
                    }
		}
		l1a_data->mirr_data[scan].mirr_impulse_enc[sample] =
		    impulse_enc[sample];
		l1a_data->mirr_data[scan].mirr_impulse_time[sample] =
		    impulse_time[sample];
	    }
	}
    } /* End of scan loop */

    if (num_times>3)
    {	/* Next, we determine empirically which range of values for the
	 * remainders identifies the affected encoder times, and also
	 * empirically determine the size of the adjustment.
	 */
	qsort(sorted_times, num_times, sizeof(sorted_times[0]), compare_double);
	cum_sum[0] = sorted_times[0][1];
	cum_sum2[0] = cum_sum[0]*cum_sum[0];
	for (low=1; low<num_times; low++)
	{
           cum_sum[low] = cum_sum[low-1] + sorted_times[low][1];
           cum_sum2[low] = cum_sum2[low-1] +
       	                   sorted_times[low][1]*sorted_times[low][1];
	}
	/* Examine all possible intervals for the affected encoder times. The
	 * number of encoder times reported in a single granule is small enough
	 * to make this brute-force approach feasible.
	 */
	for (low=0; low+2<num_times; low++)
	{
	for (high=low+1; high+2<num_times; high++)
	{
       	    if(sorted_times[low][0]+0.0175 < sorted_times[high][0] &&
		sorted_times[high][0] < sorted_times[low][0]+0.0185)
       	    {
       		njump = high-low;
       		nnorm = num_times - njump;
       		norm_sum = cum_sum[low] + cum_sum[num_times-1] - cum_sum[high];
       		jump_sum = cum_sum[high]-cum_sum[low];
       		norm_sum2 =
       		    cum_sum2[low] + cum_sum2[num_times-1] - cum_sum2[high];
       		jump_sum2 = cum_sum2[high]-cum_sum2[low];
		/* This is the standard deviation of the time counts in the
		 * jump region, plus the standard deviation of the time counts
		 * not in that region.
		 */
       		criterion = (norm_sum2-norm_sum*norm_sum/(double)nnorm) +
       		    (jump_sum2-jump_sum*jump_sum/(double)njump);
       		if (criterion < best)
       		{ /* The right jump region is the one that minimizes that sum */
       		    best = criterion;
       		    best_low = low;
       		    best_high = high;
       		    encoder_adjustment =
       			norm_sum/(double)nnorm - jump_sum/(double)njump;
       		}
       	    }
	}
	}

	/* Determine the range of value for the remainders that correspond to
	 * tht_low and best_high values.
	 */
	gap_start =
	    0.5*(sorted_times[best_low][0] + sorted_times[best_low+1][0]);
	gap_end =
	    0.5*(sorted_times[best_high][0] + sorted_times[best_high+1][0]);
	/* Because some time was dropped, the clock rate needs to be adjusted
	 * for that fact.
	 */
	factor = 1.0 - encoder_adjustment/mirr_params->packet_interval;
    }

    /* Now that we know how to identify the affected encoder times, correct
     * those times and recalculate the absolute times accordingly.
     */
    for(scan = 0; scan < l1a_data->num_scans; scan++)
    {
	impulse_time[0] = l1a_data->mirr_data[scan].mirr_impulse_time[0];
	impulse_enc[0] = l1a_data->mirr_data[scan].mirr_impulse_enc[0];
        num_impulse = 0;
	for(sample = 1; sample < MAX_IMPULSE_NUMBER; sample++)
	{
	    if (earth_encoder_times[scan][sample] ==
		    l1a_data->fill_values.raw_mir_enc ||
	        earth_encoder_times[scan][sample-1] ==
		    l1a_data->fill_values.raw_mir_enc)
		encoder_flags[sample] = BAD_DATA;
	    else
		encoder_flags[sample] = GOOD_DATA;
	    impulse_enc[sample] =
		l1a_data->mirr_data[scan].mirr_impulse_enc[sample];
	    impulse_time[sample] = impulse_time[sample-1] +
		factor*(l1a_data->mirr_data[scan].mirr_impulse_time[sample]-
		l1a_data->mirr_data[scan].mirr_impulse_time[sample-1]);
	    if(sample>1)
	    {
		phase = fmod(l1a_data->frame_data[scan].EV_start +
			l1a_data->mirr_data[scan].mirr_impulse_time[sample-1], 
			mirr_params->packet_interval);
		if (phase > gap_start && phase < gap_end)
		    impulse_time[sample] += encoder_adjustment;
		time_count[sample-1] = (impulse_time[sample] -
		    impulse_time[sample-1]) / mirr_params->t_encoder;
	    }
	}

	/* Now that the encoder time jumps have been corrected, we can
	 * re-validate using stricter limits.
	 */
	if(GEO_abs_limit_check(time_count+1, MAX_IMPULSE_NUMBER-2,
	    (double*)mirr_params->mirr_abs_limit, encoder_flags+1) != SUCCESS
	|| GEO_del_limit_check(time_count, MAX_IMPULSE_NUMBER-1,
	    mirr_params->mirr_del_limit, encoder_flags) != SUCCESS)
	{
	    sprintf(msgbuf, "in scan %d", scan);
	    modsmf(MODIS_E_DATA_SCAN, msgbuf, filefunc);
	    l1a_data->mirr_data[scan].impulse_flag = BAD_DATA; 
	}
	for(sample = 0; sample < MAX_IMPULSE_NUMBER; sample++)
	{
            if( (encoder_flags[sample] == GOOD_DATA) ||
	        (sample && encoder_flags[sample-1] == GOOD_DATA) )
	    {
                l1a_data->mirr_data[scan].mirr_impulse_enc[num_impulse] =
		    impulse_enc[sample];
                l1a_data->mirr_data[scan].mirr_impulse_time[num_impulse] =
		    impulse_time[sample];
                num_impulse++;
            }
	}
	l1a_data->mirr_data[scan].num_impulse = num_impulse;
    }
 
    return PGS_S_SUCCESS;
}
