/*-----------------------------------------------------------------------------
    File : calibrate_l1a.c

    Contents:
        calibrate_l1a   -  takes an array of level-1A raw counts and returns
                           a corresponding arry of radiance values after 
                           applying the sensor calibration 

    Other relevant files:
        cal_l1a.h       -  various #defined constants, other include files 
                                (get_cal.h, getcal_proto.h, call1a_proto.h) and 
                                also includes hdf.h
        get_cal.c       -  a higher layer of calibration input functions
        get_cal_misc.c  -  a lower layer of calibration input functions 
        tcal_l1a.c      -  a test routine to test calibrate_l1a

    Notes:
        HDF library used - HDF3.3r3p1

        Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      06/07/94    Original development
        Lakshmi Kumar	 Hughes STX	 10/30/95    Filtering scan_mod corr.
                                                     to data of dtype = SOL &
                                                     TDI.
        Lakshmi Kumar	 Hughes STX	 11/20/95    Filtering scan_mod corr.
                                                     to data of dtype = IGC
        Lakshmi Kumar	 Hughes STX	 01/25/96    Modified local variable
                                                     "rads" to global var
        Gene Eplee          GSC          05/24/96    Implemented quadratic
                                                     temporal calibration
                                                     correction (gain & offset)
        Lakshmi Kumar    Hughes STX      10/17/96    Removed cal_year, cal_day
                                                     output arguments.
                                                     Ref. V5.0 I/O specs.
        Gene Eplee          GSC          01/14/97    Reworked populating of
                                                     g_f at the calibration
                                                     knees.
        Gene Eplee          GSC          04/25/97    Modified system gain and
                                                     offset corrections to
                                                     allow calling routine to
                                                     override cal table values.
        Gene Eplee          GSC          09/09/97    Modified dark restore
                                                     subtraction to use
                                                     mean values for GAC, LAC,
                                                     and HRPT data.
-----------------------------------------------------------------------------*/

#include "calib_cal_l1a.h"
#include "calib_call1a_proto.h"

/*-----------------------------------------------------------------------------
    Function: calibrate_l1a 

    Returns: int32 (status)
        Returns a status code of 0 when successful.  Otherwise returns
        -1 - to indicate calibration file open/close error
        -2 - to indicate read error
        -3 - to indicate time error (if the given time cannot be found)
        -4 - to indicate insufficient memory error

    Description:
        The function calibrate_l1a takes an array of level-1A raw counts
        and returns a corresponding array of radiance values (level-1B data)
        after applying the sensor calibration.

    Arguments: (in calling order)
      Type       Name             I/O     Description
      ----       ----             ---     -----------
      char *     cal_path          I      calibration file path
      int16      syear             I      year of data start time
      int16      sday              I      day-of-year for data start time
      int32      smsec             I      milliseconds of the day for data 
                                          start tm as returned by get_l1a_open
      int16      eday              I      day-of-year for data end time
      int32      msec              I      millisecs of the day for data start
                                          time as returned by get_l1a_record
      char  *    dtype             I      data type flag
      int32      st_samp           I      start pixel/sample number to process 
      int32      nsamp             I      samples per scan line
      int16 *    dark_rest         I      dark restore values for all 8 bands
      float32 *  dark_mean         I      mean dark restore values for all 8
                                          bands
      int16 *    gain              I      gains for all 8 bands; as returned
                                          by get_l1a_record
      int16 *    tdi               I      input TDI for all 8 bands
      int16 *    scan_temp         I      digitized scan temperatures
      int16 *    side              I      mirror side of scan line
      int16 *    l1a_data          I      Level-1A data; Band interleave by pixel
      float32 *  l1b_data          O      Sensor calibrated radiance values
                                          corresponding to l1a_data
      struct cal_mod_struc         I      Structure for override control of
                 cal_mod                  system gain and offset

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      06/07/94    Original development
        Lakshmi Kumar	 Hughes STX	 02/07/95    Bug fix--Applied time
                                                     dependent correction.
                                                     Added cal_year & cal_day
                                                     arguments. 
                                                        (ref v4.2 I/O specs.)
        Lakshmi Kumar    Hughes STX      10/30/95    Filtering scan_mod corr.
                                                     to data of dtype = SOL &
                                                     TDI.
        Lakshmi Kumar    Hughes STX      11/20/95    Filtering scan_mod corr.
                                                     to data of dtype = IGC
        Lakshmi Kumar	 Hughes STX	 01/25/96    Modified loc var "rads"
                                                     to global var "cal_rads"
                                                     to make it accessible to
                                                     l1a_read rtn.
        Gene Eplee          GSC          05/24/96    Implemented quadratic
                                                     temporal calibration
                                                     correction (gain & offset)
        Lakshmi Kumar    Hughes STX      10/17/96    Removed cal_year, cal_day
                                                     output arguments.
                                                     Ref. V5.0 I/O specs.
        Gene Eplee          GSC          01/14/97    Reworked populating of
                                                     g_f at the calibration
                                                     knees.
        Gene Eplee          GSC          04/25/97    Modified system gain and
                                                     offset corrections to
                                                     allow calling routine to
                                                     override cal table values.
        Lakshmi Kumar	    GSC		 05/02/97    Removed non prototype 
                                                     defns.
        Gene Eplee          GSC          09/09/97    Modified dark restore
                                                     subtraction to use
                                                     mean values for GAC, LAC,
                                                     and HRPT data.
-----------------------------------------------------------------------------*/
int32 calibrate_l1a(char *cal_path, int16 syear, int16 sday, int32 smsec,
        int16 eday, int32 msec, char *dtype, int32 st_samp,
        int32 nsamp, int16 *dark_rest, float32 *dark_mean, int16 *gain,
        int16 *tdi, int16 *scan_temp, int16 side, int16 *l1a_data,
        float32 *l1b_data, cal_mod_struc *cal_mod) {

    static int32 called_get_call = 0, first_call = 1;
    static int16 pr_syear = 0, pr_sday = 0;
    static int32 pr_smsec = 0;
    static int16 pr_gain[NBANDS] = {-1, -1, -1, -1, -1, -1, -1, -1};
    static int16 pr_tdi[NBANDS] = {-1, -1, -1, -1, -1, -1, -1, -1};
    static float32 pr_cal_gain[NBANDS] = {-1, -1, -1, -1, -1, -1, -1, -1};


    int16 cal_year, cal_day; /* calibrate entry year & day       */
    int16 status, n;
    int16 ref_year, ref_day, ref_min; /* calibration model reference time */
    int16 band, knee, count, pixel; /* local indeces            */
    int16 count1, count2, tdi_flag, gain_flag; /* temporary storage fields */
    int16 l1_data; /*temp space for dark-restore corrected cnt*/
    int32 ref_jday, data_jday; /* julian reference day and data day */
    float32 slope; /*loc_slope for converting l1a cnt and rad */
    static float32 cal_gain[NBANDS]; /* calibration model system gain */
    float64 delta_min; /* calibration model timescales */
    static float64 delta_t;
    static float32 g_f[NBANDS][1024]; /*gain factors look up table        */

    /*get_cal output parameters         */
    static float32 temps[256][NBANDS]; /*temp correction co_efficients     */
    static float32 scan_mod[2][1285]; /*scan modulation corr factors      */
    static float32 mirror[2][NBANDS]; /*mirror side-0 and 1 corr factors  */
    static float64 tfactor_const[NBANDS]; /* time dependent constant term */
    static float64 tfactor_linear[NBANDS]; /* time dependent linear term */
    static float64 tfactor_quadratic[NBANDS]; /* time dependent quadratic term */
    static float32 cal_offset[NBANDS]; /* calibration model system offset */

    float32 cal_counts[BANDS][4][5];
    float32 cal_rads[BANDS][4][5];

    for (band = 0; band < NBANDS; band++) {
        if (tdi[band] < 0)
            tdi[band] = 0;
        if (tdi[band] > 255)
            tdi[band] = 255;
        if (scan_temp[band] < 0)
            scan_temp[band] = 0;
        if (scan_temp[band] > 255)
            scan_temp[band] = 255;
    }

    for (band = 0; band < NBANDS && tdi[band] == pr_tdi[band]; band++);

    if (band < NBANDS)
        tdi_flag = 1;
    else
        tdi_flag = 0;

    if (first_call || syear != pr_syear || sday != pr_sday || smsec != pr_smsec
            || tdi_flag) {

        first_call = 0;
        pr_syear = syear;
        pr_sday = sday;
        pr_smsec = smsec;
        for (band = 0; band < NBANDS; band++)
            pr_tdi[band] = tdi[band];

        if ((status = get_cal(cal_path, syear, sday, eday, msec, dtype, tdi,
                &cal_year, &cal_day, &ref_year, &ref_day, &ref_min, temps,
                scan_mod, mirror, tfactor_const, tfactor_linear,
                tfactor_quadratic, cal_offset, cal_counts, cal_rads)) < 0)
            return status;
        called_get_call = 1;
#ifdef DEBUG
        printf("\ncalled get_cal\n");
#endif

        /* Define the timescale for the calibration system gain. */
        ref_jday = jul2jday(ref_year, ref_day);
        data_jday = jul2jday(syear, sday);
        delta_min = (float64) ((data_jday - ref_jday)*1440 - ref_min);
        delta_t = delta_min + (float64) msec / 60000.0;

        /* Compute the system gain */
        for (band = 0; band < NBANDS; band++) {
            cal_gain[band] = (float32) (tfactor_const[band]
                    + delta_t * tfactor_linear[band] + delta_t * delta_t * tfactor_quadratic[band]);
            if (cal_gain[band] != pr_cal_gain[band]) {
                printf("band %d:  dt = %dmin  cal_offset = %g  cal_gain = %7.5f + %.3e * dt + %.3e * dt * dt = %g\n",
                        band, (int32_t) delta_t, cal_offset[band], tfactor_const[band], tfactor_linear[band], tfactor_quadratic[band], cal_gain[band]);
                pr_cal_gain[band] = cal_gain[band];
            }
        }

    }

    for (band = 0; band < NBANDS && gain[band] == pr_gain[band]; band++);

    if (band < NBANDS)
        gain_flag = 1;
    else
        gain_flag = 0;

    if (called_get_call || gain_flag) {
        called_get_call = 0;
        for (band = 0; band < NBANDS; band++)
            pr_gain[band] = gain[band];
        for (band = 0; band < NBANDS; band++) {
#ifdef DEBUG
            printf("TDI  band gain  knee  cnt1   cnt2      slope\n");
            printf("-----------------------------------------\n\n");
#endif
            for (knee = 1; knee <= 4; knee++) {
                n = 1;
                while (((int16) cal_counts[band][gain[band]][knee] ==
                        (int16) cal_counts[band][gain[band]][knee - n]) && n <= knee)
                    n++;
                count1 = (int16) cal_counts[band][gain[band]][knee - n] + 1;
                count2 = (int16) cal_counts[band][gain[band]][knee];
                if (knee == 1)
                    count1 = 0;
                if (knee == 4)
                    count2 = 1023;
                slope = (cal_rads[band][gain[band]][knee] -
                        cal_rads[band][gain[band]][knee - n]) /
                        (cal_counts[band][gain[band]][knee] -
                        cal_counts[band][gain[band]][knee - n]);
#ifdef DEBUG
                printf("%3d  %3d  %3d  %4d  %5d  %5d    %10.8f\n",
                        tdi[band], band, gain[band], knee, count1, count2, slope);
#endif
                for (count = count1; count <= count2; count++)
                    g_f[band][count] =
                        slope * (count - cal_counts[band][gain[band]][knee - n]) +
                    cal_rads[band][gain[band]][knee - n];
            }
        }
    }

#ifdef DEBUG
    printf("\n\n------------- g_f table values --------------\n");
    printf("\nindex b1\tb2\tb3\tb4\tb5\tb6\tb7\tb8\n");
    for (i = 0; i < 1024; i++) {
        printf("\n%4d", i);
        for (j = 0; j < NBANDS; j++)
            printf("%9.5f", g_f[j][i]);
    }

#endif


    if (st_samp < 1)
        st_samp = 1;
    for (band = 0; band < BANDS; band++) {

        /* Check for override on system gain and offset */
        /* and pass back the gain, offset used  */
        if (cal_mod) {
            if (cal_mod->flag == 1) {
                cal_gain[band] = (float32) cal_mod->gain[band];
                cal_mod->offset[band] = cal_offset[band];
            } else if (cal_mod->flag == 2) {
                cal_offset[band] = (float32) cal_mod->offset[band];
                cal_mod->gain[band] = cal_gain[band];
            } else if (cal_mod->flag == 3) {
                cal_gain[band] = (float32) cal_mod->gain[band];
                cal_offset[band] = (float32) cal_mod->offset[band];
            } else {
                cal_mod->gain[band] = cal_gain[band];
                cal_mod->offset[band] = cal_offset[band];
            }
        }

        for (pixel = st_samp - 1; pixel <= st_samp + nsamp - 2; pixel++) {

            /* Use dark restore values from each scan line for calibration data.
               Use mean dark restore values for imaging data. */
            if (!strcmp(dtype, "SOL") ||
                    !strcmp(dtype, "TDI") ||
                    !strcmp(dtype, "IGC") ||
                    !strcmp(dtype, "LUN"))
                l1_data = l1a_data[pixel * BANDS + band] - dark_rest[band];
            else
                l1_data = l1a_data[pixel * BANDS + band] - (int16) (dark_mean[band] + 0.5);

            if (l1_data < 0) l1_data = 0;
            if (l1_data > 1023) l1_data = 1023;
            l1b_data[pixel * BANDS + band] = g_f[band][l1_data] * mirror[side][band] * temps[scan_temp[band]][band];
            if (strcmp(dtype, "SOL") && strcmp(dtype, "TDI") && strcmp(dtype, "IGC"))
                l1b_data[pixel * BANDS + band] *= scan_mod[band % 2][pixel];
            l1b_data[pixel * BANDS + band] =
                    l1b_data[pixel * BANDS + band] * cal_gain[band] + cal_offset[band];
            if (l1b_data[pixel * BANDS + band] > 100)
                printf("l1b_data %g  l1a_data %d  cal_gain %g  g_f %g  mirror %g  temps %g  scan_mod %g\n", l1b_data[pixel * BANDS + band], l1a_data[pixel * BANDS + band], cal_gain[band], g_f[band][l1_data], mirror[side][band], temps[scan_temp[band]][band], scan_mod[band % 2][pixel]);

        }
    }
    return SUCCEED;
}

/*---------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------

   Function: jul2jday

   Returns type: int

   Description:

      Convert year and day-of-year pair to Julian Day. This routine uses
      the same cal2jday routine by just specifying the month to be 1 and
      the day-of-month to be day-of-year.

   Parameters: (in calling order)
      Type              Name    I/O     Description
      ----              ----    ---     -----------
      int               year    I       year
      int               yday    I       day of year [1,366]

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      Frank Chen        30-Jul-1993     Original development

------------------------------------------------------------------------------*/

int jul2jday(int year, int yday) {
    int jday;

    jday = cal2jday(year, 1, yday);
    return (jday);
}

/*------------------------------------------------------------------------------

   Function: cal2jday

   Returns type: int

   Description:

      This function converts a calendar date to the corresponding Julian
      day starting at noon on the calendar date. The algorithm used is
      from Van Flandern and Pulkkinen, Ap. J. Supplement Series 41,
      November 1979, p.400
      This will also work when month is 1 and day-of-month is specified
      as day-of-year.

   Parameters: (in calling order)
      Type              Name    I/O     Description
      ----              ----    ---     -----------
      int               year    I       year
      int               month   I       month [1,12]
      int               mday    I       day of month [1,31]

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      Fred Patt         04-Nov-1992     Original written in FORTRAN(jd)
      Frank Chen        30-Jul-1993     Rewrite in C

------------------------------------------------------------------------------*/

int cal2jday(int year, int month, int mday) {
    int jday;

    jday = 367 * year
            - 7 * (year + (month + 9) / 12) / 4
            + 275 * month / 9
            + mday
            + 1721014;
    /* additional calculation is needed only for dates outside of   */
    /* the period March 1, 1900 to Feburary 28, 2100                */
    /*
    jday = jday + 15 - 3 * ((year + (month - 9) / 7) / 100 + 1) / 4;
     */
    return (jday);
}

