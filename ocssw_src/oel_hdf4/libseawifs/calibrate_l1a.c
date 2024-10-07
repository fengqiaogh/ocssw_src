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
        Gene Eplee       SAIC GSC        05/08/98    Made calibration model reference
                                                     times static.
        Gene Eplee       SAIC GSC        05/12/98    Modified cal_mod to multiply
                                                     cal_gain rather than overwrite it.
        Gene Eplee       SAIC GSC        06/01/98    Fix gain 2 / gain 3 telemetry bit flip.
        Gene Eplee       SAIC GSC        02/08/00    Fix dtype IF-Statement comparison errors
                                                     for dark_restore and scan modulation.
        Gene Eplee       SAIC GSC        12/12/00    Convert time correction to
                                                     exponential and add time-dependent 
                                                     mirror side correction.
        Gene Eplee       SAIC            03/08/04    Converted to simultaneous
                                                     exponentials for time and
                                                     mirror sides
        Gene Eplee       SAIC            07/26/06    Added focal plane
                                                     temperatures and 
                                                     instrument electronics
                                                     temperature corrections
        Gene Eplee	SAIC	        06/27/2007   Added focal plane and
                                                     instrument electronics 
                                                     reference temperatures
-----------------------------------------------------------------------------*/

#include "cal_l1a.h"
#include "call1a_proto.h"
#include <mfhdf.h>

float cal_counts[BANDS_DIMS_1A][GAINS_DIMS_1A][KNEES_DIMS_1A]; /*digital cnts (zero-offs corrected */
float cal_rads[BANDS_DIMS_1A][GAINS_DIMS_1A][KNEES_DIMS_1A]; /*radiances corresponding to knees  */

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
      int32      nsta              I      start pixel/sample number to process 
      int32      nsamp             I      samples per scan line
      float32 *  dark_mean         I      mean dark restore values for all 8
                                          bands
      int16 *    gain              I      gains for all 8 bands; as returned
                                          by get_l1a_record
      int16 *    tdi               I      input TDI for all 8 bands
      int16 *    scan_temp         I      digitized scan temperatures
      int16 *    side              I      mirror side of scan line
      int16 *    l1a_data          I      Level-1A data; as returned by 
                                          get_l1a_record
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
       Lakshmi Kumar     Hughes STX      05/02/97    Removed non prototype defns.
        Gene Eplee          GSC          09/09/97    Modified dark restore
                                                     subtraction to use
                                                     mean values for GAC, LAC,
                                                     and HRPT data.
        Gene Eplee       SAIC GSC        05/08/98    Made calibration model reference
                                                     times static.
        Gene Eplee       SAIC GSC        05/12/98    Modified cal_mod to multiply
                                                     cal_gain rather than overwrite it.
        Gene Eplee       SAIC GSC        06/01/98    Fix gain 2 / gain 3 telemetry bit flip.
        Gene Eplee       SAIC GSC        02/08/00    Fix dtype IF-Statement comparison errors
                                                     for dark_restore and scan modulation.
        Joel Gales       Futuretech      10/02/00    Change st_samp to nta
                                                     Pass nsamp,nsta,&ninc
                                                     to get_cal
        Gene Eplee       SAIC GSC        12/12/00    Convert time correction to
                                                     exponential and add time-dependent 
                                                     mirror side correction.
        B. Franz         SAIC            06/02/02    removed dark_rest param
-----------------------------------------------------------------------------*/
int32_t calibrate_l1a(char *cal_path, int16_t syear, int16_t sday, int32_t smsec,
        int16_t eday, int32_t msec, char *dtype, int32_t nsta, int32_t ninc,
        int32_t nsamp, float *dark_mean, int16_t *gain, int16_t *tdi,
        int16_t *scan_temp, float *inst_temp, int16_t side,
        int16_t *l1a_data, float *l1b_data, cal_mod_struc *cal_mod) {

    static int32 called_get_call = 0, first_call = 1;
    static int16 pr_syear = 0, pr_sday = 0;
    static int32 pr_smsec = 0;
    static int16 pr_gain[NBANDS] = {-1, -1, -1, -1, -1, -1, -1, -1};
    static int16 pr_tdi[NBANDS] = {-1, -1, -1, -1, -1, -1, -1, -1};


    int16 cal_year, cal_day; /* calibrate entry year & day       */
    int16 status, n;
    static int16 ref_year, ref_day, ref_min; /* calibration model reference time */
    int16 band, knee, count, pixel; /* local indeces            */
    int16 gn[NBANDS]; /* gain buffer */
    int16 count1, count2, tdi_flag, gain_flag; /* temporary storage fields */
    int16 l1_data; /*temp space for dark-restore corrected cnt*/
    static int32 ref_jday; /* julian reference day */
    int32 data_jday; /* julian data day */
    float32 slope; /*loc_slope for converting l1a cnt and rad */
    static float32 g_f[NBANDS][1024]; /*gain factors look up table        */

    static float32 fp_temps[256][NBANDS]; /*focal plane temperatures     */
    static float32 scan_mod[2][1285]; /*scan modulation corr factors      */

    /* Exponential Cal Table Format */
    float32 mirror[NBANDS]; /* mirror side correction */
    float32 itemp_corr[NBANDS]; /* instrument temp correction */
    float32 fptemp_corr[NBANDS]; /* fp temp correction */
    float64 cal_gain[NBANDS]; /* calibration model system gain */
    float64 delta_day, delta_t; /* calibration model timescales */
    static float64 tfactor_const[NBANDS]; /* time dependent constant term */
    static float64 tfactor_linear_1[NBANDS]; /* time dependent linear term */
    static float64 tfactor_exponential_1[NBANDS]; /* time dependent exponential term */
    static float64 tfactor_linear_2[NBANDS]; /* time dependent linear term */
    static float64 tfactor_exponential_2[NBANDS]; /* time dependent exponential term */
    static float64 cal_offset[NBANDS]; /* calibration model system offset */
    static float64 inst_tcorr[NBANDS]; /* instrument temp corrections */
    static float64 inst_tref[NBANDS]; /* instrument reference temp */
    static float64 fp_tcorr[NBANDS]; /* fp temp corrections */
    static float64 fp_tref[NBANDS]; /* fp reference temp */
    static float64 mside1_const[NBANDS]; /* mirror side1 constant term */
    static float64 mside1_linear_1[NBANDS]; /* mirror side1 linear term */
    static float64 mside1_exponential_1[NBANDS]; /* mirror side1 exponential term */
    static float64 mside1_linear_2[NBANDS]; /* mirror side1 linear term */
    static float64 mside1_exponential_2[NBANDS]; /* mirror side1 exponential term */
    static float64 mside2_const[NBANDS]; /* mirror side2 constant term */
    static float64 mside2_linear_1[NBANDS]; /* mirror side2 linear term */
    static float64 mside2_exponential_1[NBANDS]; /* mirror side2 exponential term */
    static float64 mside2_linear_2[NBANDS]; /* mirror side2 linear term */
    static float64 mside2_exponential_2[NBANDS]; /* mirror side2 exponential term */

    for (band = 0; band < NBANDS; band++) {
        if (tdi[band] < 0)
            tdi[band] = 0;
        if (tdi[band] > 255)
            tdi[band] = 255;
        if (scan_temp[band] < 0)
            scan_temp[band] = 0;
        if (scan_temp[band] > 255)
            scan_temp[band] = 255;
        /* Fix gain 2 / gain 3 telemetry bit flip */
        switch (gain[band]) {
        case 0:
            gn[band] = 0;
            break;
        case 1:
            gn[band] = 2;
            break;
        case 2:
            gn[band] = 1;
            break;
        case 3:
            gn[band] = 3;
            break;
        }

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

        if ((status = get_cal(cal_path, syear, sday, eday, msec, nsamp, nsta,
                ninc, dtype, tdi, &cal_year, &cal_day, &ref_year,
                &ref_day, &ref_min, fp_temps, scan_mod, tfactor_const,
                tfactor_linear_1, tfactor_exponential_1, tfactor_linear_2,
                tfactor_exponential_2, cal_offset, inst_tcorr, inst_tref,
                fp_tcorr, fp_tref, mside1_const, mside1_linear_1,
                mside1_exponential_1, mside1_linear_2,
                mside1_exponential_2, mside2_const, mside2_linear_1,
                mside2_exponential_1, mside2_linear_2,
                mside2_exponential_2, cal_counts, cal_rads)) < 0)
            return status;
        /* Define the reference timescale for the calibration system gain. */
        ref_jday = jul2jday(ref_year, ref_day);
        called_get_call = 1;
    }

    for (band = 0; band < NBANDS && gn[band] == pr_gain[band]; band++);

    if (band < NBANDS)
        gain_flag = 1;
    else
        gain_flag = 0;

    if (called_get_call || gain_flag) {
        called_get_call = 0;
        for (band = 0; band < NBANDS; band++)
            pr_gain[band] = gn[band];
        for (band = 0; band < NBANDS; band++) {
            for (knee = 1; knee <= 4; knee++) {
                n = 1;
                while (((int16) cal_counts[band][gn[band]][knee] ==
                        (int16) cal_counts[band][gn[band]][knee - n]) && n <= knee)
                    n++;
                count1 = (int16) cal_counts[band][gn[band]][knee - n] + 1;
                count2 = (int16) cal_counts[band][gn[band]][knee];
                if (knee == 1)
                    count1 = 0;
                if (knee == 4)
                    count2 = 1023;
                slope = (cal_rads[band][gn[band]][knee] -
                        cal_rads[band][gn[band]][knee - n]) /
                        (cal_counts[band][gn[band]][knee] -
                        cal_counts[band][gn[band]][knee - n]);
                for (count = count1; count <= count2; count++)
                    g_f[band][count] =
                        slope * (count - cal_counts[band][gn[band]][knee - n]) +
                    cal_rads[band][gn[band]][knee - n];
            }
        }
    }

    /* Define the timescale for the calibration system gain. */
    data_jday = jul2jday(syear, sday);
    delta_day = (float64) (data_jday - ref_jday) - (float64) ref_min / 1440.0;
    delta_t = delta_day + (float64) msec / 86400000.0;

    for (band = 0; band < BANDS; band++) {

        /* Compute the system gain */
        cal_gain[band] = 1.0 / (tfactor_const[band] -
                tfactor_linear_1[band]*(1.0 - exp(-tfactor_exponential_1[band] * delta_t)) -
                tfactor_linear_2[band]*(1.0 - exp(-tfactor_exponential_2[band] * delta_t)));

        /* Check for override on system gain and offset */
        /* and pass back the gain, offset used  */
        if (cal_mod->flag == 1) {
            cal_gain[band] = cal_mod->gain[band] * cal_gain[band];
        } else if (cal_mod->flag == 2) {
            cal_offset[band] = cal_mod->offset[band];
        } else if (cal_mod->flag == 3) {
            cal_gain[band] = cal_mod->gain[band] * cal_gain[band];
            cal_offset[band] = cal_mod->offset[band];
        }            /*  return final gain and offset in use only if not externally altered */
        else {
            cal_mod->gain[band] = cal_gain[band];
            cal_mod->offset[band] = cal_offset[band];
        }

        /* Compute the mirror side correction */
        if (side == 0) mirror[band] =
                (float32) (1.0 / (mside1_const[band] -
                mside1_linear_1[band]*(1.0 - exp(-mside1_exponential_1[band] * delta_t)) -
                mside1_linear_2[band]*(1.0 - exp(-mside1_exponential_2[band] * delta_t))));
        else mirror[band] =
                (float32) (mside2_const[band] -
                mside2_linear_1[band]*(1.0 - exp(-mside2_exponential_1[band] * delta_t)) -
                mside2_linear_2[band]*(1.0 - exp(-mside2_exponential_2[band] * delta_t)));

        /* Compute the focal plane temperature correction. */
        fptemp_corr[band] = (float32) (1.0 +
                fp_tcorr[band]*(fp_temps[scan_temp[band]][band] - fp_tref[band]));

        /* Compute the instrument electronics temperature correction.
           The instrument electronics temperature is entry 15 in the 
           inst_ana  telemetry array. */
        itemp_corr[band] = (float32) (1.0 +
                inst_tcorr[band]*(inst_temp[14] - inst_tref[band]));

        for (pixel = 0; pixel < nsamp; pixel++) {
            /* Use mean dark restore values for all data. */
            l1_data = l1a_data[band * nsamp + pixel] - (int16) (dark_mean[band] + 0.5);
            if (l1_data < 0)
                l1_data = 0;
            if (l1_data > 1023)
                l1_data = 1023;
            l1b_data[band * nsamp + pixel] = g_f[band][l1_data];
            l1b_data[band * nsamp + pixel] =
                    l1b_data[band * nsamp + pixel] * mirror[band];
            l1b_data[band * nsamp + pixel] =
                    l1b_data[band * nsamp + pixel] * fptemp_corr[band];
            l1b_data[band * nsamp + pixel] =
                    l1b_data[band * nsamp + pixel] * itemp_corr[band];
            if ((strcmp(dtype, "SOL") != 0) && (strcmp(dtype, "TDI") != 0) && (strcmp(dtype, "IGC") != 0))
                l1b_data[band * nsamp + pixel] =
                    l1b_data[band * nsamp + pixel] * scan_mod[band % 2][pixel];
            l1b_data[band * nsamp + pixel] =
                    (float32) ((float64) l1b_data[band * nsamp + pixel] * cal_gain[band]
                    + cal_offset[band]);
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


