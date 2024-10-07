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

#include "cal_l1a_osmi.h"
#include "call1a_proto_osmi.h"
#include "InstStatData.h"
#include <mfhdf.h>

#define DEBUG 0  /* DEBUG= 0->no debug, 1->logic tracking, 2->data tracking */

/* function prototype */
int getCalDay(int year, int jday, int *month, int *day);

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
      char  *    dtype             I      data type ("SOL","TDI","IGC","LUN","")
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
        Lakshmi Kumar    Hughes STX      06/07/94    Original development for SeaWiFS

-----------------------------------------------------------------------------*/

int32_t calibrate_l1a_osmi(char *cal_path,
        int16_t syear,
        int16_t sday,
        int32_t smsec,
        int16_t eday,
        int32_t msec,
        char *dtype, /* "","SOL","TDI","IGC","LUN" */
        int32_t st_samp, /* 1 */
        int32_t nsamp, /* 1044 */
        int32_t fpixel, /* 0-95 */
        int16_t gain[4], /* quadrant gain settings from modified ISD */
        int16_t offset, /* quadrant offset settings from mod ISD */
        int16_t scan_temp, /* sensor temps from mod ISD */
        int16_t *l1a_data,
        float *l1b_data,
        cal_mod_struc *cal_mod) {

    static int32 first_call = 1;
    static int16 pr_syear = 0;
    static int16 pr_sday = 0;
    static int32 pr_smsec = 0;

    /*                                                                  */
    /* Electronic gains computed from EIDP pp 480-487, BAF              */
    /*    A,510,-    B,510,+    C,865,-    D,865,+                      */
    static float32 egain[NGAINS][NQUADS] = {
        {1.00000, 1.00000, 1.00000, 1.00000}, /* gain setting 0 */
        {1.54269, 1.52875, 1.52107, 1.57918}, /* gain setting 1 */
        {2.08538, 2.05750, 2.04215, 2.15835}, /* gain setting 2 */
        {2.62807, 2.58624, 2.56322, 2.73753}, /* gain setting 3 */
        {3.17076, 3.11499, 3.08430, 3.31671}, /* gain setting 4 */
        {3.71346, 3.64374, 3.60537, 3.89588}, /* gain setting 5 */
        {4.25615, 4.17248, 4.12645, 4.47506}, /* gain setting 6 */
        {4.79884, 4.70123, 4.64752, 5.05424} /* gain setting 7 */
    };


    /* get_cal output parameters */

    static int16 cal_year, cal_day; /* calibrate effective yr,day,msec   */
    static int32 cal_msec;
    static float32 eoffset; /* electronic offset conv to counts  */
    static float32 egain_dummy[NGAINS]; /* electronic gain dummy argument    */
    static float32 mirror[NBANDS]; /* mirror corr factors               */
    static float32 t_const[NBANDS]; /* time dependent constant term      */
    static float32 t_linear[NBANDS]; /* time dependent linear term        */
    static float32 t_quadratic[NBANDS]; /* time dependent quadratic term     */
    static float32 tempcorr[NBANDS]; /* temp correction coefficients      */
    static float32 slopes[NBANDS * 96]; /* sensor calibration slopes         */
    static float32 dcs[NBANDS * 96]; /* sensor calibration offsets        */
    static float32 scan_mod[1044]; /* scan modulation corr (sweep dir)  */

    /* Add static holders for instrument status data.  These data would normally be
       passed into this routine (gain,offset,scan_temp).  Because this data is only
       used here, we will add the routine to look this data up in the instrument 
       status data file as required.
     */
    static char isd_date[256] = "";
    char scan_date[256];
    int scan_month;
    int scan_day;
    int hr, min, sec;
    int of;

    /* Local working variables */

    //int16   do_scan_mod_correction;            /* local flag for calibration-mode   */
    int16 band, pixel; /* local band, pixel indices         */
    int16 field, quad, gs; /* local field, quad, gain indices   */
    int16 status; /* local return status               */
    int16 cal_jday, data_jday; /* local julian day counters         */

    float32 l1_count; /* local dark-current corrected cnt  */
    float32 l1b_rad; /* local radiance output             */

    float64 delta_msec, delta_t; /* local calibration model timescale */
    float32 cal_gain, cal_offset; /* local system cal gain & offset    */
    float32 slope, dark; /* local combined slope & dc         */


    /* === Load calibration file if required === */

    if (first_call ||
            syear != pr_syear ||
            sday != pr_sday ||
            smsec != pr_smsec) {

        first_call = 0;

        pr_syear = syear;
        pr_sday = sday;
        pr_smsec = smsec;

        if ((status = get_cal_osmi(cal_path,
                syear, sday, eday, msec,
                &cal_year, &cal_day, &cal_msec,
                &eoffset,
                &egain_dummy[0],
                &tempcorr[0],
                &mirror[0],
                &t_const[0],
                &t_linear[0],
                &t_quadratic[0],
                &slopes[0],
                &dcs[0],
                &scan_mod[0])) < 0) {
            return status;
        }

    }

    /* Format the scan date/time for isd look up YYYY/MM/DD HH:MM:SS*/
    hr = (int) (msec / 1000.0 / 3600.0);
    min = (int) (((msec / 1000.0 / 3600.0) - hr) * 60.0);
    sec = (int) (((((msec / 1000.0 / 3600.0) - hr) * 60.0) - min) * 60.0);
    getCalDay(syear, sday, &scan_month, &scan_day);
    sprintf(scan_date, "%d/%02d/%02d %02d:%02d:%02d", syear, scan_month, scan_day, hr, min, sec);


    /* If the scan date has changed, try and find an ISD record that matches. */
    /* Dong-Han : Error of ISD file, Sep. 2. 2000
    if(strcmp(scan_date, isd_date)) {

        if(FindISDRecordByGPSTime(scan_date, &isd, &records) != 0) {
                printf("Calibrate_l1a failed to find ISD record, using last known values\n");
                printf("%s\n", scan_date);
        }
        else {
     * Get the offsets, gain, and temp from the ISD record. *
                isd_offset[0] = isd.isd_osmi.offset1;
                isd_offset[1] = isd.isd_osmi.offset2;
                isd_offset[2] = isd.isd_osmi.offset3;
                isd_offset[3] = isd.isd_osmi.offset4;
                isd_gain = isd.isd_osmi.gain;
                isd_scan_temp = isd.isd_osmi.temp;
        }
                strcpy(isd_date, scan_date);
    }
     */
    strcpy(isd_date, scan_date);


    /* Determine the timescale (minutes) for the calibration system gain */
    cal_jday = jul2jday(cal_year, cal_day);
    data_jday = jul2jday(syear, sday);
    delta_msec = (float64) (data_jday - cal_jday)*86400000L - cal_msec;
    delta_t = (delta_msec + (float64) msec) / 60000.0;


    /* Determine if Scan Modulation corrections are required (not in cal modes) */
    //do_scan_mod_correction = (strcmp(dtype, "SOL") != 0) &&
    //                        (strcmp(dtype, "TDI") != 0) &&
    //                        (strcmp(dtype, "IGC") != 0);


    /* === Loop for each band of data === */

    if (st_samp < 1) st_samp = 1;

    fpixel = 95 - fpixel;

    for (band = 0; band < NBANDS; band++) {

        /* Determine the CCD quadrant (quad) based on band & frame pixel */
        field = fpixel / 48;
        if (band < 5)
            quad = 1 - field;
        else
            quad = 3 - field;

        /* Get gain setting (gs) from input gain argument */
        /* 
        gs = isd_gain;      
        of = isd_offset[quad];
         */
        gs = 4; /* Standard operational gain setting is 4 */
        of = 0;

        /* Compute the longterm system calibration gain factor for band */
        cal_gain = (float32) (t_const[band] +
                t_linear[band]*(delta_t) +
                t_quadratic[band]*(delta_t * delta_t));

        /* Compute the calibration offset */
        cal_offset = eoffset*of;


        /* Incorporate Longterm, Mirror Reflectance, Temperature Correction, 
           Sensor Response, and Electronics Gain factors into single  
           conversion slope and offset for current band and frame-pixel      */

        slope = cal_gain * /* longterm or override    */
                mirror[band]* /* mirror reflectance      */
                (1.0 - tempcorr[band] * scan_temp) * /* temperature correction  */
                slopes[band * 96 + fpixel] / /* sensor responsivity     */
                egain[gs][quad]; /* electronics gains       */

        dark = dcs[band * 96 + fpixel] + /* sensor offset           */
                cal_offset; /* electronics offset      */

        /*
        fprintf(stderr,
          "pixel= %d band= %d slope= %f gain= %f dark= %f egain[%d][%d]= %f\n",
          fpixel,band,slopes[band*96+fpixel],slope,dark,quad,gs,egain[gs][quad]);
         */


        /* === Loop for each input pixel -- convert to radiance === */
        for (pixel = st_samp - 1; pixel <= st_samp + nsamp - 2; pixel++) {

            /* Remove dark-current from input pixel count (trunc to 0..1023) */
            l1_count = (float32) l1a_data[band * nsamp + pixel] - dark;

            if (l1_count < 0L) l1_count = 0L;
            if (l1_count > 1023L) l1_count = 1023L;

            /* Convert counts to radiance */
            l1b_rad = l1_count*slope;


            /* Correct for scan modulation if not in a calibration mode */
            /*
            if (do_scan_mod_correction)
               l1b_rad = l1b_rad*scan_mod[pixel];
             */


            /* Output final L1B radiance value */
            l1b_data[band * nsamp + pixel] = l1b_rad;

        }

    }

    return SUCCEED;

}
