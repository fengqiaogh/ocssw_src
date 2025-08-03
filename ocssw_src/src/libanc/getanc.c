/*-----------------------------------------------------------------------------
    File:  getanc.c
      
    Contents:
        get_ancillary    - will open ancillary data product or climatology
                           files when required, retrieve the ancillary data, 
                           perform the interpolations, and close any open 
                           files.
        set_files	 - sets file1 and file2 from input file names depend-
                                ing upon syear, sday and eday
        set_file_dt      - sets file1, file2 and dt1 and dt2 from input file 
                           names depending upon the time of the current scane line
        ck_files_in_buf  - checks whether the requested data is already in
                           the buffers and if so, sets the rd_flag.
        read_climatology - reads climatology data file
        read_NRT         - reads ancillary NRT data file/s
        extract_data_pts - extracts 4 surrounding coords and corresponding 
                           data for the given lat/lon
        gregor		 - converts julian day into gregorian month/day
        interpolate      - initializes/sets parameters and calls dataintp rtn.
        dataintp_        - returns interpolated value for the given surrounding

    Notes:

    Other relevant files:
        anc.h  	      - various #defined constants for ancillary data, also 
                        includes hdf.h
        ancproto.h    - prototypes for ancillary data input routines
        HDFroutines.c - a lower layer of ancillary data input functions
        dataintp.f    - interpolation routine

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      04/02/93    Original development
        Lakshmi Kumar    Hughes STX      12/07/93    Modified to incorporate
                                                     HDF3.3r1 changes 
        Lakshmi Kumar    Hughes STX      01/04/94    Corrected dims bug 
        Lakshmi Kumar    Hughes STX      06/21/94    Incorporated I/O Specs 
                                                        v3.0 and v3.2 changes 
        Lakshmi Kumar	 HITC		 05/25/95    Modified to read data
                                                     from the new redesigned
                                                     ancillary data files
                                                     (Ref: I/O specs V4.3 &
                                                        product specs V2.7)
        Lakshmi Kumar	 Hughes STX	 08/10/95    Initialized DT1 & DT2 to
                                                     zero
        Lakshmi Kumar	 Hughes STX	 10/10/95    Removed PA_flag from 
                                                     input arguments to rtn
                                                     "get_ancillary".  Changed
                                                     time checking logic
                                                     (ref: V4.4 I/O specs)
        Lakshmi Kumar  	 Hughes STX      11/27/95    Changed the order of data 
                                                     pts in data_list1 & 
                                                     data_list2
        Lakshmi Kumar	 Hughes STX	 11/28/95    changed the logic of 
                                                     choosing files based on
                                                     file's start and end times
                                                     rather than its mean time
                                                     Changed extract_data_pts
                                                     logic to send 4 diff pts
                                                     always to interp rtn
        Lakshmi Kumar	 Hughes STX	 03/08/96    changed ftime structure
                                                     and moved it from hfile
                                                     to this file to avoid 
                                                     re inclusion warning msgs.
        Lakshmi Kumar    Hughes STX      04/18/96    Removed a redundant loop.
        Lakshmi Kumar    Hughes STX      08/02/96    Removed unused variables
                                                     & redefined MAX as MAXVAL
        Lakshmi Kumar  	 Hughes STX	 10/24/96    when parm_flag = 3, 
                                                     precipitable water is 
                                                     accessed instead of 
                                                     humidity (ref. V3.0 
                                                     product specs.)
        Lakshmi Kumar	 Hughes STX	 01/17/97    made changes to use range
                                                     & def values of humidity
                                                     when PA_flag is set and
                                                     parm_flag = 3, as clima-
                                                     tology files at present
                                                     contain humidity not 
                                                     precipitable water. 
        Lakshmi Kumar	 Hughes STX	 02/25/97    Made changes to allow 
                                                     access to both humidity
                                                     & precipitable water.
                                                     Accesses humidity when
                                                     param_flag = 5 and
                                                     precipitable water when
                                                     param_flag = 3.  Removed
                                                     non-prototype declarations
        Lakshmi Kumar    Hughes STX      07/07/97    Added 'get_ancillary_' to
                                                     allow Fortran interface.
        Lakshmi Kumar	 Hughes STX	 08/25/97    Incorporated Miami's 
                                                     binary search.
        Ewa Kwiatkowska  GSC		 11/12/99    Time reference for TOMS data
                                                     added.
        Don Shea         SAIC             3/20/09    moved function BinarySearch
                                                     so prototype not necessary
                                                     in ancproto.h
        Arya Akmal       SSAI             3/10/25    Check in function anc_get_time
                                                     if Start Day attribute exists
                                                     else halt execution

-----------------------------------------------------------------------------*/

#include <mfhdf.h>
#include "anc.h"
#include "ancproto.h"

#include <genutils.h>

/* Private Global Variables */
static int32 fid1, fid2, sdfid1, sdfid2;
static void *parm_buf1, *parm_buf2;
char ERR_MSG[1024];

#define FLIMIT  12

#define BS_INCR 0       /* increasing order */
#define BS_DECR 1       /* decreasing order */

static int BinarySearch(int dim, float *bufp, float val, int order,
        int16 *lower, int16 *upper) {
    int lo, hi, mid;
    int n = 0;

    lo = 0;
    hi = dim - 1;

    /* quick checks to see if value still bounded in last interval */
    /* and some additional optimizations for (just) outside last interval */
    if (*lower >= 0 && *upper >= 0) {
        int Lo, Hi;

        Lo = *lower;
        Hi = *upper;
        if (order == BS_INCR) {
            if (bufp[Lo] <= val && val <= bufp[Hi]) {
                /* bounded in last interval */
                return n;
            }
            if (val < bufp[Lo]) {
                if (Lo > 0 && bufp[Lo - 1] <= val && val <= bufp[Lo]) {
                    /* slide one index lower */
                    *lower = Lo - 1;
                    *upper = Lo;
                    return n;
                }
                hi = Lo; /* bounded above */
            } else {
                if (Hi < dim - 1 && bufp[Hi] <= val && val <= bufp[Hi + 1]) {
                    /* slide one index higher */
                    *lower = Hi;
                    *upper = Hi + 1;
                    return n;
                }
                lo = Hi; /* bounded below */
            }
        } else {
            if (bufp[Lo] >= val && val >= bufp[Hi]) {
                /* bounded in last interval */
                return n;
            }
            if (val > bufp[Lo]) {
                if (Lo > 0 && bufp[Lo - 1] >= val && val >= bufp[Lo]) {
                    /* slide one index lower */
                    *lower = Lo - 1;
                    *upper = Lo;
                    return n;
                }
                hi = Lo; /* bounded above */
            } else {
                if (Hi < dim - 1 && bufp[Hi] >= val && val >= bufp[Hi + 1]) {
                    /* slide one index higher */
                    *lower = Hi;
                    *upper = Hi + 1;
                    return n;
                }
                lo = Hi; /* bounded below */
            }
        }
    }

    /* perform standard binary search */
    if (order == BS_INCR) {
        /* ascending order list */
        /* bufp[lower] < val < bufp[upper] */
        /* or bufp[lower=upper] = val */
        if (val < bufp[0] || val > bufp[dim - 1]) {
            *lower = *upper = -1;
            return n;
        }
        for (;;) {
            mid = (hi + lo) / 2;
            n++;
            if (bufp[mid] > val)
                hi = mid;
            else
                lo = mid;
            if (lo == hi - 1) {
                *lower = lo;
                *upper = hi;
                return n;
            }
        }
    } else {
        /* descending order list */
        /* bufp[lower] > val > bufp[upper] */
        /* or bufp[lower=upper] = val */
        if (val > bufp[0] || val < bufp[dim - 1]) {
            *lower = *upper = -1;
            return n;
        }
        for (;;) {
            mid = (hi + lo) / 2;
            n++;
            if (bufp[mid] < val)
                hi = mid;
            else
                lo = mid;
            if (lo == hi - 1) {
                *lower = lo;
                *upper = hi;
                return n;
            }
        }
    }
}

/*-----------------------------------------------------------------------------
    Function:  get_ancillary

    Returns:   int32 (status)
        The return code is a negative value if any error occurs.

    Description:  
        The function get_ancillary calls appropiate functions to perform
        the open, read, and closing of the required files and to perform
        the interpolations for the given coordinates.  
        See "INTERFACE SPECIFICATIONS FOR SeaWiFS OPERATIONAL PRODUCT
        INPUT AND OUTPUT SOFTWARE" (By F. Patt, J. Firestone and M. Darzi)
        and SeaWiFS ANCILLARY DATA FORMAT (HDF) (By J. Firestone, M. Darzi)
        for details.

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        float32   lat	       I    An array of latitudes
        float32   lon          I    An array of longitudes
        int16     nsamp	       I    Count of pixels (no. of lats/lons)
        int16     syear        I    year of data start time
        int16     sday         I    day-of-year for data start time
        int16     eday         I    day-of-year for data end time
        int32     msec	       I    time in milliseconds
        char *    filename1    I    name of the ancillary data product for the
                                     nearest time preceding the scan time or 
                                     the corresponding climatology 
        char *    filename2    I    name of the ancillary data product for the
                                     nearest time following the scan time
                                     if = NULL then read file1 for climatology
        char *    filename3    I    name of the anc data product for the near-
                                     est tm following the l-1A's last scan-
                                     line's time 
        char *    anc_cor_file I    ancillary data correction file
        int16     Parm_flag    I    Flag indicating the parameter whose data
                                     are to be returned in interp: 
                                        0 - zonal wind component(meters/sec),
                                        1 - meridional wind component,
                                        2 - surface pressure (millibars)
                                        3 - precipitable water(kg m^-2)
                                        4 - total ozone (Dobson units).
                                        5 - relative humidity (percent)
        float *   interp       O    Interpolated data points
        float *   anc_unc      O    ancillary data value uncertainty
        int16 *   qcflag       O    array corresponding to interp values,     
                                     containing flags indicating the integrity
                                     of those values:	
                                        if = 0: all grid pts used in 
                                                interpolation were not suspect
                                        if = 1: suspect grid points used in
                                                interpolation


    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     04/02/93  Original development
    Lakshmi Kumar  Hughes STX  	  12/09/93  Incorporated HDF3.3r1 changes and 
                                            updated comments.
    Lakshmi Kumar  Hughes STX     01/04/94  Corrected dimensions problem 
    Lakshmi Kumar  Hughes STX     06/21/94  Incorporated I/O Specs v3.0 and 
                                            v3.2 changes 
    Lakshmi Kumar  Hughes STX	  10/10/05  Removed PA_flag as input argument
    Lakshmi Kumar  Hughes STX	  10/24/96  when parm_flag = 3, precipitable 
                                            water is accessed instead of 
                                            humidity (ref. V3.0 product specs.)
    Lakshmi Kumar  Hughes STX	  02/25/97  Added parm_flag = 5 for accessing
                                            relative humidity
----------------------------------------------------------------------------*/

int get_ancillary(float *lat, float *lon, int16_t nsamp, int16_t syear,
        int16_t sday, int16_t eday, int32_t msec, char *filename1,
        char *filename2, char *filename3,char *par_anc_file, char *anc_cor_file,
        int16_t parm_flag, float *interp, float *anc_unc,
        int16_t *qcflag) {

    int32 i, ret = 0;
    int16 PA_flag = 0, read_flag, day, month=0;
    intn toms; /* BOOLEAN value whether EP TOMS data are used */
    typedef int8 qc_buf_t[NPARMS][LAT * LON];
    static qc_buf_t *qc_buf; /* NPARMS = 6(uvppoh) */
    typedef float32 uvpph_buf_t[NPARMS - 1][LAT * LON];
    static uvpph_buf_t *uvpph_buf; /* oz has its own buf */
    static float32 *c_uvpph_lat; /* buf for clim lats */
    static float32 *c_uvpph_lon; /* buf for clim lons */
    static float32 *nrt_uvpph_lat; /* buf for nrt  lats */
    static float32 *nrt_uvpph_lon; /* buf for nrt lons  */
    static float32 *oz_lat; /* buf for oz  lats  */
    static float32 *oz_lon; /* buf for oz  lons  */
    static float32 *lat_buf;
    static float32 *lon_buf;
    typedef int16 oz_buf_t[LAT * LON];
    static oz_buf_t *oz_buf;
    static int32 oz_dims[2] = {0, 0}; /* array to hold dimensions       */
    static int32 c_uvpph_dims[2] = {0, 0}; /*array to hold dimensions       */
    static int32 nrt_uvpph_dims[2] = {0, 0}; /* array to hold dimensions     */
    int32 *dims; /* array to hold dimensions       */
    static int16 error_flag = -1; /* flag indicating prev read error */
    timech DT1, DT2; /* delta time1 and time 2 for the scan line */
    float32 lat_list[4]; /* surrounding lats of -1 time step*/
    float32 lon_list[4]; /* surrounding lons of -1 time step*/
    char file1[MAXVAL], file2[MAXVAL]; /* local variable for input files  */
    char *FUNC = "get_ancillary";

    // allocate data
    if (qc_buf == NULL) {
        qc_buf = (qc_buf_t*) allocateMemory(NFILE * sizeof (qc_buf_t), "qc_buf");
        uvpph_buf = (uvpph_buf_t*) allocateMemory(NFILE * sizeof (uvpph_buf_t), "uvpph_buf");
        c_uvpph_lat = (float32*) allocateMemory(LAT * sizeof (float32), "c_uvpph_lat");
        c_uvpph_lon = (float32*) allocateMemory(LON * sizeof (float32), "c_uvpph_lon");
        nrt_uvpph_lat = (float32*) allocateMemory(LAT * sizeof (float32), "nrt_uvpph_lat");
        nrt_uvpph_lon = (float32*) allocateMemory(LON * sizeof (float32), "nrt_uvpph_lon");
        oz_lat = (float32*) allocateMemory(LAT * sizeof (float32), "oz_lat");
        oz_lon = (float32*) allocateMemory(LON * sizeof (float32), "oz_lon");
        oz_buf = (oz_buf_t*) allocateMemory(NFILE * sizeof (oz_buf_t), "oz_buf");
    }

    /* need to init these, to avoid test of unititialized values when not TOMS */
    toms = FAIL;
    DT1.start = 0.0;
    DT1.inc = 0.0;
    DT2.start = 0.0;
    DT2.inc = 0.0;

    /*** check filenames */
    if (filename1 == NULL || filename1[0] == 0) {
        printf("****ERROR: %s: File name error: filename1 equals NULL", FUNC);
        return FAIL;
    }

    if (filename2 == NULL || filename2[0] == 0)
        PA_flag = 1;
    else
        if (filename3 == NULL || filename3[0] == 0) {
        printf("****ERROR: %s: File name error: filename3 equals NULL ", FUNC);
        return FAIL;
    }

    for (i = 0; i < nsamp; i++)
        qcflag[i] = -1;

    if (parm_flag < 0 || parm_flag > 10) {
        printf("****ERROR: %s: Invalid parm_flag input ", FUNC);
        return FAIL;
    }

    if (parm_flag == OZONE) {
        parm_buf1 = (void *) oz_buf[F1];
        parm_buf2 = (void *) oz_buf[F2];
        lat_buf = (float *) oz_lat;
        lon_buf = (float *) oz_lon;
        dims = oz_dims;
    } else {
        if (parm_flag == 5) { /* relative humidity */
            parm_buf1 = (void *) uvpph_buf[F1][4];
            parm_buf2 = (void *) uvpph_buf[F2][4];
        } else {
            parm_buf1 = (void *) uvpph_buf[F1][parm_flag];
            parm_buf2 = (void *) uvpph_buf[F2][parm_flag];
        }
        if (PA_flag) {
            lat_buf = (float *) c_uvpph_lat;
            lon_buf = (float *) c_uvpph_lon;
            dims = c_uvpph_dims;
        } else {
            lat_buf = (float *) nrt_uvpph_lat;
            lon_buf = (float *) nrt_uvpph_lon;
            dims = nrt_uvpph_dims;
        }
    }

    if (PA_flag) {
        gregor(sday, syear, &month, &day); /* convert julian day to   */
        strcpy(file1, filename1); /* gregorian month and day */
    } else
        if ((set_files(parm_flag, syear, sday, eday, msec,
            filename1, filename2, filename3, file1, file2, &DT1, &DT2, &toms)) < 0) {
        printf("\n****ERROR: %s: %s\n", FUNC, ERR_MSG);
        return FAIL;
    }

    if ((ck_files_in_buf(PA_flag, parm_flag, file1, file2, month,
            &error_flag, &read_flag)) < 0) {
        printf("\n****ERROR: %s: %s\n", FUNC, ERR_MSG);
        return FAIL;
    }

    if (read_flag) {

        if (PA_flag == 1)
            ret = read_climatology(file1, parm_flag, month,
                dims, lat_buf, lon_buf, parm_buf1);
        else
            ret = read_NRT(file1, file2, anc_cor_file, parm_flag, read_flag, dims,
                lat_buf, lon_buf, parm_buf1, parm_buf2,
                qc_buf[F1][parm_flag], qc_buf[F2][parm_flag]);
    }


    if (ret >= 0)
        extract_data_pts(PA_flag, parm_flag, DT1, DT2, nsamp, lat, lon, lat_buf,
            lon_buf, dims, qc_buf[F1][parm_flag], qc_buf[F2][parm_flag],
            parm_buf1, parm_buf2, toms, lat_list, lon_list, qcflag,
            interp, anc_unc);

    else {
        error_flag = parm_flag;
        for (i = 0; i < nsamp; i++)
            interp[i] = 0;
        printf("\n****ERROR: %s: %s\n", FUNC, ERR_MSG);
        return FAIL;
    }

    return SUCCEED;
}

intn get_ancillary_(float32 *lat, float32 *lon, int16 nsamp, int16 syear,
        int16 sday, int16 eday, int32 msec, char *filename1,
        char *filename2, char *filename3,char *par_anc_file, char *anc_cor_file,
        int16 parm_flag, float32 *interp, float32 *anc_unc,
        int16 *qcflag) {
    return get_ancillary(lat, lon, nsamp, syear, sday, eday, msec,
            filename1, filename2, filename3,par_anc_file, anc_cor_file,
            parm_flag, interp, anc_unc, qcflag);
}

/*----------------------------------------------------------------------------
    Function: check_on_TOMS 

    Returns: FAIL for applying the conventional approximation scheme;
             SUCCEED for applying EP TOMS with nearest-time approximation scheme.

    Arguments: (in calling order)
      Type         Name      I/O     Description
      ----         ----      ---     -----------
     int16      parm_flag     I      indicates wh parm it is being processed 
     char *     filename1     I      name of the ancillary data product for the
                                     nearest time preceding the scan time or 
                                     the corresponding climatology 
     char *     filename2     I      name of the ancillary data product for the
                                     nearest time following the scan time
                                     if = NULL then read file1 for climatology
     char *     filename3     I      name of the anc data product for the near-
                                     est tm following the l-1A's last scan-
                                     line's time 
     float32 *  in_lonlist    I      input longitude array
     int16      nsamp         I      number of scan line points
     float64    s_jd1         I      Julian time for the start of auxilliary filename1
     float64    s_jd2         I      Julian time for the start of auxilliary filename2
     float64    s_jd3         I      Julian time for the start of auxilliary filename3
     float64    e_jd1         I      Julian time for the end of auxilliary filename1 
     float64    e_jd2         I      Julian time for the end of auxilliary filename2
     float64    e_jd3         I      Julian time for the end of auxilliary filename3
     float64    d_jd          I      Julian time for the data scan line
     timech  *  dt1           O      Time diff bet actual time & time step 1
     timech  *  dt2    	      O      Time diff bet actual time & time step 2

    History:
      Programmer     Organization    Date      Description of change
      -------------- ------------    --------  ---------------------
      Ewa Kwiatkowska   GSC	      11/12/99  Function fully programmed
 */

int check_on_TOMS(int16_t parm_flag, char *file1, char *file2,
        char *filename1, char *filename2, char *filename3, double s_jd1,
        double s_jd2, double s_jd3, double e_jd1, double e_jd2,
        double e_jd3, double d_jd, timech *dt1, timech *dt2) {


    if (strcmp(file1, file2) == 0) {
        //printf ("\n\n check_on_TOMS failed" );
        return FAIL;
    }
    if (parm_flag != OZONE) {
        //printf ("\n\n check_on_TOMS failed" );
        return FAIL;
    };

    //length = strlen( file1 );
    //if (strcmp (((char *)(&(file1[length-10]))), "TOMS.OZONE" ) != 0) {

    if ((strstr(file1, "TOMS") == NULL) &&
            (strstr(file1, "OMI") == NULL)) {
        //printf ("\n\n check_on_TOMS detected non TOMS OZONE file" );
        return FAIL;
    };

    //length = strlen( file2 );
    //if (strcmp (((char *)(&(file2[length-10]))), "TOMS.OZONE" ) != 0) {

    if ((strstr(file2, "TOMS") == NULL) &&
            (strstr(file2, "OMI") == NULL)) {
        //printf ("\n\n check_on_TOMS detected non TOMS OZONE file" );
        return FAIL;
    };


    if ((strcmp(file1, filename1)) == 0) {

        dt1->start = d_jd - (e_jd1 + s_jd1) / 2.0;
        dt1->inc = (e_jd1 - s_jd1) / 360.0; /*358.75;*/
    } else
        if ((strcmp(file1, filename2)) == 0) {

        dt1->start = d_jd - (e_jd2 + s_jd2) / 2.0;
        dt1->inc = (e_jd2 - s_jd2) / 360.0; /*358.75;*/
    } else
        return FAIL;

    if ((strcmp(file2, filename2)) == 0) {

        dt2->start = d_jd - (e_jd2 + s_jd2) / 2.0;
        dt2->inc = (e_jd2 - s_jd2) / 360.0; /*358.75;*/
    } else
        if ((strcmp(file2, filename3)) == 0) {

        dt2->start = d_jd - (e_jd3 + s_jd3) / 2.0;
        dt2->inc = (e_jd3 - s_jd3) / 360.0; /*358.75;*/
    } else
        return FAIL;



    //printf ("\n\n Check on TOMS has been successful\n\n" ) ; exit(1);

    return SUCCEED;

}

/*----------------------------------------------------------------------------
    Function: read_climatology

    Returns:  int32 (Status) 
       The return code is FAIL (-1) if an error occurs.  Otherwise, SUCCEED (0) 

    Description:
       The function read_climatology will be called if the climatology file
       needs to be read.  It opens the file and calls appropriate functions
       to read the requested months data in to the given data buffer and
       lat/lon to the latitude and longitude buffers.  

    Arguments: (in calling order)
      Type         Name      I/O     Description
      ----         ----      ---     -----------
      char *       file1      I      Climatology file name 
      int16        parm_flag  I	     Parameter flag identifying the parameter 
                                     to be read 
      int16        month      I      data month
      int32 *      dims       O      dimensions of the data
      float32 *    lat_buf    O      buffer containing latitude data
      float32 *    lon_buf    O      buffer containing longitude data
      float32 *    parm_buf   O      buffer containing parameter data

    Modification history:
      Programmer     Organization    Date      Description of change
      -------------- ------------    --------  ---------------------
      Lakshmi Kumar  Hughes STX      04/02/93  Original development
      Lakshmi Kumar  Hughes STX      12/09/93  Incorporated HDF3.3r1 
-----------------------------------------------------------------------------*/
int read_climatology(char *file1, int16_t parm_flag, int16_t month,
        int32_t *dims, float *lat_buf,
        float *lon_buf, void *parm_buf) {
    if ((openHDF(file1, &sdfid1, &fid1)) < 0)
        return FAIL;

    if ((get_clim_data(fid1, sdfid1, parm_flag, month,
            dims, lat_buf, lon_buf, parm_buf)) < 0)
        return FAIL;

    if ((closeHDF(sdfid1, fid1)) < 0)
        return FAIL;
    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function: read_NRT

    Returns:  int32 (Status)
       The return code is FAIL (-1) if an error occurs.  Otherwise, SUCCEED (0)

    Description:
       The function read_NRT will be called if the Near Real Time file/s
       needs to be read.  It opens the file and calls appropriate functions
       to read the requested data in to the given data buffer and
       lat/lon to the latitude and longitude buffers.

    Arguments: (in calling order)
      Type         Name      I/O     Description
      ----         ----      ---     -----------
      char *       file1      I      NRT file with the nearest time preceding
                                        the scan time
      char *       file2      I      NRT file with the nearest time following
                                        the scan time
      char *       anc_cor_file  I   ancillary data correction file
      int16        parm_flag  I      Parameter flag identifying the parameter
                                     to be read
      int16        rd_flag    I      flag indicating wh file to read:
                                      1 for file1, 2 for file2 and 3 for
                                      file1 and file2
      int32 *      dims       O      dimensions of the data
      float32 *    lat_buf    O      buffer containing latitude data
      float32 *    lon_buf    O      buffer containing longitude data
      float32 *    parm_buf1  O      buffer containing parameter data of file1
      float32 *    parm_buf2  O      buffer containing parameter data of file2
      int8    *    qc_buf1    O      buffer containing qc flags of the given 
                                        parameter from file1
      int8    *    qc_buf2    O      buffer containing qc flags of the given
                                        parameter from file2

    Modification history:
      Lakshmi Kumar  Hughes STX      04/02/93  Original development
     
      Lakshmi Kumar  Hughes STX      12/09/93  Incorporated HDF3.3r1 changes
      W. Robinson, SAIC 4 Dec 2013  added code to handle anc array size 
                         mismatch for ozone (in lon dimension) but fail for 
                         other products and lat # mismatch
                         also pass in name of ancillary correction file
------------------------------------------------------------------------------*/
int read_NRT(char *file1, char *file2, char *anc_cor_file, int16_t parm_flag,
        int16_t rd_flag, int32_t *dims, float *lat_buf,
        float *lon_buf, void *parm_buf1, void *parm_buf2,
        int8_t *qc_buf1, int8_t *qc_buf2) {
    int32 i;
    static float32 *lon_buf_p4;
    static int32 lon_dim_p4;

    if (lon_buf_p4 == NULL) {
        lon_buf_p4 = (float32*) allocateMemory(LON * sizeof (float32), "lon_buf_p4");
    }

    if (rd_flag == 1 || rd_flag == 3) {
        if ((openHDF(file1, &sdfid1, &fid1)) < 0)
            return FAIL;
        if ((get_NRT_data(fid1, sdfid1, anc_cor_file, parm_flag, dims,
                lat_buf, lon_buf, parm_buf1, qc_buf1)) < 0)
            return FAIL;
        if ((closeHDF(sdfid1, fid1)) < 0)
            return FAIL;
        if (parm_flag == 4) {
            lon_dim_p4 = dims[1];
            for (i = 0; i < dims[1]; i++)
                lon_buf_p4[i] = lon_buf[i];
        }
    }


    if (rd_flag == 2 || rd_flag == 3) {
        if ((openHDF(file2, &sdfid2, &fid2)) < 0)
            return FAIL;
        if ((get_NRT_data(fid2, sdfid2, anc_cor_file, parm_flag, dims,
                lat_buf, lon_buf, parm_buf2, qc_buf2)) < 0)
            return FAIL;
        if ((closeHDF(sdfid2, fid2)) < 0)
            return FAIL;
        if ((parm_flag == 4) && (dims[1] != lon_dim_p4)) {
            /*
             *  in rare case where the 2 ozone arrays are different in lon size
             *  as between EP and AURA , make the parm_buf1 the same size as parm_buf2
             */
            printf(
                    "%s, %d - I: Detect different longitude size in ancilary data\n",
                    __FILE__, __LINE__);
            printf("   adapting file1\n");
            printf("     file1: %s\n", file1);
            printf("     file2: %s\n", file2);
            /*  resize array 1 to size of array 2 */
            resiz_anc((int16 *) parm_buf1, qc_buf1, dims[0], lon_dim_p4, dims[1],
                    lon_buf_p4, lon_buf);
        }
    }
    return SUCCEED;
}

/*-----------------------------------------------------------------------------
  Function: resiz_anc
  Returns:  None
  Description:
    resize a data and qc array in the longitude dimension to match that of
    another array.  Use simple lineas interpolation

      Type         Name        I/O     Description
      ----         ----        ---     -----------
      int16 *      data        I/O     data array to resize (size change:
                                nlat, nlon1 -> nlat, nlon2)
      int8 *       qc          I/O     qc array
      int32        nlat         I      # points in latitude
      int32        nlon1        I      # lon points in input array
      int32        nlon2        I      # lon points in final array
      float32 *    lons1        I      input longitude array
      float32 *    lons2        I      final longitude array

  W. Robinson, SAIC, 3 Dec 2013  Original development
  Note that this oonly applies to the ozone, which comes in int16 size
      
-----------------------------------------------------------------------------*/
void resiz_anc(int16_t *data, int8_t *qc, int32_t nlat, int32_t nlon1, int32_t nlon2,
        float * lons1, float * lons2) {
    int32 ix, ilon, ilat, intloc, ix1, ix2;
    float st1, del1, lonwt, wt1, wt2, lon, floc;
    int8 *cp_qc, q_fin, q1, q2;
    int16 *cp_data, dat1, dat2, dat_fin;
    /*
     *  we assume a ( -180 -> 180 ) system for both arrays
     */
    st1 = lons1[0];
    del1 = lons1[1] - st1;
    /*
     *  allocate a temporary storage for the data and qc and copy it to these 
     */
    cp_data = (int16 *) malloc(nlat * nlon1 * sizeof ( int16));
    cp_qc = (int8 *) malloc(nlat * nlon1 * sizeof ( int8));

    memcpy(cp_data, data, nlat * nlon1 * sizeof ( int16));
    memcpy(cp_qc, qc, nlat * nlon1 * sizeof ( int8));
    /*
     *  Note ix is the index of the values in data to interpolate (with ix + 1)
     *  when mapping to the final array.  lonwt is weight to apply
     */
    for (ilon = 0; ilon < nlon2; ilon++) {
        lon = lons2[ilon];
        floc = (lon - st1) / del1;
        intloc = (int) (floc + 1.) - 1; /* the '1' makes sure the lon is 
                              always ahead of what is pointed to by iloc */
        lonwt = floc - (float) intloc;
        intloc = (intloc < 0) ? (nlon1 - 1) : intloc;
        ix = intloc;
        /*
         *  do the interpolation for the row of latitudes at this longitude
         */
        ix1 = ix;
        ix2 = (ix1 == nlon1 - 1) ? 0 : ix + 1;
        wt1 = 1. - lonwt;
        wt2 = lonwt;
        for (ilat = 0; ilat < nlat; ilat++) {
            dat1 = *(cp_data + ix1 + nlon1 * ilat);
            dat2 = *(cp_data + ix2 + nlon1 * ilat);
            /*
             *  only the new qc of 20 will indicate bad qc
             */
            q1 = (*(cp_qc + ix1 + nlon1 * ilat) == 20) ? 1 : 0;
            q2 = (*(cp_qc + ix2 + nlon1 * ilat) == 20) ? 1 : 0;

            if ((q1 == 1) && (q2 == 1)) {
                dat_fin = dat1;
                q_fin = 20;
            } else if ((q1 == 0) && (q2 == 1)) {
                dat_fin = dat1;
                q_fin = 0;
            } else if ((q1 == 1) && (q2 == 0)) {
                dat_fin = dat2;
                q_fin = 0;
            } else {
                dat_fin = dat1 * wt1 + dat2 * wt2;
                q_fin = 0;
            };

            *(data + ilon + nlon2 * ilat) = dat_fin;
            *(qc + ilon + nlon2 * ilat) = q_fin;
        }
    }
    /*
     *  at the end, release workspace
     */
    free(cp_data);
    free(cp_qc);
}

/*-----------------------------------------------------------------------------
    Function: extract_data_pts

    Returns:  None

    Description:
        The function extract_data_pts extracts the four surrounding data
        points and the latitude/longitude points for the given coordinates
        and calls interpolation routine to interpolate the data points.

    Arguments: (in calling order)
      Type         Name        I/O     Description
      ----         ----        ---     -----------
      int16        PA_flag      I      product type: if 1, climatology else NRT
      int16        parm_flag    I      indicates wh parm it is being processed 
      timech   *   DT1          I      Time diff bet actual time & time step 1
      timech   *   DT2          I      Time diff bet actual time & time step 2
      int16        nsamp        I      number of scan line points
      float32 *    in_latlist   I      input latitude array 
      float32 *    in_lonlist   I      input longitude array
      float32 *    lat_bufp     I      latitude buffer
      float32 *    lon_bufp     I      longitude buffer
      int32 *      dims         I      dimensions of the data
      int8  *      qc_buf1      I      buffer containing QC flags from file1
      int8  *      qc_buf2      I      buffer containing QC flags from file2
      void  *      parm_buf1    I      buffer containing data from file1
      void  *      parm_buf2    I      buffer containing data from file2
      intn         toms;        I      BOOLEAN value whether EP TOMS data are used
      float32 *    out_lat_list O      surrounding latitude points
      float32 *    out_lon_list O      surrounding longitude points
      int16 *      qcflag       O      array corresponding to interp values,
                                        containing flags indicating the 
                                        integrity of those values.
      float32 *    interp       O      interpolated values of the data 
      float32 *    anc_unc      O      uncertainty in the ancillary value

    Modification history:
      Lakshmi Kumar  Hughes STX      04/02/93  Original development

      Lakshmi Kumar  Hughes STX      12/09/93  Modified to call this rtn only
                                                once for the given list of 
                                                lat/lons to avoid function call
                                                overhead
      Lakshmi Kumar  Hughes STX	     11/27/95  Changed the order of data pts
                                                in data_list1 and data_list2
      Lakshmi Kumar  Hughes STX      07/02/97  fixed code to extract correct 
                                                4 corner pts when input lat/lon 
                                                lies outside the lat/lon buffs 
      Ewa Kwiatkowska  GSC	     11/12/99  Added parts supporting TOMS 
                                                ancillary specification
      W Robinson, SAIC              4 Aug 2016   Pass through the ancillary 
                                                 uncertainty: anc_unc
-----------------------------------------------------------------------------*/
void extract_data_pts(int16_t PA_flag, int16_t parm_flag, timech DTime1, timech DTime2,
        int16_t nsamp, float *in_latlist, float *in_lonlist, float *lat_bufp,
        float *lon_bufp, int32_t *dims, int8_t *qc_buf1, int8_t *qc_buf2,
        void *parm_buf1, void *parm_buf2, int toms, float *out_lat_list,
        float *out_lon_list, int16_t *qcflag, float *interp, float *anc_unc) {

    int8 qc1_vals[4];
    int8 qc2_vals[4];
    int16 i, p, lat_index1 = -1, lat_index2 = -1;
    int16 lon_index1 = -1, lon_index2 = -1, index;
    int16 *OZ_list1, *OZ_list2, outside = 0;
    int16 *oz_p1, *oz_p2, loop;
    int32 int_qc;
    float32 in_lat, in_lon, *wph_p1, *wph_p2, data_list1[4], data_list2[4];
    float32 *WPH_list1, *WPH_list2, gridstep;
    float64 DT1, DT2;


    OZ_list1 = (int16 *) data_list1;
    OZ_list2 = (int16 *) data_list2;
    oz_p1 = (int16 *) parm_buf1;
    oz_p2 = (int16 *) parm_buf2;

    WPH_list1 = (float32 *) data_list1;
    WPH_list2 = (float32 *) data_list2;
    wph_p1 = (float32 *) parm_buf1;
    wph_p2 = (float32 *) parm_buf2;

    gridstep = lon_bufp[1] - lon_bufp[0];


    for (loop = 0; loop < nsamp; loop++) {

        in_lat = in_latlist[loop];
        in_lon = in_lonlist[loop];
        for (i = 0; i < 4; i++)
            qc1_vals[i] = qc2_vals[i] = -1;

        lat_index1 = lat_index2 = -1;
#if 1
        BinarySearch(dims[0], lat_bufp, in_lat, BS_DECR, &lat_index1,
                &lat_index2);
#else
        for (done = 0, i = 0; i < dims[0] && done == 0; i++) {
            if (in_lat == lat_bufp[i]) {
                lat_index1 = i;
                if (i + 1 == dims[0])
                    lat_index2 = i - 1;
                else
                    lat_index2 = i + 1;
                done = 1;
            }
            else {
                if (in_lat > lat_bufp[i]) {
                    lat_index1 = i - 1;
                    lat_index2 = i;
                    done = 1;
                }
            }
        }

#endif
        if (lat_index1 == -1 && lat_index2 == -1) {
            if (in_lat > 0)
                lat_index1 = lat_index2 = 0;
            else
                lat_index1 = lat_index2 = dims[0] - 1;

            /* lat_index2 = dims[0] - 2; */
        } else
            if (lat_index1 == -1)
            lat_index1 = 0;

        lon_index1 = lon_index2 = -1;

#if 1
        BinarySearch(dims[1], lon_bufp, in_lon, BS_INCR, &lon_index1,
                &lon_index2);
#else
        for (done = 0, i = 0; i < dims[1] && done == 0; i++) {
            if (in_lon == lon_bufp[i]) {
                lon_index1 = i;
                if (i + 1 == dims[1])
                    lon_index2 = i - 1;
                else
                    lon_index2 = i + 1;
                done = 1;
            }
            else {
                if (in_lon < lon_bufp[i]) {
                    lon_index1 = i - 1;
                    lon_index2 = i;
                    done = 1;
                }
            }
        }

#endif
        if (lon_index1 == -1 || lon_index2 == -1) { /* global  */

            if (((lon_bufp[0] + 360 - gridstep) == lon_bufp[dims[1] - 1]) || (toms == SUCCEED)) {
                lon_index1 = dims[1] - 1;
                lon_index2 = 0;
                outside = 1;
            } else {
                if (in_lon < lon_bufp[0]) { /* regional */
                    lon_index1 = 0;
                    lon_index2 = 1;
                } else {
                    lon_index1 = dims[1] - 1;
                    lon_index2 = dims[1] - 2;
                }
            }
        }

        out_lat_list[3] = lat_bufp[lat_index1];
        out_lat_list[0] = lat_bufp[lat_index1];
        out_lat_list[1] = lat_bufp[lat_index2];
        out_lat_list[2] = lat_bufp[lat_index2];

        out_lon_list[3] = lon_bufp[lon_index1];
        out_lon_list[0] = lon_bufp[lon_index2];
        out_lon_list[1] = lon_bufp[lon_index2];
        out_lon_list[2] = lon_bufp[lon_index1];

        if (parm_flag == OZONE) {

            if ((outside) && (toms == SUCCEED)) {

                if (in_lon < 0)
                    index = lon_index2;
                else
                    index = lon_index1;

                OZ_list1[3] = oz_p1[lat_index1 * dims[1] + index];
                OZ_list1[0] = oz_p1[lat_index1 * dims[1] + index];
                OZ_list1[1] = oz_p1[lat_index2 * dims[1] + index];
                OZ_list1[2] = oz_p1[lat_index2 * dims[1] + index];

                OZ_list2[3] = oz_p2[lat_index1 * dims[1] + index];
                OZ_list2[0] = oz_p2[lat_index1 * dims[1] + index];
                OZ_list2[1] = oz_p2[lat_index2 * dims[1] + index];
                OZ_list2[2] = oz_p2[lat_index2 * dims[1] + index];
            } else {

                OZ_list1[3] = oz_p1[lat_index1 * dims[1] + lon_index1];
                OZ_list1[0] = oz_p1[lat_index1 * dims[1] + lon_index2];
                OZ_list1[1] = oz_p1[lat_index2 * dims[1] + lon_index2];
                OZ_list1[2] = oz_p1[lat_index2 * dims[1] + lon_index1];

                if (PA_flag == 1) {
                    for (i = 0; i < CORNERS; i++)
                        OZ_list2[i] = 0;
                } else {
                    OZ_list2[3] = oz_p2[lat_index1 * dims[1] + lon_index1];
                    OZ_list2[0] = oz_p2[lat_index1 * dims[1] + lon_index2];
                    OZ_list2[1] = oz_p2[lat_index2 * dims[1] + lon_index2];
                    OZ_list2[2] = oz_p2[lat_index2 * dims[1] + lon_index1];
                }
            }
        } else {

            WPH_list1[3] = wph_p1[lat_index1 * dims[1] + lon_index1];
            WPH_list1[0] = wph_p1[lat_index1 * dims[1] + lon_index2];
            WPH_list1[1] = wph_p1[lat_index2 * dims[1] + lon_index2];
            WPH_list1[2] = wph_p1[lat_index2 * dims[1] + lon_index1];

            if (PA_flag == 1) {
                for (i = 0; i < CORNERS; i++)
                    WPH_list2[i] = 0;
            } else {
                WPH_list2[3] = wph_p2[lat_index1 * dims[1] + lon_index1];
                WPH_list2[0] = wph_p2[lat_index1 * dims[1] + lon_index2];
                WPH_list2[1] = wph_p2[lat_index2 * dims[1] + lon_index2];
                WPH_list2[2] = wph_p2[lat_index2 * dims[1] + lon_index1];
            }
        }

        if (PA_flag != 1) {

            if ((outside) && (toms == SUCCEED)) {

                qc1_vals[3] = qc_buf1[lat_index1 * dims[1] + index];
                qc1_vals[0] = qc_buf1[lat_index1 * dims[1] + index];
                qc1_vals[1] = qc_buf1[lat_index2 * dims[1] + index];
                qc1_vals[2] = qc_buf1[lat_index2 * dims[1] + index];

                qc2_vals[3] = qc_buf2[lat_index1 * dims[1] + index];
                qc2_vals[0] = qc_buf2[lat_index1 * dims[1] + index];
                qc2_vals[1] = qc_buf2[lat_index2 * dims[1] + index];
                qc2_vals[2] = qc_buf2[lat_index2 * dims[1] + index];
            } else {

                qc1_vals[3] = qc_buf1[lat_index1 * dims[1] + lon_index1];
                qc1_vals[0] = qc_buf1[lat_index1 * dims[1] + lon_index2];
                qc1_vals[1] = qc_buf1[lat_index2 * dims[1] + lon_index2];
                qc1_vals[2] = qc_buf1[lat_index2 * dims[1] + lon_index1];

                qc2_vals[3] = qc_buf2[lat_index1 * dims[1] + lon_index1];
                qc2_vals[0] = qc_buf2[lat_index1 * dims[1] + lon_index2];
                qc2_vals[1] = qc_buf2[lat_index2 * dims[1] + lon_index2];
                qc2_vals[2] = qc_buf2[lat_index2 * dims[1] + lon_index1];
            }

            //         qc_check = 0;
            //         for (i = 0; i < 4; i++) {
            //            if(qc1_vals[i] != 0 || qc2_vals[i] != 0)
            //              qc_check = 1;
            //          }
        }

        /*
        if (PA_flag != 1 && qc_check == 1)
           qcflag[loop] = -1;
        else 
           qcflag[loop] = 0; 
        above code is replaced by the int_qc from interpolate below */


        if (outside) {
            if (in_lon < 0) {
                for (p = 0; p < 4; p++)
                    if (out_lon_list[p] > 0)
                        out_lon_list[p] -= 360;
            } else {
                for (p = 0; p < 4; p++)
                    if (out_lon_list[p] < 0)
                        out_lon_list[p] += 360;
            }
        }


        DT1 = DTime1.start + DTime1.inc*in_lon;
        if ((toms == SUCCEED)&&(DT1 < 0.0)) DT1 *= -1.0;
        DT2 = DTime2.start + DTime2.inc*in_lon;
        if ((toms == SUCCEED)&&(DT2 < 0.0)) DT2 *= -1.0;


        outside = 0;

        interpolate(PA_flag, parm_flag, DT1, DT2, in_lat, in_lon,
                out_lat_list, out_lon_list, data_list1, data_list2,
                qc1_vals, qc2_vals, &interp[loop], &anc_unc[loop],
                &int_qc);
        qcflag[loop] = int_qc;

        /*
         */
    }
}

/*-----------------------------------------------------------------------------
    Function:  gregor

    Returns:   none

    Description:
        Converts Julian day into Gregorian Month/Day.

    Parameters:
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int16     sday         I    Julian day of syear   
        int16     syear        I    Year in which sday occurs
        int16	  day	       O    Day of Month in Year
        int16     month        O    Month in Year


    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Michael Darzi  GSC            10/89     Original development in FORTRAN
    Lakshmi Kumar  Hughes STX     04/02/93  Converted to C
    Lakshmi Kumar  Hughes STX     12/09/93  Modified the comments
-----------------------------------------------------------------------------*/
void gregor(int16_t sday, int16_t syear, int16_t *month, int16_t *day) {
    int32 sdays[2][12] = {
        {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334},
        {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335}
    };
    int16 row = 0, mth = 0, done = 0;

    if (syear % 4 == 0)
        row = 1;

    mth = 11;
    while (!done) {
        if (sday > sdays[row][mth]) {
            *day = sday - sdays[row][mth];
            done = 1;
        } else
            mth = mth - 1;
    }
    *month = mth + 1;
}

/*----------------------------------------------------------------------------
    Function: ck_files_in_buf 

    Returns:  Status

    Description:
      The function ck_files_in_buf checks whether the requested real time 
        data or the climatology data are in the buffers or not.  Sets the
        readflag to 1 if file1 needs to be read, 2 if file2 needs to be read
        and 3 if both file1 and file2 needs to be read.

    Arguments: (in calling order)
      Type         Name        I/O     Description
      ----         ----        ---     -----------
      int16        PA_flag      I      product type: 1 for climatology, else
                                        NRT 
      int16        parm_flag    I      flag indicating the parameter
      char *       f1           I      file 1     
      char *       f2           I      file 2       
      int16        month        I      data month
      int16 *      error_flag  I/O     error_flag indicating previous read
                                        errors
      int16 *      read_flag    O      flag indicating read:
                                        0 - no read, 1 - read file1,
                                        2 - read file2, 3 -read file1 and file2	
    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     04/02/93  Original development
    Lakshmi Kumar  Hughes STX     12/09/93  Moved the climatology data check
                                            for the given month from read-
                                            climatology fnc. to this fnc.
    Lakshmi Kumar  Hughes STX     04/18/96  Fixed a potential array out of
                                            bounds problem (ref. to error_flag)
    Lakshmi Kumar  Hughes STX     08/16/96  Changed const 5 to NPARMS in the
                                            if condn. below
-----------------------------------------------------------------------------*/
int32_t ck_files_in_buf(int16_t PA_flag, int16_t parm_flag, char *f1, char *f2,
        int16_t month, int16_t *error_flag, int16_t *read_flag) {
    int16 i;
    /* prev parameter files read for -t step */
    typedef char p_files_t[MAXVAL];
    static p_files_t *p_files1;
    /* prev parameter files read for +t step */
    static p_files_t *p_files2;
    static int16 *prev_mth, fst_call = 1;
    char *p_f1, *p_f2; /* pointers to prev file1 and file2 */
    void *temp_p;

    if (PA_flag == 1) {
        if (strlen(f1) <= 0)
            return FAIL;
    } else {
        if (strlen(f1) <= 0 || strlen(f2) <= 0)
            return FAIL;
    }

    if (fst_call) {
        p_files1 = (p_files_t*) allocateMemory(NPARMS * sizeof (p_files_t), "p_files1");
        p_files2 = (p_files_t*) allocateMemory(NPARMS * sizeof (p_files_t), "p_files2");
        prev_mth = (int16*) allocateMemory(NPARMS * sizeof (int16), "prev_mth");
        for (i = 0; i < NPARMS; i++)
            prev_mth[i] = -1;
        fst_call = 0;
    }

    if (*error_flag >= 0 && *error_flag < NPARMS) {
        strcpy(p_files1[*error_flag], "");
        strcpy(p_files2[*error_flag], "");
        *error_flag = -1;
    }

    p_f1 = p_files1[parm_flag];
    p_f2 = p_files2[parm_flag];

    *read_flag = 0;

    if (PA_flag == 1) {
        if ((strcmp(f1, p_f1)) != 0) {
            *read_flag = 1;
            strcpy(p_f1, f1);
        } else
            if (prev_mth[parm_flag] != month) { /* ck if data in buf is for*/
            *read_flag = 1; /* the requested month, if */
        }
        prev_mth[parm_flag] = month; /* not, set the read_flag  */
    } else {
        if (((strcmp(f1, p_f1)) == 0) && (strcmp(f2, p_f2)) == 0)
            *read_flag = 0;
        else {
            /*
             *  for the ozone, force a read of both files
             */
            if (parm_flag == 4) {
                *read_flag = 3;
                strcpy(p_f1, f1);
                strcpy(p_f2, f2);
            } else {
                if (strcmp(f1, p_f2) == 0) {
                    *read_flag = 2;
                    temp_p = parm_buf1; /* switch pointers   */
                    parm_buf1 = parm_buf2;
                    parm_buf2 = temp_p;
                    strcpy(p_f1, f2);
                } else {
                    if (strcmp(f2, p_f1) == 0) {
                        *read_flag = 1;
                        temp_p = parm_buf1;
                        parm_buf1 = parm_buf2;
                        parm_buf2 = temp_p;
                        strcpy(p_f2, f1);
                    } else {
                        if (strcmp(f2, p_f2) == 0)
                            *read_flag = 1;
                        else
                            *read_flag = 3;
                        strcpy(p_f1, f1);
                        strcpy(p_f2, f2);
                    }
                }
            }
        }
    }

    return SUCCEED;
}

/*----------------------------------------------------------------------------
    Function:  interpolate  

    Returns:   None

    Description:
        Initializes and sets some input and output parameters and calls
        FORTRAN subroutine "dataintp" to do the interpolation.

    Parameters: 
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int16	  PA_flag      I    Flag set to 1 to read climatology else
                                    to read NRT data files
        int16     parm_flag    I    indicates for wh parmeter it is being 
                                    processed 
        float32   DT1          I    Time diff bet file1 and file2 
        float32   DT2          I    Time diff bet file1 and file2
        float32   in_lat       I    input latitude
        float32   in_lon       I    input longitude
        float32 * lat_list     I    adjacent latitude pts of the given latitude 
        float32 * lon_list     I    adjacent lon pts of the given longitude
        void    * data_p1      I    adjacent data pts of file 1 
        void    * data_p2      I    adjacent data pts of file 2
        int8 *    qc1, qc2     I    the qc values for each data_p1, data_p2
        float32 * intpdata     O    interpolated data for given lat/lon   
        float32 * anc_unc      O    ancillary data uncertainty
        int32 *   int_qc       O    qc for interpolation 0 - good, 1 - 
                                    insufficient good data values to interpolate

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     03/22/93  Original development
    Lakshmi Kumar  Hughes STX     12/09/93  Modified comments
    Lakshmi Kumar  Hughes STX     04/18/96  Removed a redundant for loop
    Lakshmi Kumar  Hughes STX     03/28/97  Defined range for humidity if parm
                                            type is 5.
    W. Robinson, SAIC, 13 Dec 2013          include the qc values for the points
    W. Robinson, SAIC, 5 Jun 2014           expand valid RH range down to 0%
                                            and set change default to 80% 
                                            instead of 90%
    W. Robinson, SAIC, 4 Aug 2016           pass through an uncertainty: anc_unc
-----------------------------------------------------------------------------*/
void interpolate(int16_t PA_flag, int16_t parm_flag, double DT1, double DT2,
        float in_lat, float in_lon, float *lat_list,
        float *lon_list, void *data_p1, void *data_p2,
        int8_t *qc1, int8_t *qc2, float *intpdata,
        float *anc_unc, int32_t *int_qc) {

    float32 in_latlon[2];
    float32 range[2];
    float32 def = 0.0;
    float32 dummy[4];
    float32 *WPH_p1, *WPH_p2;
    float32 data_list1[4], data_list2[4];
    float32 dataout = -1.343, unc = 0.;
    int32 ipt = 4;
    int32 nband = 1;
    int32 intporder;
    int32 row = 1, col = 1, i, int_bad;
    int16 *OZ_p1, *OZ_p2;


    in_latlon[0] = in_lat;
    in_latlon[1] = in_lon;

    range[0] = -98.00;
    range[1] = 9999.99;

    switch (parm_flag) {
    case 0: /* wind-u and wind-v */
    case 1:
        def = 6;
        range[0] = -40;
        range[1] = 40;
        break;
    case 2: /* pressure */
        def = 1013;
        range[0] = 850;
        range[1] = 1084;
        break;
    case 3: /* precipitable water */
        def = 50;
        range[0] = 0;
        range[1] = 200;
        break;
    case 4: /* ozone */
        def = 360;
        range[0] = 80;
        range[1] = 600;
        break;
    case 5: /* relative humidity */
        def = 80;
        range[0] = 0;
        range[1] = 100;
        break;
    default:
        break;
    }

    intporder = SPATIAL_TEMPORAL;

    if (PA_flag == 1) /* climatology requested */
        DT1 = DT2 = 0;

    if (DT1 == 0 && DT2 == 0) /* if climatology or same NRT files */
        intporder = SPATIAL;
    /*
     *  we make sure that any data value that is qc flagged with 20 (= bad data)
     *  will have a value outside the good data range
     */
    if (parm_flag == OZONE) {
        OZ_p1 = (int16 *) data_p1;
        OZ_p2 = (int16 *) data_p2;
        for (i = 0; i < 4; i++) {
            data_list1[i] = (qc1[i] == 20) ? range[0] - 1 : OZ_p1[i];
            data_list2[i] = (qc2[i] == 20) ? range[0] - 1 : OZ_p2[i];
        }
    } else {
        WPH_p1 = (float32 *) data_p1;
        WPH_p2 = (float32 *) data_p2;
        for (i = 0; i < 4; i++) {
            data_list1[i] = (qc1[i] == 20) ? range[0] - 1 : WPH_p1[i];
            data_list2[i] = (qc2[i] == 20) ? range[0] - 1 : WPH_p2[i];
        }
    }
    // dataintp_(in_latlon, lat_list, lon_list, data_list1, &DT1, data_list2,
    //         &DT2, &ipt, &nband, &range, &def, &intporder, dummy,
    //         &dataout, &unc, &int_bad, &row, &col);
    dataintp(in_latlon, lat_list, lon_list, data_list1, &DT1, data_list2,
            &DT2, &ipt, &nband, range, &def, &intporder, dummy,
            &dataout, &unc, &int_bad, &row, &col);

    *intpdata = dataout;
    *anc_unc = unc;
    *int_qc = int_bad;

}

/*----------------------------------------------------------------------------
    Function:  set_files    

    Returns:   Status

    Description:
        Set_files sets file1 and file2 to input filename1 and filename2 to 
        generate the interp values if the scan-line's date and time, as
        determined by syear, sday, eday, and msec, fall between the time range
        represented by these two files; if the scan's date and time fall 
        outside their range, it will use filename2 and filname3 for generating
        interp; ignored if time fall within the range of filename1 and 
        filename2

    Parameters: 
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int16     syear        I    year of data start time
        int16     sday	       I    day-of-year for data start time
        int16     eday         I    day-of-year for data end time
        int32     msec         I    milliseconds-of-day for data start time
        char *    filename1    I    input filename1
        char *    filename2    I    input filename2
        char *    filename3    I    input filename3
        char *    file1        O    file to be used for generating interp vals 
        char *    file2        O    file to be used for generating interp vals 
        timech  * dtime1       O    Time diff bet actual scan time & time step 1    
        timech  * dtime2       O    Time diff bet actual scan time & time step 2
        intn    * toms         O    Boolean value whether EP TOMS approximation is used

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     06/21/94  Original development
    Lakshmi Kumar  HITC		  05/25/95  Added code to return error messages
    Lakshmi Kumar  Hughes STX	  10/10/95  Changed time checking logic
                                            (Ref: V4.4 I/O Specs.)
    Lakshmi Kumar  Hughes STX	  03/08/96  changed get_time interface to
                                            return files start and end times
                                            in julian days.  section of the
                                            code that was converting time to
                                            julian days is moved to get_time 
    Ewa Kwiatkowska  GSC	  11/12/99  Added parts supporting TOMS 
                                            ancillary specification 
-----------------------------------------------------------------------------*/

int32_t set_files(int16_t parm_flag, int16_t syear, int16_t sday, int16_t eday,
        int32_t msec, char *filename1, char *filename2, char *filename3,
        char *file1, char *file2, timech *dtime1, timech *dtime2, int *toms) {
    int16 dday, dyear, month, day;
    int16 hh, mm;
    int32 ss;
    float64 d_jd, s_jd1, e_jd1, s_jd2, e_jd2, s_jd3, e_jd3;
    float64 dtin[2];
    div_t quot1, quot2;
    char *FUNC = "set_files";
    static int32 prev_msec = -1;
    float64 dt1, dt2;

    if (prev_msec == -1)
        prev_msec = msec;
    /*** get start and end times of all 3 files, calculate center point time */

    if ((anc_get_time(filename1, &s_jd1, &e_jd1)) < 0)
        return FAIL;

    if (strcmp(filename2, filename1) == 0) {
        s_jd2 = s_jd1;
        e_jd2 = e_jd1;
    } else {
        if ((anc_get_time(filename2, &s_jd2, &e_jd2)) < 0)
            return FAIL;
    }

    if (strcmp(filename3, filename2) == 0) {
        s_jd3 = s_jd2;
        e_jd3 = e_jd2;
    } else {
        if ((anc_get_time(filename3, &s_jd3, &e_jd3)) < 0)
            return FAIL;
    }

    /*** convert given scan time to julian day */
    dyear = syear;
    dday = sday;
    if (sday != eday && msec < 43200000)
        dday = eday;
    if (dday < sday)
        dyear += 1;

    /**** call gregor to calculate month and day from julian day */
    gregor(dday, dyear, &month, &day);
    dtin[0] = (syear * 100 + month)*100 + day;

    quot1 = div(msec, MSECHOUR);
    hh = quot1.quot;
    quot2 = div(quot1.rem, MSECMIN);
    mm = quot2.quot;
    ss = quot2.rem;
    dtin[1] = (hh * 100 + mm)*100 + (ss / 1000.0);

    julian_(dtin, &d_jd);

    /**** check if the given time is in error */

    /****set file1 and file2 depending upon the given input date and time */

    if (s_jd1 == s_jd2 && e_jd1 == e_jd2) {
        strcpy(file1, filename2);
        strcpy(file2, filename2);
        dt1 = (d_jd - s_jd2);
        dt2 = (d_jd - s_jd2);
    } else {
        if (d_jd < s_jd1) {
            sprintf(ERR_MSG, "%s: Input time error: \nthe given time %f "
                    "(in jdays) is less than \n"
                    "%s start time of %f", FUNC, d_jd, filename1, s_jd1);
            return FAIL;
        }
        if (d_jd <= e_jd2) { /* data time is bet file1 & file2 times */
            strcpy(file1, filename1);
            strcpy(file2, filename2);
            dt1 = (d_jd - s_jd1);
            dt2 = (d_jd - s_jd2);
        } else { /* if (d_jd > s_jd2) */
            strcpy(file1, filename2);
            dt1 = (d_jd - s_jd2);
        }
    }

    if (s_jd2 == s_jd3 && e_jd2 == e_jd3) {
        strcpy(file2, filename2);
        dt2 = (d_jd - s_jd2);
    } else {
        if (d_jd > e_jd3) { /* data time falls outside file 3 ? */
            sprintf(ERR_MSG, "%s: Input time error: \nthe given time %f "
                    "(in jdays) is greater than %s \n"
                    "end time of %f", FUNC, d_jd, filename3, e_jd3);
            return FAIL;
        }
        if (d_jd >= s_jd2) { /* data time is bet file2 and file3 times? */
            strcpy(file2, filename3);
            dt2 = (d_jd - s_jd3);
        }
    }


    *toms = check_on_TOMS(parm_flag, file1, file2, filename1, filename2,
            filename3, s_jd1, s_jd2, s_jd3, e_jd1, e_jd2, e_jd3, d_jd, dtime1, dtime2);

    if (*toms == FAIL) {
        dtime1->start = dt1;
        if (dtime1->start < 0) dtime1->start *= -1.0;
        dtime2->start = dt2;
        if (dtime2->start < 0) dtime2->start *= -1.0;
        dtime1->inc = 0.0;
        dtime2->inc = 0.0;
    }

    if ((strcmp(file1, file2)) == 0)
        dtime1->start = dtime2->start = dtime1->inc = dtime2->inc = 0;


    return SUCCEED;
}

/*----------------------------------------------------------------------------
    Function:  get_time     

    Returns:   int32 status
               The return code is a negative value if any error occurs.	

    Description:
        get_time checks ftime structure to see if the given input file
        name and times exists in that structure.  If so, it reads the
        start and end times from that structure and returns the values. 
        Otherwise, it opens the given HDF file, reads global attributes 
        syear, sday, smsec, eyear, eday, and emsec.  Calls julian_ to
        convert start and end times to julian days and sets these values
        in the ftime structure, and returns.

    Parameters: 
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        char *    filename     I    file name
        float64	 *s_jd 	       O    start time in julian days
        float64  *e_jd	       O    end time in julian days

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     06/21/94  Original development
    Lakshmi Kumar  HITC		  05/22/95  Modified to access date and time
                                            from global attributes rather than
                                            reading it from the file name
    Lakshmi Kumar  Hughes STX	  10/10/95  Modified to read start and end
                                            date and time from the input file
                                            and return them in the format shown
                                            above. 
    Lakshmi Kumar  Hughes STX	  03/08/96  interface has been changed and
                                            added code to convert file times
                                            to julian days.
-----------------------------------------------------------------------------*/

int32_t anc_get_time(char *filename, double *s_jd, double *e_jd) {

    typedef struct file_time {
        char *fn;
        float64 s_jd;
        float64 e_jd;
    } ftm;
    static ftm *ftime;

    int32 i, ret, done, fid, sdfid;
    int16 syear, eyear, sday, eday, month, day;
    int32 smsec, emsec;
    int32 hh, mm, ss;
    float64 sdate, stime, edate, etime;
    float64 stin[2], etin[2];
    div_t quot1, quot2;
    static int32 file_num = 0;


    if (ftime == NULL) {
        ftime = (ftm*) allocateMemory(FLIMIT * sizeof (ftm), "ftime");
    }

    done = 0;
    for (i = 0; !done && i < FLIMIT && ftime[i].fn != NULL; i++) {
        if ((strcmp(ftime[i].fn, filename)) == 0) {
            *s_jd = ftime[i].s_jd;
            *e_jd = ftime[i].e_jd;
            done = 1;
        }
    }

    if (done)
        return SUCCEED;
    else if (i >= FLIMIT) {
        i = file_num % FLIMIT;
        free(ftime[i].fn);
    }

    ftime[i].fn = (char *) malloc(sizeof (char) * (strlen(filename) + 1));
    strcpy(ftime[i].fn, filename);

    /**** Open the given data file */
    if ((openHDF(filename, &sdfid, &fid)) < 0) return FAIL;
    file_num++;

    /**** read global attribute "Start Year" */
    ret = rdancattr(sdfid, SYEAR, (VOIDP *) & syear);
    if (ret < 0)
        return FAIL;
    if (DFNT_INT16 != ret) {
        fprintf(stderr, "\nWARNING: Datatype of %s read from %s is not int16. Processing from this point is unpredictable", SYEAR, filename);
    }

    /**** read global attribute "Start Day" */
    if ((ret = rdancattr(sdfid, SDAY, (VOIDP *) & sday)) < 0) {
        printf("Execution halted %s:%d: 'Start Day' global attribute not found in %s\n",
                    __FILE__, __LINE__, filename);
        exit(EXIT_FAILURE);
    }
    if (DFNT_INT16 != ret) {
        fprintf(stderr, "\nWARNING: Datatype of %s read from %s is not int16.  Processing from this point is unpredictable", SDAY, filename);
    }

    /**** read global attribute "Start Millisec" */
    if ((ret = rdancattr(sdfid, SMSEC, (VOIDP *) & smsec)) < 0) return FAIL;
    if (ret != DFNT_INT32) {
        fprintf(stderr, "\nWARNING: Datatype of %s read from %s is not int32.  Processing from this point is unpredictable", SMSEC, filename);
    }

    /**** read global attribute "End Year" */
    if ((ret = rdancattr(sdfid, EYEAR, (VOIDP *) & eyear)) < 0) return FAIL;
    if (ret != DFNT_INT16) {
        fprintf(stderr, "\nWARNING: Datatype of %s read from %s is not int16.  Processing from this point is unpredictable", EYEAR, filename);
    }

    /**** read global attribute "End Day" */
    if ((ret = rdancattr(sdfid, EDAY, (VOIDP *) & eday)) < 0) return FAIL;
    if (ret != DFNT_INT16) {
        fprintf(stderr, "\nWARNING: Datatype of %s read from %s is not int16.  Processing from this point is unpredictable", EDAY, filename);
    }

    /**** read global attribute "End Millisec" */
    if ((ret = rdancattr(sdfid, EMSEC, (VOIDP *) & emsec)) < 0) return FAIL;
    if (ret != DFNT_INT32) {
        fprintf(stderr, "\nWARNING: Datatype of %s read from %s is not int32.  Processing from this point is unpredictable", EMSEC, filename);
    }

    /**** call gregor to calculate month and day from julian day */
    gregor(sday, syear, &month, &day);
    sdate = (syear * 100 + month)*100 + day;

    gregor(eday, eyear, &month, &day);
    edate = ((eyear * 100 + month)*100 + day);


    /**** extract hours from msec */
    quot1 = div(smsec, MSECHOUR);
    hh = quot1.quot;
    quot2 = div(quot1.rem, MSECMIN);
    mm = quot2.quot;
    ss = quot2.rem;

    stime = ((hh * 100 + mm)*100 + (ss / 1000.0));

    quot1 = div(emsec, MSECHOUR);
    hh = quot1.quot;
    quot2 = div(quot1.rem, MSECMIN);
    mm = quot2.quot;
    ss = quot2.rem;

    etime = ((hh * 100 + mm)*100 + (ss / 1000.0));

    stin[0] = sdate;
    stin[1] = stime;
    etin[0] = edate;
    etin[1] = etime;
    julian_(stin, s_jd);
    julian_(etin, e_jd);

    ftime[i].s_jd = *s_jd;
    ftime[i].e_jd = *e_jd;

    /**** close input file */
    closeHDF(sdfid, fid);

    return SUCCEED;
}

