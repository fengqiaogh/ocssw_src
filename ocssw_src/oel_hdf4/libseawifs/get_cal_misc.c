/*-----------------------------------------------------------------------------
    File : get_cal_misc.c

    Contents:
        get_ref_time	-  reads reference time
        get_tindex	-  reads time vdata and returns appropriate index to 
                           access data
        read_parm_data  -  reads parameter data
        calc_knees	-  calculates knee1,2,3 and 4 counts and radiances
        setup_scanmod   -  sets up the scan-modulation correction factor array
                           for GAC data
        attach_vdata  	-  attaches to the requested vdata
        rdvdata		-  reads requested data from the given vdata
        sort_srad       -  sorts saturated radiances and returns ordered 
                           indices
        read_SDS        -  reads the requested SDS

    Other relevant files:
        cal.h		-  various #defined constants, TDI table, and also
                                includes hdf.h
        getcal_proto.h  -  prototypes for get_cal functions
        get_cal.c       -  a higher layer of calibration input functions

    Notes:

        Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      03/11/94    Original development
        Lakshmi Kumar    Hughes STX      06/07/94    Updated to reflect v3.1
                                                        interface spcifications
        Lakshmi Kumar    Hughes STX	 03/21/96    Corrected non-prototype
                                                     declarations
        Lakshmi Kumard Hughes STX	 03/17/97    Removed non-ANSI proto
                                                     declarations.  In-detector
                                                     offsets are redefined as
                                                     idoffs[8][16]. 
------------------------------------------------------------------------------*/

#include "get_cal.h"
#include "getcal_proto.h"
#include "hdf4utils.h"

#include <hdf.h>
#include <mfhdf.h>


/*-----------------------------------------------------------------------------
    Function: get_tindex 

    Returns: int32 (status)
        Returns status

    Description:
        The function get_ref_time reads reference date and time from the 
        input calibration table. 

    Arguments: (in calling order)
      Type       Name        I/O     Description
      ----       ----        ---     -----------
      int32      sdfid        I      SD file ID
      int16      ref_year     I      Reference Year
      int16      ref_day      I      Reference Day
      int16      ref_min      I      Reference Minute

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      05/22/96    Original development
------------------------------------------------------------------------------*/
int32_t get_ref_time(int32_t sdfid, int16_t *ref_year, int16_t *ref_day, int16_t *ref_min) {
    int32_t attrnum;

    attrnum = SDfindattr(sdfid, REFYEAR);
    if ((SDreadattr(sdfid, attrnum, ref_year)) < 0)
        return FAIL;

    attrnum = SDfindattr(sdfid, REFDAY);
    if ((SDreadattr(sdfid, attrnum, ref_day)) < 0)
        return FAIL;

    attrnum = SDfindattr(sdfid, REFMIN);
    if ((SDreadattr(sdfid, attrnum, ref_min)) < 0)
        return FAIL;

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function: get_tindex 

    Returns: int32 (index)
        On successful it returns the index of the given time entry and if
        time entry not found, returns -3.

    Description:
        The function get_tindex reads time vdata and searches for the
        given time entry.  If the given time found, it rerurns the entry
        number to access related information from slopes and parameter vdatas.

    Arguments: (in calling order)
      Type       Name        I/O     Description
      ----       ----        ---     -----------
      int32      fid          I      file ID
      int16      syear        I      year of data start time 
      int16      sday         I      day-of-year for data start time 
      int16      eday         I      day-of-year for data end time
      int32      smsec        I      milliseconds-of-day for data start time
      int16      *cal_year    O      year the cal entry was made
      int16      *cal_day     O      day of the year the cal entry was made

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      03/11/94    Original development
        Lakshmi Kumar    Hughes STX      06/07/94    Updated to reflect v3.1
                                                        interface spcifications
        Lakshmi Kumar	 Hughes STX	 02/07/94    Added code to return 
                                                     cal entry year and day
                                                     (ref to I/O specs v4.2)
         Gene Eplee       SAIC GSC        05/11/98     Fix for eyear-syear = 1
        W. Robinson,      SAIC           10/16/03   re-cast the index computation
                                                    to fix a time problem
------------------------------------------------------------------------------*/

int32_t get_tindex(int32_t fid, int16_t syear, int16_t sday, int16_t eday, int32_t msec,
        int16_t *cal_year, int16_t *cal_day) {

    int16 dyear, dday, *cal_syear, *cal_sday, *cal_eyear, *cal_eday;
    int16 *entry_year, *entry_day, ahead;
    int32_t i, *cal_smsec, *cal_emsec, vsid, elts;

    if ((vsid = attach_vdata(fid, TIME)) < 0)
        return RDERR;

    if ((elts = VSelts(vsid)) < 0)
        return RDERR;

    if ((cal_syear = (int16 *) malloc(elts * sizeof (int16))) == NULL)
        return BUFERR;

    if ((cal_eyear = (int16 *) malloc(elts * sizeof (int16))) == NULL)
        return BUFERR;

    if ((cal_sday = (int16 *) malloc(elts * sizeof (int16))) == NULL)
        return BUFERR;

    if ((cal_eday = (int16 *) malloc(elts * sizeof (int16))) == NULL)
        return BUFERR;

    if ((cal_smsec = (int32_t *) malloc(elts * sizeof (int32_t))) == NULL)
        return BUFERR;

    if ((cal_emsec = (int32_t *) malloc(elts * sizeof (int32_t))) == NULL)
        return BUFERR;

    if ((entry_year = (int16 *) malloc(elts * sizeof (int16))) == NULL)
        return BUFERR;

    if ((entry_day = (int16 *) malloc(elts * sizeof (int16))) == NULL)
        return BUFERR;

    rdvdata(vsid, SYEAR, 0, elts, (unsigned char *) cal_syear);
    rdvdata(vsid, SDAY, 0, elts, (unsigned char *) cal_sday);
    rdvdata(vsid, SMSEC, 0, elts, (unsigned char *) cal_smsec);

    rdvdata(vsid, EYEAR, 0, elts, (unsigned char *) cal_eyear);
    rdvdata(vsid, EDAY, 0, elts, (unsigned char *) cal_eday);
    rdvdata(vsid, EMSEC, 0, elts, (unsigned char *) cal_emsec);

    rdvdata(vsid, ENTRY_YEAR, 0, elts, (unsigned char *) entry_year);
    rdvdata(vsid, ENTRY_DAY, 0, elts, (unsigned char *) entry_day);

    dyear = syear;
    dday = sday;
    if (sday != eday && msec < 43200000)
        dday = eday;
    if (dday < sday)
        dyear += 1;

    for (i = elts - 1; i > 0; i--) {
        if (cal_eyear[i] == 0) /* onwards rec */ {
            if (dyear > cal_syear[i])
                break;
            if (dyear == cal_syear[i] && dday > cal_sday[i])
                break;
            if (dyear == cal_syear[i] && dday == cal_sday[i] &&
                    msec >= cal_smsec[i])
                break;
        } else /* not an onwards rec */ {
            ahead = 0;
            if (dyear > cal_syear[i])
                ahead = 1;
            else if ((dyear == cal_syear[i]) && (dday > cal_sday[i]))
                ahead = 1;
            else if ((dyear == cal_syear[i]) && (dday == cal_sday[i]) &&
                    (msec >= cal_smsec[i]))
                ahead = 1;

            if (ahead == 1) {
                if (dyear < cal_eyear[i])
                    break;
                else if ((dyear == cal_eyear[i]) && (dday < cal_eday[i]))
                    break;
                else if ((dyear == cal_eyear[i]) && (dday == cal_eday[i]) &&
                        (msec <= cal_emsec[i]))
                    break;
            }
        }
    } /* end for loop */
    *cal_year = entry_year[i];
    *cal_day = entry_day[i];

    VSdetach(vsid);
    free(cal_syear);
    free(cal_sday);
    free(cal_smsec);
    free(cal_eyear);
    free(cal_eday);
    free(cal_emsec);
    free(entry_year);
    free(entry_day);

    if (i <= 0)
        return TMERR;
    else
        return i;
}

/*-----------------------------------------------------------------------------
    Function: read_parm_data

    Returns: int32 (Status)
        On success returns 0, otherwise returns -2 indicating read error. 
    Description:
        The function  attaches to the requested vdata

    Arguments: (in calling order)
      Type       Name        I/O     Description
      ----       ----        ---     -----------
      int32      fid          I       HDF file ID
      int32      sdfid        I       HDF SD file ID
      int32      index        I       element number
      int32      idoffs[8][16] O      detector zero offset counts        
      float32    gains[8][16]  O      slopes (band*detector*gains)
      float32    fp_temps[256][8] O   fp temperatures
      float32    scan_mod[2][1285] O  scan-modulation buffer
      float64    tfactor_const[8] O        time correction constant term
      float64    tfactor_linear[8] O       time correction linear term
      float64    tfactor_exponential[8] O  time correction exponential term
      float64    cal_offset[8] O      offset term
      float64    inst_tcorr[8] O      instrument temp correction  term
      float64    fp_tcorr[8]   O      fp temp correction  term
      float64    mside1_const[8]  O   mirror side1 constant term
      float64    mside2_const[8]  O   mirror side2 constant term
      float64    mside1_linear[8] O   mirror side1 linear term
      float64    mside2_linear[8] O   mirror side2 linear term
      int16      tdi_list[256][4] O  TDI values for all 8 bands

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      03/11/94    Original development
        Lakshmi Kumar    Hughes STX      06/07/94    Updated to reflect v3.1
                                                        interface spcifications
        Gene Eplee       SAIC GSC        12/07/00    Update time correction and
                                                     mirror correction terms
        Gene Eplee       SAIC            03/09/04    Convert time correction and
                                                     mirror correction terms
                                                     to simultaneous exponentials.
        Gene Eplee       SAIC            07/26/06    Add focal plane
                                                     temperatures and
                                                     instrument electronics
                                                     temperature corrections.
        Gene Eplee	SAIC	        06/27/2007   Added focal plane and
                                                     instrument electronics 
                                                     reference temperatures
-----------------------------------------------------------------------------*/
int32_t read_parm_data(int32_t fid, int32_t sdfid, int32_t index, int32_t idoffs[8][16],
        float gains[8][16], float fp_temps[256][8],
        float scan_mod[2][1285], double *tfactor_const,
        double *tfactor_linear_1, double *tfactor_exponential_1,
        double *tfactor_linear_2, double *tfactor_exponential_2,
        double *cal_offset, double *inst_tcorr, double *inst_tref,
        double *fp_tcorr, double *fp_tref, double *mside1_const,
        double *mside1_linear_1, double *mside1_exponential_1,
        double *mside1_linear_2, double *mside1_exponential_2,
        double *mside2_const, double *mside2_linear_1,
        double *mside2_exponential_1, double *mside2_linear_2,
        double *mside2_exponential_2, int16_t tdi_list[256][4]) {

    int32_t i, slpid, parmid;
    float64 mirror_buf[8][10], tfactor_buf[8][10];

    for (i = 0; i < BANDS; i++) {
        if ((slpid = attach_vdata(fid, slp_names[i])) < 0)
            return RDERR;
        if ((rdvdata(slpid, slp_flds, index, 1, (unsigned char *) gains[i])) < 0)
            return RDERR;
        VSdetach(slpid);

        if ((parmid = attach_vdata(fid, parm_names[i])) < 0)
            return RDERR;

        if ((rdvdata(parmid, OFFSET_FLDS, index, 1,
                (unsigned char *) idoffs[i])) < 0) return RDERR;

        if ((rdvdata(parmid, TFACTOR_FLDS_NEW, index, 1,
                (unsigned char *) tfactor_buf[i])) < 0) return RDERR;

        if ((rdvdata(parmid, MIRROR_FLDS, index, 1,
                (unsigned char *) mirror_buf[i])) < 0) return RDERR;

        VSdetach(parmid);
    }

    for (i = 0; i < BANDS; i++) {
        tfactor_const[i] = tfactor_buf[i][0];
        tfactor_linear_1[i] = tfactor_buf[i][1];
        tfactor_exponential_1[i] = tfactor_buf[i][2];
        tfactor_linear_2[i] = tfactor_buf[i][3];
        tfactor_exponential_2[i] = tfactor_buf[i][4];
        cal_offset[i] = tfactor_buf[i][5];
        inst_tcorr[i] = tfactor_buf[i][6];
        inst_tref[i] = tfactor_buf[i][7];
        fp_tcorr[i] = tfactor_buf[i][8];
        fp_tref[i] = tfactor_buf[i][9];
        mside1_const[i] = mirror_buf[i][0];
        mside1_linear_1[i] = mirror_buf[i][1];
        mside1_exponential_1[i] = mirror_buf[i][2];
        mside1_linear_2[i] = mirror_buf[i][3];
        mside1_exponential_2[i] = mirror_buf[i][4];
        mside2_const[i] = mirror_buf[i][5];
        mside2_linear_1[i] = mirror_buf[i][6];
        mside2_exponential_1[i] = mirror_buf[i][7];
        mside2_linear_2[i] = mirror_buf[i][8];
        mside2_exponential_2[i] = mirror_buf[i][9];
    }

    if ((read_SDS(sdfid, TDILIST, (void *) tdi_list)) < 0)
        return RDERR;

    if ((read_SDS(sdfid, FPTEMPS, (void *) fp_temps)) < 0)
        return RDERR;

    if ((read_SDS(sdfid, SCANMOD, (void *) scan_mod)) < 0)
        return RDERR;

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function: calc_knees

    Returns: void

    Description:
        The function calc_knees calculates knee1,2,3 and 4 counts and
        radiances.  It also calculates zero offset counts.

    Arguments: (in calling order)
      Type       Name            I/O     Description
      ----       ----            ---     -----------
      int16 *    tdi              I  input TDI values for all 8 bands
      int16      tdi_list[256][4] I  TDI Detector combination table
      int32      idoffs[8][16]     I  input detector offsets
      float32    gains[8][16]     I  input gains
      float32    counts[8][4][5]  O  Digital counts corresponding to each knee
      float32    rads[8][4][5]    O  Radiances corrsponding to each knee

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      03/11/94    Original development
        Lakshmi Kumar    Hughes STX      06/07/94    Updated to reflect v3.1
                                                        interface spcifications
        W. Robinson, GSC, 13 May 97  remove double incriment of l in
                                     slopes, cnts calculation

------------------------------------------------------------------------------*/
void calc_knees(int16_t *tdi, int16_t tdi_list[256][4], int32_t idoffs[8][16],
        float gains[8][16], float counts[8][4][5],
        float rads[8][4][5]) {

    int16 dets[4];
    int32_t i, j, k, l;
    int32_t scnts[4]; /* saturation counts      */
    float32 srads[4]; /* saturation radiance    */
    float32 loc_slopes[4];
    float32 slopes[BANDS][4][4];
    int32_t cnts[BANDS][4][4];
    int32_t oindex[DETS];


    for (i = 0; i < BANDS; i++)
        for (j = 0, l = 0; j < 4; j++)
            for (k = 0; k < 4; k++) {
                /*	    slopes[i][j][k] = gains[i][l++];  */
                slopes[i][j][k] = gains[i][l];
                cnts[i][j][k] = idoffs[i][l++];
            }

    for (i = 0; i < BANDS; i++) {
        for (j = 0; j < 4; j++)
            dets[j] = tdi_list[tdi[i]][j] - 1;
        for (j = 0; j < GAINS; j++) {
            for (k = 0; k < DETS; k++) {
                scnts[k] = 1023 - cnts[i][j][dets[k]];
                srads[k] = scnts[k] * slopes[i][j][dets[k]];
                loc_slopes[k] = slopes[i][j][dets[k]];
            }

            sort_srads(srads, oindex);

            rads[i][j][0] = 0;
            for (k = 1; k < 5; k++)
                rads[i][j][k] = srads[oindex[k - 1]];

            counts[i][j][0] = 0;
            counts[i][j][1] = (scnts[oindex[0]] +
                    srads[oindex[0]] / loc_slopes[oindex[1]] +
                    srads[oindex[0]] / loc_slopes[oindex[2]] +
                    srads[oindex[0]] / loc_slopes[oindex[3]]) / 4.0;

            counts[i][j][2] = (scnts[oindex[0]] + scnts[oindex[1]] +
                    srads[oindex[1]] / loc_slopes[oindex[2]] +
                    srads[oindex[1]] / loc_slopes[oindex[3]]) / 4.0;

            counts[i][j][3] = (scnts[oindex[0]] + scnts[oindex[1]] +
                    scnts[oindex[2]] +
                    srads[oindex[2]] / loc_slopes[oindex[3]]) / 4.0;

            counts[i][j][4] = (scnts[oindex[0]] + scnts[oindex[1]] +
                    scnts[oindex[2]] + scnts[oindex[3]]) / 4.0;

        }
    }
}

/*-----------------------------------------------------------------------------
    Function: setup_scanmod

    Returns: void

    Description:
        Set up the scan-modulation correction factor array for GAC data.    
        These factors are stored in the calibration table for an entire LAC 
        scan line.

    Arguments: (in calling order)
      Type       Name        I/O     Description
      ----       ----        ---     -----------
      char *     dtype        I      data type (GAC, LAC, ...)
      float32 *  scan_mod    I/O     scan modulation correction factors

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      06/07/94    Original development
        Lakshmi Kumar    Hughes STX      03/21/96    Made it compatible for
                                                     non ANSI compilation
        Joel Gales       Futuretech      10/02/00    Use npix, nsta (start),
                                                     and ninc (stride) to 
                                                     perform subsetting and
                                                     subsampling of scan_mod
------------------------------------------------------------------------------*/

void setup_scanmod(int32_t npix, int32_t nsta, int32_t ninc,
        float scan_mod[2][1285]) {
    int32_t pixel;

    for (pixel = 0; pixel < npix; pixel++) {
        scan_mod[0][pixel] = scan_mod[0][(nsta - 1) + ninc * pixel];
        scan_mod[1][pixel] = scan_mod[1][(nsta - 1) + ninc * pixel];
    }

#ifdef DEBUG 
    printf("\n\n--------- GAC scan_mod values --------------\n");
    for (pixel = 0; pixel <= 247; pixel++) {
        printf("%4d\t%f\t%f\n", pixel, scan_mod[0][pixel], scan_mod[1][pixel]);
#endif
    }

    /*-----------------------------------------------------------------------------
        Function: sort_srad

        Returns: void

        Description:
            The function sort_srad sorts the given saturation radiances and
            returns the ordered indices.

        Arguments: (in calling order)
          Type         Name      I/O     Description
          ----         ----      ---     -----------
          float32 *    srads      I      saturation radiances
          int32  *     oindex     O      ordered indecies representing order
                                            of saturation radiances
        Notes:

        Modification history:
              Programmer     Organization      Date      Description of change
            --------------   ------------    --------    ---------------------
            Lakshmi Kumar    Hughes STX      06/07/94    Original development

    ------------------------------------------------------------------------------*/
    void sort_srads(float *srads, int32_t * oindex) {
        int32_t i, done = 0, exchange = 0, loc_index[DETS], temp_index;
        float32 loc_srads[DETS], temp;

        for (i = 0; i < DETS; i++) {
            loc_index[i] = i;
            loc_srads[i] = srads[i];
        }

        while (!done) {
            for (exchange = 0, i = 0; i < DETS - 1; i++)
                if (loc_srads[i] > loc_srads[i + 1]) {
                    exchange = 1;
                    temp = loc_srads[i];
                    temp_index = loc_index[i];
                    loc_srads[i] = loc_srads[i + 1];
                    loc_index[i] = loc_index[i + 1];
                    loc_srads[i + 1] = temp;
                    loc_index[i + 1] = temp_index;
                }
            if (!exchange)
                done = 1;
        }

        for (i = 0; i < DETS; i++)
            oindex[i] = loc_index[i];

    }



