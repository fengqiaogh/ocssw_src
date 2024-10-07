/*-----------------------------------------------------------------------------
    File : get_cal_misc.c

    Contents:
        get_ref_time	-  reads reference time
        get_index	-  reads time vdata and returns appropriate index to 
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
        Lakshmi Kumar	 Hughes STX	 03/17/97    Removed non-ANSI proto
                                                     declarations.  In-detector
                                                     offsets are redefined as
                                                     idoffs[8][16]. 
------------------------------------------------------------------------------*/
#include <hdf4utils.h>

#include "calib_get_cal.h"
#include "calib_getcal_proto.h"

static char slp_flds[] = {
    "g1d1,g1d2,g1d3,g1d4,g2d1,g2d2,g2d3,g2d4,g3d1,g3d2,g3d3,g3d4,g4d1,g4d2,g4d3,g4d4"
};

static char OFFSET_FLDS[] = {
    "g1offs1,g1offs2,g1offs3,g1offs4,g2offs1,g2offs2,g2offs3,g2offs4,g3offs1,g3offs2,g3offs3,g3offs4,g4offs1,g4offs2,g4offs3,g4offs4"
};

static char TFACTOR_FLDS[] = {
    "t_const,t_linear,t_quadratic"
};

static char CORRECTION_FLDS[] = {
    "cal_offs,mirror1,mirror2"
};

static char *slp_names[] = {
    "B1Slopes", "B2Slopes", "B3Slopes", "B4Slopes", "B5Slopes", "B6Slopes",
    "B7Slopes", "B8Slopes"
};

static char *parm_names[] = {
    "B1Parms", "B2Parms", "B3Parms", "B4Parms", "B5Parms", "B6Parms",
    "B7Parms", "B8Parms"
};

/*-----------------------------------------------------------------------------
    Function: get_index 

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
int32 get_ref_time(int32 sdfid, int16 *ref_year, int16 *ref_day, int16 *ref_min) {
    int32 attrnum;

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
    Function: get_index 

    Returns: int32 (index)
        On successful it returns the index of the given time entry and if
        time entry not found, returns -3.

    Description:
        The function get_index reads time vdata and searches for the
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
------------------------------------------------------------------------------*/

int32_t get_index(int32_t fid, int16_t syear, int16_t sday, int16_t eday, int32_t msec,
        int16_t *cal_year, int16_t *cal_day) {

    int16 dyear, dday, *cal_syear, *cal_sday, *cal_eyear, *cal_eday;
    int16 *entry_year, *entry_day;
    int32 i, *cal_smsec, *cal_emsec, vsid, elts;

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

    if ((cal_smsec = (int32 *) malloc(elts * sizeof (int32))) == NULL)
        return BUFERR;

    if ((cal_emsec = (int32 *) malloc(elts * sizeof (int32))) == NULL)
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
        if (cal_eyear[i] == 0) { /* onwards rec */
            if (dyear > cal_syear[i])
                break;
            if (dyear == cal_syear[i] && dday > cal_sday[i])
                break;
            if (dyear == cal_syear[i] && dday == cal_sday[i] &&
                    msec >= cal_smsec[i])
                break;
        } else { /* not an onwards rec */
            if (dyear > cal_syear[i]) {
                if (dyear < cal_eyear[i])
                    break;
                if (dyear == cal_eyear[i]) {
                    if (dday < cal_eday[i])
                        break;
                    if (dday == cal_eday[i] && msec <= cal_emsec[i])
                        break;
                }
            } else
                if (dyear == cal_syear[i]) {
                if (dyear < cal_eyear[i])
                    break;
                if (dyear == cal_eyear[i]) {
                    if (dday >= cal_sday[i] && dday < cal_eday[i])
                        break;
                    if (dday >= cal_sday[i] && dday == cal_eday[i] &&
                            msec <= cal_emsec[i])
                        break;
                }
            }
        }
    }

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
      int32      fid          I      HDF file ID
      int32      sdfid        I      HDF SD file ID
      int32      index        I      element number
      int32      idoffs[8][16] O      detector zero offset counts        
      float32    gains[8][16] O      slopes (band*detector*gains)
      float32    temps[256][2]O      temp and ref temp correction factors
      float32    scan_mod[2][1285]O  scan-modulation buffer
      float32    mirror[8][2] O      mirror side1 and side2 for all bands
      int16      tdi_list[256][4] O  TDI values for all 8 bands

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      03/11/94    Original development
        Lakshmi Kumar    Hughes STX      06/07/94    Updated to reflect v3.1
                                                        interface spcifications
------------------------------------------------------------------------------*/
int32 read_parm_data(int32 fid, int32 sdfid, int32 index, int32 idoffs[8][16],
        float32 gains[8][16], float32 temps[256][8],
        float32 scan_mod[2][1285], float64 *tfactor_const,
        float64 *tfactor_linear, float64 *tfactor_quadratic,
        float32 *cal_offset, float32 mirror[2][8],
        int16 tdi_list[256][4]) {

    int32 i, slpid, parmid;
    float32 parm_buf[8][3];
    float64 tfactor_buf[8][3];

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

        if ((rdvdata(parmid, TFACTOR_FLDS, index, 1,
                (unsigned char *) tfactor_buf[i])) < 0) return RDERR;

        if ((rdvdata(parmid, CORRECTION_FLDS, index, 1,
                (unsigned char *) parm_buf[i])) < 0) return RDERR;

        VSdetach(parmid);
    }

    for (i = 0; i < BANDS; i++) {
        tfactor_const[i] = tfactor_buf[i][0];
        tfactor_linear[i] = tfactor_buf[i][1];
        tfactor_quadratic[i] = tfactor_buf[i][2];
        cal_offset[i] = parm_buf[i][0];
        mirror[0][i] = parm_buf[i][1];
        mirror[1][i] = parm_buf[i][2];
    }

    if ((read_SDS(sdfid, TDILIST, (void *) tdi_list)) < 0)
        return RDERR;

    if ((read_SDS(sdfid, TEMPS, (void *) temps)) < 0)
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
void calc_knees(int16 *tdi, int16 tdi_list[256][4], int32 idoffs[8][16],
        float32 gains[8][16], float32 counts[8][4][5],
        float32 rads[8][4][5]) {

    int16 dets[4];
    int32 i, j, k, l;
    int32 scnts[4]; /* saturation counts      */
    float32 srads[4]; /* saturation radiance    */
    float32 loc_slopes[4];
    float32 slopes[BANDS][4][4];
    int32 cnts[BANDS][4][4];
    int32 oindex[DETS];


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
------------------------------------------------------------------------------*/

void setup_scanmod(char *dtype, float32 scan_mod[2][1285]) {
    int32 i, pixel;
    char loc_dtype[20];

    for (i = 0; i < strlen(dtype); i++) /* convert to lower case */
        loc_dtype[i] = tolower(dtype[i]);
    loc_dtype[i] = '\0';

    if ((strcmp(loc_dtype, GAC)) == 0)
        for (pixel = 0; pixel <= 247; pixel++) {
            scan_mod[0][pixel] = scan_mod[0][146 + 4 * pixel];
            scan_mod[1][pixel] = scan_mod[1][146 + 4 * pixel];
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
    void sort_srads(float32 *srads, int32 * oindex) {
        int32 i, done = 0, exchange = 0, loc_index[DETS], temp_index;
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

#ifdef DEBUG
        for (i = 0; i < 4; i++)
            printf("\n %d\t%8.5f ", oindex[i], loc_srads[i]);
#endif
    }
