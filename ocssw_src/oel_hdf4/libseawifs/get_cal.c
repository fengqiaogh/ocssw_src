/*-----------------------------------------------------------------------------
    File: get_cal.c

    Contents:
        get_cal         -  opens the given calibration HDF file, retrieves 
                           the sensor calibration data, calculates knee counts 
                           and radiances, closes the file and returns status 

    Other relevant files:
        get_cal.h       -  various #defined constants, and also includes hdf.h
        getcal_proto.h  -  prototypes for get_cal functions
        get_cal_misc.c  -  a lower layer of calibration input functions

    Notes:
        o A test program appears at the end of this file.  It reads the     
          given calibration file.  The input calibration file name is hard 
          coded.  To test the code, uncomment the TESTCODE section
          and compile with -DTESTCODE on the compile line.  
          To get debug output include, "-DDEBUG" on the compile line.

        Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      03/11/94    Original development
        Lakshmi Kumar    Hughes STX      06/07/94    Updated to reflect v3.1 
                                                        interface spcifications
        Lakshmi Kumar 	 Hughes STX 	 11/15/94    Modified comments
        Lakshmi Kumar	 Hughes STX	 05/22/96    Modified to read revised 
                                                     calibration table
                                                     Removed time_factor
                                                     Added output arguments
                                                     reference year, day,
                                                     min, t_const, t_linear,
                                                     t_quadratic & cal_offs.
                                                     Ref. V5.0 I/O specs.
        Lakshmi Kumar	 Hughes STX	 11/01/96    fixed output parameters 
                                                     of get_ref_time call.
        Lakshmi Kumar	 Hughes STX	 03/17/97    Chnaged idoffs[8][4] to
                                                     idoffs[8][16] and non-
                                                     prototype declarations 
                                                     have been removed
------------------------------------------------------------------------------*/

#include "get_cal.h"
#include "getcal_proto.h"
#include <mfhdf.h>

/*-----------------------------------------------------------------------------
    Function: get_cal  

    Returns: int32 (status)
        Returns a status code of 0 when successful.  Otherwise returns
        -1 - to indicate file open/close error
        -2 - to indicate read error
        -3 - to indicate time error (if the given time cannot be found)
        -4 - to indicate insufficient memory error

    Description:
        The function get_cal reads given HDF file, determines the detector 
        combination to be used, calculates knee counts and radiances.

    Arguments: (in calling order)
      Type       Name             I/O     Description
      ----       ----             ---     -----------
      char *     cal_path          I      calibration file path 
      int16      syear             I      year of data start time
      int16      sday              I      day of year for data start time
      int16      eday              I      day of year for data end time
      int32      msec              I      milliseconds of the day
      char  *    dtype             I      data type flag 
      int16 *    tdi          	   I      input TDI for all 8 bands
      int16 *    cal_year	   O      the year the calibration table entry
                                                was make
      int16 *    cal_day	   O      the day of year the calibration table
                                                entry was made
      float32    fp_temps[256][8]  O      focal plane temperatures 
      float32    scan_mod[2][1285] O      scan modulation correction factors
      float64    t_const[8]        O      time correction constant term
      float64    t_linear[8]       O      time correction linear term
      float64    t_exponential[8]  O      time correction exponential term
      float64    cal_offs[8]       O      calibration offset term
      float64    inst_tcorr[8]     O      instrument temperature correction term
      float64    fp_tcorr[8]       O      fp temperature correction term
      float64    ms1_const[8]      O      mirror side1 constant term
      float64    ms2_const[8]      O      mirror side2 constant term
      float64    ms1_linear[8]     O      mirror side1 linear term
      float64    ms2_linear[8]     O      mirror side2 linear term
      float32    counts[8][4][5]   O      Digital cnts corresponding to eh knee
      float32    rads[8][4][5]     O      Radiances corresponding to each knee

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      03/11/94    Original development
        Lakshmi Kumar    Hughes STX      06/07/94    Updated to reflect v3.0 
                                                        interface spcifications
        Lakshmi Kumar	 Hughes STX	 02/07/95    Added output arguments
                                                     cal_year and cal_day
                                                     (ref. I/O specs v4.2)
        Lakshmi Kumar	 Hughes STX	 05/22/96    Modified to read revised 
                                                     calibration table
        Lakshmi Kumar	 Hughes STX	 11/01/96    removed '&' sign from
                                                     ref_year, ref_day & ref_
                                                     minute variables.
        Joel Gales       Futuretech      10/02/00    Pass npix, nsta, & ninc
                                                     to setup_scanmod
        Gene Eplee       SAIC GSC        12/07/00    Update time correction and
                                                     mirror correction terms
        Gene Eplee       SAIC            03/09/04    Convert time correction and
                                                     mirror correction terms
                                                     to simultaneous exponentials.
        Gene Eplee       SAIC            07/26/06    Add focal plane
                                                     temperatures and 
                                                     instrument electronics
                                                     temperature corrections
        Gene Eplee	 SAIC	         06/27/2007  Added focal plane and
                                                     instrument electronics 
                                                     reference temperatures
------------------------------------------------------------------------------*/

int32_t get_cal(char *cal_path, int16_t syear, int16_t sday, int16_t eday, int32_t msec,
        int32_t npix, int32_t nsta, int32_t ninc,
        char *dtype, int16_t *tdi, int16_t *cal_year, int16_t *cal_day,
        int16_t *ref_year, int16_t *ref_day, int16_t *ref_min,
        float fp_temps[256][BANDS], float scan_mod[2][1285],
        double *t_const, double *t_linear_1, double *t_exponential_1,
        double *t_linear_2, double *t_exponential_2,
        double *cal_offs, double *inst_tcorr, double *inst_tref,
        double *fp_tcorr, double *fp_tref, double *ms1_const,
        double *ms1_linear_1, double *ms1_exponential_1,
        double *ms1_linear_2, double *ms1_exponential_2,
        double *ms2_const, double *ms2_linear_1,
        double *ms2_exponential_1, double *ms2_linear_2,
        double *ms2_exponential_2, float counts[BANDS][4][5],
        float rads[BANDS][4][5]) {
    int16 tdi_list[256][4];
    int32 fid, sdfid, status = 0, index = 0;
    int32 idoffs[BANDS][16];
    float32 gains[BANDS][16];


    /*  open given HDF file  */
    if ((fid = Hopen(cal_path, DFACC_RDONLY, 0)) < 0)
        return FAIL;
    sdfid = SDstart(cal_path, DFACC_RDONLY);
    Vstart(fid);

    /*
     * Read global attributes Reference Year, Reference Day and Reference Minute
     *  from calibration data file
     */

    if ((status = get_ref_time(sdfid, ref_year, ref_day, ref_min)) < 0)
        return status;

    if ((index = get_tindex(fid, syear, sday, eday, msec, cal_year, cal_day)) < 0)
        return index;

    if ((status = read_parm_data(fid, sdfid, index, idoffs, gains, fp_temps,
            scan_mod, t_const, t_linear_1, t_exponential_1, t_linear_2,
            t_exponential_2, cal_offs, inst_tcorr, inst_tref, fp_tcorr,
            fp_tref, ms1_const, ms1_linear_1, ms1_exponential_1,
            ms1_linear_2, ms1_exponential_2, ms2_const, ms2_linear_1,
            ms2_exponential_1, ms2_linear_2, ms2_exponential_2,
            tdi_list)) < 0)
        return status;

    calc_knees(tdi, tdi_list, idoffs, gains, counts, rads);

    setup_scanmod(npix, nsta, ninc, scan_mod);

    /* close the given HDF file */
    Vend(fid);
    if ((Hclose(fid)) < 0)
        return FAIL;
    if ((SDend(sdfid)) < 0)
        return FAIL;
    return SUCCEED;
}

