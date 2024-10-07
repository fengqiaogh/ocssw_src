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

#include "calib_get_cal.h"
#include "calib_getcal_proto.h"

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
      float32    temps[256][8]     O      temperature correction coefficients 
      float32    scan_mod[2][1285] O      scan modulation correction factors
      float32    mirror[2][8]      O      mirror side-0 & 1 correction factors  
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
------------------------------------------------------------------------------*/
int32 get_cal(char *cal_path, int16 syear, int16 sday, int16 eday, int32 msec,
        char *dtype, int16 *tdi, int16 *cal_year, int16 *cal_day,
        int16 *ref_year, int16 *ref_day, int16 *ref_min,
        float32 temps[256][BANDS], float32 scan_mod[2][1285],
        float32 mirror[2][BANDS], float64 *t_const,
        float64 *t_linear, float64 *t_quadratic,
        float32 *cal_offs, float32 counts[BANDS][4][5],
        float32 rads[BANDS][4][5]) {
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

    if ((index = get_index(fid, syear, sday, eday, msec, cal_year, cal_day)) < 0)
        return index;

    if ((status = read_parm_data(fid, sdfid, index, idoffs, gains, temps,
            scan_mod, t_const, t_linear, t_quadratic,
            cal_offs, mirror, tdi_list)) < 0)
        return status;

    calc_knees(tdi, tdi_list, idoffs, gains, counts, rads);

    setup_scanmod(dtype, scan_mod);

    /* close the given HDF file */
    Vend(fid);
    if ((Hclose(fid)) < 0)
        return FAIL;
    if ((SDend(sdfid)) < 0)
        return FAIL;
    return SUCCEED;
}
