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

#include "get_cal_osmi.h"
#include "getcal_proto_osmi.h"
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
int32_t get_cal_osmi(char *cal_path, int16_t syear, int16_t sday, int16_t eday, int32_t msec,
        int16_t *cal_year, int16_t *cal_day, int32_t *cal_msec,
        float *eoffset, float *egain, float *temp,
        float *mirror, float *t_const,
        float *t_linear, float *t_quadratic,
        float *slope, float *dc,
        float *sm) {
    int32 fid, sdfid, status = 0, index = 0;


    /*  open given HDF file  */
    if ((fid = Hopen(cal_path, DFACC_RDONLY, 0)) < 0)
        return FAIL;

    sdfid = SDstart(cal_path, DFACC_RDONLY);
    Vstart(fid);
    /*
     * Read global attributes Reference Year, Reference Day and Reference Minute
     *  from calibration data file
     */


    if ((index = get_index_osmi(fid, syear, sday, eday, msec, cal_year, cal_day, cal_msec)) < 0)
        return index;

    if ((status = read_parm_data_osmi(fid, sdfid, index, eoffset, egain, temp, mirror,
            t_const, t_linear, t_quadratic,
            slope, dc, sm)) < 0)
        return status;

    /* close the given HDF file */
    Vend(fid);
    if ((Hclose(fid)) < 0)
        return FAIL;
    if ((SDend(sdfid)) < 0)
        return FAIL;
    return SUCCEED;
}
