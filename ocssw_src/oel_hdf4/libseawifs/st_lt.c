/*-----------------------------------------------------------------------------
    Function: stray_light_corr

    Returns: int32 (status)
        Returns as 0, "not done", to indicate that the routine did not 
        return any useful information for the last call, but was filling
        its buffers with scan line data; and 1, "done", to indicate that 
        the l1b_data and sl_flag arrays for scan line sl_scan have been 
        set.

    Description:
        The function stray_light_corr calls stray_light_lac or stray_light_gac
        depending upon input data type.  If the processed scan line is 
        returned (status = DONE) checks sl_flag arry for stray light flags
        and sets corresponding l2_flag bits.

    Parameters: (in calling order)
        Type      Name         I/O   	Description
        ----      ----         ---   	-----------
        int32	*initial	I	Flag, should be set to "true" (a value
                                        of 1) if it is the fst call for a scene
        float32  Ltyp_frac	I	fraction of Ltypical for band 8
        float32  Styp_frac	I	fraction value that will be applied
                                        to knee value to calculate straylight
                                        threshold
        int32 	 nscans		I	the number of scan lines in the scene
                                        (to be processed).  Must be set to <5
                                        to process only in the along scan
        int32 	 nsamples	I	number of pixel data values in 
                                        l1b_data and sl_flag
        int32	 scan_no	I	scan line number of l1b_data input
        char 	*dtype		I 	l1b_data data type
        float32	*l1b_data      I/O	array of nsamples long containing 
                                        l1b_data of scan line scan_no
        int32	*sl_scan	O	scan line number for which the 
                                        returned l1b_data apply
        int16	*l2_flags       O	an array of nsamples long containing
                                        level2 flags of scan line scan_no. The
                                        corresponding pixel flags will be set 
                                        for affected pixels
        int32	*AS_pixels	I	Along Scan pixels (constant set to 12)

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar	 Hughes STX	 02/10/97    Added input argument 
                                                     Styp_frac & is passed to
                                                     stlt lac and gac routines
        Lakshmi Kumar	 Hughes STX	 03/21/96    Fixed code for
                                                     non ANSI compile options
        Lakshmi Kumar    Hughes STX      10/12/95    Original development
        W. Robinson, GSC                 4 Nov 97    remove hard coding of
                                                     range and fix problem
                                                     with sl_flag of > 1000
------------------------------------------------------------------------------*/

#include "st_lt.h"
#include "st_proto.h"
#include <mfhdf.h>


int32_t
stray_light_corr(int32_t *initial, float Ltyp_frac, float Styp_frac,
        int32_t nscans, int32_t nsamples, int32_t scan_no, char *dtype,
        int16_t gn, float *rads, float *l1b_data, int32_t *sl_scan,
        int16_t *l2_flags, int32_t *AS_pixels) {
    int32 status; /* Function returned value 	*/
    int32 sl_flag[MAXSAMPS]; /* Flag list identifying BTs    */
    int32 range; /* Range of pixels flagged bad  */
    int32 i, left, right;
    static float tsum = 0;
    static float t1, t2;
    div_t quot1;

    /* call appropriate correction routine according to data type */
    for (i = 0; i < nsamples; i++)
        l2_flags[i] = 0;

    t1 = (float) clock() / CLOCKS_PER_SEC;

    if ((strcmp(dtype, "GAC")) == 0)
        status = stray_light_gac(initial, Ltyp_frac, Styp_frac, nscans, nsamples,
            scan_no, gn, rads, l1b_data, sl_scan, sl_flag);
    else
        status = stray_light_lac(initial, Ltyp_frac, Styp_frac, nscans, nsamples,
            scan_no, gn, rads, l1b_data, sl_scan, sl_flag);

    t2 = (float) clock() / CLOCKS_PER_SEC;
    tsum += t2 - t1;

    /*
    if (scan_no == nscans)
        printf("\n%d\t%f\n", nscans, tsum);
     */


    /* if status = done do the following loop */
    if (status == DONE) {
        range = *AS_pixels;
        /* WDR keep the variable above     range = 14;  */
        if ((strcmp(dtype, "GAC")) == 0) {
            quot1 = div(range, 4);
            if (quot1.rem > 0)
                range = quot1.quot + 1;
            else
                range = quot1.quot;
        }
        for (i = 0; i < nsamples; i++) {
            if (sl_flag[i] < 1000) {
                if ((sl_flag[i] > -10) && (sl_flag[i] <= range))
                    l2_flags[i] = l2_flags[i] + 256;
            } else {
                left = sl_flag[i] / 1000;
                /* WDR as was:right = sl_flag[i] - 1000; */
                right = sl_flag[i] % 1000;
                if ((left <= range) || (right <= range))
                    l2_flags[i] = l2_flags[i] + 256;
            } /* end else */
        } /* end for */
    } /* end if */
    return status;
}
