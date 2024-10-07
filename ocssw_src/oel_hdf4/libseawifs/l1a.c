/*

$Header: /app/shared/RCS/irix-5.2/seawifsd/src/hdfio/L1A.V4.5/l1a.c,v 4.18 1996/01/17 16:53:52 seawifsd Exp seawifsd $
$Log: l1a.c,v $
Revision 4.18  1996/01/17 16:53:52  seawifsd
added ON_TIME()/OFF_TIME() to report execution elapse time.

Revision 4.17  1995/12/19 14:51:04  seawifsd
fixed bug in testing main program. The bug only occurs when
L1B HDF file was not created(specified). Fixes are done to
use the fid of L1A HDF file instead of that of L1B.

Revision 4.16  1995/12/07 19:21:01  seawifsd
dded support to compile with calibration and straylight correction routines.
added snode,enode,startclat,startclon,endclat,endclon and
made it IO_SPEC_V4.4 and PROD_SPEC_V2.8 compliant.
added l2_flags, checked and set HIGHLT1 flag. Created additional
SDS to store L1B data in the test main program.

Revision 4.15  1995/05/12 13:53:00  seawifsd
1. code conformed to specifications IO_SPEC_V43, PROD_SPEC_V27,
and NONPROD_SPEC_V12(no funtionality changes, just eliminate
code related to versions before IO_SPEC_V42)

Revision 4.14  1995/05/04 15:28:24  seawifsd
changed all references to 'exist.h' to 'generic.h'.

Revision 4.13  1995/04/03 16:06:36  seawifsd
modified testing part of the l1a.c for more consistent and accurate output.

Revision 4.12  1995/02/15 18:51:23  seawifsd
added additional code to implement the objects 'scsol_z', 'entry_year',
'entry_day', and 'csol_z'.
renamed macro  LAKSHMI to EXTRA_FOR_BROWSE.
specifically declared the datatype of fid[] and fid2[].
added an optional third command parameter <skip> in the testing
main program to indicate number of scan line to skip in displaying
l1a science data. Corresponding print statements were revised to
add this support and the support of TRANSPOSE_IMAGE macro.

Revision 4.11  1995/01/27 18:20:33  seawifsd
changed the order of mirror from mirror[bands][sides] to mirror[sides][bands]

Revision 4.10  1995/01/17 20:03:18  seawifsd
Jan. 17, 1994 V4.1

Revision 4.2  1995/01/17 14:34:50  seawifsd
added calling FORTRAN cdata_() before geonav_() because the changes
in MOPS geonav.f code.

Revision 4.1  1995/01/17 14:15:37  seawifsd
Jan. 9, 1994, 4.0

Revision 3.6  1994/12/15 16:06:42  seawifsd
code conformed to IO_SPEC_V41, PROD_SPEC_V24 and NONPROD_SPEC_V10(no
changes were needed from IO_SPEC_V40)

Revision 3.5  1994/12/15 16:04:57  seawifsd
made sure that the pointers passed into the i/o routines in the test
main program are all null by declaring as static.

Revision 3.4  1994/12/06 19:31:41  seawifsd
made the code comformed to IO_SPEC_V40 and PROD_SPEC_V24.
elliminated all references to IO_SPEC_V10 and related features under that
version.
made the L1A image data to be pixel interlaced(to fit the spec.) instead
of the original band interlaced.(by define TRANSPOSE_IMAGE in wifs_conf.h)

Revision 3.3  1994/11/08 18:45:55  seawifsd
Nov. 8, 1994, 3.3a3

Revision 3.3  1994/11/08 15:03:53  seawifsd
Nov. 8, 1994, 3.3a2

Revision 1.1.1.5  1994/11/03 20:24:00  frank
renamed MIN macro to MINV.

Revision 1.1.1.4  1994/11/01 16:02:28  frank
added processing of 'meta_l1a.infiles' because of the new requirement
from L1A Browse routines.

Revision 1.1.1.3  1994/10/21 19:43:51  frank
turned off some message printing by using PRINTF instead of printf.

Revision 1.1.1.2  1994/10/05 20:14:01  frank
1. updated L1A source code to fit the specification of the
"SeaWiFS OPERATIONAL ARCHIVE PRODUCT SPECIFICATIONS" v 1.0, 07/22/94,
and the "INTERFACE SPECIFICATIONS FOR SeaWiFS OPERATIONAL PRODUCT
INPUT, OUTPUT, AND SENSOR CALIBRATION SOFTWARE" v 3.3 07/12/94

Revision 1.1.1.1  1994/05/23 19:04:43  frank
variables ptr_orb_vec, ptr_l_vert, ptr_sun_ref, ptr_att_ang, ptr_sen_mat,
ptr_scan_ell, ptr_nflag, nsta, ninc, and npix are declared as PRIVATE(static)
'int32_t *longptr' was added in get_l1a_record() routine and used as the
lhs to nav.nflag so that to match the prototype.

Revision 1.2  1994/05/10 18:44:37  seawifst
May 6, 1994 version 1.2

Revision 1.1  1994/04/19 13:23:38  seawifst
Initial revision

Removed all code that was only compiled when GENERAL, USE_STRUCT,
or TESTCODE were defined. (687 lines of code removed.)
Norman Kuring		6-Nov-1996

Added code to handle "Calibration" Vgroup data.
Norman Kuring		6-Nov-1996

Removed all code that was only compiled when EVERYTHING was defined.
Removed non-prototype style function definitions.
Norman Kuring		7-Nov-1996

Put back all of the GENERAL and EVERYTHING code.  Heavy sigh.
Norman Kuring		25-Nov-1996

fix values of l1a_data and geoloc when stray light processing is done.
W. Robinson  3-Apr-1997

W. Robinson , GSC, 20 Mar 98  update the dark restore calculation to filter
unreasonable values and to take the median
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
/* for ceil used in STRAY_LIGHT_COEF macro */
#include <math.h>

// why reinvent the wheel
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#include "hdf.h"
//#include	"hdfhdr.h"
//#include	"hdfmac.h"
//#include	"usrhdr.h"
#include "usrmac.h"

//#include	"generic.h"

#include "ffm.h"
#include "cdl_object.h"
#include "hdf_object.h"
#include "SeaWiFS.h"
#include "navigation.h"
#include "l1a.h"
#include "level_1a_index.h"

#include "datatype.h"
#include "tlm.h"
#include "l1a_proto.h"

/* prototype for get_WIFSinfo() */
#include "WIFSHDF.h"

/* stray light */
#include "st_lt.h"
#include "cal_l1a.h"
#include "get_cal.h"
#ifndef NO_SET_HIGHLT1
#include "l2_flags_sw.h"
#endif /* !NO_SET_HIGHLT1 */

/* cpp_processing_start_here_please_do_not_delete_this_line */

/* Make following as default in all future release(040194)		*/


#define FIRST_GET_L1A_OPEN 1

#ifndef DEFAULT_BUFFER_BLKSIZE
#define DEFAULT_BUFFER_BLKSIZE 10
#endif

#define BLKSIZE_NOT_CHANGED 0
#define BLKSIZE_CHANGED  1

/*
   Following are for get_l1a_record routine to keep track of:
   1. current access record number
   2. total record number
 */

#define CLEAR 0
#define INIT 1
#define RECORDS 2

#ifndef NO_SET_HIGHLT1
#define FIRST_KNEE 1

/* global storage space for BAND 8 gain value to be used in stray_light	*/
short gain8;

#endif /* !NO_SET_HIGHLT1 */

#define STRAY_LIGHT_COMPLETE 0
#define STRAY_LIGHT_LATE_START 1
#define STRAY_LIGHT_EARLY_STOP 2
#define STRAY_LIGHT_RANDOM 4
#define PIX_SUB_V 4.0
#define GAC_STRAY_LIGHT_COEF(v) (ceil(roundingup = v/PIX_SUB_V))
#define LAC_STRAY_LIGHT_COEF(v) (v)
#define STRAY_LIGHT_COEF(v,dtyp) (strcmp(dtyp,"GAC")?LAC_STRAY_LIGHT_COEF(v):GAC_STRAY_LIGHT_COEF(v))

#define NTRY_WARNING 5
#define NTRY_LIMIT 20


#ifdef L1A_INIT_NEEDED

#ifdef PROTOTYPE
PRIVATE void init(void)
#else

PRIVATE void init()
#endif
{
    char *FUNC = "l1a_init";
    int i, j;

    ON_TIME();
    /* invalidate all fids						*/
    for (i = 0; i < MAX_HDF_L1AGET; i++) {
        fid[i] = -1;
        blksize_changed[i] = BLKSIZE_NOT_CHANGED;
        blksize[i] = DEFAULT_BUFFER_BLKSIZE;
        /* set to zero, indicate no call to get_l1a_record yet	*/
        blksrec[i] = 0;
        cal_table_file[i] = NULL;
        stray_light_flag[i] = 0;
        ltyp_frac_coef[i] = 0.0;
        out_band_flag[i] = 0;
        recursive_flag[i] = FALSE;
        recursive_free = 0;
        stray_light_calling_flag[i] = STRAY_LIGHT_COMPLETE;
        stray_light_scan_no[i] = 1;
    }
    OFF_TIME();
}

#endif // L1A_INIT_NEEDED

/*------------------------------------------------------------------------------

    calculate the mean and standard deviation of the 8 bands of 
    dark-restore values for each band
    W. Robinson, GSC, 18 Mar 98  change dark computation to filter
        dark values <5, > 35 and then find the median on the remainder
        Note routines d_stats, i16comp and d_sd are a part of this now

------------------------------------------------------------------------------*/
#ifdef PROTOTYPE
int dark_rest_stat(int16 *data, int nrec, float *dark_mean, float *dark_std)
#else

int dark_rest_stat(data, nrec, dark_mean, dark_std)
int16 *data; /* 8 bands of dark restore data */
int nrec; /* number of scan lines */
float *dark_mean; /* 8 bands of dark restore mean values */
float *dark_std; /* 8 bands of dark restore std values */
#endif
{
    double *ddata;
    double dmedian;
    double dark;
    int32 num_out, i, j;
    ddata = malloc(nrec * sizeof (double));

    for (i = 0; i < 8; i++) {
        for (j = 0; j < nrec; j++) {
            ddata[j] = (double) *(data + j * 8 + i);
        }
        gsl_sort(ddata, 1, nrec);
        dmedian = gsl_stats_median_from_sorted_data(ddata, 1, nrec);
        memset(ddata, 0, nrec);
        num_out = 0;
        for (j = 0; j < nrec; j++) {
            dark = (double) *(data + j * 8 + i);
            if (dark > dmedian * 0.85 && dark < dmedian * 1.15) {
                ddata[num_out] = (double) *(data + j * 8 + i);
                num_out++;
            }
        }
        *(dark_mean + i) = (float) gsl_stats_mean(ddata, 1, num_out);
        *(dark_std + i) = (float) gsl_stats_sd(ddata, 1, num_out);
    }
    free(ddata);
    return 0;

}
