#include <stdio.h>
#include "l1czcs.h"
#include "hdf.h"
#include "mfhdf.h"

int wrt_czcs_qual(int32 sdfid, int32 vid, int nlin, l1_data_struc l1_data)
/*******************************************************************

   wrt_czcs_qual

   purpose: write the CZCS file quality data to the raw data vgroup
      including:  cal_sum - errors for the line
                  cal_scan - band specific errors

   Returns type: int - FAIL if failed, else SUCCEED

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             sdfid            I      sd file id from hdf
      int32             vid              I      v-group id
      int               nlin             I      number of lines of data
      l1_data_struc     l1_data          I      structure of data arrays

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       20-Apr-2004     Original development

 *******************************************************************/ {
#define CAL_SUM_LNAM "Calibration Qcality Summary"
#define CAL_SCAN_LNAM "Calibration Quality per Scan"

    int32 create_sds(int32, char *, int32, int32, int32 *, int32, VOIDP *);
    int32 set_dim_names(int32, char *, char *, char *);
    int dims[3] = {0, 0, 0}, sdsid, iret;
    char *tot_line = "Number of Scan Lines";

    /*
     * cal_sum
     */
    dims[0] = nlin;
    dims[1] = 5;
    if ((sdsid = create_sds(sdfid, "cal_sum", DFNT_UINT8, 2,
            (int32 *) dims, vid, (VOIDP) l1_data.cal_sum)) < 0) {
        fprintf(stderr, "wrt_czcs_qual: create_sds failed on cal_sum\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, tot_line, "num_qual", NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_qual: set_dim_names failed on cal_sum\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "long_name", DFNT_CHAR,
            strlen(CAL_SUM_LNAM) + 1, (VOIDP) CAL_SUM_LNAM)), 0) {
        fprintf(stderr, "wrt_czcs_qual: SDsetattr failed on cal_sum, long_name\n");
        return FAIL;
    }
    /*
     * cal_scan
     */
    dims[1] = 6;
    if ((sdsid = create_sds(sdfid, "cal_scan", DFNT_UINT8, 2,
            (int32 *) dims, vid, (VOIDP) l1_data.cal_scan)) < 0) {
        fprintf(stderr, "wrt_czcs_qual: create_sds failed on cal_scan\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, tot_line, "Number of Bands", NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_qual: set_dim_names failed on cal_scan\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "long_name", DFNT_CHAR,
            strlen(CAL_SCAN_LNAM) + 1, (VOIDP) CAL_SCAN_LNAM)), 0) {
        fprintf(stderr,
                "wrt_czcs_qual: SDsetattr failed on cal_scan, long_name\n");
        return FAIL;
    }
    return SUCCEED;
}
