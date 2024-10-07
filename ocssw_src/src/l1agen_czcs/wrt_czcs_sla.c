#include <stdio.h>
#include "l1czcs.h"
#include "hdf.h"
#include "mfhdf.h"

int wrt_czcs_sla(int32 sdfid, int32 vid, int nlin, l1_data_struc l1_data)
/*******************************************************************

   wrt_czcs_sla

   purpose: write the CZCS file scan line attributes vgroup
      including:  msec, the line-by-line milliseconds
                  tilt, line-by-line tilt
                  slat, slon, line-by-line start lat, lon
                  clat, clon, line-by-line center lat, lon
                  elat, elon, line-by-line end lat, lon

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
      W. Robinson       14-Apr-2004     Original development

 *******************************************************************/ {
#define MSEC_LNAM "Scan-line time, milliseconds of day"
#define MSEC_VRNG "(0, 86399999)"
#define MSEC_UNITS "milliseconds"

#define TILT_LNAM "Tilt angle for scan line"
#define TILT_VRNG "(-20.1,20.1)"
#define TILT_UNITS "degrees"

#define SLAT_LNAM "Scan start-pixel latitude"
#define CLAT_LNAM "Scan center-pixel longitude"
#define ELAT_LNAM "Scan end-pixel latitude"
#define SLON_LNAM "Scan start-pixel longitude"
#define CLON_LNAM "Scan center-pixel latitude"
#define ELON_LNAM "Scan end-pixel longitude"
#define LAT_RNG "(-90.,90.)"
#define LON_RNG "(-180., 180.)"
#define NUM_SCN_LN "Number of Scan Lines"
#define PIX_PER_SCN "Pixels per Scan Line"
    int32 create_sds(int32, char *, int32, int32, int32 *, int32, VOIDP *);
    int32 set_dim_names(int32, char *, char *, char *);
    int dims[3] = {0, 0, 0}, sdsid, iret;

    /*
     * msec
     */
    dims[0] = nlin;
    if ((sdsid = create_sds(sdfid, "msec", DFNT_INT32, 1, (int32 *) dims, vid,
            (VOIDP) l1_data.msec)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on msec\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, NULL, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on msec\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "long_name", DFNT_CHAR,
            strlen(MSEC_LNAM) + 1, (VOIDP) MSEC_LNAM)), 0) {
        fprintf(stderr, "wrt_czcs_sla: SDsetattr failed on msec, long_name\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "valid_range", DFNT_CHAR,
            strlen(MSEC_VRNG) + 1, (VOIDP) MSEC_VRNG)), 0) {
        fprintf(stderr, "wrt_czcs_sla: SDsetattr failed on msec, valid_range\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "units", DFNT_CHAR,
            strlen(MSEC_UNITS) + 1, (VOIDP) MSEC_UNITS)), 0) {
        fprintf(stderr, "wrt_czcs_sla: SDsetattr failed on msec, units\n");
        return FAIL;
    }
    /*
     * tilt
     */
    if ((sdsid = create_sds(sdfid, "tilt", DFNT_FLOAT32, 1, (int32 *) dims, vid,
            (VOIDP) l1_data.tilt)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on tilt\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, NULL, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on tilt\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "long_name", DFNT_CHAR,
            strlen(TILT_LNAM) + 1, (VOIDP) TILT_LNAM)), 0) {
        fprintf(stderr, "wrt_czcs_sla: SDsetattr failed on tilt, long_name\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "valid_range", DFNT_CHAR,
            strlen(TILT_VRNG) + 1, (VOIDP) TILT_VRNG)), 0) {
        fprintf(stderr, "wrt_czcs_sla: SDsetattr failed on tilt, valid_range\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "units", DFNT_CHAR,
            strlen(TILT_UNITS) + 1, (VOIDP) TILT_UNITS)), 0) {
        fprintf(stderr, "wrt_czcs_sla: SDsetattr failed on tilt, units\n");
        return FAIL;
    }
    /*
     * slat
     */
    if ((sdsid = create_sds(sdfid, "slat", DFNT_FLOAT32, 1, (int32 *) dims, vid,
            (VOIDP) l1_data.slat)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on slat\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, NULL, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on slat\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "long_name", DFNT_CHAR,
            strlen(SLAT_LNAM) + 1, (VOIDP) SLAT_LNAM)), 0) {
        fprintf(stderr, "wrt_czcs_sla: SDsetattr failed on slat, long_name\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "valid_range", DFNT_CHAR,
            strlen(LAT_RNG) + 1, (VOIDP) LAT_RNG)), 0) {
        fprintf(stderr, "wrt_czcs_sla: SDsetattr failed on slat, valid_range\n");
        return FAIL;
    }
    /*
     * slon
     */
    if ((sdsid = create_sds(sdfid, "slon", DFNT_FLOAT32, 1, (int32 *) dims, vid,
            (VOIDP) l1_data.slon)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on slon\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, NULL, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on slon\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "long_name", DFNT_CHAR,
            strlen(SLAT_LNAM) + 1, (VOIDP) SLON_LNAM)), 0) {
        fprintf(stderr, "wrt_czcs_sla: SDsetattr failed on slon, long_name\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "valid_range", DFNT_CHAR,
            strlen(LON_RNG) + 1, (VOIDP) LON_RNG)), 0) {
        fprintf(stderr, "wrt_czcs_sla: SDsetattr failed on slon, valid_range\n");
        return FAIL;
    }
    /*
     * clat
     */
    if ((sdsid = create_sds(sdfid, "clat", DFNT_FLOAT32, 1, (int32 *) dims, vid,
            (VOIDP) l1_data.clat)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on clat\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, NULL, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on clat\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "long_name", DFNT_CHAR,
            strlen(CLAT_LNAM) + 1, (VOIDP) CLAT_LNAM)), 0) {
        fprintf(stderr, "wrt_czcs_sla: SDsetattr failed on clat, long_name\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "valid_range", DFNT_CHAR,
            strlen(LAT_RNG) + 1, (VOIDP) LAT_RNG)), 0) {
        fprintf(stderr, "wrt_czcs_sla: SDsetattr failed on clat, valid_range\n");
        return FAIL;
    }
    /*
     * clon
     */
    if ((sdsid = create_sds(sdfid, "clon", DFNT_FLOAT32, 1, (int32 *) dims, vid,
            (VOIDP) l1_data.clon)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on clon\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, NULL, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on clon\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "long_name", DFNT_CHAR,
            strlen(CLON_LNAM) + 1, (VOIDP) CLON_LNAM)), 0) {
        fprintf(stderr, "wrt_czcs_sla: SDsetattr failed on clon, long_name\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "valid_range", DFNT_CHAR,
            strlen(LON_RNG) + 1, (VOIDP) LON_RNG)), 0) {
        fprintf(stderr, "wrt_czcs_sla: SDsetattr failed on clon, valid_range\n");
        return FAIL;
    }
    /*
     * elat
     */
    if ((sdsid = create_sds(sdfid, "elat", DFNT_FLOAT32, 1, (int32 *) dims, vid,
            (VOIDP) l1_data.elat)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on elat\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, NULL, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on elat\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "long_name", DFNT_CHAR,
            strlen(ELAT_LNAM) + 1, (VOIDP) ELAT_LNAM)), 0) {
        fprintf(stderr, "wrt_czcs_sla: SDsetattr failed on elat, long_name\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "valid_range", DFNT_CHAR,
            strlen(LAT_RNG) + 1, (VOIDP) LAT_RNG)), 0) {
        fprintf(stderr, "wrt_czcs_sla: SDsetattr failed on elat, valid_range\n");
        return FAIL;
    }
    /*
     * elon
     */
    if ((sdsid = create_sds(sdfid, "elon", DFNT_FLOAT32, 1, (int32 *) dims, vid,
            (VOIDP) l1_data.elon)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on elon\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, NULL, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on elon\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "long_name", DFNT_CHAR,
            strlen(ELAT_LNAM) + 1, (VOIDP) ELON_LNAM)), 0) {
        fprintf(stderr, "wrt_czcs_sla: SDsetattr failed on elon, long_name\n");
        return FAIL;
    }
    if ((iret = SDsetattr(sdsid, "valid_range", DFNT_CHAR,
            strlen(LON_RNG) + 1, (VOIDP) LON_RNG)), 0) {
        fprintf(stderr, "wrt_czcs_sla: SDsetattr failed on elon, valid_range\n");
        return FAIL;
    }
    /*
     * for extra outputs of geometry and calibrated radiance data
     */
#ifdef GEOM_CAL
    dims[1] = 1968;
    if ((sdsid = create_sds(sdfid, "sen_zen", DFNT_FLOAT32, 2,
            (int32 *) dims, vid, (VOIDP) l1_data.sen_zen)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on sen_zen\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, PIX_PER_SCN, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on sen_zen\n");
        return FAIL;
    }

    if ((sdsid = create_sds(sdfid, "sen_az", DFNT_FLOAT32, 2,
            (int32 *) dims, vid, (VOIDP) l1_data.sen_az)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on sen_az\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, PIX_PER_SCN, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on sen_az\n");
        return FAIL;
    }

    if ((sdsid = create_sds(sdfid, "sol_zen", DFNT_FLOAT32, 2,
            (int32 *) dims, vid, (VOIDP) l1_data.sol_zen)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on sol_zen\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, PIX_PER_SCN, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on sol_zen\n");
        return FAIL;
    }

    if ((sdsid = create_sds(sdfid, "sol_az", DFNT_FLOAT32, 2,
            (int32 *) dims, vid, (VOIDP) l1_data.sol_az)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on sol_az\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, PIX_PER_SCN, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on sol_az\n");
        return FAIL;
    }

    if ((sdsid = create_sds(sdfid, "all_lat", DFNT_FLOAT32, 2,
            (int32 *) dims, vid, (VOIDP) l1_data.all_lat)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on all_lat\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, PIX_PER_SCN, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on all_lat\n");
        return FAIL;
    }

    if ((sdsid = create_sds(sdfid, "all_lon", DFNT_FLOAT32, 2,
            (int32 *) dims, vid, (VOIDP) l1_data.all_lon)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on all_lon\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, PIX_PER_SCN, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on all_lon\n");
        return FAIL;
    }

    if ((sdsid = create_sds(sdfid, "Lt_443", DFNT_FLOAT32, 2,
            (int32 *) dims, vid, (VOIDP) l1_data.Lt_443)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on Lt_443\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, PIX_PER_SCN, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on Lt_443\n");
        return FAIL;
    }

    if ((sdsid = create_sds(sdfid, "Lt_520", DFNT_FLOAT32, 2,
            (int32 *) dims, vid, (VOIDP) l1_data.Lt_520)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on Lt_520\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, PIX_PER_SCN, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on Lt_520\n");
        return FAIL;
    }

    if ((sdsid = create_sds(sdfid, "Lt_550", DFNT_FLOAT32, 2,
            (int32 *) dims, vid, (VOIDP) l1_data.Lt_550)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on Lt_550\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, PIX_PER_SCN, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on Lt_550\n");
        return FAIL;
    }

    if ((sdsid = create_sds(sdfid, "Lt_670", DFNT_FLOAT32, 2,
            (int32 *) dims, vid, (VOIDP) l1_data.Lt_670)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on Lt_670\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, PIX_PER_SCN, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on Lt_670\n");
        return FAIL;
    }

    if ((sdsid = create_sds(sdfid, "Lt_750", DFNT_FLOAT32, 2,
            (int32 *) dims, vid, (VOIDP) l1_data.Lt_750)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on Lt_750\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, PIX_PER_SCN, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on Lt_750\n");
        return FAIL;
    }

    if ((sdsid = create_sds(sdfid, "Lt_11500", DFNT_FLOAT32, 2,
            (int32 *) dims, vid, (VOIDP) l1_data.Lt_11500)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: create_sds failed on Lt_11500\n");
        return FAIL;
    }
    if ((iret = set_dim_names(sdsid, NUM_SCN_LN, PIX_PER_SCN, NULL)) < 0) {
        fprintf(stderr, "wrt_czcs_sla: set_dim_names failed on Lt_11500\n");
        return FAIL;
    }

#endif
    return SUCCEED;
}
