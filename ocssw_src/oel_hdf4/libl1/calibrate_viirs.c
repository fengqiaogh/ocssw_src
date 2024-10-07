
#include "nc4utils.h"

/**
 * 
 * @param calfile
 * @param UTC58usec
 * @param ftable
 * @return 
 */

int load_fcal_lut(char* calfile, int64_t UTC58usec, double ****ftable) {
    int fileid, varid;
    var_str_nc *vars;
    size_t ntimes, nbands, ndets, nms;
    size_t itime, iband, idet, ims;
    int64_t *times = NULL;

    /* VIIRS-SDR-F-LUT:Dimensions = "[time, band, detector, gain, mirrorside]" */
    size_t start[] = {0, 0, 0, 0, 0};
    size_t count[] = {2, 1, 1, 1, 1}; // reading high-gain only
    double flut_vals[] = {0, 0}; // before & after values
    double factor; // for linear interpolation

    /* Open file and get info for root-level variables */
    TRY_NC(__FILE__, __LINE__,
            nc_open(calfile, NC_NOWRITE, &fileid));
    vars = load_vars_nc(fileid);

    /* Populate Beginning_Time_IET array */
    /* IET = UTC*1e6 (leap-second corrected microseconds since 1958) */
    TRY_NC(__FILE__, __LINE__,
            nc_inq_varid(fileid, "Beginning_Time_IET", &varid));
    ntimes = vars[varid].dim[0].len;
    readall_var(&vars[varid]);
    times = vars[varid].data;

    /* Determine time indices bracketing specified time */
    if (UTC58usec < times[0]) {
        itime = 0;
        printf("\n- W - %s, %d: WARNING: ", __FILE__, __LINE__);
        printf("granule time of %ld ", (long int) UTC58usec);
        printf("is below F-table start time of %ld\n", (long int) times[0]);
    } else if (UTC58usec > times[ntimes - 1]) {
        itime = ntimes - 2;
        printf("\n- W - %s, %d: WARNING: ", __FILE__, __LINE__);
        printf("granule time of %ld ", (long int) UTC58usec);
        printf("is above F-table end time of %ld\n", (long int) times[ntimes - 1]);
    } else {
        for (itime = 0; itime < ntimes - 1; itime++)
            if (times[itime + 1] > UTC58usec)
                break;
    }

    /* Calculate time interpolation factor */
    factor =
            (double) (UTC58usec - times[itime]) /
            (double) (times[itime + 1] - times[itime]);

    /* Get VIIRS-SDR-F-LUT variable info */
    TRY_NC(__FILE__, __LINE__,
            nc_inq_varid(fileid, "VIIRS-SDR-F-LUT", &varid));

    /* Assume dimensions for M-bands, high-gain only: */
    nbands = 16;
    ndets = 16;
    nms = 2;
    /* if reading all bands & gains, use:
       nbands = vars[varid].dim[1].len;  // =22; covers all bands
       ndets  = vars[varid].dim[2].len;  // =32 to cover I-bands
       ngain  = vars[varid].dim[3].len;  // =3 high, low, mixed
       nms    = vars[varid].dim[4].len;  // =2 HAMside
     */

    /* Allocate space for F-cal array: [band, detector, mirrorside] */
    /*
      Modified to allocate a three-dimensional array passed back
      through a command-line parameter (ftable)

      JMG 03/08/16
     */

    *ftable = (double ***) malloc(nbands * sizeof (double**));
    for (iband = 0; iband < nbands; iband++) {
        ftable[0][iband] = (double **) malloc(ndets * sizeof (double*));
        for (idet = 0; idet < ndets; idet++)
            ftable[0][iband][idet] = (double *) malloc(nms * sizeof (double));
    }

    /* Read F-LUT values and interpolate to observation time */
    start[0] = itime;
    for (iband = 0; iband < nbands; iband++) {
        start[1] = iband + 5; // M-bands start at index 5
        for (idet = 0; idet < ndets; idet++) {
            start[2] = idet;
            for (ims = 0; ims < nms; ims++) {
                start[4] = ims;

                TRY_NC(__FILE__, __LINE__,
                        nc_get_vars_double(fileid, varid,
                        start, count, NULL, flut_vals));

                ftable[0][iband][idet][ims] =
                        (flut_vals[0] * (1.0 - factor)) +
                        (flut_vals[1] * factor);

            } // ims
        } // idet
    } // iband

    /* Close file and free memory allocated for related variables */
    //free_vars_nc(vars,2);
    TRY_NC(__FILE__, __LINE__, nc_close(fileid));

    return 0;
}