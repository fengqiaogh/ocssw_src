#include "viirs_sim_sdr.h"
#include "l12_parms.h"
#include <string.h>
#include <genutils.h>

static h5io_str file_id, *dat_ids;

int rd_rhos_scan(char *file, int npix, int iscn, int ndet, float **rhos)
/*-----------------------------------------------------------------------------
    Program:   rd_rhos_scan.c

    Description:  read a scan of reflectance data

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        char *    file          I    file of reflectance
        int       npix          I    # pixels in a line
        int       iscn          I    scan to read
        int       ndet          I    # detector lines in a scan
        float **  rhos          O    a scan's worth of reflectance

    Note that the reflectance in the file can be at wavelengths other than 
    the VIIRS wavelengths - the final reflectances will be interpolated to
    VIIRS

    Modification history:

    W. Robinson, SAIC  18 Nov 2009  Original development
    W. Robinson, SAIC  12 Mar 2010  modify for scan read

----------------------------------------------------------------------------*/ {
    static float *lrhos, *out_rhos, lam_lst[30], viirs_bands[16];
    static int first = 0, l_nlam;
    float vec_rho[20], vec_lrho[20];
    int ibnd, ipix, ilin, bndlist[MAX_BND], start[2], count[2];
    int rd_rhos_open(char *, float *, int *);
    /*
     *  for null file, clean up and exit
     */
    if (file == NULL) {
        /*  de-allocate storage, close file ids  */
        free(out_rhos);
        free(lrhos);

        for (ibnd = 0; ibnd < l_nlam; ibnd++)
            h5io_close(dat_ids + ibnd);
        h5io_close(&file_id);
        free(dat_ids);

        return 0;
    }
    /*
     *  initial time in, open the reflectance file
     */
    if (first == 0) {
        first = 1;
        /*  open file and get the band list  */
        if (rd_rhos_open(file, lam_lst, &l_nlam) != 0)
            return 1;
        /*  allocate rhos space for output */
        if ((out_rhos = (float *)
                malloc(npix * ndet * N_VNIR_BND * sizeof ( float))) == NULL) {
            printf("%s, %d: Unable to allocate rhos output storage\n",
                    __FILE__, __LINE__);
            return 1;
        }
        *rhos = out_rhos;
        /*  allocate local rhos space  */
        if ((lrhos = (float *)
                malloc(npix * ndet * l_nlam * sizeof ( float))) == NULL) {
            printf("%s, %d: Unable to allocate rhos storage\n",
                    __FILE__, __LINE__);
            return 1;
        }
        /* get the viirs band list and keep in local storage */
        if ((bnd_ix_2_sen_info("Lambda", (void *) bndlist)) < 0) {
            printf("%s, %d: failure to read sensor information\n",
                    __FILE__, __LINE__);
        }
        for (ibnd = 0; ibnd < N_VNIR_BND; ibnd++)
            *(viirs_bands + ibnd) = (float) *(bndlist + ibnd);
    }
    /*
     *  read the line's reflectance
     */
    ilin = iscn * ndet;
    start[0] = ilin;
    count[0] = ndet;
    start[1] = 0;
    count[1] = npix;
    for (ibnd = 0; ibnd < l_nlam; ibnd++) {
        if (h5io_rd_ds_slice(&dat_ids[ibnd], start, count,
                (void *) (lrhos + npix * ndet * ibnd)) != 0) {
            printf(
                    "%s, %d: Failed to read reflectance data set slice for line %d\n",
                    __FILE__, __LINE__, ilin);
            return 1;
        }
    }
    /*
     *  interpolate the values to the VIIRS bands
     */
    for (ilin = 0; ilin < ndet; ilin++) {
        for (ipix = 0; ipix < npix; ipix++) {
            for (ibnd = 0; ibnd < l_nlam; ibnd++)
                vec_lrho[ibnd] = *(lrhos + ipix + (ilin + ibnd * ndet) * npix);
            lspline(lam_lst, vec_lrho, l_nlam, viirs_bands, vec_rho, N_VNIR_BND);
            for (ibnd = 0; ibnd < N_VNIR_BND; ibnd++)
                *(out_rhos + ipix + (ilin + ibnd * ndet) * npix) = vec_rho[ibnd];
        }
    }
    return 0;
}

int rd_rhos_open(char *file, float *lam_lst, int *nlam)
/*-----------------------------------------------------------------------------
    Program:   rd_rhos_open

    Description:  Open the reflectance file

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        char *    file          I    file of reflectance
        float *   lam_lst       O    list of wavelengths in the file
        int *     nlam          O    # of wavelengths 

    Modification history:

    W. Robinson, SAIC  18 Nov 2009  Original development

----------------------------------------------------------------------------*/ {
    char *sloc, parms[20][15], parm_list[1000];
    int more_chars, ibnd;
    /*
     *  open the reflectance file
     */
    if (h5io_openr(file, 0, &file_id) != 0) {
        printf("%s, %d: Unable to open the reflectance file: %s\n", __FILE__,
                __LINE__, file);
        return 1;
    }
    /*
     * get the parm names, wavelengths, and #
     */
    if (h5io_rd_attr(&file_id, "parm_list", (void *) parm_list) != 0) {
        printf("%s, %d: Unable to read parm_list from reflectance file: %s\n",
                __FILE__, __LINE__, file);
        return 1;
    }

    more_chars = 1;
    *nlam = 0;
    while (more_chars) {
        if (*nlam == 0)
            sloc = strtok(parm_list, ",");
        else
            sloc = strtok(NULL, ",");
        if (sloc != NULL) {
            strcpy(parms[*nlam], sloc);
            sscanf(parms[*nlam], "rhos_%f", lam_lst + *nlam);
            (*nlam)++;
        } else
            more_chars = 0;
    }
    /*
     * open datasets for each parm
     */
    if ((dat_ids = (h5io_str *) malloc(*nlam * sizeof ( h5io_str)))
            == NULL) {
        printf("%s, %d: Unable to reserve h5io ids for refl datasets\n",
                __FILE__, __LINE__);
        return 1;
    }
    for (ibnd = 0; ibnd < *nlam; ibnd++) {
        if (h5io_set_ds(&file_id, parms[ibnd], (dat_ids + ibnd)) != 0) {
            printf("%s, %d: Unable to set to dataset %d in refl file %s\n",
                    __FILE__, __LINE__, ibnd, file);
            return 1;
        }
    }
    return 0;
}
