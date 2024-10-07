#include "viirs_sim_sdr.h"
/*
 *  as the hdf 4 will be hopefully temporary, just have the includes here
 */
#include "hio.h"

#include <hdf.h>
#include <mfhdf.h>

/*
 * file vset_cal.c containing 
 *   vset_cal_gain  for the gain information
 *   vset_cal_rvs  for the RVS info
 */

int vset_cal_gain(char *file, vir_gain_struc *gain)
/*-----------------------------------------------------------------------------
    Routine:   vset_cal_gain

    Description:  get the gain related coefficients for VIIRS calibration
      
    Returns type: int - 0 if good, else I/O error

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        char *    file          I    calibration file name
        vir_gain_struc * gain   O    structure with quadratic gain coeffs,
                                     dark count, and max Lt for the high gain
                                     (to cue shifting to use of low gain)

    Modification history:

    W. Robinson, SAIC  16 Sep 2010  Original development

----------------------------------------------------------------------------*/ {
    float dn_trans_pct = 0.83; /* the gain switches to low at 83% of the 4095 
                                counts for the 12 bit dn (small) values  */
    float c0, c1, c2, dark, lt_max, dn, big_dn, dn_max = 4095.;
    float min_hg_lim[MAX_BND], max_hg_lim[MAX_BND],
            min_lg_lim[MAX_BND], max_lg_lim[MAX_BND],
            min_dn_lg[MAX_BND], max_dn_lg[MAX_BND],
            min_dn_hg[MAX_BND], max_dn_hg[MAX_BND],
            min_big_dn_lg[MAX_BND], max_big_dn_lg[MAX_BND],
            min_big_dn_hg[MAX_BND], max_big_dn_hg[MAX_BND];
    int iham, idet, igain, ibnd, ext_gain = 0;
    int nham = 2, ndet = 16, ngain = 2, nbnd = 11;
    char gain_ct[] = {2, 2, 2, 2, 2, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1};
    /* for above, # gain states for each m band */
    h5io_str fid;
    /*
     *  If file is unspecified, fill general array with values from the TVAC test 
     *  which are band-averaged.  They are also just the c1 as c0 and c2 are 0.
     *  The file source is
     *  VIIRS_RSB_Band_Average_Gains.xlsx, other info:
     *  From NICST Murphy Chart RC_02_V3 (1-12-2010)
     *  Units are dn/(W/m^2/str/micron)
     *  SO to get the c1 value (W.../dn) we'll need to invert these
     *  and if c2 is in W^2/dn in future, we'll have to re-form the equations
     *  11 reflective M bands, low gain followed by high gain, 
     *  0. for bands with only high gain
     */
    float m_gains[22] = {6.4, 4.90, 4.44, 4.76, 5.15, -1., 9.97,
        -1., -1., -1., -1.,
        27.26, 23.84, 27.35, 37.13, 50.97, 109.37, 84.19,
        28.11, 41.99, 49.52, 91.49};

    /*
     *  if a file was specified, open it and get the size information
     */
    if (strcmp(file, "Unspecified") != 0) {
        ext_gain = 1;
        if (h5io_openr(file, 0, &fid) != 0) {
            printf("%s, %d, Unable to open gain calibration file: \n%s\n",
                    __FILE__, __LINE__, file);
            return 1;
        }
        if (h5io_rd_attr(&fid, "ngain", (void *) &(gain->ngain)) != 0) {
            printf("%s, %d, Unable to read ngain attr from calibration file: \n%s\n",
                    __FILE__, __LINE__, file);
            return 1;
        }
        if (h5io_rd_attr(&fid, "nham", (void *) &(gain->nham)) != 0) {
            printf("%s, %d, Unable to read nham attr from calibration file: \n%s\n",
                    __FILE__, __LINE__, file);
            return 1;
        }
        if (h5io_rd_attr(&fid, "ndet", (void *) &(gain->ndet)) != 0) {
            printf("%s, %d, Unable to read ndet attr from calibration file: \n%s\n",
                    __FILE__, __LINE__, file);
            return 1;
        }
        if (h5io_rd_attr(&fid, "nbnd", (void *) &(gain->nbnd)) != 0) {
            printf("%s, %d, Unable to read nbnd attr from calibration file: \n%s\n",
                    __FILE__, __LINE__, file);
            return 1;
        }
        if ((gain->ngain != ngain) || (gain->nham != nham) ||
                (gain->ndet != ndet) || (gain->nbnd != nbnd)) {
            printf("%s, %d, Unexpected gain array dimensions found in file: \n%s\n",
                    __FILE__, __LINE__, file);
            printf("   ngain: %d, expected: %d\n", gain->ngain, ngain);
            printf("   nham: %d, expected: %d\n", gain->nham, nham);
            printf("   ndet: %d, expected: %d\n", gain->ndet, ndet);
            printf("   nbnd: %d, expected: %d\n", gain->nbnd, nbnd);
            return 1;
        }
    } else {
        gain->nham = nham;
        gain->ndet = ndet;
        gain->ngain = ngain;
        gain->nbnd = nbnd;
    }
    /*
     *  allocate space for the coefficient arrays
     */
    if ((gain->c0 = (float *)
            malloc(nham * ndet * ngain * nbnd * sizeof ( float))) == NULL) {
        printf("%s, %d: Error allocating the gain->c0\n", __FILE__, __LINE__);
        return 1;
    }
    if ((gain->c1 = (float *)
            malloc(nham * ndet * ngain * nbnd * sizeof ( float))) == NULL) {
        printf("%s, %d: Error allocating the gain->c1\n", __FILE__, __LINE__);
        return 1;
    }
    if ((gain->c2 = (float *)
            malloc(nham * ndet * ngain * nbnd * sizeof ( float))) == NULL) {
        printf("%s, %d: Error allocating the gain->c2\n", __FILE__, __LINE__);
        return 1;
    }
    if ((gain->dark = (float *)
            malloc(nham * ndet * ngain * nbnd * sizeof ( float))) == NULL) {
        printf("%s, %d: Error allocating the gain->dark\n", __FILE__, __LINE__);
        return 1;
    }
    /*
     *  set the min, max low and high gain limits radiances (without RVS)
     *  As the hi-lo gain transition happens at a % of dn, get that instead for h-l
     */
    for (ibnd = 0; ibnd < nbnd; ibnd++) {
        min_hg_lim[ibnd] = 999.;
        max_hg_lim[ibnd] = -999.;
        min_lg_lim[ibnd] = 999.;
        max_lg_lim[ibnd] = -999.;
        min_dn_lg[ibnd] = 9999.;
        max_dn_lg[ibnd] = -9999.;
        min_dn_hg[ibnd] = 9999.;
        max_dn_hg[ibnd] = -9999.;
        min_big_dn_lg[ibnd] = 9999.;
        max_big_dn_lg[ibnd] = -9999.;
        min_big_dn_hg[ibnd] = 9999.;
        max_big_dn_hg[ibnd] = -9999.;
    }
    dn = dn_max * dn_trans_pct;
    /*
     * If no gain file specified, use the internal values, else read in
     */
    if (ext_gain == 1) {
        if (h5io_grab_ds(&fid, "c0", (void *) gain->c0) != 0) {
            printf("%s, %d, Unable to read c0 dataset from calibration file: \n%s\n",
                    __FILE__, __LINE__, file);
            return 1;
        }
        if (h5io_grab_ds(&fid, "c1", (void *) gain->c1) != 0) {
            printf("%s, %d, Unable to read c1 dataset from calibration file: \n%s\n",
                    __FILE__, __LINE__, file);
            return 1;
        }
        if (h5io_grab_ds(&fid, "c2", (void *) gain->c2) != 0) {
            printf("%s, %d, Unable to read c2 dataset from calibration file: \n%s\n",
                    __FILE__, __LINE__, file);
            return 1;
        }
        if (h5io_grab_ds(&fid, "dark", (void *) gain->dark) != 0) {
            printf(
                    "%s, %d, Unable to read dark dataset from calibration file: \n%s\n",
                    __FILE__, __LINE__, file);
            return 1;
        }
    } else {
        /*
         *  fill the coefficients from the local information
         */
        for (iham = 0; iham < nham; iham++) {
            for (idet = 0; idet < ndet; idet++) {
                for (ibnd = 0; ibnd < nbnd; ibnd++) {
                    for (igain = 0; igain < ngain; igain++) {
                        if ((gain_ct[ibnd] == 1) && (igain == 0)) continue;
                        /*  set c0, 1, 2, dark and gain from internal values */
                        c0 = 0.;
                        /* note the provided c1s are dn / W..., so invert to use  */
                        c1 = 1. / *(m_gains + ibnd + igain * nbnd);
                        c2 = 0.;
                        dark = 250. + ibnd;
                        *(gain->c0 + ibnd + nbnd * (igain + ngain *
                                (idet + ndet * iham))) = c0;
                        *(gain->c1 + ibnd + nbnd * (igain + ngain *
                                (idet + ndet * iham))) = c1;
                        *(gain->c2 + ibnd + nbnd * (igain + ngain *
                                (idet + ndet * iham))) = c2;
                        *(gain->dark + ibnd + nbnd * (igain + ngain *
                                (idet + ndet * iham))) = dark;
                    }
                }
            }
        }
    }
    /*
     *  for gathering the min, max of the limit radiance, compute here
     */
    for (iham = 0; iham < nham; iham++) {
        for (idet = 0; idet < ndet; idet++) {
            for (ibnd = 0; ibnd < nbnd; ibnd++) {
                for (igain = 0; igain < ngain; igain++) {
                    if ((gain_ct[ibnd] == 1) && (igain == 0)) continue;
                    big_dn = dn + *(gain->dark + ibnd + nbnd * (igain + ngain *
                            (idet + ndet * iham)));
                    c0 = *(gain->c0 + ibnd + nbnd * (igain + ngain *
                            (idet + ndet * iham)));
                    c1 = *(gain->c1 + ibnd + nbnd * (igain + ngain *
                            (idet + ndet * iham)));
                    c2 = *(gain->c2 + ibnd + nbnd * (igain + ngain *
                            (idet + ndet * iham)));
                    lt_max = c0 + c1 * dn + c2 * dn * dn;
                    if (igain == 0) {
                        if (lt_max < min_lg_lim[ibnd]) min_lg_lim[ibnd] = lt_max;
                        if (lt_max > max_lg_lim[ibnd]) max_lg_lim[ibnd] = lt_max;
                        if (dn < min_dn_lg[ibnd]) min_dn_lg[ibnd] = dn;
                        if (dn > max_dn_lg[ibnd]) max_dn_lg[ibnd] = dn;
                        if (big_dn < min_big_dn_lg[ibnd]) min_big_dn_lg[ibnd] = big_dn;
                        if (big_dn > max_big_dn_lg[ibnd]) max_big_dn_lg[ibnd] = big_dn;
                    } else {
                        if (lt_max < min_hg_lim[ibnd]) min_hg_lim[ibnd] = lt_max;
                        if (lt_max > max_hg_lim[ibnd]) max_hg_lim[ibnd] = lt_max;
                        if (dn < min_dn_hg[ibnd]) min_dn_hg[ibnd] = dn;
                        if (dn > max_dn_hg[ibnd]) max_dn_hg[ibnd] = dn;
                        if (big_dn < min_big_dn_hg[ibnd]) min_big_dn_hg[ibnd] = big_dn;
                        if (big_dn > max_big_dn_hg[ibnd]) max_big_dn_hg[ibnd] = big_dn;
                    }
                    if (lt_max < 0) {
                        printf("%s, %d: lt_max < 0 at ham: %d, det: %d, bnd: %d, gain:%d\n",
                                __FILE__, __LINE__, iham, idet, ibnd, igain);
                    }
                }
            }
        }
    }
    /*
     *  report the radiance, dn limits found
     */
    printf("\n\n%s, %d: info:\n", __FILE__, __LINE__);
    printf("\nRadiance, dn limit min and max per band\n");
    printf("(file: %s)\n", file);
    if (strcmp(file, "Unspecified") == 0)
        printf("Using the general gains from TVAC tests\n");
    printf("For HIGH GAIN setting\n");
    printf(
            "\nM band    Min rad    Max rad    Min dn    Max dn    Min DN    Max DN\n");
    for (ibnd = 0; ibnd < nbnd; ibnd++)
        printf("   %3d  %9.4f  %9.4f  %8.2f  %8.2f  %8.2f  %8.2f\n", ibnd + 1,
            min_hg_lim[ibnd], max_hg_lim[ibnd], min_dn_hg[ibnd], max_dn_hg[ibnd],
            min_big_dn_hg[ibnd], max_big_dn_hg[ibnd]);
    printf("\nFor LOW GAIN setting\n");
    printf(
            "\nM band    Min rad    Max rad    Min dn    Max dn    Min DN    Max DN\n");
    for (ibnd = 0; ibnd < nbnd; ibnd++)
        printf("   %3d  %9.4f  %9.4f  %8.2f  %8.2f  %8.2f  %8.2f\n", ibnd + 1,
            min_lg_lim[ibnd], max_lg_lim[ibnd], min_dn_lg[ibnd], max_dn_lg[ibnd],
            min_big_dn_lg[ibnd], max_big_dn_lg[ibnd]);
    printf("----------------------------------------\n\n");
    /*
     *  transfer the # gain states
     */
    for (ibnd = 0; ibnd < nbnd; ibnd++)
        *(gain->gain_ct + ibnd) = gain_ct[ibnd];
    return 0;
}

int vset_cal_rvs(char *file, vir_rvs_struc *rvs)
/*-----------------------------------------------------------------------------
    Routine:   vset_cal_rvs

    Description:  get the rvs related coefficients for VIIRS calibration

    Returns type: int - 0 if good, else I/O error

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        char *    file          I    RVS calibration file name up to the 
                                     band wavelength.  individual files
                                     are for each band's RVS
        vir_rvs_struc * rvs     O    structure with quadratic RVS coeffs,
                                     as a function of AOI on HAM 

    Modification history:

    W. Robinson, SAIC  16 Sep 2010  Original development
    W. Robinson, SAIC   7 Jan 2011  Change so that if 'Unspecified' file
                                    is used, a unity RVS is used

----------------------------------------------------------------------------*/ {
    int ibnd, iham, idet, bnd_lambda[16];
    int nbnd = 11, npoly = 3, nham = 2, ndet = 1;;
    char hfile[500];
    float *lcl_buf;
    hio_struct finfo;
    /*
     *  The coefficients for RVS are for a polinomial applied to the angle of
     *  incidence on the mirror AOI: RVS = a0 + a1 * aoi + a2 * aoi^2
     *
     * set up array sizes
     */
    rvs->nham = nham;
    rvs->ndet = ndet;
    rvs->nbnd = nbnd;
    /*
     *  There is a RVS for I1 at 640 nm, but for this version, I1 will be ignored.
     */
    if (bnd_ix_2_sen_info("Lambda", (void *) bnd_lambda) != 0) {
        printf("%s, %d: Error in getting lambda, Exiting\n", __FILE__, __LINE__);
        return 1;
    }
    /*
     *  although the hdf files could be read for the array size, as other 
     *  parts are hard-coded, hard-code the size as we know for these files 
     *
     *  allocate the rvs storage
     */
    if ((rvs->a0 = (float *) malloc(nham * ndet * nbnd * sizeof ( float)))
            == NULL) {
        printf("%s, %d: Error in allocating rvs->a0\n", __FILE__, __LINE__);
        return 1;
    }
    if ((rvs->a1 = (float *) malloc(nham * ndet * nbnd * sizeof ( float)))
            == NULL) {
        printf("%s, %d: Error in allocating rvs->a1\n", __FILE__, __LINE__);
        return 1;
    }
    if ((rvs->a2 = (float *) malloc(nham * ndet * nbnd * sizeof ( float)))
            == NULL) {
        printf("%s, %d: Error in allocating rvs->a2\n", __FILE__, __LINE__);
        return 1;
    }
    /*
     *  Print the input file 
     */
    printf("%s, %d, RVS file input: %s\n", __FILE__, __LINE__, file);
    /*
     *  if RVS file was unspecified, use a unit RVS
     */
    if (strcmp(file, "Unspecified") == 0) {
        printf("%s, %d: Note that unit RVS will be used\n", __FILE__, __LINE__);
        for (ibnd = 0; ibnd < nbnd; ibnd++)
            for (iham = 0; iham < nham; iham++)
                for (idet = 0; idet < ndet; idet++) {
                    *(rvs->a0 + ibnd + nbnd * (idet + ndet * iham)) = 1.;
                    *(rvs->a1 + ibnd + nbnd * (idet + ndet * iham)) = 0.;
                    *(rvs->a2 + ibnd + nbnd * (idet + ndet * iham)) = 0.;
                }
    } else {
        /*
         *  for a rvs file specified,
         *  make a read buffer for the hdf 4 i/o of the RVS
         */
        if ((lcl_buf = (float *) malloc(npoly * nham * ndet * sizeof ( float)))
                == NULL) {
            printf("%s, %d: Error in allocating lcl read buffer\n",
                    __FILE__, __LINE__);
            return 1;
        }
        /*
         *  open each hdf 4 file and read the gain 
         */
        for (ibnd = 0; ibnd < nbnd; ibnd++) {
            /*
             *  make file name
             */
            sprintf(hfile, "%s%s%d%s", file, "_", bnd_lambda[ibnd], ".hdf");

            printf("%s, %d, opening RVS file %s\n", __FILE__, __LINE__, hfile);
            /*
             *  open file
             */
            if (hio_open(hfile, DFACC_RDONLY, &finfo) != 0) {
                printf("%s, %d: Error, Unable to open the HDF file: %s\n",
                        __FILE__, __LINE__, hfile);
                return 1;
            }
            /*
             *  read the rvs
             */
            if (hio_r_sds(finfo, "rvs", DFNT_FLOAT32, (void *) lcl_buf) != 0) {
                printf("%s, %d: Error, Unable to read the 'rvs' SDS for HDF file: %s\n",
                        __FILE__, __LINE__, hfile);
                return 1;
            }
            /*
             *  read the aoi range
             */
            if (hio_r_sds(finfo, "aoirange", DFNT_FLOAT32, (void *) rvs->aoi_range)
                    != 0) {
                printf(
                        "%s, %d: Error, Unable to read the 'aoirange' SDS for HDF file: %s\n",
                        __FILE__, __LINE__, hfile);
                return 1;
            }
            printf("%s, %d: Info: aoi range is: %f  %f\n", __FILE__, __LINE__,
                    rvs->aoi_range[0], rvs->aoi_range[1]);
            /*
             *   close the file
             */
            hio_close(finfo);
            /*
             *  transfer the information into the rvs coeff storage
             */
            for (iham = 0; iham < nham; iham++) {
                for (idet = 0; idet < ndet; idet++) {
                    *(rvs->a0 + ibnd + nbnd * (idet + ndet * iham)) =
                            *(lcl_buf + idet + ndet * (iham + nham * 0));
                    *(rvs->a1 + ibnd + nbnd * (idet + ndet * iham)) =
                            *(lcl_buf + idet + ndet * (iham + nham * 1));
                    *(rvs->a2 + ibnd + nbnd * (idet + ndet * iham)) =
                            *(lcl_buf + idet + ndet * (iham + nham * 2));
                }
            }
        }
        /*
         * release local storage
         */
        free(lcl_buf);
    }
    return 0;
}
