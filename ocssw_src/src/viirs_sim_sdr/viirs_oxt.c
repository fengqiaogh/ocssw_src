#include "viirs_sim_sdr.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#define N_INT_RNG 16  /* # of inter-band ranges */
static float *ib_field[N_INT_RNG];

int viirs_oxt(ctl_struc *ctl, in_rec_struc *in_rec)
/*-----------------------------------------------------------------------------
    Program:   viirs_oxt.c

    Description:  Add the VIIRS optical crosstalk artifact to the simulated 
      radiance field

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        ctl_struc *  ctl        I    input controls
        in_rec_struc * in_rec    I/O  controls for input record reading

    Modification history:

    W. Robinson, SAIC  26 May 2010  Original development

----------------------------------------------------------------------------*/ {
    static int oxt_init = 0, n_vis_asg = 7;
    static int vis_asg[] = {1, 3, 5, 7, 10, 12, 14}, ib_status[N_INT_RNG];
    static int int_low_loc[N_INT_RNG], int_hi_loc[N_INT_RNG];
    static float cw_wav[N_INT_RNG], av_rad[N_INT_RNG];
    int irng, last_ib2, last_ib, iend, ib_targ, ib_ct, ib_pair[2];
    int in_loc, nval;
    int loc1, loc2, ilin, ipix;
    float rad0, rad1, rad2, lam0, lam1, lam2, r_bnd1, r_bnd2, rfill;

    /*
     *  do initialization
     */
    if (oxt_init == 0) {
        oxt_init = 1;

        printf("OXT parms:\n");
        printf("oxt_mode = %d\n", ctl->oxt_mode);
        printf("oxt_coef = %s\n", ctl->oxt_coef);
        printf("inter_band = %s\n", ctl->inter_band);
        /*
         *  read in the inter-band TOA radiance information
         */
        if (viirs_oxt_ib_read(ctl->inter_band, cw_wav, av_rad) != 0)
            return 1;
        /*
         *  set up the inter-band storage area
         */
        for (irng = 0; irng < N_INT_RNG; irng++) {
            if ((ib_field[irng] = (float *)
                    malloc(in_rec->npix * NDET * sizeof (float))) == NULL) {
                printf("%s, %d: unable to allocate interband field storage\n",
                        __FILE__, __LINE__);
                return 1;
            }
        }
        /*
         *  set up the inter-band status (0-unfilled, 1 filled) and the low and 
         *  high bands to use in the interpolation
         */
        for (irng = 0; irng < N_INT_RNG; irng++)
            ib_status[irng] = 0;
        for (irng = 0; irng < n_vis_asg; irng++)
            ib_status[ vis_asg[irng] ] = 1;
        /*
         *  last_ib and last_ib2 will retain locations of up to 2 of the previous 
         *  inter-bands that are filled and are updated as loop is run.  ib_pair
         *  will contain up to 2 next inter-bands that are filled.  These are used
         *  to set the best inter-bands to interpolate with in filling the 
         *  empty inter-bands
         */
        last_ib2 = -1;
        last_ib = -1;
        for (irng = 0; irng < N_INT_RNG; irng++) {
            if (ib_status[irng] == 0) {
                /*
                 *  get the next filled bands
                 */
                iend = 0;
                ib_targ = irng;
                ib_ct = 0;
                while (iend != 1) {
                    ib_targ++;
                    if (ib_targ == N_INT_RNG)
                        iend = 1;
                    else
                        if (ib_status[ib_targ] == 1) ib_pair[ib_ct++] = ib_targ;
                    if (ib_ct == 2) iend = 1;
                }
                /*
                 *  use the previously and next filled bands to get best 
                 *  interpolation pair
                 */
                if ((last_ib != -1) && (ib_ct >= 1)) { /* 2 surrounding locations, next lower and higher */
                    ib_pair[1] = ib_pair[0];
                    ib_pair[0] = last_ib;
                } else if ((ib_ct == 0) && (last_ib2 != -1)) { /* previous 2 locs */
                    ib_pair[0] = last_ib2;
                    ib_pair[1] = last_ib;
                } else if ((ib_ct == 2) && (last_ib == -1)) { /*  next 2 locations they are already in ib_pair */
                    ;
                } else { /* insufficient points to do the fill - should not happen */
                    printf("%s, %d: insufficient data to interpolate\n",
                            __FILE__, __LINE__);
                    return 1;
                }
                int_low_loc[irng] = ib_pair[0];
                int_hi_loc[irng] = ib_pair[1];
            } else { /* filled band found, maintain previous filled band locations */
                last_ib2 = last_ib;
                last_ib = irng;
            }
        }
    }
    /*
     *  insert the band radiances into the inter-band storage
     *  in_loc is any offset due to a line margin which te storage dosen't need
     */
    in_loc = in_rec->margin[0] * in_rec->npix;
    nval = sizeof (float) * in_rec->npix * NDET;
    for (irng = 0; irng < n_vis_asg; irng++) {
        memcpy(ib_field[vis_asg[irng]], in_rec->bnd_lt[irng] + in_loc, nval);
    }
    /*
     *  fill the empty inter-band values with an interpolation of surrounding
     *  radiances, conditioned by the expected TOA spectrum
     */
    for (ilin = 0; ilin < NDET; ilin++)
        for (ipix = 0; ipix < in_rec->npix; ipix++)
            for (irng = 0; irng < N_INT_RNG; irng++) {
                if (ib_status[irng] == 0) {
                    loc1 = int_low_loc[irng];
                    loc2 = int_hi_loc[irng];
                    lam1 = cw_wav[loc1];
                    lam2 = cw_wav[loc2];
                    lam0 = cw_wav[irng];
                    rad1 = av_rad[loc1];
                    rad2 = av_rad[loc2];
                    rad0 = av_rad[irng];
                    r_bnd1 = *(ib_field[loc1] + ipix + ilin * in_rec->npix);
                    r_bnd2 = *(ib_field[loc2] + ipix + ilin * in_rec->npix);
                    rfill = rad0 / (lam2 - lam1) *
                            (r_bnd1 * (lam2 - lam0) / rad1 +
                            r_bnd2 * (lam0 - lam1) / rad2);
                    *(ib_field[irng] + ipix + ilin * in_rec->npix) = rfill;
                }
            }
    /*
     *  Get and apply the coefficients to proper detector in the inter-band field
     *  and modify the radiances with it
     */
    if (viirs_oxt_comp(ctl, in_rec) != 0)
        return 1;
    /*
    printf( "PREMATURE END\n" );
    return 1;
     */
    return 0;
}

int viirs_oxt_ib_read(char *ib_file, float *cw_wav, float *av_rad)
/*-----------------------------------------------------------------------------
    routine:   viirs_oxt_read

    Description:  set up the information needed for the inter-band radiance
      determination

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        char *    ib_file       I    file containing the inter-band radiance
                                     information for typical TOA radiances
                                     in VIIRS iinterband wavelength ranges
        float *   cw_wav        O    radiance-weighted center wavelengths 
                                     of the TOA spectrum
        float *   av_rad        O    average radiance in the inter-bands

    Modification history:

    W. Robinson, SAIC  26 May 2010  Original development

----------------------------------------------------------------------------*/ {
    FILE *fp;
    int items_set[N_INT_RNG * 2], ipar, ip1, ip2;
    char line[80], name[80], value[80], parm[80];
    char *p, *p1, *p2;
    char *ib_par_nm[] = {"CW_WAVE", "AV_RAD"};
    char *ib_par_rng[] = {"[IB01]", "[M1]", "[IB12]", "[M2]", "[IB23]", "[M3]",
        "[IB34]", "[M4]", "[IB4I1]", "[I1]", "[M5]", "[IB56]", "[M6]", "[IB67]",
        "[M7]", "[IB70]"};
    float fval;
    for (ipar = 0; ipar < N_INT_RNG * 2; ipar++)
        items_set[ipar] = 0;
    /*
     *  open the file
     */
    if ((fp = fopen(ib_file, "r")) == NULL) {
        printf("%s, line %d: unable to open inter-band file: %s\n",
                __FILE__, __LINE__, ib_file);
        return 1;
    }
    /*
     *  Loop through parameters to define
     */
    while (fgets(line, 80, fp)) {
        memset(name, '\0', sizeof ( name));
        memset(value, '\0', sizeof ( value));
        /*
         * Skip comment lines, empty lines, and lines without a
         * name = value pair.
         */
        if (line[0] == '#' || line[0] == '\n')
            continue;
        if (!(p = strchr(line, '=')))
            continue;

        /*
         * Parse parameter name string
         */
        p1 = line;
        while (isspace(*p1))
            p1++;
        p2 = p - 1;
        while (isspace(*p2))
            p2--;
        strncpy(name, p1, p2 - p1 + 1);

        /*
         * Parse parameter value string
         */
        p1 = p + 1;
        while (isspace(*p1))
            p1++;
        p2 = p1;
        while (!isspace(*p2))
            p2++;
        strncpy(value, p1, p2 - p1);
        printf("int band values are:  %s = %s\n", name, value);

        /*
         *  for each parameter, set the array value
         */
        for (ipar = 0; ipar < N_INT_RNG * 2; ipar++) {
            ip1 = 1;
            if (ipar < N_INT_RNG)
                ip1 = 0;
            ip2 = ipar - ip1 * N_INT_RNG;
            sprintf(parm, "%s%s", ib_par_nm[ip1], ib_par_rng[ip2]);
            /*
             *  
             */
            if (strcmp(name, parm) == 0) {
                items_set[ipar] = 1;
                fval = (float) atof(value);
                if (ip1 == 0)
                    cw_wav[ip2] = fval;
                else
                    av_rad[ip2] = fval;
                break;
            }
        }
    }
    /*
     *  make sure all values were read from the file
     */
    for (ipar = 0; ipar < N_INT_RNG * 2; ipar++)
        if (items_set[ipar] == 0) {
            printf("%s, %d: An item in the file %s was not set\n",
                    __FILE__, __LINE__, ib_file);
            return 1;
        }
    fclose(fp);
    return 0;
}

int viirs_oxt_comp(ctl_struc *ctl, in_rec_struc *in_rec)
/*-----------------------------------------------------------------------------
    routine:   viirs_oxt_comp

    Description:  compute the optical crosstalk artifact and modify the 
      radiance field with it

   Returns type: int - 0 if good

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        ctl_struc *  ctl        I    input controls
        in_rec_struc * in_rec    I/O  controls for input record reading

    Note that ib_field - the field of inter-band radiances, is common 
    with this routine and the viirs_oxt routine

    Modification history:

    W. Robinson, SAIC  1 Jun 2010  Original development

----------------------------------------------------------------------------*/ {
    static int ifirst = 0, n_rec_bnd, n_rec_det, n_int_rng, n_snd_rng, n_snd_det;
    /* registration offset for the re-ordered sender bands:
       (where is band M? relative to M1?)
                             M1 M2  M3  M4  M5  M6  M7   I1  I2   */
    static int snd_off[] = {0, -3, -9, -6, -18, -21, -15, -12, -14};
    static float *oxt_coeff;
    h5io_str oxt_fid;
    int irlin, irpix, irbnd, ilam, isbnd, isdet, reg_off, isoff;
    float mod_rad;
    /*
     *  initialization will read in the crosstalk coefficients
     */
    if (ifirst == 0) {
        ifirst = 1;
        /*
         *  open the optical crosstalk coefficient file
         */
        if (h5io_openr(ctl->oxt_coef, 0, &oxt_fid) != 0) {
            printf("%s, %d - could not open HDF 5 Optical Xtalk file: %s\n",
                    __FILE__, __LINE__, ctl->oxt_coef);
            return 1;
        }
        /*
         *  read array sizes:
         */
        if (h5io_rd_attr(&oxt_fid, "number of receiver bands", &n_rec_bnd) != 0) {
            printf("%s, %d - failed to read Optical Xtalk attrib n_rec_bnd\n",
                    __FILE__, __LINE__);
            return 1;
        }
        if (h5io_rd_attr(&oxt_fid, "number of receiver detectors", &n_rec_det)
                != 0) {
            printf("%s, %d - failed to read Optical Xtalk attrib n_rec_det\n",
                    __FILE__, __LINE__);
            return 1;
        }
        if (h5io_rd_attr(&oxt_fid, "number of inter-band ranges", &n_int_rng)
                != 0) {
            printf("%s, %d - failed to read Optical Xtalk attrib n_int_rng\n",
                    __FILE__, __LINE__);
            return 1;
        }
        if (h5io_rd_attr(&oxt_fid, "number of sender bands", &n_snd_rng) != 0) {
            printf("%s, %d - failed to read Optical Xtalk attrib n_snd_rng\n",
                    __FILE__, __LINE__);
            return 1;
        }
        if (h5io_rd_attr(&oxt_fid, "number of sender detectors", &n_snd_det)
                != 0) {
            printf("%s, %d - failed to read Optical Xtalk attrib n_snd_det\n",
                    __FILE__, __LINE__);
            return 1;
        }
        printf("%s, %d: optical crosstalk array dimensions\n",
                __FILE__, __LINE__);
        printf("    n_rec_bnd = %d,  n_rec_det = %d, n_int_rng = %d\n",
                n_rec_bnd, n_rec_det, n_int_rng);
        printf("    n_snd_rng = %d, n_snd_det = %d\n", n_snd_rng, n_snd_det);
        /*
         *  make storage for the coefficients and read them in
         */
        if ((oxt_coeff = malloc(n_rec_bnd * n_rec_det * n_int_rng *
                n_snd_rng * n_snd_det * sizeof ( float))) == NULL) {
            printf("%s, %d: Unable to allocate optical crosstalk coeff storage\n",
                    __FILE__, __LINE__);
            return 1;
        }

        if (h5io_grab_ds(&oxt_fid, "optical crosstalk coefficients",
                (void *) oxt_coeff) != 0) {
            printf("%s, %d: Unable to read crosstalk coefficients\n",
                    __FILE__, __LINE__);
            return 1;
        }
        /*
         *  and finish with the file
         */
        if (h5io_close(&oxt_fid) != 0) {
            printf("%s, %d: Unable to close crosstalk coefficient file\n",
                    __FILE__, __LINE__);
            return 1;
        }
    }
    /*
     *  Apply the coefficients to the inter-band field and modify the radiances
     *  with the artifact
     */
    for (irlin = 0; irlin < NDET; irlin++) {
        for (irpix = 0; irpix < in_rec->npix; irpix++) {
            for (irbnd = 0; irbnd < N_VNIR_BND; irbnd++) {
                mod_rad = 0.;
                for (isdet = 0; isdet < NDET; isdet++) {
                    for (isbnd = 0; isbnd < 9; isbnd++) {
                        /*
                         *  isoff is the offset to the sender detector, including the 
                         *  registration correction required
                         *  the registration offset of the sender relative to M1, snd_off
                         *  can be shown below
                         *  sender  M1   M2   M3   M4   M5   M6   M7   I1   I2
                         *  order    0    1    2    3    4    5    6    7    8
                         *  snd_off  0,  -3,  -9,  -6, -18, -21, -15, -12, -14
                         */
                        reg_off = snd_off[ isbnd ] - snd_off[ irbnd ];
                        /*
                         *  note that if data beyond margin is needed, skip including it
                         */
                        if ((reg_off + irpix >= 0) && (reg_off + irpix < in_rec->npix)) {
                            isoff = isdet * in_rec->npix + irpix + reg_off;
                            for (ilam = 0; ilam < 16; ilam++) {
                                /*
                                 *  oxt_coeff[ irbnd, irlin, ilam, isbnd, islin ]
                                 */
                                mod_rad +=
                                        *(oxt_coeff + irbnd + N_VNIR_BND * (irlin + NDET *
                                        (ilam + 16 * (isbnd + 9 * isdet)))) *
                                        *(ib_field[ilam] + isoff);
                            }
                        }
                    }
                }
                /*
                 *  add the artifact to the radiance
                 */
                *(in_rec->bnd_lt[irbnd] + irpix +
                        in_rec->npix * (irlin + in_rec->margin[0])) += mod_rad;
            }
        }
    }
    /*
     *  scan modified, return
     */
    return 0;
}
