#include "viirs_sim_sdr.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int viirs_noise(ctl_struc *ctl, in_rec_struc *in_rec, int noise_stop)
/*-----------------------------------------------------------------------------
    Program:   viirs_noise.c

    Description:  Add the VIIRS noise artifact to the simulated 
      radiance field

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        ctl_struc *  ctl        I    input controls
        in_rec_struc * in_rec    I/O  controls for input record reading
        int       noise_stop    I    if 1, shut down the noise generation
                                     do at program end

    Modification history:

    W. Robinson, SAIC  25 May 2011  Original development

----------------------------------------------------------------------------*/ {
    int npix, nlin, ibnd, ilin, ipix;
    double rad, noise, noise_sd, maxr, rad_lim;
    const gsl_rng_type * T;
    static gsl_rng * r;
    static int noise_init = 0;
    static float a0[N_VNIR_BND], a1[N_VNIR_BND], a2[N_VNIR_BND];
    static float maxrad[N_VNIR_BND];

    /*
     *  close down if that switch is set
     */
    if (noise_stop == 1) {
        gsl_rng_free(r);
    } else {
        /*
         *  do initialization
         */
        if (noise_init == 0) {
            noise_init = 1;

            printf("%s, %d: Entering routine viirs_noise\n", __FILE__, __LINE__);
            printf("noise_coef = %s\n", ctl->noise_coef);
            /*
             *  read in the coefficients for the SNR for the bands
             */
            if (viirs_noise_coef_rd(ctl->noise_coef, a0, a1, a2, maxrad) != 0)
                return 1;
            /*
             *  set up for the random numbers needed
             */
            gsl_rng_env_setup();
            T = gsl_rng_default;
            r = gsl_rng_alloc(T);
        }
        /*
         *  for each pixel (in a line in a band) derive and apply the noise
         *  except for bad values
         */
        npix = in_rec->npix;
        nlin = in_rec->ndet_scan;
        for (ibnd = 0; ibnd < N_VNIR_BND; ibnd++) {
            maxr = *(maxrad + ibnd);
            for (ilin = 0; ilin < nlin; ilin++) {
                for (ipix = 0; ipix < npix; ipix++) {
                    if (*(in_rec->bnd_q[ibnd] + ilin * npix + ipix) == 0) {
                        rad = *(in_rec->bnd_lt[ibnd] + ilin * npix + ipix);
                        /*
                         *  limit the radiance used in the noise computation
                         *  to the maxrad for that band
                         */
                        rad_lim = (rad > maxr) ? maxr : rad;

                        noise_sd = rad_lim / (*(a0 + ibnd) + *(a1 + ibnd) * rad_lim +
                                *(a2 + ibnd) * rad_lim * rad_lim);
                        noise = gsl_ran_gaussian(r, noise_sd);
                        *(in_rec->bnd_lt[ibnd] + ilin * npix + ipix) = rad + noise;
                    }
                } /* end pix loop */
            } /* end line loop */
        } /* end band loop */
    }
    /*
     *  and end
     */
    return 0;
}

int viirs_noise_coef_rd(char *file, float *a0, float *a1, float *a2,
        float *maxrad)
/*-----------------------------------------------------------------------------
    routine:   viirs_noise_coef_rd

    Description:  read quadratic coefficients describing the SNR as a 
        f(radiance)

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        char *    file          I    text file containing the SNR coefficients
        float *   a0            O    constant term of SNR for the VIS/NIR bands
        float *   a1            O    linear term of SNR
        float *   a2            O    quadratic term of SNR
        float *   maxrad        O    max radiance to use in SNR computation

    SNR = a0[band] + a1[band] * rad + a2[band]  rad^2
    with rad = TOA radiance || maxrad, whichever is smallest

    Modification history:

    W. Robinson, SAIC  25 May 2011  Original development

----------------------------------------------------------------------------*/ {
    FILE *fp;
    int items_set[N_VNIR_BND * 4], ipar, iord, ibnd, ibnd_no, val_set;
    char line[80], name[80], value[80], parm[80];
    char *p, *p1, *p2;
    float fval;

    for (ipar = 0; ipar < N_VNIR_BND * 4; ipar++)
        items_set[ipar] = 0;
    /*
     *  open the file
     */
    if ((fp = fopen(file, "r")) == NULL) {
        printf("%s, line %d: unable to open noise coefficient file: %s\n",
                __FILE__, __LINE__, file);
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
        printf("noise coeff values are:  %s = %s\n", name, value);

        /*
         *  for each parameter, set the array value
         */
        val_set = 0;
        for (ibnd = 0; ibnd < N_VNIR_BND; ibnd++) {
            ibnd_no = ibnd + 1;
            for (iord = 0; iord < 3; iord++) {
                sprintf(parm, "a%1.1d[%1.1d]", iord, ibnd_no);

                if (strcmp(name, parm) == 0) {
                    ipar = ibnd * 4 + iord;
                    items_set[ipar] = 1;
                    val_set = 1;
                    fval = (float) atof(value);
                    if (iord == 0)
                        *(a0 + ibnd) = fval;
                    else if (iord == 1)
                        *(a1 + ibnd) = fval;
                    else if (iord == 2)
                        *(a2 + ibnd) = fval;
                    break;
                }
            }
            if (val_set == 1) break;

            sprintf(parm, "maxrad[%1.1d]", ibnd_no);
            if (strcmp(name, parm) == 0) {
                ipar = ibnd * 4 + 3;
                items_set[ipar] = 1;
                fval = (float) atof(value);
                *(maxrad + ibnd) = fval;
                break;
            }
        }
    }
    /*
     *  make sure all values were read from the file
     */
    for (ipar = 0; ipar < N_VNIR_BND * 4; ipar++)
        if (items_set[ipar] == 0) {
            printf("%s, %d: An item in the file %s was not set\n",
                    __FILE__, __LINE__, file);
            return 1;
        }
    fclose(fp);
    return 0;
}
