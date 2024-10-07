#include "viirs_sim_sdr.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

int viirs_straylt(ctl_struc *ctl, in_rec_struc *in_rec, int stray_stop)
/*-----------------------------------------------------------------------------
    Program:   viirs_straylt.c

    Description:  Add the VIIRS stray light artifact to the simulated 
      radiance field

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        ctl_struc *  ctl        I    input controls
        in_rec_struc * in_rec    I/O  controls for input record reading
        int       stray_stop    I    if 1, remove storage used for stray 
                                     light processing - do at program end

    Modification history:

    W. Robinson, SAIC  3 June 2011  Original development

----------------------------------------------------------------------------*/ {
    int ibnd, ilin, ipix, idet, lpsf, ppsf, pdat;
    static int nbnd_s, npix, nlin, npix_s, nlin_s, lin_margin;
    static vir_straylt_struc stlt;
    static int stray_init = 0;
    static float *bnd_store; /* storage for the unchanged radiance */
    float *psf_bd; /* a pointer to convol array to use for a band, detector */
    float sum_psf, sum_rad, psf_val;

    /*
     *  close down if stop switch is set
     */
    if (stray_stop == 1) {
        free(stlt.psf);
        free(bnd_store);
        printf("%s: Shutting down the stray light artifact addition\n",
                __FILE__);
    } else {
        /*
         *  do initialization
         */
        if (stray_init == 0) {
            stray_init = 1;

            printf("%s, %d: Entering routine viirs_straylt\n", __FILE__, __LINE__);
            printf("stray_tbl = %s\n", ctl->stray_tbl);
            /*
             *  read in the psf = stray light convol array, and make a
             *  holding array for 1 band of a scan
             */
            if (viirs_straylt_rd(ctl->stray_tbl, &stlt) != 0)
                return 1;

            nbnd_s = stlt.nbands; /* assume bands in psf go M1 - Mn bands */
            nlin = in_rec->ndet_scan;
            npix = in_rec->npix;
            lin_margin = in_rec->margin[0];
            npix_s = stlt.nsamp;
            nlin_s = stlt.ndet;

            if ((bnd_store = (float *) malloc(nlin * npix * sizeof ( float)))
                    == NULL) {
                printf(
                        "%s, %d: Unable to allocate band storage for stray light artifact\n",
                        __FILE__, __LINE__);
                return 1;
            }
        }

        /*
         *  go through the bands, lines, and pixels of scan to apply stray light to
         */
        for (ibnd = 0; ibnd < nbnd_s; ibnd++) {
            /*
             *  Keep a pristine copy of the original band's data
             */
            memcpy(bnd_store, in_rec->bnd_lt[ibnd], nlin * npix * sizeof (float));

            for (idet = 0; idet < nlin_s; idet++) {
                ilin = idet + lin_margin; /* line in in_rec */
                /*
                 *  the convol array can now be selected from the table for 
                 *  this band and detector
                 */
                psf_bd = stlt.psf + (npix_s * nlin_s) * (ibnd + nbnd_s * idet);

                for (ipix = 0; ipix < npix; ipix++) {
                    /*
                     *  we are at the point of collecting the correction for a sample
                     *  apply it only if that value is good
                     */
                    if (*(in_rec->bnd_q[ibnd] + ipix + npix * ilin) == 0) {
                        /*
                         *  This is where the convolution begins, but only using 
                         *  good radiance values
                         */
                        sum_psf = 0.;
                        sum_rad = 0.;

                        for (lpsf = 0; lpsf < nlin_s; lpsf++) {
                            for (ppsf = 0; ppsf < npix_s; ppsf++) {
                                /*
                                 *  compute te pixel of the radiance to apply convol to,
                                 *  account for the center sample of the psf array
                                 */
                                pdat = ipix + ppsf - stlt.csamp;
                                if ((pdat >= 0) && (pdat < npix)) {
                                    /*
                                     *  if the quality at this point is good, add the info
                                     */
                                    if (*(in_rec->bnd_q[ibnd] + ipix + ilin * npix) == 0) {
                                        psf_val = *(psf_bd + ppsf + npix_s * lpsf);
                                        sum_psf += psf_val;
                                        sum_rad += psf_val * *(bnd_store + pdat +
                                                (lpsf + lin_margin) * npix);
                                    }
                                }
                            }
                        }
                        /*
                         *  make the artifact-modified radiance value
                         */
                        *(in_rec->bnd_lt[ibnd] + ipix + ilin * npix) = sum_rad / sum_psf;
                    }
                }
            }
        }
    }
    /*
     *  and end
     */
    return 0;
}

int viirs_straylt_rd(char *file, vir_straylt_struc *stlt)
/*-----------------------------------------------------------------------------
    Routine:  viirs_straylt_rd

    Description:  read the stray light artifact coefficients
      
    Returns type: int - 0 if good, else I/O error

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        char *    file          I    stray light table psf file name
        vir_straylt_struc *  gain   O    structure with stray light 
                                     coefficients and description of it
                                     more descrip in vir_straylt_struc
                                     descrip in viirs_sim_sdr.h

    Modification history:

    W. Robinson, SAIC  3 June 2011  Original development

----------------------------------------------------------------------------*/ {
    int nrec_det = 16, ndet = 16, nbands, nsamp;

    h5io_str fid;

    /*
     *  open the stray light table file and get the size information
     */
    if (strcmp(file, "Unspecified") != 0) {
        if (h5io_openr(file, 0, &fid) != 0) {
            printf("%s, %d, Unable to open stray light table file: \n%s\n",
                    __FILE__, __LINE__, file);
            return 1;
        }
        if (h5io_rd_attr(&fid, "nbands", (void *) &(stlt->nbands)) != 0) {
            printf(
                    "%s, %d, Unable to read nbands attr from stray light file: \n%s\n",
                    __FILE__, __LINE__, file);
            return 1;
        }
        if (h5io_rd_attr(&fid, "nrec_det", (void *) &(stlt->nrec_det)) != 0) {
            printf(
                    "%s, %d, Unable to read nrec_det attr from stray light file: \n%s\n",
                    __FILE__, __LINE__, file);
            return 1;
        }
        if (h5io_rd_attr(&fid, "nsamp", (void *) &(stlt->nsamp)) != 0) {
            printf(
                    "%s, %d, Unable to read nsamp attr from stray light file: \n%s\n",
                    __FILE__, __LINE__, file);
            return 1;
        }
        if (h5io_rd_attr(&fid, "ndet", (void *) &(stlt->ndet)) != 0) {
            printf(
                    "%s, %d, Unable to read ndet attr from stray light file: \n%s\n",
                    __FILE__, __LINE__, file);
            return 1;
        }
        if (h5io_rd_attr(&fid, "csamp", (void *) &(stlt->csamp)) != 0) {
            printf(
                    "%s, %d, Unable to read csamp attr from stray light file: \n%s\n",
                    __FILE__, __LINE__, file);
            return 1;
        }

        nbands = stlt->nbands;
        nsamp = stlt->nsamp;

        if ((stlt->nrec_det != nrec_det) || (stlt->ndet != ndet)) {
            printf(
                    "%s, %d, Unexpected stray light array dimensions found in file: \n%s\n",
                    __FILE__, __LINE__, file);
            printf("   nrec_det: %d, expected: %d\n", stlt->nrec_det, nrec_det);
            printf("   ndet: %d, expected: %d\n", stlt->ndet, ndet);
            return 1;
        }
    } else {
        printf("%s, %d: Unexpectedly found Unspecified file stray light table\n",
                __FILE__, __LINE__);
        return 1;
    }
    /*
     *  allocate space for the psf array and read it in
     */
    if ((stlt->psf = (float *)
            malloc(nbands * nrec_det * nsamp * ndet * sizeof ( float))) == NULL) {
        printf("%s, %d: Error allocating the stlt->psf\n", __FILE__, __LINE__);
        return 1;
    }

    if (h5io_grab_ds(&fid, "psf", (void *) stlt->psf) != 0) {
        printf(
                "%s, %d, Unable to read psf dataset from stray light table file: \n%s\n",
                __FILE__, __LINE__, file);
        return 1;
    }

    if (h5io_close(&fid) != 0) {
        printf("%s, %d, Error closing the stray light table:\n%s\n", __FILE__,
                __LINE__, file);
        return 1;
    }
    /*
     *  return the psf and its description
     */
    return 0;
}
