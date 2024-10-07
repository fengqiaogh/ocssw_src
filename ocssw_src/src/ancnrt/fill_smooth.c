/*  hope we don't need this
#include "ancil.h"
#include "no2_avg.h"
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "ancil.h"
#include "o3_toms.h"

int fill_smooth(int16_t *arr, char *arrqc, int32_t npix, int32_t nlin)
/*******************************************************************

   fill_smooth

   purpose:  make the ozone fill seamlessly join to the day's data
      Don't extend data beyond what was filed already

   Returns type: 0 if all went well, else 1

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int16 *           arr             I/O     input array
      char *            arrqc           I       QC values assigned in 
                                                o3_toms: 0 good data, 1 - 19
                                                fill data, 20 bad = data = 0
      int32             npix            I       # pixels in the input array
      int32             nlin            I       # lines in the input array

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       12 Dec 2013     Original development

 *******************************************************************/
 {
    int32_t ntot, iv;
    float *lcl_data, *wt;
    char *lcl_qc;
    /*
     *  we will apply field_extend to the good data using a distance of 5
     *  pixels and then weight it in to the fill made before
     */
    /*
     *  allocate space to work in and transfer the good array values to a float 
     *  array and make a qc array for use in field_extend
     */
    ntot = npix * nlin;
    if (((lcl_data = (float *) malloc(ntot * sizeof ( float))) == NULL)
            || ((lcl_qc = (char *) malloc(ntot * sizeof ( char))) == NULL)
            || ((wt = (float *) malloc(npix * nlin * sizeof ( float))) == NULL)) {
        printf("%s, %d E: unable to allocate local data storage\n",
                __FILE__, __LINE__);
        return 1;
    }

    for (iv = 0; iv < ntot; iv++) {
        *(lcl_data + iv) = 0.;
        *(lcl_qc + iv) = 10;
        if (*(arrqc + iv) == 0) {
            *(lcl_data + iv) = (float) *(arr + iv);
            *(lcl_qc + iv) = 0;
        }
    }
    /*
     *  use field_extend to extend the local data and qc
     */
    if (field_extend(lcl_data, lcl_qc, npix, nlin, 7, wt) != 0)
        return 1;
    /*
     *  Now, go back and update the grid array
     */
    for (iv = 0; iv < ntot; iv++) {
        /*  if lcl_qc is 0 or unchanged, leave the input data alone */
        /*  if lcl_qc is 10, no weighting gets applied here either and no change */
        /*  in the field extend region that has fill, we combine the extended 
            field weighted by the weight with the fill field weighted by 
            1 - weight */
        if ((*(lcl_qc + iv) == 1) && (*(arrqc + iv) != 20)) {
            *(arr + iv) = *(wt + iv) * *(lcl_data + iv) +
                    (1. - *(wt + iv)) * (float) *(arr + iv);
        }
    }
    /*
     *  free the space allocated for work and return
     */
    free(lcl_data);
    free(lcl_qc);
    free(wt);

    return 0;
}

int field_extend(float *grid, char *flags, int np_grid, int nl_grid,
        int max_extend, float* weights)
/*******************************************************************

   field_extend

   purpose: extens a field of data beyond its limits

   Returns type: int - 0 if good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      float *           grid            I/O     earth covering grid containing
                                                the data field which will be
                                                extended
      char *            flags           I/O     state of grid points in grid:
                                                0 - data there, 1 - fill
                                                (extended) data
                                                10 - no data
                                                On input, expect only 0 or 1
      int               np_grid          I      # pixels in grid
      int               nl_grid          I      # lines in grid
      int               max_extend       I      maximum extent of the field
      float *           weights          O      array of weights - linearly
                                                decreasing from 1 to 0 at
                                                max_extend

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       17-Mar-2008     Original development
      W. Robinson, SAIC 21 Jun 2013     re-think algorithm some and adopt
                                      same flag meaning as elsewhere in no2_avg

 *******************************************************************/ {
    float wt_shell, *ker, wt, sum, dist, r1, r2;
    int ishell, ker_np, ker_nl, ker_cp, ker_cl, ln_grid, px_grid,
            ker_lin, ker_pix, tl_grid, tp_grid, n_shell;
    /*
     *  fill the adjacent pixels to those already filled out to the
     *  specified distance of max_extend
     */
    printf("%s: extending field with distance of %d\n", __FILE__, max_extend);
    ishell = 0;
    do {
        /*  set the shell's weight in units of 1 / ( max_extend + 1 )  */
        wt_shell = (float) (max_extend - ishell) / (float) (max_extend + 1);
        /*
         *  new distance scheme with one distance of radius at start and another
         *  at the end just linearly expanding (r2 should be at least larger
         *  than max_extend or else it could have no data to use
         */
        r1 = 5.;
        r2 = max_extend * 2.0;
        dist = r1 + (float) ishell * (r2 - r1) / (float) max_extend;
        /*
         printf( "%s, %d, Shell %d of %d, wt: %f, dist: %f\n", __FILE__, __LINE__,
           ishell, max_extend, wt_shell, dist );
         */
        /*
         *  create an averaging kernel with radius dist and value 1
         */
        if (mk_ker(dist, &ker, &ker_np, &ker_nl, &ker_cp, &ker_cl) != 0)
            return 1;
        n_shell = 0;
        /*
        printf( "wt_shell = %f, dist = %f, ker_np = %d, ker_nl = %d, ker_cp = %d, ker_cl = %d\n", wt_shell, dist, ker_np, ker_nl, ker_cp, ker_cl );
         */
        /*
         *  get the kernel weighted average of the original field and assign
         *  the weight
         */
        for (ln_grid = 0; ln_grid < nl_grid; ln_grid++) {
            for (px_grid = 0; px_grid < np_grid; px_grid++) {
                /*
                 *  process at each grid point
                 */
                if (*(flags + ln_grid * np_grid + px_grid) != 10) {
                    /*
                     *  for filled already, just set the weights to 1
                     */
                    if (ishell == 0)
                        *(weights + ln_grid * np_grid + px_grid) = 1.;
                } else /* for a point to extend into */ {
                    /*
                     *  make sure only pixels in the next shell out are considered
                     *  (ie, an adjacent pixel has data)
                     */
                    wt = 0.;
                    for (ker_lin = -1; ker_lin < 2; ker_lin++) {
                        tl_grid = ln_grid + ker_lin;
                        if ((tl_grid >= 0) && (tl_grid < nl_grid)) {
                            for (ker_pix = -1; ker_pix < 2; ker_pix++) {
                                tp_grid = px_grid + ker_pix;
                                if (tp_grid < 0)
                                    tp_grid = tp_grid + np_grid;
                                else if (tp_grid >= np_grid)
                                    tp_grid = tp_grid - np_grid;
                                if (*(flags + tl_grid * np_grid + tp_grid) <= 1)
                                    wt = 1.;
                            }
                        }
                    }
                    /*
                     *  The following pixels will be in the next shell
                     */
                    if (wt == 1.) {
                        /*
                         *  for unfilled, accumulate the weight and sum for the grid point
                         */
                        wt = 0;
                        sum = 0.;
                        for (ker_lin = 0; ker_lin < ker_nl; ker_lin++) {
                            /*
                             *  get the grid line at the start kernel line and
                             *  only compute if the line is in the grid
                             */
                            tl_grid = ln_grid + ker_lin - ker_cl;
                            if ((tl_grid >= 0) && (tl_grid < nl_grid)) {
                                for (ker_pix = 0; ker_pix < ker_np; ker_pix++) {
                                    /*  only do where the kernel is > 0  */
                                    if (*(ker + ker_lin * ker_np + ker_pix) > 0.) {
                                        /*
                                         *  wrap the pixel around
                                         */
                                        tp_grid = px_grid + ker_pix - ker_cp;
                                        if (tp_grid < 0)
                                            tp_grid = tp_grid + np_grid;
                                        else if (tp_grid >= np_grid)
                                            tp_grid = tp_grid - np_grid;
                                        /*     */
                                        if (*(flags + tl_grid * np_grid + tp_grid) == 0) {
                                            wt += *(ker + ker_lin * ker_np + ker_pix);
                                            sum += *(grid + tl_grid * np_grid + tp_grid);
                                        }
                                    }
                                } /*  end of kernel pixel loop  */
                            }
                        } /* end of kernel line loop  */
                    }
                    /*
                     *  get the extended value
                     */
                    if (wt > 0.) {
                        *(grid + ln_grid * np_grid + px_grid) = sum / wt;
                        *(flags + ln_grid * np_grid + px_grid) = 2;
                        *(weights + ln_grid * np_grid + px_grid) = wt_shell;
                        n_shell++;
                    }
                }
            } /*  end grid pixels */
        } /* end grid lines */
        free(ker);
        /*
         *  Note that I set the flags for the new shell to 2 to avoid detecting
         *  the shell as a previous shell.  Now the shell is done, set flag 2 -> 1
         */
        for (px_grid = 0; px_grid < (np_grid * nl_grid); px_grid++)
            if (*(flags + px_grid) == 2)
                *(flags + px_grid) = 1;
        ishell++;
        /* printf( "n_shell = %d\n", n_shell ); */
    } while ((ishell < max_extend) && (n_shell > 0)); /* end shell loop */
    /*
     *  return the extended field and the weights
     */
    return 0;
}

int mk_ker(float dist, float **ker, int *ker_np, int *ker_nl,
        int *ker_cp, int *ker_cl)
/*******************************************************************
   mk_ker

   purpose: create a weighting kernel of a flat value of 1.

   Returns type: int - 0 if good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      float             dist             I      maximum distance of influence
                                                of the averaging kernel
      float **          ker             I/O     kernel with weights
                                                (pass in pointer to space)
      int *             ker_np           O      # pixels in kernel
      int *             ker_nl           O      # lines in kernel
      int *             ker_cp           O      center pixel of kernel
      int *             ker_cl           O      center line of kernel

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       17-Mar-2008     Original development

 *******************************************************************/ {
    int il, ip, n_out, iret = 0;
    float rad;
    /*
     */
    n_out = (int) (dist + 1.);
    *ker_np = 2 * n_out + 1;
    *ker_nl = *ker_np;
    *ker_cp = n_out;
    *ker_cl = *ker_cp;
    /*
     *  create the kernel array
     *  it will just be a flat weight of 1 and with a circular edge
     */
    if ((*ker = (float *) malloc(*ker_np * *ker_nl * sizeof (float))) == NULL) {
        printf("%s:  unable to allocate storage for ker\n", __FILE__);
        iret = 1;
    } else {
        for (il = 0; il < *ker_nl; il++) {
            for (ip = 0; ip < *ker_np; ip++) {
                rad = sqrt(pow((float) (*ker_cl - il), 2) +
                        pow((float) (*ker_cp - ip), 2));
                *(*ker + il * *ker_np + ip) = (rad < dist) ? 1. : 0.;
            }
        }
    }
    return iret;
}
