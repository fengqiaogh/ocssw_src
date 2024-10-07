/*
 *  fuzzy_member_func.c
 *  
 *
 *  Created by Tim Moore on 9/1/09.
 *  Copyright 2009 UNH. All rights reserved.
 * 
 *  Input:
 *		rrs - a spectral rrs(0-) at specified wavelengths from a satellite
 *      urrs - the matrix of mean rrs(0-) spectra for the different optical water types
 *		y3inv - the inverted covariance matrix for the optical water types
 *		nclasses - the number of optical water types
 *		df - the number of rrs wavelengths
 *
 *		
 *		outdata - pointer to the output fuzzy memberships to the optical water types for the rrs input
 *
 */

#define NRANSI
#include "nr.h"
#include "nrutil.h"
#include <math.h>
#define DEBUG_FUZZY 0
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*

void fuzzy_func_v3(float *rrs, float **urrs, float ***y3inv, int nclasses, int nowts, int df, double *outdata)
 {

    int nmax, i, j, k;
    float *sa, *ax, *b, *y, *alphachi, **yinv, gamm, x, gln, chiz, z2;
    unsigned long msize, *ija;

    nmax = 2 * df * df + 1;

    ija = lvector(1, nmax);
    ax = vector(1, df);
    b = vector(1, df);
    sa = vector(1, nmax);

    yinv = matrix(1, df, 1, df);
    y = vector(1, df);
    alphachi = vector(1, nclasses);


    if (DEBUG_FUZZY) {
        printf("\ninside fuzzy_func...");
        printf("\nnclasses = %d", nclasses);
        printf("\ndf = %d", df);
    }

    for (i = 1; i < nclasses + 1; i++) {
        for (j = 1; j < df + 1; j++) {
            y[j] = rrs[j - 1] - urrs[j][i];
            if (DEBUG_FUZZY) printf("\n%d %d %f %f %f", i, j, y[j], rrs[j - 1], urrs[j][i]);
            for (k = 1; k < df + 1; k++) {
                yinv[j][k] = y3inv[i][j][k];
            }

        }

        sprsin(yinv, df, 0.5, nmax, sa, ija);
        msize = ija[1] - 2;
        sprsax(sa, ija, y, b, msize);
        z2 = 0;
        for (j = 1; j < df + 1; j++) {
            z2 = z2 + y[j] * b[j];
        }
        x = z2 / 2.0;
        chiz = ((float) df) / 2.0;
        if (DEBUG_FUZZY) printf("\nz2=%f", z2);

        /*      
        if (x < 0) {   // baf
          free_vector(y,1,df+1);
          free_vector(alphachi,1,nclasses);
          free_matrix(yinv,1,df,1,df);

          free_vector(sa,1,nmax);
          free_vector(b,1,df);
          free_vector(ax,1,df);
          free_lvector(ija,1,nmax);
          return(1);
        }
         */

        if (x <= (chiz + 1.0)) {
            gser(&gamm, chiz, x, &gln);
            alphachi[i] = 1.0 - gamm;
        } else {
            gcf(&gamm, chiz, x, &gln);
            alphachi[i] = gamm;
        }

        if (DEBUG_FUZZY) printf("\ngamm=%f", gamm);

    }

    if (DEBUG_FUZZY) {
        for (i = 0; i < nclasses; i++) {
            printf("\nalpha[%d]=%f", i, alphachi[i + 1]);
        }
    }

    for (i = 0; i < nowts - 1; i++) {
        *outdata++ = (double) alphachi[i + 1];
    }

    *outdata++ = alphachi[9] + alphachi[10] + alphachi[11] + alphachi[12] + alphachi[13] + alphachi[14] + alphachi[15] + alphachi[16];

    free_vector(y, 1, df + 1);
    free_vector(alphachi, 1, nclasses);
    free_matrix(yinv, 1, df, 1, df);

    free_vector(sa, 1, nmax);
    free_vector(b, 1, df);
    free_vector(ax, 1, df);
    free_lvector(ija, 1, nmax);

    if (DEBUG_FUZZY) printf("\nleaving fuzzy_member_func...");


}
