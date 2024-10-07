/*
 *  covariance_inversion.c
 *  
 *
 *  Created by Tim Moore on 12/16/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "nr.h"
#include "nrutil.h"
#include <math.h>
#define DEBUG_COVAR 0
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

void covariance_inversion(float *rrs_cov, int nclasses, int df, float ***y3inv) {

    int i, j, k, ii, jj, *indx;
    float d, **a, *col;

    void ludcmp();

    indx = ivector(1, df);
    col = vector(1, df);
    a = matrix(1, df, 1, df);

    if (DEBUG_COVAR) {
        printf("\nptr addr =%p x0=%f", rrs_cov, *(rrs_cov));
        printf("\nnclasses = %d", nclasses);
        printf("\ndf = %d", df);
    }

    /******* get inversion of covariance matrix ***********/

    for (i = 0; i < nclasses; i++) {
        for (j = 0; j < df; j++) {
            for (k = 0; k < df; k++) a[j + 1][k + 1] = *rrs_cov++;
        }

        ludcmp(a, df, indx, &d);
        for (jj = 1; jj <= df; jj++) {
            for (ii = 1; ii <= df; ii++) col[ii] = 0.0;
            col[jj] = 1.0;
            lubksb(a, df, indx, col);
            for (ii = 1; ii <= df; ii++) y3inv[i + 1][ii][jj] = col[ii];
        }

    }

    free_ivector(indx, 1, df);
    free_vector(col, 1, df);
    free_matrix(a, 1, df, 1, df);

}
