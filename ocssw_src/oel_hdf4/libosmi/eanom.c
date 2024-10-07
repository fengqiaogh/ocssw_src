/*
 *----------------------------------------------------------------------
 * @(#) eanom.c		1.0	02 Apr 98	<shc>
 * Copyright (c) 1993, CSIRO Division of Oceanography
 * Copyright (c) 1998, Datron Transco Inc.
 *----------------------------------------------------------------------
 *
 * eanom --
 *
 *	Mean anomaly to eccentric anomaly iteration.
 *
 * Results:
 *
 *	Calculate eccentric anomaly from mean anomaly using modified Newton
 *	iteration.  No angle normalisation is done here in order to retain
 *	the notion of the mean and eccentric anomalies increasing monotonically
 *	with time.  Inputs are mean anomaly (MANOM) and orbit eccentricity
 *	(ECC), result is the corresponding eccentric anomaly to absolute
 *	precision EPS.
 *
 *	On error, ERRSTR is set to point to an error message and  the
 *	value FP_ERRVAL is returned, since many systems do not have NaN
 *	definitions.
 *
 * Side effects:
 *	None.
 *
 * History:
 *  30 Mar 93  <shc>
 *	Converted from FORTRAN to C.
 *
 *  12 Dec 94  <shc>
 *	Added errmsg stuff for integration with Python
 *
 *  02 Apr 98  <shc>
 *	Added dual EPS limits and increased MAXIT from 10 to 100 to
 *	reduce incidence of convergence failures.
 *
 *  13 Apr 98  <shc>
 *	Changed to use same convergence scheme as ktosgp4(): iterate
 *	to convergence, then do another couple of steps for good measure.
 *
 *----------------------------------------------------------------------
 */

#include <stdio.h>
#include "orbit.h"

#define MAXIT 100  /* Maximum numer of iterations */
#define EPS 1.0e-12  /* Error goal */

double
eanom(manom, ecc, errmsg)
double manom, ecc;
char **errmsg;
{
    int i, converged;
    double delta, ea, func, ean;

    if (ecc < 0.0 || ecc > 1.0) {
        if (errmsg != NULL)
            *errmsg = "eanom - non-elliptic orbit";
        return FP_ERRVAL;
    }

    ea = manom;
    for (i = converged = 0; i < MAXIT && converged < 2; i++) {
        func = ea - ecc * sin(ea) - manom;
        ean = ea - func / (1.0 - ecc * cos(ea - func / 2.0));
        delta = ean - ea;
        if (fabs(delta) < EPS) {
            converged += 1;
        }
        ea = ean;
    }

    if (converged) {
        return (ea);
    } else {
        if (errmsg != NULL)
            *errmsg = "eanom: eccentric anomaly iteration "
                "failed to converge";
        return FP_ERRVAL;
    }
}
