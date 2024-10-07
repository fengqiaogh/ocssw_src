/* ============================================================== */
/* Module:  swl0_utils.c                                          */
/* Purpose: Utilities for manipulating seawifs L0 data.           */
/* Author:  B.A. Franz, General Scences Corp., 9/97               */
/* ============================================================== */

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <math.h>
#include "swl0_proto.h"

#define TIMELEN 25

/* -------------------------------------------------------------- */
/* scid2mnftype() - decoms the spacecraft ID to minor frame type  */

/* -------------------------------------------------------------- */
INT16 scid2mnftype(INT16 scid[]) {
    return (scid[1] & 0x000F);
}


/* -------------------------------------------------------------- */
/* scid2mnftype() - decoms the spacecraft ID to minor frame numb  */

/* -------------------------------------------------------------- */
INT16 scid2mnfnum(INT16 scid[]) {
    return ( (scid[0] & 0x0180) >> 7);
}


/* -------------------------------------------------------------- */
/* ttag2unix() - converts seawifs timetags to secs since 1/1/70   */

/* -------------------------------------------------------------- */
FLOAT64 ttag2unix(INT16 ttag[]) {
    FLOAT64 sec = 7.2688320e8; /* secs from 1/1/70 to 1/13/93  */
    INT32 day = 0;
    INT32 msec = 0;

    day = ttag[0] << 3 | ttag[1] >> 7;
    msec = (ttag[1] & 0x7F) << 20 | (ttag[2] << 10) | ttag[3];
    sec += ((double) day) * 86400.0 + ((double) msec) / 1000.0;

    return sec;
}


/* -------------------------------------------------------------- */
/* unix2timeStr() - converts secs since 1/1/70 to timedate string */

/* -------------------------------------------------------------- */
char *unix2timeStr(FLOAT64 usec) {
    struct tm *trec;
    time_t utime = (time_t) usec;
    static char timeStr[TIMELEN];

    trec = gmtime(&utime);
    strftime(timeStr, TIMELEN, "%b %d %Y %H:%M:%S:000", trec);

    return (timeStr);
}


/* -------------------------------------------------------------- */
/* ttag2ydmsec() - converts seawifs timetags to year, day, msec   */

/* -------------------------------------------------------------- */
void ttag2ydmsec(INT16 ttag[], INT16 *year, INT16 *day, INT32 *msec) {
    INT32 jday; /* julian day number            */
    INT32 rday; /* reference day number         */

    *msec = (ttag[1] & 0x7F) << 20 | (ttag[2] << 10) | ttag[3];
    jday = (ttag[0] << 3 | ttag[1] >> 7) + 2449001;

    /* The following was taken from jdate.f of F. S. Patt         */
    rday = jday - 2415020; /* Days since January 0, 1900 */
    *year = 4 * rday / 1461; /* Years since 1900           */

    *day = rday - 1461 * (*year - 1) / 4 - 365;
    *year += 1900;

    return;
}


/* -------------------------------------------------------------- */
/* timeError() - returns > 0 if any time error flag is set.       */

/* -------------------------------------------------------------- */
BYTE timeError(swl0indx *indx, INT32 irec) {
    return ( indx->rec[irec].tRanError ||
            indx->rec[irec].tSeqError ||
            indx->rec[irec].tDifError);
}


/* -------------------------------------------------------------- */
/* timeSeqError() - returns > 0 if timetags not increasing. Only  */
/*                  identifies single-frame sequence errors.      */

/* -------------------------------------------------------------- */
BYTE timeSeqError(swl0indx *indx, INT32 irec) {
    INT32 i1 = MAX(irec - 1, 0);
    INT32 i2 = irec;
    INT32 i3 = MIN(irec + 1, indx->nrecs - 1);

    if (irec < 0 || irec >= indx->nrecs)
        return (0);

    /* Use times that were previously unflagged, if possible. */
    while (i1 > 0 && timeError(indx, i1)) i1--;
    while (i3 < indx->nrecs - 1 && timeError(indx, i3)) i3++;

    /* Don't attempt test if frames are already flagged */
    if (timeError(indx, i1) || timeError(indx, i2) || timeError(indx, i3))
        return (0);

    /* Frame type must be consistent for this test. */
    if (indx->type != HRPT &&
            (indx->rec[i2].mnftype != indx->rec[i1].mnftype ||
            indx->rec[i2].mnftype != indx->rec[i3].mnftype)) {
        return (0);
    }

    /* Check for large deviations at start and end frames */
    if ((irec == 0) &&
            (indx->rec[i3].time - indx->rec[i2].time > DELSCENEGAC))
        return (1);
    if ((irec == indx->nrecs - 1) &&
            (indx->rec[i2].time - indx->rec[i1].time > DELSCENEGAC))
        return (1);
    if ((irec == 0) &&
            (indx->rec[i3].time - indx->rec[i2].time <= 0))
        return (1);
    if ((irec == indx->nrecs - 1) &&
            (indx->rec[i2].time - indx->rec[i1].time <= 0))
        return (1);

    /* Bounding times must be sequential for test to be valid. We don't */
    /* to flag the wrong frame.                                         */
    if (indx->rec[i3].time < indx->rec[i1].time)
        return (0);

    /* General case */
    if ((indx->rec[i2].time < indx->rec[i1].time) ||
            (indx->rec[i2].time > indx->rec[i3].time) ||
            ((i1 != i2) && (indx->rec[i2].time == indx->rec[i1].time)))

        return (1);

    return (0);
}


/* -------------------------------------------------------------- */
/* timeContiguous() - returns 1 if timetags are contiguous        */

/* -------------------------------------------------------------- */
BYTE timeContiguous(swl0indx *indx, INT32 irec) {
    INT32 i1;
    INT32 i2;
    INT32 i3;
    FLOAT64 DelTime;
    FLOAT64 tdiff1;
    FLOAT64 tdiff2;

    if (irec < 0 || irec >= indx->nrecs)
        return (0);

    /* Determine appropriate center and bounding indices */
    if (irec == 0) {
        i1 = 0;
        i2 = MIN(1, indx->nrecs - 1);
        i3 = MIN(2, indx->nrecs - 1);
    } else if (irec == indx->nrecs - 1) {
        i3 = indx->nrecs - 1;
        i2 = MAX(i3 - 1, 0);
        i1 = MAX(i3 - 2, 0);
    } else {
        i1 = MAX(irec - 1, 0);
        i2 = irec;
        i3 = MIN(irec + 1, indx->nrecs - 1);
    }

    /* Frame must be (thus far) free of timing errors */
    if (timeError(indx, i1) || timeError(indx, i2) || timeError(indx, i3)) {
        printf("Previous timing error(s) at frame %d\n", irec);
        return (0);
    }

    /* Frame SCID field must be valid */
    if (indx->rec[i1].scidError ||
            indx->rec[i2].scidError ||
            indx->rec[i3].scidError) {
        printf("Previous SCID error(s) at frame %d\n", irec);
        return (0);
    }

    /* Frame type must be consistent as well. */
    if (indx->rec[i2].mnftype != indx->rec[i1].mnftype ||
            indx->rec[i2].mnftype != indx->rec[i3].mnftype) {
        printf("Inconsistent frame types at frame %d\n", irec);
        return (0);
    }

    /* Determine expected time difference between frames */
    if (indx->rec[i2].mnftype == GACTYPE)
        DelTime = DTGAC;
    else
        DelTime = DTLAC;

    /* Compute actual differences */
    tdiff1 = indx->rec[i2].time - indx->rec[i1].time;
    tdiff2 = indx->rec[i3].time - indx->rec[i2].time;

    /* Check center time relative to bounding frames */
    if (tdiff1 > DelTime - CLOCKVAR &&
            tdiff1 < DelTime + CLOCKVAR &&
            tdiff2 > DelTime - CLOCKVAR &&
            tdiff2 < DelTime + CLOCKVAR)
        return (1);
    else {
        printf("%d %lf %lf\n", irec, tdiff1, tdiff2);
        return (0);
    }
}


/* -------------------------------------------------------------- */
/* timeConsistent() - returns 0 if timetag differences are not    */
/*                    multiples of expected frame rate            */
/* Note that this will not detect shifts, as it looks for timetag */
/* which is inconsistent with the preceeding and following times  */

/* -------------------------------------------------------------- */
BYTE timeConsistent(swl0indx *indx, INT32 irec) {
    INT32 i1;
    INT32 i2;
    INT32 i3;
    FLOAT64 DelTime;
    FLOAT64 DelScene;
    FLOAT64 tdiff1;
    FLOAT64 tdiff2;
    FLOAT64 terror1;
    FLOAT64 terror2;

    if (irec < 0 || irec >= indx->nrecs)
        return (1);

    /* Determine appropriate center and bounding indices */
    if (irec == 0) {
        i1 = 0;
        i2 = MIN(1, indx->nrecs - 1);
        i3 = MIN(2, indx->nrecs - 1);
    } else if (irec == indx->nrecs - 1) {
        i3 = indx->nrecs - 1;
        i2 = MAX(i3 - 1, 0);
        i1 = MAX(i3 - 2, 0);
    } else {
        i1 = MAX(irec - 1, 0);
        i2 = irec;
        i3 = MIN(irec + 1, indx->nrecs - 1);
    }

    /* Avoid time difference test across frame type changes. */
    if (indx->rec[i2].mnftype == GACTYPE &&
            (indx->rec[i3].mnftype != indx->rec[i1].mnftype)) {
        return (1);
    }

    /* Determine expected time difference between frames */
    if (indx->rec[i2].mnftype == GACTYPE) {
        DelTime = DTGAC;
        DelScene = DELSCENEGAC;
    } else {
        DelTime = DTLAC;
        DelScene = DELSCENELAC;
    }

    /* Compute actual differences */
    tdiff1 = indx->rec[i2].time - indx->rec[i1].time;
    tdiff2 = indx->rec[i3].time - indx->rec[i2].time;

    /* Avoid consistency check across scene breaks, GAC/LAC files only */
    if (indx->type == GAC && (tdiff1 > DelScene || tdiff2 > DelScene))
        return (1);

    /* Compute residual error */
    terror1 = fmod(tdiff1, DelTime);
    terror2 = fmod(tdiff2, DelTime);
    terror1 = ABS(terror1);
    terror2 = ABS(terror2);
    if (terror1 > DelTime / 2) terror1 = ABS(terror1 - DelTime);
    if (terror2 > DelTime / 2) terror2 = ABS(terror2 - DelTime);

    /* Check error relative to allowed clock variance. We require time */
    /* error on both sides of the center frame, as we want to avoid    */
    /* flagging both frames across a time shift.                       */

    if (tdiff1 == 0.0 || (terror1 > CLOCKVAR && terror2 > CLOCKVAR)) {
        printf("Time Difference Error at Frame %d: %lf %lf %lf\n",
                irec, tdiff1, terror1, terror2);
        return (0);
    } else {
        return (1);
    }
}


/* -------------------------------------------------------------- */
/* timeShifted() - returns 0 if timetag differences are not       */
/*                    multiples of expected frame rate            */

/* -------------------------------------------------------------- */
BYTE timeShifted(swl0indx *indx, INT32 irec, FLOAT64 *shiftval) {
    INT32 i1;
    INT32 i2;
    INT32 i3;
    FLOAT64 DelTime;
    FLOAT64 DelScene;
    FLOAT64 tdiff1;
    FLOAT64 tdiff2;
    FLOAT64 terror1;
    FLOAT64 terror2;

    *shiftval = 0.0;

    if (irec < 0 || irec >= indx->nrecs)
        return (1);

    /* Determine appropriate center and bounding indices */
    if (irec == 0) {
        i1 = 0;
        i2 = MIN(1, indx->nrecs - 1);
        i3 = MIN(2, indx->nrecs - 1);
    } else if (irec == indx->nrecs - 1) {
        i3 = indx->nrecs - 1;
        i2 = MAX(i3 - 1, 0);
        i1 = MAX(i3 - 2, 0);
    } else {
        i1 = MAX(irec - 1, 0);
        i2 = irec;
        i3 = MIN(irec + 1, indx->nrecs - 1);
    }

    /* Avoid time difference test across frame type changes. */
    if (indx->rec[i2].mnftype == GACTYPE &&
            (indx->rec[i3].mnftype != indx->rec[i1].mnftype)) {
        return (0);
    }

    /* Determine expected time difference between frames */
    if (indx->rec[i2].mnftype == GACTYPE) {
        DelTime = DTGAC;
        DelScene = DELSCENEGAC;
    } else {
        DelTime = DTLAC;
        DelScene = DELSCENELAC;
    }

    /* Compute actual differences */
    tdiff1 = indx->rec[i2].time - indx->rec[i1].time;
    tdiff2 = indx->rec[i3].time - indx->rec[i2].time;

    /* Avoid check across scene breaks, GAC/LAC files only */
    if (indx->type == GAC && (tdiff1 > DelScene || tdiff2 > DelScene))
        return (0);

    /* Compute residual error */
    terror1 = fmod(tdiff1, DelTime);
    if (tdiff1 < DelTime)
        terror1 -= DelTime;

    terror2 = fmod(tdiff2, DelTime);
    if (tdiff2 < DelTime)
        terror2 -= DelTime;


    /* Check error relative to allowed clock variance. We require time */
    /* shift preceeding center frame, and no timeshift following.      */

    if (fabs(terror1) > CLOCKVAR && fabs(terror2) <= CLOCKVAR) {
        *shiftval = terror1;
        printf("Time Shift at Frame %d: %lf of %lf secs\n",
                irec, tdiff1, *shiftval);
        return (1);
    } else {
        return (0);
    }
}


/* -------------------------------------------------------------- */
/* sohHdrError() - returns > 0 if bytes do not match SOH header   */

/* -------------------------------------------------------------- */
BYTE sohHdrError(BYTE hdr[]) {
    BYTE status = 0;
    BYTE msghdr[6] = {225, 1, 193, 28, 0, 3};
    int i;

    for (i = 0; i < 6; i++)
        if (hdr[i + 2] != msghdr[i])
            status = 1;

    /*
    if (status)
      printf("Bad SOH header: %3d %3d %3d %3d %3d %3d \n",
             hdr[2],hdr[3],hdr[4],hdr[5],hdr[6],hdr[7]);
     */

    return (status);
}


/* -------------------------------------------------------------- */
/* bitError() - computes bit errors based on image data           */
/*                                                                */
/* Returns 1 if any errors found.                                 */

/* -------------------------------------------------------------- */
BYTE bitError(BYTE mnf[], INT32 *totbits, INT32 *toterrs) {
    INT32 numbits;
    INT32 numerrs;

    startBitError(mnf, &numbits, &numerrs);
    *totbits = numbits;
    *toterrs = numerrs;

    stopBitError(mnf, &numbits, &numerrs);
    *totbits += numbits;
    *toterrs += numerrs;

    if (*toterrs > 0)
        return (1);
    else
        return (0);
}


/* -------------------------------------------------------------- */
/* startBitError() - computes bit errors based on image data      */
/*                                                                */
/* Uses the known pattern of the start-pixel field to compute the */
/* number of bits which deviate from expectation.  There is one   */
/* start-pixel in LAC and 5 in GAC, so the test is more           */
/* comprehensive for GAC frames. Returns 1 if any errors found.   */

/* -------------------------------------------------------------- */
BYTE startBitError(BYTE mnf[], INT32 *numbits, INT32 *numerrs) {
    INT16 scid[2];
    INT16 mnftype;
    BYTE mask[] = {3, 255, 0, 0, 3, 255, 0, 0, 3, 255, 0, 0, 3, 255};
    BYTE npix;
    INT16 offset = 0; /* byte offset to first start-pixel */
    BYTE m, x, z;
    int i, j;

    memcpy(scid, &mnf[O_SCID], sizeof (scid));
    if (endianess() == 1)
        swapc_bytes((char *) scid, 2, 2);
    mnftype = scid2mnftype(scid);

    if (mnftype == GACTYPE) {
        npix = 5;
        *numbits = 350;
        offset = 806;
    } else if (mnftype == LACTYPE) {
        npix = 1;
        *numbits = 70;
        offset = 894;
    } else {
        npix = 0;
        *numbits = 0;
    }
    *numerrs = 0;

    for (i = 0; i < npix; i++) {
        for (j = 0; j < 14; j++) {
            x = mnf[offset + (4032 * i) + j];
            m = mask[j];
            if (x != m) {
                z = x ^ m;
                while (z) {
                    if (z & 0x1)
                        (*numerrs)++;
                    z = z >> 1;
                }
            }
        }
    }

    if (*numerrs > 0)
        return (1);
    else
        return (0);
}


/* -------------------------------------------------------------- */
/* stopBitError() - computes bit errors based on image data       */
/*                                                                */
/* Uses the known pattern of the stop-pixel field to compute the  */
/* number of bits which deviate from expectation.  There is one   */
/* start-pixel in LAC and 5 in GAC, so the test is more           */
/* comprehensive for GAC frames. Returns 1 if any errors found.   */

/* -------------------------------------------------------------- */
BYTE stopBitError(BYTE mnf[], INT32 *numbits, INT32 *numerrs) {
    INT16 scid[2];
    INT16 mnftype;
    BYTE mask[] = {3, 255, 3, 255, 0, 0, 0, 0, 3, 255, 3, 255, 0, 0};
    BYTE npix;
    INT16 offset = 0; /* byte offset to first start-pixel */
    BYTE m, x, z;
    int i, j;

    memcpy(scid, &mnf[O_SCID], sizeof (scid));
    if (endianess() == 1)
        swapc_bytes((char *) scid, 2, 2);
    mnftype = scid2mnftype(scid);

    if (mnftype == GACTYPE) {
        npix = 5;
        *numbits = 350;
        offset = 806 + 2000 * 2;
    } else if (mnftype == LACTYPE) {
        npix = 1;
        *numbits = 70;
        offset = 894 + 10296 * 2;
    } else {
        npix = 0;
        *numbits = 0;
    }
    *numerrs = 0;

    for (i = 0; i < npix; i++) {
        for (j = 0; j < 14; j++) {
            x = mnf[offset + (4032 * i) + j];
            m = mask[j];
            if (x != m) {
                z = x ^ m;
                while (z) {
                    if (z & 0x1)
                        (*numerrs)++;
                    z = z >> 1;
                }
            }
        }
    }

    if (*numerrs > 0)
        return (1);
    else
        return (0);
}

/* -------------------------------------------------------------- */
/* pixVariance() - computes pixel-to-pixel change                 */
/*                                                                */

/* -------------------------------------------------------------- */
INT32 pixVariance(BYTE mnf[]) {
    INT16 scid[2];
    INT16 mnftype;
    INT16 npix;
    INT16 offset;
    INT16 i, j;
    INT16 *data;
    INT32 var = 0;
    INT16 x, lastx = 0, dx;

    memcpy(scid, &mnf[O_SCID], sizeof (scid));
    if (endianess() == 1)
        swapc_bytes((char *) scid, 2, 2);
    mnftype = scid2mnftype(scid);

    if (mnftype == GACTYPE) {
        return (var);
    } else if (mnftype == LACTYPE) {
        npix = 1285;
        offset = 926;
        data = (INT16 *) & mnf[offset];
    } else {
        return (var);
    }

    for (j = 0; j < NBANDS; j++) {
        for (i = 0; i < npix; i++) {
            x = data[i * NBANDS + j];
            if (endianess() == 1)
                swapc_bytes((char *) &x, 2, 1);
            if (i > 0) {
                dx = x - lastx;
                var = var + (INT32) abs(dx);
            }
            lastx = x;
        }
    }

    return (var);
}









