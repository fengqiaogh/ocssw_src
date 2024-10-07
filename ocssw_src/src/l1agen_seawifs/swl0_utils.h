#ifndef _SWL0_UTILS_H
#define _SWL0_UTILS_H

#include "swl0_types.h"
#include "swl0_struc.h"

INT16 scid2mnftype(INT16 scid[]);
INT16 scid2mnfnum(INT16 scid[]);
void ttag2ydmsec(INT16 ttag[], INT16 *year, INT16 *day, INT32 *msec);
FLOAT64 ttag2unix(INT16 ttag[]);
char *unix2timeStr(FLOAT64 usec);
BYTE timeError(swl0indx *indx, INT32 irec);
BYTE timeSeqError(swl0indx *indx, INT32 irec);
BYTE timeContiguous(swl0indx *indx, INT32 irec);
BYTE timeConsistent(swl0indx *indx, INT32 irec);
BYTE timeShifted(swl0indx *indx, INT32 irec, FLOAT64 *shiftval);
BYTE sohHdrError(BYTE hdr[]);
BYTE startBitError(BYTE mnf[], INT32 *numbits, INT32 *numerrs);
BYTE stopBitError(BYTE mnf[], INT32 *numbits, INT32 *numerrs);
BYTE bitError(BYTE mnf[], INT32 *numbits, INT32 *numerrs);
INT32 pixVariance(BYTE mnf[]);

/* Macro definitions */
#ifndef MAX
#define MAX(A,B)    ((A) > (B) ? (A) : (B))  /* Greater of (A,B) */ 
#endif

#ifndef MIN
#define MIN(A,B)    ((A) < (B) ? (A) : (B))  /* Lesser  of (A,B) */
#endif

#define ABS(A)      ((A) > 0 ? (A) : -(A))   /* Absolute Value   */ 



#endif
