
#ifndef _L12_PARMS_H
#define _L12_PARMS_H

#include <sensorDefs.h>
#include <l1.h>

#define PROGRAM    "l2gen"

/* #define NBANDS        16 */
#define NQMIN          3
#define NQMAX        500
#define FILTMAX      200
#define MAX_OFILES    10
#define MAX_IFILES  1024
#define NSSTFLAGS     16
#define NGIOPFLAGS    16
#define NINPRODS       3
#define NQSSTFLAGS     5

#define MAXAERMOD     100
#define AERWHITE        0
#define AERRH          -1
#define AERRHNIR       -2
#define FIXMODPAIR     -4
#define FIXMODPAIRNIR  -5
#define FIXANGSTROM    -6
#define FIXANGSTROMNIR -7
#define FIXAOT         -8
#define AERRHSWIR      -9
#define AERRHMUMM      -16
#define AERRHMSEPS     -17
#define AERRHSM        -18
#define AERNULL        -99

#define DEFAULT_CHL     0
#define CHL_MIN      0.00
#define CHL_MAX     100.0
#define AOT_MIN      0.00
#define AOT_MAX       1.0

#define DEM_WIDTH  43200
#define DEM_HEIGHT 21600

#define NOBRDF    0     /* brdf  */
#define FRESNSEN  1     /* bit 1 */
#define FRESNSOL  2     /* bit 2 */
#define FOQMOREL  4     /* bit 3 */
#define DTBRDF    8     /* bit 4 */
#define QMOREL   16     /* bit 5 */

#define O3_BIT     1
#define CO2_BIT    2
#define NO2_BIT    4
#define H2O_BIT    8
#define ATREM_BIT 16
#define GAS_TRANS_TBL_BIT 32
#define CO_BIT 64
#define CH4_BIT 128
#define N2O_BIT 256

#define IOPNONE   0
#define IOPCARDER 1
#define IOPGSM    2
#define IOPQAA    3
#define IOPPML    4
#define IOPNIWA   5
#define IOPLAS    6
#define IOPGIOP   7
#define IOPSWIM   8
#define IOPDEFAULT IOPQAA

#define QAABLEND  0
#define QAA555    1
#define QAA640    2

#define NOMATCH_ERROR  110
#define FILESIZE_ERROR 111

#define STDPR 1013.25

#define DAYSCENE       0
#define NIGHTSCENE     1
#define DAYANDNIGHT    2
#define UNKNOWNSCENE   3

#define SWN  0
#define SWA  1
#define SWBB 2

#endif
