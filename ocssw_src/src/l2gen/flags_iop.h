#ifndef _FLAGS_IOP_H
#define _FLAGS_IOP_H

#include <stdint.h>
#include "l12_parms.h"

/* flag bit settings */
#define IOPF_ISMASKED 0x0001   // pixel masked
#define IOPF_FAILED   0x0002   // algorithm signals failure
#define IOPF_MAXITER  0x0004   // max iterations reached
#define IOPF_BADRRS   0x0008   // insufficient valid Rrs
#define IOPF_NAN      0x0010   // inversion returned NAN
#define IOPF_RRSDIFF  0x0020   // mean Rrs diff exceeded threshold
#define IOPF_ALO      0x0040   // retrieved a below threshold
#define IOPF_AHI      0x0080   // retrieved a above threshold
#define IOPF_APHLO    0x0100   // retrieved aph below threshold
#define IOPF_APHHI    0x0200   // retrieved aph above threshold
#define IOPF_ADGLO    0x0400   // retrieved adg below threshold
#define IOPF_ADGHI    0x0800   // retrieved adg above threshold
#define IOPF_BBLO     0x1000   // retrieved bb below threshold
#define IOPF_BBHI     0x2000   // retrieved abb above threshold
#define IOPF_BBPLO    0x4000   // retrieved bbp below threshold
#define IOPF_BBPHI    0x8000   // retrieved bbp above threshold

#ifdef __cplusplus
extern "C" {
#endif

static const char * const giop_flag_lname[NGIOPFLAGS] = {"ISMASKED",
    "FAILED",
    "MAXITER",
    "BADRRS",
    "NAN",
    "RRSDIFF",
    "ALO",
    "AHI",
    "APHLO",
    "APHHI",
    "ADGLO",
    "ADGHI",
    "BBLO",
    "BBHI",
    "BBPLO",
    "BBPHI"};

typedef struct iopflagctl_struc {
    int32_t nwave;
    int32_t maxiter;

    float *a_lo;
    float *a_hi;
    float *a_on;

    float *aph_lo;
    float *aph_hi;
    float *aph_on;

    float *adg_lo;
    float *adg_hi;
    float *adg_on;

    float *bb_lo;
    float *bb_hi;
    float *bb_on;

    float *bbp_lo;
    float *bbp_hi;
    float *bbp_on;

} iopfstr;

void set_iop_flag(float wave[], int32_t nwave,
        float a[], float aph[], float adg[],
        float bb[], float bbp[], int16_t *flag);

#ifdef __cplusplus
}
#endif

#endif
