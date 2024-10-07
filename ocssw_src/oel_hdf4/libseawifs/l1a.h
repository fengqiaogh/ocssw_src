
#ifndef L1A_H
#define L1A_H

#if 0
#ifdef GENERAL
#define get_l1a_open get_l1a_open_all
#define get_l1a_record get_l1a_record_all
#define get_l1a_close get_l1a_close_all
#define put_l1a_open put_l1a_open_all
#define put_l1a_record put_l1a_record_all
#define put_l1a_close put_l1a_close_all
#define set_l1aget_buffer_blksize set_l1aget_buffer_blksize_all
#define set_l1aput_buffer_blksize set_l1aput_buffer_blksize_all
#endif /* GENERAL */

#include "navigation.h"
#include "level_1a_index.h"

#ifndef MAX_HDF
#define MAX_HDF 16
#endif

#endif

#ifndef MAX_HDF_L1AGET
#define MAX_HDF_L1AGET 4
#endif

#include "level_1a_index.h"

/* out of band correction */
#define OOB_OFF   0
#define OOB_DEFAULT_METHOD 1
#define OXYGEN_CORR_FACTOR 1.12f
extern int out_band_corr(float *radiances, float oxygen_factor, int nsample);

/* calibration modification structure definition */
typedef struct cal_mod_def {
    int flag; /* use flag: 0 - use cal file values
                                   1 - use input gain, cal file offset
                                   2 - use cal file gain, input offset
                                   3 - use input gain and cal  */
    double gain[8]; /* calibration gain */
    double offset[8]; /* calibration offset */
} cal_mod_struc;

#endif /* L1A_H */
