#ifndef CALL1A_H
#define CALL1A_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "hdf.h"
#include "mfhdf.h"

typedef struct cal_mod_def {
    int flag; /* use flag: 0 - use cal file values
                                   1 - use input gain, cal file offset
                                   2 - use cal file gain, input offset
                                   3 - use input gain and cal  */
    double gain[8]; /* calibration gain */
    double offset[8]; /* calibration offset */
} cal_mod_struc;
/*
#ifndef L1A_H
#include "l1a_modif.h"
#endif
 */
#include "calib_get_cal.h"
#include "calib_call1a_proto.h"
#include "calib_getcal_proto.h"

#define  NBANDS 8
#define  NGAINS 4

#endif /*CALL1A_H*/
