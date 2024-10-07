#ifndef CALL1A_H
#define CALL1A_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifndef L1A_H
#include "l1a.h"
#endif
#include "get_cal.h"
#include "call1a_proto.h"
#include "getcal_proto.h"

#define  NBANDS 8
#define  NGAINS 4

extern float cal_counts[BANDS_DIMS_1A][GAINS_DIMS_1A][KNEES_DIMS_1A]; /*digital cnts (zero-offs corrected */
extern float cal_rads[BANDS_DIMS_1A][GAINS_DIMS_1A][KNEES_DIMS_1A]; /*radiances corresponding to knees  */

#endif /*CALL1A_H*/
