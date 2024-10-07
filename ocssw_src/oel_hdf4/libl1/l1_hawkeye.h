#ifndef L1_HAWKEYE_H
#define L1_HAWKEYE_H

#ifdef __cplusplus
extern "C" {
#endif
    
#include <stdint.h>
#include "l1.h"

typedef struct hawkeye_t {
    int exposureID;
    float roll_offset;
    float time_offset;
} hawkeye_t;

// Call sequence
// openl1_hawkeye      > qc_hawkeye_CCD_T      > nan_wmean                     // once
// readl1_hawkeye      > read_cal_hawkeye                                      // once
// readl1_hawkeye      > interp_hawkeye_CCD_T  > prep_for_interp_double        // per scan
// readl1_hawkeye      > calibrate_hawkeye     > prep_for_interp_double        // per scan

// Open, read, close
int openl1_hawkeye(filehandle *file);
int readl1_hawkeye(filehandle *file, int32_t recnum, l1str *l1rec);
int closel1_hawkeye(filehandle *file);

    
#ifdef __cplusplus
}
#endif
    
#endif
