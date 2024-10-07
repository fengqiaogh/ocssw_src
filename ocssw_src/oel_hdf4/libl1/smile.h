#ifndef  SMILE_H
#define  SMILE_H

#include "l1.h"

#ifdef __cplusplus
extern "C" {
#endif

void smile_init(int num_bands, int num_detectors, const char* bandinfo_filename,
        float* detectorWLs, float* detectorE0s);
void radcor(l1str *l1rec, int32_t ip, int32_t land, int32_t escorrected);

#ifdef __cplusplus
}
#endif

#endif
