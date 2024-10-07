#ifndef ST_PROTO_H
#define ST_PROTO_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
int32_t stray_light_gac(int32_t *initial, float Ltyp_frac, float Styp_frac,
        int32_t  nscans, int32_t  nsamples, int32_t scan_no, int16_t gn,
        float *rads, float *l1b_data, int32_t *sl_scan,
        int32_t *sl_flag);

int32_t stray_light_lac(int32_t *initial, float Ltyp_frac, float Styp_frac,
        int32_t nscans, int32_t nsamples, int32_t scan_no, int16_t gn,
        float *rads, float *l1b_data, int32_t *sl_scan,
        int32_t *sl_flag);

int32_t stray_light_corr(int32_t *initial, float Ltyp_frac, float Styp_frac,
        int32_t nscans, int32_t nsamples, int32_t scan_no, char *dtype,
        int16_t gn, float *rads, float *l1b_data, int32_t *sl_scan,
        int16_t *l2_flags, int32_t *AS_pixels);

#ifdef __cplusplus
}
#endif


#endif  /* ST_PROTO_H */

