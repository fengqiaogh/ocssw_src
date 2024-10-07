#ifndef CALPROTOI_H /* avoid re-inclusion */
#define CALPROTOI_H

#include "l1a.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void read_caltable(char *cal_path);

int32_t l1b_rad(int syear, int sday, int32_t smsec, int32_t msec,
        char *dtype, int32_t nsta, int32_t ninc, int32_t npix,
        float *dark_mean, short *gain, short *tdi,
        short *scan_temp, float *inst_temp, int mside,
        short *l1a_data, float *l1b_data, cal_mod_struc *cal_mod);

int32_t get_cal(char *cal_path, int16_t syear, int16_t sday, int16_t eday,
        int32_t msec, int32_t npix, int32_t nsta, int32_t ninc,
        char *dtype, int16_t *tdi, int16_t *cal_year, int16_t *cal_day,
        int16_t *ref_year, int16_t *ref_day, int16_t *ref_min,
        float fp_temps[256][8], float scan_mod[2][1285],
        double *tfactor_const, double *tfactor_linear_1,
        double *tfactor_exponential_1, double *tfactor_linear_2,
        double *tfactor_exponential_2, double *cal_offset,
        double *inst_tcorr, double *inst_tref, double *fp_tcorr,
        double *fp_terf, double *mside1_const,
        double *mside1_linear_1, double *mside1_exponential_1,
        double *mside1_linear_2, double *mside1_exponential_2,
        double *mside2_const, double *mside2_linear_1,
        double *mside2_exponential_1, double *mside2_linear_2,
        double *mside2_exponential_2, float counts[8][4][5],
        float rads[8][4][5]);

int32_t read_parm_data(int32_t fid, int32_t sdfid, int32_t index, int32_t idoffs[8][16],
        float gains[8][16], float fp_temps[256][8],
        float scan_mod[2][1285], double *tfactor_const,
        double *tfactor_linear_1, double *tfactor_exponential_1,
        double *tfactor_linear_2, double *tfactor_exponential_2,
        double *cal_offset, double *inst_tcorr,
        double *inst_tref, double *fp_tcorr, double *fp_tref,
        double *mside1_const, double *mside1_linear_1,
        double *mside1_exponential_1, double *mside1_linear_2,
        double *mside1_exponential_2, double *mside2_const,
        double *mside2_linear_1, double *mside2_exponential_1,
        double *mside2_linear_2, double *mside2_exponential_2,
        int16_t tdi_list[256][4]);

void calc_knees(int16_t *tdi, int16_t tdi_list[256][4], int32_t idoffs[8][16],
        float gains[8][16], float counts[8][4][5],
        float rads[8][4][5]);

//int32_t read_time_data(int32_t fid, int16_t *year, int16_t *day, int32_t *msec,
//        int32_t *elts);

int32_t get_index(int32_t fid, int16_t syear, int16_t sday, int16_t eday, int32_t msec,
        int16_t *cal_year, int16_t *cal_day);

void setup_scanmod(int32_t npix, int32_t nsta, int32_t ninc,
        float scan_mod[2][1285]);

void sort_srads(float *srads, int32_t *oindex);

int32_t get_ref_time(int32_t sdfid,
        int16_t *ref_year, int16_t *ref_day, int16_t *ref_min);


int32_t get_tindex(int32_t fid, int16_t syear, int16_t sday, int16_t eday, int32_t msec,
        int16_t *cal_year, int16_t *cal_day);

#ifdef __cplusplus
}
#endif

#endif  /* CALPROTOI_H */
