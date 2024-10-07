#ifndef STATPROTO_H
#define STATPROTO_H
#include "usrmac.h"
extern int32 read_cntldata
PROTO((char *control_file, char *fsttim, cntl1_str *gn1,
        cntl1_str *gn2, cntl1_str *zero, cntl2_str *l1hicnt,
        cntl2_str *l1locnt, float32 *nav_thresh1,
        float32 *nav_thresh2, float32 *l1tilt_thresh,
        float32 *pct_noise_thresh, float32 *pct_encrypt_thresh,
        thr_ctl_def *thr_ctl, int16 *rpt_negtim));

extern int32 l1file
PROTO((int32 sdfid, int32 *nsamp, int32 *nscans, int16 *dtynum));

extern int32 chk_gn
PROTO((int32 sdfid, cntl1_str *gn1, cntl1_str *gn2, int16 dtynum,
        int32 nsamps, int32 nscans));

extern int32 chk_zero
PROTO((int32 sdfid, cntl1_str *zero_str, int32 nscans, int32 nsamp));

extern void stat_exit
PROTO((int status));

extern int32 rdattr
PROTO((int32 sdfid, char *attr_name, void *buf));

extern int32 rdslice
PROTO((int32 sdfid, char *name, int32 *start, int32 *edge, void *buf));

extern int32 chk_count
PROTO((int32 sdfid, int32 nscans, int32 nsamp, int16 dtynum,
        cntl2_str *l1hicnt, cntl2_str *l1locnt, int *spike_cnt,
        float *line_sd));

extern void get_hicnt
PROTO((int32 nrec, int32 nsamp, int32 nbands, int16 *databuf,
        cntl2_str *l1hicnt, int32 *hicnt));

extern void get_lowcnt
PROTO((int32 nrec, int32 nsamp, int32 nbands, int16 *databuf,
        cntl2_str *l1hicnt, int32 *lowcnt));

extern int32 chk_nav
PROTO((int32 sdfid, int32 nscans, float32 nav_thresh1,
        float32 nav_thresh2, int32 ntilts, int16 tilt_ranges[20][2],
        int16 tilt_flags[20], int16 dtynum, int16 rpt_negtim));

extern int32 chk_tilt
PROTO((int32 sdfid, int16 dtynum, float32 l1tilt_thresh,
        int32 *ntilts, int16 tilt_ranges[20][2], int16 tilt_flags[20]));

extern int geovex_(float32 *orb_vec, float32 *sen_mat, float32 *scan_ell,
        float32 *sun_ref, float32 *v0);

extern void ck_trng(char *);

extern void anal_noise(int32, int32, int32, int32, int32, int16 *, int *,
        float *, int32_t *, int32_t *);

extern void rpt_noise(int32, int16, int32, int32, int *, float *,
        float, float);

extern void chk_inst_ana(int32, int32, thr_ctl_def);

extern void chk_gainv(int32, int16, int32, thr_ctl_def);

extern void chk_tdiv(int32, int16, int32, thr_ctl_def);
#endif /* STATPROTO_H */
