#ifndef STATPROTO_H
#define STATPROTO_H

extern int32_t read_cntldata(char *cntl_file, cntl_str *databinchk,
        cntl_str *grschk, cntl_str *statchk, int32_t *wtchk,
        clim_str *climchk);

extern int32_t l3file(int32_t sdfid, int32_t c_sdfid, int32_t *nbins, int32_t *c_nbins,
        char *ptype);

extern int32_t chk_databin(int32_t sdfid, cntl_str *databinchk);

extern int32_t l3data_chk(char *clim_file, int32_t fid, int32_t nbins, int32_t c_nbins,
        int32_t c_fid, cntl_str *grschk, cntl_str *statchk,
        clim_str *climchk);

extern void stat_exit(int status);

extern int32_t chk_weight(int32_t fid, int32_t nbins);

extern int32_t get_wtcnt(int32_t vsid, char *vsname, int32_t start, int32_t elts,
        float *wt_buf);

extern double calc_sd(double xmean, float wts, int16_t nseg, float sumxx);

extern void pr_grs_results(cntl_str *grschk, int32_t nbins, double *value1,
        double *value2);

extern void pr_stat_results(cntl_str *statchk, int32_t nbins, double *sumxbar,
        double *sumsd, double *, double *);

extern void pr_clim_results(char *clim_file, clim_str *climchk, double npts,
        double *num1Hstd, double *num1Lstd, double *num2Hstd,
        double *num2Lstd, double *num3Hstd, double *num3Lstd);

extern int32_t pr_error(char *label, int32_t nvals, int32_t nvals_read);

extern int32_t rdattr(int32_t sdfid, char *attr_name, void *buf);

extern int32_t get_vsid(int32_t fid, char *vs_name);

#endif /* STATPROTO_H */
