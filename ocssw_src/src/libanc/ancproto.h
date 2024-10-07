#ifndef ANCPROTO_H     /* avoid re-inclusion */
#define ANCPROTO_H

#include <stdint.h>

typedef struct time_change {
    double start;
    double inc;
} timech;


int get_ancillary(float *in_lat, float *in_lon, int16_t cnt, int16_t syear,
        int16_t sday, int16_t eday, int32_t time, char *filename1,
        char *filename2, char *filename3,char *par_anc_file, char *anc_cor_file,
        int16_t parm_flag, float *interp, float *anc_unc,
        int16_t *qcflag);

int32_t set_files(int16_t parm_flag, int16_t syear, int16_t sday, int16_t eday,
        int32_t msec, char *filename1, char *filename2, char *filename3,
        char *file1, char *file2, timech *dtime1, timech *dtime2, int *toms);

int32_t ck_files_in_buf(int16_t PA_flag, int16_t parm_flag, char *f1, char *f2,
        int16_t month, int16_t *error_flag, int16_t *read_flag);

int read_climatology(char *file1, int16_t parm_flag, int16_t month,
        int32_t *dims, float *lat_buf,
        float *lon_buf, void *parm_buf);

int check_on_TOMS(int16_t parm_flag, char *file1, char *file2,
        char *filename1, char *filename2, char *filename3, double s_jd1,
        double s_jd2, double s_jd3, double e_jd1, double e_jd2,
        double e_jd3, double d_jd, timech *dt1, timech *dt2);

int read_NRT(char *file1, char *file2, char *anc_cor_file,
        int16_t parm_flag, int16_t read_flag, int32_t *dims, float *lat_buf,
        float *lon_buf, void *parm_buf1, void *parm_buf2,
        int8_t *qc_buf1, int8_t *qc_buf2);

void resiz_anc(int16_t *data, int8_t *qc, int32_t nlat, int32_t nlon1, int32_t nlon2,
        float * lons1, float * lons2);

void extract_data_pts(int16_t PA_flag, int16_t parm_flag, timech DTime1, timech DTime2,
        int16_t nsamp, float *in_latlist, float *in_lonlist,
        float *lat_bufp, float *lon_bufp, int32_t *dims, int8_t *qc_buf1,
        int8_t *qc_buf2, void *parm_buf1, void *parm_buf2, int toms,
        float *out_lat_list, float *out_lon_list, int16_t *qcflag,
        float *interp, float *anc_unc);

void gregor(int16_t jday, int16_t year, int16_t *month, int16_t *day);

void interpolate(int16_t PA_flag, int16_t parm_flag, double DT1, double DT2,
        float in_lat, float in_lon, float *lat_list,
        float *lon_list, void *data_p1, void *data_p2,
        int8_t *qc1, int8_t *qc2, float *intpdata,
        float *anc_unc, int32_t *int_qc);

int get_clim_data(int32_t fid, int32_t sdfid, int16_t parm_flag,
        int16_t month, int32_t *dims,
        float *lat_buf, float *lon_buf, void *parm_buf);

int get_NRT_data(int32_t fid, int32_t sdfid, char *anc_cor_file,
        int16_t parm_flag, int32_t *dims, float *lat_buf,
        float *lon_buf, void *parm_buf, int8_t *qc_buf);

int openHDF(char *infile, int32_t *sdfid, int32_t *fid);

int rdlatlon(int32_t fid, int32_t *dims, float *lat_buf,
        float *lon_buf);

int get_refs(int32_t fid, int32_t sdfid, int32_t vid, char *parm_label,
        char *QC_label, int32_t *geom_id, int32_t *parm_sdsid,
        int32_t *QC_sdsid, int32_t *dims);

int32_t anc_get_time(char *filename, double *s_jd, double *e_jd);

int closeHDF(int32_t sdfid, int32_t fid);

#endif /* ANCPROTO_H */
