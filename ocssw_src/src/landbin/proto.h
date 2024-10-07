#include "mfhdf.h"

void chand(float xphi, float xmuv, float xmus, float *xtau, float *rhoray, double *trup, double *trdown, int nbands, char *process);
int compute_l1b(int32 sd_id, int iline, float *l1b_data[], char *calibtable, char newfile);
float csalbr(float xtau);
void make_psurf(signed short *psurf, char *file, int Nlines, int Npixels, float psurf0);
int get_attributes(int32 sd_id, char *prefix, int *Nlines, int *Npixels, int *startpix, int *subsampling, int *year, int *doy, int *orbno, char *title);
int get_calib_sds(int32 sd_id, int16 *dark_rest, int16 *gain, int16 *tdi, int16 *scan_temp, int16 *side, int32 *msec);
int get_navig_sds(int32 sd_id, float32 *sen_mat, float32 *orb_vec, float32 *sun_ref, float32 *scan_ell, int32 *msec);
int get_navig_sds_line(int32 sd_id, int iline, float32 *sen_mat, float32 *orb_vec, float32 *sun_ref, float32 *scan_ell, char newfile);
int read_attr(int32 sd_id, char *attrname, void *attr);
int32 open_sds(int32 sd_id, char *sds_name);
int read_sds(int32 sd_id, char *sds_name, void *dataP);
int read_sds_block(int32 sds_id, int iline, int buflines, void *dataP);
int write_sds(int32 sd_id, char *SDSname, void *data, int startline, int buflines, int Nl, int Np, int32 numtype, void *fillvalue, float64 scale_factor, char verbose, char gzip);
