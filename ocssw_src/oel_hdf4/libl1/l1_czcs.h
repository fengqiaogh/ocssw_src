/*
 *  W. Robinson, SAIC, 10 Dec 2004  new for CZCS
 */
#ifndef L1_CZCSW_H
#define L1_CZCSW_H

#define   NREC_IN_BUF             1
#define   NBND_CZCS               5
#define   POS_ERR_THRESH          2000.  /* orbit position error tolerance */

int openl1_czcs(filehandle *file);
int readl1_czcs(filehandle *file, int32_t recnum, l1str *l1rec);
int closel1_czcs(filehandle *file);

int czcs_ring(int gain, float *lt750, char *ring_sat, l1str *l1rec);
int get_czcscal(char *file, int orbit, int16_t year, int16_t day, int msec, short l1acnt[], float slope750, float intercept750, int16_t igain, float l1brads[]);
int cz_posll_2_satang(float *, int, float *, float *, float *, float *);
void matrix_mult(double[3], double[3][3], double[3]);
void cross_prod(double *, double *, double *);

#endif
