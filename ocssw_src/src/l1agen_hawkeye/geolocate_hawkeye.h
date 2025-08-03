#ifndef _GEOLOCATE_HAWKEYE_H_
#define _GEOLOCATE_HAWKEYE_H_

#include <genutils.h>

typedef float quat_array[4];
typedef float orb_array[3];

float constexpr focal_length = 45.184;

int j2000_to_ecr(int32_t iyr, int32_t idy, double sec, double ecmat[3][3]);
int j2000_to_mod(int32_t iyr, int32_t idy, double sec, double j2mod[3][3]);
int get_nut(int32_t iyr, int32_t idy, double xnut[3][3]);
int get_ut1(int32_t iyr, int32_t idy, double &ut1utc);
int ephparms(double t, double &xls, double &gs, double &xlm, double &omega);
int nutate(double t, double xls, double gs, double xlm, double omega,
           double &dpsi, double &eps, double &epsm);
int gha2000(int32_t iyr, double day, double &gha);
int mtoq(double rm[3][3], double q[4]);
int qprod(double q1[4], float q2[4], double q3[4]);
int qprod(float q1[4], float q2[4], float q3[4]);
int orb_interp(size_t n_SC_rec, size_t sdim,
               double *torb, orb_array *p, orb_array *v,
               double *time, orb_array *posi, orb_array *veli);
int q_interp(size_t n_SC_rec, size_t sdim, double *tq, quat_array *q,
             double *time, quat_array *qi);
int l_sun(size_t sdim, int32_t iyr, int32_t iday,
          double *sec, orb_array *sunr);
int sun2000(size_t sdim, int32_t iyr, int32_t idy,
            double *sec, orb_array *sun);
int qtom(float quat[4], double rm[3][3]);
int scan_ell(float p[3], double sm[3][3], double coef[10]);
int uni_geonav(float pos[3], float vel[3], double smat[3][3], double coef[10],
               float sunr[3], orb_array *pview, size_t npix,
               float *xlat, float *xlon, short *solz, short *sola,
               short *senz, short *sena, short *range);
int mat2rpy(float pos[3], float vel[3], double smat[3][3], float rpy[3], double om[3][3]);
int euler(float a[3], double xm[3][3]);

#endif  // _GEOLOCATE_HAWKEYE_H_
