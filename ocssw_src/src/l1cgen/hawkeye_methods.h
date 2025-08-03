#ifndef _HAWKEYE_METHODS_H_
#define _HAWKEYE_METHODS_H_

#include <string.h>
#include <string>

typedef float quat_array[4];
typedef double quat_array2[4];
typedef double orb_array2[3];
typedef double orb_array[3];

int orb_interp2(size_t n_SC_rec, size_t sdim, double *torb, orb_array2 *p, orb_array2 *v, double *time,
                orb_array2 *posi, orb_array2 *veli);
int j2000_to_ecr(int32_t iyr, int32_t idy, double sec, double ecmat[3][3]);
int j2000_to_mod(int32_t iyr, int32_t idy, double sec, double j2mod[3][3]);
int get_nut(int32_t iyr, int32_t idy, double xnut[3][3]);
int get_ut1(int32_t iyr, int32_t idy, double &ut1utc);
int ephparms(double t, double &xls, double &gs, double &xlm, double &omega);
int nutate(double t, double xls, double gs, double xlm, double omega, double &dpsi, double &eps,
           double &epsm);
int gha2000(int32_t iyr, double day, double &gha);

int expandEnvVar(std::string *sValue);

int orb_to_latlon(size_t ix_swt_ini,size_t ix_swt_end,size_t num_gridlines, int nbinx, double *orb_time_tot,orb_array2 *p, orb_array2 *v, double mgv1, double *tmgv1, double *tmgvf, float **lat_gd,float **lon_gd, float **alt,int FirsTerrain, bool swathCrossesUtcDay);
void cross_product_double2(double vector_a[], double vector_b[], double temp[]);
double cross_product_norm_double2(double vector_a[], double vector_b[]);
int interp_gap(size_t n_orb_rec, int *geogap,double *torb, double *latorb, double *lonorb);
#endif  // _GEOLOCATE_HAWKEYE_H_
