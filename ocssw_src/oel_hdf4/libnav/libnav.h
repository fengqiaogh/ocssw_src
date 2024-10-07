#ifndef LIBNAV_H
#define LIBNAV_H

/*
 * Function prototypes for routines defined in src/libnav
 */

#ifdef  __cplusplus
extern "C" {
#endif

    int crossp_(float *v1, float *v2, float *v3);
    void sunangs_( int *year, int *day, float *gmt, float *lon, float *lat,
                   float *sunz, float *suna);

    void get_zenaz(float *pos, float lon, float lat, float *zenith, float *azimuth);
    
    void compute_alpha(float lon[], float lat[],
                       float senz[], float sena[],
                       double mnorm[3], int npix, float alpha[]);

    void l_sun_(int *iyr, int *idoy, double *sec, float sunr[3], float *rs);

    void ocorient_(float *pos, float *vel, float *att, float (*)[3], float *coef);

    double esdist_(int *year, int *day, int *msec);

    void cdata_();


  void nav_get_vel(float ilat, float mlat, float ilon, float mlon, float *vel);
  void nav_get_pos(float lat, float lon, float alt, float *p);
  void nav_gd_orient(float *pos, float *vel, float *att, float *smat[], float *coef);
  void nav_get_geonav(float *sunr, float *pos, float *pview, float *coef, float *smat[], float *xlon, float *xlat, float *solz, float *sola, float *senz, float *sena);
  double angular_distance(double lat1, double lon1, double lat2, double lon2);
  double deg2rad(double deg);
  double rad2deg(double rad);


  
#ifdef  __cplusplus
}
#endif


#endif /* LIBNAV_H */
