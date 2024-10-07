/* this file is a partial copy of profile_utils.c/.h from UW-M 1km CT code */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "profile_management.h"

#define GRAV_ACCEL_EARTH	9.80665			/* m / s^2	*/
#define SPEC_GAS_CONST_AIR	287.05			/* J / kg K	*/
#define SPEC_GAS_CONST_H2O	461.51			/* J / kg K	*/

/*******************************************************************************
 /
 *******************************************************************************/
double lin_int(double x1, double x2, double x, double y1, double y2) {
	
	return (x - x1) / (x2 - x1) * (y2 - y1);
}


/*******************************************************************************
/
*******************************************************************************/

int profile_to_101_(double *p, double *t, double *w, int *nn, double *nlat,
                   double *pp, double *tt, double *ww, int *is_o3) {
     int i;

     int iy;
     int iz;
     int nlx = 101;
     
     int i_str;
     int i_end;

     double anum;
     double aden;

     double rlogp;

     double delt;
     double delw;

     double wmin;
     double wmax;

     double pb[2];
     double tb[2];
     double wb[2];
	 
	 int n;
	 double lat;
	 
	 n = *nn;
	 lat = *nlat;

     if (n < 2) {
          printf("number of levels must be atleast 2\n");
          return -1;
     }

	if (*is_o3 == 1) 
		wmin = 0.;
	else
		wmin = 0.003;


     /*make_profile_101(pp);*/
     int_levels_pp(p, t, w, n, pp, tt, ww, nlx, &i_str, &i_end);

     if (i_str >= 35) {
          printf("ERROR: temperature profile doesn't go high enough\n");
          return -1;
     }

     for (i = 0; i < i_str; ++i) {
          tt[i] = -1.;
          ww[i] = wmin;
     }

     extem101_64_(tt, &lat);
     
     iz = n - 1;

     if (p[iz] < pp[nlx-1]) {
          pb[0] = p[iz];
          tb[0] = t[iz];
          wb[0] = w[iz];

          pb[1] = pp[nlx-1];

          iy = n - 2;
/*
          while (p[iy] >= p[iz])
               --iy;
*/
          anum = log(pb[1] / pb[0]);
          aden = log(pb[1] / p[iy]);

          rlogp = anum / aden;

          delt = t[iz] - t[iy];
          tb[1] = t[iz] + delt * rlogp;

          delw = w[iz] - w[iy];
          wb[1] = w[iz] + delw * rlogp;

          if (wb[1] < wmin)
               wb[1] = wmin;
          else {
               wmax = get_satmix(pb[1], tb[1]);

               if (wb[1] > wmax)
                    wb[1] = wmax;
          }

          i = i_end + 1;

          int_levels_pp(pb, tb, wb, 2, pp+i, tt+i, ww+i, nlx-i, &i_str, &i_end);
     }

     return 0;
}

float get_satmix(float P, float T) {


/*! After Rogers and Yau.*/


    const float c = 4187.;
	const float cpv = 1870. ;
	const float L0 = 2.501e6 ;
	const float T0 = 273.;

	float L, esat, ws;
  
	if (T < 235.) 
		L = 2.83e6;   /*! latent heat of ice, doesnt change very much */
    else
		L = L0 - (c - cpv) * (T - T0); /* ! latent heat of water */

	esat = 6.11657*exp((L/461.51)*(1/273. - 1/T) ); /*! saturation vapor pressure in mb */

	ws = esat * 0.622 / (P - esat) * 1000.; /*! saturation mixing ratio in g/kg */

	return ws;

}




/*******************************************************************************
/
*******************************************************************************/
int make_profile_101_(double *p) {

     int i;

     double l;

     double a = -1.550789414500298e-04;
     double b = -5.593654380586063e-02;
     double c =  7.451622227151780e+00;

     l = 101;
     for (i = 0; i < 101; ++i) {
            p[i] = pow((a*l*l + b*l + c), (7./2.));
            l = l - 1;
     }

     return 0;
}

/*******************************************************************************
/ 
*******************************************************************************/
int int_levels_pp(double *p,  double *t,  double *w, int n,
                  double *pp, double *tt, double *ww, int nn, int *i_str, int *i_end) {

     int i;
     int ii;

     double dl;

     double slope_t;
     double slope_w;

     if (n < 2) {
          *i_str =  0;
          *i_end = -1;

          return 0;
     }

     for (i = 0; i < nn; ++i) {
          if (pp[i] >= p[0])
               break;

          tt[i] = 0.;
          ww[i] = 0.;
     }

     *i_str = i;

     ii = 1;

     dl    = log(p[ii] / p[ii-1]);
     slope_t = (t[ii] - t[ii-1]) / dl;
     slope_w = (w[ii] - w[ii-1]) / dl;

     for (     ; i < nn;    ) {
          if (pp[i] < p[ii]) {
               dl      = log(pp[i] / p[ii-1]);
               tt[i] = t[ii-1] + slope_t * dl;
               ww[i] = w[ii-1] + slope_w * dl;

               ++i;
          }
          else {
               ++ii;
               if (ii >= n) {
                    if (p[ii-1] == pp[i]) {
                         tt[i] = t[ii-1];
                         ww[i] = w[ii-1];

                         ++i;
                    }

                    break;
               }

               dl    = log(p[ii] / p[ii-1]);
               slope_t = (t[ii] - t[ii-1]) / dl;
               slope_w = (w[ii] - w[ii-1]) / dl;
          }
     }

     *i_end = i - 1;

     for (     ; i < nn; ++i) {
          tt[i] = 0.;
          ww[i] = 0.;
     }

     return 0;
}

/*******************************************************************************
/
*******************************************************************************/
int height_profile_(double *p, double *t, double *w, double *z, int *nn, double *np0) {

     int i;

     int i_up;
     int i_dn;

     double g;

     double Rd;
     double Rv;

     double epsilon;

     double t0;
     double w0;

     double z_cur;

     double p_last;
     double t_last;
     double w_last;

     double P;
     double T;
     double W;

     double e;

     double rho;
	 
	 int n;
	 double p0;
	 
	 n = *nn-1; /* we're calling from fortran */
	 p0 = *np0;

     g = GRAV_ACCEL_EARTH;

     Rd = SPEC_GAS_CONST_AIR;
     Rv = SPEC_GAS_CONST_H2O;

     epsilon = Rd / Rv;


     for (i = 0; i < n; ++i) {
          if (p[i] > p0)
               break;
     }

     if (i == 0) {
          i_up = i;
          i_dn = i+1;
     }
     else {
          i_up = i-1;
          i_dn = i;
     }

     t0 = t[i_up] + lin_int(p[i_up], p[i_dn], p0, t[i_up], t[i_dn]);
     w0 = w[i_up] + lin_int(p[i_up], p[i_dn], p0, w[i_up], w[i_dn]);

     z_cur = 0;

     p_last = p0;
     t_last = t0;
     w_last = w0;
     for (i = i_up+1; i >= 0; --i) {
          P = (p_last + p[i]) / 2. * 100.;
          T = (t_last + t[i]) / 2.;
          W = (w_last + w[i]) / 2. / 1000.;

          e = W / (W + epsilon) * P;

          rho = (P-e) / (Rd*T) + e / (Rv*T);

          z_cur += (p_last-p[i])*100. / (g*rho);

          z[i] = z_cur;

          p_last = p[i];
          t_last = t[i];
          w_last = w[i];
     }


     z_cur = 0;

     p_last = p0;
     t_last = t0;
     w_last = w0;
     for (i = i_dn; i < n; ++i) {
          P = (p_last + p[i]) / 2. * 100.;
          T = (t_last + t[i]) / 2.;
          W = (w_last + w[i]) / 2. / 1000.;

          e = W / (W + epsilon) * P;

          rho = (P-e) / (Rd*T) + e / (Rv*T);

          z_cur += (p_last-p[i])*100. / (g*rho);

          z[i] = z_cur;

          p_last = p[i];
          t_last = t[i];
          w_last = w[i];
     }
     
     
     
     return 0;
}
