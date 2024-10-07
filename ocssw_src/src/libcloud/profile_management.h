#ifndef _PROFILE_MANAGEMENT_H_
#define _PROFILE_MANAGEMENT_H_

int profile_to_101_(double *p, double *t, double *w, int *n, double *lat,
                   double *pp, double *tt, double *ww, int *is_o3) ;
				   
int make_profile_101_(double *p) ;
	
int int_levels_pp(double *p,  double *t,  double *w, int n,
                  double *pp, double *tt, double *ww, int nn, int *i_str, int *i_end);				   
							   			   

int height_profile_(double *p, double *t, double *w, double *, int *, double *);

float get_satmix(float, float);
double lin_int(double, double, double, double, double);

void extem101_64_(double *tt, double *slat);

#endif