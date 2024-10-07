#ifndef _XCAL_H
#define _XCAL_H

#include "l1.h"

#define XTNTIME 300
#define XTNDET   40
#define XTNMSIDE  2
#define XTNORDER  6
#define XRVS      0
#define XM12      1
#define XM13      2

double *get_xcal(l1str *l1rec, int type, int bandnum);
double *get_fpm_xcal(char *fpm_file); //added by Sudipta to support FPM based band correction

//double *xcal_modis(l1str *l1rec, int type, int bandnum);

#endif
