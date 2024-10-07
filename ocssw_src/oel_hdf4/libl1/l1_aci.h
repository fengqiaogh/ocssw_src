#ifndef  _L1_ACI_H
#define  _L1_ACI_H

/* noaa avhrr satellites */
#define NO06  6    /* the 2nd character is a capital oh, the 3rd is a zero */   
#define NO07  7                                                                 
#define NO08  8                                                                 
#define NO09  9                                                                 
#define NO10 10                                                                 
#define NO11 11                                                                 
#define NO12 12                                                                 
#define NO14 14                                                                 
#define NO15 15                                                                 
#define NO16 16                                                                 
#define NO17 17                                                                 
#define NO18 18                                                                 
#define NO19 19                                                                 

#include <stdint.h>
#include "l1.h"

int closel1_aci(filehandle *l1file);
int openl1_aci(filehandle *l1file);
int readl1_aci(filehandle *l1file, int32_t recnum, l1str *l1rec);
const char* xsatid2name(int xsatid);
int satname2xsatid(const char* satname);


#endif
