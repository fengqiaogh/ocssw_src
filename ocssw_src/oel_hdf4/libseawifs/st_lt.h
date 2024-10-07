#ifndef ST_LT_H
#define ST_LT_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "st_proto.h"
#include "cf.h"

#define BT       0
#define AT       -1
#define DIAG     -2
#define BLANK    -10
#define NO       1 
#define RIGHT    2 
#define LEFT     3 

#ifdef TRUE
#undef TRUE
#endif

#define TRUE     1
#define FALSE    0
#define DONE     1
#define NOTDONE  0
#define MAXLINES        5
#define GACLINES 3
#define MAXSAMPS        1285
#define GACSAMPS 248
#define NBANDS  8
#define GAINS  4
#define KNEES  5 

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B)) 

#endif /* ST_LT_H */
