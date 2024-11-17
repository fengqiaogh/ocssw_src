#ifndef NAV_H_
#define NAV_H_
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "stdbool.h"
#include "math.h"
#include <stdint.h>
#define DSIGN(A, B) B >= 0 ? fabs(A) : -fabs(A)
const static int32_t maxlin = 16000;
const static double pi = 3.141592653589793e0, radeg = 57.29577951e0, re_const = 6378.137e0, rem = 6371.e0,
                    f_const = 1.e0 / 298.257e0, omf2 = (1.e0 - f_const) * (1.e0 - f_const),
                    omegae = 7.29211585494e-5;  //   common /gconst/pi,radeg,re,rem,f,omf2,omegae

/**
 * @brief
 * computes the norm of a vector
 * @param vec vector
 * @param n dim
 * @return double
 */
double norm(const double* vec, size_t n);

/**
 * @brief
 * computes the square of the  norm of a vector
 * @param vec vector
 * @param n dim
 * @return double
 */
double square(const double* vec, size_t n);

/**
 * @brief
 * computes the scalar/dot product of two vectors
 * @param vec1 vector1
 * @param vec2 vector1
 * @param n dim
 * @return double
 */
double dot(const double* vec1, const double* vec2, size_t n);


/**
 * @brief Description:

Returns remainder calculated as:

     A - (INT(A / P) * P)
P must not be zero.
 *
 * @param a INTEGER or REAL; scalar; INTENT(IN).
 * @param p INTEGER or REAL; scalar; INTENT(IN).
 * @return double
 */
double dmod(double a, double p);


/**
 * @brief c    This function converts a calendar date to the corresponding Julian
day starting at noon on the calendar date.  The algorithm used is
from Van Flandern and Pulkkinen, Ap. J. Supplement Series 41,
November 1979, p. 400.  Written by Frederick S. Patt, GSC, November 4, 1992
 *
 * @param i Year - e.g. 1970; I
 * @param j Month - (1-12): I
 * @param k Day  - (1-31) : I
 * @return int32_t  Julian day; O
 */
int32_t jd(int32_t i, int32_t j, int32_t k);


/**
 * @param tin     TIN(1) = CALENDER DATE, YYYYMMDD.
       (2) = CALENDER TIME, HHMMSS.SSSS
 * @return  TOUT   = JULIAN DATE (DAYS)
 */
double julian(double tin[2]);

/**
c       Compute days since January 0, 1900
 * @param jd julian date
 * @param i output year
 * @param k output day of year (0-365/366)
 */
void jdate(int jd, int *i, int *k);


/**
 * @brief Construct a new ddate objectC
  C       This routine computes the calendar date corresponding to
  C       a given Julian day.  This code was brazenly copied from
  C       a routine written by Myron Shear for CSC on Julian Day 1.
  C       
  C       ARGUMENT        TYPE    I/O     DESCRIPTION     
  C       __________________________________________________________
  C        JD             I*4      I      Julian Day (reference Jan 1, 4713 BC)
  C        I              I*4      O      Year 
  C        J              I*4      O      Month   
  C        K              I*4      0      Day of Month
  C
 * 
 * @param jd 
 * @param i 
 * @param j 
 * @param k 
 */
void jddate(int jd, int *i, int *j, int *k);


/**
 * @brief c -------------------------------------------------------------
  c Subroutine ymdhms2jul
  c
  c Convers from Year, Month, Day of Month, Hour, Minute, Second
  c to Julian time.
  c
  c BA Franz, GSC, 1/97
 * 
 * @param year input year
 * @param month input month
 * @param day input day
 * @param hour input hour
 * @param minute input minute
 * @param sec input second
 * @param jul output julian date
 */
void ymdhms2jul(int32_t year,int32_t month, int32_t day,int32_t hour, int32_t minute, double sec, double *jul);
#endif