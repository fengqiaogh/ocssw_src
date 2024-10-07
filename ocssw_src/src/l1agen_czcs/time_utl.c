#include <sys/types.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "time_utl.h"

/*
  NOTE that the $SWFAPP/QCdisp/time_utl.c has the fortran code for 
    jddate
    jdate
    jd
  from which the algorithms for the c code was taken
     W. Robinson, SAIC, 23 Sep 2005

 */
int jd_to_yydoy(int jd, int *yy, int *doy)
/*******************************************************************

   jd_to_yydoy

   purpose: convert julian day to year and day of year

   Returns type: int - 0

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               jd               I      julian day (ref to 4xxx BC)
      int *             yy               O      4 digit year
      int *             doy              O      day of the year

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 23-Sep-2005     development from ftn code

 *******************************************************************/ {
    int day1900;
    /*
     *  Compute days since January 0, 1900 and years since 1900
     */
    day1900 = jd - 2415020;

    *yy = 4 * day1900 / 1461;
    /*
     *  Compute day-of-year and Add first two digits of year
     */
    *doy = day1900 - 1461 * (*yy - 1) / 4 - 365;
    *yy = *yy + 1900;
    /*
     *  and end
     */
    return 0;
}

int jd_to_ydoy(int jd, int *y, int *doy)
/*******************************************************************

   jd_to_yydoy

   purpose: convert julian day to year and day of year

   Returns type: int - 0

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               jd               I      julian day (ref to 4xxx BC)
      int *             y                O      2 digit year
      int *             doy              O      day of the year

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 23-Sep-2005     development from ftn code

 *******************************************************************/ {
    int yy;
    jd_to_yydoy(jd, &yy, doy);
    *y = yy - 1900;
    return 0;
}

int yymd_to_jd(int yy, int m, int d)
/*******************************************************************

   yymd_to_jd

   purpose: convert year, month, day to julian day

   Returns type: int - julian day (ref 4xxx BC)

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               yy               I      4 digit year
      int               m                I      month 1 - 12
      int               d                I      day of month 1 - 31

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 23-Sep-2005     development from ftn version

 *******************************************************************/ {
    return ( 367 * yy -
            7 * (yy + (m + 9) / 12) / 4 +
            275 * m / 9 +
            d + 1721014);
}

int yydoy_to_md70(int yy, int doy, int *m, int *d)
/*******************************************************************
  
   yydoy_to_md70
  
   purpose: convert year and day of year into month and day
    NOTE *** only good for dates post 1970
    see yydoy_to_md for more conventional but universal conversion

   Returns type: int - 0 no meaning

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               yy               I      4 digit year
      int               doy              I      day of year
      int *             m                O      month
      int *             d                O      day of month

 *******************************************************************/ {
    struct tm mytm, *out_tm;
    time_t mytime;

    /*
     *  place the day of year in as day of month and set month to 0,
     *  mktime and localtime will set the month and day of month properly
     */
    mytm.tm_year = yy - 1900; /* Years after 1900 */
    mytm.tm_mon = 0;
    mytm.tm_mday = doy;
    mytm.tm_sec = 0;
    mytm.tm_min = 0;
    mytm.tm_hour = 0;
    mytm.tm_isdst = 0;

    mytime = mktime(&mytm);
    out_tm = localtime(&mytime);

    *m = out_tm->tm_mon + 1;
    *d = out_tm->tm_mday;
    return 0;
}

int yydoy_to_jd(int yy, int doy)
/*******************************************************************

   yydoy_to_jd

   purpose: convert 4 digit year, and day of year to julian day

   Returns type: int - julian day

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               yy               I      2 digit year in 1900 - 1999
      int               doy              I      day of year

 *******************************************************************/ {
    int m, d;
    yydoy_to_md(yy, doy, &m, &d);
    return yymd_to_jd(yy, m, d);
}

int ydoy_to_jd(int y, int doy)
/*******************************************************************

   ydoy_to_jd

   purpose: convert 2 digit year, and day of year to julian day

   Returns type: int - julian day

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               y                I      4 digit year
      int               doy              I      day of year

 *******************************************************************/ {
    int yy;

    yy = y + 1900;
    return yydoy_to_jd(yy, doy);
}

int yydoy_to_md(int yy, int doy, int *m, int *d)
/*******************************************************************

   yydoy_to_md

   purpose: convert year and day of year into month and day

   Returns type: int - 0 no meaning

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               yy               I      4 digit year
      int               doy              I      day of year
      int *             m                O      month
      int *             d                O      day of month

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC  7-Nov-2005     development from ftn gregor.f

 *******************************************************************/ {
    int daymon[2][12] = {
        { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334},
        { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335}
    };
    int ly;

    ly = leap(yy);
    *m = 11;
    while (doy <= daymon[ly][*m])
        *m -= 1;
    *d = doy - daymon[ly][*m];
    *m += 1;
    return 0;
}

int leap(int yy)
/*******************************************************************

   leap

   purpose: find if the year is a leap year

   Returns type: int - 0 if not leap year, 1 if leap year

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               yy               I      4 digit year

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC  7-Nov-2005     development from leap.f

 *******************************************************************/ {
    return ( (yy % 400 == 0) ||
            ((yy % 4 == 0) && (yy % 100 != 0))) ? 1 : 0;
}
