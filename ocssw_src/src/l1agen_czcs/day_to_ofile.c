#include <stdio.h>
#include <stdlib.h>
#include "time_utl.h"

char *day_to_ofile(int day, char *floc)
/*******************************************************************

   day_to_ofile

   purpose: find the orbit file name for the julian day

   Returns type: char * - file name

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               day              I      julian day
      char *            floc             I      file partial path
                                                  
   NOTE that the ful file name is gotten using the year and day out of day
      and the partial file path this way:
   name = <floc>YYYY/SMMR_ORB_YYDDD.dat
     <floc> is partial path
     YYYY, YY are 4 and 2 digit years
     DDD is the julian day

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 22-Sep-2005     Original development

 *******************************************************************/ {
    static char fname[500], *cptr;
    int year, doy;
    /*
     *  get the year, day of year and create the file name
     */
    jd_to_ydoy(day, &year, &doy);
    sprintf(fname, "%s19%2d/SMMR_ORB_%2d%03d.dat",
            floc, year, year, doy);
    cptr = fname;
    return cptr;
}
