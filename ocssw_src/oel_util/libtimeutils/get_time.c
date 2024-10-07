/*----------------------------------------------------------------------------
  Function: get_time 
  
  Returns: None
  
  Description:
  The function get_time accesses system current time (GMT) and passes the  
  concatenated string of digits for year, day-of-year, hours, minutes,
  seconds, and fraction of seconds in the format of YYYYDDDHHMMSSFFF
  to the calling routine.
  
  Arguments: (in calling order)
  Type       Name             I/O     Description
  ----       ----             ---     -----------
  char *     pr_time           O      Processing time 
  
  Notes:
  
  Modification history:
  Programmer     Organization      Date      Description of change
  --------------   ------------    --------    ---------------------
  Lakshmi Kumar    Hughes STX      07/14/94    Original development
  ---------------------------------------------------------------------------*/
#include <time.h>
#include <stdio.h>
#include <timeutils.h>

void get_time(char *pr_time) {
    struct tm *tp;
    time_t ptime;

    ptime = time(&ptime);
    tp = gmtime(&ptime);

    sprintf(pr_time, "%4d%03d%02d%02d%02d%03d",
            tp->tm_year + 1900,
            tp->tm_yday + 1,
            tp->tm_hour,
            tp->tm_min,
            tp->tm_sec,
            0);
}








