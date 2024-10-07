#include <stdio.h>

int time_str(short year, short jday, int millisec, char *time)
/*******************************************************************

   time_str

   purpose: Make the time string from year, day and millisec

   Returns type: int 0 if OK, -1 if error writing to the time string

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      short             year             I      year
      short             jday             I      day of year
      int               millisec         I      time of day in milliseconds
      char *            time             O      string with time in form:
                                                YYYYDDDHHMMSSFFF Note that
                                                time should be created with 
                                                at least 17 characters

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       13-Apr-2004     Original development

 *******************************************************************/ {
    int sec, msec, min, hr;

    /*
     *  Break apart the millisecs to hours, mins, secs and fraction
     */
    sec = millisec / 1000;
    msec = millisec - sec * 1000;
    min = sec / 60;
    sec = sec - min * 60;
    hr = min / 60;
    min = min - hr * 60;
    /*
     * fill the string
     */
    if ((sprintf(time, "%04d%03d%02d%02d%02d%03d", year, jday, hr,
            min, sec, msec)) != 16)
        return -1;
    return 0;
}
