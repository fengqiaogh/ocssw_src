#define         IS_LEAP_YEAR(y)         ( (!((y)%4) && (y)%100) || !((y)%400) )

int day2mday(int year, int day_of_year, int *month, int *day_of_month)
/*------------------------------------------------------------------------------
   Function: day2mday

   Returns type: 

   Description: Given a year and a day of that year, determine what the month
                and day of that month are.

                This function returns 1 upon successful completion and 0
                otherwise.

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               year            I       all 4 digits (e.g. 1985 or 2001)
      int               day_of_year     I       1 - 365 (or 366 for a leap year)
      unsigned char **  month           O       3 characters plus a terminator
      int *             day_of_month    O       # of day in month (1st day is 1)
   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      Norman Kuring     10-Feb-1993     Original development
      W. Robinson       25-Oct-1996     modify to output integer month

------------------------------------------------------------------------------*/ {

    static char days_per_month[] = {
        31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
    };

    int date;
    int month_index;

    if (IS_LEAP_YEAR(year))
        days_per_month[1] = 29;
    else
        days_per_month[1] = 28;

    for (month_index = 0, date = day_of_year;
            date > days_per_month[month_index];
            date -= days_per_month[month_index++])
        ;

    if (month_index > 11)
        return ( -1);

    *month = month_index + 1;
    *day_of_month = date;

    return ( 0);

}
