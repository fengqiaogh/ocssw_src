#include "l1io.h"

int anc_daymon(char *fname, int *jday, char *monstr)
/*******************************************************************

   anc_daymon

   purpose: find the ancillary file start day and convert that to the month

   Returns type: int - return status: 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            fname            I      name of file to open
      int *             jday             O      julian start day
      char *            monstr           O      month the start is in

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       17-Jun-1997     Original development

 *******************************************************************/
 {
    l1info_struct hdf_info;
    int month, day, year;
    char *mon_list[] = {"January", "February", "March", "April", "May",
        "June", "July", "August", "September", "October",
        "November", "December"};
    int32 count = 1, n_type = DFNT_INT16;
    int16 i16;

    /*
     *  open the file and read the julian day
     */
    if (open_hdf(fname, &hdf_info) != 0) {
        printf("anc_daymon: open_hdf error\n");
        return -1;
    }
    if (read_g_attr(hdf_info, "Start Day", &n_type, &count, (void *) &i16) != 0) {
        printf("anc_daymon: read_g_attr (day) error\n");
        return -1;
    }
    *jday = i16;
    if (read_g_attr(hdf_info, "Start Year", &n_type, &count, (void *) &i16) != 0) {
        printf("anc_daymon: read_g_attr (year) error\n");
        return -1;
    }
    year = i16;
    if (*jday < 1 || *jday > 366) {
        printf("anc_daymon: ancillary file julian day out of bounds: %d\n",
                *jday);
        if (*jday < 1) *jday = 1;
        if (*jday > 366) *jday = 365;
    }

    l1io_close(hdf_info);

    /*
     *  get the month now
     */
    if (day2mday(year, *jday, &month, &day) != 0) {
        printf("anc_daymon: trouble with day2mday\n");
        return -1;
    }

    strcpy(monstr, mon_list[month - 1]);
    return 0;
}
