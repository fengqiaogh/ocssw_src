#include <string.h>
#include <libgen.h>
#include "l1stat.h"
#include "l1stat_proto.h"
extern int32 stat_status;
extern char bad_stat_str[320];
#define N_TRNG 20

void ck_trng(char *file)
/*******************************************************************

   ck_trng

   purpose: check the file name to see if it is in one of any bad
            time ranges, fail the file and print out the reason

   Returns type: void - none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            file             I      level-1 file name

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       24-Jan-1998     Original development
      W. Robinson       14-May-1998     Add new time for nav off due to
                                        GPS error
      W. Robinson       2-Jul-1998      Add new time for nav off due to
                                        clock error associated with a hit
                                        on the S/C
      W. Robinson       19-Oct-1998     new time range to auto fail for
                                        GPS clock problem
      W. Robinson       17-May-1999     GPS reset caused some nav probs,
                                        add another fail period
      W. Robinson       28-Mar-2000     GPS outage affects 1 GAC orbit 
                                        by 36 km, so flag it
      W. Robinson       7-Apr-2000      GPS time error added to the 
                                        checking
      W. Robinson       19-Sep-2000     1 sec GPS error for new time range
                                        added
      W. Robinson       13-Nov-2000     GPS error for 1 GAC orbit was added
      W. Robinson       17-Jan-2001     A time tag problem for part of an orbit
                                        on day 2001 011 was set to be flagged
                                        Add 348 and 352 at this time also
      W. Robinson       19-Apr-2001     GPS error for part of 1 GAC orbit
      W. Robinson       21-Jun-2001     add another mandatory fail for a 
                                        1 sec time slip
      W. Robinson       6-Sep-2001      flag a GPS reset causing 1 sec time 
                                        error for just HRPTs after the GAC
      W. Robinson       18-sep-2001     flag GAC and time period in day 256
      W. Robinson       6-nov-2001      report a mnemonic and code as std

 *******************************************************************/ {
    char filefrag[15];
    char stime[N_TRNG][100] = {"S1997292202120", "S1997295010921",
        "S1998133173600", "S1998181170800", "S1998289233900",
        "S1999084031100", "S19991341625", "S20000841825",
        "S20000850602", "S20002591908", "S20003180900", "S20003041710",
        "S20003481922", "S20003520725", "S20010111215", "S20011090325",
        "S20011670840", "S20011940812", "S20012481847", "S20012560303"};
    char etime[N_TRNG][100] = {"S1997292233904", "S1997297023851",
        "S1998133175600", "S1998182161400", "S1998289234551",
        "S19990841407", "S19991341653", "S20000841910",
        "S20000851815", "S20002591923", "S20003180958", "S20003041735",
        "S20003481934", "S20003520737", "S20010111225", "S20011090345",
        "S20011670911", "S20011940852", "S20012481902", "S20012560319"};
    char tcomnt[N_TRNG][100] ={"Time stamp problem: large along-track nav errors",
        "Period of large along-track nav errors",
        "GPS problem: time error causing > 6 km along-track nav errors",
        "Spacecraft Clock problem, causing a large along-track nav error",
        "GPS problem: time error causing 1 sec = 6 km along-track nav errors",
        "GPS problem causing variable pitch, roll errors",
        "GPS problem: time error causing some nav offsets",
        "GPS outage: possible 36 km N-S nav offset in this time range",
        "GPS time off, causing navigation error in whole period",
        "GPS 1 sec error, causing along track error in navigation",
        "GPS time error, causing navigation error in whole period",
        "HRPT requires special L1 processing for gain 2 in this period",
        "GPS time error, causing navigation error in whole period",
        "GPS time error, causing navigation error in whole period",
        "GPS time error, causing navigation error in whole period",
        "GPS time error, causing navigation error in whole period",
        "GPS time slip for this period requiring failure of HRPT",
        "GPS time slip for this period requiring failure of HRPT",
        "GPS reset, causing navigation error for part of the scene",
        "GPS reset, causing navigation error for part of the scene"};
    int i, icode;
    char str[12];

    /*
     *  This is very easy, if the time in the file name is between the
     *  ranges declared in stime and etime, write the comment and fail the file
     */
    strncpy(filefrag, basename(file), 14);
    filefrag[14] = '\0';

    icode = 0;
    for (i = 0; i < N_TRNG; i++) {
        if (strncmp(stime[i], filefrag, 14) <= 0 &&
                strncmp(etime[i], filefrag, 14) >= 0) {
            /* update the general status to product failure */
            stat_status = stat_status | 2;

            /* print the range and the reason for failure */
            /* Also, note that manual inspection is necessary */
            printf("\n\n********* ck_trng: Product is in a mandatory product failure time range\n");
            printf("   (The product requires interactive QC)\n");
            printf("\nFile: '%s'\nis in the mandatory failure time range: '%s' to '%s'\n\n",
                    file, stime[i], etime[i]);
            printf("Reason is:\n  '%s'\n\n", tcomnt[i]);
            icode = 1;
            sprintf(str, "TRNG_CHK ");
            if (strlen(bad_stat_str) <= 300)
                strcat(bad_stat_str, str);
            break;
        }
    }
    /*
     *  print the std mnemonic
     */
    printf("\n\nFor file in time range check...\n");
    printf("    Name  code\n");
    printf("--------  ----\n");
    printf("TRNG_CHK  %4d\n\n", icode);
}
