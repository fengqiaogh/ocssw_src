/* -------------------------------------------------------------- */
/* getorbnum() - compute seastar orbit number for given time      */
/*                                                                */
/* Synopsis:                                                      */
/*                                                                */
/*   orbnum = getorbnum( time );                                  */
/*                                                                */
/*   INT32    orbnum  - orbit number                              */
/*   FLOAT64  time    - time in seconds since 1/1/70              */
/*                                                                */
/* Description:                                                   */
/*                                                                */
/*   Given a table of epoch, orbital period, and orbit number,    */
/*   where the epoch is time in seconds at the start of a new     */
/*   orbit number, just divide time since last epoch by orbital   */
/*   period to get orbit number.  This function is not intended   */
/*   to be highly accurate as to the exact time of orbit number   */
/*   change, but if the last table entry is within a few months   */
/*   of the input time, and the orbit does not suddenly degrade,  */
/*   the transition time should be accurate to the minute.        */
/*                                                                */
/*   New table entries should be added every few months.  A       */
/*   utility to generate them is called SWorbnum.                 */
/*                                                                */
/* Modification Hisory:                                           */
/*                                                                */
/*   03 Nov 1997, B. A. Franz, Initial Devalopment                */
/*                                                                */
/* -------------------------------------------------------------- */

#include "swl0_proto.h"
#include <string.h>
#include <stdio.h>

typedef struct orbtab_struct {
    INT16 year;
    INT16 day;
    FLOAT32 sec;
    FLOAT64 time;
    FLOAT32 period;
    INT32 orbnum;
} orbtab_str;

INT32 getorbnum(FLOAT64 usec) {
    static int firstCall = 1;
    static orbtab_str orbtab[MAXORBTAB];
    static INT32 ntab = 0;

    FLOAT64 delsec;
    INT32 orbnum;
    INT32 itab;

    if (firstCall) {

        FILE *fp = NULL;
        char *tmp_str;
        char filename[1024] = "";
        char line[80];

        INT16 year;
        INT16 day;
        FLOAT32 sec;
        FLOAT64 time;
        FLOAT32 period;

        firstCall = 0;

        if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
            printf("OCDATAROOT environment variable is not defined.\n");
            exit(1);
        }
        strcpy(filename, tmp_str);
        strcat(filename, "/seawifs/nav/orbit_table.dat");

        printf("\nLoading orbit number table from %s\n", filename);

        if ((fp = fopen(filename, "r")) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: unable to open %s for reading\n",
                    __FILE__, __LINE__, filename);
            exit(1);
        }

        itab = 0;
        while (fgets(line, 80, fp)) {
            sscanf(line, "%hd %hd %f %lf %f %d",
                    &year, &day, &sec, &time, &period, &orbnum);

            if (itab > (MAXORBTAB - 1)) {
                fprintf(stderr,
                        "-E- %s line %d: orbit table exceeded allocation\n",
                        __FILE__, __LINE__);
                exit(1);
            }

            if (year < 1997 || year > 2020)
                continue;

            orbtab[itab].year = year;
            orbtab[itab].day = day;
            orbtab[itab].sec = sec;
            orbtab[itab].time = time;
            orbtab[itab].period = period;
            orbtab[itab].orbnum = orbnum;
            itab++;
        }

        ntab = itab;
        fclose(fp);

        printf("%d orbit records loaded\n", ntab);
    }


    for (itab = 1; itab < ntab; itab++) {
        if (usec < orbtab[itab].time) {
            itab--;
            break;
        }
    }

    if (itab == ntab) itab--;

    printf("Using orbit number extrapolated from %d %d\n\n",
            orbtab[itab].year, orbtab[itab].day);



    delsec = usec - orbtab[itab].time;
    orbnum = orbtab[itab].orbnum + (int32_t) (delsec / orbtab[itab].period);

    if (delsec < 0.0)
        orbnum--;

    if (delsec / 86400. > 180)
        printf("-W- Orbit table is out of date, closest entry is %lf days old.\n",
            delsec / 86400.);

    return (orbnum);
}


#ifdef GETORBIT

void main(int argc, char *argv[]) {
    INT16 year;
    INT16 day;
    FLOAT64 sec;
    FLOAT64 usec;
    INT32 orbnum;
    INT32 orbnum2;

    year = atoi(argv[1]);
    day = atoi(argv[2]);

    orbnum2 = 0;

    for (sec = 0.0; sec < 86400.0; sec += 0.1) {
        usec = yds2unix(year, day, sec);
        orbnum = getorbnum(usec);
        if (orbnum != orbnum2) {
            orbnum2 = orbnum;
            printf("%d %d %lf %d\n", year, day, sec, orbnum);
        }
    }

}
#endif

