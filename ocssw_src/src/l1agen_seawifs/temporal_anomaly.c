/* ============================================================== */
/* Module:  timetable.c                                           */
/* Purpose: generate and search temporal anomaly table            */
/* Author:  B.A. Franz, General Scences Corp., 9/2001             */
/* ============================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>

#include "swl0_proto.h"
#include <timeutils.h>

#define MAXTTABLE 10000

static timeTab timeTable[MAXTTABLE];
static int32_t nTab = 0;

int compTimeTable(timeTab *rec1, timeTab *rec2) {
    int retval = strcmp(rec1->time1, rec2->time1);

    if (retval == 0)
        retval = (int) rec1->type - rec2->type;

    return (retval);
}

INT32 read_timetable() {
    FILE *fp = NULL;
    char *filename;
    int status = 1;

    /*									*/
    /* Open file for reading						*/
    /*									*/
    if ((filename = getenv("TEMPORAL_ANOMALY")) == NULL) {
        printf("-W- %s: TEMPORAL_ANOMALY env variable undefined.\n", __FILE__);
        return (FATAL_ERROR);
    }

    if ((fp = fopen(filename, "r")) == NULL) {
        fprintf(stderr,
                "-W- %s line %d: unable to open %s for reading\n",
                __FILE__, __LINE__, filename);
        return (FATAL_ERROR);
    }

    while (status != EOF && nTab < MAXTTABLE) {

        status = fscanf(fp, "%16s %6d %16s %6d %6d\n",
                timeTable[nTab].time1,
                &timeTable[nTab].type,
                timeTable[nTab].time2,
                &timeTable[nTab].action,
                &timeTable[nTab].delmsec);

        if (status != EOF) nTab++;

    }

    fclose(fp);
    return (0);
}

INT32 update_timetable(FILE *fp, timeTab timeTable) {
    if (fp != NULL)
        fprintf(fp, "%16s %6d %16s %6d %6d\n",
            timeTable.time1,
            timeTable.type,
            timeTable.time2,
            timeTable.action,
            timeTable.delmsec);

    return (0);
}

INT32 write_timetable(char *filename) {
    FILE *fp = NULL;
    INT32 itab = 0;

    /*									*/
    /* Open file for writing						*/
    /*									*/
    if ((fp = fopen(filename, "w")) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: unable to open %s for writing\n",
                __FILE__, __LINE__, filename);
        return (FATAL_ERROR);
    }


    for (itab = 0; itab < nTab; itab++)
        update_timetable(fp, timeTable[itab]);


    fclose(fp);
    return (0);
}

INT32 locate_temporal_anomalies(swl0indx *indx, char *updfile) {
    INT32 srec = indx->srec;
    INT32 erec = indx->erec;
    INT32 irec = srec - 1;
    INT32 i1, i2;
    INT32 uniq[MAXTTABLE];
    INT32 delmsec;
    FLOAT64 delsec;
    FILE *fp = NULL;

    /* Load any previously defined timing anomaly records */
    read_timetable();

    /* If no update file specified, skip search for new anomalies */
    if (updfile == NULL || indx->type != GAC)
        return (0);

    /* Open update file for writing */
    if ((fp = fopen(updfile, "w")) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: unable to open %s for writing\n",
                __FILE__, __LINE__, updfile);
        exit(FATAL_ERROR);
    }

    /* Search for the common 1-sec by 30-sec time shift */
    while (irec < erec - 10) {

        irec++;

        if (indx->rec[irec].mnftype != GACTYPE)
            continue;

        /* Check for the 1-sec shift with 9-frame interval */
        /* ----------------------------------------------- */

        i1 = irec;
        for (i2 = irec + 8; i2 <= irec + 10; i2++) {

            if (indx->rec[i1].tShfError && indx->rec[i2].tShfError &&
                    fabs(fabs(indx->rec[i1].timeShift) - 1.0) < CLOCKVAR &&
                    fabs(indx->rec[i1].timeShift + indx->rec[i2].timeShift) < 2 * CLOCKVAR) {

                delsec = indx->rec[i1].timeShift;
                delmsec = (int32_t) (delsec * 1000);

                printf("Correctable timeshift located at frame %d of %d msecs\n",
                        i1, delmsec);

                /* set GAC frames to correct time shift */
                strcpy(timeTable[nTab].time1, unix2ydhmsf(indx->rec[i1].time, 'G'));
                strcpy(timeTable[nTab].time2, unix2ydhmsf(indx->rec[i2].time, 'G'));
                timeTable[nTab].delmsec = delmsec;
                timeTable[nTab].type = GACTYPE;
                timeTable[nTab].action = 1; /* fix time   */
                nTab++;
                if (nTab >= MAXTTABLE) {
                    fprintf(stderr,
                            "-E- %s line %d: temporal anomaly table size exceeded: %d\n",
                            __FILE__, __LINE__, nTab);
                    exit(FATAL_ERROR);
                }
                update_timetable(fp, timeTable[nTab - 1]);


                /* Set LAC frames for nav failure at leading boundary of time shift */
                strcpy(timeTable[nTab].time1, unix2ydhmsf(indx->rec[i1].time - delsec - 19.5 * DTLAC, 'G'));
                strcpy(timeTable[nTab].time2, unix2ydhmsf(indx->rec[i1].time, 'G'));
                timeTable[nTab].delmsec = (int32_t) 0;
                timeTable[nTab].type = LACTYPE;
                timeTable[nTab].action = 2; /* nav fail */
                nTab++;
                if (nTab >= MAXTTABLE) {
                    fprintf(stderr,
                            "-E- %s line %d: temporal anomaly table size exceeded: %d\n",
                            __FILE__, __LINE__, nTab);
                    exit(FATAL_ERROR);
                }
                update_timetable(fp, timeTable[nTab - 1]);

                /* Set LAC frames for time correction for n-1 GAC frames of time shift */
                /* Note: we can't predict where in the last frame the shift occurred.  */
                strcpy(timeTable[nTab].time1, unix2ydhmsf(indx->rec[i1].time, 'G'));
                strcpy(timeTable[nTab].time2, unix2ydhmsf(indx->rec[i2].time - delsec - 19.5 * DTLAC, 'G'));
                timeTable[nTab].delmsec = delmsec;
                timeTable[nTab].type = LACTYPE;
                timeTable[nTab].action = 1; /* fix time */
                nTab++;
                if (nTab >= MAXTTABLE) {
                    fprintf(stderr,
                            "-E- %s line %d: temporal anomaly table size exceeded: %d\n",
                            __FILE__, __LINE__, nTab);
                    exit(FATAL_ERROR);
                }
                update_timetable(fp, timeTable[nTab - 1]);

                /* Set LAC frames for nav failure within last GAC frame of time shift */
                strcpy(timeTable[nTab].time1, unix2ydhmsf(indx->rec[i2].time - delsec - 19.5 * DTLAC, 'G'));
                strcpy(timeTable[nTab].time2, unix2ydhmsf(indx->rec[i2].time, 'G'));
                timeTable[nTab].delmsec = (int32_t) 0;
                timeTable[nTab].type = LACTYPE;
                timeTable[nTab].action = 2; /* nav fail */
                nTab++;
                if (nTab >= MAXTTABLE) {
                    fprintf(stderr,
                            "-E- %s line %d: temporal anomaly table size exceeded: %d\n",
                            __FILE__, __LINE__, nTab);
                    exit(FATAL_ERROR);
                }
                update_timetable(fp, timeTable[nTab - 1]);

                /* Advance irec beyond shift and exit shift-match search */
                irec = i2 + 1;
                break;

            }
        } /* end for over shift pairs */

    } /* end while over irec */

    /* Ensure the time table array is sorted and unique */
    if (nTab > 0) {

        qsort(timeTable, nTab, sizeof (timeTab),
                (int (*)(const void *, const void *)) compTimeTable);

        i1 = 0;
        i2 = 0;
        uniq[i2] = i1;
        for (i1 = 1; i1 < nTab; i1++)
            if (strcmp(timeTable[i1].time1, timeTable[i1 - 1].time1) != 0 ||
                    timeTable[i1].type != timeTable[i1 - 1].type) {
                i2++;
                uniq[i2] = i1;
            }

        nTab = i2 + 1;
        for (i1 = 1; i1 < nTab; i1++)
            if (uniq[i1] != i1)
                memcpy(&timeTable[i1], &timeTable[uniq[i1]], sizeof (timeTab));
    }

    /*       
    write_timetable(timefile);
     */

    if (fp != NULL) fclose(fp);

    return (0);
}

timeTab *temporal_anomaly(INT16 mnftype, INT16 year, INT16 day, INT32 msec) {
    timeTab *trec = NULL;
    int32_t itab = 0;
    FLOAT64 sec = ((double) msec + 0.5) / 1000.0;
    char ztime[17];

    strcpy(ztime, unix2ydhmsf(yds2unix(year, day, sec), 'G'));

    for (itab = 0; itab < nTab; itab++) {

        if (timeTable[itab].type != mnftype)
            continue;

        if (strcmp(ztime, timeTable[itab].time2) < 0) {
            if (strcmp(ztime, timeTable[itab].time1) >= 0) {
                trec = &timeTable[itab];
            }
            break;
        }
    }

    return (trec);
}
