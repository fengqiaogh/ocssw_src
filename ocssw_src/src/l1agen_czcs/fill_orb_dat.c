#include "l1czcs.h"
#include <stdlib.h>
#include "time_utl.h"

typedef struct orb_struct_def {
    int n_good; /* # of good points */
    int flag[4]; /* data value flag: -1 - no measurement, 1 measurement */
    double sec[4]; /* seconds attached to record */
    double pos[4][3]; /* position [point #][x, y, or z] */
    double vel[4][3]; /* velocity */
    double pos_err[4]; /* position max error */
} orb_str;

#define s_per_day 24 * 60 * 60
#define ms_per_day s_per_day * 1000

/* function prototypes */
int get_orb_dat(int year, int doy, int st_msec, int en_msec, orb_str *orb);
int int_orb_dat(orb_str orb, int nlines, int *msec, float *orb_vec, float *pos_err);

int rd_smmr_orb(char *sfile, int *nrec, int *syear, int *sday,
        double **orbvec, double **time, float **pos_err);

void asap_int2_(int32_t *nstp, double *tsap, double *asap, double *pos_erri,
        int32_t *ngps, double *gpsec, double *vecs, float *pos_erro);

int fill_orb_dat(l1_data_struc *l1_data, gattr_struc *gattr)
/*******************************************************************

   fill_orb_dat

   purpose: fill the orbit vector with Nimbus 7 data

   Returns type: int - 0 if all went well

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      l1_data_struc *   l1_data          I      level-1 data struct
      gattr_struc *     gattr           I/O     global attribute struct

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 2 Dec 2005      Original development

 *******************************************************************/
 {
    int en_msec;
    orb_str orb;
    /*
     *  Get the orbit data for the period covered by the granule
     */
    en_msec = (gattr->end_day != gattr->start_day) ?
            gattr->end_msec + ms_per_day : gattr->end_msec;
    if (get_orb_dat(gattr->start_year, gattr->start_day, gattr->start_msec,
            en_msec, &orb) != 0)
        return 1;

    /*
     *  interpolate the orbit points to each line of the czcs data
     */
    int_orb_dat(orb, gattr->scan_lines, l1_data->msec, l1_data->orb_vec,
            l1_data->pos_err);
    return 0;
}

int get_orb_dat(int year, int doy, int st_msec, int en_msec, orb_str *orb)
/*******************************************************************

   get_orb_dat

   purpose: extract orbit data needed for a czcs file

   Returns type: int - 0 if no problems

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               year             I      4 digit year of data
      int               doy              I      day of the year
      int               st_msec          I      start millisec of czcs data
      int               en_msec          I      end millisec of czcs data,
                                                relative to start day
      orb_struc *       orb              O      orbit data for at most 4 samples

   Note that as a CZCS granule is at most slightly longer than 2 minutes and
   orbit info is at 1 min intervals, at most 4 samples are required to 
   interpolate a granule

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC  4-Nov-2005     Original development

 *******************************************************************/ {
    char *floc, *fname;
    int i, jd, st_min, en_min, ix_interp, ix_tim;
    int gap_lim = 5;
    int nrec0, nrec1, nrec2, year0, year1, year2, day0, day1, day2;
    double *orb0, *orb1, *orb2, *tmin0, *tmin1, *tmin2;
    float *pos_err0, *pos_err1, *pos_err2;

    /*  defs needed for operation */
    char *day_to_ofile(int day, char *floc);
    /*
     *  initialize some items
     */
    for (i = 0; i < 4; i++)
        orb->flag[i] = -1;
    orb->n_good = 0;
    /*
     *  read in orbit data for the year / day
     */
    jd = yydoy_to_jd(year, doy);
    if ((floc = getenv("OCDATAROOT")) == NULL) {
        printf("%s: environment variable OCDATAROOT undefined\n", __FILE__);
        exit(1);
    }
    strcat(floc, "/czcs/nav/");

    fname = day_to_ofile(jd, floc);
    if (rd_smmr_orb(fname, &nrec1, &year1, &day1, &orb1, &tmin1,
            &pos_err1) != 0) {
        printf(
                "%s: Primary Nimbus 7 orbit file: %s\nwas not found but should exist\n",
                __FILE__, fname);
        return -1;
    }
    /*
     *  if no main day entries, abandon the efort
     */
    if (nrec1 <= 0) {
        orb->n_good = 0;
        return 0;
    }
    /*
     *  get the surrounding minute values from the start, end czcs times
     * and handle minor (> -100 min) times, probably never happen
     */
    st_min = (st_msec + 6000000) / (1000 * 60) - 100;
    en_min = ((en_msec + 59999) / 1000) / 60;
    ix_interp = 0;
    /*
     *  get data from previous day if needed
     */
    if (*(tmin1) > st_min) {
        /*
         *  read previous day
         */
        fname = day_to_ofile((jd - 1), floc);
        if (rd_smmr_orb(fname, &nrec0, &year0, &day0, &orb0, &tmin0,
                &pos_err0) == 0) {
            if (nrec0 > 0) {
                if ((st_min - (*(tmin0 + nrec0 - 1) - 1440)) > gap_lim) {
                    /* we can't add this in as it is beyond interpolation limit */
                    orb->flag[ ix_interp ] = -1;
                } else {
                    orb->flag[ ix_interp ] = 1;
                    orb->sec[ ix_interp ] = *(tmin0 + nrec0 - 1) - 1440;
                    orb->pos_err[ ix_interp ] = *(pos_err0 + nrec0 - 1);
                    for (i = 0; i < 3; i++) {
                        orb->pos[ix_interp][i] = *(orb0 + (nrec0 - 1) * 6 + i);
                        orb->vel[ix_interp][i] = *(orb0 + (nrec0 - 1) * 6 + 3 + i);
                    }
                    ix_interp++;
                }
                /*
                 *  remove storage for the prev day
                 */
                free(orb0);
                free(tmin0);
                free(pos_err0);
            }
        }
    }
    /*
     *  get back to the current day
     *  Note that we will start collecting samples starting at start min - 
     *  gap_lim, but discard them if closer times are in file.
     */
    ix_tim = 0;
    while ((*(tmin1 + ix_tim) < (en_min + gap_lim)) &&
            (ix_tim < nrec1) && (ix_interp < 4)) {
        if (*(tmin1 + ix_tim) >= st_min - gap_lim) {
            /* discard a time found before start if a closer one exists */
            if ((ix_interp > 0) && *(tmin1 + ix_tim) <= st_min)
                ix_interp--;
            orb->flag[ ix_interp ] = 1;
            orb->sec[ ix_interp ] = *(tmin1 + ix_tim);
            orb->pos_err[ ix_interp ] = *(pos_err1 + ix_tim);
            for (i = 0; i < 3; i++) {
                orb->pos[ix_interp][i] = *(orb1 + ix_tim * 6 + i);
                orb->vel[ix_interp][i] = *(orb1 + ix_tim * 6 + 3 + i);
            }
            ix_interp++;
        }
        ix_tim++;
    }
    /*
     *  free storage
     */
    free(orb1);
    free(tmin1);
    free(pos_err1);
    /*
     *  possibly, orbit information from the following day may also be needed
     */
    if (ix_interp < 4) {
        if (*(tmin1 + nrec1 - 1) <= en_min) {
            fname = day_to_ofile((jd + 1), floc);
            if (rd_smmr_orb(fname, &nrec2, &year2, &day2, &orb2, &tmin2,
                    &pos_err2) == 0) {
                if (nrec2 > 0) {
                    if ((*tmin2 + 1440 - en_min) <= gap_lim) {
                        ix_tim = 0;
                        while (((*(tmin2 + ix_tim) + 1440) < (en_min + gap_lim)) &&
                                (ix_tim < nrec2) && (ix_interp < 4)) {
                            if ((*(tmin2 + ix_tim) + 1440) >= st_min) {
                                orb->flag[ ix_interp ] = 1;
                                orb->sec[ ix_interp ] = *(tmin2 + ix_tim) + 1440;
                                orb->pos_err[ ix_interp ] = *(pos_err2 + ix_tim);
                                for (i = 0; i < 3; i++) {
                                    orb->pos[ix_interp][i] = *(orb2 + ix_tim * 6 + i);
                                    orb->vel[ix_interp][i] = *(orb2 + ix_tim * 6 + 3 + i);
                                }
                                ix_interp++;
                            }
                            ix_tim++;
                        }
                    }
                    /*
                     *  free storage
                     */
                    free(orb2);
                    free(tmin2);
                    free(pos_err2);
                }
            }
        }
    }
    orb->n_good = ix_interp;

    /*
     *  convert the orb->sec from minutes to seconds
     */
    for (i = 0; i < 4; i++) {
        orb->sec[i] = orb->sec[i] * 60.;
    }

    /*
     *  and end
     */
    return 0;
}

int int_orb_dat(orb_str orb, int nlines, int *msec, float *orb_vec,
        float *pos_err)
/*******************************************************************

   int_orb_dat

   purpose: interpolate the Nimbus 7 orbit information to each line of CZCS

   Returns type: int - 0 if all is OK

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      orb_str           orb              I      structure containing the 
                                                nimbus 7 orbit data available
                                                at the CZCS data time range
      int               nlines           I      # lines of CZCS data
      int *             msec             I      time of each CZCS line in
                                                msec
      float *           orb_vec          O      orbit data for each line to fill
      float *           pos_err          O      orbit position error estimate:
                                                < 0, no orbit data available
                                                0, no error or original smmr
                                                orbit data
                                                > 0, orbit modelled data with
                                                this error estimate in m.
   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC  2-Dec-2005     Original development

 *******************************************************************/ {
    int start_sec, i, j;
    double *sec, asap[24], *lcl_orb;
    /*
     *  at this time, use the asap_int.f routine, slightly modified to interpolate
     */
    /*
     *  we must adapt some inputs for asap
     *  msec to seconds
     *  and orbit position, velocity to just an orbit vector
     */
    start_sec = *msec / 1000 - 1;

    if ((sec = (double *) malloc(nlines * sizeof ( double))) == NULL) {
        printf("%s: malloc failed for seconds array\n", __FILE__);
        return 1;
    }
    if ((lcl_orb = (double *) malloc(nlines * 6 * sizeof ( double))) == NULL) {
        printf("%s: malloc failed for lcl_orb array\n", __FILE__);
        return 1;
    }
    for (i = 0; i < nlines; i++) {
        *(sec + i) = (double) *(msec + i) / 1000.;
        if (*(sec + i) < start_sec)
            *(sec + i) = *(sec + i) + ms_per_day;
    }

    /*
     *  with 1 or less good points, the entire czcs file is designated to 
     *  not have orbit data
     */
    if (orb.n_good <= 1) {
        for (i = 0; i < nlines; i++) {
            *(pos_err + i) = -1.;
        }
    } else {
        /*
         *  with good orbit to interpolate, proceed
         *  get the orbit info set up correctly for asap routine
         */
        for (i = 0; i < orb.n_good; i++) {
            for (j = 0; j < 3; j++) {
                *(asap + j + i * 6) = orb.pos[i][j];
                *(asap + j + 3 + i * 6) = orb.vel[i][j];
            }
        }
        /*
         *  call the interpolation routine
         */
        asap_int2_(&(orb.n_good), orb.sec, asap, orb.pos_err, &nlines,
                sec, lcl_orb, pos_err);
        /*
         *  place position part of orbit in orb_vec
         */
        for (i = 0; i < nlines; i++) {
            for (j = 0; j < 3; j++) {
                *(orb_vec + j + 3 * i) = *(lcl_orb + j + 6 * i);
            }
        }
    }
    /*
     *  free the seconds and local orbit
     */
    free(lcl_orb);
    free(sec);
    /*
     *  all done
     */
    return 0;
}
