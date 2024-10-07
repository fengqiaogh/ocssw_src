#include "viirs_sim_sdr.h"
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#define jul_day_1958 2436205
#define sec_per_day 86400

int rd_geo_init(ctl_struc *ctl, sdr_info_struc *sdr_info,
        in_rec_struc *in_rec)
/*-----------------------------------------------------------------------------
    Routine:   rd_geo_init

    Description:  prepare the geolocation data file for use
       The attributes  and small datasets are read and the file is set so 
       datasets can then be read

   Returns type: int - 0 if good

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        ctl_struc * ctl         I    input controls, use geolocation file
                                     name and ctrate time (optionally)
        sdr_info_struc * sdr_info I/O  general SDR information
        in_rec_struc * in_rec  I/O   input file information

    Modification history:

    W. Robinson, SAIC  20 Nov 2008  Original development
    W. Robinson, SAIC  02 Oct 2009  accomodate the sensor attitude, 
                     position and velocity
    W. Robinson, SAIC  10 Mar 2010  go to scan-oriented processing and use 
                     more scan options (agregated, un-aggregated, margin)

----------------------------------------------------------------------------*/ {
    int month, dom, jul_day, hour, min, sec, h5_stat;
    int16_t lcl_margin[2];
    time_t sec_time;
    struct tm *lcltime, init_tim, *new_tim;
    double st_sec;
    long extra_sec;
    /*
     *  Open the HDF 5 file
     */
    if (h5io_openr(ctl->in_geo_file, 0, &(in_rec->geo_fid)) != 0) {
        printf("%s - could not open HDF 5 file: %s\n",
                __FILE__, ctl->in_geo_file);
        return 1;
    }
    /*
     *  attributes may exist for the other scan formats, if not, it is standard
     */
    if ((h5_stat = h5io_attr_exist(&(in_rec->geo_fid), "scn_fmt"))
            == 0) {
        printf("%s, %d - Geo file indicates detailed scan format information\n",
                __FILE__, __LINE__);
        if (h5io_rd_attr(&(in_rec->geo_fid), "scn_fmt", &(in_rec->scn_fmt))
                != 0) {
            printf("%s, %d - failed to read the scn_fmt attribute\n",
                    __FILE__, __LINE__);
            return 1;
        }
        if ((in_rec->scn_fmt < 0) || (in_rec->scn_fmt > 2)) {
            printf("%s, %d - scn_fmt attribute outside valid range\n",
                    __FILE__, __LINE__);
            return 1;
        }
        if (h5io_rd_attr(&(in_rec->geo_fid), "ndet", &(in_rec->ndet_scan))
                != 0) {
            printf("%s, %d - failed to read the ndet attribute\n",
                    __FILE__, __LINE__);
            printf("        (after reading scn_fmt attribute)\n");
            return 1;
        }
        if (h5io_rd_attr(&(in_rec->geo_fid), "margin", lcl_margin) != 0) {
            printf("%s, %d - failed to read the margin attribute\n",
                    __FILE__, __LINE__);
            printf("        (after reading scn_fmt attribute)\n");
            return 1;
        }
        in_rec->margin[0] = lcl_margin[0];
        in_rec->margin[1] = lcl_margin[1];
    } else {
        in_rec->scn_fmt = 0;
        in_rec->margin[0] = 0;
        in_rec->margin[1] = 0;
        in_rec->ndet_scan = NDET;
    }
    printf("%s, %d - scn_fmt = %d, ndet_scan: %d, margin = %d  %d\n",
            __FILE__, __LINE__, in_rec->scn_fmt, in_rec->ndet_scan, in_rec->margin[0],
            in_rec->margin[1]);
    /*
     *  read the normal attributes 
     */
    if (h5io_rd_attr(&(in_rec->geo_fid), "year", &(sdr_info->year))
            != 0) {
        printf("%s - failed to read the year attribute\n", __FILE__);
        return 1;
    }
    if (h5io_rd_attr(&(in_rec->geo_fid), "day", &(sdr_info->day))
            != 0) {
        printf("%s - failed to read the day attribute\n", __FILE__);
        return 1;
    }
    if (h5io_rd_attr(&(in_rec->geo_fid), "start_sec",
            &(sdr_info->start_sec)) != 0) {
        printf("%s - failed to read the start_sec attribute\n", __FILE__);
        return 1;
    }
    if (h5io_rd_attr(&(in_rec->geo_fid), "npix", &(in_rec->npix))
            != 0) {
        printf("%s - failed to read the npix attribute\n", __FILE__);
        return 1;
    }
    if (h5io_rd_attr(&(in_rec->geo_fid), "nlin", &(in_rec->nlin))
            != 0) {
        printf("%s - failed to read the nlin attribute\n", __FILE__);
        return 1;
    }
    if (h5io_rd_attr(&(in_rec->geo_fid), "nscan", &(in_rec->nscan))
            != 0) {
        printf("%s - failed to read the nscan attribute\n", __FILE__);
        return 1;
    }
    /*
     *  read the small datasets: position, velocity, scan time, and sensor 
     *  attitude The large, pixel-oriented datasets will be transfered a 
     *  line at a  time
     */
    if ((sdr_info->geo_pos = (float *) malloc(in_rec->nscan * 3 *
            sizeof (float))) == NULL) {
        printf("%s:%d - failed to allocate space for geo_pos dataset\n", __FILE__,
                __LINE__);
        return 1;
    }
    if (h5io_grab_ds(&(in_rec->geo_fid), "pos", (void *) sdr_info->geo_pos)
            != 0) {
        printf("%s:%d - failed to grab the geo_pos dataset\n", __FILE__,
                __LINE__);
        return 1;
    }

    if ((sdr_info->geo_vel = (float *) malloc(in_rec->nscan * 3 *
            sizeof (float))) == NULL) {
        printf("%s:%d - failed to allocate space for geo_vel dataset\n", __FILE__,
                __LINE__);
        return 1;
    }
    if (h5io_grab_ds(&(in_rec->geo_fid), "vel", (void *) sdr_info->geo_vel)
            != 0) {
        printf("%s:%d - failed to grab the geo_vel dataset\n", __FILE__,
                __LINE__);
        return 1;
    }

    if ((sdr_info->scan_time =
            (double *) malloc(in_rec->nscan * sizeof (double)))
            == NULL) {
        printf("%s:%d - failed to allocate space for scan_time dataset\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_grab_ds(&(in_rec->geo_fid), "scan_time",
            (void *) sdr_info->scan_time) != 0) {
        printf("%s:%d - failed to grab the scan_time dataset\n", __FILE__,
                __LINE__);
        return 1;
    }

    if ((sdr_info->geo_att = (float *) malloc(in_rec->nscan * 3 *
            sizeof (float))) == NULL) {
        printf("%s:%d - failed to allocate space for geo_att dataset\n", __FILE__,
                __LINE__);
        return 1;
    }
    if (h5io_grab_ds(&(in_rec->geo_fid), "att_arr", (void *) sdr_info->geo_att)
            != 0) {
        printf("%s:%d - failed to grab the geo_att dataset\n", __FILE__,
                __LINE__);
        return 1;
    }

    printf("%s: year, day, start_sec:  %d  %d  %f\n", __FILE__, sdr_info->year,
            sdr_info->day, sdr_info->start_sec);
    printf("%s: npix, nlin, nscan:  %d  %d  %d\n", __FILE__, in_rec->npix,
            in_rec->nlin, in_rec->nscan);
    /*
     *  create the start of the day in a VIIRS used unit - microseconds
     *  past 1/1/1958
     */
    day2mday(sdr_info->year, sdr_info->day, &month, &dom);
    jul_day = jd_c(sdr_info->year, month, dom);
    sdr_info->t58_day = 1.e6 * sec_per_day * ((int64) jul_day - jul_day_1958);
    /*
     *  set the start and end in the 1958 reference too
     */
    sdr_info->st_58_t = sdr_info->t58_day + sdr_info->start_sec * 1.e6;
    sdr_info->en_58_t = sdr_info->st_58_t + SEC_PER_SCAN * in_rec->nscan * 1.e6;
    /*
     *  set up the create dates and times needed for filling the dataset here
     *  create date, time in formats YYYYMMDD, HHMMSS.SSSSSSZ
     *  start end, date, time can only be completely done after trading the geo
     *  time
     */
    if (strcmp(ctl->cre_time, "Unspecified") == 0) {
        time(&sec_time);
        lcltime = localtime(&sec_time);
        sprintf(sdr_info->cre_date, "%4.4d%2.2d%2.2d", lcltime->tm_year + 1900,
                lcltime->tm_mon + 1, lcltime->tm_mday);

        sprintf(sdr_info->cre_time, "%2.2d%2.2d%2.2d.000000Z",
                lcltime->tm_hour, lcltime->tm_min, lcltime->tm_sec);
    } else {
        strncpy(sdr_info->cre_date, ctl->cre_time, 8);
        sdr_info->cre_date[8] = 0;
        strncpy(sdr_info->cre_time, ctl->cre_time + 9, 13);
        sdr_info->cre_time[13] = 'Z';
        sdr_info->cre_time[14] = 0;
    }
    /*
     *  For convenience in creating output file names, create a string 
     *  with the start year, day, hour, min, sec
     */
    hour = (int) (sdr_info->start_sec / 3600.);
    min = (int) (sdr_info->start_sec / 60.) % 60;
    sec = (int) (sdr_info->start_sec) % 60;
    sprintf(sdr_info->ofile_base, "%4.4d%3.3d%2.2d%2.2d%2.2d",
            sdr_info->year, sdr_info->day, hour, min, sec);
    /*
     *  set the start and end times 
     *  to avoid daylight savings complications, set to GMT
     */
    if (putenv("TZ=GMT 0") != 0) {
        printf("%s, %d: GMT time zone switch error\n", __FILE__, __LINE__);
        return 1;
    }
    day2mday(sdr_info->year, sdr_info->day, &month, &dom);
    sprintf(sdr_info->st_date, "%4.4d%2.2d%2.2d", sdr_info->year,
            month, dom);
    st_sec = sdr_info->start_sec;
    init_tim.tm_year = sdr_info->year - 1900;
    init_tim.tm_mday = sdr_info->day;
    init_tim.tm_mon = 0;
    init_tim.tm_hour = 0;
    init_tim.tm_min = 0;
    init_tim.tm_sec = (int) st_sec;
    init_tim.tm_isdst = 0;
    extra_sec = 1.e6 * (st_sec - init_tim.tm_sec);
    sec_time = mktime(&init_tim);
    new_tim = localtime(&sec_time);
    sprintf(sdr_info->st_time, "%2.2d%2.2d%2.2d.%6.6ldZ", new_tim->tm_hour,
            new_tim->tm_min, new_tim->tm_sec, extra_sec);

    init_tim.tm_hour = 0;
    init_tim.tm_min = 0;
    init_tim.tm_sec = (int) (st_sec + (in_rec->nscan - 1) * SEC_PER_SCAN);
    /*  the extra_sec is the fractional seconds that don't get in mktime et al */
    extra_sec = 1.e6 * (st_sec + (in_rec->nscan - 1) * SEC_PER_SCAN
            - init_tim.tm_sec);

    sec_time = mktime(&init_tim);
    new_tim = localtime(&sec_time);

    if (unsetenv("TZ") != 0) {
        printf("%s, %d: GMT time zone unset error\n", __FILE__, __LINE__);
        return 1;
    }

    sprintf(sdr_info->en_date, "%4.4d%2.2d%2.2d", new_tim->tm_year + 1900,
            new_tim->tm_mon + 1, new_tim->tm_mday);
    sprintf(sdr_info->en_time, "%2.2d%2.2d%2.2d.%6.6ldZ", new_tim->tm_hour,
            new_tim->tm_min, new_tim->tm_sec, extra_sec);
    /*
     */
    return 0;
}
