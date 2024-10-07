/* -------------------------------------------------------------- */
/* mkmeta() - creates the meta file associated with each L1A file */
/*                                                                */
/* Synopsis:                                                      */
/*                                                                */
/*   int status =  mkmeta(char *dir, swl0scene  *scene)           */
/*                                                                */
/*   status     error status (-1 = error)                         */
/*   dir        output directory                                  */
/*   scene      fully loaded scene structure                      */
/*                                                                */
/* Description:                                                   */
/*                                                                */
/*   Creates a file in dir with using the L1A naming convention   */
/*   and an additional .meta suffix.  This function assumes the   */
/*   scene structure has been fully loaded, including navigation  */
/*   fields, and the L1A file has been created.  The L1A file     */
/*   must exist for the file_size field of the meta data to be    */
/*   correct.                                                     */
/*                                                                */
/* Modification History:                                          */
/*                                                                */
/*   30 Oct 1997, B. A. Franz, Original development               */
/*                                                                */
/* -------------------------------------------------------------- */

#include <stdio.h>
#include <time.h>
#include <libgen.h>
#include <stdlib.h>
#include <sys/param.h>
#include "swl0_parms.h"
#include "swl0_proto.h"
#include "swl1_hdf.h"
#include <timeutils.h>

int mkmeta(const char *metaFile, const char *l1aFile, swl0scene *scene, swl0ctl *l0ctl) {

    BYTE dataType;
    FILE *fp;
    INT16 year;
    INT16 day;
    FLOAT64 sec;
    char fullname[MAXPATHLEN];
    char *tmppath;
    int navstat = 0;

    if ((fp = fopen(metaFile, "w")) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: error opening %s for writing.\n",
                __FILE__, __LINE__, metaFile);
        return (1);
    }

    unix2yds(scene->stime, &year, &day, &sec);

    fprintf(fp, "%s_version=%s\n", l0ctl->progname, L01VERSION);

    tmppath = strdup(l1aFile);
    fprintf(fp, "l1_filename=%s\n", basename(tmppath));
    free(tmppath);

    tmppath = strdup(l1aFile);
    fprintf(fp, "l1_pathname=%s\n", realpath(dirname(tmppath), fullname));
    free(tmppath);

    fprintf(fp, "file_size=%d\n", filesize(l1aFile));

    if (scene->type == HRPT)
        dataType = 16;
    else
        dataType = scene->mnftype;
    fprintf(fp, "datatype=%s\n", DTypeString(dataType));

    if (scene->lower_left_lat < -90.0)
        navstat = 1;

    fprintf(fp, "format=LEVEL_1\n");
    fprintf(fp, "scan_line_count=%d\n", scene->nscan);
    fprintf(fp, "filled_scan_count=%d\n", 0);
    fprintf(fp, "start_time=%s\n", unix2timeStr(scene->stime));
    fprintf(fp, "stop_time=%s\n", unix2timeStr(scene->etime));
    fprintf(fp, "orbit_number=%d\n", scene->orbnum);

    fprintf(fp, "navstat=%d\n", navstat);
    fprintf(fp, "lower_left_lat=%.2f\n", scene->lower_left_lat);
    fprintf(fp, "lower_left_lon=%.2f\n", scene->lower_left_lon);
    fprintf(fp, "lower_right_lat=%.2f\n", scene->lower_right_lat);
    fprintf(fp, "lower_right_lon=%.2f\n", scene->lower_right_lon);
    fprintf(fp, "upper_left_lat=%.2f\n", scene->upper_left_lat);
    fprintf(fp, "upper_left_lon=%.2f\n", scene->upper_left_lon);
    fprintf(fp, "upper_right_lat=%.2f\n", scene->upper_right_lat);
    fprintf(fp, "upper_right_lon=%.2f\n", scene->upper_right_lon);

    tmppath = strdup(scene->l0file);
    fprintf(fp, "raw_filename=%s\n", basename(scene->l0file));
    free(tmppath);

    tmppath = strdup(scene->l0file);
    fprintf(fp, "raw_pathname=%s\n", realpath(dirname(tmppath), fullname));
    free(tmppath);

    fprintf(fp, "dataday1=%4d%03d\n", year, day);
    fprintf(fp, "dataday2=%d\n", 0);
    fprintf(fp, "day_night_scene=D\n");

    fclose(fp);

    return (1);
}


