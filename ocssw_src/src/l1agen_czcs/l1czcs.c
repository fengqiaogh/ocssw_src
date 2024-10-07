/*-----------------------------------------------------------------------------
    Program:   l1czcs

    Description:
               Create CZCS L1-A HDF file from a CZCS CRT full-resolution
               L1 files. 

    Arguments: 
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       argc          I    count of command line args - 3 always
        char *[]  argv          I    command line arguments:
                                     [1] = input file name
                                     [2] = output path

    Modification history:

    W. Robinson, SAIC  5 Mar 2004   derived from seadas l1aczcs.  Improved
                                    to have additional parameters for improved
                                    processing and a switch to a c routine for 
                                    reading the CRTT dataset

----------------------------------------------------------------------------*/

#include <unistd.h>
#include <libgen.h>
#include <string.h>
#include <stdio.h>
#include <sys/utsname.h>

#include "l1czcs.h"

#include <timeutils.h>
#include <passthebuck.h>
#include <GetStationInfo.h>

#define PROGRAM "l1czcs"
#define VERSION "1.1"
#define RET_GOOD 0
#define RET_ERR 1
#define RET_WARN 2

#define MAXLINE 3000

int main(int argc, char *argv[]) {
    /*
     *  for the library functions, define here
     */

    char ifile[132], ofile[256], opath[132];
    char *fbase;

    gattr_struc gattr;
    l1_data_struc l1_data;

    char seadas_vs[64], prog_version[64];

    int i, ilen;

    memset(ifile, 0, sizeof (ifile));

    if (argc == 3) {
        strcpy(ifile, argv[1]);
        strcpy(opath, argv[2]);
    } else {
        fbase = basename(argv[0]);
        printf("%s %s (%s %s)\n",
                fbase, VERSION, __DATE__, __TIME__);
        printf("Usage: %s ifile opath\n\n", fbase);
        printf("Return 0 if OK, 1 if file cannot be created, 2 if no navigation\n");
        return RET_ERR;
    }

    if (access(ifile, F_OK) || access(ifile, R_OK)) {
        printf("%s - Input file '%s' does not exist or cannot read\n",
                argv[0], ifile);
        return RET_ERR;
    }

    /*
     *  Read in all the data from the CRTT file
     */
    if (read_crtt(ifile, &gattr, &l1_data) != 0) {
        printf("%s - Error in reading or processing the data\n", argv[0]);
        return RET_ERR;
    }
    /*
     *  clean up the data (eliminate bad msec lines)
     */
    if (cz_clean(&gattr, &l1_data) != 0) {
        printf("%s - Dataset contains no good data lines\n", argv[0]);
        return RET_ERR;
    }

    StationInfo stationInfo;
    if (GetStationInfo(NULL, &stationInfo) == LIFE_IS_GOOD) {
        strcpy(gattr.datacenter, stationInfo.data_center);
        strcpy(gattr.stn_name, stationInfo.station_name);
        gattr.stn_lat = stationInfo.station_latitude;
        gattr.stn_lon = stationInfo.station_longitude;
    }

    sprintf(prog_version, "%s %s", PROGRAM, VERSION);
    struct utsname osname;
    uname(&osname);
    sprintf(gattr.soft_id, "%s, %s, %s %s", seadas_vs, prog_version, osname.sysname,
            osname.release);

    strcpy(gattr.datatype, "LAC");
    gattr.lac_pixl_start_no = 1;
    gattr.lac_pixl_subsample = 1;

    get_time((char*) &gattr.process_time);

    strcpy(gattr.proc_ctl, argv[0]);
    for (i = 1; i < argc; i++) {
        strcat(gattr.proc_ctl, " ");
        strcat(gattr.proc_ctl, argv[i]);
    }

    strcpy(gattr.input_files, basename(ifile));

    ilen = strlen(opath);
    if (opath[ ilen - 1 ] == '/')
        sprintf(ofile, "%sC%13.13s.L1A_LAC", opath, gattr.start_time);
    else
        sprintf(ofile, "%s/C%13.13s.L1A_LAC", opath, gattr.start_time);

    if (!access(ofile, F_OK)) {
        printf("%s - Output file '%s' already existed\n", argv[0], ofile);
        return RET_ERR;
    }

    printf("\n Writing data to output file...\n");
    if (czcs_l1_write(ofile, l1_data, gattr)
            == 0)
        printf(" output file created\n");
    else {
        printf(" output file not created\n");
        return RET_ERR;
    }
    /*
     *  if a file can be made but it is unsuitable for further processing,
     *  report a status of 1.  conditions are:
     *  - no navigation (ILT code = 0, or lat, lon limits all =
     */
    if ((gattr.ilt_flags) == 0 || (gattr.limits[0] == gattr.limits[1])
            || (gattr.limits[2] == gattr.limits[3])) {
        printf("%s - Navigation is unavailable for this granule\n", argv[0]);
        return RET_ERR;
    } else
        return RET_GOOD;
}
