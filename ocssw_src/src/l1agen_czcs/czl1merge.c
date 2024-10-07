#include <libgen.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <sys/utsname.h>

#include "l1czcs.h"

#include <passthebuck.h>
#include <GetStationInfo.h>
#include <timeutils.h>

#define MAXFILES 80
#define PROGRAM "czl1merge"
#define VERSION "1.1"

int main(int argc, char *argv[])
/*******************************************************************

   czl1merge

   purpose: merge several CZCS L1 files together into 1 file

   Returns type: int - 0 if all good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               argc             I      count of command line 
                                                args
      char *            argv[]           I      command line arguments
                                                They are:
                                                [1] input L1 list file
                                                [2] merged file name or
                                                just the path for it

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       9 Aug 2004      Original development

 *******************************************************************/ {
    struct stat file_stat;
    l1_data_struc l1_data; /* output L1 data storage */
    gattr_struc gattr; /* and global attributes */
    mstr_struc mstr; /* master storage of time and quality */
    timqual_struc init_info; /* time, quality from candidate file */
    int n_mstr, i, mstr_last, in_ptr, mstr_ptr, nfiles, ifil, new_start,
            new_end, rstat, n_ctl_pt, mtch_st, olap_end, nscans, init = 1,
            orb_st_day, num_used[MAXFILES], scans_per_file[50];
    int32_t mmsec, imsec;
    char infiles[MAXFILES][FILENAME_MAX]; /* list of input files */
    short file_used[MAXFILES]; /* 1 if the input L1 was used to make output */
    char outfile[FILENAME_MAX], outpath[FILENAME_MAX];
    char inflist[FILENAME_MAX];
    char seadas_vs[64], prog_version[64];
    /*
     *  for the library functions, define here
     */
    int get_swf_def_meta1();
    int gen_soft_id();
    int get_current_time();

    for (i = 0; i < MAXFILES; i++) {
        file_used[i] = 0;
        num_used[i] = 0;
    }
    olap_end = 0;
    n_ctl_pt = 0;
    /*
     *  get inputs
     */
    if (argc == 2) {
        strcpy(inflist, argv[1]);
        strcpy(outpath, ".");
    } else if (argc == 3) {
        strcpy(inflist, argv[1]);
        strcpy(outpath, argv[2]);
    } else
        usage(basename(argv[0]));
    /*
     *  get the start times from the file list and order by start time
     */
    nfiles = read_file_list(inflist, (char**) infiles, 80);
    if (nfiles <= 0) {
        printf("-E- %s: error reading list file %s\n", argv[0], inflist);
        exit(-1);
    }

    /*
     *  inflist can be sorted with qsort and thus, end up
     *  in the correct order (Use strcmp).  
     */
    qsort((char *) infiles, nfiles, sizeof ( *infiles),
            (int (*)(const void*, const void*))strcmp);

    printf("\nordered list of input L1 CZCS files:\n");
    printf("numb   name\n");
    for (ifil = 0; ifil < nfiles; ifil++)
        printf("%4d   %s\n", ifil, infiles[ifil]);
    /*
     *  allocate the storage for the master list of time and quality
     */
    n_mstr = nfiles * 970;
    mstr.msec = (int32_t *) malloc(n_mstr * sizeof ( int32_t));
    mstr.exist = (short *) calloc(n_mstr, sizeof ( short));
    mstr.qual = (short *) malloc(n_mstr * sizeof ( short));
    mstr.ds_num = (short *) malloc(n_mstr * sizeof ( short));
    mstr.out_scan = (int32_t *) malloc(n_mstr * sizeof ( int32_t));
    mstr.scan = (short *) malloc(n_mstr * sizeof ( short));
    mstr.in_msec = (int32_t *) malloc(n_mstr * sizeof ( int32_t));
    mstr.in_exist = (short *) malloc(n_mstr * sizeof ( short));
    mstr.in_qual = (short *) malloc(n_mstr * sizeof ( short));
    mstr.in_scan = (short *) malloc(n_mstr * sizeof ( short));
    /*
     *  loop through the input datasets and construct the master list 
     *  of the best scans
     */
    mstr_last = -1;
    orb_st_day = -1; /* start day of the orbit, to get msec in mstr good */

    for (ifil = 0; ifil < nfiles; ifil++) {
        printf("czl1merge: processing input file %d: %s\n", ifil, infiles[ifil]);
        new_start = -1;
        new_end = -1;

        if ((rstat = cztimqual(infiles[ifil], &init_info, &orb_st_day)) == 0) {
            scans_per_file[ifil] = init_info.nscan;
            /*
             *  There is data to use from this file.  if the first usable set,
             *  or new data is beyond the current data, just put into the 
             *  list as the current set.
             */
            n_ctl_pt = init_info.n_ctl_pt;

            if ((mstr_last == -1) || (mstr.msec[mstr_last] < init_info.msec[0])) {
                fill_mstr(&mstr_last, &mstr, &init_info, ifil, 0,
                        init_info.nscan);
            } else {
                /*
                 *  There is some overlap.  back up in the current msec to before 
                 *  the new data time
                 */
                mstr_ptr = mstr_last;
                mtch_st = 0;
                while (mtch_st == 0) {
                    if (mstr.msec[mstr_ptr] <= init_info.msec[0])
                        mtch_st = 1;
                    else
                        mstr_ptr--;
                }
                printf("czl1merge: new file matches at index %d, current end at %d\n",
                        mstr_ptr, mstr_last);
                /*
                 *  now, mstr_ptr points to the overlap start.  Flip through the 
                 *  msec in mtch and init_info to get matches (approximate in gaps)
                 */
                olap_end = 0;
                in_ptr = 0;
                while (olap_end == 0) {
                    mmsec = mstr.msec[mstr_ptr];
                    imsec = init_info.msec[in_ptr];
                    if (imsec == mmsec) {
                        /*
                         *  times match.  place the new data in the master list in 
                         *  the '_in' arrays
                         */
                        if (new_start == -1)
                            new_start = mstr_ptr;
                        new_end = mstr_ptr;

                        mstr.in_msec[mstr_ptr] = init_info.msec[in_ptr];
                        mstr.in_exist[mstr_ptr] = 1;
                        mstr.in_qual[mstr_ptr] = init_info.qual[in_ptr];
                        mstr.in_scan[mstr_ptr] = in_ptr;

                        in_ptr++;
                        mstr_ptr++;
                    } else if ((mstr.exist[mstr_ptr] == 0) &&
                            (imsec - mmsec < 40) && (imsec - mmsec > -40)) {
                        /*
                         *  incoming time matches with a rough predicted time in a gap.
                         *  do as above
                         */
                        if (new_start == -1)
                            new_start = mstr_ptr;
                        new_end = mstr_ptr;

                        mstr.in_msec[mstr_ptr] = init_info.msec[in_ptr];
                        mstr.in_exist[mstr_ptr] = 1;
                        mstr.in_qual[mstr_ptr] = init_info.qual[in_ptr];
                        mstr.in_scan[mstr_ptr] = in_ptr;

                        in_ptr++;
                        mstr_ptr++;
                    } else if (imsec > mmsec) {
                        mstr.in_exist[mstr_ptr] = 0;
                        mstr_ptr++;
                        printf("czl1merge: gap found in new data vs current at mstr_ptr: %d\n",
                                mstr_ptr);
                    } else if (mmsec > imsec)
                        in_ptr++;

                    /*
                     * set the end condition if we exhausted the current or input data
                     */
                    if ((mstr_ptr > mstr_last) || (in_ptr >= init_info.nscan))
                        olap_end = 1;
                } /* end loop to deal with ovelap filling of mstr */
                /*
                 *  The overlap is filled.  resolve the overlap
                 */
                printf("czl1merge: resolving overlap from index %d to %d\n",
                        new_start, new_end);
                olap_resolve(&mstr, mstr_last, new_start, new_end, ifil);
                /*
                 *  for remaining scans, if any, just fill the master array
                 */
                if (in_ptr < init_info.nscan) {
                    fill_mstr(&mstr_last, &mstr, &init_info, ifil, in_ptr,
                            init_info.nscan);
                }
            }
        } else
            exit(-1);

        /*
         *  free the arrays under init_info
         */
        free(init_info.qual);
        free(init_info.msec);
    }

    printf("\nFinished set-up of best data to merge, Outputting data\n\n");

    if (mstr_last == -1) {
        printf("No suitable data exists for this set of files\n");
        exit(-1);
    }
    /*
     *  the best scans are described in mstr.  Get the scans from the 
     *  contributing files to make the output file
     *
     *  first count the real # scans and find which files contribute
     *  note below that mstr_last is pointer to last entry, not the #
     *  of entrys, thus the <= in the for loop
     */
    nscans = 0;
    for (i = 0; i <= mstr_last; i++) {
        if (mstr.exist[i] != 0) {
            mstr.out_scan[i] = nscans;
            nscans++;
            file_used[ mstr.ds_num[i] ] = 1;
            num_used[ mstr.ds_num[i] ]++;
        }
    }
    /*
     *  set storage in the data structure
     */
    cz_dat_alloc(nscans, n_ctl_pt, 0, &l1_data);

    /*
     *  loop through the datasets and use cz_mov_scn() to place the scan lines 
     *  using mstr as guide
     */
    init = 1;
    for (i = 0; i < nfiles; i++)
        if (file_used[i] != 0) {
            if (
                    cz_mov_scn(init, i, infiles[i], &mstr, mstr_last, nscans,
                    &gattr, &l1_data) != 0) exit(-1);
            if (init == 1) init = 0;
        }

    gattr.scan_lines = nscans;
    /*
     *  update lat, lon metadata
     */
    cz_ll_upd(&l1_data, &gattr);
    /*
     *  adjust the metadata to reflect the new span of data
     */
    cz_meta_adj(&l1_data, &gattr);
    /*
     *  create output file name using SWl0merge methodology
     */
    if (stat(outpath, &file_stat) == 0 && S_ISDIR(file_stat.st_mode)) {
        strcat(outpath, "/");
        strcpy(outfile, basename(infiles[0]));
        strcpy(outfile + 19, "MLAC");
        strcat(outpath, outfile);
    }

    strcpy(gattr.proc_ctl, argv[0]);
    for (i = 1; i < argc; i++) {
        strcat(gattr.proc_ctl, " ");
        strcat(gattr.proc_ctl, argv[i]);
    }

    strcpy(gattr.input_files, basename(infiles[0]));
    if (nfiles > 1)
        for (i = 1; i < nfiles; i++) {
            strcat(gattr.input_files, " ");
            strcat(gattr.input_files, basename(infiles[i]));
        }
    /*
     *  set up some attrib info 
     */
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

    strcpy(gattr.datatype, "MLAC");
    get_time((char*) &gattr.process_time);
    /*
     *  call the czcs write routine to write out the gattr and data structs
     */
    printf("\n Writing data to output file...\n");
    if (czcs_l1_write(outpath, l1_data, gattr) == 0) {
        printf("\n\nCreated output file: %s \n", outpath);
        printf("With %d lines\n", nscans);
        printf("\nContributions from input files:\n\n");
        printf("seq #               file name  # lines   # lines contributed\n");
        for (ifil = 0; ifil < nfiles; ifil++) {
            strcpy(outfile, basename(infiles[ifil]));
            i = ifil + 1;
            printf("%5d  %22s    %5d    %5d\n", i, outfile, scans_per_file[ifil],
                    num_used[ifil]);
        }
    } else {
        printf(" output file not created\n");
        exit(-1);
    }
    /*
     *  close the output and end
     */
    exit(0);
}

void usage(char *file) {
    printf("%s %s (%s %s)\n",
            file, VERSION, __DATE__, __TIME__);
    printf("\nUsage: %s input-CZCS-listfile  [output-file-name or dir]\n",
            file);

    exit(-1);
}
