/* ==================================================================== */
/*                                                                      */
/* SWl01 - converts seawifs level-0 data to level-1                     */
/*                                                                      */
/* Synopsis:							  	*/
/*									*/
/*     SWl01 input-L0-filename [output-L1-directory]			*/
/*									*/
/* Description:                                                         */
/*                                                                      */
/*     The program reads the input level-0 file and generates one or    */
/*     more level-1a hdf-formatted output files.  One file will be      */
/*     generated for each contiguous scene, and for HRPT data there is  */
/*     only one scene.  For stored GAC, there is one scene for each     */
/*     orbit, and stored LAC scenes are broken when the data type       */
/*     changes or the time between frames is greater than n seconds.    */
/*                                                                      */
/*     The output file(s) will be named Syyyydddhhmmss.L1A_*, where the */
/*     trailing suffix will indicate the HRPT station ID, or the stored */
/*     data type.  Each output L1A file will be accompanied by a meta   */
/*     data file with the same name and a ".meta" extension.            */
/*                                                                      */
/*     The processing steps are as follows:                             */
/*                                                                      */
/*     1 Read through input file and generate index of frame content    */
/*       and quality, as well as scene break points.                    */
/*                                                                      */
/*     2 Read selected frames to accumulate raw GPS telemetry, and      */
/*       fit GPS to orbit model to produce filtered orbit vectors for   */
/*       the time period over which the input data extends.             */
/*                                                                      */
/*     3 For each scene, read the state-of-health and instrument tlm    */
/*       from the "good" quality frames, and generate navigation for    */
/*       every scan.                                                    */
/*                                                                      */
/*     4 For each scene, read image data from the "good" quality frames,*/
/*       combine with the navigation information, and write to L1A file.*/
/*                                                                      */
/* Written By:                                                          */
/*                                                                      */
/*     Bryan A. Franz 							*/
/*     General Sciences Corp.                                           */
/*     27 September 1997                                                */
/*									*/
/* =====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <libgen.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "swl0_proto.h"
#include "swl1_hdf.h"
#include "swl1_netcdf.h"

#define CMD_ARGS "mhgt:b:s:e:d:n:f:c"

void usage(char *file) {
    printf("%s %s (%s %s)\n",
            file, L01VERSION, __DATE__, __TIME__);
    if (endianess() == 1)
        printf("      (byte-swapping on)\n");
    printf("\nUsage: %s [-%s] level0_filename [output-file-name or dir]\n",
            file, CMD_ARGS);
    printf("       -m : turn-off metadata file generation.\n");
    printf("       -h : force hrpt processing.\n");
    printf("       -g : force GAC processing.\n");
    printf("       -t n : limit timetags to +/- n-days around data start time.\n");
    printf("       -b n : maximum number of bit errors allowed per frame (def=50).\n");
    printf("       -s n : fix gain setting for HRPT (def=0, -1 to determine from telemetry).\n");
    printf("       -d n : stop-time delta (seconds).\n");
    printf("       -e file : timing anomaly update filename, GAC only. (def=NULL)\n");
    printf("                 If not provided, no timing error search will be done.\n");
    printf("       -n environment : runtime environment: seadas or sdps (def=seadas).\n");
    printf("       -f file : station info file (def=$HRPT_STATION_IDENTIFICATION_FILE).\n");
    printf("       -c : generate a NetCDF file instead of HDF4.\n");

    exit(FATAL_ERROR);
}


/* -------------------------------------------------------------------- */
/*                            main                                      */

/* -------------------------------------------------------------------- */
int main(int argc, char* argv[]) {
    swl0indx indx; /* L0 file content and quality */
    INT16 nscenes = 0; /* Number of scenes in L0 file */
    swl1rec l1rec[5]; /* Array of L1 records         */
    INT16 nl1rec; /* Number of L1 recs (1 or 5)  */
    INT16 is; /* Scene number                */
    INT16 iframe; /* Minor frame number          */
    INT16 irec; /* Record number               */
    INT32 sceneScanNum; /* Scan number within scene    */

    INT32 norb; /* Number of orbit input recs  */
    orbit_sType orbit; /* Processed orbit data        */
    INT32 nframes; /* Number of input recs for nav*/
    navctl_sType navctl; /* Navigation control info     */
    navqc_sType navqc; /* Navigation quality info     */
    INT32 nlines; /* Number of navblk recs       */
    tilt_states_sType tiltblk; /* Tilt states per scene       */
    FLOAT32 xnodel; /* Longitude of dscending node */
    INT32 tnode; /* Time of node crossing       */

    swl0ctl l0ctl; /* Control str for L0 indexing */
    char *l0file = NULL; /* Input L0 file path          */
    char l1file[FILENAME_MAX]; /* Output L1A file path        */
    char *outpath = "."; /* Output file dir or path     */
    char *metadir; /* Output metafile directory   */
    char metafile[FILENAME_MAX]; /* Output metafile path     */
    char *outdir; /* Output file dir             */
    char *filename; /* Output L1A filename         */
    unsigned char dataType; /* Flag for file naming        */
    struct stat file_stat; /* File status buffer          */
    char *tmppath; /* Temporary file path         */
    int haveFilename = 0; /* If set, don't create fileame*/
    int wantMeta = 1; /* Generate meta data by def   */
    char proccon[1024] = ""; /* Store command-line sequence */
    char proclog[1024] = ""; /* Store processing log        */
    char timefile[FILENAME_MAX]; /* time shift filename      */
    char *updfile = NULL;
    char stationInfoFile[FILENAME_MAX]; /* station filename  */
    int wantNetCDF = 0;

    int i;
    INT32 status = 0;

    /* These are declared static to force solaris to allocate the full  */
    /* space at load time.  Otherwise, the code core dumps at startup   */

    static input_sType orbinp[MAXFRAMES]; /* Input for orbit filtering  */
    static input_sType navinp[MAXFRAMES]; /* Input for navigation       */
    static navblk_sType navblk[MAXFRAMES]; /* Processed navigation data  */
    static swl0scene scene[MAXSCENES]; /* Scene index structure      */

    /* Parameters for getopt() */
    extern int opterr;
    extern int optind;
    extern char *optarg;
    int c;

    if (argc > 1 && !strcmp(argv[1], "--version")) {
        printf("%s %s (compiled %s, %s)\n",
                argv[0], L01VERSION, __DATE__, __TIME__);
        if (endianess() == 1)
            printf("      (byte-swapping on)\n");
        exit(0);
    }

    /*									*/
    /* Load telemetry decomm commons used by swfnav library 		*/
    /*									*/

    /* MDM Oct. 15, 2004 */
    /* due to the linker not including these block data routines, I've put
       them inline in the one place each of them is used.. therefore they
       are now completely eliminated and Bryan's calls are no longer req'd */
    /* acs_block_(); */
    /* ins_block_(); */

    /*									*/
    /* Enable line buffering for better log-file monitoring		*/
    /*									*/
    setvbuf(stdout, NULL, _IOLBF, 0);

    /*									*/
    /* Log command-line sequence to process-control string              */
    /*									*/
    strcpy(proccon, argv[0]);
    for (i = 1; i < argc; i++) {
        strcat(proccon, " ");
        strcat(proccon, argv[i]);
    }


    /*									  */
    /* Initialize control structure                                       */
    /*									  */
    l0ctl.fileType = -1; /* Get filetype from header           */
    l0ctl.timerangeFactor = 0; /* Get timerange limits from header   */
    l0ctl.maxBitErrors = 50; /* Max bit errors before tossing frame*/
    l0ctl.gainSetting = 0; /* Forced gain setting, for HRPT only */
    l0ctl.stopTimeDelta = 1800; /* Stop-time delta (seconds)          */
    l0ctl.env = SEADAS; /* Environment flag                   */
    l0ctl.progname = argv[0]; /* Program name for meta             */
    l0ctl.stationInfoFile = NULL; /* File containing station meta info */


    /*									*/
    /* Process command-line arguments                                   */
    /*									*/
    while ((c = getopt(argc, argv, CMD_ARGS)) != EOF) {
        switch (c) {
        case 'm':
            wantMeta = 0;
            break;
        case 'b':
            l0ctl.maxBitErrors = atoi(optarg);
            break;
        case 's':
            l0ctl.gainSetting = atoi(optarg);
            break;
        case 't':
            l0ctl.timerangeFactor = atoi(optarg);
            break;
        case 'd':
            l0ctl.stopTimeDelta = atoi(optarg);
            break;
        case 'h':
            l0ctl.fileType = HRPT;
            break;
        case 'g':
            l0ctl.fileType = GAC;
            break;
        case 'e':
            strncpy(timefile, optarg, FILENAME_MAX);
            updfile = timefile;
            break;
        case 'n':
            if (strcmp(optarg, "sdps") == 0)
                l0ctl.env = SDPS;
            else
                l0ctl.env = SEADAS;
            break;
        case 'f':
            strncpy(stationInfoFile, optarg, FILENAME_MAX);
            l0ctl.stationInfoFile = stationInfoFile;
            break;
        case 'c':
            wantNetCDF = 1;
            break;
        default:
            usage(argv[0]);
            break;
        }
    }
    switch (argc - optind + 1) {
    case 3:
        l0file = argv[optind + 0];
        outpath = argv[optind + 1];
        break;
    case 2:
        l0file = argv[optind + 0];
        break;
    default:
        usage(argv[0]);
        break;
    }


    printf("\nBegin %s Version %s Processing for %s using the ",
            l0ctl.progname, L01VERSION, l0file);
    if (l0ctl.env == SEADAS)
        printf("SeaDAS environment\n\n");
    else
        printf("SDPS environment\n\n");


    /*									*/
    /* Extract output directory from full (or partial) path.            */
    /*									*/
    if (stat(outpath, &file_stat) == 0 && S_ISDIR(file_stat.st_mode)) {
        /* output path contains output directory */
        outdir = strdup(outpath);
    } else {
        /* output path contains full directory and filename */
        haveFilename = 1;
        tmppath = strdup(outpath);
        outdir = strdup(dirname(tmppath));
        free(tmppath);
    }

    /* Set metadata output dir */
    if ((metadir = getenv("L1A_META_DIR")) == NULL) {
        metadir = strdup(outdir);
    }

    /*									*/
    /* Get L0 content and quality index					*/
    /*	                                                                */
    printf("\nGenerating Level-0 file index ...\n");
    nframes = getl0indx(l0file, &l0ctl, &indx);
    if (nframes <= 0) {
        printf("-E- %s: no frames found in %s\n", argv[0], l0file);
        exit(FATAL_ERROR);
    }
    locate_temporal_anomalies(&indx, updfile);
    printindx(&indx);


    /*									*/
    /* Create scene index						*/
    /*									*/
    printf("\nGenerating Level-0 scene index ...\n");
    nscenes = getl0scene(&indx, scene);
    if (nscenes < 1) {
        printf("-E- %s: no valid scenes found in %s\n", argv[0], l0file);
        exit(FATAL_ERROR);
    }
    printscene(nscenes, scene);


    /*									*/
    /* Gather telemetry for orbit determination	and compute filtered   	*/
    /* orbit vectors.							*/
    /*									*/
    printf("\nGathering raw GPS orbit data ...\n");
    norb = getorbdata(&indx, orbinp);
    printf("Found %d frames with valid GPS data.\n", norb);
    if (norb < MINORBVEC) {
        printf("-E- %s: insufficient GPS telemetry.\n", argv[0]);
        exit(FATAL_ERROR);
    }
    printf("\nGenerating filtered GPS orbit vectors ...\n");
    initnav_(orbinp, &norb, &navctl, &navqc, &orbit, &status);
    if (status != 0) {
        printf("-E- %s: unable to process orbit vectors.\n", argv[0]);
        exit(FATAL_ERROR);
    }

    /*									*/
    /* Process each scene	   				        */
    /*									*/
    for (is = 0; is < nscenes; is++) {

        printf("\nProcessing scene %d ...\n", is);
        
        /*                                                              */
        /* Construct output filename, or extract from output path.      */
        /*                                                              */
        if (!haveFilename) {
            if (scene[is].type == HRPT)
                dataType = 16;
            else
                dataType = scene[is].mnftype;
            filename = L1aFilename_netcdf(&l0ctl, scene[is].stime, dataType);
            sprintf(l1file, "%s/%s", outpath, filename);
        } else {
            strcpy(l1file, outpath);
            filename = strdup(basename(outpath));
        }
        printf("\nOutput filename is %s\n", filename);
        printf("Output file dir: %s\n", outdir);
        printf("Metafile dir:    %s\n\n", metadir);

        /*								*/
        /* Load navigation input structure, skip scene if error occurs  */
        /*								*/
        printf("Loading navigation data\n");
        nframes = getnavdata(&scene[is], navinp);
        if (nframes != scene[is].nrec) {
            printf("-E- %s: error reading nav data for scene %d\n",
                    argv[0], is);
            status++;
            continue;
        }

        /*								*/
        /* Compute scene navigation.                                    */
        /*								*/
        printf("Navigating scene\n");

        swfnav_(&navqc, &navctl, navinp, &nframes, &orbit, &nlines,
                navblk, &tiltblk, &xnodel, &tnode);

        if (getl0scene_nav(xnodel, tnode, navblk, &tiltblk, &scene[is]) != 0) {
            printf("-E- %s: no valid navigation for scene %d\n",
                    argv[0], is);
            status++;
            continue;
        }
        printnav(&scene[is]);


        /*								*/
        /* Write L1A file for this scene                                */
        /*								*/

        if (wantNetCDF) {
            if (CreateL1aFile_netcdf(l1file, &scene[is], proccon, proclog, &l0ctl) != 0) {
            printf("-E- %s: error creating L1A file for scene %d\n",
                    argv[0], is);
            status++;
            continue;
            }
        } else {
            if (CreateL1aFile(l1file, &scene[is], proccon, proclog, &l0ctl) != 0) {
            printf("-E- %s: error creating L1A file for scene %d\n",
                    argv[0], is);
            status++;
            continue;
            }
        }
        

        sceneScanNum = 0;
        for (iframe = 0; iframe < scene[is].nrec; iframe++) {

            if ((iframe % 500) == 0)
                printf("Writing frame %5d of scene %2d\n", iframe, is);

            nl1rec = getl1rec(iframe, &scene[is], &l0ctl,
                    navinp, navblk, &tiltblk, l1rec);

            for (irec = 0; irec < nl1rec; irec++) {
                addL1Metrics(sceneScanNum, &l1rec[irec]);

                // toggle which scan data to use, netcdf or hdf4
                if (wantNetCDF) {
                    if (WriteScanData_netcdf(sceneScanNum, &l1rec[irec]) != 0) {
                        printf("-E- %s: error writing to L1A file for scene %d\n",
                                argv[0], is);
                        break;
                    }   
                } else {
                    if (WriteScanData(sceneScanNum, &l1rec[irec]) != 0) {
                        printf("-E- %s: error writing to L1A file for scene %d\n",
                                argv[0], is);
                        break;
                    }
                }
                
                sceneScanNum++;
            }

            /* If an error occurs, end processing for this scene        */
            if (irec != nl1rec) {
                status++;
                break;
            }
        }
        if (wantNetCDF)
            CloseL1aFile_netcdf(getL1Metrics());
        else
            CloseL1aFile(getL1Metrics());

        /* Write associated metadata file */
        if (wantMeta) {
            sprintf(metafile, "%s/%s.meta", metadir, filename);
            mkmeta(metafile, l1file, &scene[is], &l0ctl);
        }
    }

    printf("\nNumber of scenes processed = %d\n", nscenes);
    printf("Number of scenes failed    = %d\n\n", status);

    exit(status);

}





