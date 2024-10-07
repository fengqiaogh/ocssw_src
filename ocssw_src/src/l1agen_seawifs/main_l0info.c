/* ==================================================================== */
/*                                                                      */
/* SWl0info - checks level-0 file structure, frame content, and errors  */
/*                                                                      */
/* Synopsis:							  	*/
/*                                                                      */
/*     SWl0info L0-filename                                             */
/*									*/
/* Description:                                                         */
/*                                                                      */
/* Reads the level-0 file and reports statistics on frame types and     */
/* frame errors.  If the file can not be opened or does not appear to   */
/* be a seawifs level-0 file (based on the header), the error status    */
/* will be set to -1, else it will be set to the number of frames read. */
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
#include <time.h>

#include "swl0_proto.h"

#define CMD_ARGS "h"

/* ----------------------------------------------------------------------- */
/*                            main                                         */

/* ----------------------------------------------------------------------- */
int main(int argc, char* argv[]) {
    swl0indx indx;
    swl0scene scene[MAXSCENES];
    INT32 nframes;
    INT32 nscenes = 0;
    swl0ctl l0ctl; /* Control str for L0 indexing */

    /* Parameters for getopt() */
    extern int opterr;
    extern int optind;
    extern char *optarg;
    int c;

    l0ctl.fileType = -1; /* Get filetype from header           */

    /*									*/
    /* Process command-line arguments                                   */
    /*									*/
    while ((c = getopt(argc, argv, CMD_ARGS)) != EOF) {
        switch (c) {
        case 'h':
            l0ctl.fileType = HRPT;
            break;
        case 'g':
            l0ctl.fileType = GAC;
            break;
        default:
            printf("%s %s (%s %s)\n",
                    "l0info_seawifs", L01VERSION, __DATE__, __TIME__);
            printf("Usage: %s [-h -g] level0_filename\n", argv[0]);
            exit(FATAL_ERROR);
            break;
        }
    }
    switch (argc - optind + 1) {
    case 2:
        break;
    default:
        printf("%s %s (%s %s)\n",
                "l0info_seawifs", L01VERSION, __DATE__, __TIME__);
        printf("Usage: %s [-h -g] level0_filename\n", argv[0]);
        exit(FATAL_ERROR);
        break;
    }

    /*									*/
    /* Initialize control structure                                     */
    /*									*/
    l0ctl.timerangeFactor = 0; /* Get timerange limits from header   */
    l0ctl.maxBitErrors = 50; /* Max bit errors before tossing frame*/
    l0ctl.gainSetting = 0; /* Forced gain setting, for HRPT only */
    l0ctl.stopTimeDelta = 2700; /* Stop-time delta (seconds)          */

    nframes = getl0indx(argv[optind], &l0ctl, &indx);
    if (nframes > 0) {
        printindx(&indx);
        nscenes = getl0scene(&indx, scene);
        printscene(nscenes, scene);
    }

    if (nscenes < 1)
        exit(1);
    else
        exit(nframes);
}




