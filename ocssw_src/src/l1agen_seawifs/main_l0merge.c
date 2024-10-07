/* ==================================================================== */
/*                                                                      */
/* SWl0merge - merges multiple L0 files to one optimal file             */
/*                                                                      */
/* synopsisl:							  	*/
/*									*/
/*     SWl0merge input-L0-listfile [output-L0-file]			*/
/*									*/
/* Description:                                                         */
/*                                                                      */
/* Takes as input a list of Level-0 filenames and writes a merged       */
/* Level-0 file.  Frames of common type and time are resolved by        */
/* quality tests based on bit-error counts. Only the best available     */
/* frame is retained in the merged output.                              */
/*                                                                      */
/* Currently, the program can merge up to 100 files, but this limit can */
/* be easily increased.                                                 */
/*                                                                      */
/* Written By:                                                          */
/*                                                                      */
/*     Bryan A. Franz 							*/
/*     General Sciences Corp.                                           */
/*     13 Septempber 2002                                               */
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
#include <time.h>

#include "swl0_proto.h"
#define MAXFILES 100
#define CMD_ARGS "gt:"

/* Function to read a list of filenames from a specified file */
INT32 read_file_list(char *filename, char files[][FILENAME_MAX], INT32 maxfiles) {
    FILE *fp = NULL;
    int nfiles = 0;

    if ((fp = fopen(filename, "r")) == NULL) {
        fprintf(stderr,
                "-W- %s line %d: unable to open %s for reading\n",
                __FILE__, __LINE__, filename);
        return (1);
    }

    while (nfiles < maxfiles && (fscanf(fp, "%s\n", files[nfiles]) != EOF))
        if (strlen(files[nfiles]) > 1)
            nfiles++;

    fclose(fp);

    return (nfiles);
}

/* qsort comparison function for the merged index */
int compmindx(swl0mindx *rec1, swl0mindx *rec2) {
    int retval = 0;

    /* frame type */
    if (rec1->qual->mnftype > rec2->qual->mnftype) retval = 1;
    if (rec1->qual->mnftype < rec2->qual->mnftype) retval = -1;

    /* frame time within type (+/- 1/60 sec) */
    if (retval == 0) {
        if (rec1->qual->time > (rec2->qual->time + DTLAC / 10)) retval = 1;
        if (rec1->qual->time < (rec2->qual->time - DTLAC / 10)) retval = -1;
    }

    /* quality within time */
    if (retval == 0) {
        if (rec1->qual->pixVariance > rec2->qual->pixVariance) retval = 1;
        if (rec1->qual->pixVariance < rec2->qual->pixVariance) retval = -1;
    }
    if (retval == 0) {
        if (rec1->qual->errBits > rec2->qual->errBits) retval = 1;
        if (rec1->qual->errBits < rec2->qual->errBits) retval = -1;
    }
    if (retval == 0) {
        if (rec1->qual->sgaError > rec2->qual->sgaError) retval = 1;
        if (rec1->qual->sgaError < rec2->qual->sgaError) retval = -1;
    }
    if (retval == 0) {
        if (rec1->qual->sacError > rec2->qual->sacError) retval = 1;
        if (rec1->qual->sacError < rec2->qual->sacError) retval = -1;
    }
    if (retval == 0) {
        if (rec1->qual->saaError > rec2->qual->saaError) retval = 1;
        if (rec1->qual->saaError < rec2->qual->saaError) retval = -1;
    }

    /* filename within quality (LAC files take precedence) */
    if (retval == 0) {
        if (strstr(rec1->filename, "LAC") != NULL)
            retval = -1;
        else if (strstr(rec2->filename, "LAC") != NULL)
            retval = 1;
        else if (strstr(rec1->filename, "HNSG") != NULL)
            retval = -1;
        else if (strstr(rec2->filename, "HNSG") != NULL)
            retval = 1;
        else if (strstr(rec1->filename, "HDUN") != NULL)
            retval = -1;
        else if (strstr(rec2->filename, "HDUN") != NULL)
            retval = 1;
        else
            retval = strcmp(rec1->filename, rec2->filename);
    }

    return (retval);
}

/* builds a seawifs-style filename from time and datatype */
char *ml0filename(double time, unsigned char dataType) {

    static char filename[23]; /* "Syyyydddhhmmss.L0_tttt\0" */
    struct tm *t;
    time_t itime;
    double rint(double);
    char suffix[5] = "MLAC";

    if (dataType == GACTYPE) strcpy(suffix, "MGAC");

    itime = (time_t) rint(time); /* Round to nearest second. */
    t = gmtime(&itime);

    sprintf(
            filename,
            "S%4d%03d%02d%02d%02d.L0_%.4s",
            t->tm_year + 1900,
            t->tm_yday + 1,
            t->tm_hour,
            t->tm_min,
            t->tm_sec,
            suffix
            );

    return (filename);
}

void usage(char *file) {
    printf("%s %s (%s %s)\n",
            file, L01VERSION, __DATE__, __TIME__);
    if (endianess() == 1)
        printf("      (byte-swapping on)\n");
    printf("\nUsage: %s [-%s] input-L0-listfile  [output-file-name or dir]\n",
            file, CMD_ARGS);
    printf("       -g : force GAC processing (allows mixed frame types).\n");
    printf("       -t n : limit timetags to +/- n-days around data start time. (GAC)\n");

    exit(FATAL_ERROR);
}


/* -------------------------------------------------------------------- */
/*                            main                                      */

/* -------------------------------------------------------------------- */
int main(int argc, char* argv[]) {
    char *l0listfile = NULL; /* Input L0 list file          */
    char outpath[FILENAME_MAX]; /* Output file dir or path  */
    char l0files[MAXFILES][FILENAME_MAX]; /* Input L0 files */
    char *filename; /* Output L0 filename          */
    unsigned char dataType = LACTYPE; /* Flag for file naming        */
    struct stat file_stat; /* File status buffer          */

    INT32 irec; /* Record (frame) number       */
    INT32 nfiles = 0; /* Number of input L0 files    */
    INT32 nframes = 0; /* Number of frames in L0 file */
    INT32 nrecs = 0; /* Number of merged frames     */
    INT32 ntotal = 0; /* Total frames from all files */

    INT32 filecnt[MAXFILES]; /* Count select frames by file */

    swl0ctl l0ctl; /* Control str for L0 indexing */
    swl0indx indx[MAXFILES]; /* L0 file content and quality */
    swl0mindx *mindx; /* Merged L0 file index        */


    FILE *ifp = NULL; /* Input L0 file pointer       */
    FILE *ofp = NULL; /* Output merged file pointer  */

    swl0hdr hdr; /* Level-0 file header         */
    BYTE mnf[L0LEN]; /* Raw minor frame             */
    INT32 fileBytePos; /* File position of L0 frame   */

    FLOAT64 lasttime;
    INT32 lastfilenum;
    INT32 i, j;


    /* Parameters for getopt() */
    extern int opterr;
    extern int optind;
    extern char *optarg;
    int c;


    /*									*/
    /* Enable line buffering for better log-file monitoring		*/
    /*									*/
    setlinebuf(stdout);


    /*									*/
    /* Initialize control structure                                     */
    /*									*/
    l0ctl.fileType = HRPT; /* Assume HRPT (ignore non-LAC frames)*/
    l0ctl.timerangeFactor = 0; /* Get timerange limits from hdr (GAC)*/
    l0ctl.maxBitErrors = 500; /* No effect on SWl0merge             */
    l0ctl.gainSetting = 0; /* No effect on SWl0merge             */
    l0ctl.stopTimeDelta = 3600; /* Stop-time delta (seconds)          */


    /*									*/
    /* Process command-line arguments                                   */
    /*									*/
    while ((c = getopt(argc, argv, CMD_ARGS)) != EOF) {
        switch (c) {
        case 't':
            l0ctl.timerangeFactor = atoi(optarg);
            break;
        case 'g':
            l0ctl.fileType = GAC;
            dataType = GACTYPE;
            break;
        default:
            usage(argv[0]);
            break;
        }
    }
    switch (argc - optind + 1) {
    case 3:
        l0listfile = argv[optind + 0];
        strcpy(outpath, argv[optind + 1]);
        break;
    case 2:
        l0listfile = argv[optind + 0];
        strcpy(outpath, ".");
        break;
    default:
        usage(argv[0]);
        break;
    }


    /*									*/
    /* Build L0 file list						*/
    /*									*/
    nfiles = read_file_list(l0listfile, l0files, MAXFILES);
    if (nfiles <= 0) {
        printf("-E- %s: error reading list file %s\n", argv[0], l0listfile);
        exit(FATAL_ERROR);
    }
    printf("\nMerging %d files.\n", nfiles);
    for (i = 0; i < nfiles; i++) printf("File %d: %s\n", i + 1, l0files[i]);
    printf("\n");


    /*									*/
    /* Loop over each L0 file and create content and quality index	*/
    /*									*/
    for (i = 0; i < nfiles; i++) {
        printf("\nGenerating Level-0 file index for %s\n", l0files[i]);
        filecnt[i] = 0;
        nframes = getl0indx(l0files[i], &l0ctl, &indx[i]);
        if (nframes <= 0) {
            printf("-E- %s: no frames found in %s\n", argv[0], l0files[i]);
            exit(FATAL_ERROR);
        }
        nrecs += nframes;
    }


    /*									*/
    /* Allocate space for merged index of all files, all frames		*/
    /*									*/
    if ((mindx = (swl0mindx *) calloc(nrecs, sizeof (swl0mindx))) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: error allocating memory for merge index array\n",
                __FILE__, __LINE__);
        return (1);
    }


    /*									*/
    /* Concatenate frame quality indices into merged index   		*/
    /*									*/
    irec = 0;
    for (i = 0; i < nfiles; i++) {
        printindx(&indx[i]);
        for (j = indx[i].srec; j <= indx[i].erec; j++) {
            if (indx[i].rec[j].tRanError == 0 &&
                    indx[i].rec[j].tSeqError == 0 &&
                    indx[i].rec[j].tDifError == 0 &&
                    indx[i].rec[j].tShfError == 0 &&
                    indx[i].rec[j].scidError == 0) {

                strcpy(mindx[irec].filename, l0files[i]);
                mindx[irec].filenum = i;
                mindx[irec].framenum = j;
                mindx[irec].qual = &(indx[i].rec[j]);
                irec++;

            }
        }
    }
    ntotal = nrecs = irec;


    /*									*/
    /* Sort merged index by time, and quality within time		*/
    /*									*/
    qsort(mindx, nrecs, sizeof (swl0mindx),
            (int (*)(const void *, const void *)) compmindx);


    /*									*/
    /* Reduce merged index to unique frame times. Assumes frames have   */
    /* sorted into time order, with "best" quality frame preceeding any */
    /* lower quality frames. Thus, first frame of each time is selected.*/
    /*									*/
    lasttime = 0.0;
    irec = 0;
    for (i = 0; i < nrecs; i++) {
        if (fabs(mindx[i].qual->time - lasttime) > DTLAC / 10.) {
            lasttime = mindx[i].qual->time;
            mindx[irec] = mindx[i];
            irec++;
        }
    }
    nrecs = irec;


    /*                                                              */
    /* Construct output filename, or extract from output path.      */
    /*                                                              */
    if (stat(outpath, &file_stat) == 0 && S_ISDIR(file_stat.st_mode)) {
        filename = ml0filename(mindx[0].qual->time, dataType);
        strcat(outpath, "/");
        strcat(outpath, filename);
    }
    printf("\nOutput filename is %s\n", outpath);


    /*                                                              */
    /* Open merged L0 file for writing.                             */
    /*                                                              */
    if ((ofp = fopen(outpath, "w")) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: unable to open %s for writing\n",
                __FILE__, __LINE__, outpath);
        exit(FATAL_ERROR);
    }


    /*                                                              */
    /* Construct and write merged L0 file header.                   */
    /*                                                              */
    strcpy((char *) hdr.id, "CWIF");
    hdr.numlines = nrecs;
    if (dataType == GACTYPE)
        hdr.type = 0;
    else
        hdr.type = 1;
    hdr.startTime = mindx[0].qual->time;
    hdr.stopTime = mindx[nrecs - 1].qual->time;
    if (endianess() == 1) {
        swapc_bytes((char *) &hdr.startTime, 4, 1);
        swapc_bytes((char *) &hdr.stopTime, 4, 1);
    }
    fwrite(&hdr, 1, 512, ofp);


    /*                                                              */
    /* Read selected frames, and write to merged output L0 file.    */
    /*                                                              */
    lastfilenum = -1;
    for (irec = 0; irec < nrecs; irec++) {

        /* Open L0 file of selected frame, if neccessary */
        if (mindx[irec].filenum != lastfilenum) {
            if (ifp != NULL) fclose(ifp);
            if ((ifp = fopen(mindx[irec].filename, "r")) == NULL) {
                fprintf(stderr,
                        "-E- %s line %d: unable to open %s for reading\n",
                        __FILE__, __LINE__, mindx[irec].filename);
                exit(FATAL_ERROR);
            }
            lastfilenum = mindx[irec].filenum;
        }

        /* Position file pointer to start of selected frame */
        fileBytePos = sizeof (swl0hdr) + mindx[irec].framenum * (L0LEN);
        if (fseek(ifp, fileBytePos, SEEK_SET) != 0) {
            fprintf(stderr,
                    "-E- %s line %d: error reading %s\n",
                    __FILE__, __LINE__, mindx[irec].filename);
            exit(FATAL_ERROR);
        }

        /* Read frame */
        if (fread(mnf, L0LEN, 1, ifp) != 1) {
            fprintf(stderr,
                    "-E- %s line %d: error reading %s\n",
                    __FILE__, __LINE__, mindx[irec].filename);
            exit(FATAL_ERROR);
        }

        /* Write frame to merged output */
        if (fwrite(mnf, L0LEN, 1, ofp) != 1) {
            fprintf(stderr,
                    "-E- %s line %d: error writing frame %d to %s\n",
                    __FILE__, __LINE__, irec, outpath);
            exit(FATAL_ERROR);
        }

        filecnt[mindx[irec].filenum]++;

    }

    fclose(ofp);
    fclose(ifp);

    printf("\n\n%d frames selected from %d frames in %d files\n", nrecs, ntotal, nfiles);
    printf("Compression ratio is %f\n", 1.0 * nrecs / ntotal);
    for (i = 0; i < nfiles; i++)
        printf("%10d frames selected from %s\n", filecnt[i], l0files[i]);

    exit(0);
}





