/* ==================================================================== */
/*                                                                      */
/* cut_elements - clip time-range of elements.dat file of SWl01.        */
/*                                                                      */
/* Synopsis:							  	*/
/*									*/
/*     cut_elements [-n num_days] infile outfile		        */
/*									*/
/* Description:                                                         */
/*                                                                      */
/*                                                                      */
/* Written By:                                                          */
/*                                                                      */
/*     Bryan A. Franz 							*/
/*     General Sciences Corp.                                           */
/*     8 March 2000                                                     */
/*									*/
/* =====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include "genutils.h"
#include "elements.h"
#include <string.h>
#include <timeutils.h>

#define CMD_ARGS "y:d:n:"

void usage(char *file) {
    printf("Usage: %s [-%s] input-elements-file output-elements-file\n",
            file, CMD_ARGS);
    printf("       -n n : number of days to include. [default=1]\n");

    exit(1);
}


/* -------------------------------------------------------------------- */
/*                            main                                      */

/* -------------------------------------------------------------------- */
int main(int argc, char* argv[]) {

    INT32 ndays = 1; /* Number of days to retain    */
    char *ifile = NULL; /* Input elements file         */
    char *ofile = NULL; /* Output elements file        */
    FILE *ifp;
    FILE *ofp;
    hdrstr ihdr;
    hdrstr ohdr;
    elmstr selm;
    elmstr elm;
    INT32 pos;
    FLOAT64 utime;
    FLOAT64 stime;
    INT32 irec;
    INT32 orec = 0;
    INT32 srec = -1;

    INT32 nrecs;
    INT32 year;
    INT32 day;
    FLOAT64 sec;

    /* Parameters for getopt() */
    extern int opterr;
    extern int optind;
    extern char *optarg;
    int c;


    /*									*/
    /* Process command-line arguments                                   */
    /*									*/
    while ((c = getopt(argc, argv, CMD_ARGS)) != EOF) {
        switch (c) {
        case 'n':
            ndays = atoi(optarg);
            break;
        default:
            usage(argv[0]);
            break;
        }
    }
    switch (argc - optind + 1) {
    case 3:
        ifile = argv[optind + 0];
        ofile = argv[optind + 1];
        break;
    default:
        usage(argv[0]);
        break;
    }

    /*									*/
    /* Open file for reading						*/
    /*									*/
    if ((ifp = fopen(ifile, "r")) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: unable to open %s for reading\n",
                __FILE__, __LINE__, ifile);
        exit(1);
    }

    /*									*/
    /* Open file for writing						*/
    /*									*/
    if ((ofp = fopen(ofile, "w")) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: unable to open %s for writing\n",
                __FILE__, __LINE__, ofile);
        exit(1);
    }

    /*									*/
    /* Read input header                                        	*/
    /*									*/
    if (fread(&ihdr, sizeof (ihdr), 1, ifp) != 1) {
        fprintf(stderr,
                "-E- %s line %d: error reading %s\n",
                __FILE__, __LINE__, ifile);
        exit(1);
    }

    /*									*/
    /* Now, we read the last record, determine the end date, and from   */
    /* that and ndays we compute the start date.                        */
    /*									*/
    nrecs = ihdr.nrecs;
    if (endianess() == 1)
        swapc_bytes((char *) &nrecs, 4, 1);
    pos = sizeof (elm)*(nrecs - 1);
    if (fseek(ifp, pos, SEEK_SET) != 0) {
        fprintf(stderr,
                "-E- %s line %d: error reading %s\n",
                __FILE__, __LINE__, ifile);
        exit(1);
    }
    if (fread(&elm, sizeof (elm), 1, ifp) != 1) {
        fprintf(stderr,
                "-E- %s line %d: error reading %s\n",
                __FILE__, __LINE__, ifile);
        exit(1);
    }

    year = elm.year;
    day = elm.day;
    sec = elm.sec;
    if (endianess() == 1) {
        swapc_bytes((char *) &year, 4, 1);
        swapc_bytes((char *) &day, 4, 1);
        swapc_bytes((char *) &sec, 8, 1);
    }

    utime = yds2unix(year, day, sec);
    stime = utime - ndays * 3600.0 * 24.0;

    /*									*/
    /* Write a dummy header, then read input and write desired output	*/
    /*									*/
    if (fwrite(&ohdr, sizeof (ohdr), 1, ofp) != 1) {
        fprintf(stderr,
                "-E- %s line %d: error writing %s\n",
                __FILE__, __LINE__, ofile);
        exit(1);
    }
    orec++;

    if (fseek(ifp, sizeof (ihdr), SEEK_SET) != 0) {
        fprintf(stderr,
                "-E- %s line %d: error reading %s\n",
                __FILE__, __LINE__, ifile);
        exit(1);
    }

    for (irec = 1; irec < nrecs; irec++) {
        if (fread(&elm, sizeof (elm), 1, ifp) != 1) {
            fprintf(stderr,
                    "-E- %s line %d: error reading %s\n",
                    __FILE__, __LINE__, ifile);
            exit(1);
        }

        year = elm.year;
        day = elm.day;
        sec = elm.sec;
        if (endianess() == 1) {
            swapc_bytes((char *) &year, 4, 1);
            swapc_bytes((char *) &day, 4, 1);
            swapc_bytes((char *) &sec, 8, 1);
        }
        utime = yds2unix(year, day, sec);
        if (utime > stime) {
            if (fwrite(&elm, sizeof (elm), 1, ofp) != 1) {
                fprintf(stderr,
                        "-E- %s line %d: error writing %s\n",
                        __FILE__, __LINE__, ofile);
                exit(1);
            }
            orec++;
            if (srec < 0) {
                srec = irec;
                memcpy(&selm, &elm, sizeof (elm));
            }
        }
    }

    /*									*/
    /* Write final header                                       	*/
    /*									*/
    nrecs = orec;
    if (endianess() == 1) {
        swapc_bytes((char *) &orec, 4, 1);
    }
    ohdr.nrecs = orec;
    if (fseek(ofp, 0, SEEK_SET) != 0) {
        fprintf(stderr,
                "-E- %s line %d: error reading %s\n",
                __FILE__, __LINE__, ifile);
        exit(1);
    }
    if (fwrite(&ohdr, sizeof (ohdr), 1, ofp) != 1) {
        fprintf(stderr,
                "-E- %s line %d: error writing %s\n",
                __FILE__, __LINE__, ofile);
        exit(1);
    }


    printf("A total of %d records written to %s\n", nrecs, ofile);
    year = selm.year;
    day = selm.day;
    sec = selm.sec;
    if (endianess() == 1) {
        swapc_bytes((char *) &year, 4, 1);
        swapc_bytes((char *) &day, 4, 1);
        swapc_bytes((char *) &sec, 8, 1);
    }
    printf("Start year, day, sec = %d, %d, %lf\n",
            year, day, sec);
    year = elm.year;
    day = elm.day;
    sec = elm.sec;
    if (endianess() == 1) {
        swapc_bytes((char *) &year, 4, 1);
        swapc_bytes((char *) &day, 4, 1);
        swapc_bytes((char *) &sec, 8, 1);
    }
    printf("End   year, day, sec = %d, %d, %lf\n",
            year, day, sec);

    exit(0);

}
