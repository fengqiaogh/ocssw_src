/* ============================================================================	*/
/*                                                                           	*/
/* raw2L0 - converts raw (10-bit, 13860 byte) SeaWifs frames to Level-0      	*/
/*                                                                           	*/
/* Synopsis:									*/
/*										*/
/*	raw2L0 -f input-raw-filename -o output-L0-filename -l reclen		*/
/*	raw2L0 < input-raw-filename > output-L0-filename -l reclen		*/
/*										*/
/* Description:                                                              	*/
/*                                                                           	*/
/*     The program reads raw SeaWiFS minor frames in 10-bit, 13860-byte	format,	*/
/*     from standard input, and writes the corresponding L0 file to standard	*/
/*     output.  Thus, it can be used as a filter. It will handle HRPT frames	*/
/*     as well, provided they are truncated to 13860 bytes from there original	*/
/*     13862.5-byte form (or, better yet, the 20-bit pad is removed).		*/
/*                                                                           	*/
/*     Code now has option to define record length. Try 13864 for some HRPT.    */
/*										*/
/* See Also:									*/
/*										*/
/*	neb2raw									*/
/*										*/
/*										*/
/* Written By:                                                               	*/
/*                                                                           	*/
/*     Bryan A. Franz 								*/
/*     General Sciences Corp.                                                	*/
/*     4 June 1996                                                           	*/
/*										*/
/*     Original version was blowup.f by Fred Patt				*/
/*                                                                           	*/
/* ============================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "genutils.h"

typedef unsigned char BYTE;
typedef short int INT16;
typedef int32_t INT32;

#define CMD_ARGS        "f:o:l:s:"                /* Valid commandline options    */

#define MAXRECLEN       13864
#define INRECLEN 13860 /* # Bytes in raw minor frame			*/
#define OUTRECLEN 21504 /* # Bytes in level-0 record			*/
#define OUTHDRLEN 512 /* # Bites in level-0 header			*/


void bit10toi2(BYTE bit10[], INT16 ints[]);

int main(int argc, char* argv[]) {
    BYTE inrec [MAXRECLEN];
    BYTE outrec[OUTRECLEN];
    BYTE outhdr[OUTHDRLEN];
    INT16 intbuf[4];

    FILE *infp = stdin;
    FILE *outfp = stdout;

    INT32 inRecNum = 0L;
    INT32 outRecNum = 0L;
    INT32 i;

    INT32 reclen = INRECLEN;
    INT32 startrec = 1;

    extern int opterr; /* used by getopt()     */
    extern int optind; /* used by getopt()     */
    extern char *optarg; /* used by getopt()     */
    int c; /* used by getopt()     */

    /*                                                                          */
    /* Process command-line inputs                                              */
    /*                                                                          */
    while ((c = getopt(argc, argv, CMD_ARGS)) != EOF) {
        switch (c) {
        case 'f':
            if ((infp = fopen(optarg, "r")) == NULL) {
                fprintf(stderr, "%s: Error opening %s for reading.\n", argv[0], optarg);
                exit(1);
            }
            break;
        case 'o':
            if ((outfp = fopen(optarg, "w")) == NULL) {
                fprintf(stderr, "%s: Error opening %s for writing.\n", argv[0], optarg);
                exit(1);
            }
            break;
        case 'l':
            reclen = atoi(optarg);
            break;
        case 's':
            startrec = atoi(optarg);
            break;
        default:
            fprintf(stderr, "Usage: %s [-f input-nebula-filename]\n", argv[0]);
            fprintf(stderr, "          [-o output-raw-filename]\n");
            fprintf(stderr, "          [-l reclen (13860 or 13864)]\n");
            fprintf(stderr, "          [-s start rec number (>1)]\n");
            break;
        }
    }

    /*										*/
    /* Initialize								*/
    /*										*/
    memset(outrec, 0, OUTRECLEN);
    memset(outhdr, 0, OUTHDRLEN);

    /*										*/
    /* Write L0 file header to standard out					*/
    /*										*/
    memcpy(outhdr, "CWIF", 5);
    if (fwrite(outhdr, OUTHDRLEN, 1, outfp) != 1) {
        printf("Error writing output file at record %d\n", outRecNum);
        exit(1);
    }


    /*										*/
    /* Read through end of input file						*/
    /*										*/
    while (fread(inrec, reclen, 1, infp) == 1) {

        inRecNum++;

        if (inRecNum < startrec) continue;

        /*									*/
        /* Copy S/C ID								*/
        /*									*/
        bit10toi2(&inrec[5], intbuf);
        memcpy(&outrec[3], &intbuf[2], 4);

        /*									*/
        /* Copy Timetag								*/
        /*									*/
        bit10toi2(&inrec[10], intbuf);
        memcpy(&outrec[7], intbuf, 8);

        /*									*/
        /* Copy SOH block (no conversion required)				*/
        /*									*/
        memcpy(&outrec[15], &inrec[15], 775);

        /*									*/
        /* Copy rest of frame (Inst., Scan, G&TDI)				*/
        /*									*/
        for (i = 0; i < 2589; i++) {
            bit10toi2(&inrec[790 + i * 5], intbuf);
            memcpy(&outrec[790 + i * 8], intbuf, 8);
        }

        /*									*/
        /* Write L0 record to standard output					*/
        /*									*/
        if (fwrite(outrec, OUTRECLEN, 1, outfp) != 1) {
            printf("Error writing output file at record %d\n", outRecNum);
            exit(1);
        } else
            outRecNum++;

    }


    /*										*/
    /* Check for error on input stream						*/
    /*										*/
    if (ferror(infp)) {
        fprintf(stderr, "Error reading input file at record %d\n", inRecNum);
        fclose(infp);
        exit(1);
    }

    /*										*/
    /* Normal Termination							*/
    /*										*/
    fprintf(stderr, "End of file reached.\n");
    fprintf(stderr, "Number of  input records = %d\n", inRecNum);
    fprintf(stderr, "Number of output records = %d\n", outRecNum);
    fclose(infp);
    fclose(outfp);

    return (0);

}


/* ---------------------------------------------------------------------------- */
/* Unpacks a 5-element byte array into a 4-element I*2 array, using the 10 	*/
/* least significant bits of each integer.                             		*/

/* ---------------------------------------------------------------------------- */
void bit10toi2(BYTE bit10[], INT16 ints[]) {
    ints[0] = (((INT16) bit10[1] & 0xC0) >> 6) + (((INT16) bit10[0] & 0xFF) << 2);
    ints[1] = (((INT16) bit10[2] & 0xF0) >> 4) + (((INT16) bit10[1] & 0x3F) << 4);
    ints[2] = (((INT16) bit10[3] & 0xFC) >> 2) + (((INT16) bit10[2] & 0x0F) << 6);
    ints[3] = (((INT16) bit10[4] & 0xFF) >> 0) + (((INT16) bit10[3] & 0x03) << 8);

    if (endianess() == 1) {
        swapc_bytes((char *) ints, 2, 4);
    }
}
