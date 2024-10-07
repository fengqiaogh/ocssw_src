/* ============================================================================ */
/*                                                                           	*/
/* neb2raw - converts OSC nebula files to native S/C frames (13860-byte,10-bit)	*/
/*                                                                           	*/
/* Synopsis:                                                                   	*/
/*                                                                           	*/
/*	neb2raw [ -f nebula-file ] [ -o raw-file ] [-b backorbit-file ]		*/
/*										*/
/*	nebula-file: 	name of input nebula file to read			*/
/*			(default is to read from standard input, see warning)  	*/
/*	raw-file:	name of output file for raw S/C frames			*/
/*			(default is to write to standard output)		*/
/*	backorbit-file:	name of output file for back-orbit SOH messages 	*/
/*			(default is to discard back-orbit records)		*/
/*                                                                           	*/
/* Description:									*/
/* 	 	                                                             	*/
/* 	A SeaWiFS data file in OSC's Nebula format is provided to the program	*/
/*	via standard input or the -f option flag.  The file is searched, record	*/
/*	by record, to find any of three types of data records:			*/
/* 	 	                                                             	*/
/*	1) HRPT records, which are stored as single nebula records consisting	*/
/*	   of: the 40-byte nebula header and 13862-bytes of the S/C frame.  An	*/
/*	   HRPT frame is actually 13862.5 bytes, but the last half-byte gets	*/
/*	   truncated by the OSC front-end.  There is no useful data longward of */
/*	   byte 13735 (the 20-bit pad), so this trunctaion doesn't matter. This	*/
/*	   program extracts the first 13860 bytes of the HRPT record, and 	*/
/*	   writes it either to standard output, or to a named file provided via	*/
/*	   the -o option.							*/
/*										*/
/*	2) Stored GAC/LAC records, which are divided over four 3540-byte nebula */
/*	   records. The first record in the sequence contains a 40-byte nebula	*/
/*	   header, an 8-byte HDLC header, a 2-byte length field, and 3490 bytes */
/*	   of a recorded S/C frame.  The next two records contain a 40-byte neb */
/*	   header, an 8-byte HDLC header, and 3492 bytes of the S/C frame.  The */
/*	   final record in the sequence has the 40-byte neb header, the 8-byte  */
/*	   HDLC header, the remaining 3386 bytes of the S/C frame, and 103 	*/
/*	   bytes of padding.  The program checks for proper sequence order, but	*/
/*	   it can not distinguish between the middle two records.  The raw 	*/
/*	   frames are reconstructed from the nebula record and written to 	*/
/*	   standard output, or to a named file provided with the -o option. 	*/
/*	   Out-of-order sequences, which are most likely associated with a 	*/
/*	   corrupted HDLC header, are discarded. 				*/
/* 	 	                                                             	*/
/*	3) Logger (back-orbit) records, which are stored in 3540-byte nebula	*/
/*	   records. Each record consists of a 40-byte neb header and an 8-byte 	*/
/*	   HDLC header, followed by a series of variable-length spacecraft 	*/
/* 	   messages originating	from the GIM or SCM memory subsystem.  Each of	*/
/*	   these messages are just wrappers around another spacecraft message.	*/
/*	   The program reads the 3540-byte nebula records and writes a 3502	*/
/*	   byte record containing the HDLC header and the 3492 bytes that     	*/
/*	   follow it. The last two bytes are a placeholder to account for the 	*/
/*	   CRC bytes that are stored in the SeaSpace version of the back-orbit	*/
/*  	   file, but that are not available from the nebula input file.	 	*/
/*  									 	*/
/*	Generally, nebula files from OSC will contain either type 1 records, or */
/*	type 2 and 3 records, but any combination can be handled.  Nebula 	*/
/*	records which do not meet one of these three types are probably just	*/
/*	status messages from OSC's high-speed front-end.			*/
/* 	 	                                                             	*/
/*										*/
/* Warnings:									*/
/*										*/
/*	Due to the random-seek capability required to step through the nebula	*/
/*	records, the program can not read from a pipe.				*/
/*										*/
/*										*/
/* See Also:									*/
/*										*/
/* 	raw2L0									*/
/*                                                                           	*/
/*                                                                           	*/
/* Written By:                                                               	*/
/*                                                                           	*/
/*     Bryan A. Franz                                                        	*/
/*     General Sciences Corp.                                                	*/
/*     4 June 1996                                                           	*/
/*                                                                           	*/
/* ============================================================================ */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "loghdr.h"


#define CMD_ARGS        "f:o:b:"         /* Valid commandline options	*/
#define  HRPTRECLEN 13902   /* Length of HRPT nebula record */
#define  STORRECLEN 3540   /* Length of recorded data rec  */
#define  FRAMELEN 13860   /* Bytes in spacecraft mnf	*/
#define SIZEOFBUF       (16*10124)  /* Data buffer length		*/
#define SIZEOFHDR sizeof(struct LOGHDR) /* Nebula header length		*/
#define SIZEOFHDLC 8   /* HDLC header length		*/
#define GACID  80   /* GAC ID in HDLC header	*/
#define LACID  96   /* LAC ID in HDLC header	*/
#define BRBID  64   /* Back-orbit ID in HDLC header	*/
#define IDBYTE  3   /* ID byte location HDLC header	*/
#define SQBYTE  5   /* Sequence # location HDLC hdr	*/



int issoh(int addr);
int extract_loggermsgs(FILE *outfp, BYTE *buf, int buflen);



/* ---------------------------------------------------------------------------- */
/* main										*/

/* ---------------------------------------------------------------------------- */
int main(int argc, char* argv[]) {
    FILE *infp = stdin; /* Input file pointer	*/
    FILE *outfp = stdout; /* Output file pointer	*/
    FILE *brbfp = NULL; /* Back-orbit file pntr	*/

    struct LOGHDR hdrbuf; /* Nebula header buffer	*/
    struct LOGHDR *hdr = &hdrbuf; /* Pointer to neb hdr 	*/
    unsigned char hdlc[SIZEOFHDLC]; /* HDLC header buffer	*/
    unsigned char buf[SIZEOFBUF]; /* Data buffer		*/

    int32_t inRecNum = 0L; /* # Input neb records	*/
    int32_t outRecNum = 0L; /* # Output mnf records */
    int32_t brbRecNum = 0L; /* # Output brb records */

    int32_t inPos = 0L; /* Byte pos in neb file */
    int swap = 0; /* Byte swapping flag	*/
    int first = 1; /* First neb rec flag	*/
    int partnum = 0; /* Frame sequence flag	*/

    extern int opterr; /* used by getopt()     */
    extern int optind; /* used by getopt()     */
    extern char *optarg; /* used by getopt()     */
    int c; /* used by getopt()     */

    /*										*/
    /* Process command-line inputs						*/
    /*										*/
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
        case 'b':
            if ((brbfp = fopen(optarg, "w")) == NULL) {
                fprintf(stderr, "%s: Error opening %s for writing.\n", argv[0], optarg);
                exit(1);
            }
            break;
        default:
            fprintf(stderr, "Usage: %s [-f input-nebula-filename]\n", argv[0]);
            fprintf(stderr, "          [-o output-raw-filename]\n");
            fprintf(stderr, "          [-b backorbit-filename]\n");
            break;
        }
    }


    /*										*/
    /* Read through end of input file						*/
    /*										*/
    while (fread(hdr, SIZEOFHDR, 1, infp) == 1) {

        /*									*/
        /* If this is the first Nebula record processed, determine swapping for */
        /* Nebula macros.							*/
        /*									*/
        if (first) {
            first = 0;
            switch (hdr->headerVersion) {
            case 0x0100:
                swap = 1;
                break;
            case 0x0001:
                swap = 0;
                break;
            default:
                swap = 0;
                break;
            }
        }

        /*									*/
        /* Do some quality checking on the Nebula header			*/
        /*									*/
        if (HDR_WORD(hdr->headerVersion) != 1) {
            fprintf(stderr, "Bad nebula header at record %d\n", inRecNum);
            fclose(infp);
            exit(5);
        }

        /*									*/
        /* If this is the record we want, read the data				*/
        /*									*/
        if (HDR_LONG(hdr->blockSize) == HRPTRECLEN) {

            /*									*/
            /* HRPT Record, just copy data from input to output			*/
            /*									*/
            if (fread(buf, FRAMELEN, 1, infp) != 1) {
                fprintf(stderr, "Error reading input file at record %d\n", inRecNum);
                exit(1);
            }
            if (fwrite(buf, FRAMELEN, 1, outfp) != 1) {
                fprintf(stderr, "Error writing output file at record %d\n", outRecNum);
                exit(1);
            } else
                outRecNum++;

        } else if (HDR_LONG(hdr->blockSize) == STORRECLEN) {

            /*									*/
            /* Possibly part of a stored GAC or LAC record, or it maybe		*/
            /* some back-orbit telemetry messages. Examine the HDLC header.	*/
            /*									*/
            if (fread(hdlc, SIZEOFHDLC, 1, infp) != 1) {
                fprintf(stderr, "Error reading input file at record %d\n", inRecNum);
                exit(1);
            }
            if (hdlc[IDBYTE] == GACID || hdlc[IDBYTE] == LACID) {

                /*								*/
                /* Part of stored LAC or GAC, but which part.			*/
                /*								*/
                if (partnum == 0 && hdlc[SQBYTE] == 2) {

                    /*								*/
                    /* Beginning of sequence, read into beginning of buffer.	*/
                    /*								*/
                    partnum++;
                    fseek(infp, 2, SEEK_CUR); /* Skip msg length field */
                    if (fread(&buf[0], 1, 3490, infp) != 3490) {
                        fprintf(stderr, "Error reading input file at record %d\n", inRecNum);
                        exit(1);
                    }

                } else if (partnum == 1 && hdlc[SQBYTE] == 0) {

                    /*								*/
                    /* Middle segment (assumed to be second)			*/
                    /*								*/
                    partnum++;
                    if (fread(&buf[3490], 1, 3492, infp) != 3492) {
                        fprintf(stderr, "Error reading input file at record %d\n", inRecNum);
                        exit(1);
                    }

                } else if (partnum == 2 && hdlc[SQBYTE] == 0) {

                    /*								*/
                    /* Middle segment (assumed to be third)			*/
                    /*								*/
                    partnum++;
                    if (fread(&buf[6982], 1, 3492, infp) != 3492) {
                        fprintf(stderr, "Error reading input file at record %d\n", inRecNum);
                        exit(1);
                    }

                } else if (hdlc[SQBYTE] == 1) {

                    /*								*/
                    /* Last segment, write output record if sequence O.K.	*/
                    /*								*/
                    partnum++;
                    if (partnum != 4) {
                        fprintf(stderr, "Sequence error at record %d\n", inRecNum);
                    } else {
                        if (fread(&buf[10474], 1, 3386, infp) != 3386) {
                            fprintf(stderr, "Error reading input file at record %d\n", inRecNum);
                            exit(1);
                        }
                        if (fwrite(buf, 1, 13860, outfp) != 13860) {
                            printf("Error writing output file at record %d\n", outRecNum);
                            exit(1);
                        } else
                            outRecNum++;
                    }
                    partnum = 0;

                } else
                    fprintf(stderr, "Sequence error at record %d\n", inRecNum);

            } else if (hdlc[IDBYTE] == BRBID && brbfp != NULL) {

                if (fwrite(hdlc, SIZEOFHDLC, 1, brbfp) != 1) {
                    fprintf(stderr, "Error writing backorbit file at record %d\n", inRecNum);
                    exit(1);
                }

                /*								*/
                /* Back-orbit telemetry. Copy data to buffer and write to the   */
                /* back-orbit file.		                         	*/
                /*								*/
                if (fread(buf, 1, 3494, infp) != 3494) {
                    printf("Error reading input file at record %d\n", inRecNum);
                    exit(1);
                }
                if (fwrite(buf, 1, 3494, brbfp) != 3494) {
                    fprintf(stderr, "Error writing backorbit file at record %d\n", inRecNum);
                    exit(1);
                }

                brbRecNum++;

            }
        }

        /*									*/
        /* Update input file position to start of next Nebula record		*/
        /*								 	*/
        inRecNum++;
        inPos += HDR_LONG(hdr->blockSize);
        fseek(infp, inPos, SEEK_SET);

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
    /* Normal Termination.							*/
    /*										*/
    fprintf(stderr, "End of file reached.\n");
    fprintf(stderr, "Number of  input records      = %d\n", inRecNum);
    fprintf(stderr, "Number of output records      = %d\n", outRecNum);
    fprintf(stderr, "Number of back-orbit records  = %d\n", brbRecNum);
    //    if (infp  != NULL && infp  != stdin ) fclose(infp );
    //    if (outfp != NULL && outfp != stdout) fclose(outfp);
    //    if (brbfp != NULL) fclose(brbfp);

    return (0);

}


/* ---------------------------------------------------------------------------- */
/* extract_loggermsgs() - write any SOH messages found in the input buffer to	*/
/*			  the output file. Returns number of messages written.  */
/*										*/
/* The input buffer is assummed to contain a stream of SeaStar memory logger	*/
/* messages from the GIM or SCM memory modules.  Each message functions as a	*/
/* wrapper around another message, and it is some of those innermost messages	*/
/* which this function extracts.  Specifically, it looks for SOH subsytem msgs	*/
/* and writes them to the output file.						*/
/*										*/
/* Written By: BA Franz, GSC, June 1996						*/
/*										*/

/* ---------------------------------------------------------------------------- */
int extract_loggermsgs(FILE *outfp, /* Output file pointer			*/
        BYTE *buf, /* Buffer holding stream of GIM/SCM msg	*/
        int buflen) /* Buffer length			*/ {
    int i = 0;
    short int addr = 0; /* Spacecraft message source address		*/
    short int len = 1; /* GIM/SCM S/C message length			*/
    short int newlen = 0; /* SOH-subsystem S/C message length		*/
    int cnt = 0; /* Number of SOH messages found			*/

    /*										*/
    /* Step-through GIM/SCM messages, looking inside for SOH subsystem messages */
    /*										*/
    while (i < buflen && len > 0) {
        memcpy(&len, &buf[i], 2);
        if (len > 0) {
            memcpy(&addr, &buf[i + 24], 2);
            if (issoh(addr)) {

                /*                                                             */
                /* Need to adjust length field in s/c msg header.  The logger  */
                /* includes the trailing four of the GIM/SCM wrapper here, but */
                /* we are extracting the message from that wrapper so that it  */
                /* will look no different than a recorded SOH message.         */
                /*                                                             */
                memcpy(&newlen, &buf[i + 10], 2);
                newlen -= 4;
                memcpy(&buf[i + 10], &newlen, 2);

                /*                                                             */
                /* Write S/C message to output file and count it               */
                /*                                                             */
                if (fwrite(&buf[i + 10], 1, len - 14, outfp) != len - 14)
                    return ( -1);
                else
                    cnt++;
            }
        }
        i += len;
    }

    return (cnt);
}



/* ---------------------------------------------------------------------------- */
/* issoh() - returns true if the input address is associated with an SOH subsys	*/

/* ---------------------------------------------------------------------------- */
int issoh(int addr) {
    switch (addr) {
    case 4416: /* SGA	 */
    case 4384: /* SAC	 */
    case 4368: /* SAA	 */
    case 9217: /* RXS-1 */
    case 9218: /* RXS-2 */
    case 9985: /* TXL	 */
    case 9729: /* TXS	 */
    case 4609: /* SMU	 */
    case 8705: /* EPS	 */
    case 1281: /* FDR-1 */
    case 1282: /* FDR-2 */
    case 8544: /* GUA	 */
    case 8464: /* GSC	 */
    case 272: /* PSA	 */
    case 320: /* PFA	 */
    case 4544: /* SHM	 */
    case 8640: /* GHM	 */
        return (1);
    default:
        return (0);
    }
}




