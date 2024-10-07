/*----------------------------------------------------------------------
Copyright (C) 2000,  Space Science and Engineering Center, University
of Wisconsin-Madison, Madison WI.
      
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <PGS_TD.h>

#define buffsize_mb 0.5
#define primary_hdr_size 6
#define OUTBUFLEN 384

#define O_STRING "nf:"

#define USAGE "write_constructor_file (%s %s)\n\
\
Usage: %s [-n] [-f outfile] infile\n\
\n\
Note:\n\
If outfile is not specified, it will be infile.constr\n\
\n\
-n : do not generate a constructor file but print statistics\n\
"

void print_usage(char *prog) {
    printf(USAGE, __DATE__, __TIME__, prog);
}

int main(int argc, char *argv[]) {
    int n_packets = 0;
    int read_cnt;
    long strm_pos = 0;
    int packet_length = 0;
    int buf_pos = 0;
    int buf_pos_last = 0;
    int bytes_left;
    int last_pkt_in_file;
    int last_pkt_in_buffer;
    int found_start_time;
    int buffsize;
    int n_night = 0;
    int n_day = 0;
    int c;
    int fopt = 0;
    int nopt = 0;
    int errflag = 0;

    static int pkt_len_off = 4;
    static int time_off = primary_hdr_size;
    static int time_cnt = 8;
    static int L0_true = 0;
    static int L0_false = 1;
    static int day_pkt_size = 642;
    static int night_pkt_size = 276;
    static unsigned char start_time_tag[8];
    static unsigned char stop_time_tag[8];
    char *cptr = (char *) &n_packets;

    double taitime_start;
    double taitime_stop;

    unsigned char *buffer, outbuf[OUTBUFLEN];

    char outfile[FILENAME_MAX];
    char *infile;
    int len;

    FILE *stream, *outfp = NULL;

    /* For getopt() */
    extern int optind;
    extern char *optarg;


    /*
     ** Process command-line options
     */
    while ((c = getopt(argc, argv, O_STRING)) != -1) {
        switch (c) {
        case 'n':
            nopt++;
            break;
        case 'f':
            len = snprintf(outfile, sizeof (outfile), "%s", optarg);
            if (strlen(outfile) == len) {
                fopt++;
            } else {
                printf("specified output file name too long!\n");
                errflag++;
            }
            break;
        default:
            errflag++;
            break;
            ;
        }
    }


    if (optind < argc && !errflag) {
        printf("write_constructor_file: Version as of %s %s\n", __DATE__, __TIME__);

        /* Set input file */
        infile = argv[argc - 1];

        /* If -n not given and no output file specified, create its name from input file */
        if (!fopt && !nopt) {
            len = snprintf(outfile, FILENAME_MAX, "%s.constr", infile);
            if (strlen(outfile) != len) {
                printf("derived output file name is too long!\n");
                errflag++;
            }
        }

        /* Open input file */
        stream = fopen(infile, "r");

        if (stream != NULL) {
            /* Open output file if -n not given */
            if (!nopt)
                outfp = fopen(outfile, "w");

            if (outfp != NULL || nopt) {

                fseek(stream, (long) 0, SEEK_SET);

                buffsize = buffsize_mb * 1000000;
                buffer = (unsigned char *) malloc(buffsize * sizeof (unsigned char));
                if (buffer == NULL) {
                    fseek(stream, (long) 0, SEEK_SET);
                    return 1;
                }

                last_pkt_in_file = L0_false;
                found_start_time = L0_false;

                while (last_pkt_in_file == L0_false) {
                    last_pkt_in_buffer = L0_false;

                    read_cnt = fread(buffer, sizeof (char), buffsize, stream);
                    //printf("read_cnt: %d  buffsize: %d\n", read_cnt,buffsize);

                    buf_pos = 0;
                    while (last_pkt_in_buffer == L0_false) {
                        bytes_left = read_cnt - buf_pos;
                        if (bytes_left < day_pkt_size) {
                            if (bytes_left == 0) {
                                last_pkt_in_buffer = L0_true;
                                continue;
                            } else if ((bytes_left == night_pkt_size) ||
                                    (bytes_left == 2 * night_pkt_size)) {
                            } else if (feof(stream) == 0) {
                                last_pkt_in_buffer = L0_true;
                                fseek(stream, strm_pos, SEEK_SET);
                                continue;
                            } else {
                                printf("Corrupt record found at end of file.\n");
                                free(buffer);
                                fseek(stream, (long) 0, SEEK_SET);
                                return 3;
                            }
                        }

                        packet_length = buffer[buf_pos + (pkt_len_off)]*256 +
                                buffer[buf_pos + pkt_len_off + 1];
                        packet_length += 1;
                        packet_length += primary_hdr_size;

                        if ((packet_length != 642) && (packet_length != 276)) {
                            //        utc_str[27] = '\0';
                            if (n_packets > 0) {
                                printf("Packet with invalid length found: (%d).\n", packet_length);
                            }
                            free(buffer);
                            fseek(stream, (long) 0, SEEK_SET);
                            return 3;
                        }

                        if (packet_length == 642) n_day++;
                        if (packet_length == 276) n_night++;

                        if (found_start_time == L0_false) {
                            memcpy(start_time_tag, &buffer[buf_pos + time_off], time_cnt);
                            found_start_time = L0_true;
                        } else {
                            memcpy(stop_time_tag, &buffer[buf_pos_last + time_off], time_cnt);
                            PGS_TD_EOSAMtoTAI(stop_time_tag, &taitime_stop);
                            PGS_TD_TAItoUTC(taitime_stop, (char*) outbuf);
                            /*printf("stoptime =%s\n", outbuf);*/
                        }

                        strm_pos += (long) packet_length;
                        buf_pos_last = buf_pos;
                        buf_pos += packet_length;
                        n_packets++;
                        if ((n_packets % 100000) == 0) {
                            printf("%10d packets read\n", n_packets);
                        }
                    } /* end while ( last_pkt_in_buffer == L0_false ) */

                    if (feof(stream) != 0) {
                        last_pkt_in_file = L0_true;
                    }

                    if (ferror(stream) != 0) {
                        free(buffer);
                        fseek(stream, (long) 0, SEEK_SET);
                        return 2;
                    }

                } /* end while ( last_pkt_in_file == L0_false ) */

                free(buffer);

                printf("%d Total packets read\n", n_packets);

                printf("%d day packets\n", n_day);
                printf("%d night packets\n", n_night);


                /* If -n option not given, write constructor record */
                if (outfp != NULL) {

                    printf("Writing constuctor record for %s to %s\n", infile, outfile);

                    memset(outbuf, 0, OUTBUFLEN);

                    fseek(outfp, 0, SEEK_SET);

                    memset(&outbuf[0x33], 1, 1);

                    memcpy(&outbuf[0x50], start_time_tag, 8);
                    memcpy(&outbuf[0x58], stop_time_tag, 8);

                    memcpy(&outbuf[0x16c], start_time_tag, 8);
                    memcpy(&outbuf[0x174], stop_time_tag, 8);

                    memcpy(&outbuf[0x74], cptr + 3, 1);
                    memcpy(&outbuf[0x75], cptr + 2, 1);
                    memcpy(&outbuf[0x76], cptr + 1, 1);
                    memcpy(&outbuf[0x77], cptr + 0, 1);

                    memset(&outbuf[0x93], 1, 1);
                    memset(&outbuf[0xa3], 1, 1);
                    memset(&outbuf[0xf7], 2, 1);

                    memset(&outbuf[0x167], 1, 1);

                    if (fwrite(outbuf, OUTBUFLEN, 1, outfp) == 0) {
                        printf("fwrite failed: %s\n", strerror(errno));
                        errflag++;
                    }

                    /* Close output file */
                    fclose(outfp);

                }

                /* Close input file */
                fclose(stream);

                PGS_TD_EOSAMtoTAI(start_time_tag, &taitime_start);
                PGS_TD_TAItoUTC(taitime_start, (char*) outbuf);
                printf("starttime=%s\n", outbuf);

                PGS_TD_EOSAMtoTAI(stop_time_tag, &taitime_stop);
                PGS_TD_TAItoUTC(taitime_stop, (char*) outbuf);
                printf("stoptime =%s\n", outbuf);

                printf("granule length =%f\n", taitime_stop - taitime_start);
            } else {
                /* Failed to open output file */
                printf("failed to open output file: %s\n", strerror(errno));
                errflag++;
            }
        } else {
            /* Failed to open input file */
            printf("failed to open input file: %s\n", strerror(errno));
            errflag++;
        }
    } else {
        print_usage(argv[0]);
    }

    return errflag;
}
