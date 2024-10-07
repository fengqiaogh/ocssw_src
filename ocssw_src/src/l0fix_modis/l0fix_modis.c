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
#include <PGS_TD.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define buffsize_mb 0.5
#define primary_hdr_size 6

int main(int argc, char *argv[]) {
    int n_packets = 0;
    int read_cnt;
    long strm_pos = 0;
    int packet_length = 0;
    int buf_pos = 0;
    int bytes_left;
    int last_pkt_in_file;
    int last_pkt_in_buffer;
    int found_start_time;
    int buffsize;
    int write;
    int n_night = 0;
    int n_day = 0;
    struct stat filestat;

    static int pkt_len_off = 4;
    static int time_off = primary_hdr_size;
    static int time_cnt = 8;
    static int L0_true = 0;
    static int L0_false = 1;
    static int day_pkt_size = 642;
    static int night_pkt_size = 276;
    static unsigned char time_tag[8];

    double taitime_start;
    double taitime_stop;
    double taitime_tag;
    double taitime_good;
    double eps = 1e-6;
    size_t status;

    unsigned char *buffer;
    FILE *stream, *stream_w = NULL;

    /*
      1) Get start and stop times

      ./fix_L0 L0file -1 -1

      2) fix L0

      ./fix_L0 L0file starttime stoptime outputL0
     */

    printf("fix_L0: Version as of 07/16/07\n\n");

    stream = fopen(argv[1], "r");
    if (stream == NULL) {
        printf("%s not found.\n", argv[1]);
        exit(-1);
    }
    fseek(stream, (long) 0, SEEK_SET);

    taitime_start = atof(argv[2]);
    taitime_stop = atof(argv[3]);

    if (taitime_stop == -1) {
        taitime_stop = 1e30;
    } else {
        stream_w = fopen(argv[4], "wb");
    }

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
                    free(buffer);

                    printf("%d Total packets read\n", n_packets);
                    printf("taitime_stop: %f\n", taitime_good);

                    if (stream_w != NULL) fclose(stream_w);
                    fclose(stream);

                    if (atof(argv[2]) != -1) {
                        stat(argv[4], &filestat);
                        if (n_day * 642 + n_night * 276 != filestat.st_size) {
                            printf("output filesize discrepency.\n");
                            printf("Computed filesize: %d\n", n_day * 642 + n_night * 276);
                            printf("Actual filesize:   %d\n", (int) filestat.st_size);
                            return 3;
                        }
                    }

                    return 0;
                }
            }

            packet_length = buffer[buf_pos + (pkt_len_off)]*256 +
                    buffer[buf_pos + pkt_len_off + 1];
            packet_length += 1;
            packet_length += primary_hdr_size;

            write = 1;
            if ((packet_length != 642) && (packet_length != 276)) {
                if (n_packets > 0) {
                    printf("Packet with invalid length found %d %d\n",
                            n_packets, packet_length);
                    write = 0;
                }
            }


            memcpy(time_tag, &buffer[buf_pos + time_off], time_cnt);
            PGS_TD_EOSAMtoTAI(time_tag, &taitime_tag);

            //      if (n_packets == 300000) printf("%f\n", taitime_tag);

            if (taitime_tag < taitime_start - eps || taitime_tag > taitime_stop + eps) {
                printf("Bad L0 packet: %d %d\n",
                        taitime_tag < taitime_start,
                        taitime_tag > taitime_stop);
                printf("Bad L0 packet: %f %f %f\n",
                        taitime_start, taitime_tag, taitime_stop);
                write = 0;
            }

            if (taitime_start == -1 && found_start_time == L0_false) {
                printf("taitime_start: %f\n", taitime_tag);
                found_start_time = L0_true;
            }

            if (write && stream_w != NULL) {
                status = fwrite(&buffer[buf_pos], packet_length, 1, stream_w);
                if (status != 1) {
                    exit(4);
                }
            }

            if (write)
                taitime_good = taitime_tag;

            strm_pos += (long) packet_length;
            buf_pos += packet_length;
            n_packets++;

            if (write) {
                if (packet_length == 642) n_day++;
                if (packet_length == 276) n_night++;
            }

            if ((n_packets % 100000) == 0) {
                printf("%10d packets read\n", n_packets);
            }

        }

        if (feof(stream) != 0) {
            last_pkt_in_file = L0_true;
        }

        if (ferror(stream) != 0) {
            free(buffer);
            fseek(stream, (long) 0, SEEK_SET);
            return 2;
        }
    }
    free(buffer);

    printf("%d Total packets read\n", n_packets);
    printf("taitime_stop: %f\n", taitime_good);

    if (stream_w != NULL) fclose(stream_w);
    fclose(stream);

    if (atof(argv[2]) != -1) {
        stat(argv[4], &filestat);
        if (n_day * 642 + n_night * 276 != filestat.st_size) {
            printf("output filesize discrepency.\n");
            printf("Computed filesize: %d\n", n_day * 642 + n_night * 276);
            printf("Actual filesize:   %d\n", (int) filestat.st_size);
            return 3;
        }
    }

    return 0;
}
