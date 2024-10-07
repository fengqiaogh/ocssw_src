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

#define buffsize_mb 0.5
#define primary_hdr_size 6

#define basename(s) (strrchr((s), '/') == NULL ? (s) : strrchr((s), '/') + 1)

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
    int packet_written = 0;
    int write;

    static int pkt_len_off = 4;
    static int time_off = primary_hdr_size;
    static int time_cnt = 8;
    static int L0_true = 0;
    static int L0_false = 1;
    static int day_pkt_size = 642;
    static int night_pkt_size = 276;
    static unsigned char start_time_tag[8];
    static unsigned char stop_time_tag[8];
    static unsigned char time_tag[8];
    char *cptr3;

    double taitime_stop;
    double taitime0;
    double taitime;
    double gran_length = 300;

    char utc_str[28], utc_str2[28];
    unsigned char *buffer, outname[384];
    FILE *stream, *outfp;



    printf("%s %s (%s %s)\n\n", argv[0], "1.x", __DATE__, __TIME__);

    if (argc == 1) {
        printf("USAGE: %s MODIS_L0_PDS_file [granule_length] [zulu_start_time]\n\n", argv[0]);
        printf("granule_length: Time in seconds of each generated granule. (default = 300 sec)\n\n");
        printf("zulu_start_time: Time of the first packet included in the granule(s). All packets\n");
        printf("                 before this time will be ignored. All packets after this time\n");
        printf("                 will be included in the output L0 granule(s). If zulu_start_time\n");
        printf("                 is not set, the first output granule's filename will be set to the\n");
        printf("                 rounded 5 minute interval preceding the time of the first packet.\n");
        printf("                 zulu_start_time example: 2006-06-12T16:50:00.00000\n");
        printf("                 (default = start time of first packet in L0 granule)\n\n");

        exit(1);
    }

    stream = fopen(argv[1], "r");
    if (stream == NULL) {
        printf("%s not found.\n", argv[1]);
        exit(1);
    }
    fseek(stream, (long) 0, SEEK_SET);

    buffsize = buffsize_mb * 1000000;
    buffer = (unsigned char *) malloc(buffsize * sizeof (unsigned char));
    if (buffer == NULL) {
        fseek(stream, (long) 0, SEEK_SET);
        return -1;
    }

    if (argc >= 3) gran_length = atof(argv[2]);
    printf("Granule Length: %f\n", gran_length);

    read_cnt = fread(buffer, sizeof (char), time_off + time_cnt, stream);
    PGS_TD_EOSAMtoTAI(&buffer[time_off], &taitime0);
    PGS_TD_TAItoUTC(taitime0, utc_str2);
    strcpy(&utc_str2[15], "0:00.000000Z");
    fseek(stream, (long) 0, SEEK_SET);

    if (argc == 4) strcpy(utc_str2, argv[3]);
    printf("UTC start time: %s\n", utc_str2);

    PGS_TD_UTCtoTAI(utc_str2, &taitime0);

    cptr3 = basename(argv[1]);
    strcpy((char*) outname, cptr3);
    /*
      cptr3++;
      if (*cptr3 == '1')
        strcpy(outname, "MOD00.P");
      else if (*cptr3 == '0')
        strcpy(outname, "MOD00.A");
      else {
        printf("Satellite Type cannot be determined\n");
        return -1;
      }
     */

    PGS_TD_ASCIItime_AtoB(utc_str2, utc_str);
    utc_str[4] = 0;
    utc_str[8] = 0;
    utc_str[11] = 0;
    utc_str[14] = 0;
    strcat((char*) outname, "_");
    strcat((char*) outname, &utc_str[0]);
    strcat((char*) outname, &utc_str[5]);
    strcat((char*) outname, &utc_str[9]);
    strcat((char*) outname, &utc_str[12]);
    /* strcat(outname, ".pds"); */
    outfp = fopen((char*) outname, "a");
    fseek(outfp, 0, SEEK_SET);

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
                } else {
                    last_pkt_in_buffer = L0_true;
                    fseek(stream, strm_pos, SEEK_SET);
                    continue;
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
                /*
                free( buffer );
                fseek( stream, (long) 0, SEEK_SET);
                return -3;
                 */
            }

            if (packet_length == 642) n_day++;
            if (packet_length == 276) n_night++;

            if (found_start_time == L0_false) {
                memcpy(start_time_tag, &buffer[buf_pos + time_off], time_cnt);
                found_start_time = L0_true;
            } else {
                memcpy(stop_time_tag, &buffer[buf_pos_last + time_off], time_cnt);
                PGS_TD_EOSAMtoTAI(stop_time_tag, &taitime_stop);
                /*	PGS_TD_TAItoUTC( taitime_stop, outname);*/
                /*printf("stoptime =%s\n", outname);*/
            }

            /* write packet if within time period */
            memcpy(time_tag, &buffer[buf_pos + time_off], time_cnt);
            PGS_TD_EOSAMtoTAI(time_tag, &taitime);
            if (taitime >= taitime0 && taitime < taitime0 + gran_length) {
                if (write) {
                    fwrite(&buffer[buf_pos], packet_length, 1, outfp);
                    packet_written = 1;
                }
            }

            if (taitime >= taitime0 + gran_length) {
                if (packet_written == 1) printf("Writing L0 file: %s\n", outname);
                fclose(outfp);

                taitime0 += gran_length;

                cptr3 = basename(argv[1]);
                strcpy((char*) outname, cptr3);
                /*
                        cptr3++;
                        if (*cptr3 == '1')
                          strcpy(outname, "MOD00.P");
                        else if (*cptr3 == '0')
                          strcpy(outname, "MOD00.A");
                        else {
                          printf("Satellite Type cannot be determined\n");
                          return -1;
                        }
                 */

                PGS_TD_TAItoUTC(taitime0 + 1, utc_str2);
                PGS_TD_ASCIItime_AtoB(utc_str2, utc_str);
                utc_str[4] = 0;
                utc_str[8] = 0;
                utc_str[11] = 0;
                utc_str[14] = 0;
                strcat((char*) outname, "_");
                strcat((char*) outname, &utc_str[0]);
                strcat((char*) outname, &utc_str[5]);
                strcat((char*) outname, &utc_str[9]);
                strcat((char*) outname, &utc_str[12]);
                /* strcat(outname, ".pds"); */
                outfp = fopen((char*) outname, "a");
                fseek(outfp, 0, SEEK_SET);

                packet_written = 0;
            }

            strm_pos += (long) packet_length;
            buf_pos_last = buf_pos;
            buf_pos += packet_length;
            n_packets++;
            if ((n_packets % 100000) == 0) {
                printf("%10d packets read %f %f\n", n_packets, taitime, taitime0);
            }
        }

        if (feof(stream) != 0) {
            last_pkt_in_file = L0_true;
        }

        if (ferror(stream) != 0) {
            free(buffer);
            fseek(stream, (long) 0, SEEK_SET);
            return -2;
        }
    }

    free(buffer);

    fclose(stream);
    fclose(outfp);

    return 0;
}
