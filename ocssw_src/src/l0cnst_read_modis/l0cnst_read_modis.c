#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <PGS_TD.h>

#define buffsize_mb 0.5
#define primary_hdr_size 6

int main(int argc, char *argv[]) {
    int tptr;

    //static int pkt_len_off = 4;
    //static int time_cnt = 8;
    //static int day_pkt_size = 642;
    //static int night_pkt_size = 276;
    static unsigned char start_time_tag[8];
    static unsigned char stop_time_tag[8];

    double taitime_start;
    double taitime_stop;

    unsigned char inbuf[384];
    FILE *stream;

    printf("read_constructor_file: Version as of 06/18/04\n\n");

    stream = fopen(argv[1], "r");
    if (stream == NULL) {
        printf("%s not found.\n", argv[1]);
        exit(-1);
    }
    fseek(stream, (long) 0, SEEK_SET);


    fread(inbuf, sizeof (char), 384, stream);

    tptr = 64 + 16 * inbuf[51];
    memcpy(start_time_tag, &inbuf[tptr], 8);
    memcpy(stop_time_tag, &inbuf[tptr + 8], 8);

    fclose(stream);

    PGS_TD_EOSAMtoTAI(start_time_tag, &taitime_start);
    PGS_TD_EOSAMtoTAI(stop_time_tag, &taitime_stop);
    printf("starttimeTAI=%f\n", taitime_start);
    printf("stoptimeTAI =%f\n", taitime_stop);

    PGS_TD_TAItoUTC(taitime_start, (char*) inbuf);
    printf("starttime=%s\n", inbuf);
    PGS_TD_TAItoUTC(taitime_stop, (char*) inbuf);
    printf("stoptime =%s\n", inbuf);

    printf("granule length =%f\n", (taitime_stop - taitime_start) / 60.0);


    return 0;
}
