/*
This program reads a SeaWiFS level-1A HDF file and writes out the
data in SeaWiFS level-0 format.
The program also writes out a text file containing the station
identification information that is needed by SWl01.  The name
of this text file is "<output levl-0 name>.station_info" .

Norman Kuring		11-Dec-1997
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mfhdf.h"
#include "genutils.h"

#define VERSION "1.0"
#define USAGE "Usage: %s SeaWiFS_level_1a_filename output_filename\n"

#define READ_GLBL_ATTR(nam,ptr) {                                           \
  if(SDreadattr(sd_id,SDfindattr(sd_id,(nam)),(VOIDP)(ptr))){               \
    fprintf(stderr,                                                         \
    "-E- %s line %d: Could not get global attribute, %s, from file, %s.\n", \
    __FILE__,__LINE__,(nam),argv[1]);                                       \
    exit(EXIT_FAILURE);                                                     \
  }                                                                         \
}

#define READ_SDS(nam,ptr,s0,s1,s2,e0,e1,e2) {                           \
  int32 start[3];                                                       \
  int32 edge[3];                                                        \
  edge[0]=(e0); edge[1]=(e1); edge[2]=(e2);                             \
  start[0]=(s0); start[1]=(s1); start[2]=(s2);                          \
  if(SDreaddata(SDselect(sd_id, SDnametoindex(sd_id, (nam))),           \
  start, NULL, edge, (VOIDP)(ptr)) == FAIL){                            \
    fprintf(stderr,"-E- %s line %d: Could not read SDS, %s.\n",         \
    __FILE__,__LINE__,(nam));                                           \
    exit(EXIT_FAILURE);                                                 \
  }                                                                     \
}

#define MALLOC(ptr,typ,num) {                                           \
  (ptr) = (typ *)malloc((num) * sizeof(typ));                           \
  if((ptr) == NULL){                                                    \
    fprintf(stderr,"-E- %s line %d: Memory allocation failure.\n",      \
    __FILE__,__LINE__);                                                 \
    exit(EXIT_FAILURE);                                                 \
  }                                                                     \
}

#define MINOR_FRAME_SIZE 21504

double doubleutime(int16, int16, int32);
void get_frame_data(int32, int32, int32, int16 *, int16 *, uint8 *,
        int16 *, int16 *, int16 *, int16 *, int16 *, int16 *);

int main(int argc, char *argv[]) {

    FILE *fh;
    int32 sd_id;
    int16 syear, sday, eyear, eday;
    int32 smsec, emsec, npix, nscan, s;
    uint8 *s_flags;
    int is_gac;
    time_t aos, los;
    uint16 bit_error_rate_denominator;
    uint8 bit_error_rate_numerator;
    uint32 ber_numerator = 0, ber_denominator = 0;
    static uint8 frame[MINOR_FRAME_SIZE];
    uint8 *fp;
    int16 sc_id[2], sc_ttag[4], inst_tlm[44];
    int16 start_syn[8], stop_syn[8], dark_rest[8], gain_tdi[8];
    int16 *l1a_data;
    uint8 sc_soh[775];
    static uint8 header[512];
    uint8 dtype[32], center[1024], station[1024];
    uint8 is_hrpt;
    int base;
    char *stationfile;
    static char code[4];
    float32 lat, lon;

    if (argc != 3) {
        fprintf(stderr, "%s %s (%s %s)\n", argv[0], VERSION, __DATE__, __TIME__);
        fprintf(stderr, USAGE, argv[0]);
        exit(EXIT_FAILURE);
    }

    /* Make sure the input file exists. */
    if ((fh = fopen(argv[1], "rb")) == NULL) {
        fprintf(stderr, "-E- %s line %d: Could not open file, %s .  ",
                __FILE__, __LINE__, argv[1]);
        perror(NULL);
        exit(EXIT_FAILURE);
    } else {
        fclose(fh);
    }

    /* Make sure the output file does not already exist. */
    if ((fh = fopen(argv[2], "rb")) != NULL) {
        fprintf(stderr, "-E- %s line %d: File, %s , already exists.\n",
                __FILE__, __LINE__, argv[2]);
        exit(EXIT_FAILURE);
    }

    /* Compose a station info filename and make sure it doesn't exist. */
    MALLOC(stationfile, char, strlen(argv[2]) + 14);
    strcpy(stationfile, argv[2]);
    strcat(stationfile, ".station_info");
    if ((fh = fopen(stationfile, "rb")) != NULL) {
        fprintf(stderr, "-E- %s line %d: File, %s , already exists.\n",
                __FILE__, __LINE__, stationfile);
        exit(EXIT_FAILURE);
    }

    /* Open the HDF input file */
    sd_id = SDstart(argv[1], DFACC_RDONLY);
    if (sd_id == FAIL) {
        fprintf(stderr, "-E- %s line %d: SDstart(%s, %d) failed.\n",
                __FILE__, __LINE__, argv[1], DFACC_RDONLY);
        exit(EXIT_FAILURE);
    }

    /* Read some of the level-1A global attributes. */
    READ_GLBL_ATTR("Data Type", dtype);
    READ_GLBL_ATTR("Start Year", &syear);
    READ_GLBL_ATTR("Start Day", &sday);
    READ_GLBL_ATTR("Start Millisec", &smsec);
    READ_GLBL_ATTR("End Year", &eyear);
    READ_GLBL_ATTR("End Day", &eday);
    READ_GLBL_ATTR("End Millisec", &emsec);
    READ_GLBL_ATTR("Pixels per Scan Line", &npix);
    READ_GLBL_ATTR("Number of Scan Lines", &nscan);
    READ_GLBL_ATTR("Data Center", center);
    READ_GLBL_ATTR("Station Name", station);
    READ_GLBL_ATTR("Station Latitude", &lat);
    READ_GLBL_ATTR("Station Longitude", &lon);

    /* Reformat the start and end times. */
    aos = (time_t) floor(doubleutime(syear, sday, smsec));
    los = (time_t) ceil(doubleutime(eyear, eday, emsec));

    /* Is it HRPT data? */
    is_hrpt = strcmp((char*) dtype, "HRPT") ? 0 : 1;

    /* Is it GAC data? */
    switch (npix) {
    case 248: is_gac = 1;
        break;
    case 1285: is_gac = 0;
        break;
    default:
        fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "\"Pixels per Scan Line\" global attribute of file, ");
        fprintf(stderr, "%s , has the unexpected value, %d.\n", argv[1], npix);
        exit(EXIT_FAILURE);
    }

    /* Consistency check. */
    if (is_hrpt && is_gac) {
        fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "\"Pixels per Scan Line\" (= %d) ", npix);
        fprintf(stderr, "is inconsistent with \"Data Type\" (= \"%s\") ", dtype);
        fprintf(stderr, "in file, %s .\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    /* Get a 3-letter code from the input filename if this is HRPT data. */
    if (is_hrpt) {
        char *p;
        p = argv[1] + strlen(argv[1]) - 3;
        strcpy(code, p);
    }

    /* Write out the station information. */
    if ((fh = fopen(stationfile, "wb")) == NULL) {
        fprintf(stderr, "-E- %s line %d: Could not create file, %s . ",
                __FILE__, __LINE__, stationfile);
        perror(NULL);
        exit(EXIT_FAILURE);
    }
    if (is_hrpt)
        fprintf(fh, "Code: %s\n", code);
    fprintf(fh, "Data Center: %s\n", center);
    fprintf(fh, "Station Name: %s\n", station);
    fprintf(fh, "Station Latitude: %f\n", lat);
    fprintf(fh, "Station Longitude: %f\n", lon);
    fclose(fh); /* Close the station info file. */

    /*
    Read the s_flags SDS.  I need this for bit error rates and
    GAC scan line positions within a minor frame.
     */
    MALLOC(s_flags, uint8, nscan * 4);
    READ_SDS("s_flags", s_flags, 0, 0, 0, nscan, 4, 1);

    /*
    Due to some confusion about whether the GAC-line-number byte is
    zero-based or one-based, different level-1A files may store that
    byte in different ways.  So, I make a guess at the scheme used
    to write the current input file by checking the very first GAC
    line number stored.
     */
    if (s_flags[2] == 0 || s_flags[2] == 5 || s_flags[2] == 10)
        base = 0;
    else
        base = 1;

    /* Sum up the bit error rate values. */
    if (is_gac) {
        int mfnumber;
        int prevmfnumber;
        int gacposition;
        int prevposition;

        for (s = 0; s < nscan; s++) {
            mfnumber = (s_flags[4 * s + 2] - base) / 5;
            gacposition = (s_flags[4 * s + 2] - base) % 5;
            if (s == 0 || gacposition <= prevposition || mfnumber != prevmfnumber) {
                ber_numerator += s_flags[4 * s];
                ber_denominator += 5 * s_flags[4 * s + 3];
            }
            prevmfnumber = mfnumber;
            prevposition = gacposition;
        }
    } else {
        for (s = 0; s < nscan; s++) {
            ber_numerator += s_flags[4 * s];
            ber_denominator += 5 * s_flags[4 * s + 3];
        }
    }

    /* Fill in the header. */
    strcpy((char*) header, "CWIF");

    if (endianess() == 1) {
        swapc_bytes((char *) &ber_denominator, 4, 1);
        swapc_bytes((char *) &ber_numerator, 4, 1);
        swapc_bytes((char *) &aos, 4, 1);
        swapc_bytes((char *) &los, 4, 1);
    }

    header[5] = is_hrpt;
    memcpy(&header[316], &ber_denominator, 4);
    memcpy(&header[320], &ber_numerator, 4);
    memcpy(&header[332], &aos, 4);
    memcpy(&header[336], &los, 4);
    header[511] = 1; /* Signal that this L0 file was generated from L1A. */

    /* Allocate some memory for the level-1a image data. */
    MALLOC(l1a_data, int16, npix * 8);

    /* Create the output file and write out the header. */
    if ((fh = fopen(argv[2], "wb")) == NULL) {
        fprintf(stderr, "-E- %s line %d: Could not create file, %s . ",
                __FILE__, __LINE__, argv[2]);
        perror(NULL);
        exit(EXIT_FAILURE);
    }
    fwrite(header, 1, 512, fh);

    if (is_gac) {
        int mfnumber;
        int prevmfnumber;
        int gacposition;
        int prevposition = 5;
        int mf_is_empty = 1;

        for (s = 0; s < nscan; s++) {
            mfnumber = (s_flags[4 * s + 2] - base) / 5;
            gacposition = (s_flags[4 * s + 2] - base) % 5;
            if ((gacposition <= prevposition || mfnumber != prevmfnumber) && s != 0) {
                /* Write out the previous minor frame. */
                fwrite(frame, 1, MINOR_FRAME_SIZE, fh);
                /* Zero out minor frame buffer. */
                memset(frame, 0, MINOR_FRAME_SIZE);
                /* Signal that the minor fame is not yet initialized. */
                mf_is_empty = 1;
            }
            prevmfnumber = mfnumber;
            prevposition = gacposition;

            get_frame_data(sd_id, npix, s, sc_id, sc_ttag, sc_soh, inst_tlm,
                    l1a_data, start_syn, stop_syn, dark_rest, gain_tdi);

            /*
            Copy the various pieces to the appropriate places in the minor frame.
            Only one of the five scans in a GAC minor frame is used to fill in
            the bit error rate, the spacecraft ID, the spacecraft time tag, and
            the state-of-health telemetry.  The rest of the pieces are copied
            from each scan line.
             */
            if (mf_is_empty) {
                bit_error_rate_denominator = 5 * s_flags[4 * s + 3];
                bit_error_rate_numerator = s_flags[4 * s ];
                fp = frame;
                if (endianess() == 1) {
                    swapc_bytes((char *) &bit_error_rate_denominator, 2, 1);
                    swapc_bytes((char *) sc_id, 2, 2);
                    swapc_bytes((char *) sc_ttag, 2, 4);
                }
                memcpy(fp, &bit_error_rate_denominator, 2);
                fp += 2;
                memcpy(fp, &bit_error_rate_numerator, 1);
                fp++;
                memcpy(fp, sc_id, 4);
                fp += 4;
                memcpy(fp, sc_ttag, 8);
                fp += 8;
                memcpy(fp, sc_soh, 775);
                mf_is_empty = 0;
            }
            fp = frame + 790 + gacposition * ((8 + 8 + 8 + 248 * 8 + 8)*2);
            if (endianess() == 1) {
                swapc_bytes((char *) gain_tdi, 2, 8);
                swapc_bytes((char *) start_syn, 2, 8);
                swapc_bytes((char *) dark_rest, 2, 8);
                swapc_bytes((char *) l1a_data, 2, 1984);
                swapc_bytes((char *) stop_syn, 2, 8);
                swapc_bytes((char *) inst_tlm, 2, 44);
            }
            memcpy(fp, gain_tdi, 16);
            fp += 16;
            memcpy(fp, start_syn, 16);
            fp += 16;
            memcpy(fp, dark_rest, 16);
            fp += 16;
            memcpy(fp, l1a_data, 3968);
            fp += 3968;
            memcpy(fp, stop_syn, 16);
            fp = frame + 790 + 5 * (8 + 8 + 8 + 248 * 8 + 8)*2 + gacposition * (44 * 2);
            memcpy(fp, inst_tlm, 88);
        }

        /* Write out the final minor frame of GAC data. */
        fwrite(frame, 1, MINOR_FRAME_SIZE, fh);

    } else { /* non-GAC data */
        for (s = 0; s < nscan; s++) {
            get_frame_data(sd_id, npix, s, sc_id, sc_ttag, sc_soh, inst_tlm,
                    l1a_data, start_syn, stop_syn, dark_rest, gain_tdi);
            bit_error_rate_denominator = 5 * s_flags[4 * s + 3];
            bit_error_rate_numerator = s_flags[4 * s ];
            fp = frame;
            if (endianess() == 1) {
                swapc_bytes((char *) &bit_error_rate_denominator, 2, 1);
                swapc_bytes((char *) sc_id, 2, 2);
                swapc_bytes((char *) sc_ttag, 2, 4);
                swapc_bytes((char *) inst_tlm, 2, 44);
                swapc_bytes((char *) gain_tdi, 2, 8);
                swapc_bytes((char *) start_syn, 2, 8);
                swapc_bytes((char *) dark_rest, 2, 8);
                swapc_bytes((char *) l1a_data, 2, 10280);
                swapc_bytes((char *) stop_syn, 2, 8);
            }
            memcpy(fp, &bit_error_rate_denominator, 2);
            fp += 2;
            memcpy(fp, &bit_error_rate_numerator, 1);
            fp++;
            memcpy(fp, sc_id, 4);
            fp += 4;
            memcpy(fp, sc_ttag, 8);
            fp += 8;
            memcpy(fp, sc_soh, 775);
            fp += 775;
            memcpy(fp, inst_tlm, 88);
            fp += 88;
            memcpy(fp, gain_tdi, 16);
            fp += 16;
            memcpy(fp, start_syn, 16);
            fp += 16;
            memcpy(fp, dark_rest, 16);
            fp += 16;
            memcpy(fp, l1a_data, 20560);
            fp += 20560;
            memcpy(fp, stop_syn, 16);
            /* Write out the minor frame. */
            fwrite(frame, 1, MINOR_FRAME_SIZE, fh);
            /* Zero out minor frame buffer. */
            memset(frame, 0, MINOR_FRAME_SIZE);
        }
    }
    fclose(fh); /* Close the level-0 file. */
    exit(0);
}





/* Leap year definition from Kernighan and Ritchie page 37 */
#define         IS_LEAP_YEAR(y)         ( (!((y)%4) && (y)%100) || !((y)%400) )
#define         EPOCH_YEAR              1970
#define         SECONDS_PER_YEAR        31536000
#define         SECONDS_PER_LEAP_YEAR   31622400

double doubleutime(int16 year, int16 day, int32 msec) {

    double utime = 0;
    int y;

    for (y = EPOCH_YEAR; y < year; y++)
        utime += IS_LEAP_YEAR(y) ? SECONDS_PER_LEAP_YEAR : SECONDS_PER_YEAR;
    utime += (day - 1) * 86400;
    utime += msec / 1000.0;
    return (utime);
}

void
get_frame_data(
        int32 sd_id, /* This variable is part of the READ_SDS macro. */
        int32 npix,
        int32 scan,
        int16 *sc_id,
        int16 *sc_ttag,
        uint8 *sc_soh,
        int16 *inst_tlm,
        int16 *l1a_data,
        int16 *start_syn,
        int16 *stop_syn,
        int16 *dark_rest,
        int16 *gain_tdi
        ) {
    int16 gain[8], tdi[8];
    int b;

    READ_SDS("sc_id", sc_id, scan, 0, 0, 1, 2, 1);
    READ_SDS("sc_ttag", sc_ttag, scan, 0, 0, 1, 4, 1);
    READ_SDS("sc_soh", sc_soh, scan, 0, 0, 1, 775, 1);
    READ_SDS("inst_tlm", inst_tlm, scan, 0, 0, 1, 44, 1);
    READ_SDS("l1a_data", l1a_data, scan, 0, 0, 1, npix, 8);
    READ_SDS("start_syn", start_syn, scan, 0, 0, 1, 8, 1);
    READ_SDS("stop_syn", stop_syn, scan, 0, 0, 1, 8, 1);
    READ_SDS("dark_rest", dark_rest, scan, 0, 0, 1, 8, 1);
    READ_SDS("gain", gain, scan, 0, 0, 1, 8, 1);
    READ_SDS("tdi", tdi, scan, 0, 0, 1, 8, 1);

    for (b = 0; b < 8; b++) {
        gain_tdi[b] = ((gain[b] << 8) & 0x0300) | (tdi[b] & 0x00FF);
    }
}
