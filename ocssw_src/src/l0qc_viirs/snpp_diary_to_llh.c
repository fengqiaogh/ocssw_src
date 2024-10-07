// procedure to read an SNPP S/C diary file and generate an LLH 
// ephemeris file in the Aquarius CODS format.

// 	Arguments
//     
//     	Name    Type 	I/O 	Description
//     	----	---- 	--- 	-----------
//       sfile   string   I      S/C diary packet file name
//		 ofile	 string   O      LLH output file name

// Liang Hong, Jan 14, 2016, Ported from IDL
// V0.1, Feb 1, 2016
// V0.2, Dec 21, 2018, Liang Hong, added error catch when altitude is < 100 km

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <math.h>
#include <libgen.h>
#include <timeutils.h>
#include "snpp_diary_to_llh.h"

#define VERSION "0.2"
#define MIN_ALT 100
#define FILL_VALUE -9999

int main(int argc, char *argv[]) {
    printf("snpp_diary_to_llh Version: %s (%s %s)\n", VERSION, __DATE__, __TIME__);
    int c, i, j;
    int nResampled1mn; // count only the packtes that are closest to 1 minute time tag
    char *sfile = NULL;
    char *ofile = NULL;

    int32_t jd0, iyr, iday;
    int16_t mon, idm, iy, ih, mn, mm, dd;
    double usec, secm;
    double otime[MAX_RECORDS];
    double orb[6 * MAX_RECORDS];
    double atime[MAX_RECORDS];
    double quat[4 * MAX_RECORDS];
    double lon[MAX_RECORDS];
    double lat[MAX_RECORDS];
    double alt[MAX_RECORDS];
    int32_t jd[MAX_RECORDS];
    int result = 1;
    int nErr_val_out_of_range = 0;

    while ((c = getopt(argc, argv, "s:o:")) != -1)
        switch (c) {
        case 's':
            sfile = optarg;
            break;
        case 'o':
            ofile = optarg;
            break;
        case '?':
            if (optopt == 's')
                fprintf(stderr, "Option -s requires an argument.\n");
            else if (optopt == 'o')
                fprintf(stderr, "Option -o requires an argument.\n");
            else
                fprintf(stderr,
                    "Unknown option character `\\x%x'.\n",
                    optopt);
            return 1;
        default:
            abort();
        }

    FILE *infile, *outfile;

    infile = fopen(sfile, "r");
    if (infile == NULL) {
        printf("Unable to open input S/C diary packet file.\n");
        printf("usage: snpp_diary_to_llh -s S/C_diary_packet_file -o llh_output_file\n");
        return 1;
    }
    printf("Opened S/C diary packet file %s\n", sfile);
    struct stat st;
    stat(sfile, &st);
    int nscd = st.st_size / 71;
    uint8_t scdpkts[71 * nscd];
    fread(scdpkts, 1, 71 * nscd, infile);
    fclose(infile);


    // Get start year and day
    //jd0 = swap_endian(fix(dstore(6:7,0),0)) + 2436205 ;Days since 1/1/1958	
    jd0 = (scdpkts[0 * 71 + 6] << 8) + (scdpkts[0 * 71 + 7]) + 2436205;
    jdate(jd0, &iyr, &iday);

    // Unpack and convert diary data
    result = convert_diary(nscd, scdpkts, otime, orb, atime, quat);
    if (result != 0) {
        printf("diary data conversion problem!\n");
        return 1;
    }

    // Find records closest to the minute from the right. e.g. use 01:02:02 in stead of 01:01:59
    nResampled1mn = 0;
    double new_otime[MAX_RECORDS], new_orb[6 * MAX_RECORDS], min_dotime, tmp_dotime, curr1m, omm;
    min_dotime = 6000; // otime - round_otime_to_closest_1_minute, in seconds
    curr1m = round(otime[0] / 60.0)*60.0;
    new_otime[0] = otime[0];
    for (j = 0; j < 6; j++) new_orb[j + 0 * 6] = orb[j + 0 * 6];

    for (i = 1; i < nscd; i++) {
        // omm = round(otime[i]/60.0)*60.0;
        // tmp_dotime = fabs(otime[i] - omm);
        omm = floor(otime[i] / 60.0)*60.0;
        tmp_dotime = otime[i] - omm;
        if (omm > curr1m) {
            curr1m = omm;
            min_dotime = tmp_dotime;
            nResampled1mn++;
            new_otime[nResampled1mn] = otime[i];
            for (j = 0; j < 6; j++) new_orb[j + nResampled1mn * 6] = orb[j + i * 6] / 1000.0;
        } else {
            if (tmp_dotime < min_dotime) {
                new_otime[nResampled1mn] = otime[i];
                for (j = 0; j < 6; j++) new_orb[j + nResampled1mn * 6] = orb[j + i * 6] / 1000.0;
                min_dotime = tmp_dotime;
            }
        }
    }

    result = orb2lla(nResampled1mn, new_orb, lon, lat, alt);
    if (result != 0) {
        printf("orbit to lon, lat and alt conversion problem!\n");
        return 1;
    }

    // Check for day rollover
    for (i = 0; i < nResampled1mn; i++) {
        jd[i] = jday(iyr, 1, iday);
        jd[i] = jd[i] + (int32_t) (new_otime[i] / 86400.0);
        new_otime[i] = new_otime[i] - floor(new_otime[i] / 86400.0)*86400.0;
    }

    // Generate output file name
    jdate(jd[nResampled1mn - 1], &iyr, &iday);
    yd2md((int16_t) iyr, (int16_t) iday, &mon, &idm);
    //sod2hms,sec,ih,mn,secm
    usec = yds2unix((int16_t) iyr, (int16_t) iday, new_otime[nResampled1mn - 1]);
    unix2ymdhms(usec, &iy, &mm, &dd, &ih, &mn, &secm);
    //sprintf(cfile,"GSFC_%4d%02d%02d_%02d%02d%02d_SNPP_ORBEPHEM_GEOD_LLH_O.TXT", iyr, mon, idm, ih, mn, (int)secm);
    //printf("cfile=%s\n",cfile);

	// no valid packets in input file
	if (nResampled1mn<1) {
		printf("Corrupted input file, no valid packets!\n");
        return 1;
	}
	
	// open llh output file	
    outfile = fopen(ofile, "w");
    if (outfile == NULL) {
        printf("Unable to open output file.\n");
        return 1;
    }
    fprintf(outfile, "basename=%s\n", basename(sfile));
    //printf("basename = %s\n",basename(sfile));
    
    // Output header record
    fprintf(outfile, "*HEADER SNPP    %d LLH    UTC    WGS84    PRO    37849 11061A    GSFC    %d/%02d/%02d %02d:%02d:%02d\n", nscd, iyr, mon, idm, ih, mn, (int) secm);

    // Output data records
    for (i = 0; i < nResampled1mn; i++) {
        jdate(jd[i], &iyr, &iday);
        yd2md((int16_t) iyr, (int16_t) iday, &mon, &idm);
        usec = yds2unix(iyr, iday, new_otime[i]);
        unix2ymdhms(usec, &iy, &mm, &dd, &ih, &mn, &secm);
        if (alt[i]<MIN_ALT) {
        	lat[i] = FILL_VALUE;
        	lon[i] = FILL_VALUE;
        	alt[i] = FILL_VALUE;			// fill value for bogus altitude 
        	nErr_val_out_of_range++;
        }
        fprintf(outfile, "%d/%02d/%02d %02d:%02d:%011.8f	%15.5f	%15.5f	%15.3f\n", iyr, mon, idm, ih, mn, secm, lat[i], lon[i], alt[i]);
    }
    fclose(outfile);

    printf("Total LLH output: %d\n", nResampled1mn);
    if (nErr_val_out_of_range>0) printf("Records with value out of range: %d\n", nErr_val_out_of_range);

    return 0;
}
