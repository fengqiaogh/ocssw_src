// procedure to perform quality control on a VIIRS Level-0 packet
// file.  The output is written to STDOUT.

//	Arguments
//     
//     	Name    Type 	I/O 	Description
//     	----	---- 	--- 	-----------
//       vfile   string   I      VIIRS packet file name
//		 rfile	 string	  O		  quality check output file

// Liang Hong, July 22, 2015, Ported from IDL
// V0.1, Aug. 4, 2015
// Liang Hong, Aug. 9, 2016; V0.2; added JPSS-1 VIIRS support

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <timeutils.h>
#include <libgen.h>
#include "l0qc_viirs.h"

#define VERSION "0.2"

int main(int argc, char *argv[]) {
    printf("l0qc_viirs Version: %s (%s %s)\n", VERSION, __DATE__, __TIME__);
    int c, i;
    char *vfile = NULL;
    char *rfile = NULL;
    int bReportFile = 0;

    int imissingpackets = 0;
    int imodetransition = 0;
    int iSameLenCalPacket = 0;
    int itimegaps = 0;
    int ierr = 0;
    int isc = 0;
    int ndy = 0;
    int nnt = 0;
    int ndb = 0;
    int modep = -1;
    int npkts;
    int lpd = -1;
    int tmplp = -1;
    int nSameLenCalPacket = 0;

    while ((c = getopt(argc, argv, "f:r:")) != -1)
        switch (c) {
        case 'f':
            vfile = optarg;
            break;
        case 'r':
            rfile = optarg;
            bReportFile = 1;
            break;
        case '?':
            if (optopt == 'f')
                fprintf(stderr, "Option -f requires an argument.\n");
            else if (optopt == 'r')
                fprintf(stderr, "Option -r requires an argument.\n");
            else
                fprintf(stderr,
                    "Unknown option character `\\x%x'.\n",
                    optopt);
            return 1;
        default:
            abort();
        }


    FILE *infile, *outfile;
    long int endfile = 0;
    infile = fopen(vfile, "r");

    if (infile == NULL) {
        printf("Unable to open input VIIRS file.\n");
        printf("usage: l0qc_viirs -f VIIRS_PDS_input_file -r QC_report_file\n");
        return 1;
    }

    // seek to end of file, find file size and store it to endfile
    fseek(infile, 0, SEEK_END);
    endfile = ftell(infile);
    fseek(infile, 0, SEEK_SET);

    // open report file
    if (bReportFile) {
        outfile = fopen(rfile, "a+");
        if (outfile == NULL) {
            printf("Unable to open output report file.\n");
            return 1;
        }
        fprintf(outfile, "basename=%s\n", basename(vfile));
    }

    printf("Opened VIIRS packet file %s\n", vfile);

    // Initialize reading of VIIRS packet file, 
    //  get first engineering packet
    uint8_t epacket[10000];
    int len_packet;
    init_packet_read(infile, epacket, &len_packet, &endfile);

    // Get packet time
    int32_t pyear, pday, eyear, eday, syear, sday;
    int16_t iyr, ih, mn, mm, dd;
    double sec, stim1, stim2, ptime, etime;
    double usec; // unix seconds
    get_viirs_packet_time(epacket, &pyear, &eyear, &syear, &pday, &eday, &sday, &stim1, &ptime, &etime);
    int32_t jd1 = jday(syear, 1, sday);
    usec = yds2unix(syear, sday, stim1);
    //        iyr = (int16_t) syear;
    //        idy = (int16_t) sday;

    unix2ymdhms(usec, &iyr, &mm, &dd, &ih, &mn, &sec);
    sec = floor(sec * 1000) / 1000; // drop off digits after millisecond
    printf("Start date and time: syear=%d, sday=%d, ih=%d, mn=%d, sec=%06.3f\n", syear, sday, ih, mn, sec);
    if (bReportFile) {
        fprintf(outfile, "start_time=%4d-%02d-%02dT%02d:%02d:%06.3fZ\n", syear, mm, dd, ih, mn, sec);
    }

    // While there are packets in the L0 file
    uint8_t *pbuffer = (uint8_t*) malloc(32000 * 1026);
    char *errstr = (char*) malloc(10000 * 100); // string array to store errors

    while (endfile) {
        if ((isc % 500) == 0) printf("Reading scan %d\n", isc);
        // pbuffer stores current isc scan
        // epacket stores next isc scan
        read_viirs_scan_packets(infile, epacket, &len_packet, pbuffer, &npkts, &endfile);

        // Get sensor mode
        int mode = -1;
        int mpkts = 0;
        int ip;
        for (ip = 0; ip < npkts; ip++) {
            if ((pbuffer[0 + ip * 32000] == 11) && (pbuffer[1 + ip * 32000] != 58)) {
                //for (ip=1;ip<npkts;ip++) {
                //if (pbuffer[0+ip*32000] == 11) {
                if (mode == -1) {
                    mode = pbuffer[46 + ip * 32000];
                    // If day mode
                    if (mode == 4) {
                        mpkts = 479;
                        ndy++;
                    }
                    // If night mode
                    if (mode == 5) {
                        mpkts = 244;
                        nnt++;
                    }
                }

                // Need to check for M11, DNB MGS and LGS packets at night
                if (mode == 5) {
                    if (pbuffer[1 + ip * 32000] == 42) mpkts += 17; //M11
                    if (pbuffer[1 + ip * 32000] == 54) mpkts += 17; //DNB MGS
                    if (pbuffer[1 + ip * 32000] == 55) mpkts += 17; //DNB LGS
                    if (pbuffer[1 + ip * 32000] == 59) mpkts += 17; //DNB HGA
                    if (pbuffer[1 + ip * 32000] == 60) mpkts += 17; //DNB HGB
                }
            }
        }
        if ((mode == 5) && (mpkts >= 278)) ndb++;

        // Check for missing packets
        if (npkts != mpkts) {
            usec = yds2unix(syear, sday, stim1);
            unix2ymdhms(usec, &iyr, &mm, &dd, &ih, &mn, &sec);
            char tmpstr[100];
            sprintf(tmpstr, "%d packets in scan %d mode %d at %02d:%02d:%02d", npkts, isc, mode, ih, mn, (int) sec);
            strcpy(&errstr[ierr * 100], tmpstr);
            printf("%s\n", tmpstr);
            imissingpackets++;
            ierr++;
        }

        if (mode != modep) {
            usec = yds2unix(syear, sday, stim1);
            unix2ymdhms(usec, &iyr, &mm, &dd, &ih, &mn, &sec);
            char tmpstr[100];
            sprintf(tmpstr, "Scan %d transition to mode %d at %02d:%02d:%02d", isc, mode, ih, mn, (int) sec);
            strcpy(&errstr[ierr * 100], tmpstr);
            printf("%s\n", tmpstr);
            imodetransition++;
            //ierr++;  // mode change is informational, not an error
        }

        // Check for same length cal packets
        if (mode == 4) {
            nSameLenCalPacket = 0;
            lpd = -1;
            for (ip = 0; ip < npkts; ip++) {
                if (pbuffer[1 + ip * 32000] == 57) {
                    tmplp = (int) pbuffer[4 + ip * 32000]*256 + (int) pbuffer[5 + ip * 32000] + 7;
                    if (lpd == -1) {
                        lpd = tmplp;
                    } else {
                        if (lpd == tmplp) {
                            nSameLenCalPacket++;
                        }
                    }
                }
                if (nSameLenCalPacket > 1) {
                    usec = yds2unix(syear, sday, stim1);
                    unix2ymdhms(usec, &iyr, &mm, &dd, &ih, &mn, &sec);
                    char tmpstr[100];
                    sprintf(tmpstr, "Same length cal packets in scan %d mode %d at %02d:%02d:%02d", isc, mode, ih, mn, (int) sec);
                    strcpy(&errstr[ierr * 100], tmpstr);
                    printf("%s\n", tmpstr);
                    iSameLenCalPacket++;
                    ierr++;
                    break;
                }
            }
        }

        // If not last scan, check time difference between current scan and next scan
        get_viirs_packet_time(epacket, &pyear, &eyear, &syear, &pday, &eday, &sday, &stim2, &ptime, &etime);
        int32_t jd2 = jday(syear, 1, sday);
        double dtime = stim2 - stim1 + 864 * (jd2 - jd1);
        if (dtime > 1.8) {
            usec = yds2unix(syear, sday, stim2);
            unix2ymdhms(usec, &iyr, &mm, &dd, &ih, &mn, &sec);
            char tmpstr[100];
            sprintf(tmpstr, "Time difference %f at scan %d, at %02d:%02d:%02d", dtime, isc, ih, mn, (int) sec);
            strcpy(&errstr[ierr * 100], tmpstr);
            printf("%s\n", tmpstr);
            itimegaps++;
            ierr++;
        }
        jd1 = jd2;
        stim1 = stim2;
        isc++;
        modep = mode;
    }

    fclose(infile);
    usec = yds2unix(syear, sday, stim2);
    unix2ymdhms(usec, &iyr, &mm, &dd, &ih, &mn, &sec);
    sec = floor(sec * 1000) / 1000; // drop off digits after millisecond
    printf("End date and time year = %d, sday = %d, ih = %d, mn = %d, sec = %06.3f\n", syear, sday, ih, mn, sec);
    printf("%d scans processed\n", isc);
    printf("%d Day-mode scans\n", ndy);
    printf("%d Night-mode scans\n", nnt);
    printf("%d DNB scans\n", ndb);
    printf("QC complete for %s\n", vfile);

    //doy2mmdd(syear, sday, &mm, &dd);
    if (bReportFile) {
        fprintf(outfile, "stop_time=%4d-%02d-%02dT%02d:%02d:%06.3fZ\n", syear, mm, dd, ih, mn, sec);
        fprintf(outfile, "total_scans=%d\n", isc);
        fprintf(outfile, "daymode_scans=%d\n", ndy);
        fprintf(outfile, "nightmode_scans=%d\n", nnt);
        fprintf(outfile, "dnb_scans=%d\n", ndb);
        fprintf(outfile, "mode_transitions=%d\n", imodetransition);
        fprintf(outfile, "errors=%d\n", ierr);
        fprintf(outfile, "incomplete_scans=%d\n", imissingpackets);
        fprintf(outfile, "same_len_cal_packet_scans=%d\n", iSameLenCalPacket);
        fprintf(outfile, "timegap_scans=%d\n", itimegaps);
        for (i = 0; i < ierr; i++) {
            fprintf(outfile, "err_%d=%.*s\n", i + 1, 100, errstr + i * 100);
        }
        fclose(outfile);
    }

    free(pbuffer);
    free(errstr);
    return 0;
}



