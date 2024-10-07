// Ported from IDL procedure to perform quality control on SNPP spacecraft packet data
// files.  The output is written to STDOUT.

//	Arguments
//     
//     	Name    Type 	I/O 	Description
//     	----	---- 	--- 	-----------
//       sfile   string   I      S/C diary packet file name
//       afile   string   I      ADCS telemetry packet file name
//       bfile   string   I      Bus-critical telemetry packet file name
//       gfile   string   I      GPS telemetry packet file name (J2 only)
//		 rfile	 string	  O		 quality check output file

// Liang Hong, July 22, 2015
// V0.1, Aug. 4, 2015
// Liang Hong, Aug. 9, 2016; V 0.2; added support for JPSS-1 VIIRS
// Liang Hong, Dec. 30, 2020; V1.0.0;
//   -- err # from index 1; 
//   -- GPS telemetry input and QC added (J2)
//   -- check for the J2 S/C ID and set packet sizes.
//   -- sample rates for the diary and ADCS data are 10Hz for J2 vs. 1Hz for SNPP and J1, so the delta time limits are modified
// Liang Hong, Aug. 23, 2021; V1.1.0;
//   -- read s/c diary data by each packet instead of whole file


#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <libgen.h>
#include <timeutils.h>

#define VERSION "1.1.0_2021-08-23"

int main(int argc, char *argv[]) {
    printf("scqc_viirs Version: %s (%s %s)\n", VERSION, __DATE__, __TIME__);
    char *sfile, *afile, *bfile, *gfile, *rfile;
    int c;
    int bReportFile = 0;
    int plat = -1;  // platform identifier. 
    int scdpktsize = 71; 
    int apktsize;
    int bpktsize;
    int gpktsize;

    while ((c = getopt(argc, argv, "s:a:b:g:r:")) != -1)
        switch (c) {
        case 's':
            sfile = optarg;
            break;
        case 'a':
            afile = optarg;
            break;
        case 'b':
            bfile = optarg;
            break;
        case 'g':
            gfile = optarg;
            break;
        case 'r':
            rfile = optarg;
            bReportFile = 1;
            break;
        case '?':
            if (optopt == 's')
                fprintf(stderr, "Option -s requires an argument.\n");
            else if (optopt == 'a')
                fprintf(stderr, "Option -a requires an argument.\n");
            else if (optopt == 'b')
                fprintf(stderr, "Option -b requires an argument.\n");
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
    if (sfile == NULL || afile == NULL || bfile == NULL || rfile == NULL) {
        printf("unable to open input files.\n");
        printf("usage: scqc_viirs -s S/C_diary_file -a ADCS_file -b S/C_bus_file [-g GPS_file] -r QC_report_file\n");
        return 1;
    }

    FILE *infile, *infile_a, *infile_b, *infile_g, *outfile;
    infile = fopen(sfile, "r");
    printf("Opened S/C diary packet file %s\n", sfile);
    struct stat st;
    stat(sfile, &st);
    int nscd = st.st_size / scdpktsize;
    uint8_t scdpkts[scdpktsize];
    fread(scdpkts, 1, scdpktsize, infile);
    //fclose(infile);
    // convert_diary,nscd,scdpkts,iyrsc,idysc,otime,orb,atime,quat

    // Check spacecraft ID
    if (scdpkts[14] == 157) {
        plat = 0;
        apktsize = 355;
        bpktsize = 207;
		// scanp = 1.7793
    }
    if (scdpkts[14] == 159) {
        plat = 1;
        apktsize = 393;
        bpktsize = 183;
		//   scanp 1.7864
    }
    if (scdpkts[14] == 177) {
		plat = 2;
		apktsize = 493;
		bpktsize = 212;
		gpktsize = 374;
		//   scanp 1.7864

		// quit with error message if JPSS-2 is being processed without gps file as input
		if (gfile == NULL) {
			printf("missing gps file for JPSS-2 (P177).\n");
			printf("usage: scqc_viirs -s S/C_diary_file -a ADCS_file -b S/C_bus_file -g GPS_file -r QC_report_file\n");
			return 1;
		}
    }

    infile_a = fopen(afile, "r");
    printf("Opened ADCS telemetry packet file %s\n", afile);
    stat(afile, &st);
    int nadc = st.st_size / apktsize;
    uint8_t adcpkts[apktsize];

    infile_b = fopen(bfile, "r");
    printf("Opened bus-critical telemetry packet file %s\n", bfile);
    stat(bfile, &st);
    int nbus = st.st_size / bpktsize;
    uint8_t buspkts[bpktsize];
    
    int ngps = 0;
    if (plat >= 2) {
		infile_g = fopen(gfile, "r");
		printf("Opened GPS telemetry packet file %s\n\n", gfile);
		stat(gfile, &st);
		ngps = st.st_size / gpktsize;
    }
    
    // open report file
    if (bReportFile) {
        outfile = fopen(rfile, "a+");
        if (outfile == NULL) {
            printf("Unable to open output report file.\n");
            return 1;
        }
        fprintf(outfile, "scdiary_basename=%s\n", basename(sfile));
        fprintf(outfile, "adcs_basename=%s\n", basename(afile));
        fprintf(outfile, "btlm_basename=%s\n", basename(bfile));
        if (plat >= 2)
            fprintf(outfile, "gps_basename=%s\n", basename(gfile));
    }

    double scdtime[nscd], adctime[nadc], bustime[nbus], gpstime[ngps], delt[2], dtime;

    // Get start year and day

    int32_t jd0, iyr, iday;
    jd0 = (scdpkts[6] << 8) + (scdpkts[7]) + 2436205;
    jdate(jd0, &iyr, &iday);

    // Loop through packets and extract times
    // S/C diary
    int16_t iy, ih, mn, mm, dd;
    int32_t ccsds_iy, idy, jd;
    int i, m;
    int ierr = 0;
    int itimegaps = 0;
    double sec, secd, usec;
    uint8_t cctime[8];
    char *errstr = (char*) malloc(10000 * 100); // string array to store errors

	if (plat <= 1) {
		delt[0] = 0.9;
		delt[1] = 1.1; 
	} else {
		delt[0] = 0.09;
		delt[1] = 0.11;
	}
	
	for (m = 0; m < 8; m++) cctime[m] = scdpkts[m + 6];
	ccsds_to_yds(cctime, &ccsds_iy, &idy, &secd);
	iy = (int16_t) ccsds_iy;
	jd = jday((int32_t) iy, 1, idy);
	scdtime[0] = secd + (jd - jd0)*86400;
    //sod2hms((int)secd,&ih,&mn,&sec);
	usec = yds2unix(iy, idy, secd);
	unix2ymdhms(usec, &iy, &mm, &dd, &ih, &mn, &sec);
	printf("S/C diary start date and time, year = %d, day = %d, hour = %d, min = %d, sec = %d\n", iy, idy, ih, mn, (int) sec);
	if (bReportFile) {
		sec = floor(sec * 1000) / 1000;
		fprintf(outfile, "scdiary_start_time=%4d-%02d-%02dT%02d:%02d:%06.3fZ\n", iy, mm, dd, ih, mn, sec);
	}
	
    for (i = 1; i < nscd; i++) {
        fread(scdpkts, 1, scdpktsize, infile);
        for (m = 0; m < 8; m++) cctime[m] = scdpkts[m + 6];
        ccsds_to_yds(cctime, &ccsds_iy, &idy, &secd);
        iy = (int16_t) ccsds_iy;
        jd = jday((int32_t) iy, 1, idy);
        scdtime[i] = secd + (jd - jd0)*86400;
		dtime = scdtime[i] - scdtime[i - 1];
		if ((dtime < delt[0]) || (dtime > delt[1])) {
			//sod2hms((int)secd,&ih,&mn,&sec);
			usec = yds2unix(iy, idy, secd);
			unix2ymdhms(usec, &iy, &mm, &dd, &ih, &mn, &sec);
			char tmpstr[100];
			sprintf(tmpstr, "Time difference %f at record %d, at%02d:%02d:%02d", dtime, i, ih, mn, (int) sec);
			strcpy(&errstr[ierr * 100], tmpstr);
			printf("%s\n", tmpstr);
			itimegaps++;
			ierr++;
		}
    }
    fclose(infile);
    
    //sod2hms((int)secd,&ih,&mn,&sec);
    usec = yds2unix(iy, idy, secd);
    unix2ymdhms(usec, &iy, &mm, &dd, &ih, &mn, &sec);
    printf("S/C diary end date and time, year = %d, day = %d, hour = %d, min = %d, sec = %d\n", iy, idy, ih, mn, (int) sec);
    printf("%d records processed\n\n", nscd);
    if (bReportFile) {
        sec = floor(sec * 1000) / 1000;
        fprintf(outfile, "scdiary_stop_time=%4d-%02d-%02dT%02d:%02d:%06.3fZ\n", iy, mm, dd, ih, mn, sec);
        fprintf(outfile, "scdiary_records=%d\n", nscd);
        fprintf(outfile, "scdiary_errors=%d\n", ierr);
        fprintf(outfile, "scdiary_timegap_records=%d\n", itimegaps);
        for (i = 0; i < ierr; i++) {
            fprintf(outfile, "scdiary_err_%d=%.*s\n", i+1, 100, errstr + i * 100);
        }
    }

    // ADCS
    ierr = 0;
    itimegaps = 0;
    free(errstr);
    errstr = (char*) malloc(10000 * 100);
    for (i = 0; i < nadc; i++) {
    	fread(adcpkts, 1, apktsize, infile_a);
        for (m = 0; m < 8; m++) cctime[m] = adcpkts[m + 6];
        ccsds_to_yds(cctime, &ccsds_iy, &idy, &secd);
        iy = (int16_t) ccsds_iy;
        jd = jday((int32_t) iy, 1, idy);
        adctime[i] = secd + (jd - jd0)*86400;
        if (i == 0) {
            //sod2hms((int)secd,&ih,&mn,&sec);
            usec = yds2unix(iy, idy, secd);
            unix2ymdhms(usec, &iy, &mm, &dd, &ih, &mn, &sec);
            printf("ADCS start date and time, year = %d, day = %d, hour = %d, min = %d, sec = %d\n", iy, idy, ih, mn, (int) sec);
            if (bReportFile) {
                sec = floor(sec * 1000) / 1000;
                fprintf(outfile, "adcs_start_time=%4d-%02d-%02dT%02d:%02d:%06.3fZ\n", iy, mm, dd, ih, mn, sec);
            }
        } else {
            dtime = adctime[i] - adctime[i - 1];
            if ((dtime < delt[0]) || (dtime > delt[1])) {
                //sod2hms((int)secd,&ih,&mn,&sec);
                usec = yds2unix(iy, idy, secd);
                unix2ymdhms(usec, &iy, &mm, &dd, &ih, &mn, &sec);
                char tmpstr[100];
                sprintf(tmpstr, "Time difference %f at record %d, at%02d:%02d:%02d", dtime, i, ih, mn, (int) sec);
                strcpy(&errstr[ierr * 100], tmpstr);
                printf("%s\n", tmpstr);
                itimegaps++;
                ierr++;
            }
        }
    }
    fclose(infile_a);
    
    //sod2hms((int)secd,&ih,&mn,&sec);
    usec = yds2unix(iy, idy, secd);
    unix2ymdhms(usec, &iy, &mm, &dd, &ih, &mn, &sec);
    printf("ADCS end date and time, year = %d, day = %d, hour = %d, min = %d, sec = %d\n", iy, idy, ih, mn, (int) sec);
    printf("%d records processed\n\n", nadc);
    if (bReportFile) {
        sec = floor(sec * 1000) / 1000;
        fprintf(outfile, "adcs_stop_time=%4d-%02d-%02dT%02d:%02d:%06.3fZ\n", iy, mm, dd, ih, mn, sec);
        fprintf(outfile, "adcs_records=%d\n", nadc);
        fprintf(outfile, "adcs_errors=%d\n", ierr);
        fprintf(outfile, "adcs_timegap_records=%d\n", itimegaps);
        for (i = 0; i < ierr; i++) {
            fprintf(outfile, "adcs_err_%d=%.*s\n", i+1, 100, errstr + i * 100);
        }
    }

    // Bus-critical
    ierr = 0;
    itimegaps = 0;
    free(errstr);
    errstr = (char*) malloc(10000 * 100);
    for (i = 0; i < nbus; i++) {
        fread(buspkts, 1, bpktsize, infile_b);
        for (m = 0; m < 8; m++) cctime[m] = buspkts[m + 6];
        ccsds_to_yds(cctime, &ccsds_iy, &idy, &secd);
        iy = (int16_t) ccsds_iy;
        jd = jday((int32_t) iy, 1, idy);
        bustime[i] = secd + (jd - jd0)*86400;
        if (i == 0) {
            //sod2hms((int)secd,&ih,&mn,&sec);
            usec = yds2unix(iy, idy, secd);
            unix2ymdhms(usec, &iy, &mm, &dd, &ih, &mn, &sec);
            printf("Bus-critical start date and time, year = %d, day = %d, hour = %d, min = %d, sec = %d\n", iy, idy, ih, mn, (int) sec);
            if (bReportFile) {
                sec = floor(sec * 1000) / 1000;
                fprintf(outfile, "btlm_start_time=%4d-%02d-%02dT%02d:%02d:%06.3fZ\n", iy, mm, dd, ih, mn, sec);
            }
        } else {
            dtime = bustime[i] - bustime[i - 1];
            if ((dtime < 0.9) || (dtime > 1.1)) {
                //sod2hms((int)secd,&ih,&mn,&sec);
                usec = yds2unix(iy, idy, secd);
                unix2ymdhms(usec, &iy, &mm, &dd, &ih, &mn, &sec);
                char tmpstr[100];
                sprintf(tmpstr, "Time difference %f at record %d, at%02d:%02d:%02d", dtime, i, ih, mn, (int) sec);
                strcpy(&errstr[ierr * 100], tmpstr);
                printf("%s\n", tmpstr);
                itimegaps++;
                ierr++;
            }
        }
    }
    fclose(infile_b);
    
    //sod2hms((int)secd,&ih,&mn,&sec);
    usec = yds2unix(iy, idy, secd);
    unix2ymdhms(usec, &iy, &mm, &dd, &ih, &mn, &sec);
    printf("Bus-critical end date and time, year = %d, day = %d, hour = %d, min = %d, sec = %d\n", iy, idy, ih, mn, (int) sec);
    printf("%d records processed\n\n", nbus);
    if (bReportFile) {
        sec = floor(sec * 1000) / 1000;
        fprintf(outfile, "btlm_stop_time=%4d-%02d-%02dT%02d:%02d:%06.3fZ\n", iy, mm, dd, ih, mn, sec);
        fprintf(outfile, "btlm_records=%d\n", nbus);
        fprintf(outfile, "btlm_errors=%d\n", ierr);
        fprintf(outfile, "btlm_timegap_records=%d\n", itimegaps);
        for (i = 0; i < ierr; i++) {
            fprintf(outfile, "btlm_err_%d=%.*s\n", i+1, 100, errstr + i * 100);
        }
    }
    
    // GPS
    if (plat >= 2 && gfile != NULL) {
    	uint8_t gpspkts[gpktsize];		
        ierr = 0;
    	itimegaps = 0;
    	free(errstr);
    	errstr = (char*) malloc(10000 * 100);
    	for (i = 0; i < ngps; i++) {
    		fread(gpspkts, 1, gpktsize, infile_g);
    		for (m = 0; m < 8; m++) cctime[m] = gpspkts[m + 6];
    		ccsds_to_yds(cctime, &ccsds_iy, &idy, &secd);
    		iy = (int16_t) ccsds_iy;
    		jd = jday((int32_t) iy, 1, idy);
    		gpstime[i] = secd + (jd - jd0)*86400;
			if (i == 0) {
				//sod2hms((int)secd,&ih,&mn,&sec);
				usec = yds2unix(iy, idy, secd);
				unix2ymdhms(usec, &iy, &mm, &dd, &ih, &mn, &sec);
				printf("GPS start date and time, year = %d, day = %d, hour = %d, min = %d, sec = %d\n", iy, idy, ih, mn, (int) sec);
				if (bReportFile) {
					sec = floor(sec * 1000) / 1000;
					fprintf(outfile, "gps_start_time=%4d-%02d-%02dT%02d:%02d:%06.3fZ\n", iy, mm, dd, ih, mn, sec);
				}
			} else {
				dtime = gpstime[i] - gpstime[i - 1];
				if ((dtime < 0.9) || (dtime > 1.1)) {
					//sod2hms((int)secd,&ih,&mn,&sec);
					usec = yds2unix(iy, idy, secd);
					unix2ymdhms(usec, &iy, &mm, &dd, &ih, &mn, &sec);
					char tmpstr[100];
					sprintf(tmpstr, "Time difference %f at record %d, at%02d:%02d:%02d", dtime, i, ih, mn, (int) sec);
					strcpy(&errstr[ierr * 100], tmpstr);
					printf("%s\n", tmpstr);
					itimegaps++;
					ierr++;
				}
			}
        }
        fclose(infile_g);
        
        //sod2hms((int)secd,&ih,&mn,&sec);
		usec = yds2unix(iy, idy, secd);
		unix2ymdhms(usec, &iy, &mm, &dd, &ih, &mn, &sec);
		printf("GPS end date and time, year = %d, day = %d, hour = %d, min = %d, sec = %d\n", iy, idy, ih, mn, (int) sec);
		printf("%d records processed\n\n", ngps);
		if (bReportFile) {
			sec = floor(sec * 1000) / 1000;
			fprintf(outfile, "gps_stop_time=%4d-%02d-%02dT%02d:%02d:%06.3fZ\n", iy, mm, dd, ih, mn, sec);
			fprintf(outfile, "gps_records=%d\n", ngps);
			fprintf(outfile, "gps_errors=%d\n", ierr);
			fprintf(outfile, "gps_timegap_records=%d\n", itimegaps);
			for (i = 0; i < ierr; i++) {
				fprintf(outfile, "gps_err_%d=%.*s\n", i+1, 100, errstr + i * 100);
			}
		}
				
    } // if (plat >= 2)

	if (bReportFile) fclose(outfile);
    free(errstr);
    return 0;
}

