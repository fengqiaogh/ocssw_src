#include "l1stat.h"
#include "l1stat_proto.h"
#include <libgen.h>

#define SET 1

char err_msg[1024];
int32 stat_status = 0; /* status of statistical checking: 0 all good,
        		   	   1 program problem, 2 statistical problem, 
		 	   	   3 both problems  */
char bad_stat_str[320]; /* summary of all mnemonics that were bad to report */

int main(int argc, char *argv[])
/*******************************************************************

   l1stat_chk

   purpose: open the SeaWiFS dataset and check the data values against
                given thresholds

   Returns type: none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               argc            I       count of input args
      char *            argv[]          I       args: [1] hdf file
                                                [2] controls file

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      L. Kumar          04-Apr-1996     Original development
      W. Robinson, GSC  31-Jul- 1998    add code to test for the amount
                                        of noise and trouble in telemetry 
                                        SDs and report if the problem
                                        is too large
      W. Robinson, GSC  29-Dec-1999     add capability for doing classified
                                        land and water comparisons to mapped
                                        values (start of changes)
      W. Robinson, GSC  17-Jan-2001     add checking of (-) time changes for 
                                        GAC files
      W. Robinson, SAIC 24 Oct 2001     add controls for the noise and 
                                  encryption % in file for reporting an error
      W. Robinson, SAIC 25 Jun 2002     add time dependence control to neg
                                  time, inst_ana, tdi and gain value checks

 *******************************************************************/ {
    int16 tilt_ranges[20][2], tilt_flags[20], dtynum, rpt_negtim;
    int32 fid, sdfid;
    int32 i, nsamp, nscans, ntilts;
    float32 nav_thresh1 = -999, nav_thresh2 = -999, l1tilt_thresh = -999,
            pct_noise_thresh, pct_encrypt_thresh;
    char fsttim[14];
    cntl1_str gn1[8], gn2[8], zero[8];
    cntl2_str l1hicnt[8], l1locnt[8];
    thr_ctl_def thr_ctl;
    int *spike_cnt;
    float *line_sd;


    bad_stat_str[0] = (char) 0;
    /*
     *  Check input arguments
     */

    if (argc != 3) {
        printf("\n******* Usage: l1stat_chk <hdf file> <control file> \n");
        printf("\t\t<hdf file> is level 1 data file\n");
        printf("\t\t<control file> is the file specifying thresholds\n");
        printf("\nOn exit, $status will be set to: \n\t0 - no problems,"
                "\n\t1 - program error, ");
        printf("\n\t2 - data problem, \n\t3 - both program and data problem\n");
        stat_status = stat_status | 1;
        stat_exit(stat_status);
    }

    printf("\n\n# Statistical Check of Dataset '%s'\n\n", argv[1]);

    /*
     *  Initialize input structures
     */

    for (i = 0; i < 8; i++) {
        gn1[i].band = gn2[i].band = zero[i].band = 0;
        gn1[i].threshold = gn2[i].threshold = zero[i].threshold = 0;
        l1hicnt[i].band = l1locnt[i].band = 0;
        l1hicnt[i].err_thresh = l1hicnt[i].cnt_thresh = 0;
        l1locnt[i].err_thresh = l1locnt[i].cnt_thresh = 0;
    }

    /*
     *   hopen will open the file and read data descriptor blocks
     *   to memory
     */

    if ((fid = Hopen(argv[1], DFACC_RDONLY, 0)) < 0) {
        printf("****** L1stat_chk: Failed on the Hopen of \n\t'%s'\n", argv[1]);
        stat_status = 1;
        stat_exit(stat_status);
    }

    /*  
     *	SDstart opens the hdf interface and initiates SD interface
     */

    if ((sdfid = SDstart(argv[1], DFACC_RDONLY)) < 0) {
        printf("******* L1stat_chk: Failure at SDstart of \n'%s'\n", argv[1]);
        Hclose(fid);
        stat_status = 1;
        stat_exit(stat_status);
    }
    /*
     * get the file start time for time-dependent checks
     */
    strncpy(fsttim, basename(argv[1]) + 1, 13);

    /*
     *  Read control file and save the information.  If the return value
     *  is -1, it means, some prolblem with opening control data file or
     *  with some control format, close HDF file, and exit with status set
     *  to 1 indicating program error.
     */

    if ((read_cntldata(argv[2], fsttim, gn1, gn2, zero, l1hicnt,
            l1locnt, &nav_thresh1, &nav_thresh2,
            &l1tilt_thresh, &pct_noise_thresh,
            &pct_encrypt_thresh, &thr_ctl, &rpt_negtim)) < 0) {
        Hclose(fid);
        stat_exit(stat_status);
    }

    /*
     *  Verify if given input file is a level 1 file
     */

    if ((l1file(sdfid, &nsamp, &nscans, &dtynum)) < 0) {
        Hclose(fid);
        printf("******* %s", err_msg);
        stat_exit(stat_status);
    }

    /*
     *  Verify Level 1 and Level 2 gain 
     */
    if ((dtynum == GAC) || (dtynum == LAC) || (dtynum == HRPT)) {
        if ((chk_gn(sdfid, gn1, gn2, dtynum, nsamp, nscans)) < 0)
            printf("\n******* %s", err_msg);
    }

    /*
     *  Verify Zero pixels 
     */
    if ((dtynum == GAC) || (dtynum == LAC) || (dtynum == HRPT)) {
        if ((chk_zero(sdfid, zero, nscans, nsamp)) < 0)
            printf("\n******* %s", err_msg);
    }

    /* 
     *  Verify High and Low count of pixels
     */
    if ((spike_cnt = (int *) calloc(8 * nscans, sizeof ( int))) == NULL) {
        printf("\n******** l1stat_chk: failure to allocate spike_cnt storage\n");
        Hclose(fid);
        stat_status = stat_status | 1;
        stat_exit(stat_status);
    }
    if ((line_sd = (float *) malloc(8 * nscans * sizeof ( float))) == NULL) {
        printf("\n******** l1stat_chk: failure to allocate line_sd storage\n");
        Hclose(fid);
        stat_status = stat_status | 1;
        stat_exit(stat_status);
    }
    if ((chk_count(sdfid, nscans, nsamp, dtynum, l1hicnt, l1locnt, spike_cnt, line_sd)) < 0)
        printf("\n******** %s", err_msg);

    /*
     *  Verify Tilt behavior to make sure tilt change lasts < 20 secs.
     */
    if ((l1tilt_thresh != -999) &&
            ((dtynum == GAC) || (dtynum == LAC) || (dtynum == HRPT))) {
        if ((chk_tilt(sdfid, dtynum, l1tilt_thresh, &ntilts, tilt_ranges,
                tilt_flags)) < 0)
            printf("\n******** %s", err_msg);
    }

    /*
     *  Verify Navigation discontinuity if requested
     */
    if ((nav_thresh1 != -999 && nav_thresh2 != -999) &&
            ((dtynum == GAC) || (dtynum == LAC) || (dtynum == HRPT))) {
        if ((chk_nav(sdfid, nscans, nav_thresh1, nav_thresh2, ntilts,
                tilt_ranges, tilt_flags, dtynum, rpt_negtim)) < 0)
            printf("\n********* %s", err_msg);
    }

    /*
     *  check file to be outside restricted time ranges
     */
    if ((thr_ctl.trng_chk_do == 1) &&
            ((dtynum == GAC) || (dtynum == LAC) || (dtynum == HRPT)))
        ck_trng(argv[1]);

    /*
     *  get other noise indicators and do a final report
     */
    rpt_noise(sdfid, dtynum, nscans, nsamp, spike_cnt, line_sd,
            pct_noise_thresh, pct_encrypt_thresh);
    /*
     *  check the instrument analog telemetry to find any unusual departures
     */
    if (dtynum == GAC) chk_inst_ana(sdfid, nscans, thr_ctl);
    /*
     *  Check the gain setting for correctness
     */
    chk_gainv(sdfid, dtynum, nscans, thr_ctl);
    /*
     *  Check the TDI settings
     */
    chk_tdiv(sdfid, dtynum, nscans, thr_ctl);

    Hclose(fid);
    stat_exit(stat_status);

    return 0;
}

/****************************************************************************
   read_cntldata

   purpose: Opens given control file and reads in the mnemonics and thresholds
                given.  If any error occurs, prints out an error message,
                sets stat_status field and returns a negative 1.	

   Returns type: integer

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char  *           cntl_file        I     	control file path/name 
      char  *           fsttim           I      input file start time, format
                                                YYYYDDDHHMMSS
      cntl1_str *	gn1		 O	structure to hold gain1 control
                                                data
      cntl1_str *	gn1		 O	structure to hold gain2 control
                                                data
      cntl1_str *	zero		 O      struct to hold Zero pix controls
      cntl2_str *	l1hicnt		 O	struct to hold hight count
                                                control infomration
      cntl2_str *	l1locnt		 O	struct to hold low count 
                                                control data
      float32 * 	nav_thresh1 	 O	threshold 1 of L1NAVDISC
      float32 *		nav_thresh2 	 O	threshold 2 of L1NAVDISC	
      float32 *		l1tilt_thresh	 O	threshold for L1TILT
      float32 *         pct_noise_thresh O      limit for % of lines with 
                                                noise for reporting
      float32 *         pct_encrypt_thresh O    limit for % of lines with 
                                                encryption for reporting
      thr_ctl_def *     thr_ctl          O      threshold controls for the 
                                                instrument analog values, the
                                                gain check values
      int16 *           rpt_negtim       O      if = 1, fail file if negative
                                                time step happens
       
   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      L. Kumar          04-Apr-1996    Original development
      W. Robinson, SAIC 25-Jun-2002    get the controls for time-dependent 
                                       group checks of neg delta time, inst_ana
                                       tdi and gain checks
      W. Robinson, SAIC 25-Mar-2010    set up bnd as a null terminated 
                                       array for safety

 ****************************************************************************/
int32 read_cntldata(char *cntl_file, char *fsttim, cntl1_str *gn1,
        cntl1_str *gn2, cntl1_str *zero, cntl2_str *l1hicnt,
        cntl2_str *l1locnt, float32 *nav_thresh1,
        float32 *nav_thresh2, float32 *l1tilt_thresh,
        float32 *pct_noise_thresh, float32 *pct_encrypt_thresh,
        thr_ctl_def *thr_ctl, int16 *rpt_negtim) {
    FILE *fid;
    char line[501], str[25], bnd[2], tstr[14];
    int32 i, band, nflds = 0, ival;
    float32 threshold, cnt_thresh, err_thresh, f1, f2, f3;

    printf("\n# Reading control file '%s'\n\n", cntl_file);

    /*
     *  Open the file
     */
    if ((fid = fopen(cntl_file, "r")) == NULL) {
        printf("********read_cntldata: Unable to open control file:\n\t'%s'\n",
                cntl_file);
        stat_status = stat_status | 1;
        return FAIL;
    }

    /*
     *  Initialize all input structrues
     */
    thr_ctl->trng_chk_do = 0;
    for (i = 0; i < 8; i++) {
        gn1[i].band = gn2[i].band = 0;
        gn1[i].threshold = gn2[i].threshold = 0;
        thr_ctl->gainv_chk_do[i] = 0;
    }
    for (i = 0; i < 32; i++)
        thr_ctl->inst_ana_do[i] = 0;

    thr_ctl->rpt_tdi_vchk = 1; /* default is to report occurences */
    thr_ctl->rpt_gainv_chk = 1; /* default is to report occurences */
    thr_ctl->rpt_inst_ana = 1; /* default is to report occurences */
    *rpt_negtim = 1; /* default is to report neg time occurences */

    while (fgets(line, 500, fid) != NULL) {
        /*  Ignore beginning blanks and tabs */
        for (i = 0; i < 500 && (line[i] == ' ' || line[i] == '\t'); i++);

        /* If first character is '#' sign, then treat the line as comment */
        if (i < 500 && line[i] == '#') {
            printf("%s", line);
            continue;
        }

        /* If not comment, check if it is GAIN 1 or GAIN 2 info */
        printf("#%s", line);
        if (strncmp(&line[i], "L1GAIN", 6) == 0) {
            if ((nflds = sscanf(line, "%s %f", str, &threshold)) != 2) {
                printf("\n*********read_cntldata: expecting 2 fields");
                printf(" but read %d field(s)", nflds);
                stat_status = stat_status | 1;
                return FAIL;
            }
            bnd[0] = str[7];
            bnd[1] = 0;
            band = atoi(bnd);
            if (band < 1 || band > 8) {
                printf("********read_cntldata: Error in band number.");
                printf("  Band # read = %d", band);
                stat_status = stat_status | 1;
                return FAIL;
            }
            if (str[6] == '1') {
                gn1[band - 1].band = 1;
                gn1[band - 1].threshold = threshold;
            } else if (str[6] == '2') {
                gn2[band - 1].band = 1;
                gn2[band - 1].threshold = threshold;
            }
        } else if (strncmp(&line[i], "L1ZERO", 6) == 0) {
            if ((nflds = sscanf(line, "%s %f", str, &threshold)) != 2) {
                printf("*********read_cntldata: expecting 2 fields");
                printf(" but read %d field(s)", nflds);
                stat_status = stat_status | 1;
                return FAIL;
            }
            bnd[0] = str[6];
            band = atoi(bnd);
            if (band < 1 || band > 8) {
                printf("********read_cntldata: Error in band number.");
                printf("  Band # read = %d", band);
                stat_status = stat_status | 1;
                return FAIL;
            }
            zero[band - 1].band = 1;
            zero[band - 1].threshold = threshold;
        } else if (strncmp(&line[i], "L1HICOUNT", 9) == 0) {
            if ((nflds = sscanf(line, "%s %f %f", str, &err_thresh, &cnt_thresh))
                    != 3) {
                printf("*********read_cntldata: expecting 3 fields");
                printf(" but read %d field(s)", nflds);
                stat_status = stat_status | 1;
                return FAIL;
            }
            bnd[0] = str[9];
            band = atoi(bnd);
            if (band < 1 || band > 8) {
                printf("*********read_cntldata: Error in band number.");
                printf("  Band # read = %d", band);
                stat_status = stat_status | 1;
                return FAIL;
            }
            l1hicnt[band - 1].band = 1;
            l1hicnt[band - 1].err_thresh = err_thresh;
            l1hicnt[band - 1].cnt_thresh = cnt_thresh;
        } else if (strncmp(&line[i], "L1LOWCOUNT", 10) == 0) {
            if ((nflds = sscanf(line, "%s %f %f", str, &err_thresh, &cnt_thresh))
                    != 3) {
                printf("*********read_cntldata: expecting 3 fields");
                printf(" but read %d field(s)", nflds);
                stat_status = stat_status | 1;
                return FAIL;
            }
            bnd[0] = str[10];
            band = atoi(bnd);
            if (band < 1 || band > 8) {
                printf("*********read_cntldata: Error in band number.");
                printf("  Band # read = %d", band);
                stat_status = stat_status | 1;
                return FAIL;
            }
            l1locnt[band - 1].band = 1;
            l1locnt[band - 1].err_thresh = err_thresh;
            l1locnt[band - 1].cnt_thresh = cnt_thresh;
        } else if (strncmp(&line[i], "L1TILT", 6) == 0) {
            if ((nflds = sscanf(line, "%s %f", str, l1tilt_thresh)) != 2) {
                printf("*********read_cntldata: expecting 2 fields");
                printf(" but read %d field(s)", nflds);
                stat_status = stat_status | 1;
                return FAIL;
            }
        }            /*
      *   WDR new check of time limit on neg time check
      */
        else if (strncmp(&line[i], "TLIM_NEGTIM", 11) == 0) {
            if ((nflds = sscanf(line, "%s %s", str, tstr)) != 2) {
                printf("***** read_cntldata: expecting 2 fields for TLIM_NEGTIM\n");
                printf(" but read %d field(s)", nflds);
                stat_status = stat_status | 1;
                return FAIL;
            }
            *rpt_negtim = 1;
            if (strncmp(fsttim, tstr, 13) < 0) *rpt_negtim = 0;
        } else if (strncmp(&line[i], "L1NAVDISC", 9) == 0) {
            if ((nflds = sscanf(line, "%s %f %f", str, nav_thresh1, nav_thresh2))
                    != 3) {
                printf("*********read_cntldata: expecting 3 fields");
                printf(" but read %d field(s)", nflds);
                stat_status = stat_status | 1;
                return FAIL;
            }
        } else if (strncmp(&line[i], "NOISE", 5) == 0) {
            if ((nflds = sscanf(line, "%s %f", str, pct_noise_thresh)) != 2) {
                printf("*********read_cntldata: NOISE: expecting 2 fields");
                printf(" but read %d field(s)", nflds);
                stat_status = stat_status | 1;
                return FAIL;
            }
        } else if (strncmp(&line[i], "ENCRYPT", 7) == 0) {
            if ((nflds = sscanf(line, "%s %f", str, pct_encrypt_thresh)) != 2) {
                printf("*********read_cntldata: ENCRYPT: expecting 2 fields");
                printf(" but read %d field(s)", nflds);
                stat_status = stat_status | 1;
                return FAIL;
            }
        } else if (strncmp(&line[i], "TRNG_CHK", 8) == 0) {
            if ((nflds = sscanf(line, "%s", str)) != 1) {
                printf("*********read_cntldata: TRNG_CHK: expecting 1 field");
                printf(" but read %d field(s)", nflds);
                stat_status = stat_status | 1;
                return FAIL;
            }
            thr_ctl->trng_chk_do = 1;
        }            /*
      *  WDR new time limit check for inst_ana
      */
        else if (strncmp(&line[i], "TLIM_INST_ANA", 13) == 0) {
            if ((nflds = sscanf(line, "%s %s", str, tstr)) != 2) {
                printf("***** read_cntldata: expecting 2 fields for TLIM_INST_ANA\n");
                printf(" but read %d field(s)", nflds);
                stat_status = stat_status | 1;
                return FAIL;
            }
            thr_ctl->rpt_inst_ana = 1;
            if (strncmp(fsttim, tstr, 13) < 0) thr_ctl->rpt_inst_ana = 0;
        } else if (strncmp(&line[i], "INST_ANA", 8) == 0) {
            if ((nflds = sscanf(line,
                    "INST_ANA%2d %f %f %f", &ival, &f1, &f2, &f3)) != 4) {
                printf("********read_cntldata: Error in INST_ANA line decode\n");
                printf("# fields read: %d\n", nflds);
                printf("Line is:\n%s\n\n", line);
                stat_status = stat_status | 1;
                return FAIL;
            }
            /*
             *  for the inst_ana, get the thresholds
             */
            if ((ival < 1) || (ival > 32)) {
                printf("********read_cntldata: Error in INST_ANA item #\n");
                printf("It must be from 1 - 32 and the entered value was %d\n",
                        ival);
                stat_status = stat_status | 1;
                return FAIL;
            }
            thr_ctl->inst_ana_do[ ival - 1 ] = 1;
            thr_ctl->inst_ana_lo[ ival - 1 ] = f1;
            thr_ctl->inst_ana_hi[ ival - 1 ] = f2;
            thr_ctl->inst_ana_pct[ ival - 1 ] = f3;
        }            /*
      *  WDR new time limit check for gain value
      */
        else if (strncmp(&line[i], "TLIM_GAINV_CHK", 14) == 0) {
            if ((nflds = sscanf(line, "%s %s", str, tstr)) != 2) {
                printf("***** read_cntldata: expecting 2 fields for TLIM_GAINV_CHK\n");
                printf(" but read %d field(s)", nflds);
                stat_status = stat_status | 1;
                return FAIL;
            }
            thr_ctl->rpt_gainv_chk = 1;
            if (strncmp(fsttim, tstr, 13) < 0) thr_ctl->rpt_gainv_chk = 0;
        } else if (strncmp(&line[i], "GAINV_CHK", 9) == 0) {
            if ((nflds = sscanf(line, "GAINV_CHK%1d %f", &ival, &f1)) != 2) {
                printf("********read_cntldata: Error in GAINV_CHK line decode\n");
                printf("# fields read: %d\n", nflds);
                printf("Line is:\n%s\n\n", line);
                stat_status = stat_status | 1;
                return FAIL;
            }
            /*
             *  for the gainv_chk, get the % error acceptable
             */
            if ((ival < 1) || (ival > 8)) {
                printf("********read_cntldata: Error in GAINV_CHK item #\n");
                printf("It must be from 1 - 8 and the entered value was %d\n",
                        ival);
                stat_status = stat_status | 1;
                return FAIL;
            }
            thr_ctl->gainv_chk_do[ ival - 1 ] = 1;
            thr_ctl->gainv_chk_pct[ ival - 1 ] = f1;
        }            /*
      *  WDR new time limit check for gain value
      */
        else if (strncmp(&line[i], "TLIM_TDIV_CHK", 13) == 0) {
            if ((nflds = sscanf(line, "%s %s", str, tstr)) != 2) {
                printf("***** read_cntldata: expecting 2 fields for TLIM_TDIV_CHK\n");
                printf(" but read %d field(s)", nflds);
                stat_status = stat_status | 1;
                return FAIL;
            }
            thr_ctl->rpt_tdi_vchk = 1;
            if (strncmp(fsttim, tstr, 13) < 0) thr_ctl->rpt_tdi_vchk = 0;
        } else if (strncmp(&line[i], "TDIV_CHK", 8) == 0) {
            if ((nflds = sscanf(line, "TDIV_CHK%1d %f", &ival, &f1)) != 2) {
                printf("********read_cntldata: Error in TDIV_CHK line decode\n");
                printf("# fields read: %d\n", nflds);
                printf("Line is:\n%s\n\n", line);
                stat_status = stat_status | 1;
                return FAIL;
            }
            /*
             *  for the tdiv_chk, get the % error acceptable
             */
            if ((ival < 1) || (ival > 8)) {
                printf("********read_cntldata: Error in TDIV_CHK item #\n");
                printf("It must be from 1 - 8 and the entered value was %d\n",
                        ival);
                stat_status = stat_status | 1;
                return FAIL;
            }
            thr_ctl->tdiv_chk_do[ ival - 1 ] = 1;
            thr_ctl->tdiv_chk_pct[ ival - 1 ] = f1;
        } else {
            printf("*********read_cntldata: Cannot recognize the mnemonic read");
            stat_status = stat_status | 1;
            return FAIL;
        }
    }
    return SUCCEED;
}

/****************************************************************************
   l1file

   purpose: Verifies if the given input file is level 1 data file

   Returns type: integer

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             sdfid           I       SD interface ID
      int32 *		nsamp		O	number of samples/pixels
      int32 *		nscans		O	number of scan lines
      int16 *		dtynum		O  	data type number 
                                        representation: GAC  0, LAC  1,
                                        HRPT 2, LUN  3, SOL  4, IGC  5, TDI  6

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      L. Kumar          04-Apr-1996     Original development

 ****************************************************************************/
int32 l1file(int32 sdfid, int32 *nsamp, int32 *nscans, int16 *dtynum) {
    char title[1024];
    char dtype[500];

    /*
     *  Read title to verify if the given input file is level 1A data file
     */

    if ((rdattr(sdfid, TITLE, &title)) < 0)
        return FAIL;

    if (strcmp(title, "SeaWiFS Level-1A Data") != 0) {
        sprintf(err_msg, "l1file: Data file is not level 1A file");
        stat_status = stat_status | 1;
        return FAIL;
    }

    if ((rdattr(sdfid, NSAMP, (VOIDP *) nsamp)) < 0)
        return FAIL;

    if ((rdattr(sdfid, NSCANS, (VOIDP *) nscans)) < 0)
        return FAIL;

    if ((rdattr(sdfid, DTYPE, dtype)) < 0)
        return FAIL;
    if (strncmp(dtype, "GAC", 3) == 0) {
        *dtynum = (int16) 0;
    } else if (strncmp(dtype, "LAC", 3) == 0) {
        *dtynum = (int16) 1;
    } else if (strncmp(dtype, "HRPT", 4) == 0) {
        *dtynum = (int16) 2;
    } else if (strncmp(dtype, "LUN", 3) == 0) {
        *dtynum = (int16) 3;
    } else if (strncmp(dtype, "SOL", 3) == 0) {
        *dtynum = (int16) 4;
    } else if (strncmp(dtype, "IGC", 3) == 0) {
        *dtynum = (int16) 5;
    } else if (strncmp(dtype, "TDI", 3) == 0) {
        *dtynum = (int16) 6;
    } else {
        printf("Data type of '%s' is unknown.  Exiting\n", dtype);
        return FAIL;
    }

    return SUCCEED;
}

/****************************************************************************
   chk_gn 

   purpose: Verifies gain 1 and 2 saturated percentage against the given
                thresholds for requested bands

   Returns type: integer 

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             sdfid           I       SD interface ID
      cntl1_str	*	gn1 		I       structure containing gain1
                                                band numbers and thresholds
      cntl1_str  *	gn2		I  	structure containing gain2
                                                band numbers and thresholds
      int16 	   	dtynum		I	data type number
      int32		nsamp		I 	number of samples
      int32	 *	nscans		I 	number of scanlines	
      
   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      L. Kumar          04-Apr-1996     Original development
      W. Robinson       31-Dec-1997     switch to checking % of good 
                                        gains for image data

 ****************************************************************************/
int32 chk_gn(int32 sdfid, cntl1_str *gn1, cntl1_str *gn2, int16 dtynum,
        int32 nsamp, int32 nscans) {
    int32 i, gn1sat[8], gn2sat[8], gn1unsat[8], gn2unsat[8];
    int32 failcode = 0;
    float64 gn1_val[8], gn2_val[8];
    float pct_good_gain;
    char str[12];

    printf("\n#Checking Level 1 and Level 2 gains ....");

    /*
     *  Following few lines Initializes value buffers 
     *  Reads gain1 and gain 2 saturated and non-saturated pixel values from
     *  global attributes of input data file
     */
    for (i = 0; i < 8; i++)
        gn1_val[i] = gn2_val[i] = 0;

    if ((rdattr(sdfid, GN1SAT, (VOIDP *) gn1sat)) < 0)
        return FAIL;

    if ((rdattr(sdfid, GN2SAT, (VOIDP *) gn2sat)) < 0)
        return FAIL;

    if ((rdattr(sdfid, GN1UNSAT, (VOIDP *) gn1unsat)) < 0)
        return FAIL;

    if ((rdattr(sdfid, GN2UNSAT, (VOIDP *) gn2unsat)) < 0)
        return FAIL;


    /* 
     * If gain 1 is > 0, calculate saturated percentage
     * else, set value to -1, indicating success
     */
    for (i = 0; i < 8; i++) {
        if ((gn1sat[i] + gn1unsat[i]) > 0)
            gn1_val[i] = ((float64) gn1sat[i] / (gn1sat[i] + gn1unsat[i])) * 100;
        else
            gn1_val[i] = -1;

        if ((gn2sat[i] + gn2unsat[i]) > 0)
            gn2_val[i] = ((float64) gn2sat[i] / (gn2sat[i] + gn2unsat[i]))*100;
        else
            gn2_val[i] = -1;
    }

    /*
     *  Check to see if all the saturated and unsaturated pixel count matches
     *  total number of pixels for each band
     */
    if ((dtynum == GAC) || (dtynum == LAC) || (dtynum == HRPT))
        for (i = 0; i < 8; i++) {
            pct_good_gain =
                    (float) (gn1sat[i] + gn1unsat[i] + gn2sat[i] + gn2unsat[i]) /
                    (float) (nsamp * nscans) * 100.;
            if (pct_good_gain < 100.) {
                if (pct_good_gain >= 90.) /* only note the underperforming gain */ {
                    printf("\nNote that percentage of good gains (%f) for band %d is < 100%%", pct_good_gain, i);
                } else {
                    printf("\n********chk_gn: Percentage of good gains (%f) for band %d is < 90 %%, an ERROR condition", pct_good_gain, i);
                    stat_status = stat_status | 2;
                    sprintf(str, "L1GAIN%1d ", (i + 1));
                    if (strlen(bad_stat_str) <= 300)
                        strcat(bad_stat_str, str);
                }
            }
        }

    /*  
     *  Print headers for l1 and l2 gains output messages
     */

    printf("\n\n#Mnemonic  Code\t  Value\t   --\tError_threshold");
    printf("\n#-----------------------------------------------\n");
    /* 
     * If gain 1 is > 0, calculate saturated percentage
     * else, set value to -1, indicating success
     */
    for (i = 0; i < 8; i++) {
        failcode = 0;
        if (gn1[i].band == 1) {
            if (gn1_val[i] != -1 && gn1_val[i] > gn1[i].threshold) {
                failcode = 1;
                stat_status = stat_status | 2;
                sprintf(str, "L1GAIN1%1d ", (i + 1));
                if (strlen(bad_stat_str) <= 300)
                    strcat(bad_stat_str, str);
            }
            printf("\nL1GAIN1%d   %d\t%f   --\t%f", i + 1, failcode, gn1_val[i],
                    gn1[i].threshold);
        }
    }

    printf("\n");
    for (i = 0; i < 8; i++) {
        failcode = 0;
        if (gn2[i].band == 1) {
            if (gn2_val[i] != -1 && gn2_val[i] > gn2[i].threshold) {
                failcode = 1;
                stat_status = stat_status | 2;
                sprintf(str, "L1GAIN2%1d ", (i + 1));
                if (strlen(bad_stat_str) <= 300)
                    strcat(bad_stat_str, str);
            }
            printf("\nL1GAIN2%d   %d\t%f   --\t%f", i + 1, failcode, gn2_val[i],
                    gn2[i].threshold);
        }
    }

#ifdef DEBUG
    for (i = 0; i < 8; i++)
        printf("\n gn1sat = %d\t gn1unsat = %d\t gn1_val = %f\n", gn1sat[i],
            gn1unsat[i], gn1_val[i]);
    for (i = 0; i < 8; i++)
        printf("\n gn2sat = %d\t gn2unsat = %d\t gn2_val = %f\n", gn2sat[i],
            gn2unsat[i], gn2_val[i]);
#endif /* DEBUG */

    return SUCCEED;
}

/****************************************************************************
   chk_zero

   purpose: Verifies Zero pixels percentage against the given
                threshold for requested band/s

   Returns type: integer 

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             sdfid           I       SD interface ID
      cntl1_str  *       zero_str        I       structure containing zero
                                                pixel band nos., and thresholds
      int32		nscans		I	Number of scan lines
      int32		nsamp		I	Number of pixels

   Modification history:

      Programmer        Date            Description of change
      ----------        ----------      ---------------------
      L. Kumar          04-Apr-1996     Original development

 ****************************************************************************/
int32 chk_zero(int32 sdfid, cntl1_str *zero_str, int32 nscans, int32 nsamp) {
    int32 i, zero_pix[8];
    int32 failcode = 0;
    float64 zero_val[8];
    char str[12];

    printf("\n\n#Checking Zero pixels percentage ....");

    /*
     *  Following few lines Initializes value buffers 
     *  Reads zero pixel values from global attribute 'Zero Pixels'
     *  of input data file
     */
    for (i = 0; i < 8; i++)
        zero_val[i] = 0;

    if ((rdattr(sdfid, ZEROPIX, (VOIDP *) zero_pix)) < 0)
        return FAIL;

    for (i = 0; i < 8; i++)
        zero_val[i] = ((float64) zero_pix[i] / (nscans * nsamp)) * 100.0;

    /*  
     *  Print headers for l1 and l2 gains output messages
     */

    printf("\n\n#Mnemonic  Code\t  Value\t   --\tError_threshold");
    printf("\n#-----------------------------------------------\n");

    for (i = 0; i < 8; i++) {
        failcode = 0;
        if (zero_str[i].band == 1) {
            if (zero_val[i] > zero_str[i].threshold) {
                failcode = 1;
                stat_status = stat_status | 2;
                sprintf(str, "L1ZERO%1d ", (i + 1));
                if (strlen(bad_stat_str) <= 300)
                    strcat(bad_stat_str, str);
            }
            printf("\nL1ZERO%d   %d\t%f   --\t%f", i + 1, failcode, zero_val[i],
                    zero_str[i].threshold);
        }
    }

#ifdef DEBUG
    for (i = 0; i < 8; i++)
        printf("\n zero_pix = %d\t zero_val = %f\n", zero_pix[i], zero_val[i]);
#endif /* DEBUG */

    return SUCCEED;
}

/****************************************************************************
   chk_count

   purpose: Verifies pixels high and low count percentages against the given
                error threshold and count threshold

   Returns type: integer 

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             sdfid           I       SD interface ID
      int32		nscans		I	Number of scan lines
      int32		nsamp		I	Number of pixels
      int16             dtynum          I       data type number
      cntl2_str  *      l1hicnt         I       struct containing high count 
                                                pixels error & count thresholds 
      cntl2_str  *      l1locnt         I       struct containing low count 
                                                pixels error & count thresholds 
      int *             spike_cnt      I/O      # lines X # bands array with
                                                the spike noise count per line
      float *           line_sd        I/O      # lines X # bands std 
                                                deviation for each line, band

   Modification history:

      Programmer        Date            Description of change
      ----------        ----------      ---------------------
      L. Kumar          04-Apr-1996     Original development
      W. Robinson, GSC  31-Jul-1998     perform an analysis of noise while the 
                                        data is available
      W. Robinson, GSC  29-oct-2001     flag only the GAC, LAC, HRPT for the 
                                        LOCOUNT

 ****************************************************************************/
int32 chk_count(int32 sdfid, int32 nscans, int32 nsamp, int16 dtynum,
        cntl2_str *l1hicnt, cntl2_str *l1locnt, int *spike_cnt, float *line_sd) {
    int16 *i16buf;
    int32 i, rec, bnd, failcode = 0;
    int32 nbands = 8, nrec, rdrecs = 256, recsleft;
    int32 start[3], edge[3];
    int32 hicnt[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    int32 lowcnt[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    float32 value[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    int32_t cnt_coin_jmp[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    int32_t jmp_hist[1024];
    char str[12];

    int DETAIL_PRT = 0; /* set to 1 to print noise detail info 0 to not print */

    for (i = 0; i < 1024; i++)
        jmp_hist[i] = 0;
    nrec = rdrecs;

    printf("\n\n#Checking Low and High count percentages ....");

    /*
     *  Allocate buffer space for reading l1a_data from given input hdf file
     */
    if ((i16buf = (int16 *) calloc(nrec * nsamp * nbands, sizeof (int16))) == NULL) {
        sprintf(err_msg,
                "chk_count: cannot allocate memory for reading l1a_data");
        return FAIL;
    }

    /*
     * It has been observed that it takes quite some time to read a large chunk
     * of data from a HDF file at once.  In order to make it time efficient
     * the l1a_data is being read in smaller chunks containing "nrec" scans.  
     * (see the defn. of nrec above in the variable declaration section)
     * Following assinment stmts. initializes start and end dimensions of the
     * first data slab that will be read.  The for loop reads in data and
     * calls get_hicnt and get_lowcnt to accumulate high and low counts of 
     * level 1 data.
     */

    start[0] = start[1] = start[2] = 0;
    edge[0] = nrec;
    edge[1] = nsamp;
    edge[2] = nbands;

    for (rec = 0; rec < nscans; rec += rdrecs) {
        if ((recsleft = nscans - rec) < rdrecs)
            nrec = recsleft;
        if (nrec > 0) {
            start[0] = rec;
            edge[0] = nrec;
            if ((rdslice(sdfid, L1ADATA, start, edge, (VOIDP) i16buf)) < 0) {
                stat_status = stat_status | 1;
                return FAIL;
            }
            get_hicnt(nrec, nsamp, nbands, i16buf, l1hicnt, hicnt);
            get_lowcnt(nrec, nsamp, nbands, i16buf, l1locnt, lowcnt);
            /*
             *  for the lines read, collect the noise estimates
             */
            anal_noise(rec, nrec, nscans, nsamp, nbands, i16buf, spike_cnt, line_sd, cnt_coin_jmp, jmp_hist);
        }
    }
    /*
     *  report the grand total info for the co-incident jumps and the histogram
     */

    if (DETAIL_PRT) {
        printf("Grand total of co-incident jumps\n");
        printf("Total # pixels: %d\n", nscans * nsamp);
        for (rec = 0; rec < 8; rec++)
            printf("#co-incidences: %d,  count: %d\n",
                (rec + 1), *(cnt_coin_jmp + rec));

        printf("\nhistogram (size 1024) of jump frequency");
        for (rec = 0; rec < 1024; rec++) {
            printf("%10d", *(jmp_hist + rec));
            if (rec % 8 == 7) printf("\n");
        }

        /*
         * temporarily report the spike and sd info here
         */

        printf("\n\nline-by-line jump count   Standard deviation\n");
        printf("line   1     2     3     4     5     6     7     8       1       2       3       4       5       6       7       8\n");
        for (rec = 0; rec < nscans; rec++) {
            printf("%5d", rec);
            for (bnd = 0; bnd < nbands; bnd++) {
                printf(" %4d", *(spike_cnt + bnd + nbands * rec));
            }
            for (bnd = 0; bnd < nbands; bnd++) {
                printf(" %7.2f", *(line_sd + bnd + nbands * rec));
            }
            printf("\n");
        }
    }

    /*
     *  Following for loops calculates the percentage of high and low counts
     *  verifies against the given hight and low error thresholds and outputs
     *  the message
     */

    printf("\n\n#Mnemonic  Code\t  Value\t   --\tErr_thresh\tCount_thresh");
    printf("\n#-----------------------------------------------------------\n");

    for (i = 0; i < nbands; i++) {
        failcode = 0;
        if (l1hicnt[i].band == 1) {
            value[i] = (hicnt[i] / (nscans * nsamp * 1.0)) * 100.0;
            if (value[i] > l1hicnt[i].err_thresh) {
                failcode = 1;
                stat_status = stat_status | 2;
                sprintf(str, "L1HICOUNT%1d ", (i + 1));
                if (strlen(bad_stat_str) <= 300)
                    strcat(bad_stat_str, str);
            }
            printf("\nL1HICOUNT%d   %d\t%f   --\t%f\t%f",
                    i + 1, failcode, value[i], l1hicnt[i].err_thresh,
                    l1hicnt[i].cnt_thresh);
        }
    }
    printf("\n");
    for (i = 0; i < nbands; i++) {
        failcode = 0;
        if (l1locnt[i].band == 1) {
            value[i] = (lowcnt[i] / (nscans * nsamp * 1.0))*100.00;
            /* only flag if it is an image dataset for the low count */
            if (((dtynum == GAC) || (dtynum == LAC) || (dtynum == HRPT)) &&
                    (value[i] > l1locnt[i].err_thresh)) {
                failcode = 1;
                stat_status = stat_status | 2;
                sprintf(str, "L1LOWCOUNT%1d ", (i + 1));
                if (strlen(bad_stat_str) <= 300)
                    strcat(bad_stat_str, str);
            }
            if ((dtynum == GAC) || (dtynum == LAC) || (dtynum == HRPT)) {
                printf("\nL1LOWCOUNT%d   %d\t%f   --\t%f\t%f",
                        i + 1, failcode, value[i], l1locnt[i].err_thresh,
                        l1locnt[i].cnt_thresh);
            } else {
                printf("\nL1LOWCOUNT%d  N/A\t%f   --\t%f\t%f",
                        i + 1, value[i], l1locnt[i].err_thresh,
                        l1locnt[i].cnt_thresh);
            }
        }
    }


    free(i16buf);

    return SUCCEED;
}

/****************************************************************************
   get_hicnt

   purpose: Accumulates high count pixels for each requested band

   Returns type: void

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             nrec            I       Number of records
      int32             nscans          I       Number of scan lines
      int32             nsamp           I       Number of pixels
      int16 	 *	databuf	     	I	data buffer
      cntl2_str  *      l1hicnt         I       struct containing high count
                                                pixels error & count thresholds
      cntl2_str  *      hicnt         	I       Array of 8 to hold high counts
                                                of each band

   Modification history:

      Programmer        Date            Description of change
      ----------        ----------      ---------------------
      L. Kumar          04-Apr-1996     Original development

 ****************************************************************************/
void
#ifdef PROTOTYPE
get_hicnt(int32 nrec, int32 nsamp, int32 nbands, int16 *databuf,
        cntl2_str *l1hicnt, int32 *hicnt)
#else
get_hicnt(nrec, nsamp, nbands, databuf, l1hicnt, hicnt)
int32 nrec, nsamp, nbands, *hicnt;
int16 *databuf;
cntl2_str *l1hicnt;
#endif
{
    int32 i;

    if (l1hicnt[BAND1].band)
        for (i = 0; i < nrec * nsamp * nbands; i += 8)
            if (databuf[i] > l1hicnt[BAND1].cnt_thresh)
                hicnt[BAND1]++;

    if (l1hicnt[BAND2].band)
        for (i = 1; i < nrec * nsamp * nbands; i += 8)
            if (databuf[i] > l1hicnt[BAND2].cnt_thresh)
                hicnt[BAND2]++;

    if (l1hicnt[BAND3].band)
        for (i = 2; i < nrec * nsamp * nbands; i += 8)
            if (databuf[i] > l1hicnt[BAND3].cnt_thresh)
                hicnt[BAND3]++;

    if (l1hicnt[BAND4].band)
        for (i = 3; i < nrec * nsamp * nbands; i += 8)
            if (databuf[i] > l1hicnt[BAND4].cnt_thresh)
                hicnt[BAND4]++;

    if (l1hicnt[BAND5].band)
        for (i = 4; i < nrec * nsamp * nbands; i += 8)
            if (databuf[i] > l1hicnt[BAND5].cnt_thresh)
                hicnt[BAND5]++;

    if (l1hicnt[BAND6].band)
        for (i = 5; i < nrec * nsamp * nbands; i += 8)
            if (databuf[i] > l1hicnt[BAND6].cnt_thresh)
                hicnt[BAND6]++;

    if (l1hicnt[BAND7].band)
        for (i = 6; i < nrec * nsamp * nbands; i += 8)
            if (databuf[i] > l1hicnt[BAND7].cnt_thresh)
                hicnt[BAND7]++;

    if (l1hicnt[BAND8].band)
        for (i = 7; i < nrec * nsamp * nbands; i += 8)
            if (databuf[i] > l1hicnt[BAND8].cnt_thresh)
                hicnt[BAND8]++;

}

/****************************************************************************
   get_lowcnt 

   purpose: Accumulates low count pixels for each band

   Returns type: void

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             nrec            I       Number of records
      int32             nscans          I       Number of scan lines
      int32             nsamp           I       Number of pixels
      int16 	 *	databuf	     	I	data buffer
      cntl2_str  *      l1locnt         I       struct containing low count
                                                pixels error & count thresholds
      cntl2_str  *      lowcnt        	I       Array of 8 to hold high counts
                                                of each band

   Modification history:

      Programmer        Date            Description of change
      ----------        ----------      ---------------------
      L. Kumar          04-Apr-1996     Original development

 ****************************************************************************/
void
#ifdef PROTOTYPE
get_lowcnt(int32 nrec, int32 nsamp, int32 nbands, int16 *databuf,
        cntl2_str *l1locnt, int32 *lowcnt)
#else
get_lowcnt(nrec, nsamp, nbands, databuf, l1locnt, lowcnt)
int32 nrec, nsamp, nbands, *lowcnt;
int16 *databuf;
cntl2_str *l1locnt;
#endif
{
    int32 i;

    if (l1locnt[BAND1].band)
        for (i = 0; i < nrec * nsamp * nbands; i += 8)
            if (databuf[i] < l1locnt[BAND1].cnt_thresh)
                lowcnt[BAND1]++;

    if (l1locnt[BAND2].band)
        for (i = 1; i < nrec * nsamp * nbands; i += 8)
            if (databuf[i] < l1locnt[BAND2].cnt_thresh)
                lowcnt[BAND2]++;

    if (l1locnt[BAND3].band)
        for (i = 2; i < nrec * nsamp * nbands; i += 8)
            if (databuf[i] < l1locnt[BAND3].cnt_thresh)
                lowcnt[BAND3]++;

    if (l1locnt[BAND4].band)
        for (i = 3; i < nrec * nsamp * nbands; i += 8)
            if (databuf[i] < l1locnt[BAND4].cnt_thresh)
                lowcnt[BAND4]++;

    if (l1locnt[BAND5].band)
        for (i = 4; i < nrec * nsamp * nbands; i += 8)
            if (databuf[i] < l1locnt[BAND5].cnt_thresh)
                lowcnt[BAND5]++;

    if (l1locnt[BAND6].band)
        for (i = 5; i < nrec * nsamp * nbands; i += 8)
            if (databuf[i] < l1locnt[BAND6].cnt_thresh)
                lowcnt[BAND6]++;

    if (l1locnt[BAND7].band)
        for (i = 6; i < nrec * nsamp * nbands; i += 8)
            if (databuf[i] < l1locnt[BAND7].cnt_thresh)
                lowcnt[BAND7]++;

    if (l1locnt[BAND8].band)
        for (i = 7; i < nrec * nsamp * nbands; i += 8)
            if (databuf[i] < l1locnt[BAND8].cnt_thresh)
                lowcnt[BAND8]++;
}

/****************************************************************************
   chk_nav 

   purpose: Verifies navigation discontinuity.  Following gives the algorithm:
                Take center pixel of every line, find the latitude, longitude
                via navigation.  Convert it to a vector from earth center and 
                compure the angular difference between this vector and the one 
                for the previous line.  At first, exclude regions of tilt 
                change as the difference would be larger.  the angular speed
                will be computed using the angular difference and the 
                line time tags.  If the speed is within tolerences,
                the test passes.  This is different from the previous
                km. difference because now, the speed is the same for GAC 
                and LAC.

   Returns type: int32

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             sdfid           I       SD interface id
      int32		nscans		I	Number of scanlines
      float32           nav_thresh1     I       input threshold 1
      float32           nav_thresh2     I       input threshold 2
      int32 *		ntilts		I	number of tilt states
      int16 		tilt_ranges[20[2] I    	scan line ranges of tilt states
      int16		tilt_flags[20]	I	flags corresponding to tilt
                                                states
      int16             dtynum          I       data type number
      int16             rpt_negtim      I       =1 to set the fail condition
                                                for a (-) delta time found (for
                                                time-dependent checking)

   Modification history:

      Programmer        Date            Description of change
      ----------        ----------      ---------------------
      L. Kumar          04-Apr-1996     Original development
      W. Robinson, GSC, 2- Dec-1997     change from line-line km distance
                                        to line-line angular speed.
      W. Robinson, GSC, 31-Dec-1997     remember the highest, lowest nav
                                        rates and report them
      H. Qi, SSAI,      18-Feb-1999     remember the invalid total number 
                                        for each nflag and report them
      W. Robinson, GSC, 17-Jan-2001     turn on (-) delta time check for
                                        GAC files
      W. Robinson, SAIC  4-Oct-2001     avoid nav checks on lines flagged
                                        for unfixed navigation
      W. Robinson, SAIC 25-Jun-2002     add time-dependent checking for 
                                        (-) time changes and make slight mod
                                        to navflag checking to conform with
                                        the L0-L1 action
      W. Robinson, SAIC 20 Aug 2003     check for a + time shift of from 20 to
                                        2000 msec and report in the (-) shift 
                                        in L1NAVDISC.  This should find the 
                                        + steps for GPS resets

 ****************************************************************************/
int32 chk_nav(int32 sdfid, int32 nscans, float32 nav_thresh1,
        float32 nav_thresh2, int32 ntilts, int16 tilt_ranges[20][2],
        int16 tilt_flags[20], int16 dtynum, int16 rpt_negtim) {
    int16 good_scans[20][2];
    int32 scan, start[3], edge[3], failcode = 0;
    int32 err_thresh = 0, err_thresh1 = 0, err_thresh2 = 0;
    int32 i, fst, t_cnt = 0;
    int32 *nflag;
    int32 suncounter = 0, tiltcounter = 0, telcounter = 0;
    int32 timecounter = 0, earthcounter = 0, failcounter = 0;
    int32 fixtimecount = 0, interpcount = 0, warncount = 0;
    int32 tiltchgcount = 0;
    float32 *orb_vec, *sun_ref, *sen_mat, *scan_ell;
    float32 value, v0[3], v1[3];
    int32_t t0, t1, tstart, delt, trend_dev, *msec;
    double v0abs, v1abs, delang, d_value;
    float lo_nav_rate, hi_nav_rate, lo_nav_rate_delang,
            hi_nav_rate_delang, lo_nav_rate_delt, hi_nav_rate_delt;
    int lo_nav_rate_loc, hi_nav_rate_loc;
    char str[12];

    printf("\n\n#Checking for Navigation discontinuity ....");

    for (i = 0; i < 20; i++)
        good_scans[i][0] = good_scans[i][1] = 0;

    /* 
     *  Allocate space for the navigation buffers
     */

    if ((orb_vec = (float32 *) calloc(nscans * 3, sizeof (float32))) == NULL) {
        sprintf(err_msg, "chk_nav: cannot allocate memory for reading orb_vec");
        return FAIL;
    }

    if ((sun_ref = (float32 *) calloc(nscans * 3, sizeof (float32))) == NULL) {
        sprintf(err_msg, "chk_nav: cannot allocate memory for reading sun_ref");
        return FAIL;
    }

    if ((sen_mat = (float32 *) calloc(nscans * 3 * 3, sizeof (float32))) == NULL) {
        sprintf(err_msg, "chk_nav: cannot allocate memory for reading sen_mat");
        return FAIL;
    }

    if ((scan_ell = (float32 *) calloc(nscans * 6, sizeof (float32))) == NULL) {
        sprintf(err_msg, "chk_nav: cannot allocate memory for reading scan_ell");
        return FAIL;
    }

    if ((msec = (int32_t *) calloc(nscans, sizeof (int32_t))) == NULL) {
        sprintf(err_msg, "chk_nav: cannot allocate memory for reading msec");
        return FAIL;
    }
    if ((nflag = (int32 *) calloc(nscans * 8, sizeof (int32))) == NULL) {
        sprintf(err_msg, "chk_nav: cannot allocate memory for reading nflag");
        return FAIL;
    }

    /*
     *  Set start and end dims for reading orb_vec sds
     */
    start[0] = start[1] = start[2] = 0;
    edge[0] = nscans;
    edge[1] = 3;
    edge[2] = 0;
    if ((rdslice(sdfid, ORBVEC, start, edge, (VOIDP) orb_vec)) < 0) {
        stat_status = stat_status | 1;
        return FAIL;
    }

    /*
     *  Set start and end dims and read sun_ref sds 
     */
    start[0] = start[1] = start[2] = 0;
    edge[0] = nscans;
    edge[1] = 3;
    edge[2] = 0;
    if ((rdslice(sdfid, SUNREF, start, edge, (VOIDP) sun_ref)) < 0) {
        stat_status = stat_status | 1;
        return FAIL;
    }

    /*
     *  Set start and end dims and read sen_mat sds
     */
    start[0] = start[1] = start[2] = 0;
    edge[0] = nscans;
    edge[1] = 3;
    edge[2] = 3;
    if ((rdslice(sdfid, SENMAT, start, edge, (VOIDP) sen_mat)) < 0) {
        stat_status = stat_status | 1;
        return FAIL;
    }

    /*
     *  Set start and end dims and read scan_ell sds
     */
    start[0] = start[1] = start[2] = 0;
    edge[0] = nscans;
    edge[1] = 6;
    edge[2] = 0;
    if ((rdslice(sdfid, SCANELL, start, edge, (VOIDP) scan_ell)) < 0) {
        stat_status = stat_status | 1;
        return FAIL;
    }

    /*
     *  Set start and end dims and read msec sds
     */
    start[0] = start[1] = start[2] = 0;
    edge[0] = nscans;
    edge[1] = 0;
    edge[2] = 0;
    if ((rdslice(sdfid, "msec", start, edge, (VOIDP) msec)) < 0) {
        stat_status = stat_status | 1;
        return FAIL;
    }

    /*
     *  Set start and end dims and read nflag sds
     */
    start[0] = start[1] = start[2] = 0;
    edge[0] = nscans;
    edge[1] = 8;
    edge[2] = 0;
    if ((rdslice(sdfid, NFLAG, start, edge, (VOIDP) nflag)) < 0) {
        stat_status = stat_status | 1;
        return FAIL;
    }

    /* 
     *  loop through tilt_ranges and exclude all the lines that has tilt status
     *  equal to 3, indicating the tilt change.  
     *  SO, good_scans contains only line ranges that we can do the nav test on
     *  and t_cnt is the # of those ranges.
     */
    for (t_cnt = 0, i = 0; i < ntilts; i++) {
        if (tilt_flags[i] == 1 || tilt_flags[i] == 2) {
            good_scans[t_cnt][0] = tilt_ranges[i][0] - 1;
            good_scans[t_cnt][1] = tilt_ranges[i][1] - 1;
            t_cnt++;
        }
    }

    /*
     *  Initially, we were just going to call geovex to get the vectors
     *  out to a point for 2 lines and get the vector difference and 
     *  check that against the thresholds.  However, now there can be missing
     *  scans and the distance may change.
     *
     *  So, now, we find the angle between the 2 vectors using the 
     *  relation:  cos( angle ) = ( a * b ) / ( |a| |b| ) for vectors a, b
     *  The time difference is provided by the msec sds (make sure you
     *  account for crossing a day boundary, add 86400 if time < start).
     *  The final result is converted to milliradians / sec and for
     *  the SeaWiFS orbit of 98.9 min, this should be 1.059 mrad / sec
     *  so the thresholds should be around that value.
     */
    lo_nav_rate = 9999.;
    hi_nav_rate = -9999.;
    lo_nav_rate_loc = hi_nav_rate_loc = -1;
    lo_nav_rate_delang = lo_nav_rate_delt = 0;
    hi_nav_rate_delang = hi_nav_rate_delt = 0;

    for (i = 0; i < t_cnt; i++) /* do the next for each good line range */ {
        /*  loop through the lines in that line range */
        for (fst = 0, scan = good_scans[i][0]; scan <= good_scans[i][1]; scan++) {
            /*  
             *  deal with the valid navigation line and skip the invalid 
             *  navigation line. also, avoid any lines where the time code 
             *  was flagged but unfixed
             */
            if ((*(nflag + scan * 8) == 0) && (*(nflag + scan * 8 + 7) == 0)) {

                /*
                 *  get the vector for the line and the time 
                 */
                geovex_((orb_vec + scan * 3), (sen_mat + scan * 9),
                        (scan_ell + scan * 6), (sun_ref + scan * 3), v0);
                t0 = *(msec + scan);

                if (fst > 0) {
                    if (t0 < tstart) t0 = t0 + 86400000;
                    delt = t0 - t1;
                    /*
                     *  Do a (-) delta time check only for GAC and report only if
                     *  there is a change at least .5 - 3. sec off of the normal
                     *  GAC change of .666 sec.  So, only flag for 
                     *  -2.340 < delta t < 0.
                     *  Also, for any forward deviations of 20 to 2000 msec = 686
                     *  to 2666 msec (< 20 msec and it is negligable, 4000 msec 
                     *  and higher is a frame drop)
                     */
                    if ((dtynum == GAC) &&
                            (((delt <= 0) && (delt > -2340)) ||
                            ((delt >= 686) && (delt <= 2666)))) {
                        trend_dev = delt - 666;
                        if (rpt_negtim == 1)failcode = 1;
                        printf("\n\n\nTIME_ERROR  scan # [%d-%d]   ", scan, scan + 1);
                        printf("Delta time = %d,  Trend deviation = %d\n\n\n", delt, trend_dev);
                        printf("NOTE that if this time is not in a mandatory product\n");
                        printf("failure time range, it may need to be put in one.\n");
                        printf("For post-repro 3 data this denotes a 30 second time\n");
                        printf("shift that wasn't detected or a GPS reset\n");
                        printf("Check around this time period for navigation offsets\n");
                        printf("using GAC, LAC and HRPT stations.\n");
                        printf("Also, inform QC and navigation manager of this problem\n");
                        printf("(This may require an update to l1stat_chk)\n\n\n");
                    }
                    v0abs = sqrt(v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2]);
                    v1abs = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);

                    /*  this may be poor for the small angles we should encounter
                     delang = acos( (v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2] ) /
                                   ( v0abs * v1abs ) );
                        SO, try method where we normalize both radius vectors, find
                        the magnitude of the vector difference d_value, and then find the
                        angle with: sin( 1/2 angle ) = 1/2 d_value
                     */
                    d_value = sqrt(pow((v0[0] / v0abs - v1[0] / v1abs), 2) +
                            pow((v0[1] / v0abs - v1[1] / v1abs), 2) +
                            pow((v0[2] / v0abs - v1[2] / v1abs), 2));

                    delang = 2. * asin(d_value / 2.);

                    value = (delang * 1.e6) / (float) delt;

                    /*
                     *  remember the lowest, highest rates
                     */
                    if (value < lo_nav_rate) {
                        lo_nav_rate = value;
                        lo_nav_rate_loc = scan;
                        lo_nav_rate_delang = delang;
                        lo_nav_rate_delt = delt;
                    }
                    if (value > hi_nav_rate) {
                        hi_nav_rate = value;
                        hi_nav_rate_loc = scan;
                        hi_nav_rate_delang = delang;
                        hi_nav_rate_delt = delt;
                    }

                    if (value < nav_thresh1) {
                        failcode = 1;
                        err_thresh++;
                        err_thresh1++;
                        printf("\nNAVERR   scan # [%d-%d]   ", scan, scan + 1);
                        printf("value = %f, threshold = %f", value, nav_thresh1);
                        printf(" Kdelang: %f, delt: %d", 1000. * delang, delt);
                    } else if (value > nav_thresh2) {
                        failcode = 1;
                        err_thresh++;
                        err_thresh2++;
                        printf("\nNAVERR   scan # [%d-%d]   ", scan, scan + 1);
                        printf("value = %f, threshold = %f", value, nav_thresh2);
                        printf(" Kdelang: %f, delt: %d", 1000. * delang, delt);
                    }
                } else {
                    tstart = t0; /* first time thru, define start time */
                }
                fst++;
                v1[0] = v0[0];
                v1[1] = v0[1];
                v1[2] = v0[2];
                t1 = t0;
            }
            if (*(nflag + scan * 8) == 1) failcounter++;
            if (*(nflag + scan * 8 + 1) == 1) interpcount++;
            if (*(nflag + scan * 8 + 2) == 1) suncounter++;
            if (*(nflag + scan * 8 + 3) == 1) earthcounter++;
            if (*(nflag + scan * 8 + 4) == 1) telcounter++;
            if (*(nflag + scan * 8 + 5) == 1) timecounter++;
            if (*(nflag + scan * 8 + 5) == 2) fixtimecount++;
            if (*(nflag + scan * 8 + 6) == 1) tiltcounter++;
            if (*(nflag + scan * 8 + 6) == 2) tiltchgcount++;
            if (*(nflag + scan * 8 + 7) == 1) warncount++;
        }
    }

    printf("\n\n# --- Cumulative navigation flag (nflag) settings ---\n");
    printf("Description                              # of occurences\n");
    printf("--------------------------------------   ---------------\n");
    printf("navigation failure      (nflag(1) = 1)   %d\n", failcounter);
    printf("nav interp problem      (nflag(2) = 1)   %d\n", interpcount);
    printf("sun sensor out          (nflag(3) = 1)   %d\n", suncounter);
    printf("earth sensor out        (nflag(4) = 1)   %d\n", earthcounter);
    printf("attitude uncertainty    (nflag(5) = 1)   %d\n", telcounter);
    printf("time code problem       (nflag(6) = 1)   %d\n", timecounter);
    printf("time code problem fixed (nflag(6) = 2)   %d\n", fixtimecount);
    printf("tilt info bad           (nflag(7) = 1)   %d\n", tiltcounter);
    printf("tilt changing           (nflag(7) = 2)   %d\n", tiltchgcount);
    printf("navigation warning      (nflag(8) = 1)   %d\n\n", warncount);

    printf("\n\n\n#Mnemonic  Code  tot # outside\tthresh_1   thresh_2  # below 1  # below 2");
    printf("\n                 thresholds ");
    printf("\n#--------------------------------------------------------------------------\n");

    printf("\nL1NAVDISC   %d\t    %d\t       %f\t  %f\t%d\t %d", failcode,
            err_thresh, nav_thresh1, nav_thresh2, err_thresh1, err_thresh2);

    /*
     *  report the low and high rates
     */
    printf("\n\nRecord of lowest, highest rates found ( ideal rate: 1.059)\n");
    printf("        rate (mrad / sec)  line location  delang (mrad) delt (msec)\n");
    printf("LOWRATE:%-12.4f        %-10d     %-12.4f  %-12.4f\n", lo_nav_rate, lo_nav_rate_loc, lo_nav_rate_delang * 1000., lo_nav_rate_delt);
    printf("HIRATE: %-12.4f        %-10d     %-12.4f  %-12.4f\n", hi_nav_rate, hi_nav_rate_loc, hi_nav_rate_delang * 1000., hi_nav_rate_delt);


    if (failcode) {
        stat_status = stat_status | 2;
        sprintf(str, "L1NAVDISC ");
        if (strlen(bad_stat_str) <= 300)
            strcat(bad_stat_str, str);
    }

    free(orb_vec);
    free(sun_ref);
    free(scan_ell);
    free(sen_mat);
    free(msec);
    free(nflag);
    return SUCCEED;
}

/****************************************************************************
   chk_tilt

   purpose: Calculates time required by tilt changes and verifies to see
                if it is less than the given threshold value. 

   Returns type: int32

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32		sdfid		I	SD interface id
      int16 		dtynum    	I  	data type number
      float32 		l1tilt_thresh	I 	input threshold for tiltchange
      int32 *		ntilts		O	number tilts
      int16 		tilt_ranges[20][2] O    scan line ranges of tilt states
      int16		tilt_flags[20]  O	flag corresponding to each 
                                                tilt states

   Modification history:

      Programmer        Date            Description of change
      ----------        ----------      ---------------------
      L. Kumar          04-Apr-1996     Original development

 ****************************************************************************/
int32 chk_tilt(int32 sdfid, int16 dtynum, float32 l1tilt_thresh, int32 *ntilts,
        int16 tilt_ranges[20][2], int16 *tilt_flags) {
    int32 i, failcode = 0, start[3], edge[3];
    float32 tilt_value = 0;
    char str[12];

    printf("\n\n#Checking Tilt behavior....");

    /* 
     *  Set start and end dims of ntilts sds and read ntilts fr given data file
     */
    start[0] = start[1] = start[2] = 0;
    edge[0] = 1;
    edge[1] = edge[2] = 0;
    if ((rdslice(sdfid, NTILTS, start, edge, (VOIDP) ntilts)) < 0) {
        stat_status = stat_status | 1;
        return FAIL;
    }

    /* 
     *  Set dimensions and read 'tilt_flags' from the given data file
     */
    edge[0] = 20;
    edge[1] = edge[2] = 0;
    if ((rdslice(sdfid, TILT_FLAGS, start, edge, (VOIDP) tilt_flags)) < 0) {
        stat_status = stat_status | 1;
        return FAIL;
    }

    /* 
     *  Set dimensions of 'tilt_ranges' sds and read it from the data file
     */
    edge[0] = 20;
    edge[1] = 2;
    edge[2] = 0;
    if ((rdslice(sdfid, TILT_RANGES, start, edge, (VOIDP) tilt_ranges)) < 0) {
        stat_status = stat_status | 1;
        return FAIL;
    }

    printf("\n\n#Mnemonic  Code\t  Value\t   --\tError_threshold");
    printf("\n#-----------------------------------------------\n");

    for (i = 0; i < *ntilts; i++) {
        if (tilt_flags[i] == 3) {
            failcode = 0;
            tilt_value = (tilt_ranges[i][1] - tilt_ranges[i][0] + 1);
            if (dtynum == GAC)
                tilt_value = (tilt_value * 2.0 / 3.0);
            else
                tilt_value = (tilt_value / 6.0);
            if (tilt_value > l1tilt_thresh)
                failcode = 1;
        }
    }
    printf("\nL1TILT     %d\t%f\t%f", failcode, tilt_value, l1tilt_thresh);
    if (failcode) {
        stat_status = stat_status | 2;
        sprintf(str, "L1TILT ");
        if (strlen(bad_stat_str) <= 300)
            strcat(bad_stat_str, str);
    }

    return SUCCEED;
}

void stat_exit(int status)
/*******************************************************************

   stat_exit

   purpose: provide a common exit from the stat_check program and
            error summary info.

   Returns type: none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               status          I       statistical check status

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       30-jun-1995     Original development

 *******************************************************************/ {
    if (status == 0)
        printf("\n\nSuccessful check, no statistical errors\n");
    else {
        if (status & 1)
            printf("\n\nFailure of statistical check due to program error\n");
        if (status & 2) {
            printf("\n\nFailure of statistical check due to data error\n\n");
            printf("Summary of mnemonics with bad status:\n");
            printf("%s\n", bad_stat_str);
        }
    }
    exit(status);
}
