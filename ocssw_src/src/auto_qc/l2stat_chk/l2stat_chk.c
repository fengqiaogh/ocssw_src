/*
 *  W. Robinson, GSC, 13 Dec 1999 update for version 4.0 - the new format
 *                                from MSL12
 */
#include <time.h>
#include "l2stat.h"
#include "l2lists.h"
#include "l2stat_proto.h"

#define MAX_NAME 255
#define SET 1

char err_msg[1024];
static int32 stat_status = 0; /* status of statistical checking: 0 all good,
        		   	   1 program problem, 2 statistical problem, 
		 	   	   3 both problems  */

int main(int argc, char *argv[])
/*******************************************************************

   l2stat_chk

   purpose: open the SeaWiFS level 2 dataset and checks the data values 
                against given thresholds

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
      L. Kumar          06-Jun-1996     Removed '#ifdef PROTOTYPE' defns. and
                                        '\t' (tabs) from the printf statements
      W. Robinson, GSC  5-Nov-1997      add extra argument for extra flags to 
                                        exclude
      W. Robinson, GSC 10-Dec-1999      update for version 4.0 of L2 format,
                                        (different parameters and 32 bit flags)
      W. Robinson, SAIC, 27jun2002      update for chlor_a stored as a float
                                        for REPRO 4
 *******************************************************************/ {
    char dtype[255], extra_flags[200];
    int32 fid, sdfid;
    int32 nsamp, nscans;
    cntl_str grschk[NPARAMS], statchk[NPARAMS];
    flag_str flgchk[NFLAGS];

    clock_t val = 0, val2 = 0;

    val = (float) clock() / CLOCKS_PER_SEC;

    /*
     *  Check input arguments
     */

    /* WDR chg:if (argc != 3) { */
    if (argc < 3 || argc > 4) {
        printf("\n******* Usage: l2stat_chk <hdf file> <control file> [<flg>]\n");
        printf("   <hdf file> is level 2 data file\n");
        printf("   <control file> is the file specifying thresholds\n");
        printf("   <flg> is any extra flag names to exclude (optional)\n");
        printf("\nOn exit, $status will be set to: \n  0 - no problems,"
                "\n  1 - program error, ");
        printf("\n  2 - data problem, \n  3 - both program and data problem\n");
        stat_status = stat_status | 1;
        stat_exit(stat_status);
    }

    printf("\n\n# Statistical Check of Dataset '%s'\n\n", argv[1]);

    /*
     *   hopen will open the file and read data descriptor blocks
     *   to memory
     */

    if ((fid = Hopen(argv[1], DFACC_RDONLY, 0)) < 0) {
        printf("****** l2stat_chk: Failed on the Hopen of \n  '%s'\n", argv[1]);
        stat_status = 1;
        stat_exit(stat_status);
    }

    /*  
     *	SDstart opens the hdf interface and initiates SD interface
     */

    if ((sdfid = SDstart(argv[1], DFACC_RDONLY)) < 0) {
        printf("******* l2stat_chk: Failure at SDstart of \n'%s'\n", argv[1]);
        Hclose(fid);
        stat_status = 1;
        stat_exit(stat_status);
    }

    /*
     *  Read control file and save the information.  If the return value
     *  is -1, it means, some prolblem with opening control data file or
     *  with some control format, close HDF file, and exit with status set
     *  to 1 indicating program error.
     */

    if ((read_cntldata(argv[2], grschk, statchk, flgchk)) < 0) {
        Hclose(fid);
        stat_exit(stat_status);
    }

    /*
     *  Verify if given input file is a level 2 file
     */

    if ((l2file(sdfid, &nsamp, &nscans, dtype)) < 0) {
        Hclose(fid);
        printf("******* %s", err_msg);
        stat_exit(stat_status);
    }

    /*
     *  before the gross check, set up the extra flag list
     *  depending on a 4th argument
     */
    if (argc == 4) {
        strcpy(extra_flags, argv[3]);
    } else
        strcpy(extra_flags, "none");

    /*
     *  Do a gross value check and general statistics check
     */

    if ((chk_grsstat(sdfid, nscans, nsamp, grschk, statchk, extra_flags)) < 0)
        printf("\n******* %s", err_msg);

    /* 
     *  Verify Flag percentage
     */

    if ((chk_flg(sdfid, flgchk)) < 0)
        printf("\n******** %s", err_msg);

    /*
     *  Free all the allocated resources
     */

    Hclose(fid);

    val2 = (float) clock() / CLOCKS_PER_SEC;
    printf("\n\n# Time took for %s = %d secs\n", argv[0], (int) (val2 - val));

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
      cntl_str *	grschk		 O	structure to hold control
                                                data related to gross checking
      cntl_str *	statchk		 O	structure to hold control
                                                data related to general stat
                                                check
      flag_str *	flgchk		 O	structure to hold control
                                                data related to flag percentage
       
   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      L. Kumar          04-Apr-1996    Original development

 ****************************************************************************/
int32 read_cntldata(char *cntl_file, cntl_str *grschk, cntl_str *statchk,
        flag_str *flgchk) {
    FILE *fid;
    char line[501], str[25];
    int32 i, param, flag, nflds = 0;
    float32 err_thresh, low_thresh, high_thresh;

    printf("\n# Reading control file '%s'\n\n", cntl_file);

    /*
     *  Open the file
     */
    if ((fid = fopen(cntl_file, "r")) == NULL) {
        printf("********read_cntldata: Unable to open control file:\n  '%s'\n",
                cntl_file);
        stat_status = stat_status | 1;
        return FAIL;
    }

    /*
     *  Initialize all input structrues
     */
    for (i = 0; i < NPARAMS; i++) {
        grschk[i].param = statchk[i].param = 0;
        grschk[i].err_thresh = statchk[i].err_thresh = 0;
        grschk[i].low_thresh = statchk[i].low_thresh = 0;
        grschk[i].high_thresh = statchk[i].high_thresh = 0;
    }

    for (i = 0; i < NFLAGS; i++)
        flgchk[i].flag = flgchk[i].err_low_thresh = flgchk[i].err_high_thresh = 0;

    /*
     *  Read control file one line at a time
     */
    while (fgets(line, 500, fid) != NULL) {
        /*  Ignore beginning blanks and tabs */
        for (i = 0; i < 500 && (line[i] == ' ' || line[i] == '\t'); i++);

        /* If first character is '#' sign, then treat the line as comment */
        if (i < 500 && line[i] == '#') {
            printf("%s", line);
            continue;
        }

        /* If not comment, check if it is Gross CHK info */
        printf("#%s", line);
        if (strncmp(&line[i], "L2GCHK", 6) == 0) {
            if ((nflds = sscanf(line, "%s %f %f %f", str, &err_thresh,
                    &low_thresh, &high_thresh)) != 4) {
                printf("\n*********read_cntldata: expecting 4 fields");
                printf(" but read %d field(s)", nflds);
                stat_status = stat_status | 1;
                return FAIL;
            }
            /*
                     prm[0] = str[6];
             */
            param = atoi(&str[6]);
            if (param < 1 || param > NPARAMS) {
                printf("********read_cntldata: Error in parameter number.");
                printf("  Parameter # read = %d", param);
                stat_status = stat_status | 1;
                return FAIL;
            }
            grschk[param - 1].param = 1;
            grschk[param - 1].err_thresh = err_thresh;
            grschk[param - 1].low_thresh = low_thresh;
            grschk[param - 1].high_thresh = high_thresh;
        } else if (strncmp(&line[i], "L2STAT", 6) == 0) {
            if ((nflds = sscanf(line, "%s", str)) != 1) {
                printf("*********read_cntldata: expecting 1 field");
                printf(" but read %d field(s)", nflds);
                stat_status = stat_status | 1;
                return FAIL;
            }
            param = atoi(&str[6]);
            if (param < 1 || param > NPARAMS) {
                printf("********read_cntldata: Error in param number.");
                printf("  Parameter # read = %d", param);
                stat_status = stat_status | 1;
                return FAIL;
            }
            statchk[param - 1].param = 1;
        } else if (strncmp(&line[i], "L2FLGCK", 7) == 0) {
            if ((nflds = sscanf(line, "%s %f %f", str, &low_thresh, &high_thresh))
                    != 3) {
                printf("*********read_cntldata: expecting 3 fields");
                printf(" but read %d field(s)", nflds);
                stat_status = stat_status | 1;
                return FAIL;
            }
            flag = atoi(&str[7]);
            if (flag < 1 || flag > NFLAGS) {
                printf("*********read_cntldata: Error in flag number.");
                printf("  Flag # read = %d", flag);
                stat_status = stat_status | 1;
                return FAIL;
            }
            flgchk[flag - 1].flag = 1;
            flgchk[flag - 1].err_low_thresh = low_thresh;
            flgchk[flag - 1].err_high_thresh = high_thresh;
        }
    }
    return SUCCEED;
}

/****************************************************************************
   l2file

   purpose: Verifies if the given input file is level 2 data file

   Returns type: integer

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             sdfid           I       SD interface ID
      int32 *		nsamp		O	number of samples/pixels
      int32 *		nscans		O	number of scan lines
      char  *		dtype		O  	data type (GAC, LAC, ...)

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      L. Kumar          04-Apr-1996     Original development

 ****************************************************************************/
int32 l2file(int32 sdfid, int32 *nsamp, int32 *nscans, char *dtype) {
    char title[1024];

    /*
     *  Read title to verify if the given input file is level 1A data file
     */

    if ((rdattr(sdfid, TITLE, &title)) < 0)
        return FAIL;

    if (strcmp(title, "SeaWiFS Level-2 Data") != 0) {
        sprintf(err_msg, "l2file: Data file is not level 2 file");
        stat_status = stat_status | 1;
        return FAIL;
    }

    if ((rdattr(sdfid, NSAMP, (VOIDP *) nsamp)) < 0)
        return FAIL;

    printf("\n# Number of samples = %d", *nsamp);

    if ((rdattr(sdfid, NSCANS, (VOIDP *) nscans)) < 0)
        return FAIL;

    printf("\n# Number of scanlines = %d", *nscans);

    if ((rdattr(sdfid, DTYPE, dtype)) < 0)
        return FAIL;

    printf("\n# Data type = %s\n", dtype);

    return SUCCEED;
}

/****************************************************************************
   chk_grsstat 

   purpose: Verifies Gross value and statistical checks against the given 
                thresholds for requested parameters

   Returns type: integer 

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             sdfid           I       SD interface ID
      int32		nscans		I	Number of scanlines
      int32		nsamp		I	Number of pixels
      cntl_str	*	grschk 		I       structure containing Gross
                                                value control thresholds
      int16 *		buf		I	buffer for reading data
      char *            extra_flags     I       user-requested extra flags
                                                to use in excluding data
      
   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      L. Kumar          04-Apr-1996     Original development
      W. Robinson, GSC  5-Nov-1997      add the extra_flags argument
      W. Robinson, GSC  13-Dec-1999     use 32 bit flag array and tau as a 
                                        16 bit value
      W. Robinson, SAIC, 27jun2002      update for chlor_a stored as a float
                                        for REPRO 4

 ****************************************************************************/
int32 chk_grsstat(int32 sdfid, int32 nscans, int32 nsamp, cntl_str *grschk,
        cntl_str *statchk, char *extra_flags) {
    uint8 *i8buf;
    uint16 mask;
    int32 *l2_flags;
    int16 *buf, flag[NPARAMS];
    int32 i, p, npix = (nsamp * nscans);
    int32 edge[3], start[3];
    int32 non_masked_pixels, low_cnt = 0, high_cnt = 0;
    int32 minloc[NPARAMS], maxloc[NPARAMS];
    int32 failcode = 0;
    float32 min[NPARAMS], max[NPARAMS], *flt32buf;
    float64 geoval, value1[NPARAMS], value2[NPARAMS];
    float64 xbar[NPARAMS], sd[NPARAMS], sum, sumsq;

    printf("\n#Checking Gross/Statistical values ....\n");

    for (i = 0; i < NPARAMS; i++) {
        value1[i] = value2[i] = flag[i] = 0;
        xbar[i] = sd[i] = 0;
        min[i] = max[i] = minloc[i] = maxloc[i] = 0;
    }

    /*
     *  Allocate buffer for reading a line of data with 4 byte values 
     *  for 4byte chlor_a values in repro4
     */

    if ((buf = (int16 *) calloc(npix * 2, sizeof (int16))) == NULL) {
        sprintf(err_msg, "chk_grsstat: cannot allocate memory for reading data");
        stat_status = stat_status | 1;
        return FAIL;
    }

    /*
     *  Allocate buffer for reading l2_flags and set start and end dims for 
     *  reading l2_flags
     */

    if ((l2_flags = (int32 *) calloc(npix, sizeof (int32))) == NULL) {
        sprintf(err_msg, "chk_grsstat: cannot allocate memory for reading data");
        stat_status = stat_status | 1;
        free(buf);
        return FAIL;
    }

    start[0] = start[1] = start[2] = 0;
    edge[0] = nscans;
    edge[1] = nsamp;
    edge[2] = 0;
    if ((rdslice(sdfid, L2FLAGS, start, edge, (VOIDP) l2_flags)) < 0) {
        stat_status = stat_status | 1;
        free(buf);
        free(l2_flags);
        return FAIL;
    }

    /* 
     * Call set_mask for setting mask, and find total number of non masked values
     */

    if (set_mask(sdfid, extra_flags, &mask) < 0) {
        stat_status = stat_status | 1;
        free(buf);
        free(l2_flags);
        return FAIL;
    }

    start[0] = start[1] = start[2] = 0;
    edge[0] = nscans;
    edge[1] = nsamp;
    edge[2] = 0;

    for (p = 0; p < NPARAMS; p++)
        if (grschk[p].param || statchk[p].param) {
            printf("\n#Parm-%d  %s", p + 1, param_name[p]);
            printf("   Slope = %f, Intercept = %f", slope[p], intercept[p]);
            low_cnt = high_cnt = 0;
            min[p] = 9999;
            max[p] = -9999;
            minloc[p] = maxloc[p] = 0;
            non_masked_pixels = sum = sumsq = 0;
            if ((rdslice(sdfid, param_name[p], start, edge, (VOIDP) buf)) < 0) {
                stat_status = stat_status | 1;
                free(buf);
                free(l2_flags);
                return FAIL;
            }
            if ((strcmp(param_name[p], "eps_78")) == 0) {
                i8buf = (uint8 *) buf;
                for (i = 0; i < npix; i++)
                    if ((l2_flags[i] & mask) == 0) {
                        geoval = (float64) i8buf[i] * slope[p] + intercept[p];
                        if (geoval < grschk[p].low_thresh)
                            low_cnt++;
                        if (geoval > grschk[p].high_thresh)
                            high_cnt++;
                        sum += geoval;
                        sumsq += geoval * geoval;
                        non_masked_pixels++;
                        if (geoval < min[p]) {
                            min[p] = geoval;
                            minloc[p] = i;
                        }
                        if (geoval > max[p]) {
                            max[p] = geoval;
                            maxloc[p] = i;
                        }
                    }
            } else if ((strcmp(param_name[p], "chlor_a")) == 0) {
                flt32buf = (float32 *) buf;
                for (i = 0; i < npix; i++)
                    if ((l2_flags[i] & mask) == 0) {
                        geoval = (float64) flt32buf[i] * slope[p] + intercept[p];
                        if (geoval < grschk[p].low_thresh)
                            low_cnt++;
                        if (geoval > grschk[p].high_thresh)
                            high_cnt++;
                        sum += geoval;
                        sumsq += geoval * geoval;
                        non_masked_pixels++;
                        if (geoval < min[p]) {
                            min[p] = geoval;
                            minloc[p] = i;
                        }
                        if (geoval > max[p]) {
                            max[p] = geoval;
                            maxloc[p] = i;
                        }
                    }
            } else {
                for (i = 0; i < npix; i++)
                    if ((l2_flags[i] & mask) == 0) {
                        geoval = (float64) buf[i] * slope[p] + intercept[p];
                        if (geoval < grschk[p].low_thresh)
                            low_cnt++;
                        if (geoval > grschk[p].high_thresh)
                            high_cnt++;
                        sum += geoval;
                        sumsq += geoval * geoval;
                        non_masked_pixels++;
                        if (geoval < min[p]) {
                            min[p] = geoval;
                            minloc[p] = i;
                        }
                        if (geoval > max[p]) {
                            max[p] = geoval;
                            maxloc[p] = i;
                        }
                    }
            }
            /* 
             * calculate percentage of non-masked values < given low_threshold and
             * values > given high threshold
             */
            if (low_cnt > 0)
                value1[p] = (low_cnt / (float64) (non_masked_pixels) * 100.0);
            if (high_cnt > 0)
                value2[p] = (high_cnt / (float) (non_masked_pixels) * 100.0);
            if (value2[p] > grschk[p].err_thresh)
                flag[p] = -1;
            if (value1[p] > grschk[p].err_thresh)
                flag[p] = 1;

            /* 
             *  Gather statistical information and output the info later
             */
            if (non_masked_pixels) {
                xbar[p] = sum = sum / non_masked_pixels;
                sumsq /= non_masked_pixels;
                if ((sumsq - sum * sum) <= 0)
                    sd[p] = 0;
                else
                    sd[p] = sqrt(sumsq - sum * sum);
            } else
                xbar[p] = sd[p] = 0;
        }

    /*  
     *  Print headers for gross check output messages
     */

    printf("\n\n#Mnemonic  Flag  %%<Low   %%>High   Err_thr   low_thr   high_thr");
    printf("\n#-------------------------------------------------------------\n");

    for (p = 0; p < NPARAMS; p++)
        if (grschk[p].param) {
            printf("\nL2GCHK%d   %d  %f  %f  %f  %f  %f", p + 1, flag[p], value1[p],
                    value2[p], grschk[p].err_thresh, grschk[p].low_thresh,
                    grschk[p].high_thresh);
            if (flag[p] != 0)
                failcode = 1;
        }

    /*  
     *  Print headers for statistical check output messages
     */

    printf("\n\n");
    printf("#Mnemonic   non      xbar        sd      min      min_loc     max    max_loc");
    printf("\n#       masked pixels                             (x,y)               (x,y)");
    printf("\n#--------------------------------------------------------------------------\n");

    for (p = 0; p < NPARAMS; p++)
        if (statchk[p].param)
            printf("\nL2STAT%d   %d   %f   %f  %f  %d,%d  %f  %d,%d",
                p + 1, non_masked_pixels, xbar[p], sd[p], min[p], minloc[p] % nsamp,
                minloc[p] / nsamp, max[p], maxloc[p] % nsamp, maxloc[p] / nsamp);

    if (failcode)
        stat_status = stat_status | 2;

    free(l2_flags);
    free(buf);
    return SUCCEED;
}

/****************************************************************************
   set_mask

   purpose: Sets mask field based on masknames stored in the datafile

   Returns type: 2 byte integer <0 is fail

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             sdfid           I       SD interface ID
      char *            extra_flags     I       string of extra flags to mask
      uint16 *          mask            O       mask returned

   Modification history:

      Programmer        Date            Description of change
      ----------        ----------      ---------------------
      L. Kumar          04-Apr-1996     Original development
      W. Robinson, GSC  5-Nov-1997      add extra_flags argument

 ****************************************************************************/
int16 set_mask(int32 sdfid, char *extra_flags, uint16 *mask) {
    int32 nt, count, i;
    char *masknames;

    /*
     *  Reads Mask Names global attribute 
     */
    *mask = 0;

    if ((getattrsz(sdfid, MASKNAMES, &nt, &count)) < 0)
        return FAIL;

    if ((masknames = (char *) malloc(sizeof (char)*count + 1)) == NULL) {
        sprintf(err_msg, "set_mask: malloc error -- while alloc space for %s",
                MASKNAMES);
        return FAIL;
    }

    if ((rdattr(sdfid, MASKNAMES, masknames)) < 0)
        return FAIL;

    if (masknames != NULL)
        for (i = 0; i < NFLAGS; i++) {
            if ((strstr(masknames, flag_names[i]) != NULL) ||
                    (strstr(extra_flags, flag_names[i]) != NULL))
                *mask = *mask + (int16) pow(2, i);
        }

    printf("\n#Masknames: %s\n", masknames);
    printf("#(note user-designated flags are: '%s')\n", extra_flags);

    free(masknames);
    return 0;
}

/****************************************************************************
   chk_flag

   purpose: Sets mask field based on masknames stored in the datafile

   Returns type: 2 byte integer 

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             sdfid           I       SD interface ID

   Modification history:

      Programmer        Date            Description of change
      ----------        ----------      ---------------------
      L. Kumar          04-Apr-1996     Original development

 ****************************************************************************/
int32 chk_flg(int32 sdfid, flag_str *flgchk) {
    int32 i, flag = 0;
    float32 perct_flags[NFLAGS];

    printf("\n\n#Checking Flag Percentages ....");

    if ((rdattr(sdfid, PERCENTFLAGS, perct_flags)) < 0)
        return FAIL;

    printf("\n\n#Mnemonic    flag       flag %%       err_thr_low      err_thr_high");
    printf("\n#------------------------------------------------------------------\n");

    for (i = 0; i < NFLAGS; i++)
        if (flgchk[i].flag) {
            flag = 0;
            if (perct_flags[i] < flgchk[i].err_low_thresh)
                flag = -1;
            if (perct_flags[i] > flgchk[i].err_high_thresh)
                flag = 1;
            if (flag != 0)
                stat_status = stat_status | 2;
            printf("\nL2FLGCK%d     %d        %f        %f        %f", i + 1, flag,
                    perct_flags[i], flgchk[i].err_low_thresh,
                    flgchk[i].err_high_thresh);
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
        if (status & 2)
            printf("\n\nFailure of statistical check due to data error\n");
    }
    exit(status);
}

/*-----------------------------------------------------------------------------
    Function:  rdattr

    Returns:   int32 (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        The function rdattr reads the requested global attribute

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     sdfid        I    ID req to access HDF SDS interface
        char  *   attr_name    I    attribute name
        void  *   buf         I/O   pointer to data buffer

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     11/07/94  Original development

----------------------------------------------------------------------------*/
int32 rdattr(int32 sdfid, char *attr_name, void *buf) {
    int32 attrnum;

    if ((attrnum = SDfindattr(sdfid, attr_name)) < 0) {
        sprintf(err_msg, "rdattr: Failure in SDfindattr while trying to read %s",
                attr_name);
        stat_status = stat_status | 1;
        return FAIL;
    }

    if ((SDreadattr(sdfid, attrnum, buf)) < 0) {
        sprintf(err_msg, "rdattr: Failure in SDreadattr while trying to read %s",
                attr_name);
        stat_status = stat_status | 1;
        return FAIL;
    }

    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  rdslice

    Returns:   int32 (sdsid)
        The return code is a negative value if any error occurs, otherwise,
        returns sdsid.

    Description:
        The function rdslice reads requested slice of data from the
        given named dataset

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     sdfid         I    ID req to access HDF SDS interface
        char      *name         I   SDS name
        int32     *start        I   start data dimension
        int32     *edge         I   no. of values to be read
        void      *buf          O   SDS data buffer

    NOTE:

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     11/02/94  Original development
----------------------------------------------------------------------------*/
int32 rdslice(int32 sdfid, char *name, int32 *start, int32 *edge, void *buf) {
    int32 index, sdsid, rank, num_type, nattrs;
    char sdsname[MAX_NAME];


    if ((index = SDnametoindex(sdfid, name)) < 0) {
        sprintf(err_msg, "rdslice: SDnametoindex failed for sds \"%s\" ", name);
        return FAIL;
    }
    if ((sdsid = SDselect(sdfid, index)) < 0) {
        sprintf(err_msg, "rdslice: SDselect failed for sds \"%s\" ", name);
        return FAIL;
    }

    if (edge[0] == 0 && edge[1] == 0 && edge[2] == 0)
        if ((SDgetinfo(sdsid, sdsname, &rank, edge, &num_type, &nattrs)) < 0) {
            sprintf(err_msg, "rdslice: SDgetinfo failed for sds \"%s\" ", name);
            return FAIL;
        }

    if ((SDreaddata(sdsid, start, NULL, edge, buf)) < 0) {
        sprintf(err_msg,
                "rdslice: SDreaddata error while reading \"%s\" ", name);
        return FAIL;
    }


    SDendaccess(sdsid);
    return SUCCEED;
}

/*-----------------------------------------------------------------------------
    Function:  getattrsz

    Returns:   int32 (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        The function getattrsz passes the requested global attribute's
        number type (data type) and the number of values

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32     sdfid        I    ID req to access HDF SDS interface
        char  *   attr_name    I    attribute name
        int32 *   nt           O    HDF data type
        int32 *   count        O    number of values in the specified attribute

    Modification history:
    Programmer     Organization   Date      Description of change
    -------------- ------------   --------  ---------------------
    Lakshmi Kumar  Hughes STX     11/07/94  Original development
----------------------------------------------------------------------------*/
int32 getattrsz(int32 id, char *attr_name, int32 *nt, int32 *count) {
    int32 attrnum;
    char name[MAX_NAME];

    attrnum = SDfindattr(id, attr_name);
    if ((SDattrinfo(id, attrnum, name, nt, count)) < 0) {
        sprintf(err_msg, "getattrsz: SDattrinfo failed for attribute - %s\n",
                attr_name);
        return FAIL;
    }
    return SUCCEED;
}

