#include <time.h>
#include <math.h>
#include <hdf4utils.h>
#include "l3stat.h"
#include "l3lists.h"
#include "l3stat_proto.h"
#include "hist_proto.h"

#include <hdf.h>
#include <mfhdf.h>

static int32_t stat_status = 0; /* status of statistical checking: 0 all good,
                                   1 program problem, 2 statistical problem, 
                                   3 both problems  */

#define NEG_FLAG -9.1E6

int main(int argc, char *argv[])
/*******************************************************************

   l3stat_chk

   purpose: open the SeaWiFS level 3 dataset and checks the data values 
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
      L. Kumar          07/19/1996     Original development
      W. Robinson       24 Nov 1997    Make the climatology optional

 *******************************************************************/ {
    char ptype[MAXVAL];
    int32_t fid, sdfid, c_fid, c_sdfid;
    int32_t nbins, c_nbins, wtchk;
    cntl_str databinchk[1], grschk[NPARAMS], statchk[NPARAMS];
    clim_str climchk[NFLAGS];

    clock_t val = 0, val3 = 0;

    val = (float) clock() / CLOCKS_PER_SEC;
    c_fid = 0;
    c_sdfid = 0;
    c_nbins = 0;
    nbins = 0;

    /*
     *  Check input arguments
     */

    if (argc != 4 && argc != 3) {
        printf("\n***** Usage: l3stat_chk <l3 bin file> <control file> ");
        printf("[<climatology bin file>]\n\n");
        printf("         <l3 bin file> is level 3 bin data file\n");
        printf("         <control file> is the file specifying thresholds\n");
        printf("         [<climatology bin file>] full name of optional clim bin file ");
        printf("\n\nOn exit, $status will be set to: \n\n  0 - no problems,"
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
        printf("****** l3stat_chk: Failed on the Hopen of \n  '%s'\n", argv[1]);
        stat_status = 1;
        stat_exit(stat_status);
    }

    if (argc == 4) {
        if ((c_fid = Hopen(argv[3], DFACC_RDONLY, 0)) < 0) {
            printf("****** l3stat_chk: Failed on the Hopen of \n  '%s'\n", argv[3]);
            stat_status = 1;
            stat_exit(stat_status);
        }
    }

    /* 
     *  Vstart initializes vgroup interface
     */

    Vstart(fid);
    if (c_fid != 0) Vstart(c_fid);

    /*  
     *   SDstart opens the hdf interface and initiates SD interface
     */

    if ((sdfid = SDstart(argv[1], DFACC_RDONLY)) < 0) {
        printf("******* l3stat_chk: Failure at SDstart of \n'%s'\n", argv[1]);
        Hclose(fid);
        stat_status = 1;
        stat_exit(stat_status);
    }

    if (c_fid != 0) {
        if ((c_sdfid = SDstart(argv[3], DFACC_RDONLY)) < 0) {
            printf("******* l3stat_chk: Failure at SDstart of \n'%s'\n", argv[3]);
            Hclose(fid);
            Hclose(c_fid);
            stat_status = 1;
            stat_exit(stat_status);
        }
    }

    /*
     *  Read control file and save the information.  If the return value
     *  is -1, it means, some prolblem with opening control data file or
     *  with some control format, close HDF file, and exit with status set
     *  to 1 indicating program error.
     */

    if ((read_cntldata(argv[2], databinchk, grschk, statchk, &wtchk, climchk))
            < 0) {
        Hclose(fid);
        stat_exit(stat_status);
    }

    /*
     *  Verify if given input file is a level 3 file
     */

    if (l3file(sdfid, c_sdfid, &nbins, &c_nbins, ptype) < 0) {
        Hclose(fid);
        stat_exit(stat_status);
    }

    /*
     *  Check percent data bins against the given thresholds
     */

    if (databinchk[0].param && (strcmp(ptype, "scene") != 0))
        if ((chk_databin(sdfid, databinchk)) < 0) {
            Hclose(fid);
            if (c_fid != 0) Hclose(c_fid);
            stat_exit(stat_status);
        }

    /*
     *  Do a weight check
     */

    if (wtchk)
        if (chk_weight(fid, nbins) < 0) {
            Hclose(fid);
            if (c_fid != 0) Hclose(c_fid);
            stat_exit(stat_status);
        }

    /*
     *  initialize the histogram accumulation
     */
    h_init(argv[1]);

    /*
     *  Do a gross value, general statistics and climatology check
     */

    if ((l3data_chk(argv[3], fid, c_fid, nbins, c_nbins, grschk, statchk,
            climchk)) < 0)
        stat_exit(stat_status);

    /*
     *  Free allocated resources and close files 
     */

    SDend(sdfid);
    Vend(fid);
    Hclose(fid);

    if (c_fid != 0) {
        SDend(c_sdfid);
        Vend(c_fid);
        Hclose(c_fid);
    }

    val3 = (float) clock() / CLOCKS_PER_SEC;
    printf("\n\n# Time took for %s = %d secs\n", argv[0], (int) (val3 - val));

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
      char  *           cntl_file        I      control file path/name 
      cntl_str *        databinchk       O      structure to hold databin
                                                check thresholds 
      cntl_str *        grschk           O      structure to hold control
                                                data related to gross checking
      cntl_str *        statchk          O      structure to hold control
                                                data related to general stat
                                                check
      flag_str *        climchk          O      structure to hold climatology
                                                check thresholds
       
   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      L. Kumar          07/19/1996     Original development

 ****************************************************************************/
int32_t read_cntldata(char *cntl_file, cntl_str *databinchk, cntl_str *grschk,
        cntl_str *statchk, int32_t *wtchk, clim_str *climchk) {
    FILE *fid;
    char line[501], str[25];
    int32_t i, param, nflds = 0;

    printf("\n# Reading control file '%s'\n\n", cntl_file);

    /*
     *  Open the file
     */
    if ((fid = fopen(cntl_file, "r")) == NULL) {
        printf("\n****read_cntldata: Unable to open control file:\n  '%s'\n",
                cntl_file);
        stat_status = stat_status | 1;
        return FAIL;
    }

    /*
     *  Initialize all input structrues
     */

    *wtchk = 0;
    databinchk[0].param = databinchk[0].err_thresh = 0;
    databinchk[0].low_thresh = databinchk[0].high_thresh = 0;
    for (i = 0; i < NPARAMS; i++) {
        grschk[i].param = statchk[i].param = 0;
        grschk[i].err_thresh = statchk[i].err_thresh = 0;
        grschk[i].low_thresh = statchk[i].low_thresh = 0;
        grschk[i].high_thresh = statchk[i].high_thresh = 0;
        climchk[i].param = 0;
        climchk[i].thresh1H = climchk[i].thresh2H = climchk[i].thresh3H = 0;
        climchk[i].thresh1L = climchk[i].thresh2L = climchk[i].thresh3L = 0;
    }

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
        if (strncmp(&line[i], "L3DAT", 5) == 0) {
            if ((nflds = sscanf(line, "%s %f %f", str, &databinchk[0].low_thresh,
                    &databinchk[0].high_thresh)) != 3)
                return (pr_error(str, 2, nflds - 1));
            databinchk[0].param = 1;
        } else if (strncmp(&line[i], "L3WT", 4) == 0)
            *wtchk = 1;
        else if (strncmp(&line[i], "L3GCHK", 6) == 0) {
            param = atoi(&line[6]);
            if (param < 1 || param > NPARAMS) {
                printf("********read_cntldata: Error in parameter number.");
                printf("  Parameter # read = %d", param);
                stat_status = stat_status | 1;
                return FAIL;
            }
            if ((nflds = sscanf(line, "%s %f %f %f", str,
                    &grschk[param - 1].err_thresh, &grschk[param - 1].low_thresh,
                    &grschk[param - 1].high_thresh)) != 4)
                return (pr_error(str, 3, nflds - 1));
            grschk[param - 1].param = 1;
        } else if (strncmp(&line[i], "L3STAT", 6) == 0) {
            param = atoi(&line[6]);
            if (param < 1 || param > NPARAMS) {
                printf("********read_cntldata: Error in param number.");
                printf("  Parameter # read = %d", param);
                stat_status = stat_status | 1;
                return FAIL;
            }
            statchk[param - 1].param = 1;
        } else if (strncmp(&line[i], "L3CCHK", 6) == 0) {
            param = atoi(&line[6]);
            if (param < 1 || param > NPARAMS) {
                printf("********read_cntldata: Error in param number.");
                printf("  Parameter # read = %d", param);
                stat_status = stat_status | 1;
                return FAIL;
            }
            if ((nflds = sscanf(line, "%s %f %f %f %f %f %f", str,
                    &climchk[param - 1].thresh1L, &climchk[param - 1].thresh1H,
                    &climchk[param - 1].thresh2L, &climchk[param - 1].thresh2H,
                    &climchk[param - 1].thresh3L, &climchk[param - 1].thresh3H)) < 7)
                return (pr_error(str, 6, nflds - 1));
            climchk[param - 1].param = 1;
        }
    }
    fclose(fid);
    return SUCCEED;
}

/****************************************************************************
   l3file

   purpose: Verifies if the given input file is level 3 bin data file

   Returns type: integer

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             sdfid           I       SD interface ID
      int32  *          nbins           I       No. of bins containing data

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      L. Kumar          07/19/1996      Original development

 ****************************************************************************/
int32_t l3file(int32_t sdfid, int32_t c_sdfid, int32_t *nbins, int32_t *c_nbins,
        char *ptype) {
    char title[1024];

    /*
     *  Read title to verify if the given input file is level 3 bin data file
     */

    if ((rdattr(sdfid, TITLE, &title)) < 0)
        return FAIL;

    if (strcmp(title, "SeaWiFS Level-3 Binned Data") != 0) {
        printf("l3file: Data file is not level 3 bin file");
        stat_status = stat_status | 1;
        return FAIL;
    }

    if ((rdattr(sdfid, "Data Bins", (VOIDP) nbins)) < 0)
        return FAIL;

    if (c_sdfid != 0) {
        if ((rdattr(c_sdfid, "Data Bins", (VOIDP) c_nbins)) < 0)
            return FAIL;
    }

    if ((rdattr(sdfid, "Product Type", ptype)) < 0)
        return FAIL;

    printf("\n# Data containing bins = %d", *nbins);
    printf("\n\n# Data containing bins in climatology file = %d", *c_nbins);

    return SUCCEED;
}

/****************************************************************************
   chk_databin 

   purpose:  Reads 'Percent Data Bins' global attribute value and verifies 
              it against the given low and high error thresholds and outputs
              the results. 

   Returns type: integer 

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32_t             sdfid           I       SD interface ID
      cntl_str  *       databinchk      I       structure containing bin
                                                percentage thresholds
   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      L. Kumar          07/19/1996      Original development

 ****************************************************************************/
int32_t chk_databin(int32_t sdfid, cntl_str *databinchk) {
    int32_t flag = 0;
    float32 pctbins = 0;

    if ((rdattr(sdfid, "Percent Data Bins", &pctbins)) < 0)
        return FAIL;

    if (pctbins > databinchk[0].high_thresh)
        flag = 1;
    else if (pctbins < databinchk[0].low_thresh)
        flag = -1;

    /*
     *  Print headers for Data bin percentage check results
     */

    printf("\n\n#Mnemonic  Flag  %%data bins   Err_low_thr   Err_high_thr");
    printf("\n#----------------------------------------------------\n");

    printf("\nL3DAT %8d %12.6f %12.6f %12.6f", flag, pctbins,
            databinchk[0].low_thresh, databinchk[0].high_thresh);

    if (flag)
        stat_status = stat_status | 2;

    return SUCCEED;
}

/****************************************************************************
   chk_weight 

   purpose:  Reads 'Percent Data Bins' global attribute value and verifies 
              it against the given low and high error thresholds and outputs
              the results. 

   Returns type: integer 

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32_t             sdfid           I       SD interface ID
      cntl_str  *       databinchk      I       structure containing bin
                                                percentage thresholds
   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      L. Kumar          07/19/1996      Original development

 ****************************************************************************/
int32_t chk_weight(int32_t fid, int32_t nbins) {
    int32_t vsref, vsid, val = 0, flag = 0;
    int32_t start = 0, elts = 0, ret;
    float32 wt_buf[BUFSZ];

    /*  Find Vdata with name "BinList"   */
    if ((vsref = VSfind(fid, "BinList")) < 0) {
        printf("\nchk_weight: VSfind failed for vdata 'BinList'");
        stat_status = stat_status | 1;
        return FAIL;
    }

    /* Attach to a vdata using the given vsid */
    if ((vsid = VSattach(fid, vsref, "r")) < 0) {
        printf("\nchk_weight: VSattach failed for vdata 'BinList'");
        stat_status = stat_status | 1;
        return FAIL;
    }

    start = 0;
    elts = BUFSZ; /* read in BUFSZ elements at a time */
    if (nbins < elts)
        elts = nbins;
    else {
        for (start = 0; start + elts < nbins; start += elts) {
            if ((ret = get_wtcnt(vsid, "weights", start, elts, wt_buf)) < 0)
                return FAIL;
            else
                val += ret;
        }

        if ((elts = (nbins - start)) > 0) {
            if ((ret = get_wtcnt(vsid, "weights", start, elts, wt_buf)) < 0)
                return FAIL;
            else
                val += ret;
        }
    }

    VSdetach(vsid);

    if (val)
        flag = 1;

    /*
     *  Print headers for weight check 
     */

    printf("\n\n#Mnemonic  Flag   Value  ");
    printf("\n#-----------------------\n");
    printf("\nL3WT        %d       %d", flag, val);

    if (flag) {
        printf("\n\n****Suspending further check..");
        printf(" %d bins have been found with weight set to 0.", val);
        stat_status = stat_status | 2;
        return FAIL;
    }
    else
        return SUCCEED;
}

/****************************************************************************
   l3data_chk

   purpose: 

   Returns type: integer

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char      *       clim_file       I       Climatology file used
      int32_t             fid             I       SD interface ID
      int32_t             c_fid           I       SD interface ID for climatology
                                                file
      cntl_str  *       grschk          I       structure containing Gross
                                                value control thresholds
      cntl_str  *       statchk         I       struct containing statistical
                                                value control thresholds
      clim_str  *       climchk         I       structure containing thresholds
                                                for climatology checks

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      L. Kumar          07/19/1996      Original development

 ****************************************************************************/

/*
 ** data      -  float array of size 'BUFSZ'(1000) for reading parameter
 **               data
 ** wts       -  float array of size 'BUFSZ' for reading weights 
 ** elts      -  number of elements (records)
 ** sumxbar   -  an array(NAPRAMS) of type double for accumulating individual
 **              mean values
 ** sumsd     -  an array(NPARAMS) of type double for accumulating individual
 **              standard deviation values
 */

int32_t l3data_chk(char *clim_file, int32_t fid, int32_t c_fid, int32_t nbins,
        int32_t c_nbins, cntl_str *grschk, cntl_str *statchk,
        clim_str *climchk) {
    int32_t i, ci, p, j, done, c_done, rd_flag, rd_c_flag;
    int32_t st, c_st, elts, c_elts;
    int32_t prev_c_st = -1;
    int32_t binlist_id, c_binlist_id, proceed = 0;
    int32_t param_vsid[NPARAMS], c_param_vsid[NPARAMS];
    int32_t binno[BUFSZ], c_binno[BUFSZ];
    float32 data[BUFSZ][2];
    float32 c_data[NPARAMS][BUFSZ][2];
    float32 wts[BUFSZ], c_wts[BUFSZ];
    int16 scenes[BUFSZ], c_scenes[BUFSZ];
    float64 xbar[NPARAMS][BUFSZ];
    float64 c_xbar = 0, c_sd = 0;
    float64 npts = 0;
    float64 loparm[NPARAMS] = {1.e10, 1.e10, 1.e10, 1.e10, 1.e10, 1.e10,
        1.e10, 1.e10, 1.e10, 1.e10, 1.e10};
    float64 hiparm[NPARAMS] = {-1.e10, -1.e10, -1.e10, -1.e10, -1.e10,
        -1.e10, -1.e10, -1.e10, -1.e10, -1.e10, -1.e10};
    float64 negflag[NPARAMS] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float64 num1Lstd[NPARAMS] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float64 num1Hstd[NPARAMS] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float64 num2Lstd[NPARAMS] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float64 num2Hstd[NPARAMS] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float64 num3Lstd[NPARAMS] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float64 num3Hstd[NPARAMS] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float64 value1[NPARAMS] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float64 value2[NPARAMS] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float64 sumxbar[NPARAMS] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float64 sumsd[NPARAMS] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    printf("\n\n#  Doing gross, statistical and climatology checks\n");

    /*  
     * return to calling routine, if there are no gross, statistical or
     * climatological checks requested
     */
    for (p = 0; p < NPARAMS; p++)
        if (grschk[p].param || statchk[p].param || climchk[p].param)
            proceed = 1;

    if (!proceed)
        return SUCCEED;

    /*
     * get access id for accessing binlist and parameter vdatas of both
     * input binfile and climatology bin file
     */

    binlist_id = get_vsid(fid, "BinList");
    if (c_fid != 0)c_binlist_id = get_vsid(c_fid, "BinList");

    for (p = 0; p < NPARAMS; p++) {
        if ((param_vsid[p] = get_vsid(fid, param_name[p])) < 0)
            return FAIL;

        if (c_fid != 0) {
            if ((c_param_vsid[p] = get_vsid(c_fid, param_name[p])) < 0)
                return FAIL;
        }
    }

    /* 
     * Set initial values for the following flags and variables, 
     * do the loop until all the bins of input file have been read
     * and verified for gross and statistical checks
     */
    rd_flag = rd_c_flag = 1;

    st = c_st = 0; /* this tracks the index of the bin (out of total #
                        in the file) that is at the startof current buffer */
    elts = c_elts = BUFSZ; /* # elements that were read to the buffer */

    if (nbins < BUFSZ)
        elts = nbins;

    if (c_nbins < BUFSZ)
        c_elts = c_nbins;

    done = c_done = 0;
    while (!done) {
        if (rd_flag) {
            /* Read in next set of bin #s weights and # scenes from the file */
            rd_flag = 0;
            i = 0; /* index into binno buffer */
            if ((rdvdata(binlist_id, "bin_num", st, elts,
                    (unsigned char *) binno)) < 0)
                return FAIL;
            if ((rdvdata(binlist_id, "weights", st, elts,
                    (unsigned char *) wts)) < 0)
                return FAIL;
            if ((rdvdata(binlist_id, "nscenes", st, elts,
                    (unsigned char *) scenes)) < 0)
                return FAIL;

            for (p = 0; p < NPARAMS; p++) /* Also, read in each parameter */ {
                if ((rdvdata(param_vsid[p], param_flds[p], st, elts,
                        (unsigned char *) data)) < 0)
                    return FAIL;

                /* For all the elements in the buffer, compute the mean */
                for (j = 0; j < elts; j++) {
                    xbar[p][j] = data[j][0] / wts[j];
                    sumxbar[p] += xbar[p][j];
                    sumsd[p] += xbar[p][j] * xbar[p][j];
                    /*
                     * get the low and high
                     */
                    if (xbar[p][j] < loparm[p]) loparm[p] = xbar[p][j];
                    if (xbar[p][j] > hiparm[p]) hiparm[p] = xbar[p][j];
                    /*
                     * add to the histograms
                     */
                    h_accum(xbar[p][j], p);

                    if (xbar[p][j] < grschk[p].low_thresh)
                        value1[p]++;
                    if (xbar[p][j] > grschk[p].high_thresh) {
                        value2[p]++;
#ifdef DEBUG
                        printf("\n#indx= %d,sum= %f,sumxx= %f,wts= %f,mean= %f",
                                j, data[j][0], data[j][1], wts[j], xbar[p][j]);
#endif
                    }
                }
            }
        }
        if (c_fid != 0) { /* WDR added climatology omission if */
            /* read the bin #s for the climatology if necessary (if any of the
               climatology data is needed for this buffer, it will be read 
               further below if required */
            if (rd_c_flag) {
                rd_c_flag = 0;
                ci = 0; /* index into clim. binno buffer */
                if ((rdvdata(c_binlist_id, "bin_num", c_st, c_elts,
                        (unsigned char *) c_binno)) < 0)
                    return FAIL;
            }

            while (i < elts && ci < c_elts) {
                if (binno[i] < c_binno[ci])
                    i++;
                else if (binno[i] > c_binno[ci])
                    ci++;
                else {
                    npts++;
                    if (prev_c_st != c_st) {
                        /*  The rest of the climatology data is read for these
                            bins so that comparisons can be done  */
                        if ((rdvdata(c_binlist_id, "weights", c_st, c_elts,
                                (unsigned char *) c_wts)) < 0)
                            return FAIL;
                        if ((rdvdata(c_binlist_id, "nscenes", c_st, c_elts,
                                (unsigned char *) c_scenes)) < 0)
                            return FAIL;
                        for (p = 0; p < NPARAMS; p++)
                            if ((rdvdata(c_param_vsid[p], param_flds[p], c_st, c_elts,
                                    (unsigned char *) c_data[p])) < 0)
                                return FAIL;
                        prev_c_st = c_st;
                    }
                    for (p = 0; p < NPARAMS; p++) {
                        /* get the mean, std deviation for climatology and check the 
                           file mean against it */
                        c_xbar = c_data[p][ci][0] / c_wts[ci];
                        c_sd = calc_sd(c_xbar, c_wts[ci], c_scenes[ci],
                                c_data[p][ci][1]);

                        if (c_sd == NEG_FLAG)
                            negflag[p]++;
                        else {
                            if (xbar[p][i] < (c_xbar - 3 * c_sd)) {
                                num3Lstd[p]++;
                                num2Lstd[p]++;
                                num1Lstd[p]++;
                            } else if (xbar[p][i] < (c_xbar - 2 * c_sd)) {
                                num2Lstd[p]++;
                                num1Lstd[p]++;
                            } else if (xbar[p][i] < (c_xbar - c_sd))
                                num1Lstd[p]++;
                            if (xbar[p][i] > (c_xbar + 3 * c_sd)) {
                                num3Hstd[p]++;
                                num2Hstd[p]++;
                                num1Hstd[p]++;
                            } else if (xbar[p][i] > (c_xbar + 2 * c_sd)) {
                                num2Hstd[p]++;
                                num1Hstd[p]++;
                            } else if (xbar[p][i] > (c_xbar + c_sd))
                                num1Hstd[p]++;
                        }
                    }
                    i++;
                    ci++;
                }
            }

            /*  if the end of the file data or climatology data buffer is hit,
                set the proper flag so that the next stretch of data is read in */

            if (ci == c_elts && c_elts != 0) {
                c_st += c_elts;
                if (c_st + c_elts > c_nbins)
                    c_elts = c_nbins - c_st;
                if (c_elts > 0)
                    rd_c_flag = 1;
                else
                    c_done = 1;
            }
        }/* WDR added climatology omission if and add else condition below */
        else {
            i = elts;
        }

        if ((i == elts && elts != 0) || c_done) {
            st += elts;
            if (st + elts > nbins)
                elts = nbins - st;
            if (elts > 0)
                rd_flag = 1;
            else
                done = 1;
        }
    }
    VSdetach(binlist_id);
    if (c_fid != 0) VSdetach(c_binlist_id);

    pr_grs_results(grschk, nbins, value1, value2);
    pr_stat_results(statchk, nbins, sumxbar, sumsd, loparm, hiparm);

    /*
     *  output the histograms
     */
    h_out(nbins, sumxbar, sumsd, loparm, hiparm);

    for (p = 0; p < NPARAMS; p++) {
        VSdetach(param_vsid[p]);
        if (c_fid != 0) VSdetach(c_param_vsid[p]);
    }

    if (c_fid != 0) pr_clim_results(clim_file, climchk, npts, num1Lstd, num1Hstd,
            num2Lstd, num2Hstd, num3Lstd, num3Hstd);

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

    Returns:   int32_t (status)
        The return code is a negative value if any error occurs, otherwise,
        returns 0.

    Description:
        The function rdattr reads the requested global attribute

    Parameters: (in calling order)
        Type      Name        I/O   Description
        ----      ----        ---   -----------
        int32_t     sdfid        I    ID req to access HDF SDS interface
        char  *   attr_name    I    attribute name
        void  *   buf         I/O   pointer to data buffer

    Modification history:
    Programmer     Organization   Date        Description of change
    -------------- ------------   --------    ---------------------
    Lakshmi Kumar  Hughes STX     11/07/94    original development

----------------------------------------------------------------------------*/
int32_t rdattr(int32_t sdfid, char *attr_name, void *buf) {
    int32_t attrnum;

    if ((attrnum = SDfindattr(sdfid, attr_name)) < 0) {
        printf("\n****rdattr: Failure in SDfindattr while trying to read %s",
                attr_name);
        stat_status = stat_status | 1;
        return FAIL;
    }

    if ((SDreadattr(sdfid, attrnum, buf)) < 0) {
        printf("\n****rdattr: Failure in SDreadattr while trying to read %s",
                attr_name);
        stat_status = stat_status | 1;
        return FAIL;
    }

    return SUCCEED;
}

/*----------------------------------------------------------------------------
    Function: pr_error

    Description:
        The function prints an error message and returns a negetive 1 to the
        calling routine
-----------------------------------------------------------------------------*/
int32_t pr_error(char *label, int32_t nvals, int32_t nvals_read) {
    printf("\n*****Expecting %d value(s) for - %s, but, read %d value(s). \n",
            nvals, label, nvals_read);
    stat_status = stat_status | 1;
    return FAIL;
}

/*-----------------------------------------------------------------------------
    Function: get_wtcnt

    Returns: int32_t(no. of wts that are less than or equal to zero)

    Description:
        The function get_wtcnt reads weights from BinIndex vdata starting 
        at position indicated by start, up to elements indicated by 'elts' 
        and returns a count of weights that are less than or equal to zero.

    Arguments: (in calling order)
      Type         Name      I/O     Description
      ----         ----      ---     -----------
      int32_t        vsid       I      access id for accessing BinIndex vdata
      char *       fld_name   I      field name (weights)
      int32_t        start      I      start element position
      int32_t        elts       I      number of elements to read
      float32 *    wt_buf     O      buffer to read the data

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      07/19/1996  Original development
------------------------------------------------------------------------------*/
int32_t get_wtcnt(int32_t vsid, char *fld_name, int32_t start, int32_t elts,
        float32 *wt_buf) {
    int32_t val = 0, i;

    if (rdvdata(vsid, fld_name, start, elts, (unsigned char *) wt_buf) < 0)
        return FAIL;

    for (i = 0; i < elts; i++)
        if (wt_buf[i] <= 0)
            val++;

    return val;
}

/*----------------------------------------------------------------------------
    Function: get_vsid 

    Returns: int32_t (vdata id)

    Description:
        The function get_vsid locates the given named vdata, sets read
         access to it and returns the access id to the calling routine.

    Arguments: (in calling order)
      Type         Name      I/O     Description
      ----         ----      ---     -----------
      int32_t        fid        I      File access ID
      char *       vs_name    I      vdata name

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      07/19/1996  Original development

-----------------------------------------------------------------------------*/
int32_t get_vsid(int32_t fid, char *vs_name) {

    int32_t vsid, vsref;

    /* Get vdata reference number for the given vdata name */
    if ((vsref = VSfind(fid, vs_name)) < 0) {
        printf("\nget_vsid: VSfind failed for vdata '%s'", vs_name);
        stat_status = stat_status | 1;
        return FAIL;
    }

    /* Attach to a vdata using the given vdata reference number */
    if ((vsid = VSattach(fid, vsref, "r")) < 0) {
        printf("\nget_vsid: VSattach failed for vdata '%s'", vs_name);
        stat_status = stat_status | 1;
        return FAIL;
    }

    return vsid;
}

/*----------------------------------------------------------------------------
    Function: calc_sd

    Returns: none 

    Description:

    Arguments: (in calling order)
      Type         Name      I/O     Description
      ----         ----      ---     -----------

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      07/19/1996  Original development

-----------------------------------------------------------------------------*/
float64 calc_sd(float64 xmean, float32 wts, int16 nseg, float32 sumxx) {
    float64 ww, sd;

    ww = wts * wts;
    if (ww > nseg) {
        sd = (sumxx / wts) - xmean * xmean;
        if (sd < 0)
            sd = 0;
        else
            sd = sqrt(sd * ww / (ww - nseg));
    } else
        sd = NEG_FLAG;

    return sd;
}

/*----------------------------------------------------------------------------
    Function: pr_grs_results

    Returns: none 

    Description:
        The function pr_grs_results outputs the gross check results 

    Arguments: (in calling order)
      Type         Name      I/O     Description
      ----         ----      ---     -----------
      cntl_str   * grschk     I      strcutre containing gross check thresholds
      int32_t        nbins      I      number of bins present in input file
      float64    * value1     I      number of data pts (mean) < low thresh
      float64    * value2     I      number of data pts (mean) > high thresh

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      07/19/1996  Original development

-----------------------------------------------------------------------------*/
void pr_grs_results(cntl_str *grschk, int32_t nbins, float64 *value1,
        float64 *value2) {
    int32_t p, flag, failcode = 0;
    float64 v1, v2;

    /*
     *  Print headers for gross check output messages
     */

    printf("\n\n#Mnemonic  Flag    %%<Low       %%>High       Err_thr     ");
    printf("low_thr     high_thr");
    printf("\n#-----------------------------------------------------------");
    printf("-------------------\n");

    for (p = 0; p < NPARAMS; p++)
        if (grschk[p].param) {
            flag = 0;
            v1 = value1[p] / nbins * 100;
            v2 = value2[p] / nbins * 100;
            if (v1 > grschk[p].err_thresh)
                flag = 1;
            else if (v2 > grschk[p].err_thresh)
                flag = -1;

            printf("\nL3GCHK%d %5d %11.6f %12.6f %12.6f %12.6f %12.6f",
                    p + 1, flag, v1, v2, grschk[p].err_thresh,
                    grschk[p].low_thresh, grschk[p].high_thresh);
            if (flag != 0)
                failcode = 1;
        }
    if (failcode)
        stat_status = stat_status | 2;
}

/*----------------------------------------------------------------------------
    Function: pr_stat_results

    Returns: none 

    Description:
        The function pr_stat_results calculates mean and std. deviation 
        of the sumx's and sumsxx's and outputs the results

    Arguments: (in calling order)
      Type         Name      I/O     Description
      ----         ----      ---     -----------
      cntl_str   * statchk    I      structure containing statistical check
                                      thresholds
      int32_t        nbins      I      number of bins present
      float64    * sumx      I      array of sz 12, containing sum of sumx's
      float64    * sumxx     I      array of sz 12, containing sum of sumxx's
      float64    * loparm    I      lowest value encountered
      float64    * hiparm    I      highest value encountered

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      07/19/1996  Original development

-----------------------------------------------------------------------------*/
void pr_stat_results(cntl_str *statchk, int32_t nbins, float64 *sumxbar,
        float64 *sumsd, float64 *loparm, float64 *hiparm) {
    int32_t p, flag = 0;
    float64 sum, sumsq, xbar, sd;

    /*
     *  Print headers for statistical check output messages
     */

    printf("\n\n");
    printf("#Mnemonic  flag   nbins       xbar         sd   Low   High ");
    printf("\n#----------------------------------------------------\n");

    for (p = 0; p < NPARAMS; p++)
        if (statchk[p].param) {
            sum = sumxbar[p] / nbins;
            sumsq = sumsd[p] / nbins;
            xbar = sum;
            sd = sqrt(sumsq - (sum * sum));

            printf("\nL3STAT%d %5d %7d %12.9f %12.9f %f  %f", p + 1, flag, nbins,
                    xbar, sd, loparm[p], hiparm[p]);
        }
}

/*----------------------------------------------------------------------------
    Function: pr_clim_results

    Returns: none 

    Description:
        The function pr_clim_results calculates percentage outside 1, 2,
        and 3 standard deviation and outputs the results

    Arguments: (in calling order)
      Type         Name      I/O     Description
      ----         ----      ---     -----------
      char       * clim_file  I      name of climatology file used
      clim_str   * climchk    I      structure containing climatology check
                                     thresholds
      int32_t        npts       I      number of input bins that matched clim.
                                     bins
      float64    * num1Hstd    I      number of bins outside 1 STD
      float64    * num1Lstd    I      number of bins outside 1 STD
      float64    * num2Hstd    I      number of bins outside 2 STD
      float64    * num2Lstd    I      number of bins outside 2 STD
      float64    * num3Hstd    I      number of bins outside 3 STD
      float64    * num3Lstd    I      number of bins outside 3 STD

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      07/19/1996  Original development
        Lakshmi Kumar    Hughes STX      06/16/1997  Added clim_file argument
-----------------------------------------------------------------------------*/
void pr_clim_results(char *clim_file, clim_str *climchk, float64 npts,
        float64 *num1Hstd, float64 *num1Lstd, float64 *num2Hstd,
        float64 *num2Lstd, float64 *num3Hstd, float64 *num3Lstd) {
    int32_t p, failcode = 0;
    int32_t flag1_H, flag2_H, flag3_H;
    int32_t flag1_L, flag2_L, flag3_L;
    float64 pct1Hstd[NPARAMS], pct2Hstd[NPARAMS], pct3Hstd[NPARAMS];
    float64 pct1Lstd[NPARAMS], pct2Lstd[NPARAMS], pct3Lstd[NPARAMS];

    for (p = 0; p < NPARAMS; p++) {
        pct1Lstd[p] = pct2Lstd[p] = pct3Lstd[p] = 0;
        pct1Hstd[p] = pct2Hstd[p] = pct3Hstd[p] = 0;
    }

    /*
     *  Print headers for statistical check output messages
     */

    printf("\n\nCLIM_FILE_USED: %s", clim_file);
    printf("\n\nCLIM_NPTS: %f\n", npts);
    printf("\n\n");
    printf("#Mnemonic  flag     #std dev   %%outside    %%thresh ");
    printf("\n#----------------------------------------------------\n");

    for (p = 0; p < NPARAMS; p++)
        if (climchk[p].param) {
            if (npts > 0) {
                pct1Hstd[p] = 100 * num1Hstd[p] / npts;
                pct1Lstd[p] = 100 * num1Lstd[p] / npts;
                pct2Hstd[p] = 100 * num2Hstd[p] / npts;
                pct2Lstd[p] = 100 * num2Lstd[p] / npts;
                pct3Hstd[p] = 100 * num3Hstd[p] / npts;
                pct3Lstd[p] = 100 * num3Lstd[p] / npts;
            }

            flag1_H = flag1_L = flag2_H = flag2_L = flag3_H = flag3_L = 0;

            if (pct3Lstd[p] > climchk[p].thresh3L)
                flag3_L = 1;
            if (pct2Lstd[p] > climchk[p].thresh2L)
                flag2_L = 1;
            if (pct1Lstd[p] > climchk[p].thresh1L)
                flag1_L = 1;

            if (pct3Hstd[p] > climchk[p].thresh3H)
                flag3_H = 1;
            if (pct2Hstd[p] > climchk[p].thresh2H)
                flag2_H = 1;
            if (pct1Hstd[p] > climchk[p].thresh1H)
                flag1_H = 1;

            if (flag1_H || flag1_L || flag2_H || flag2_L || flag3_H || flag3_L)
                failcode = 1;

            printf("\n\nL3CCHK%d %6d %10s %12.6lf %12.6f",
                    p + 1, flag1_L, "LO_1", pct1Lstd[p], climchk[p].thresh1L);
            printf("\nL3CCHK%d %6d %10s %12.6lf %12.6f",
                    p + 1, flag1_H, "HI_1", pct1Hstd[p], climchk[p].thresh1H);
            printf("\nL3CCHK%d %6d %10s %12.6f %12.6f",
                    p + 1, flag2_L, "LO_2", pct2Lstd[p], climchk[p].thresh2L);
            printf("\nL3CCHK%d %6d %10s %12.6f %12.6f",
                    p + 1, flag2_H, "HI_2", pct2Hstd[p], climchk[p].thresh2H);
            printf("\nL3CCHK%d %6d %10s %12.6f %12.6f",
                    p + 1, flag3_L, "LO_3", pct3Lstd[p], climchk[p].thresh3L);
            printf("\nL3CCHK%d %6d %10s %12.6f %12.6f",
                    p + 1, flag3_H, "HI_3", pct3Hstd[p], climchk[p].thresh3H);
        }
    if (failcode)
        stat_status = stat_status | 2;
}


