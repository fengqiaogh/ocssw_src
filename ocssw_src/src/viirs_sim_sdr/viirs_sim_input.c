#include "viirs_sim_sdr.h"
#include <stdio.h>
#include <assert.h>
#include "clo.h"
#include <ctype.h>

int viirs_sim_input(int argc, char *argv[], ctl_struc *ctl)
/*-----------------------------------------------------------------------------
    Routine:   viirs_sim_input

    Description:  get inputs for viirs_sim_sdr

    Returns:  int 0 if all good

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       argc          I    count of command line args 
        char *[]  argv          I    command line arguments to be processed
                                     by the clo library

    Modification history:

    W. Robinson, SAIC  21 Sep 2010  Original development

----------------------------------------------------------------------------*/ {
    int errflg, opt_count, opt_count2, verr;
    char tmpstr[FILENAME_MAX], filetmp[FILENAME_MAX];
    int i, nbands, ar_count, orb_geo_use = 0, out_loc_use = 0;
    clo_optionList_t* list;
    float *fArray;
    char localsuite[FILENAME_MAX];
    char prog_name[] = "viirs_sim_sdr", sensor_id[] = "viirsn", *dataroot;

    if ((dataroot = getenv("OCDATAROOT")) == NULL) {
        printf("%s, %d -E- OCDATAROOT environment variable is not defined.\n",
                __FILE__, __LINE__);
        return (-1);
    }
    /*
     *  follow the clo.h and start it up
     */
    sprintf(tmpstr, "viirs_sim_sdr version 1.0 (%s %s)", __DATE__, __TIME__);
    clo_setVersion(tmpstr);
    list = clo_createList();

    assert(list);

    /*
     *  fill some control values
     */
    ctl->vic_cal_chg = 0;
    for (i = 0; i < MAX_BND; i++) {
        ctl->gain[i] = 1.;
        ctl->offset[i] = 0;
    }
    /*
     *  set default arguments and add callback for a par file
     */
    strcpy(tmpstr, "Input VIIRS geolocation file (required)");
    clo_addOption(list, "geofile", CLO_TYPE_IFILE, NULL, tmpstr);

    strcpy(tmpstr, "path to place output SDR files (required)");
    clo_addOption(list, "opath", CLO_TYPE_OFILE, NULL, tmpstr);

    ctl->l2_use = 0;
    strcpy(tmpstr, "level-2 file containing TOA radiances");
    clo_addOption(list, "il2file", CLO_TYPE_IFILE, "Unspecified", tmpstr);

    ctl->rhos_use = 0;
    strcpy(tmpstr, "Land reflectance data file");
    clo_addOption(list, "iland_refl_file", CLO_TYPE_IFILE, "Unspecified",
            tmpstr);

    strcpy(tmpstr, "reflectance replacement options:\n");
    strcat(tmpstr,
            "\t0: fill only where there are no valid TOA values\n");
    strcat(tmpstr, "\t1: replace all points with reflectance taken to TOA");
    clo_addOption(list, "refl_fil_opt", CLO_TYPE_INT, "0", tmpstr);

    strcpy(tmpstr, "output scan format, resultant format is\n");
    strcat(tmpstr, "\tminimum of input format and oscn_fmt_opt:\n");
    strcat(tmpstr, "\t0 - make aggregated\n");
    strcat(tmpstr, "\t1 - make unaggregated\n");
    strcat(tmpstr, "\t2 - same as input format");
    clo_addOption(list, "oscn_fmt_opt", CLO_TYPE_INT, "2", tmpstr);

    strcpy(tmpstr, "M-band create option:\n");
    strcat(tmpstr, "\t0 - make the Vis, NIR bands (M1 - 7)\n");
    strcat(tmpstr, "\t1 - make all M bands (M1 - 16)");
    clo_addOption(list, "band_opt", CLO_TYPE_INT, "0", tmpstr);

    strcpy(tmpstr, "optical crosstalk option:\n");
    strcat(tmpstr, "\t0 - no added cross talk\n");
    strcat(tmpstr, "\t1 - add optical crosstalk (OXT) artifact");
    clo_addOption(list, "oxt_opt", CLO_TYPE_INT, "0", tmpstr);

    strcpy(filetmp, "$VIIRS_SIM_DATA/oxt_coeff_default.h5\n");
    sprintf(tmpstr, "optical crosstalk influence coefficients");
    clo_addOption(list, "oxt_coef_file", CLO_TYPE_IFILE, filetmp, tmpstr);

    strcpy(filetmp, "$VIIRS_SIM_DATA/inter_band_default.dat");
    sprintf(tmpstr, "Inter-band data table\n");
    clo_addOption(list, "inter_band_file", CLO_TYPE_IFILE, filetmp, tmpstr);

    strcpy(tmpstr, "metadata output file (not implemented yet)");
    clo_addOption(list, "metafile", CLO_TYPE_OFILE, "Unspecified", tmpstr);

    strcpy(tmpstr,
            "bow tie deletion option (default, on to do bow tie delete)");
    clo_addOption(list, "bow_tie_opt", CLO_TYPE_BOOL, "on", tmpstr);

    strcpy(tmpstr, "Calibration gains to remove from TOA radiances");
    clo_addOption(list, "gain", CLO_TYPE_FLOAT, "1.,1.", tmpstr);

    strcpy(tmpstr, "Calibration offsets to remove from TOA radiances");
    clo_addOption(list, "offset", CLO_TYPE_FLOAT, "0.,0.", tmpstr);

    strcpy(tmpstr, "optional create time in form YYYYMMDD/HHMMSS.FFFFFF\n");
    strcat(tmpstr, "\tif not specified, the time the program ran");
    clo_addOption(list, "cre_time", CLO_TYPE_STRING, "Unspecified", tmpstr);

    strcpy(tmpstr, "SDR file name create option:\n");
    strcat(tmpstr, "\t0 - make standard VIIRS format\n");
    strcat(tmpstr, "\t1 - omit the create time part (_cYYYY...)\n");
    clo_addOption(list, "fname_opt", CLO_TYPE_INT, "0", tmpstr);

    strcpy(tmpstr, "SDR overwrite option:\n");
    strcat(tmpstr, "\t0 - do not permit overwrite\n");
    strcat(tmpstr, "\t1 - allow overwrite to happen\n");
    clo_addOption(list, "sdr_overwrite", CLO_TYPE_INT, "0", tmpstr);

    /*
     *  the count calibration (decalibration) controls
     */
    ctl->count_cal_opt = 0;
    strcpy(tmpstr, "count calibration option:\n");
    strcat(tmpstr, "\t0 - no count calibration done\n");
    strcat(tmpstr, "\t1 - perform de-cal and cal (mandatory for exect X-talk\n");
    strcat(tmpstr, "\t2 - (1), but also integerize counts");
    clo_addOption(list, "count_cal_opt", CLO_TYPE_INT, "0", tmpstr);

    strcpy(tmpstr, "Count calibration gain file");
    strcat(tmpstr, "\t(if omitted, an internal TVAC general set is used)");
    clo_addOption(list, "count_cal_gain_file", CLO_TYPE_IFILE, "Unspecified",
            tmpstr);

    strcpy(tmpstr, "Count calibration rvs file leader (pre-lambda)");
    strcat(tmpstr, "\t(if omitted, unity rvs is used)");
    clo_addOption(list, "count_cal_rvs_file", CLO_TYPE_IFILE, "Unspecified",
            tmpstr);

    strcpy(tmpstr, "Count de-calibration gain file\n");
    strcat(tmpstr, "\t(if omitted, an internal TVAC general set is used)");
    clo_addOption(list, "count_decal_gain_file", CLO_TYPE_IFILE, "Unspecified",
            tmpstr);

    strcpy(tmpstr, "Count de-calibration rvs file (pre-lambda)\n");
    strcat(tmpstr, "\t(if omitted, unity rvs is used)");
    clo_addOption(list, "count_decal_rvs_file", CLO_TYPE_IFILE, "Unspecified",
            tmpstr);

    ctl->count_dark_opt = 0;
    strcpy(tmpstr, "count calibration dark subtract option:\n");
    strcat(tmpstr,
            "\t0 - do not perform dark subtract in decal / cal process\n");
    strcat(tmpstr, "\t1 - do the dark subtract in decal / cal process\n");
    clo_addOption(list, "count_dark_opt", CLO_TYPE_INT, "0", tmpstr);

    strcpy(tmpstr, "product suite string for loading\n");
    strcat(tmpstr, "        suite-specific defaults");
    clo_addOption(list, "suite", CLO_TYPE_STRING, "EMPTY", tmpstr);

    ctl->ext_opt = 0;
    strcpy(tmpstr, "Electronic crosstalk processing option:\n");
    strcat(tmpstr, "\t0 - no crosstalk applied\n");
    strcat(tmpstr, "\t1 - Apply electronic crosstalk\n");
    clo_addOption(list, "ext_opt", CLO_TYPE_INT, "0", tmpstr);

    strcpy(tmpstr, "Electronic crosstalk coefficient file\n");
    strcat(tmpstr,
            "\t(if omitted, a trivial crosstalk, zero influence, is applied)\n");
    clo_addOption(list, "ext_coeff_file", CLO_TYPE_IFILE, "Unspecified",
            tmpstr);

    strcpy(tmpstr, "SDR origin designation\n");
    strcat(tmpstr,
            "\t(will appear in attributes and file name origin location)\n");
    clo_addOption(list, "id_origin", CLO_TYPE_STRING, "gsfc", tmpstr);

    strcpy(tmpstr, "SDR domain designation\n");
    strcat(tmpstr,
            "\t(will appear in attributes and file name domain location)\n");
    clo_addOption(list, "id_domain", CLO_TYPE_STRING, "SCI", tmpstr);

    /*
     *  for the noise, control if used and have file of noise coeffs
     */
    ctl->noise_mode = 0;
    strcpy(tmpstr, "Noise addition option: 0 to not use 1 to apply\n");
    clo_addOption(list, "noise_opt", CLO_TYPE_INT, "0", tmpstr);

    strcpy(tmpstr, "Noise coefficient file name\n");
    strcat(tmpstr, "\t(if omitted or Unspecified, no noise will be added\n");
    clo_addOption(list, "noise_coeff", CLO_TYPE_IFILE, "Unspecified", tmpstr);
    /*
     *  for stray light artifact, control if used and have file of coeffs
     */
    ctl->stray_opt = 0;
    strcpy(tmpstr, "Stray light addition option: 0 to not use 1 to apply\n");
    clo_addOption(list, "stray_opt", CLO_TYPE_INT, "0", tmpstr);

    strcpy(tmpstr, "Stray light coefficient file name\n");
    strcat(tmpstr, "\t(if omitted or Unspecified, no stray light is added\n");
    clo_addOption(list, "stray_tbl", CLO_TYPE_IFILE, "Unspecified", tmpstr);
    /*
     *  All code to get data from the default file, arguments, suite file...
     *
     *
     */
    /*
     *  get the option info and re-get so the command line set takes precedence
     */
    clo_setEnableDumpOptions(0);
    opt_count = clo_getNumOptions(list);
    /*
     *  go through the option heirarchy
     */
    clo_readArgs(list, argc, argv);

    localsuite[0] = '\0';
    if (clo_isSet(list, "suite")) {
        strcpy(localsuite, clo_getString(list, "suite"));
    }
    sprintf(tmpstr, "%s/%s/%s_defaults.par", dataroot, sensor_id, prog_name);
    clo_readFile(list, tmpstr);

    if (localsuite[0] == '\0')
        strcpy(localsuite, clo_getString(list, "suite"));
    sprintf(tmpstr, "%s/%s/%s_defaults_%s.par", dataroot, sensor_id, prog_name,
            localsuite);
    clo_readFile(list, tmpstr);

    clo_setEnableDumpOptions(1);
    clo_readArgs(list, argc, argv);
    /*
     *  end heirarchy, check for any new options
     */
    opt_count2 = clo_getNumOptions(list);
    if (opt_count != opt_count2) {
        printf("\n\nError: initial # options: %d, final # options: %d\n",
                opt_count, opt_count2);
        printf("%s, %d: Error, a mis-spelled or non-existant option has\n",
                __FILE__, __LINE__);
        printf("\tbeen entered.  The extra options should be at the end of\n");
        printf("\tthe following option list:\n");
        clo_printUsage(list);
        exit(1);
    }
    /*
     *  End all code to get data from the default file, arguments, suite file...
     *
     *
     */
    /*
     *  gather the inputs to the ctl struct and check input validity
     */
    errflg = 0;

    if (clo_isSet(list, "geofile")) {
        orb_geo_use = 1;
        strcpy(ctl->in_geo_file, clo_getString(list, "geofile"));
    }

    if (clo_isSet(list, "opath")) {
        out_loc_use = 1;
        strcpy(ctl->out_loc, clo_getString(list, "opath"));
    }

    if (clo_isSet(list, "il2file")) ctl->l2_use = 1;
    strcpy(ctl->l2_file, clo_getString(list, "il2file"));

    if (clo_isSet(list, "iland_refl_file")) ctl->rhos_use = 1;
    strcpy(ctl->rhos_file, clo_getString(list, "iland_refl_file"));

    ctl->rhos_opt = clo_getInt(list, "refl_fil_opt");
    if ((ctl->rhos_opt < 0) || (ctl->rhos_opt > 1)) {
        printf("refl_fil_opt must be 0 or 1\n");
        errflg++;
    }

    ctl->out_scn_fmt = clo_getInt(list, "oscn_fmt_opt");
    if ((ctl->out_scn_fmt < 0) || (ctl->out_scn_fmt > 2)) {
        printf("out_scn_fmt must be 0 to 2\n");
        errflg++;
    }

    ctl->make_m = clo_getInt(list, "band_opt");
    if ((ctl->make_m < 0) || (ctl->make_m > 1)) {
        printf("make_m must be 0 or 1\n");
        errflg++;
    }
    if (ctl->make_m == 0)
        nbands = 7;
    else
        nbands = 16;

    ctl->oxt_mode = clo_getInt(list, "oxt_opt");
    if ((ctl->oxt_mode < 0) || (ctl->oxt_mode > 1)) {
        printf("oxt_mode must be 0 to 1\n");
        errflg++;
    }

    strcpy(ctl->oxt_coef, clo_getString(list, "oxt_coef_file"));

    strcpy(ctl->inter_band, clo_getString(list, "inter_band_file"));

    ctl->meta_use = (clo_isSet(list, "metafile")) ? 1 : 0;
    strcpy(ctl->meta_file, clo_getString(list, "metafile"));

    ctl->bowtie_opt = (clo_getBool(list, "bow_tie_opt")) ? 1 : 0;

    ctl->sdr_overwrite = clo_getInt(list, "sdr_overwrite");
    if ((ctl->sdr_overwrite != 0) && (ctl->sdr_overwrite != 1)) {
        printf("sdr_overwrite must be 0 or 1\n");
        errflg++;
    }

    ctl->fname_opt = clo_getInt(list, "fname_opt");
    if ((ctl->fname_opt != 0) && (ctl->fname_opt != 1)) {
        printf("fname_opt must be 0 or 1\n");
        errflg++;
    }

    if (clo_isSet(list, "gain")) {
        ctl->vic_cal_chg = 1;
        fArray = clo_getFloats(list, "gain", &ar_count);
        if (ar_count != nbands) {
            printf("%s, %d: Error, number of gain elements must be %d\n",
                    __FILE__, __LINE__, nbands);
            errflg++;
        } else {
            for (i = 0; i < nbands; i++)
                ctl->gain[i] = fArray[i];
        }
    }

    if (clo_isSet(list, "offset")) {
        ctl->vic_cal_chg = 1;
        fArray = clo_getFloats(list, "offset", &ar_count);
        if (ar_count != nbands) {
            printf("%s, %d: Error, number of offset elements must be %d\n",
                    __FILE__, __LINE__, nbands);
            errflg++;
        } else {
            for (i = 0; i < nbands; i++)
                ctl->offset[i] = fArray[i];
        }
    }

    strcpy(ctl->cre_time, clo_getString(list, "cre_time"));
    if (strcmp(ctl->cre_time, "Unspecified") != 0) {
        verr = 0;
        for (i = 0; i < 8; i++)
            if (!isdigit(ctl->cre_time[i])) verr = 1;
        if (ctl->cre_time[8] != '/') verr = 1;
        for (i = 9; i < 15; i++)
            if (!isdigit(ctl->cre_time[i])) verr = 1;
        if (ctl->cre_time[15] != '.') verr = 1;
        for (i = 16; i < 22; i++)
            if (!isdigit(ctl->cre_time[i])) verr = 1;

        if (verr == 1) {
            printf("%s, %d: Error, cre_time has improper format\n",
                    __FILE__, __LINE__);
            errflg++;
        }
    }
    ctl->count_cal_opt = clo_getInt(list, "count_cal_opt");
    if ((ctl->count_cal_opt < 0) || (ctl->count_cal_opt > 2)) {
        printf("count_cal_opt must be 0, 1 or 2\n");
        errflg++;
    }

    strcpy(ctl->count_cal_gain_file, clo_getString(list,
            "count_cal_gain_file"));
    strcpy(ctl->count_cal_rvs_file, clo_getString(list,
            "count_cal_rvs_file"));
    strcpy(ctl->count_decal_gain_file, clo_getString(list,
            "count_decal_gain_file"));
    strcpy(ctl->count_decal_rvs_file, clo_getString(list,
            "count_decal_rvs_file"));

    ctl->ext_opt = clo_getInt(list, "ext_opt");
    if ((ctl->ext_opt < 0) || (ctl->ext_opt > 1)) {
        printf("ext_opt must be 0 or 1\n");
        errflg++;
    }

    strcpy(ctl->ext_coeff_file, clo_getString(list, "ext_coeff_file"));

    strcpy(ctl->id_domain, clo_getString(list, "id_domain"));
    strcpy(ctl->id_origin, clo_getString(list, "id_origin"));

    ctl->noise_mode = clo_getInt(list, "noise_opt");
    strcpy(ctl->noise_coef, clo_getString(list, "noise_coeff"));
    if (strcmp(ctl->noise_coef, "Unspecified") == 0)
        ctl->noise_mode = 0;

    ctl->stray_opt = clo_getInt(list, "stray_opt");
    strcpy(ctl->stray_tbl, clo_getString(list, "stray_tbl"));
    if (strcmp(ctl->stray_tbl, "Unspecified") == 0)
        ctl->stray_opt = 0;

    ctl->any_artifact = 0;

    if (!errflg) {
        /*
         *  make sure the minimal inputs were provided
         */
        if ((orb_geo_use == 0) || (out_loc_use == 0)) {
            printf("Both geofile and opath must be set\n");
            errflg++;
        }
        /*
         *  if electronic crosstalk is enabled, make sure the cal / decal is done
         */
        if ((ctl->ext_opt != 0) && (ctl->count_cal_opt == 0))
            ctl->count_cal_opt = 1;
        /*
         *  if any artifact addition or decal / cal is done, set general switch
         */
        if ((ctl->oxt_mode != 0) || (ctl->count_cal_opt > 0) ||
                (ctl->noise_mode != 0) || (ctl->stray_opt != 0))
            ctl->any_artifact = 1;
    }
    if (errflg) {
        clo_printUsage(list);
        exit(1);
    }
    return (0);
}
