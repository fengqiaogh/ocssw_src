#include "l1.h"

#include <ctype.h>
#include <dirent.h>
#include <stdlib.h>
#include <unistd.h>

// global L1 input pointer
l1_input_t *l1_input;

void l1_input_init() {
    int i;

    l1_input = allocateMemory(sizeof(l1_input_t), "l1_input");

    l1_input->calfile[0] = 0;
    l1_input->xcal_file[0] = 0;
    l1_input->btfile[0] = 0;

    strcpy(l1_input->pversion, "Unspecified");
    l1_input->input_parms[0] = 0;
    l1_input->input_files[0] = 0;
    l1_input->rad_opt = 1;
    l1_input->geom_per_band = 0;
    l1_input->xcal_nwave = 0;
    l1_input->xcal_opt = NULL;
    l1_input->xcal_wave = NULL;
    l1_input->resolution = -1;
    l1_input->newavhrrcal = 0;
    l1_input->sl_pixl = -1;
    l1_input->sl_frac = 0.25;
    for (i = 0; i < 10; i++) l1_input->ch22detcor[i] = 1.0; /* for modis: Ltir = Ltir / detcor */
    for (i = 0; i < 10; i++) l1_input->ch23detcor[i] = 1.0;
    for (i = 0; i < 2; i++) l1_input->cirrus_thresh[i] = -1;
    l1_input->albedo = -1.0;
    l1_input->cloud_wave = 865.0;
    l1_input->cloud_eps = -1.0;
    l1_input->glint = 0.005;
    l1_input->extreme_glint=0.03;
    l1_input->sunzen = 75.0;
    l1_input->satzen = 60.0;
    l1_input->hipol = 0.50;
    l1_input->gain = NULL;
    l1_input->offset = NULL;

    l1_input->outband_opt = 99;

    l1_input->spixl = 1;
    l1_input->epixl = -1;
    l1_input->dpixl = 1;
    l1_input->sline = 1;
    l1_input->eline = -1;
    l1_input->dline = 1;
    
    l1_input->evalmask = 0;
    l1_input->landmask = 1;
    l1_input->bathmask = 0;
    l1_input->cloudmask = 0;
    l1_input->glintmask = 0;
    l1_input->sunzenmask = 0;
    l1_input->satzenmask = 0;
    l1_input->hiltmask = 0;
    l1_input->stlightmask = 0;
}

void l1_input_delete(l1_input_t *input) {
    free(input->xcal_opt);
    input->xcal_opt = NULL;
    free(input->xcal_wave);
    input->xcal_wave = NULL;
    free(input->gain);
    input->gain = NULL;
    free(input->offset);
    input->offset = NULL;
}

void l1_add_options(clo_optionList_t* list) {
    char tmpStr[2048];
    clo_option_t *option;

    strcpy(tmpStr, "processing version string");
    clo_addOption(list, "pversion", CLO_TYPE_STRING, "Unspecified", tmpStr);

    strcpy(tmpStr, "radiation correction option (sensor-specific)\n");
    strcat(tmpStr, "        0: no correction\n");
    strcat(tmpStr, "        1: apply MERIS Smile correction");
    clo_addOption(list, "rad_opt", CLO_TYPE_INT, "0", tmpStr);

    strcpy(tmpStr, "VIIRS L1A calibration parameter file name (VIIRS only)");
    clo_addOption(list, "viirscalparfile", CLO_TYPE_IFILE, NULL, tmpStr);

    clo_addOption(list, "calfile", CLO_TYPE_IFILE, NULL, "system calibration file");

    strcpy(tmpStr, "geometry per band option:\n");
    strcat(tmpStr, "        0: use nominal viewing geometry - same for all bands\n");
    strcat(tmpStr, "        1: use band-specific viewing geometry (if available)\n");
    clo_addOption(list, "geom_per_band", CLO_TYPE_BOOL, "0", tmpStr);


    clo_addOption(list, "xcalfile", CLO_TYPE_IFILE, NULL, "cross-calibration file");

    strcpy(tmpStr, "cross-calibration option (sensor-specific) comma separated\n");
    strcat(tmpStr, "        list of option values, 1 per band, with bands listed in xcal_wave.\n");
    strcat(tmpStr, "        3: apply cross-calibration corrections (polarization and rvs)\n");
    strcat(tmpStr, "        2: apply cross-calibration polarization corrections\n");
    strcat(tmpStr, "        1: apply cross-calibration rvs corrections\n");
    strcat(tmpStr, "        0: no correction");
    clo_addOption(list, "xcal_opt", CLO_TYPE_INT, NULL, tmpStr);

    strcpy(tmpStr, "wavelengths at which to apply cross-calibration.  Comma\n");
    strcat(tmpStr, "        separated list of sensor wavelength values associated with xcal_opt.");
    clo_addOption(list, "xcal_wave", CLO_TYPE_FLOAT, NULL, tmpStr);

    clo_addOption(list, "btfile", CLO_TYPE_IFILE, NULL, "IR brightness temperature file");

    strcpy(tmpStr, "processing resolution (MODIS only)\n");
    strcat(tmpStr, "       -1: standard ocean 1km processing\n");
    strcat(tmpStr, "     1000: 1km resolution including aggregated 250 and 500m land bands\n");
    strcat(tmpStr, "      500: 500m resolution including aggregated 250 land bands and\n");
    strcat(tmpStr, "           replication for lower resolution bands\n");
    strcat(tmpStr, "      250: 250m resolution with replication for lower resolution bands");
    clo_addOption(list, "resolution", CLO_TYPE_INT, "-1", tmpStr);

    clo_addOption(list, "newavhrrcal", CLO_TYPE_INT, "0", "=1 for new noaa-16 calibration");

    strcpy(tmpStr, "\n        Channel 22 detector corrections (MODIS only)");
    clo_addOption(list, "ch22detcor", CLO_TYPE_FLOAT, "[1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]", tmpStr);

    strcpy(tmpStr, "\n        Channel 23 detector corrections (MODIS only)");
    clo_addOption(list, "ch23detcor", CLO_TYPE_FLOAT, "[1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]", tmpStr);

    clo_addOption(list, "sl_pixl", CLO_TYPE_INT, "-1", "SeaWiFS only, number of LAC pixels for\n        straylight flagging");
    clo_addOption(list, "sl_frac", CLO_TYPE_FLOAT, "0.25", "SeaWiFS only, straylight fractional\n        threshold on Ltypical");

    strcpy(tmpStr, "out-of-band correction for water-leaving\n");
    strcat(tmpStr, "        radiances\n");
    strcat(tmpStr, "        2: On (default for MODIS, SeaWiFS, OCTS)\n");
    strcat(tmpStr, "        0: Off (default for MOS, OSMI)");
    clo_addOption(list, "outband_opt", CLO_TYPE_INT, "99", tmpStr);

    strcpy(tmpStr, "evaluation bitmask\n");
    strcat(tmpStr, "          0: standard processing\n");
    strcat(tmpStr, "          1: init to old aerosol models\n");
    strcat(tmpStr, "          2: enables MODIS and MERIS cloud Mask for HABS\n");
    strcat(tmpStr, "         16: enables MODIS cirrus mask\n");
    strcat(tmpStr, "         32: use test sensor info file\n");
    strcat(tmpStr, "         64: use test rayleigh tables\n");
    strcat(tmpStr, "        128: use test aerosol tables\n");
    strcat(tmpStr, "        256: use test polarization tables\n");
    strcat(tmpStr, "       1024: mask modis mirror-side 1 (navfail)\n");
    strcat(tmpStr, "       2048: mask modis mirror-side 2 (navfail)\n");
    strcat(tmpStr, "       4096: don't apply 'cold-only' or equatorial aerosol tests for SST\n");
    strcat(tmpStr, "       8192: use alt sensor info file in eval\n");
    strcat(tmpStr, "      32768: enables spherical path geom for dtran");
    clo_addOption(list, "eval", CLO_TYPE_INT, "0", tmpStr);

    clo_addOption(list, "maskland", CLO_TYPE_BOOL, "on", "land mask option");
    clo_addOption(list, "maskbath", CLO_TYPE_BOOL, "off", "shallow water mask option");
    clo_addOption(list, "maskcloud", CLO_TYPE_BOOL, "on", "cloud mask option");
    clo_addOption(list, "maskglint", CLO_TYPE_BOOL, "off", "glint mask option");
    clo_addOption(list, "masksunzen", CLO_TYPE_BOOL, "off", "large sun zenith angle mask option");
    clo_addOption(list, "masksatzen", CLO_TYPE_BOOL, "off", "large satellite zenith angle mask option");
    clo_addOption(list, "maskhilt", CLO_TYPE_BOOL, "on", "high Lt masking");
    clo_addOption(list, "maskstlight", CLO_TYPE_BOOL, "on", "stray light masking");

    clo_addOption(list, "sunzen", CLO_TYPE_FLOAT, "75.0", "sun zenith angle threshold in deg.");
    clo_addOption(list, "satzen", CLO_TYPE_FLOAT, "60.0", "satellite zenith angle threshold");
    clo_addOption(list, "hipol", CLO_TYPE_FLOAT, "0.5", "threshold on degree-of-polarization to set\n        HIPOL flag");

    option = clo_addOption(list, "glint_thresh", CLO_TYPE_FLOAT, "0.005", "high sun glint threshold");
    clo_addOptionAlias(option, "glint");

    option = clo_addOption(list, "extreme_glint", CLO_TYPE_FLOAT, "0.03", "extreme sun glint threshold");

    option = clo_addOption(list, "cloud_thresh", CLO_TYPE_FLOAT, "0.027", "cloud reflectance\n        threshold");
    clo_addOptionAlias(option, "albedo");
    clo_addOption(list, "cloud_wave", CLO_TYPE_FLOAT, "865.0", "wavelength of cloud reflectance test");
    clo_addOption(list, "cloud_eps", CLO_TYPE_FLOAT, "-1.0", "cloud reflectance ratio threshold\n        (-1.0=disabled)");

    clo_addOption(list, "cloud_mask_file", CLO_TYPE_IFILE, NULL, "cloud mask file");

    strcpy(tmpStr, "cloud mask file variable name\n");
    strcat(tmpStr, "          0: cloud_flag\n");
    strcat(tmpStr, "          1: cloud_flag_dilated");
    clo_addOption(list, "cloud_mask_opt", CLO_TYPE_INT, "0", tmpStr);

    clo_addOption(list, "gain", CLO_TYPE_FLOAT, NULL, "calibration gain multiplier");
    clo_addOption(list, "offset", CLO_TYPE_FLOAT, NULL, "calibration offset adjustment");

    clo_addOption(list, "spixl", CLO_TYPE_INT, "1", "start pixel number");
    clo_addOption(list, "epixl", CLO_TYPE_INT, "-1", "end pixel number (-1=the last pixel)");
    clo_addOption(list, "dpixl", CLO_TYPE_INT, "1", "pixel sub-sampling interval");
    clo_addOption(list, "sline", CLO_TYPE_INT, "1", "start line number");
    clo_addOption(list, "eline", CLO_TYPE_INT, "-1", "end line number (-1=the last line)");
    clo_addOption(list, "dline", CLO_TYPE_INT, "1", "line sub-sampling interval");

}


void l1_read_default_files(clo_optionList_t *list, filehandle *l1file, const char *ifile) {
    const char *l1_defaults_prefix = "msl12";

    char *dataRoot;
    char tmpStr[FILENAME_MAX];


    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        printf("-E- OCDATAROOT environment variable is not defined.\n");
        exit(EXIT_FAILURE);
    }

    strcpy(l1file->name, ifile);
    file_format format = getFormat(l1file->name);
    if (format.type == FT_INVALID) {
        printf("-E- %s Line %d: Could not find type for file %s.\n", __FILE__, __LINE__, ifile);
        exit(EXIT_FAILURE);
    }
    l1file->format = format.type;
    l1file->sensorID = format.sensor_id;
    l1file->subsensorID = format.subsensor_id;
    l1file->nbands = rdsensorinfo(l1file->sensorID, 0, "Nbands", NULL);


    // load l1 defaults
    //sprintf(tmpStr, "%s/common/%s_defaults.par", dataRoot, l1_defaults_prefix);
    //if (want_verbose)
    //    printf("Loading default parameters from %s\n", tmpStr);
    //clo_readFile(list, tmpStr);


    // load the sensor specific defaults file
    sprintf(tmpStr, "%s/%s/%s_defaults.par", dataRoot,
            sensorId2SensorDir(l1file->sensorID), l1_defaults_prefix);
    if (want_verbose)
        printf("Loading default parameters for %s from %s\n",
            sensorId2SensorName(l1file->sensorID), tmpStr);
    clo_readFile(list, tmpStr);

    // load the sub-sensor specific defaults file
    if (l1file->subsensorID >= 0) {
        sprintf(tmpStr, "%s/%s/%s/%s_defaults.par", dataRoot,
                sensorId2SensorDir(l1file->sensorID),
                subsensorId2SubsensorDir(l1file->subsensorID), l1_defaults_prefix);
        if (access(tmpStr, R_OK) != -1) {
            if (want_verbose)
                printf("Loading default sub-sensor parameters for %s from %s\n",
                    sensorId2SensorName(l1file->sensorID), tmpStr);
            clo_readFile(list, tmpStr);
        }
    } // if sub-sensor

}


void l1_load_options(clo_optionList_t *list, filehandle *l1file) {

    clo_option_t *option;
    int count;
    int i;
    char tmp_file[FILENAME_MAX];
    char *strVal;
    int *iArray;
    float *fArray;

    l1_input->xcal_opt = calloc_nbandsi32t(l1file->nbands, l1_input->xcal_opt, 0);
    l1_input->xcal_wave = calloc_nbandsf(l1file->nbands, l1_input->xcal_wave, -1.0);
    l1_input->gain = calloc_nbandsf(l1file->nbands, l1_input->gain, 1.0);
    l1_input->offset = calloc_nbandsf(l1file->nbands, l1_input->offset, 0.0);

    strVal = clo_getString(list, "pversion");
    strcpy(l1_input->pversion, strVal);

    option = clo_findOption(list, "viirscalparfile");
    if (clo_isOptionSet(option)) {
        strVal = clo_getOptionString(option);
        parse_file_name(strVal, tmp_file);
        strcpy(l1_input->viirscalparfile, tmp_file);
    }

    l1_input->rad_opt = clo_getInt(list, "rad_opt");

    strVal = clo_getString(list, "calfile");
    if(strVal && strVal[0]) {
        parse_file_name(strVal, tmp_file);
        strcpy(l1_input->calfile, tmp_file);
    }

    l1_input->geom_per_band = clo_getBool(list, "geom_per_band");

    option = clo_findOption(list, "xcalfile");
    if (clo_isOptionSet(option)) {
        strVal = clo_getOptionString(option);
        parse_file_name(strVal, tmp_file);
        strcpy(l1_input->xcal_file, tmp_file);
    } else {
        // Look for default xcalfile if not previously provided
        char *varRoot;
        char xcaldir[FILENAME_MAX];
        if ((varRoot = getenv("OCVARROOT")) == NULL) {
            printf("-E- %s, %d: OCVARROOT env variable undefined.\n", __FILE__,
                    __LINE__);
            exit(EXIT_FAILURE);
        }
        strcpy(xcaldir, varRoot);
        char *lcsensor = strdup(sensorId2SensorName(l1file->sensorID));
        for (i = 0; lcsensor[i]; i++) {
            lcsensor[i] = tolower(lcsensor[i]);
        }
        strcat(xcaldir, "/");
        strcat(xcaldir, lcsensor);
        strcat(xcaldir, "/xcal/OPER");
        free(lcsensor);

        DIR *dir;
        struct dirent *ent;
        if ((dir = opendir(xcaldir)) != NULL) {
            /* print all the files and directories within directory */
            while ((ent = readdir(dir)) != NULL) {
                if (strncmp(ent->d_name, "xcal_", 5) == 0)
                    break;
            }
            if(ent) {
                char *xcalfile_prefix = strdup(ent->d_name);
                closedir(dir);
                char *tmpptr = strrchr(xcalfile_prefix, '_');
                *tmpptr = 0;
                strcat(xcaldir, "/");
                strcat(xcaldir, xcalfile_prefix);
                free(xcalfile_prefix);
                strcpy(l1_input->xcal_file, xcaldir);
            }
        }
    }

    option = clo_findOption(list, "xcal_opt");
    if (clo_isOptionSet(option)) {
        iArray = clo_getOptionInts(option, &count);
        l1_input->xcal_nwave = count;
        if (count > l1file->nbands) {
            printf("-E- number of xcal_opt elements (%d) must be %d or less\n", count, l1file->nbands);
            exit(1);
        }
        for (i = 0; i < count; i++)
            l1_input->xcal_opt[i] = iArray[i];
    }

    option = clo_findOption(list, "xcal_wave");
    if (clo_isOptionSet(option)) {
        fArray = clo_getOptionFloats(option, &count);
        if (count > l1file->nbands) {
            printf("-E- number of xcal_wave elements (%d) must be %d or less\n", count, l1file->nbands);
            exit(1);
        }
        if (count != l1_input->xcal_nwave) {
            printf("-W- Number of xcal_wave elements (%d) should be equal to xcal_opt number elements (%d)\n", count, l1_input->xcal_nwave);
        }
        l1_input->xcal_nwave = count;
        for (i = 0; i < count; i++)
            l1_input->xcal_wave[i] = fArray[i];
    }

    option = clo_findOption(list, "btfile");
    if (clo_isOptionSet(option)) {
        strVal = clo_getOptionString(option);
        parse_file_name(strVal, tmp_file);
        strcpy(l1_input->btfile, tmp_file);
    }

    l1_input->resolution = clo_getInt(list, "resolution");
    if (l1_input->resolution == -1000)
        l1_input->resolution = 1000;


    l1_input->newavhrrcal = clo_getInt(list, "newavhrrcal");

    fArray = clo_getFloats(list, "ch22detcor", &count);
    if (count != 10) {
        printf("-E- number of ch22detcor elements must be 10 \n");
        exit(1);
    }
    for (i = 0; i < count; i++)
        l1_input->ch22detcor[i] = fArray[i];

    fArray = clo_getFloats(list, "ch23detcor", &count);
    if (count != 10) {
        printf("-E- number of ch23detcor elements must be 10 \n");
        exit(1);
    }
    for (i = 0; i < count; i++)
        l1_input->ch23detcor[i] = fArray[i];

    l1_input->sl_pixl = clo_getInt(list, "sl_pixl");
    l1_input->sl_frac = clo_getFloat(list, "sl_frac");
    l1_input->outband_opt = clo_getInt(list, "outband_opt");
    l1_input->evalmask = clo_getInt(list, "eval");
    l1_input->landmask = clo_getBool(list, "maskland");
    l1_input->bathmask = clo_getBool(list, "maskbath");
    l1_input->cloudmask = clo_getBool(list, "maskcloud");
    l1_input->glintmask = clo_getBool(list, "maskglint");
    l1_input->sunzenmask = clo_getBool(list, "masksunzen");
    l1_input->satzenmask = clo_getBool(list, "masksatzen");
    l1_input->hiltmask = clo_getBool(list, "maskhilt");
    l1_input->stlightmask = clo_getBool(list, "maskstlight");
    l1_input->sunzen = clo_getFloat(list, "sunzen");
    l1_input->satzen = clo_getFloat(list, "satzen");
    l1_input->hipol = clo_getFloat(list, "hipol");
    l1_input->glint = clo_getFloat(list, "glint_thresh");
    l1_input->extreme_glint = clo_getFloat(list, "extreme_glint");
    l1_input->albedo = clo_getFloat(list, "cloud_thresh");
    l1_input->cloud_wave = clo_getFloat(list, "cloud_wave");
    l1_input->cloud_eps = clo_getFloat(list, "cloud_eps");

    option = clo_findOption(list, "cloud_mask_file");
    if (clo_isOptionSet(option)) {
        strVal = clo_getOptionString(option);
        parse_file_name(strVal, tmp_file);
        strcpy(l1_input->cld_msk_file, tmp_file);
    }
    l1_input->cloud_mask_opt = clo_getInt(list, "cloud_mask_opt");

    option = clo_findOption(list, "gain");
    if (clo_isOptionSet(option)) {
        fArray = clo_getOptionFloats(option, &count);
        if (count != l1file->nbands) {
            printf("-E- number of gain elements (%d) must be equal to number of bands (%d)\n", count, l1file->nbands);
            exit(1);
        }
        for (i = 0; i < count; i++)
            l1_input->gain[i] = fArray[i];
    }

    option = clo_findOption(list, "offset");
    if (clo_isOptionSet(option)) {
        fArray = clo_getOptionFloats(option, &count);
        if (count != l1file->nbands) {
            printf("-E- number of offset elements (%d) must be equal to nu,ber of bands (%d)\n", count, l1file->nbands);
            exit(1);
        }
        for (i = 0; i < count; i++)
            l1_input->offset[i] = fArray[i];
    }

    l1_input->spixl = clo_getInt(list, "spixl");
    l1_input->epixl = clo_getInt(list, "epixl");
    l1_input->dpixl = clo_getInt(list, "dpixl");
    l1_input->sline = clo_getInt(list, "sline");
    l1_input->eline = clo_getInt(list, "eline");
    l1_input->dline = clo_getInt(list, "dline");

}



void l1_get_input_params(filehandle *l1file, char *input_parms) {
    int i;
    char str_buf[FILENAME_MAX];

    sprintf(str_buf, "pversion = %s", l1_input->pversion);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "viirscalparfile = %s", l1_input->viirscalparfile);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "rad_opt = %3d", l1_input->rad_opt);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "calfile = %s", l1_input->calfile);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "geom_per_band = %3d", l1_input->geom_per_band);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "xcalfile = %s", l1_input->xcal_file);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "xcal_opt = %3d", l1_input->xcal_opt[0]);
    strcat(input_parms, str_buf);
    for (i = 1; i < l1_input->xcal_nwave; i++) {
        sprintf(str_buf, ", %3d", l1_input->xcal_opt[i]);
        strcat(input_parms, str_buf);
    }

    strcat(input_parms, "\n");
    sprintf(str_buf, "xcal_wave = %8.4f", l1_input->xcal_wave[0]);
    strcat(input_parms, str_buf);
    for (i = 1; i < l1_input->xcal_nwave; i++) {
        sprintf(str_buf, ", %8.4f", l1_input->xcal_wave[i]);
        strcat(input_parms, str_buf);
    }
    strcat(input_parms, "\n");

    sprintf(str_buf, "btfile = %s", l1_input->btfile);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "resolution = %3d", l1_input->resolution);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "newavhrrcal = %d", l1_input->newavhrrcal);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "ch22detcor = %9.6f", l1_input->ch22detcor[0]);
    strcat(input_parms, str_buf);
    for (i = 1; i < 10; i++) {
        sprintf(str_buf, ", %9.6f", l1_input->ch22detcor[i]);
        strcat(input_parms, str_buf);
    }
    strcat(input_parms, "\n");

    sprintf(str_buf, "ch23detcor = %9.6f", l1_input->ch23detcor[0]);
    strcat(input_parms, str_buf);
    for (i = 1; i < 10; i++) {
        sprintf(str_buf, ", %9.6f", l1_input->ch23detcor[i]);
        strcat(input_parms, str_buf);
    }
    strcat(input_parms, "\n");

    sprintf(str_buf, "sl_pixl = %3d", l1_input->sl_pixl);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "sl_frac = %8.4f", l1_input->sl_frac);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "outband_opt = %3d", l1_input->outband_opt);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "eval = %5d", l1_input->evalmask);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "maskland = %2d", l1_input->landmask);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "maskbath = %2d", l1_input->bathmask);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "maskcloud = %2d", l1_input->cloudmask);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "maskglint = %2d", l1_input->glintmask);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "masksunzen = %2d", l1_input->sunzenmask);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "masksatzen = %2d", l1_input->satzenmask);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "maskhilt = %2d", l1_input->hiltmask);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "maskstlight = %2d", l1_input->stlightmask);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "sunzen = %8.3f", l1_input->sunzen);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "satzen = %8.3f", l1_input->satzen);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "hipol = %8.3f", l1_input->hipol);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "glint_thresh = %8.3f", l1_input->glint);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "extreme_glint = %8.3f", l1_input->extreme_glint);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "cloud_thresh = %8.3f", l1_input->albedo);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "cloud_wave = %8.3f", l1_input->cloud_wave);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "cloud_eps = %8.3f", l1_input->cloud_eps);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "cloud_mask_file= %s", l1_input->cld_msk_file);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "gain = %8.4f", l1_input->gain[0]);
    strcat(input_parms, str_buf);
    for (i = 1; i < l1file->nbands; i++) {
        sprintf(str_buf, ", %8.4f", l1_input->gain[i]);
        strcat(input_parms, str_buf);
    }
    strcat(input_parms, "\n");

    sprintf(str_buf, "offset = %8.5f", l1_input->offset[0]);
    strcat(input_parms, str_buf);
    for (i = 1; i < l1file->nbands; i++) {
        sprintf(str_buf, ", %8.5f", l1_input->offset[i]);
        strcat(input_parms, str_buf);
    }
    strcat(input_parms, "\n");

    sprintf(str_buf, "spixl = %5d", l1_input->spixl);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "epixl = %5d", l1_input->epixl);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "dpixl = %5d", l1_input->dpixl);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "sline = %5d", l1_input->sline);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "eline = %5d", l1_input->eline);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");

    sprintf(str_buf, "dline = %5d", l1_input->dline);
    strcat(input_parms, str_buf);
    strcat(input_parms, "\n");





}



void l1_get_input_files(filehandle *l1file, char *input_files) {
    char *tmp_str;
    char str_buf[FILENAME_MAX];

    if (l1_input->viirscalparfile[0]) {
        tmp_str = strrchr(l1_input->viirscalparfile, '/');
        tmp_str = (tmp_str == 0x0) ? l1_input->viirscalparfile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input_files, str_buf);
    }

    if (l1_input->calfile[0]) {
        tmp_str = strrchr(l1_input->calfile, '/');
        tmp_str = (tmp_str == 0x0) ? l1_input->calfile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input_files, str_buf);
    }

    if (l1_input->xcal_file[0]) {
        tmp_str = strrchr(l1_input->xcal_file, '/');
        tmp_str = (tmp_str == 0x0) ? l1_input->xcal_file : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input_files, str_buf);
    }

    if (l1_input->btfile[0]) {
        tmp_str = strrchr(l1_input->btfile, '/');
        tmp_str = (tmp_str == 0x0) ? l1_input->btfile : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input_files, str_buf);
    }

    if (l1_input->cld_msk_file[0]) {
        tmp_str = strrchr(l1_input->cld_msk_file, '/');
        tmp_str = (tmp_str == 0x0) ? l1_input->cld_msk_file : tmp_str + 1;
        sprintf(str_buf, ",%s", tmp_str);
        strcat(input_files, str_buf);
    }
}
