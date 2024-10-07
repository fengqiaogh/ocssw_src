/* ============================================================================ */
/* module l1_viirs_nc.c - functions to read VIIRS L1A for MSL12                 */
/* Written By:  Joel M. Gales GSFC Futuretech, Sep. 2011.                       */
/*                                                                              */
/* ============================================================================ */

// #include <stdbool.h>
#include <nc4utils.h>
#include "libnav.h"
#include "calibrate_viirs.h"
#include <Calibrate_Viirs_Connector.h>
#include <productInfo.h>
#include <libnav.h>
#include "l1.h"
#include "l2_flags.h"

#define MBAND_NUM_DETECTORS 16

static int32_t prevScan = -1;

static int scanLineGrp;
static int scanStartTimeId;
static int HAMSideId;

static int geoFileId;
static int geoNavigationGrp;
static int geoGeolocationGrp;
static int geoScanLineGrp;
static int lonId, latId, senzId, senaId, solzId, solaId, angId, posId, velId, scanQualityId, pixelQualityId;

static double ***f_cal_corr = NULL; /* f table correction [band][det][ms] */

static float *Fobar; // reflectance to radiance conversion factors
static int extract_pixel_start = 0;
static int extract_pixel_stop = 0;

static short *tmpShort;
static unsigned char *tmpByte;
static int nscan, nline, npix;

static double starttime;
static double lastvalidtime;
static int lastvalidscan = 0;
static double time_interval;

static float latGeoFillValue = -999.9;
static float lonGeoFillValue = -999.9;
static short senzGeoFillValue = -32768;
static short senaGeoFillValue = -32768;
static short solzGeoFillValue = -32768;
static short solaGeoFillValue = -32768;

static float latL2FillValue = -999.0;
static float lonL2FillValue = -999.0;
static float senzL2FillValue = -32767;
static float senaL2FillValue = -32767;
static float solzL2FillValue = -32767;
static float solaL2FillValue = -32767;

static double scan_start_tai = -999;

int openl1_viirs_l1a(filehandle * file) {
    char *fltime;

    size_t tmpSize;
    int fileID, dimid, status;
    size_t att_len; // Change to size_t JMG 08/12/13
    int orbit_number;

    // Open the netcdf4 input file
    status = nc_open(file->name, NC_NOWRITE, &fileID);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
                __FILE__, __LINE__, file->name);
        return (1);
    }

    // Get pixel and scan dimensions
    status = nc_get_att_int(fileID, NC_GLOBAL, "number_of_filled_scans", &nscan);
    check_err(status, __LINE__, __FILE__);
    nline = nscan * MBAND_NUM_DETECTORS;

    status = nc_inq_dimid(fileID, "Mband_pixels", &dimid);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_dimlen(fileID, dimid, &tmpSize);
    check_err(status, __LINE__, __FILE__);
    npix = tmpSize;

    nc_type vr_type; /* attribute type */
    size_t vr_len; /* attribute length */
    if ((nc_inq_att(fileID, NC_GLOBAL, "extract_pixel_start", &vr_type, &vr_len) == NC_NOERR)) {
        status = nc_get_att_int(fileID, NC_GLOBAL, "extract_pixel_start", &extract_pixel_start);
        check_err(status, __LINE__, __FILE__);
        extract_pixel_start--; // Attribute is one-based
        status = nc_get_att_int(fileID, NC_GLOBAL, "extract_pixel_stop", &extract_pixel_stop);
        check_err(status, __LINE__, __FILE__);
        extract_pixel_stop--; // Attribute is one-based
        if (npix != (extract_pixel_stop - extract_pixel_start + 1)) {
            printf("-E- Problem with the extracted L1A file pixel dimension.\n");
            printf("    npix(%d), extract_pixel_stop(%d), extract_pixel_start(%d) do not work together.\n",
                    npix, extract_pixel_stop, extract_pixel_start);
            exit(EXIT_FAILURE);
        }
    }

    if (want_verbose) {
        printf("VIIRS L1A Npix  :%d Nlines:%d\n", npix, nline);
    } // want_verbose

    // get start and end time
    status = nc_inq_attlen(fileID, NC_GLOBAL, "time_coverage_start", &att_len);
    check_err(status, __LINE__, __FILE__);

    /* allocate required space before retrieving values */
    fltime = (char *) malloc(att_len + 1); /* + 1 for trailing null */

    /* get attribute values */
    status = nc_get_att_text(fileID, NC_GLOBAL, "time_coverage_start", fltime);
    check_err(status, __LINE__, __FILE__);
    fltime[att_len] = '\0';

    starttime = lastvalidtime = isodate2unix(fltime);
    lastvalidscan = 0;
    free(fltime);

    status = nc_inq_attlen(fileID, NC_GLOBAL, "time_coverage_end", &att_len);
    check_err(status, __LINE__, __FILE__);

    /* allocate required space before retrieving values */
    fltime = (char *) malloc(att_len + 1); /* + 1 for trailing null */

    /* get attribute values */
    status = nc_get_att_text(fileID, NC_GLOBAL, "time_coverage_end", fltime);
    check_err(status, __LINE__, __FILE__);
    fltime[att_len] = '\0';

    double stoptime = isodate2unix(fltime);
    free(fltime);

    time_interval = (stoptime - starttime) / (nscan - 1); /* in sec */

    if ((nc_inq_att(fileID, NC_GLOBAL, "orbit_number", &vr_type, &vr_len) == NC_NOERR)) {
        status = nc_get_att_int(fileID, NC_GLOBAL, "orbit_number", &orbit_number);
        check_err(status, __LINE__, __FILE__);
    } else {
        status = nc_get_att_int(fileID, NC_GLOBAL, "OrbitNumber", &orbit_number);
        check_err(status, __LINE__, __FILE__);
    }

    file->sd_id = fileID;
    file->nbands = 10;
    file->npix = npix;
    file->nscan = nline;
    file->ndets = MBAND_NUM_DETECTORS;
    file->terrain_corrected = 1; // presumed.
    file->orbit_number = orbit_number;
    strcpy(file->spatialResolution, "750 m");

    rdsensorinfo(file->sensorID, l1_input->evalmask,
            "Fobar", (void **) &Fobar);

    if (want_verbose)
        printf("file->nbands = %d\n", (int) file->nbands);

    status = nc_inq_ncid(file->sd_id, "scan_line_attributes", &scanLineGrp);
    check_err(status, __LINE__, __FILE__);

    // get tai93
    status = nc_inq_varid(scanLineGrp, "scan_start_time", &scanStartTimeId);
    check_err(status, __LINE__, __FILE__);

    // get mirror side
    status = nc_inq_varid(scanLineGrp, "HAM_side", &HAMSideId);
    check_err(status, __LINE__, __FILE__);

    // Setup geofile pointers
    if (file->geofile && file->geofile[0]) {

        status = nc_open(file->geofile, NC_NOWRITE, &geoFileId);
        if (status != NC_NOERR) {
            printf("-E- Could not open GEO file \"%s\"\n", file->geofile);
            exit(EXIT_FAILURE);
        }

        status = nc_inq_grp_ncid(geoFileId, "geolocation_data", &geoGeolocationGrp);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "longitude", &lonId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_var_fill(geoGeolocationGrp, lonId, NULL, &lonGeoFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "latitude", &latId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_var_fill(geoGeolocationGrp, latId, NULL, &latGeoFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "sensor_zenith", &senzId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_var_fill(geoGeolocationGrp, senzId, NULL, &senzGeoFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "sensor_azimuth", &senaId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_var_fill(geoGeolocationGrp, senaId, NULL, &senaGeoFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "solar_zenith", &solzId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_var_fill(geoGeolocationGrp, solzId, NULL, &solzGeoFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "solar_azimuth", &solaId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_var_fill(geoGeolocationGrp, solaId, NULL, &solaGeoFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "quality_flag", &pixelQualityId);
        check_err(status, __LINE__, __FILE__);

        status = nc_inq_grp_ncid(geoFileId, "navigation_data", &geoNavigationGrp);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoNavigationGrp, "att_ang_mid", &angId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoNavigationGrp, "orb_pos_ev_mid", &posId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoNavigationGrp, "orb_vel_ev_mid", &velId);
        check_err(status, __LINE__, __FILE__);

        status = nc_inq_grp_ncid(geoFileId, "scan_line_attributes", &geoScanLineGrp);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoScanLineGrp, "scan_quality", &scanQualityId);
        check_err(status, __LINE__, __FILE__);
    } // geofile

    // Setup the fill values for the geo products
    productInfo_t* info = allocateProductInfo();
    status = findProductInfo("lat", VIIRSN, info);
    if (status)
        latL2FillValue = info->fillValue;
    status = findProductInfo("lon", VIIRSN, info);
    if (status)
        lonL2FillValue = info->fillValue;
    status = findProductInfo("sena", VIIRSN, info);
    if (status)
        senaL2FillValue = info->fillValue;
    status = findProductInfo("senz", VIIRSN, info);
    if (status)
        senzL2FillValue = info->fillValue;
    status = findProductInfo("sola", VIIRSN, info);
    if (status)
        solaL2FillValue = info->fillValue;
    status = findProductInfo("solz", VIIRSN, info);
    if (status)
        solzL2FillValue = info->fillValue;
    freeProductInfo(info);

    return (LIFE_IS_GOOD);
}

int readl1_viirs_l1a(filehandle *file, int32_t line, l1str *l1rec) {
    static int firstCall = 1;

    int32_t ip, ib, ipb;
    int i;
    double f_corr;

    int status; //, ncid, varid;
    size_t start[] = { 0, 0, 0 };
    size_t count[] = { 1, 1, 1 };

    int32_t scan = line / MBAND_NUM_DETECTORS;

    for (ip = 0; ip < npix; ip++) {
        l1rec->pixnum[ip] = ip + extract_pixel_start;
    }

    if (firstCall) {
        firstCall = 0;

        // init calibration only if a geo file was given
        if (file->geofile && file->geofile[0]) {

            /*----- Calibration LUT -----*/
            if (l1_input->calfile[0]) {
                double grantime = starttime + 378691200.0;  // unix (1970) -> UTC58
                grantime *= 1000000.0; // convert to IET
                load_fcal_lut(l1_input->calfile, (int64_t) grantime, &f_cal_corr);
            }

            // Initialize L1A calibration
            VcstViirsCal_initialize(file->name, l1_input->viirscalparfile);

            tmpShort = (short *) malloc(npix * sizeof(short));
            tmpByte = (unsigned char *) malloc(npix);
        }
    }

    //    l1rec->sensorID = file->sensorID;
    l1rec->npix = file->npix;

    //      set time for this scan - if scan_start_time value not properly set, estimate from scene start time.
    if (scan != prevScan) {
        start[0] = scan;
        status = nc_get_var1_double(scanLineGrp, scanStartTimeId, start, &scan_start_tai);
        check_err(status, __LINE__, __FILE__);
    }

    if (scan_start_tai > 0) {
        lastvalidtime = tai93_to_unix(scan_start_tai);
        lastvalidscan = line;
        l1rec->scantime = lastvalidtime;
    } else {
        l1rec->scantime = lastvalidtime + (time_interval * (line - lastvalidscan));
    }

    //------------------------------------
    // if there is no geo file just return
    // This is used for l1info with only a L1A file and no GEO
    //-------------------------------------
    if (!file->geofile || !file->geofile[0]) {
        return 0;
    }

    // first check the scan quality flag
    // 1   SCE_side_A_B
    // 2   SCE_side_invalid
    // 4   Sector_rotation
    // 8   Encoder_degraded
    // 16  SAA
    // 32  Solar_eclipse
    // 64  Lunar_eclipse
    // 128 HAM_side
    //
    // Sector_rotation
    short scanQualityWarnMask = 2 | 8 | 128;
    short scanQualityFailMask = 4;
    static short scanQualityFlag = 0;

    if (scan != prevScan) {
        start[0] = scan;
        status = nc_get_var1_short(geoScanLineGrp, scanQualityId, start, &scanQualityFlag);
        check_err(status, __LINE__, __FILE__);
    }
    if (scanQualityFlag & scanQualityFailMask) {
        for (ip = 0; ip < npix; ip++)
            l1rec->flags[ip] |= NAVFAIL;
        return 0;
    }
    if (scanQualityFlag & scanQualityWarnMask) {
        for (ip = 0; ip < npix; ip++)
            l1rec->flags[ip] |= NAVWARN;
    }

    static unsigned char HAMSideVal = 0;
    if (scan != prevScan) {
        start[0] = scan;
        status = nc_get_var1_uchar(scanLineGrp, HAMSideId, start, &HAMSideVal);
        check_err(status, __LINE__, __FILE__);
    }
    l1rec->mside = HAMSideVal;

    // set up to read all pixels of the line.
    start[0] = line;
    start[1] = 0;
    count[0] = 1;
    count[1] = npix; // read all pixels

    status = nc_get_vara_float(geoGeolocationGrp, latId, start, count, l1rec->lat);
    check_err(status, __LINE__, __FILE__);
    for (i = 0; i < npix; i++)
        if (l1rec->lat[i] == latGeoFillValue)
            l1rec->lat[i] = latL2FillValue;

    status = nc_get_vara_float(geoGeolocationGrp, lonId, start, count, l1rec->lon);
    check_err(status, __LINE__, __FILE__);
    for (i = 0; i < npix; i++)
        if (l1rec->lon[i] == lonGeoFillValue)
            l1rec->lon[i] = lonL2FillValue;

    status = nc_get_vara_short(geoGeolocationGrp, solzId, start, count, tmpShort);
    check_err(status, __LINE__, __FILE__);
    for (i = 0; i < npix; i++)
        if (tmpShort[i] == solzGeoFillValue)
            l1rec->solz[i] = solzL2FillValue;
        else
            l1rec->solz[i] = tmpShort[i] * 0.01;

    status = nc_get_vara_short(geoGeolocationGrp, solaId, start, count, tmpShort);
    check_err(status, __LINE__, __FILE__);
    for (i = 0; i < npix; i++)
        if (tmpShort[i] == solaGeoFillValue)
            l1rec->sola[i] = solaL2FillValue;
        else
            l1rec->sola[i] = tmpShort[i] * 0.01;

    status = nc_get_vara_short(geoGeolocationGrp, senzId, start, count, tmpShort);
    check_err(status, __LINE__, __FILE__);
    for (i = 0; i < npix; i++)
        if (tmpShort[i] == senzGeoFillValue)
            l1rec->senz[i] = senzL2FillValue;
        else
            l1rec->senz[i] = tmpShort[i] * 0.01;

    status = nc_get_vara_short(geoGeolocationGrp, senaId, start, count, tmpShort);
    check_err(status, __LINE__, __FILE__);
    for (i = 0; i < npix; i++)
        if (tmpShort[i] == senaGeoFillValue)
            l1rec->sena[i] = senaL2FillValue;
        else
            l1rec->sena[i] = tmpShort[i] * 0.01;

    /* Load Angles */
    float ang[3]; // degrees
    float pos[3]; // km
    float vel[3]; // km/sec
    size_t s3[] = { scan, 0 };
    size_t c3[] = { 1, 3 };
    status = nc_get_vara_float(geoNavigationGrp, angId, s3, c3, ang);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(geoNavigationGrp, posId, s3, c3, pos);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(geoNavigationGrp, velId, s3, c3, vel);
    check_err(status, __LINE__, __FILE__);
    for (i = 0; i < 3; i++) {
        pos[i] /= 1000.; // m   -> km
        vel[i] /= 1000.; // m/s -> km/s
    }

    /* Compute polarization rotation angles */
    float sen_mat[3][3], coeff[10];
    double mnorm[3];
    ocorient_(pos, vel, ang, sen_mat, coeff);
    for (i = 0; i < 3; i++)
        mnorm[i] = sen_mat[i][0];
    compute_alpha(l1rec->lon, l1rec->lat,
            l1rec->senz, l1rec->sena,
            mnorm, l1rec->npix, l1rec->alpha);

    /* Check pixel values */
    status = nc_get_vara_uchar(geoGeolocationGrp, pixelQualityId, start, count, tmpByte);
    check_err(status, __LINE__, __FILE__);
    // 1 Input_invalid
    // 2 Pointing_bad
    // 4 Terrain_bad
    unsigned char qualityFailMask = 1 | 2;
    unsigned char qualityWarnMask = 4;
    for (i = 0; i < npix; i++) {
        if (tmpByte[i] & qualityFailMask)
            l1rec->flags[i] |= NAVFAIL;
        if (tmpByte[i] & qualityWarnMask)
            l1rec->flags[i] |= NAVWARN;
    }

    /* Earth-sun distance correction for this scan */
    static double esdist = -999.9;
    if (scan != prevScan) {
	int16_t year, day;
	double dsec;
	unix2yds(l1rec->scantime, &year, &day, &dsec);
	int32_t yr = (int32_t) year;
	int32_t dy = (int32_t) day;
	int32_t msec = (int32_t) (dsec * 1000.0);
        esdist = esdist_(&yr, &dy, &msec);
    }
    l1rec->fsol = pow(1.0 / esdist, 2);

    int nbands = 16; //, nRSBbands = 10, nCIRbands = 1, nTEBbands = 5;

    // read in calibrated L1B data
    static float *l1bptrs[16] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
        0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};

    static char *bandType[16] = {"RSB", "RSB", "RSB", "RSB", "RSB", "RSB",
        "RSB", "RSB", "CIR", "RSB", "RSB", "TEB",
        "TEB", "TEB", "TEB", "TEB"};

    // Note: l1bptrs arrays are 3200 pixels wide
    int oldVerbose = want_verbose;
    want_verbose = 0;
    VcstViirsCal_calibrateMOD(line, nbands, l1bptrs);
    want_verbose = oldVerbose;

    l1rec->detnum = line % file->ndets;

    int irsb = 0, iteb = 0;
    int l1bptrs_scan_idx;
    for (ib = 0; ib < nbands; ib++) {

        /* get specific f table cal correction  */
        f_corr = (f_cal_corr == NULL) ? 1.0
                : f_cal_corr[ib][l1rec->detnum][l1rec->mside];

        if (strcmp(bandType[ib], "TEB") == 0) {

            for (ip = 0; ip < npix; ip++) {
                ipb = ip * NBANDSIR + iteb;
                l1rec->Ltir[ipb] = 0;

                l1bptrs_scan_idx = l1rec->detnum * 3200 + ip + extract_pixel_start;

                if (l1bptrs[ib][l1bptrs_scan_idx] != -32767) {
                    l1rec->Ltir[ipb] = l1bptrs[ib][l1bptrs_scan_idx] / 10.0;

                    /* Apply F-factor */
                    l1rec->Ltir[ipb] *= f_corr;
                }

            }
            iteb++;

        } else if (strcmp(bandType[ib], "CIR") == 0) {

            for (ip = 0; ip < npix; ip++) {

                l1bptrs_scan_idx = l1rec->detnum * 3200 + ip + extract_pixel_start;

                if (l1bptrs[ib][l1bptrs_scan_idx] != -32767) {
                    l1rec->rho_cirrus[ip] = l1bptrs[ib][l1bptrs_scan_idx];

                    /* Normalize reflectance by solar zenith angle */
                    l1rec->rho_cirrus[ip] /= cos(l1rec->solz[ip] / RADEG);

                    /* Apply F-factor */
                    l1rec->rho_cirrus[ip] *= f_corr;
                }

            }

        } else if (strcmp(bandType[ib], "RSB") == 0) {

            l1rec->Fo[irsb] = Fobar[irsb] * l1rec->fsol;

            // copy to Lt record.
            for (ip = 0; ip < npix; ip++) {
                ipb = ip * l1rec->l1file->nbands + irsb;

                l1bptrs_scan_idx = l1rec->detnum * 3200 + ip + extract_pixel_start;

                if (l1bptrs[ib][l1bptrs_scan_idx] != -32767) {
                    l1rec->Lt[ipb] = l1bptrs[ib][l1bptrs_scan_idx];

                    /* convert from reflectance to radiance */
                    l1rec->Lt[ipb] *= l1rec->Fo[irsb] / PI;

                    /* Apply F-factor */
                    l1rec->Lt[ipb] *= f_corr;
                }

            }
            irsb++;
        } // if RSB

    } // for ib

    radiance2bt(l1rec, -1); // calculate brightness temperature

    for (ip = 0; ip < npix; ip++) {
        flag_bowtie_deleted(l1rec, ip, extract_pixel_start);
    }

    prevScan = scan;
    return (LIFE_IS_GOOD);
}

int readl1_lonlat_viirs_nc(filehandle *file, int32_t line, l1str *l1rec) {
    int32_t ip;
    int status;
    size_t start[] = { 0, 0, 0 };
    size_t count[] = { 1, 1, 1 };

    int32_t scan = line / MBAND_NUM_DETECTORS;

    if (!file->geofile || !file->geofile[0]) {
        printf("-E- Geolocation file needs to be set\n");
        exit(EXIT_FAILURE);
    }

    //    l1rec->sensorID = file->sensorID;
    l1rec->npix = file->npix;

    //      set time for this scan - if scan_start_time value not properly set, estimate from scene start time.
    if (scan != prevScan) {
        start[0] = scan;
        status = nc_get_var1_double(scanLineGrp, scanStartTimeId, start, &scan_start_tai);
        check_err(status, __LINE__, __FILE__);
    }

    if (scan_start_tai > 0) {
        lastvalidtime = tai93_to_unix(scan_start_tai);
        lastvalidscan = line;
        l1rec->scantime = lastvalidtime;
    } else {
        l1rec->scantime = lastvalidtime + (time_interval * (line - lastvalidscan));
    }

    // set up to read all pixels of the line.
    start[0] = line;
    start[1] = 0;
    count[0] = 1;
    count[1] = npix; // read all pixels

    status = nc_get_vara_float(geoGeolocationGrp, latId, start, count, l1rec->lat);
    check_err(status, __LINE__, __FILE__);

    status = nc_get_vara_float(geoGeolocationGrp, lonId, start, count, l1rec->lon);
    check_err(status, __LINE__, __FILE__);

    prevScan = scan;

    for (ip = 0; ip < npix; ip++) {
        // convert the fill values to the proper value
        if (l1rec->lat[ip] == latGeoFillValue)
            l1rec->lat[ip] = latL2FillValue;
        if (l1rec->lon[ip] == lonGeoFillValue)
            l1rec->lon[ip] = lonL2FillValue;
    }

    return (LIFE_IS_GOOD);
}

int closel1_viirs_l1a(filehandle *file) {
    int status;

    status = nc_close(file->sd_id);
    check_err(status, __LINE__, __FILE__);

    if (file->geofile && file->geofile[0]) {
        status = nc_close(geoFileId);
        check_err(status, __LINE__, __FILE__);

        free(tmpShort);
        free(tmpByte);
    }

    return (LIFE_IS_GOOD);
}

