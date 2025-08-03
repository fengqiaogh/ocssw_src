/**
 * @file l1b_viirs_nc.c
 * @brief Reads VIIRS Level 1B and geolocation data from NetCDF4 files
 * @author Gwyn Fireman
 * @author Ocean Ecology Laboratory, NASA Goddard Space Flight Center
 */
#include "l1_viirs_l1b.h"
#include "l1.h"

#include "nc4utils.h"
#include <libnav.h>
#include <float.h>
#include <math.h>

typedef unsigned short ushort;

/***********************************************************************
                             Calibration
 ***********************************************************************/
static double ***f_cal_corr = NULL; /* f table correction [band][det][ms] */

/***********************************************************************
                            Initialization
 ***********************************************************************/
static float *Fobar; // reflectance to radiance conversion factors
static int extract_pixel_start = 0;
static int extract_pixel_stop = 0;

/* Variables to load for each input file and group */
typedef struct {
    size_t index;
    char *name;
} varlist;

/*----- GEO file -----*/

/* Scan Line attributes, to be read on initialization.
   All values are dimensioned as [number_of_scans]. */
static const char* SCN_GRP = "scan_line_attributes";

enum scn_var {
    SCN_STIME, /*< double scan_start_time */
    //SCN_ETIME, /*< double scan_end_time   */
    //SCN_MTIME, /*< double ev_mid_time     */
    SCN_MSIDE, /*< ubyte  HAM_side        */
    SCN_MODE, /*< ubyte  sensor_mode     */
    SCN_QUAL, /*< short  scan_quality    */
    NVARS_SCN
};
static const varlist VARLIST_SCN[] = {
    { SCN_STIME, "scan_start_time"},
    //{ SCN_ETIME, "scan_end_time"   },
    //{ SCN_MTIME, "ev_mid_time"     },
    { SCN_MSIDE, "HAM_side"},
    { SCN_MODE, "sensor_mode"},
    { SCN_QUAL, "scan_quality"}
};
static var_str_nc *scn[NVARS_SCN];

/* Navigation data, to be read on initialization.
   All values are float; dimensioned [number_of_scans,3] except as noted. */
static const char* NAV_GRP = "navigation_data";

enum nav_var {
    //NAV_QUAT,  /*< att_quat_ev [number_of_scans, quaternion_elements] */
    NAV_ANG, /*< att_ang     */
    NAV_POS, /*< orb_pos_ev  */
    NAV_VEL, /*< orb_vel_ev  */
    //NAV_SOLJ, /*< solar_j2000 */
    //NAV_SOLV, /*< solar_inst  */
    //NAV_SOLD, /*< earth_sun_distance [number_of_scans] */
    //NAV_LUNJ,  /*< lunar_j2000 */
    //NAV_LUNV,  /*< lunar_inst  */
    //NAV_LUND,  /*< earth_moon_distance [number_of_scans] */
    NVARS_NAV
};
static const varlist VARLIST_NAV[] = {
    //{ NAV_QUAT, "att_quat_ev"         },
    { NAV_ANG, "att_ang"}, // (roll, pitch, yaw: degrees)
    { NAV_POS, "orb_pos_ev"}, // (ECR: meters)
    { NAV_VEL, "orb_vel_ev"}, // (ECR: meters/second)
    //{ NAV_SOLJ, "solar_j2000"},
    //{ NAV_SOLV, "solar_inst"},
    //{ NAV_SOLD, "earth_sun_distance"}, // AU
    //{ NAV_LUNJ, "lunar_j2000"         },
    //{ NAV_LUNV, "lunar_inst"          },
    //{ NAV_LUND, "earth_moon_distance" }
};
static var_str_nc *nav[NVARS_NAV]; /* read on initialization */

/* Geolocation data, to be read line-by-line.
   All values are dimensioned as [number_of_lines, number_of_pixels]. */
static const char* GEO_GRP = "geolocation_data";

enum geo_var {
    GEO_LAT, /*< float latitude        */
    GEO_LON, /*< float longitude       */
    GEO_HGT, /*< short height         (-> float) */
    //GEO_RNG,     /*< short range          (-> float) */
    GEO_SENA, /*< short sensor_azimuth (-> float) */
    GEO_SENZ, /*< short sensor_zenith  (-> float) */
    GEO_SOLA, /*< short solar_azimuth  (-> float) */
    GEO_SOLZ, /*< short solar_zenith   (-> float) */
    //GEO_MASK,    /*< ubyte land_water_mask */
    GEO_QUAL, /*< ubyte quality_flag    */
    NVARS_GEO
};
static const varlist VARLIST_GEO[] = {
    { GEO_LAT, "latitude"},
    { GEO_LON, "longitude"},
    { GEO_HGT, "height"},
    //{ GEO_RNG,  "range"           },
    { GEO_SENA, "sensor_azimuth"},
    { GEO_SENZ, "sensor_zenith"},
    { GEO_SOLA, "solar_azimuth"},
    { GEO_SOLZ, "solar_zenith"},
    //{ GEO_MASK, "land_water_mask" },
    { GEO_QUAL, "quality_flag"}
};
static var_str_nc *geo[NVARS_GEO];

/*----- L1B file -----*/

/* Scan Line attributes, to be read on initialization.
   All values are dimensioned as [number_of_scans]. */
enum l1bscn_var {
    L1BSCN_STIME, /*< double scan_start_time */
    //L1BSCN_ETIME, /*< double scan_end_time   */
    //L1BSCN_MTIME, /*< double ev_mid_time     */
    L1BSCN_MODE, /*< ubyte  scan_state_flags    */
    L1BSCN_QUAL, /*< ubyte  scan_quality_flags  */
    NVARS_L1BSCN
};
static const varlist VARLIST_L1BSCN[] = {
    { L1BSCN_STIME, "scan_start_time"},
    //{ L1BSCN_ETIME, "scan_end_time"   },
    //{ L1BSCN_MTIME, "ev_mid_time"     },
    { L1BSCN_MODE, "scan_state_flags"},
    { L1BSCN_QUAL, "scan_quality_flags"}
};
static var_str_nc *l1bscn[NVARS_L1BSCN];

/* Radiometric data, to be read line-by-line.
   All values are dimensioned as [number_of_lines, number_of_pixels]. */
static const char* L1B_GRP = "observation_data";

enum bandtypes {
    RSB, TEB, CIR
};

static const varlist VARLIST_L1B[] = {
    /*** Reflective Solar Bands (aka VIS/SWIR) ***/
    { RSB, "M01"}, //   410 nm
    { RSB, "M02"}, //   443 nm
    { RSB, "M03"}, //   486 nm (blue)
    { RSB, "M04"}, //   551 nm (green)
    { RSB, "M05"}, //   671 nm (red)
    { RSB, "M06"}, //   745 nm
    { RSB, "M07"}, //   862 nm
    { RSB, "M08"}, //  1238 nm
    { CIR, "M09"}, //  1378 nm (cirrus)
    { RSB, "M10"}, //  1601 nm
    { RSB, "M11"}, //  2257 nm
    /*** Thermal Emissive Bands ***/
    { TEB, "M12"}, //  3700 nm
    { TEB, "M13"}, //  4050 nm
    { TEB, "M14"}, //  8550 nm
    { TEB, "M15"}, // 10763 nm
    { TEB, "M16"} // 12013 nm
};
#define MAXBANDS 16
static var_str_nc *l1b[MAXBANDS]; /* ushort -> float */
static var_str_nc *l1bq[MAXBANDS]; /* ubyte; name = band+"_quality_flags" */

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

/*---------------------------------------------------------------------*/

/* Info for any VIIRS input file */

typedef struct {
    int32_t id; /*< NetCDF4 file ID */
    char file[FILENAME_MAX]; /*< file path */
    char title[FILENAME_MAX];
    size_t nscans;
    size_t nlines;
    size_t npixls;
    int orbit_number;
    char start_time[25];
    char end_time[25];
    double *scan_start_time;
} viirs_file;

void print_viirs_file(const viirs_file info) {
    printf("file:\t%s\n", info.file);
    printf("title:\t%s\n", info.title);
    printf("time range:\t%s %s\n", info.start_time, info.end_time);
    printf("orbit: \t%d\n", info.orbit_number);
    printf("nscans =\t%d\n", (int) info.nscans);
    printf("nlines =\t%d\n", (int) info.nlines);
    printf("npixls =\t%d\n", (int) info.npixls);
}

int init_viirs_file(const char filename[FILENAME_MAX],
                    viirs_file *info) {
    int32_t grpid, varid, dimid;
    size_t attlen; // text attribute length
    int status;

    /* open file */
    bzero(info, sizeof(*info));
    strcpy(info->file, filename);
    status = nc_open(info->file, NC_NOWRITE, &info->id);
    if(status != NC_NOERR) {
        printf("-E- init_viirs_file - Could not open netCDF file \"%s\"\n", info->file);
        exit(EXIT_FAILURE);
    }
    
    
    /* load dimensions */
    nc_inq_dimid(info->id, "number_of_scans", &dimid);
    nc_inq_dimlen(info->id, dimid, &info->nscans);
    nc_inq_dimid(info->id, "number_of_lines", &dimid);
    nc_inq_dimlen(info->id, dimid, &info->nlines);
    nc_inq_dimid(info->id, "number_of_pixels", &dimid);
    nc_inq_dimlen(info->id, dimid, &info->npixls);

    /* load global attributes */
    nc_get_att_text(info->id, NC_GLOBAL,
            "title", (char*) &info->title);
    nc_inq_attlen(info->id, NC_GLOBAL, "title", &attlen);
    info->title[attlen] = '\0'; // null terminate

    nc_get_att_text(info->id, NC_GLOBAL,
            "time_coverage_start", (char*) &info->start_time);
    nc_inq_attlen(info->id, NC_GLOBAL, "time_coverage_start", &attlen);
    info->start_time[attlen] = '\0';

    nc_get_att_text(info->id, NC_GLOBAL,
            "time_coverage_end", (char*) &info->end_time);
    nc_inq_attlen(info->id, NC_GLOBAL, "time_coverage_end", &attlen);
    info->end_time[attlen] = '\0';

    // yeah, borrowed this from the example code on UCAR's website...
    nc_type vr_type; /* attribute type */
    size_t vr_len; /* attribute length */
    if ((nc_inq_att(info->id, NC_GLOBAL, "orbit_number", &vr_type, &vr_len) == NC_NOERR)) {
        nc_get_att_int(info->id, NC_GLOBAL,
                "orbit_number", &info->orbit_number);
    } else {
        nc_get_att_int(info->id, NC_GLOBAL,
                "OrbitNumber", &info->orbit_number);
    }

    if ((nc_inq_att(info->id, NC_GLOBAL, "extract_pixel_start", &vr_type, &vr_len) == NC_NOERR)) {
        nc_get_att_int(info->id, NC_GLOBAL,
                "extract_pixel_start", &extract_pixel_start);
        extract_pixel_start--; // Attribute is one-based
        nc_get_att_int(info->id, NC_GLOBAL,
                "extract_pixel_stop", &extract_pixel_stop);
        extract_pixel_stop--; // Attribute is one-based
    }

    /* load start time for each scan */
    nc_inq_grp_ncid(info->id, "scan_line_attributes", &grpid);
    nc_inq_varid(grpid, "scan_start_time", &varid);
    TRYMEM(__FILE__, __LINE__,
            (info->scan_start_time = malloc(info->nscans *
                    sizeof(*info->scan_start_time))));
    nc_get_var_double(grpid, varid, info->scan_start_time);

    return SUCCESS;
}

/*---------------------------------------------------------------------*/

/* Info for VIIRS L1B and GEO files */

int init_viirs_l1bfile(viirs_file l1binfo) {
    size_t i;
    const char* title =
            "VIIRS M-band Reflected Solar Band and Thermal Emissive Band Data";
    grp_str_nc file = { 0 };
    char qaname[FILENAME_MAX + 1];

    /* Verify file type from global attributes */
    if (strcmp(l1binfo.title, title) != 0) {
        fprintf(stderr,
                "Input file %s has unexpected product title!\n"
                "Expected:\t\"%s\"\n"
                "Actual:  \t\"%s\"\n",
                l1binfo.file, title, l1binfo.title);
        return -1;
    }

    /* Determine file contents */
    file.id = l1binfo.id;
    load_grp_nc(&file);

    /* Populate array holding info and data for L1BSCN variables */
    for (i = 0; i < NVARS_L1BSCN; i++) {
        l1bscn[i] = find_var_byname_nc(file, VARLIST_L1BSCN[i].name, SCN_GRP);
        readall_var(l1bscn[i]);
    }

    /* For each L1B variable, */
    for (i = 0; i < MAXBANDS; i++) {

        /* Populate array holding info for L1B variables */
        l1b[i] = find_var_byname_nc(file, VARLIST_L1B[i].name, L1B_GRP);
        // returns NULL if var not found.

        /* also load QA flag info (unscaled) */
        sprintf(qaname, "%s%s", VARLIST_L1B[i].name, "_quality_flags");
        l1bq[i] = find_var_byname_nc(file, qaname, L1B_GRP);

    }

    /* TO DO: free unused memory allocated to grp_str_nc file. */
    return SUCCESS;
}

int init_viirs_geofile(viirs_file geoinfo) {
    size_t i;
    const char* title = "VIIRS M-band Geolocation Data";
    grp_str_nc file = { 0 };

    /* Verify file type from global attributes */
    if (strcmp(geoinfo.title, title) != 0) {
        fprintf(stderr,
                "Input file %s has unexpected product title!\n"
                "Expected:\t\"%s\"\n"
                "Actual:  \t\"%s\"\n",
                geoinfo.file, title, geoinfo.title);
        return -1;
    }

    /* Determine file contents */
    file.id = geoinfo.id;
    load_grp_nc(&file);

    /* Populate array holding info for GEO variables */
    for (i = 0; i < NVARS_GEO; i++) {
        geo[i] = find_var_byname_nc(file, VARLIST_GEO[i].name, GEO_GRP);
    }

    TRY_NC(__FILE__, __LINE__,
            nc_get_att_float(geo[GEO_LAT]->grpid, geo[GEO_LAT]->id, "_FillValue", &latGeoFillValue));
    TRY_NC(__FILE__, __LINE__,
            nc_get_att_float(geo[GEO_LON]->grpid, geo[GEO_LON]->id, "_FillValue", &lonGeoFillValue));
    TRY_NC(__FILE__, __LINE__,
            nc_get_att_short(geo[GEO_SENZ]->grpid, geo[GEO_SENZ]->id, "_FillValue", &senzGeoFillValue));
    TRY_NC(__FILE__, __LINE__,
            nc_get_att_short(geo[GEO_SENA]->grpid, geo[GEO_SENA]->id, "_FillValue", &senaGeoFillValue));
    TRY_NC(__FILE__, __LINE__,
            nc_get_att_short(geo[GEO_SOLZ]->grpid, geo[GEO_SOLZ]->id, "_FillValue", &solzGeoFillValue));
    TRY_NC(__FILE__, __LINE__,
            nc_get_att_short(geo[GEO_SOLA]->grpid, geo[GEO_SOLA]->id, "_FillValue", &solaGeoFillValue));

    // grab the fill values for L2 products
    int status;
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

    /* Populate array holding info and data for SCN variables */
    for (i = 0; i < NVARS_SCN; i++) {
        scn[i] = find_var_byname_nc(file, VARLIST_SCN[i].name, SCN_GRP);
        readall_var(scn[i]);
    }

    /* Populate array holding info and data for NAV variables */
    for (i = 0; i < NVARS_NAV; i++) {
        nav[i] = find_var_byname_nc(file, VARLIST_NAV[i].name, NAV_GRP);
        if (!nav[i]) {
            char *tmpStr = malloc(strlen(VARLIST_NAV[i].name) + 5);
            sprintf(tmpStr, "%s_mid", VARLIST_NAV[i].name);
            nav[i] = find_var_byname_nc(file, tmpStr, NAV_GRP);
            if (!nav[i]) {
                printf("-E- init_viirs_geofile - could not find variable %s in the geo file\n", VARLIST_NAV[i].name);
                exit(EXIT_FAILURE);
            }
            free(tmpStr);
        }
        readall_var(nav[i]);
    }

    /* TO DO: free unused memory allocated to grp_str_nc file. */
    return SUCCESS;
}

/*---------------------------------------------------------------------*/

/* Open and validate L1B and GEO files */

int openl1_viirs_l1b(filehandle *l1file) {
    viirs_file l1binfo, geoinfo;
    int iscan;

    // make sure a GEO file was passed in
    if(!l1file->geofile || !l1file->geofile[0]) {
        printf("-E- openl1b_viirs_nc - VIIRS L1B files require a GEO file.\n");
        exit(EXIT_FAILURE);
    }
    
    /*----- Initialize L1B and GEO info -----*/
    init_viirs_file(l1file->name, &l1binfo);
    init_viirs_file(l1file->geofile, &geoinfo);

    /*----- Populate filehandle structure -----*/
    l1file->ndets = l1binfo.nlines / l1binfo.nscans;
    l1file->nscan = l1binfo.nlines;
    l1file->npix = l1binfo.npixls;
    l1file->terrain_corrected = 1; // presumed.
    l1file->orbit_number = l1binfo.orbit_number;
    strcpy(l1file->spatialResolution, "750 m");

    /*----- Check that input files are compatible -----*/

    /* dimensions */
    if (
            (l1binfo.nscans != geoinfo.nscans) ||
            (l1binfo.nlines != geoinfo.nlines) ||
            (l1binfo.npixls != geoinfo.npixls)) {
        fprintf(stderr, "Geometry mismatch!\n");
        fprintf(stderr, "L1B: nscans = %3d, nlines = %4d, npixls = %4d\n",
                (int) l1binfo.nscans, (int) l1binfo.nlines, (int) l1binfo.npixls);
        fprintf(stderr, "GEO: nscans = %3d, nlines = %4d, npixls = %4d\n",
                (int) geoinfo.nscans, (int) geoinfo.nlines, (int) geoinfo.npixls);
        return -1;
    }

    /* time and orbit */
    if (
            (l1binfo.orbit_number != geoinfo.orbit_number) ||
            (strcmp(l1binfo.start_time, geoinfo.start_time) != 0) ||
            (strcmp(l1binfo.end_time, geoinfo.end_time) != 0)) {
        fprintf(stderr, "Time coverage mismatch!\n");
        fprintf(stderr, "L1B: Orbit %6d, %s - %s\n",
                l1binfo.orbit_number, l1binfo.start_time, l1binfo.end_time);
        fprintf(stderr, "GEO: Orbit %6d, %s - %s\n",
                geoinfo.orbit_number, geoinfo.start_time, geoinfo.end_time);
        return -1;
    }

    /* scan start times */
    for (iscan = 0; iscan < l1binfo.nscans; iscan++)
        if (fabs(l1binfo.scan_start_time[iscan] -
                geoinfo.scan_start_time[iscan]) > DBL_EPSILON) {
            fprintf(stderr, "Scan time mismatch!\n");
            fprintf(stderr, "Scan %d\tL1B=%f\tGEO=%f\tdiff=%f\n",
                    iscan,
                    l1binfo.scan_start_time[iscan],
                    geoinfo.scan_start_time[iscan],
                    l1binfo.scan_start_time[iscan] -
                    geoinfo.scan_start_time[iscan]);
            return -1;
        }

    /*----- Load L1B and GEO info into global variables -----*/
    init_viirs_l1bfile(l1binfo); // l1b, l1bq
    init_viirs_geofile(geoinfo); // scn, nav, geo

    /*----- Ancillary Data -----*/
    // TODO: remove this?
    rdsensorinfo(l1file->sensorID, l1_input->evalmask,
            "Fobar", (void **) &Fobar);

    free(l1binfo.scan_start_time); // free memory allocated
    free(geoinfo.scan_start_time); //  by init_viirs_file()
    return SUCCESS;
}

/***********************************************************************
                              Read Lines
 ***********************************************************************/

int read_var_1line(var_str_nc *var, size_t iline) {
    size_t typesize;
    size_t npixl = 0;
    size_t start[2] = {0, 0}; // assume dims [nline,npixl]
    size_t count[2] = {1, 1};

    if (var->ndims == 2) {

        /* determine # bytes to read */
        npixl = var->dim[1].len;
        TRY_NC(__FILE__, __LINE__,
                nc_inq_type(var->id, var->type, NULL, &typesize));

        /* allocate one line as needed (first call) */
        if (var->data == NULL) {
            TRYMEM(__FILE__, __LINE__,
                    (var->data = calloc(npixl, typesize)));
            var->dim[0].len = 1; // to represent data stored
        }
        /* read data */
        start[0] = iline; // read one line
        count[1] = npixl; // read all pixels
        TRY_NC(__FILE__, __LINE__,
                nc_get_vara(var->grpid, var->id, start, count, var->data));
    }
    return npixl;
}

int scale_short(var_str_nc *var, float* dest) {
    /* cheating here: assume scaling is always short to float. */

    size_t ipix;
    size_t npix = var->dim[1].len;
    short tmpval, fillval;
    short minval, maxval;
    float scale, offset;

    /* load scaling factors */
    TRY_NC(__FILE__, __LINE__,
            nc_get_att_short(var->grpid, var->id, "_FillValue", &fillval));
    TRY_NC(__FILE__, __LINE__,
            nc_get_att_short(var->grpid, var->id, "valid_min", &minval));
    TRY_NC(__FILE__, __LINE__,
            nc_get_att_short(var->grpid, var->id, "valid_max", &maxval));
    TRY_NC(__FILE__, __LINE__,
            nc_get_att_float(var->grpid, var->id, "scale_factor", &scale));
    TRY_NC(__FILE__, __LINE__,
            nc_get_att_float(var->grpid, var->id, "add_offset", &offset));

    for (ipix = 0; ipix < npix; ipix++) {
        tmpval = ((short *) var->data)[ipix];
        if ((tmpval == fillval) ||
                (tmpval < minval) ||
                (tmpval > maxval))
            dest[ipix] = (float) tmpval;
        else
            dest[ipix] = scale * (float) tmpval + offset;
    }

    return SUCCESS;
}

int scale_l1bvals(l1str *l1rec) {

    uint32_t tmpval, fillval;
    uint32_t minval, maxval;
    float scale, offset;
    double f_corr;

    size_t ipb, ipix, npix;
    size_t iband;
    size_t irsb = 0;
    size_t iteb = 0;
    enum bandtypes bandtype;
    //char flag; // ubyte

    /* Initializations */
    npix = l1rec->npix;

    /* Step through all bands */
    for (iband = 0; iband < MAXBANDS; iband++) {
        bandtype = VARLIST_L1B[iband].index;

        if (l1b[iband] == NULL) {
            if (bandtype == TEB) iteb++;
            if (bandtype == RSB) irsb++;
            continue; // skip rest of processing for missing band
        }

        TRY_NC(__FILE__, __LINE__,
                nc_get_att_uint(l1b[iband]->grpid, l1b[iband]->id,
                "_FillValue", &fillval));
        TRY_NC(__FILE__, __LINE__,
                nc_get_att_uint(l1b[iband]->grpid, l1b[iband]->id,
                "valid_min", &minval));
        TRY_NC(__FILE__, __LINE__,
                nc_get_att_uint(l1b[iband]->grpid, l1b[iband]->id,
                "valid_max", &maxval));
        TRY_NC(__FILE__, __LINE__,
                nc_get_att_float(l1b[iband]->grpid, l1b[iband]->id,
                "scale_factor", &scale));
        TRY_NC(__FILE__, __LINE__,
                nc_get_att_float(l1b[iband]->grpid, l1b[iband]->id,
                "add_offset", &offset));

        /* get specific f table cal correction  */
        f_corr = (f_cal_corr == NULL) ? 1.0
                : f_cal_corr[iband][l1rec->detnum][l1rec->mside];

        /*** Thermal Emissive Bands ***/
        if (bandtype == TEB) {
            for (ipix = 0; ipix < npix; ipix++) {
                //ipb = ipix * l1rec->l1file->nbandsir + iteb; // band varies fastest
                ipb = ipix * NBANDSIR + iteb; // band varies fastest
                l1rec->Ltir[ipb] = 0; // init to fill value

                /* skip out-of-bounds values */
                if (strcmp(VARLIST_L1B[iband].name, "M13") == 0) {
                    tmpval = ((uint32_t *) l1b[iband]->data)[ipix];
                } else {
                    tmpval = ((ushort *) l1b[iband]->data)[ipix];
                }
                if ((tmpval == fillval) ||
                        (tmpval < minval) ||
                        (tmpval > maxval))
                    continue;

                /* apply radiance scaling */
                l1rec->Ltir[ipb] = scale * (float) tmpval + offset;

                /* Convert radiance from W/m^2/um/sr to mW/cm^2/um/sr */
                l1rec->Ltir[ipb] /= 10.0;

                /* Apply F-factor */
                l1rec->Ltir[ipb] *= f_corr;

            } // ipix
            iteb++;
        } // TEB

        /*** Cirrus Band ***/
        else if (bandtype == CIR) {
            for (ipix = 0; ipix < npix; ipix++) {
                l1rec->rho_cirrus[ipix] = BAD_FLT; // init to fill value

                /* skip out-of-bounds values */
                tmpval = ((ushort *) l1b[iband]->data)[ipix];
                if ((tmpval == fillval) ||
                        (tmpval < minval) ||
                        (tmpval > maxval))
                    continue;

                /* apply reflectance scaling */
                l1rec->rho_cirrus[ipix] = scale * (float) tmpval + offset;

                /* Normalize reflectance by solar zenith angle */
                l1rec->rho_cirrus[ipix] /= cos(l1rec->solz[ipix] / OEL_RADEG);

                /* Apply F-factor */
                l1rec->rho_cirrus[ipix] *= f_corr;

            } // ipix
        } // CIR

        /*** Reflective Solar Bands ***/
        else { // if (bandtype == RSB)

            l1rec->Fo[irsb] = Fobar[irsb] * l1rec->fsol;

            for (ipix = 0; ipix < npix; ipix++) {
                ipb = ipix * l1rec->l1file->nbands + irsb; // band varies fastest
                l1rec->Lt[ipb] = BAD_FLT; // init to fill value

                /* skip night data for visible bands */
                if (l1rec->solz[ipix] > SOLZNIGHT)
                    continue;

                /* skip out-of-bounds values */
                tmpval = ((ushort *) l1b[iband]->data)[ipix];
                if ((tmpval == fillval) ||
                        (tmpval < minval) ||
                        (tmpval > maxval))
                    continue;

                /* apply reflectance scaling */
                l1rec->Lt[ipb] = scale * (float) tmpval + offset;

                /* convert from reflectance to radiance */
                l1rec->Lt[ipb] *= l1rec->Fo[irsb] / OEL_PI;

                /* Apply F-factor */
                l1rec->Lt[ipb] *= f_corr;

                /* TO DO: flag any suspect values */
                /* TO DO: flag hilt */
                //flag = ((char*)l1bq[iband]->data)[ipix];
                /* if ((flag & 4) && (irsb < 13)) ; */

            } // ipix
            irsb++;
        } // RSB

    } // iband

    return SUCCESS;
}

int readl1_viirs_l1b(filehandle *l1file,
                     const int32_t iline,
                     l1str *l1rec,
                     int lonlat) {
    float nbytes = l1file->npix * sizeof(float);
    size_t ivar;
    size_t ipix;
    size_t iscan;
    size_t i;

    float ang[3]; // degrees
    float pos[3]; // km
    float vel[3]; // km/sec
    float sen_mat[3][3], coeff[10]; // for ocorient & mnorm
    double mnorm[3];

    /*----- Basic info -----*/
    //    l1rec->sensorID = l1file->sensorID;
    l1rec->npix = l1file->npix;
    l1rec->detnum = iline % l1file->ndets;

    /*----- Scan-line values -----*/
    iscan = iline / l1file->ndets;

    /* Scan Start Time */
    double TAI93sec = ((double*) scn[SCN_STIME]->data)[iscan];
    if (TAI93sec < 0.0) {
        l1rec->scantime = BAD_FLT;
    } else {
        l1rec->scantime = tai93_to_unix(TAI93sec);
    }

    /* Mirror Side */
    l1rec->mside = (int32_t) ((char*) scn[SCN_MSIDE]->data)[iscan];

    /* Earth-sun distance correction for this scan */
    int16_t year, day;
    double dsec;
    unix2yds(l1rec->scantime, &year, &day, &dsec);
    int32_t yr = (int32_t) year;
    int32_t dy = (int32_t) day;
    int32_t msec = (int32_t) (dsec * 1000.0);
    double esdist = esdist_(&yr, &dy, &msec);
    l1rec->fsol = pow(1.0 / esdist, 2);

    /* TO DO: check SCN_MODE, SCN_QUAL for both GEO and L1B */
    int8_t scanqual = ((int8_t*) l1bscn[L1BSCN_QUAL]->data)[iscan];

    if (scanqual & 1) {
        l1file->sv_with_moon = 1;
    }

    // if lonlat mode, read only lon, lat and solz
    if (lonlat) {
        read_var_1line(geo[GEO_LAT], iline);
        read_var_1line(geo[GEO_LON], iline);
        read_var_1line(geo[GEO_SOLZ], iline);
    } 
    // else, read all the variables
    else {
        /*----- Geolocation swath values -----*/
        for (ivar = 0; ivar < NVARS_GEO; ivar++)
            read_var_1line(geo[ivar], iline);
    }

    // copy lon, lat and solz into l1rec
    memcpy(l1rec->lat, geo[GEO_LAT]->data, nbytes);
    memcpy(l1rec->lon, geo[GEO_LON]->data, nbytes);
    scale_short(geo[GEO_SOLZ], l1rec->solz);

    // stop reading if lonlat mode, dont need anything else
    if (lonlat)
        return LIFE_IS_GOOD;

    scale_short(geo[GEO_HGT], l1rec->height);
    scale_short(geo[GEO_SOLA], l1rec->sola);
    scale_short(geo[GEO_SENZ], l1rec->senz);
    scale_short(geo[GEO_SENA], l1rec->sena);

    /* Load Angles */
    for (i = 0; i < 3; i++) {
        ang[i] = ((float*) nav[NAV_ANG]->data)[iscan * 3 + i]; // degrees
        pos[i] = ((float*) nav[NAV_POS]->data)[iscan * 3 + i] / 1000.; // m   -> km
        vel[i] = ((float*) nav[NAV_VEL]->data)[iscan * 3 + i] / 1000.; // m/s -> km/s
    }

    /* Check for non-nominal roll, pitch, or yaw */
    /* badatt = */
    /*     (fabs(ang[0]) > MAX_ATTERR) || */
    /*     (fabs(ang[1]) > MAX_ATTERR) || */
    /*     (fabs(ang[2]) > MAX_ATTERR); */

    /* Compute polarization rotation angles */
    ocorient_(pos, vel, ang, sen_mat, coeff);
    for (i = 0; i < 3; i++)
        mnorm[i] = sen_mat[i][0];
    compute_alpha(l1rec->lon, l1rec->lat,
            l1rec->senz, l1rec->sena,
            mnorm, l1rec->npix, l1rec->alpha);

    //check for moon in spaceview port
    //l1file->sv_with_moon = 1;  // as in l1_viirs_h5.c

    /*----- Radiometric swath values -----*/
    for (ivar = 0; ivar < MAXBANDS; ivar++) {
        if (l1b[ivar] == NULL) continue; // skip missing band
        read_var_1line(l1b[ivar], iline);
        read_var_1line(l1bq[ivar], iline);
    }
    scale_l1bvals(l1rec); // get reflectance & radiance
    radiance2bt(l1rec, -1); // calculate brightness temperature

    /*----- Check pixel values -----*/
    for (ipix = 0; ipix < l1rec->npix; ipix++) {
        l1rec->pixnum[ipix] = ipix + extract_pixel_start;
        if (l1rec->scantime < 0)
            l1rec->flags[ipix] |= NAVFAIL;
        if (scanqual & 4)
            l1rec->flags[ipix] |= NAVFAIL;

        // 1 Input_invalid
        // 2 Pointing_bad
        // 4 Terrain_bad
        // check for Input_invalid or Pointing_bad (==3)
        if (((unsigned char*) (geo[GEO_QUAL]->data))[ipix] & 3)
            l1rec->flags[ipix] |= NAVFAIL;

        // check for Terrain_bad (==4)
        if (((unsigned char*) (geo[GEO_QUAL]->data))[ipix] & 4)
            l1rec->flags[ipix] |= NAVWARN;

        flag_bowtie_deleted(l1rec, ipix, extract_pixel_start);
        //l1rec->navwarn[ipix] |= (l1file->sv_with_moon == 1);

        // convert the fill values to the proper value
        if (l1rec->lat[ipix] == latGeoFillValue)
            l1rec->lat[ipix] = latL2FillValue;
        if (l1rec->lon[ipix] == lonGeoFillValue)
            l1rec->lon[ipix] = lonL2FillValue;
        if (l1rec->sena[ipix] == senaGeoFillValue)
            l1rec->sena[ipix] = senaL2FillValue;
        if (l1rec->senz[ipix] == senzGeoFillValue)
            l1rec->senz[ipix] = senzL2FillValue;
        if (l1rec->sola[ipix] == solaGeoFillValue)
            l1rec->sola[ipix] = solaL2FillValue;
        if (l1rec->solz[ipix] == solzGeoFillValue)
            l1rec->solz[ipix] = solzL2FillValue;
    }

    return SUCCESS;
}

void flag_bowtie_deleted(l1str *l1rec, size_t ipix, int extract_offset) {
    int pix = ipix + extract_offset;
    if (pix < AGZONE1 || pix >= AGZONE5) {
        if (l1rec->detnum < 2 || l1rec->detnum > 13) {
            l1rec->flags[ipix] |= BOWTIEDEL;
        }
    } else if ((pix >= AGZONE1 && pix < AGZONE2)
            || (pix >= AGZONE4 && pix < AGZONE5)) {
        if (l1rec->detnum == 0 || l1rec->detnum == 15) {
            l1rec->flags[ipix] |= BOWTIEDEL;
        }
    }
}

/***********************************************************************
                              Cleanup
 ***********************************************************************/

int closel1_viirs_l1b() {
    /* int32_t i; */

    /* Free memory allocated for global variables */
    /* for (i = 0; i < NVARS_SCN; i++) */
    /*     free_vars_nc(scn[i],1) ; */
    /* for (i = 0; i < NVARS_NAV; i++) */
    /*     free_vars_nc(nav[i],1) ; */
    /* for (i = 0; i < NVARS_GEO; i++) */
    /*     free_vars_nc(geo[i],1) ; */
    /* for (i = 0; i < NVARS_L1BSCN; i++) */
    /*     free_vars_nc(l1bscn[i],1) ; */
    /* for (i = 0; i < MAXBANDS; i++) { */
    /*     free_vars_nc(l1b[i], 1) ; */
    /*     free_vars_nc(l1bq[i],1) ; */
    /* } */

    /* Free any memory allocated for input files */
    /* End NetCDF access */

    /* TO DO: make grp_str_nc l1bfile & geofile global?
     Then use free_grp_nc & nc_close.  Won't need free_vars_nc calls above,
     since they point to same var_str_nc memory. */

    return SUCCESS;
}
