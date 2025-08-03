/* ============================================================================ */
/* module l1_hmodis_hdf.c - functions to read MODIS HIRES L1B for MSL12         */
/*                                                                              */
/* Written By:  B. Franz, NASA/SIMBIOS, January 2003.                           */
/* Modified By: J. Gales, NASA/OBPG,    August 2005.                            */
/* Conversion for 16-band MODIS: B. Franz, J. Gales, 2006.                      */
/* Completely restructured by G. Fireman in 2014.                               */
/*                                                                              */
/* ============================================================================ */

#define _XOPEN_SOURCE_EXTENDED /* for strdup() */
#include "libnav.h"
#include "hdf4utils.h"
#include "l1.h"
#include <libgen.h>
#include <hdf.h>
#include <mfhdf.h>
#undef _XOPEN_SOURCE_EXTENDED

#define SCAN_TIME_INTERVAL  1.4771

/***********************************************************************/
#define TRYMEM(file,line,memstat) {                                     \
        if (memstat == NULL) {                                          \
            fprintf(stderr,                                             \
                    "-E- %s line %d: Memory allocation error.\n",       \
                    file, line);                                        \
            exit(1); }                                                  \
    }
/***********************************************************************/

/* SDSs of interest in each type of input file */
typedef struct {
    int32_t index;
    int32_t scandim;
    char *name;
} sdslist;

/* Positioning info to preload from MODIS Geolocation file */
enum REF_SDS {
    GEO_TAISEC, /**< EV start time   */
    GEO_MSIDE, /**< Mirror side     */
    GEO_ANGLES, /**< attitude_angles */
    GEO_MNORM, /**< T_inst2ECR      */
    REF_NUM_SDS /**< last entry = number of SDSs */
};
static const sdslist SDSLIST_REF[] = {
    { GEO_TAISEC, 0, "EV start time"},
    { GEO_MSIDE, 0, "Mirror side"},
    { GEO_ANGLES, 0, "attitude_angles"},
    { GEO_MNORM, 0, "T_inst2ECR"},
    { 0, 0, NULL}
};

/* MODIS Geolocation swath SDSs */
enum GEO_SDS {
    GEO_LON, /**< Longitude     */
    GEO_LAT, /**< Latitude      */
    GEO_HGT, /**< Height        */
    GEO_SOLZ, /**< SolarZenith   */
    GEO_SOLA, /**< SolarAzimuth  */
    GEO_SENZ, /**< SensorZenith  */
    GEO_SENA, /**< SensorAzimuth */
    GEO_NUM_SDS /**< last entry = number of SDSs */
};
static const sdslist SDSLIST_GEO[] = {
    { GEO_LON, 0, "Longitude"}, /* float32 */
    { GEO_LAT, 0, "Latitude"}, /* float32 */
    { GEO_HGT, 0, "Height"}, /* int16   */
    { GEO_SOLZ, 0, "SolarZenith"}, /* int16   */
    { GEO_SOLA, 0, "SolarAzimuth"}, /* int16   */
    { GEO_SENZ, 0, "SensorZenith"}, /* int16   */
    { GEO_SENA, 0, "SensorAzimuth"}, /* int16   */
    { 0, 0, NULL}
};

enum GEO_COEFFS {
    GEO_SCALE, /**< scale_factor */
    GEO_NUM_COEFFS /**< last entry = number of scale attributes */
};
static const char* geo_coeff[] = {/* index with enum GEO_COEFFS */
    "scale_factor"
};

/* MODIS L1B Band groupings */
enum L1B_SDS {
    RSB_250, /**< Bands 1 & 2         */
    RSB_500, /**< Bands 3 - 7         */
    RSB_1KM, /**< all other reflective bands */
    CIR_1KM, /**< Band 26             */
    TEB_1KM, /**< all thermal bands   */
    L1B_NUM_SDS /**< last entry = number of SDSs */
};
static sdslist SDSLIST_1KM[] = {// 1KM file contents
    { RSB_250, 1, "EV_250_Aggr1km_RefSB"},
    { RSB_500, 1, "EV_500_Aggr1km_RefSB"},
    { RSB_1KM, 1, "EV_1KM_RefSB"},
    { CIR_1KM, 0, "EV_Band26"}, // 1st dim=scan, not band
    { TEB_1KM, 1, "EV_1KM_Emissive"},
    { 0, 0, NULL}
};
static sdslist SDSLIST_HKM[] = {// HKM file contents
    { RSB_250, 1, "EV_250_Aggr500_RefSB"},
    { RSB_500, 1, "EV_500_RefSB"},
    { 0, 0, NULL}
};
static sdslist SDSLIST_QKM[] = {// QKM file
    { RSB_250, 1, "EV_250_RefSB"},
    { 0, 0, NULL}
};

enum L1B_COEFFS {
    REFL_SCALE, /**< reflectance_scales  */
    REFL_OFFSET, /**< reflectance_offsets */
    RAD_SCALE, /**< radiance_scales     */
    RAD_OFFSET, /**< radiance_offsets    */
    L1B_NUM_COEFFS /**< last entry = number of scale attributes */
};
static const char* l1b_coeff[] = {/* index with enum L1B_COEFFS */
    "reflectance_scales",
    "reflectance_offsets",
    "radiance_scales",
    "radiance_offsets"
};

/***********************************************************************/

/* Structure holding info for each MODIS input file */
typedef struct {
    int32_t id; /**< HDF4 file ID */
    char file[FILENAME_MAX]; /**< file path */
    char shortname[FILENAME_MAX]; /**< M[OY]D03 or M[OY]D02[1HQ]KM */
    double starttime; /**< granule start time from metadata, in unixtime */
    int16_t resolution; /**< native resolution: 1000, 500 or 250 meters */
    int32_t nlines; /**< (nscans*ndets) */
    int32_t npixls; /**< (usually 1354)*1000/resolution */
    int32_t n_sds; /**< number of science data sets loaded */
    sds_struct *sds; /**< array of data set pointers */
} modis_file;

/***********************************************************************/

/* Structure holding info needed for each swath-based SDS */
typedef struct {
    sds_struct sds; /**< a single SDS = 1 or more bands */
    int16_t resolution; /**< native resolution: 1000, 500 or 250 meters */
    int32_t nbands; /**< number of bands in this SDS */
    int32_t nscans; /**< 1KM scans; usually ~203 */
    int32_t ndets; /**< 10, 20 or 40 detectors */
    int32_t nframes; /**< 1KM frames; usually 1354 */
    int32_t nsf; /**< number of subframes = ndets/10 */
    int32_t scandim; /**< index of scan dimension */
} modis_sds;

/*
  Notes:
  All  GEO arrays are dimensioned [nscans*ndets, nframes]
  Most L1B arrays are dimensioned [nbands, nscans*ndets, nframes*nsf]
  where:
  nscans and nframes are at 1KM resolution
  nsf = #subframes = 1000/resolution
  ndets = 10*nsubframes

  Band 26 drops the 1st dimension.
  Geolocation is inherently 1KM resolution: nbands==1, ndets==10, nsf==1.

  ndets and scandim are needed for reading data one scan at a time.

  There is some redundancy between resolution, ndets, and nsf,
  but it's helpful to have them all defined.
 */

/***********************************************************************/
/* Global variables */

static modis_file file_geo, file_1km, file_hkm, file_qkm;
static modis_sds ref[REF_NUM_SDS], geo[GEO_NUM_SDS], l1b[L1B_NUM_SDS];

/***********************************************************************/

/* Substitute provided string for LAC/1KM part of filename */
int modpath_1km(const char *oldpath, const char* newchars, char* newpath) {
    int32_t status;
    char *ptr = NULL;
    char *tmpbuf = NULL;
    char *filename = NULL;

    /* directory name */
    tmpbuf = strdup(oldpath);
    strcpy(newpath, dirname(tmpbuf));
    strcat(newpath, "/");

    /* filename up to substring */
    tmpbuf = strdup(oldpath);
    filename = basename(tmpbuf);
    status = ((ptr = strstr(filename, "_MODIS.")) == NULL); 
    if (status) {
        status = (((ptr = strstr(filename, "LAC")) == NULL) &&
            ((ptr = strstr(filename, "1KM")) == NULL));  
    
        if (status) {
            fprintf(stderr,
                    "Input file %s does not follow the standard name convention;"
                    " can't find \"%s\" file.\n",
                    oldpath, newchars);
            newpath[0] = '\0';
            return status;
        }
        
        strncat(newpath, filename, ptr - filename);

        /* add new substring and the rest of the filename */
        strcat(newpath, newchars);
        strcat(newpath, ptr + strlen(newchars));
    } else {
        strncat(newpath, filename, ptr - filename + 6);

        /* add new substring and the rest of the filename */
        strcat(newpath, "_");
        strcat(newpath, newchars);
        strcat(newpath, ptr + strlen(newchars) + 3);
    }
    free(tmpbuf);

    return status;
}

/***********************************************************************/

int open_modis_l1bfile(modis_file *l1bfile) {
    int32_t status = FAIL;
    int32_t i;
    char result[FILENAME_MAX];
    sdslist *file_content;

    /* Open input file */
    status = ((l1bfile->id = SDstart(l1bfile->file, DFACC_RDONLY)) == FAIL);
    if (status) {
        fopen_warn(l1bfile->file, __FILE__, __LINE__);
        return status;
    }

    /*----- Read Global Attributes to determine file type -----*/
    get_hdfeos_meta(l1bfile->id, "CoreMetadata.0", "SHORTNAME", result);
    if ((strcmp(result, "MYD021KM") == 0) ||
            (strcmp(result, "MOD021KM") == 0)) { // 1KM
        l1bfile->resolution = 1000;
        file_content = SDSLIST_1KM;
    } else if ((strcmp(result, "MYD02HKM") == 0) ||
            (strcmp(result, "MOD02HKM") == 0)) { // HKM
        l1bfile->resolution = 500;
        file_content = SDSLIST_HKM;
    } else if ((strcmp(result, "MYD02QKM") == 0) ||
            (strcmp(result, "MOD02QKM") == 0)) { // QKM
        l1bfile->resolution = 250;
        file_content = SDSLIST_QKM;
    } else {
        fprintf(stderr,
                "Input file %s has type %s; "
                "does not contain MODIS calibrated radiances.\n",
                l1bfile->file, result);
        return FAIL;
    }
    strcpy(l1bfile->shortname, result);

    /*----- Data set info -----*/
    l1bfile->n_sds = 0;
    while (file_content[l1bfile->n_sds].name)
        l1bfile->n_sds++;
    TRYMEM(__FILE__, __LINE__,
            (l1bfile->sds = calloc(l1bfile->n_sds, sizeof (sds_struct))));
    for (i = 0; i < l1bfile->n_sds; i++)
        init_sds_byname(l1bfile->id, file_content[i].name, &(l1bfile->sds[i]));

    /*----- Native dimensions -----*/
    l1bfile->nlines = l1bfile->sds[RSB_250].dimlen[1]; //(bands,NSCANS*NDETS,nframes*nsf)
    l1bfile->npixls = l1bfile->sds[RSB_250].dimlen[2]; //(bands,nscans*ndets,NFRAMES*NSF)

    /* Note: Every MODIS L1B file contains a 250m-band SDS, aggregated for
       coarser resolutions, so we can use [RSB_250] to reference native
       dimensions for every file. Later we'll use [RSB_250] to index the
       highest-resolution SDS, to specify output dimensions. */

    return SUCCESS;
}

int init_l1b(const char filename[FILENAME_MAX], int32_t *max_resolution) {
    int32_t i, j;
    int32_t status = FAIL;

    /* Start fresh */
    memset(&file_1km, 0, sizeof (modis_file));
    memset(&file_hkm, 0, sizeof (modis_file));
    memset(&file_qkm, 0, sizeof (modis_file));

    /* Initialize L1B 1KM file info */
    strcpy(file_1km.file, filename);
    status = open_modis_l1bfile(&file_1km);
    if (status) {
        fprintf(stderr,
                "Error reading %s; please specify a valid Level 1B file.\n",
                filename);
        return status;
    }

    /* Input L1B file must have 1KM resolution */
    status = (file_1km.resolution != 1000);
    if (status) {
        fprintf(stderr,
                "Input file %s has %dm resolution; "
                "please specify a 1km-resolution Level 1B file.\n",
                filename, file_1km.resolution);
        return status;
    }

    /* Open other L1B files as needed for specified resolution */
    switch (*max_resolution) {

    case 250: // QKM
        status = modpath_1km(filename, "QKM", file_qkm.file);
        if (status)
            return status;
        status = open_modis_l1bfile(&file_qkm);
        if (status) {
            fprintf(stderr, "File not found: %s\n", file_qkm.file);
            fprintf(stderr,
                    "Processing at %im resolution requires a QKM file to be "
                    "present in the same directory as the 1KM L1B file.\n",
                    *max_resolution);
            return status;
        }
        /* no break: also load HKM */

    case 500: // HKM
        status = modpath_1km(filename, "HKM", file_hkm.file);
        if (status)
            return status;
        status = open_modis_l1bfile(&file_hkm);
        if (status) {
            fprintf(stderr, "File not found: %s\n", file_hkm.file);
            fprintf(stderr,
                    "Processing at %im resolution requires a HKM file to be "
                    "present in the same directory as the 1KM L1B file.\n",
                    *max_resolution);
            return status;
        }
        break;

    case -1: // default to 1KM
        *max_resolution = 1000;
        /* no break: proceed as for 1KM */

    case 1000: // 1KM
        /* file_1km has already been initialized */
        break;

    default:
        fprintf(stderr, "Invalid resolution %i; ", *max_resolution);
        fprintf(stderr, "must be 1000, 500 or 250 m.\n");
        status = FAIL;
        exit(status);
        break;
    }

    /* Populate global info for each L1B SDS */
    for (i = 0; i < L1B_NUM_SDS; i++) {
        l1b[i].sds = file_1km.sds[i]; // initialize to 1KM inputs
        l1b[i].resolution = file_1km.resolution;
        l1b[i].nscans = file_1km.nlines / 10;
        l1b[i].nframes = file_1km.npixls;
        l1b[i].scandim = SDSLIST_1KM[i].scandim;
    }
    if (*max_resolution < 1000) {
        l1b[RSB_500].sds = file_hkm.sds[RSB_500]; // HKM or QKM
        l1b[RSB_500].resolution = file_hkm.resolution;

        l1b[RSB_250].sds = file_hkm.sds[RSB_250]; // HKM only
        l1b[RSB_250].resolution = file_hkm.resolution;
    }
    if (*max_resolution < 500) {
        l1b[RSB_250].sds = file_qkm.sds[RSB_250]; // QKM only
        l1b[RSB_250].resolution = file_qkm.resolution;
    }

    /* Derived values */
    for (i = 0; i < L1B_NUM_SDS; i++) {
        l1b[i].nbands =
                (l1b[i].scandim == 0) ? 1 // Band 26 only
                : l1b[i].sds.dimlen[0];
        l1b[i].nsf = 1000 / l1b[i].resolution;
        l1b[i].ndets = 10 * l1b[i].nsf;

        /* Load scaling coefficient attribute(s) */
        TRYMEM(__FILE__, __LINE__,
                (l1b[i].sds.atts = calloc(L1B_NUM_COEFFS, sizeof (att_struct))));
        l1b[i].sds.natts = L1B_NUM_COEFFS;
        for (j = 0; j < L1B_NUM_COEFFS; j++)
            load_att_byname(l1b[i].sds.id, l1b_coeff[j], &l1b[i].sds.atts[j]);
    }

    return status;
}

/***********************************************************************/

int open_modis_geofile(modis_file *geofile) {
    int32_t status = FAIL;
    int32_t i;
    char result[FILENAME_MAX];

    /* Open input file */
    status = ((geofile->id = SDstart(geofile->file, DFACC_RDONLY)) == FAIL);
    if (status) {
        fopen_warn(geofile->file, __FILE__, __LINE__);
        return status;
    }

    /*----- Read Global Attributes to determine file type -----*/
    get_hdfeos_meta(geofile->id, "CoreMetadata.0", "SHORTNAME", result);
    if ((strcmp(result, "MYD03") == 0) ||
            (strcmp(result, "MOD03") == 0)) { // GEO
    } else {
        fprintf(stderr, "Input file %s has type %s; "
                "does not contain MODIS geolocation.\n",
                geofile->file, result);
        return FAIL;
    }
    strcpy(geofile->shortname, result);

    /*----- Data set info -----*/
    geofile->n_sds = GEO_NUM_SDS;
    TRYMEM(__FILE__, __LINE__,
            (geofile->sds = calloc(GEO_NUM_SDS, sizeof (sds_struct))));
    for (i = 0; i < GEO_NUM_SDS; i++)
        init_sds_byname(geofile->id, SDSLIST_GEO[i].name, &(geofile->sds[i]));

    /*----- Native dimensions -----*/
    geofile->nlines = geofile->sds[GEO_LAT].dimlen[0];
    geofile->npixls = geofile->sds[GEO_LAT].dimlen[1];

    return SUCCESS;
}

int init_geo(const char filename[FILENAME_MAX]) {
    int32_t i, j;
    int32_t status = FAIL;

    /* Start fresh */
    memset(&file_geo, 0, sizeof (modis_file));

    /* Initialize GEO file info */
    strcpy(file_geo.file, filename);
    status = open_modis_geofile(&file_geo);
    if (status) {
        fprintf(stderr,
                "Error reading %s; please specify a valid geolocation file.\n",
                filename);
        return status;
    }

    /* Populate global info for each GEO SDS */
    for (i = 0; i < GEO_NUM_SDS; i++) {
        geo[i].sds = file_geo.sds[i];
        geo[i].resolution = 1000;
        geo[i].nbands = 1;
        geo[i].nscans = file_geo.nlines / 10;
        geo[i].ndets = 10;
        geo[i].nframes = file_geo.npixls;
        geo[i].nsf = 1;
        geo[i].scandim = SDSLIST_GEO[i].scandim;

        /* Load scaling coefficient attribute(s) */
        TRYMEM(__FILE__, __LINE__,
                (geo[i].sds.atts = calloc(GEO_NUM_COEFFS, sizeof (att_struct))));
        geo[i].sds.natts = GEO_NUM_COEFFS;
        for (j = 0; j < GEO_NUM_COEFFS; j++)
            load_att_byname(geo[i].sds.id, geo_coeff[j], &geo[i].sds.atts[j]);
    }

    /*----- Read preloaded variables -----*/
    for (i = 0; i < REF_NUM_SDS; i++) {
        init_sds_byname(file_geo.id, SDSLIST_REF[i].name, &(ref[i].sds));
        readall_sds(&(ref[i].sds));
    }

    return status;
}

int read_sds_1scan(sds_struct *sds,
        int32_t iscan,
        int32_t scandim,
        int32_t per_scan) {
    int32_t i;
    int32_t status = FAIL;
    int32_t nvals = 1;
    int32_t *start, *edges;

    /* don't bother if SDS handle is invalid */
    if (sds->id != FAIL) {

        /* set array indices */
        TRYMEM(__FILE__, __LINE__,
                (start = malloc(sds->ndims * sizeof (int32_t))));
        TRYMEM(__FILE__, __LINE__,
                (edges = malloc(sds->ndims * sizeof (int32_t))));
        for (i = 0; i < sds->ndims; i++) {
            if (i == scandim) {
                start[i] = iscan * per_scan;
                edges[i] = per_scan;
            } else {
                start[i] = 0;
                edges[i] = sds->dimlen[i];
            }
            nvals *= edges[i];
        }

        /* allocate one scan as needed (first call) */
        if (sds->data == NULL)
            TRYMEM(__FILE__, __LINE__,
                (sds->data = malloc(nvals * DFKNTsize(sds->ntype))));

        /* read data */
        status = SDreaddata(sds->id, start, NULL, edges, sds->data);
        free(start);
        free(edges);
    }
    if (status != FAIL)
        return nvals;
    else
        return status;
}

/***********************************************************************/
/*  Everything below this line is specific to OCSSW implementation.    */
/***********************************************************************/
#define MIN_BAD_SI  65500
#define MAX_ATTERR  0.017453293 /* radians, ~= 1 degree */
#define LOOPI for (i = 0; i < nvals; i++)

/* Global variables from l1_hmodis_hdf.c */
static int32_t spixl = 0; /* first valid pixel of extracted granule */
static float *Fobar;

/* RSMAS corrections to Terra IR radiances (20,22,23,27,28,29,31,32) */
/* latband v6 sst coefficients were estimated from Bt's from which these */
/* corrections have been applied, so we have to keep using them. */
static float radoff[][10] = {
    {-0.000972, 0.014200, 0.000022, 0.000238, -0.000003,
        0.003760, 0.002337, 0.000084, 0.008848, 0.005050},
    {-0.002300, 0.000600, 0.000236, -0.000280, -0.000003,
        0.002798, 0.003496, 0.018035, 0.002942, 0.001787},
    {-0.003833, 0.003657, 0.002118, 0.000798, -0.000003,
        0.000208, 0.000399, 0.000553, 0.000258, 0.021128},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {-0.000423, -0.000242, -0.000330, -0.000065, -0.000001,
        -0.000006, 0.000064, 0.000398, 0.000362, 0.000322},
    {-0.000433, -0.000246, -0.000222, -0.000148, -0.000001,
        0.000068, 0.000221, 0.000221, 0.000437, 0.000303}
};
static float mfact[] = {0.0, 0.0, 0.0005, 0.0, 0.0, 0.0, 0.0028, 0.00278};

/**
 *
 */
typedef struct {
    int16_t sdsindex; /**< which SDS contains the band */
    int16_t bandindex; /**< index of band within SDS */
    int16_t flagindex; /**< starting index within "all-detector" arrays */
    int16_t wavelength; /**< (nm) */
    char *bandname; /**< for clarity */
} bandinfo;

/* Bands used by l2gen, in expected order */
static const bandinfo BANDLIST[] = {
    /*** Reflective Solar Bands ***/
    { RSB_1KM, 0, 180, 412, "8"}, //  0
    { RSB_1KM, 1, 190, 443, "9"}, //  1
    { RSB_500, 0, 80, 469, "3"}, //  2
    { RSB_1KM, 2, 200, 488, "10"}, //  3
    { RSB_1KM, 3, 210, 531, "11"}, //  4
    { RSB_1KM, 4, 220, 547, "12"}, //  5
    { RSB_500, 1, 100, 555, "4"}, //  6
    { RSB_250, 0, 0, 645, "1"}, //  7
    { RSB_1KM, 5, 230, 667, "13lo"}, //  8
    { RSB_1KM, 7, 250, 678, "14lo"}, //  9
    { RSB_1KM, 9, 270, 748, "15"}, // 10
    { RSB_250, 1, 40, 859, "2"}, // 11
    { RSB_1KM, 10, 280, 869, "16"}, // 12
    { RSB_500, 2, 120, 1240, "5"}, // 13
    { RSB_500, 3, 140, 1640, "6"}, // 14
    { RSB_500, 4, 160, 2130, "7"}, // 15
    /*** Thermal Emissive Bands ***/
    { TEB_1KM, 0, 330, 3750, "20"}, // 16
    { TEB_1KM, 2, 350, 3959, "22"}, // 17
    { TEB_1KM, 3, 360, 4050, "23"}, // 18
    { TEB_1KM, 6, 390, 6715, "27"}, // 19
    { TEB_1KM, 7, 400, 7325, "28"}, // 20
    { TEB_1KM, 8, 410, 8550, "29"}, // 21
    { TEB_1KM, 10, 430, 11000, "31"}, // 22
    { TEB_1KM, 11, 440, 12000, "32"}, // 23
    /*** Cirrus Band ***/
    { CIR_1KM, 0, 320, 1375, "26"} // 24
};

/***********************************************************************/

/**
 *
 * @param[in,out] l1file
 * @return
 */
int openl1_modis(filehandle *l1file) {
    char geometa[EOSMETALEN] = "";
    char result[EOSMETALEN] = "";
    char datetime[32] = "";
    int32_t resolution;
    int32_t intval;
    uint32_t gflags[8];

    /*----- Initialize L1B and GEO info -----*/
    if (init_l1b(l1file->name, &l1_input->resolution))
        exit(EXIT_FAILURE);
    if (init_geo(l1file->geofile))
        exit(EXIT_FAILURE);

    /* Make sure this is not an ocean-band subsetted L1A file */
    int i;
    for (i = 0; i < file_1km.n_sds; i++) {
        if (strcmp(file_1km.sds[i].name, "EV_250_Aggr1km_RefSB") == 0) {
            if (SDfindattr(file_1km.sds[i].id, "Rescaled Ocean R") != FAIL) {
                fprintf(stderr, "\n"
                        "This L1B file contains only the ocean band subset.\n"
                        "Processing is no longer supported; exiting.\n");
                exit(EXIT_FAILURE);
            }
            break;
        }
    }

    /*----- Populate filehandle structure -----*/
    l1file->ndets = l1b[RSB_250].ndets;
    l1file->nscan = l1b[RSB_250].nscans * l1b[RSB_250].ndets;
    l1file->npix = l1b[RSB_250].nframes * l1b[RSB_250].nsf;

    /* Capture processing resolution */
    resolution = l1b[RSB_250].resolution;
    if (resolution == 1000)
        strcpy(l1file->spatialResolution, "1 km");
    else
        sprintf(l1file->spatialResolution, "%d m", resolution);

    if (want_verbose) {
        printf("Processing at %d meter resolution.\n", resolution);
        printf("    1000-meter file: %s\n", file_1km.file);
        if (resolution < 1000) {
            printf("     500-meter file: %s\n", file_hkm.file);
            if (resolution < 500)
                printf("     250-meter file: %s\n", file_qkm.file);
        }
        printf("\n");
    }

    /*----- Read Level 1B metadata -----*/
    get_hdfeos_meta(file_1km.id, "CoreMetadata.0", "INPUTPOINTER", result);
    strcpy(l1_input->calfile, result);

    /*----- Read Geolocation metadata -----*/
    if (read_att(file_geo.id, "CoreMetadata.0", geometa)) {
        fprintf(stderr,
                "-E- %s line %d: Error reading CoreMetadata.0 from %s.\n",
                __FILE__, __LINE__, file_geo.file);
        exit(EXIT_FAILURE);
    }

    /* Orbit info */
    parse_odl(geometa, "ORBITNUMBER", result);
    if (strlen(result) > 0)
        l1file->orbit_number = (int32_t) atoi(result);
    parse_odl(geometa, "EQUATORCROSSINGLONGITUDE", result);
    if (strlen(result) > 0)
        l1file->orbit_node_lon = atof(result);

    parse_odl(geometa, "EQUATORCROSSINGDATE", result);
    strcpy(datetime, result);
    strcat(datetime, "T");
    parse_odl(geometa, "EQUATORCROSSINGTIME", result);
    strcat(datetime, result);
    l1file->node_crossing_time = isodate2unix(datetime);

    /* Granule start time */
    parse_odl(geometa, "RANGEBEGINNINGDATE", result);
    strcpy(datetime, result);
    strcat(datetime, "T");
    parse_odl(geometa, "RANGEBEGINNINGTIME", result);
    strcat(datetime, result);
    file_geo.starttime = isodate2unix(datetime);

    /* Terrain height */
    if (read_att(file_geo.id, "Cumulated gflags", gflags))
        fprintf(stderr, "-W- %s line %d: Error reading gflags from %s.\n",
            __FILE__, __LINE__, file_geo.file);
    else
        l1file->terrain_corrected = (((int32_t*) gflags)[5] == 0);

    /*** Adjust dimensions for extracted granules ***/
    /*
      When Level 1 MODIS granules are extracted, the full scan
      width is retained but pixels outside the ROI are set to fill
      value.  We must keep track of pixel offset and count in
      order to load the ROI correctly.  Line offset and count are
      handled in the main() loop.
     */
    if (read_att(file_1km.id, "Extract Pixel Count", &intval) == 0
            && intval > 0) {
        int32_t sline, nline, npixl;
        npixl = intval * l1b[RSB_250].nsf;
        read_att(file_1km.id, "Extract Pixel Offset", &intval);
        spixl = intval * l1b[RSB_250].nsf;
        read_att(file_1km.id, "Extract Line Offset", &intval);
        sline = intval * l1b[RSB_250].ndets;
        if (want_verbose) {
            nline = l1file->nscan;
            fprintf(stdout,
                    "File was generated from L1A extract of size %d x %d\n",
                    npixl, nline);
            fprintf(stdout, "  Pixels %d - %d\n", spixl + 1, spixl + npixl);
            fprintf(stdout, "  Lines  %d - %d\n", sline + 1, sline + nline);
        }
        l1file->npix = npixl;
    }

    /* Read reflectance to radiance conversion factors */
    rdsensorinfo(l1file->sensorID, l1_input->evalmask, "Fobar", (void **) &Fobar);

    return SUCCESS;
}

/***********************************************************************/

/**  Holds an entire scan of data, interpolated to output resolution.
     Geolocation values are scaled, and have data type same as in l1rec. */
static struct {
    int32_t iscan; /**< scan number (0-based) */
    int32_t nvals; /**< number of values per variable per scan */
    int32_t npix; /**< number of values per variable per line */
    int32_t nbands; /**< number of bands */

    /* read from input files */
    double taisec; /**< EV start time */
    int32_t mside; /**< Mirror side */
    double angles[3]; /**< attitude_angles */
    double mnorm[3]; /**< T_inst2ECR (need only 1st 3 values) */
    float *lon; /**< Longitude */
    float *lat; /**< Latitude */
    float *hgt; /**< Height */
    float *solz; /**< SolarZenith */
    float *sola; /**< SolarAzimuth */
    float *senz; /**< SensorZenith */
    float *sena; /**< SensorAzimuth */
    float *allbands; /**< Scaled Integers for all bands */
} scan;

/**
 * Allocates memory for global structure "scan", and initializes some constants
 * @param[in] nvals number of values per variable per scan
 * @return
 */
int alloc_scan(const int32_t nvals) {

    /* set stack variables */
    memset(&scan, 0, sizeof (scan));
    scan.mside = -1;
    scan.nvals = nvals;
    scan.npix = nvals / l1b[RSB_250].ndets;
    scan.nbands = sizeof (BANDLIST) / sizeof (BANDLIST[0]); // works on stack vars only

    /* allocate memory for dynamic variables */
    TRYMEM(__FILE__, __LINE__, (scan.lon = malloc(nvals * sizeof (*scan.lon))));
    TRYMEM(__FILE__, __LINE__, (scan.lat = malloc(nvals * sizeof (*scan.lat))));
    TRYMEM(__FILE__, __LINE__, (scan.hgt = malloc(nvals * sizeof (*scan.hgt))));
    TRYMEM(__FILE__, __LINE__, (scan.solz = malloc(nvals * sizeof (*scan.solz))));
    TRYMEM(__FILE__, __LINE__, (scan.sola = malloc(nvals * sizeof (*scan.sola))));
    TRYMEM(__FILE__, __LINE__, (scan.senz = malloc(nvals * sizeof (*scan.senz))));
    TRYMEM(__FILE__, __LINE__, (scan.sena = malloc(nvals * sizeof (*scan.sena))));
    TRYMEM(__FILE__, __LINE__,
            (scan.allbands = malloc(scan.nbands * nvals * sizeof (*scan.allbands))));

    return SUCCESS;
}

/***********************************************************************/

typedef struct {
    int16_t nrad; /**< number of radiance levels */
    int16_t nsf; /**< number of subframes at highest resolution */
    double *SI_val; /**< radiance level converted to scaled integer */
    double *corr; /**< subframe correction factors */
} sfcorr_table;

void read_sfcorrtables(sfcorr_table *sfcorr, int32_t sensorID) {
    char *dataroot;
    char file[FILENAME_MAX] = "";
    char sensorName[32];
    FILE *fp = NULL;
    char line[81];
    int32_t iband, isf, nsf, irad, nrad;
    size_t isds, ical;
    att_struct *atts;
    double scale, offset;
    double radiance, SI, radiance2, SI2, rad_corr[4];

    /* calibration directory path */
    if ((dataroot = getenv("OCDATAROOT")) == NULL) {
        printf("OCDATAROOT environment variable is not defined.\n");
        exit(1);
    }

    /* for each band, */
    /* get max possible subframes */
    if (want_verbose)
        printf("\nLoading subframe correction tables:\n");
    for (iband = 0; iband < scan.nbands; iband++) {
        nrad = 0;
        isds = BANDLIST[iband].sdsindex;
        if (isds == RSB_250)
            nsf = 4;
        else if (isds == RSB_500)
            nsf = 2;
        else
            nsf = 1;

        /* if a hi-res band, */
        if (nsf > 1) {

            /* get radiance factors */
            atts = l1b[isds].sds.atts;
            ical = BANDLIST[iband].bandindex;
            scale = (double) ((float*) atts[RAD_SCALE].data)[ical];
            offset = (double) ((float*) atts[RAD_OFFSET].data)[ical];

            strcpy(sensorName, sensorId2SensorName(sensorID));
            lowcase(sensorName);
            
            /* open subframe correction file */
            sprintf(file, "%s/%s/%s/cal/subframecor_%s_%d.dat",
                    dataroot,
                    sensorId2SensorDir(sensorID),
                    subsensorId2SubsensorDir(sensorId2SubsensorId(sensorID)),
                    sensorName,
                    BANDLIST[iband].wavelength);
            if (want_verbose)
                printf("    %s\n", file);
            if ((fp = fopen(file, "r")) == NULL) {
                fprintf(stderr,
                        "-E- %s line %d: unable to open %s for reading\n",
                        __FILE__, __LINE__, file);
                exit(1);
            }

            /* find number of radiance levels and allocate space */
            while (fgets(line, 80, fp))
                if (strncmp(line, "Number of radiance levels", 25) == 0) {
                    sscanf(line, "Number of radiance levels: %d", &nrad);
                    TRYMEM(__FILE__, __LINE__,
                            (sfcorr[iband].SI_val =
                            malloc(nrad * sizeof (*sfcorr[iband].SI_val))));
                    TRYMEM(__FILE__, __LINE__,
                            (sfcorr[iband].corr =
                            malloc(nrad * 4 * sizeof (*sfcorr[iband].corr))));
                    break;
                }
            while (fgets(line, 80, fp)) /* skip over rest of header */
                if (strncmp(line, "/end_header", 11) == 0)
                    break;

            /* read radiances and correction factors */
            for (irad = 0; irad < nrad; irad++) {
                if (fgets(line, 80, fp) == NULL) {
                    fprintf(stderr,
                            "-E- %s line %d: File %s contains only %d data records\n",
                            __FILE__, __LINE__, file, irad);
                    exit(1);
                }
                sscanf(line, "%lf %lf %lf %lf %lf", &radiance,
                        &rad_corr[0], &rad_corr[1], &rad_corr[2], &rad_corr[3]);

                /* Scale reference radiances and corrections to SI
                   space so we don't have to scale each pixel later.
                   Also invert factor so we later multiply instead of
                   divide.
                 */
                SI = radiance / scale + offset;
                sfcorr[iband].SI_val[irad] = SI;
                for (isf = 0; isf < 4; isf++) {
                    radiance2 = radiance / rad_corr[isf];
                    SI2 = radiance2 / scale + offset;
                    sfcorr[iband].corr[4 * irad + isf] = SI2 / SI;
                }

            }
            fclose(fp);

        } // hi-res band
        sfcorr[iband].nrad = nrad;
        sfcorr[iband].nsf = nsf;

    } // iband

    if (want_verbose)
        printf("Subframe destriping corrections enabled.\n\n");
}

void subframe_correction(const int32_t sensorID,
        const int32_t iband,
        const modis_sds mds,
        double *data) {
    /* Applies subframe destriping to an entire full-width scan at once,
       so we don't need to correct for any extract offset. */
    sfcorr_table lut;
    int32_t i;
    int32_t idet, iframe, isf, ip;
    double x; // input value (scaled radiance)
    double y; // interpolated subframe correction

    static sfcorr_table *sfcorr = NULL;

    /* Initialize subframe correction tables */
    if (sfcorr == NULL) {
        TRYMEM(__FILE__, __LINE__,
                (sfcorr = malloc(scan.nbands * sizeof (*sfcorr))));
        read_sfcorrtables(sfcorr, sensorID);
    }
    lut = sfcorr[iband];

    /* Destripe hi-res bands at native resolution only */
    if ((lut.nsf > 1) && (lut.nrad > 0) && (lut.nsf == mds.nsf)) {
        for (idet = 0; idet < mds.ndets; idet++)
            for (iframe = 0; iframe < mds.nframes; iframe++)
                for (isf = 0; isf < mds.nsf; isf++) {
                    ip = (idet * mds.nframes + iframe) * mds.nsf + isf;

                    if (data[ip] < MIN_BAD_SI) {
                        x = data[ip];

                        /* find closest lookup table indices */
                        /* TO DO: implement more efficient lookup */
                        for (i = 0; i < lut.nrad - 1; i++) {
                            if (x < lut.SI_val[i])
                                break;
                        }

                        /* use end value if out-of-bounds */
                        if (i == 0 || i == (lut.nrad - 1))
                            y = lut.corr[4 * i + isf];

                            /* otherwise, interpolate between two closest values */
                        else {
                            double x0 = lut.SI_val[i - 1];
                            double x1 = lut.SI_val[i];
                            double y0 = lut.corr[4 * (i - 1) + isf];
                            double y1 = lut.corr[4 * i + isf];
                            y = y0 + (x - x0) * (y1 - y0) / (x1 - x0);
                        }

                        /* Apply the destriping correction */
                        data[ip] *= (float) y;
                    }
                } // isf, iframe, idet loops
    } // eligible band
}

/***********************************************************************/

typedef struct {
    int32_t ndets; /**< number of detectors at native resolution */
    int32_t *prev; /**< previous valid detector for each detector */
    int32_t *next; /**< next valid detector for each detector */
} fill_table;

void init_detfill(fill_table *fill) {
    int32_t ndets; // number of detectors
    int32_t isds, iband, idet, i0, x0, x1;
    char dead[EOSMETALEN] = "";

    /* Read "Dead Detector List" global attribute */
    if (read_att(file_1km.id, "Dead Detector List", dead)) {
        fprintf(stderr,
                "-E- %s line %d: Error reading \"Dead Detector List\" from %s.\n",
                __FILE__, __LINE__, file_1km.file);
        exit(EXIT_FAILURE);
    }

    /* Load values into dead detector table */
    for (iband = 0; iband < scan.nbands; iband++) {
        i0 = BANDLIST[iband].flagindex;

        /* Find number of detectors at native resolution */
        isds = BANDLIST[iband].sdsindex;
        if (isds == RSB_250)
            ndets = 40;
        else if (isds == RSB_500)
            ndets = 20;
        else
            ndets = 10;
        fill[iband].ndets = ndets;

        /* Allocate space for arrays */
        TRYMEM(__FILE__, __LINE__,
                (fill[iband].prev = malloc(ndets * sizeof (*fill->prev))));
        TRYMEM(__FILE__, __LINE__,
                (fill[iband].next = malloc(ndets * sizeof (*fill->next))));

        /* Find previous and next good detectors */
        /*
          These arrays hold indices for the nearest previous and next
          valid detector for each detector.  Since detector arrays are
          small, opt for code clarity over efficiency.

          Special values:
          prev=idet=next: detector is alive
          prev = -1   : no previous good detector
          next = ndets: no next good detector
         */
        for (idet = 0; idet < ndets; idet++) {
            x0 = idet; // previous
            while ((x0 > -1) && dead[i0 + x0]) x0--;
            fill[iband].prev[idet] = x0;
            x1 = idet; // next
            while ((x1 < ndets) && dead[i0 + x1]) x1++;
            fill[iband].next[idet] = x1;
        }

    } // iband loop
}

void fill_dead_detectors(const int32_t iband, const modis_sds mds, double *data) {

    int32_t idet, ndets, ipix, npix;
    size_t nbytes; // bytes per line (detector)

    int32_t *prev, *next; // previous and next valid detector for each det
    int32_t x, x0, x1; // index into data array for current, prev, next det
    double y0, y1; // previous and next valid value for current detector
    double dx, dy; // change in det, value over valid interval

    static fill_table *fill = NULL;

    /* Initialize previous and next valid detector array table */
    if (fill == NULL) {
        TRYMEM(__FILE__, __LINE__, (fill = malloc(scan.nbands * sizeof (*fill))));
        init_detfill(fill);
    }

    /* Fill dead detectors at native resolution only */
    ndets = fill[iband].ndets;
    if (ndets == mds.ndets) {

        prev = fill[iband].prev;
        next = fill[iband].next;
        npix = mds.nframes * mds.nsf; // # values per line (detector)
        nbytes = npix * sizeof (*data);

        /* Substitute data values as appropriate for each detector: */
        for (idet = 0; idet < ndets; idet++) {
            if (prev[idet] != next[idet]) { // continue for bad dets only

                /* Between two good detectors: interpolate */
                if ((prev[idet] > -1) && (next[idet] < ndets)) {
                    dx = next[idet] - prev[idet];

                    for (ipix = 0; ipix < npix; ipix++) {
                        x = idet * npix + ipix;
                        x0 = prev[idet] * npix + ipix;
                        x1 = next[idet] * npix + ipix;
                        y0 = data[x0];
                        y1 = data[x1];

                        /* if both values are valid, interpolate */
                        if ((y0 < MIN_BAD_SI) && (y1 < MIN_BAD_SI)) {
                            dy = y1 - y0;
                            data[x] = y0 + (idet - prev[idet]) * (dy / dx);
                        }
                            /* otherwise copy adjacent valid value */
                        else if (y0 < MIN_BAD_SI) data[x] = data[x0];
                        else if (y1 < MIN_BAD_SI) data[x] = data[x1];

                            /* if none, copy maximum adjacent flag */
                        else
                            data[x] = MAX(data[x0], data[x1]);
                    }
                }
                    /* No previous good detector: copy next value or flag */
                else if ((prev[idet] == -1) && (next[idet] < ndets))
                    memcpy(&data[idet * npix], &data[next[idet] * npix], nbytes);

                    /* No next good detector: copy previous value or flag */
                else if ((prev[idet] > -1) && (next[idet] == ndets))
                    memcpy(&data[idet * npix], &data[prev[idet] * npix], nbytes);

                /* Default case: no good detectors = no change. */

            } // dead detector
        } // idet loop
    } // band is at native resolution
}

/***********************************************************************/

double interp_bilin(double val[4], double dx, double dy) {
    return (1 - dy) * (val[0] * (1 - dx) + val[1] * dx)
            + dy * (val[2] * (1 - dx) + val[3] * dx);
}

double interp_modis_EVdata(double val[4], double dx, double dy) {
    int32_t i, ngood = 0;
    double sum_val = 0;
    double max_val = 0;
    double result;

    /* interpolation: account for flagged values */
    for (i = 0; i < 4; i++) {
        if (val[i] < MIN_BAD_SI) {
            sum_val = sum_val + val[i];
            ngood++;
        }
        if (val[i] > max_val)
            max_val = val[i];
    }

    /* If we have four non-flagged samples, bilinearly interpolate. */
    /* If all samples are flagged, set to the max flag value.       */
    /* If we have 1-3 samples, use average of unflagged values.     */
    switch (ngood) {
    case 4:
        result = interp_bilin(val, dx, dy);
        break;
    case 0:
        result = max_val;
        break;
    default:
        /* TO DO: should be WEIGHTED avg of good values */
        result = sum_val / ngood;
        break;
    }
    return result;
}

double interp_modis_Longitude(double val[4], double dx, double dy) {
    int32_t i;
    double min_val, max_val;

    /* longitude interpolation: account for dateline */
    min_val = max_val = val[0];
    for (i = 1; i < 4; i++) {
        if (val[i] < min_val)
            min_val = val[i];
        else if (val[i] > max_val)
            max_val = val[i];
    }
    if (max_val > 180.0 + min_val) {
        for (i = 0; i < 4; i++) {
            if (val[i] < 0.0)
                val[i] += 360.0;
        }
    }
    return interp_bilin(val, dx, dy);
}

double (*interp_modis_var)(double val[4], double dx, double dy);

int modis_interp(const modis_sds s,
        const int32_t sdsband,
        const int32_t newres,
        double *olddata,
        double *newdata) {

    int32_t ny = s.ndets;
    int32_t nx = s.nframes * s.nsf;
    int32_t nvals = ny * nx;
    int32_t oldres = s.resolution;

    int32_t nxnew, nynew;
    int32_t ixnew, iynew;

    int32_t x0, y0, x1, y1; // 4 points in old array
    double fx, fy;
    double dx, dy;
    double val[4];
    int32_t i;
    double factor, xoffset, yoffset;
    int32_t maxi = 2; /* 1 duplicates end pixel; 2 extrapolates past end */

    int32_t sdsoffset = nvals * sdsband; // offset within SDS
    short alloc_old;

    /* Cast input data to double precision */
    alloc_old = (olddata == NULL);
    if (alloc_old) {
        TRYMEM(__FILE__, __LINE__, (olddata = malloc(nvals * sizeof (*olddata))));

        switch (s.sds.ntype) {
        case DFNT_FLOAT32:
            LOOPI olddata[i] = (double) ((float*) s.sds.data)[sdsoffset + i];
            break;
        case DFNT_INT16:
            LOOPI olddata[i] = (double) ((int16_t*) s.sds.data)[sdsoffset + i];
            break;
        case DFNT_UINT16:
            LOOPI olddata[i] = (double) ((uint16_t*) s.sds.data)[sdsoffset + i];
            break;
        default:
            fprintf(stderr,
                    "-E- %s line %d: Cannot interpolate type code: %d\n",
                    __FILE__, __LINE__, s.sds.ntype);
            exit(EXIT_FAILURE);
            break;
        }
    } // alloc_old

    /* Set interp function pointer according to SDS */
    if (strncmp(s.sds.name, "EV_", 3) == 0)
        interp_modis_var = interp_modis_EVdata;
    else if (strcmp(s.sds.name, "Longitude") == 0)
        interp_modis_var = interp_modis_Longitude;
    else
        interp_modis_var = interp_bilin;

    /* Calculate offset between different pixel resolutions */
    factor = (double) newres / (double) oldres;
    xoffset = 0; // frame centers are aligned
    yoffset = 0.5 * (1.0 - factor); // detector edges are aligned
    /* yoffset = 0;    // for testing only */

    /* Size of new array */
    nxnew = nx * oldres / newres;
    nynew = ny * oldres / newres;

    /* Calculate index into original data */
    for (iynew = 0; iynew < nynew; iynew++) {
        fy = factor * iynew - yoffset;
        y0 = MIN(MAX((int32_t) floor(fy), 0), ny - maxi);
        y1 = MIN(y0 + 1, ny - 1);
        dy = fy - (double) y0;

        for (ixnew = 0; ixnew < nxnew; ixnew++) {
            fx = factor * ixnew - xoffset;
            x0 = MIN(MAX((int32_t) floor(fx), 0), nx - maxi);
            x1 = MIN(x0 + 1, nx - 1);
            dx = fx - (double) x0;

            /* Do bilinear interpolation */
            val[0] = olddata[nx * y0 + x0];
            val[1] = olddata[nx * y0 + x1];
            val[2] = olddata[nx * y1 + x0];
            val[3] = olddata[nx * y1 + x1];
            newdata[nxnew * iynew + ixnew] = interp_modis_var(val, dx, dy);
        }
    }

    if (alloc_old)
        free(olddata);

    return SUCCESS;
}

/***********************************************************************/

/**
 * load values for one scan
 * @param[in] iscan
 * @return
 */
int load_modis_scan(const int32_t iscan,
        const int32_t resolution,
        int32_t sensorID) {
    modis_sds s;
    int32_t isds, iband;
    int32_t i, i1, i2, nold;
    static int32_t nvals = 0; /* pixels per variable at output resolution */
    double *olddata = NULL; /* at native resolution */
    double *newdata = NULL; /* at output resolution */

    /* Allocate space to hold all data for one scan */
    if (nvals == 0) {
        nvals = l1b[RSB_250].ndets * l1b[RSB_250].nframes * l1b[RSB_250].nsf;
        alloc_scan(nvals);
    }

    /* Allocate space to hold interpolation buffer for one variable */
    TRYMEM(__FILE__, __LINE__, (newdata = malloc(nvals * sizeof (*newdata))));

    /* ----- Geolocation values ----- */
    for (isds = 0; isds < GEO_NUM_SDS; isds++)
        read_sds_1scan(&geo[isds].sds,
            iscan,
            geo[isds].scandim,
            geo[isds].ndets);

    if (geo[GEO_LON].resolution == resolution) { /* native resolution */
        /* memcpy works only when src and dst have same data type */
        memcpy(scan.lon, geo[GEO_LON].sds.data, nvals * sizeof (float));
        memcpy(scan.lat, geo[GEO_LAT].sds.data, nvals * sizeof (float));

        /* must copy each element when datatypes differ */
        LOOPI{
            scan.hgt[i] = (float) ((int16_t *) geo[GEO_HGT].sds.data)[i];
            scan.solz[i] = (float) ((int16_t *) geo[GEO_SOLZ].sds.data)[i] * 0.01;
            scan.sola[i] = (float) ((int16_t *) geo[GEO_SOLA].sds.data)[i] * 0.01;
            scan.senz[i] = (float) ((int16_t *) geo[GEO_SENZ].sds.data)[i] * 0.01;
            scan.sena[i] = (float) ((int16_t *) geo[GEO_SENA].sds.data)[i] * 0.01;}

    } else { /* must interpolate */
        modis_interp(geo[GEO_LON], 0, resolution, NULL, newdata);
        LOOPI scan.lon[i] = (float) newdata[i];
        modis_interp(geo[GEO_LAT], 0, resolution, NULL, newdata);
        LOOPI scan.lat[i] = (float) newdata[i];
        modis_interp(geo[GEO_HGT], 0, resolution, NULL, newdata);
        LOOPI scan.hgt[i] = (float) newdata[i];
        modis_interp(geo[GEO_SOLZ], 0, resolution, NULL, newdata);
        LOOPI scan.solz[i] = (float) newdata[i] * 0.01;
        modis_interp(geo[GEO_SOLA], 0, resolution, NULL, newdata);
        LOOPI scan.sola[i] = (float) newdata[i] * 0.01;
        modis_interp(geo[GEO_SENZ], 0, resolution, NULL, newdata);
        LOOPI scan.senz[i] = (float) newdata[i] * 0.01;
        modis_interp(geo[GEO_SENA], 0, resolution, NULL, newdata);
        LOOPI scan.sena[i] = (float) newdata[i] * 0.01;
    }

    /* ----- Radiometric values ----- */
    for (isds = 0; isds < L1B_NUM_SDS; isds++)
        read_sds_1scan(&l1b[isds].sds,
            iscan,
            l1b[isds].scandim,
            l1b[isds].ndets);

    for (iband = 0; iband < scan.nbands; iband++) {
        s = l1b[BANDLIST[iband].sdsindex];
        nold = s.ndets * s.nframes * s.nsf;
        i1 = nold * BANDLIST[iband].bandindex; // offset within SDS
        i2 = nvals * iband; // offset within scan.allbands

        /* Cast SI values to double precision, in native spatial resolution */
        TRYMEM(__FILE__, __LINE__, (olddata = malloc(nold * sizeof (*olddata))));
        for (i = 0; i < nold; i++)
            olddata[i] = (double) ((uint16_t*) s.sds.data)[i1 + i];

        /* Apply subframe correction to destripe BEFORE interpolation */
        if ((sensorID == MODISA) && (resolution < 1000))
            subframe_correction(sensorID, iband, s, olddata);

        /* Fill any dead detectors with neighboring values */
        fill_dead_detectors(iband, s, olddata);

        /* Interpolate to higher resolution, as needed */
        if (s.resolution != resolution) {
            modis_interp(l1b[BANDLIST[iband].sdsindex],
                    BANDLIST[iband].bandindex,
                    resolution, olddata, newdata);

            /* Copy adjusted data into scan structure */
            LOOPI scan.allbands[i2++] = (float) newdata[i]; // interpolated
        } else
            LOOPI scan.allbands[i2++] = (float) olddata[i]; // native resolution

        free(olddata);
    }

    /* ----- Store preloaded values for later ----- */
    scan.iscan = iscan;
    scan.taisec = ((double*) ref[GEO_TAISEC].sds.data)[iscan];
    scan.mside = (int32_t) ((int16_t*) ref[GEO_MSIDE].sds.data)[iscan];
    for (i = 0; i < 3; i++) {
        scan.angles[i] = ((double*) ref[GEO_ANGLES].sds.data)[3 * iscan + i];
        scan.mnorm[i] = ((double*) ref[GEO_MNORM].sds.data)[3 * (3 * iscan + i)];
    }
    free(newdata);
    return SUCCESS;
}

/**
 * set all scan-based values in l1rec
 * @param[in] iscan
 * @param[in,out] l1rec
 * @return
 */
int set_l1rec_scanvals(const int32_t iscan, l1str *l1rec) {
    double esdist;
    static double lastTAIsec = -1;
    int16_t badatt, badmside;
    double TAIsec, dsec;
    int16_t year, day;
    int32_t i;

    /* Mirror Side */
    l1rec->mside = scan.mside;
    badmside = (l1rec->mside < 0);
    if (badmside) /* Fix mirror side if set to -1 */
        l1rec->mside = 0;

    /* Check for non-nominal roll, pitch, or yaw */
    badatt =
            (fabs(scan.angles[0]) > MAX_ATTERR) ||
            (fabs(scan.angles[1]) > MAX_ATTERR) ||
            (fabs(scan.angles[2]) > MAX_ATTERR);

    /* Flag each pixel of scan for geolocation errors */
    for (i = 0; i < l1rec->npix; i++) {
        l1rec->pixnum[i] = spixl + i;
        l1rec->navfail[i] = badmside;
        l1rec->navwarn[i] = badatt;
    }

    /* Scan Start Time */
    TAIsec = scan.taisec;
    if (TAIsec < 0) { /* Workaround for scan time errors */
        fprintf(stderr, "-W- %s: Bad time for scan %d; ", __FILE__, iscan);
        if (iscan == 0) {
            fprintf(stderr, "using granule start time.\n");
            TAIsec = unix_to_tai93(file_geo.starttime);
        } else {
            fprintf(stderr, "extrapolating from last good scan.\n");
            TAIsec = lastTAIsec + SCAN_TIME_INTERVAL;
        }
    }
    l1rec->scantime = tai93_to_unix(TAIsec);
    lastTAIsec = TAIsec;

    /* Earth-sun distance correction for this scan */
    unix2yds(l1rec->scantime, &year, &day, &dsec);
    int32_t msec = (int32_t) (dsec * 1000.0);
    int32_t yr = (int32_t) year;
    int32_t dy = (int32_t) day;
    esdist = esdist_(&yr, &dy, &msec);
    l1rec->fsol = pow(1.0 / esdist, 2);

    return SUCCESS;
}

/**
 *
 * @param[in,out] l1file
 * @param[in] line
 * @param[in,out] l1rec
 * @return
 */
int readl1_modis(filehandle *l1file,
        const int32_t line,
        l1str *l1rec, int lonlat) {
    static int32_t prev_scan = -1;
    int32_t iband, iscan, idet, ndets;
    size_t ip, ipb;
    size_t detoffset, dataptr;
    size_t nbytes;
    size_t irsb, iteb, ical;
    att_struct *atts;
    float scale, offset, SI_val, value;

    /* Scan and detector numbers for this line */
    ndets = l1b[RSB_250].ndets;
    iscan = line / ndets;
    idet = line % ndets;

    /* Load all values for next scan as needed */
    if (iscan != prev_scan) {
        prev_scan = iscan;
        load_modis_scan(iscan, l1_input->resolution, l1file->sensorID);
    }

    /* The l1rec structure is reinitialized */
    /* each line, so always reload basic info. */
    l1rec->npix = l1file->npix;
    l1rec->detnum = (int32_t) idet;
    set_l1rec_scanvals(iscan, l1rec);
    unix_to_tai93(file_geo.starttime);

    /* "scan" structure data pointer offset due */
    /* to line number and extract start pixel */
    detoffset = idet * scan.npix + spixl;

    /* ----- Geolocation values ----- */
    nbytes = l1rec->npix * sizeof (float);
    memcpy(l1rec->lon, &scan.lon[detoffset], nbytes);
    memcpy(l1rec->lat, &scan.lat[detoffset], nbytes);
    memcpy(l1rec->solz, &scan.solz[detoffset], nbytes);
    
    // lonlat only needs time, lon, lat and solar zenith
    if (lonlat)
        return LIFE_IS_GOOD;
        
    memcpy(l1rec->height, &scan.hgt[detoffset], nbytes);
    memcpy(l1rec->solz, &scan.solz[detoffset], nbytes);
    memcpy(l1rec->sola, &scan.sola[detoffset], nbytes);
    memcpy(l1rec->senz, &scan.senz[detoffset], nbytes);
    memcpy(l1rec->sena, &scan.sena[detoffset], nbytes);

    /* Don't need to check for navigation errors
       - it's done by readl1() */

    /* Compute polarization frame rotation angles */
    compute_alpha(l1rec->lon, l1rec->lat,
            l1rec->senz, l1rec->sena,
            scan.mnorm, l1rec->npix, l1rec->alpha);

    /* ----- Radiometric values ----- */
    irsb = iteb = 0;
    for (iband = 0; iband < scan.nbands; iband++) {

        /* calibration info */
        atts = l1b[BANDLIST[iband].sdsindex].sds.atts;
        ical = BANDLIST[iband].bandindex;

        /* offset within scan.allbands */
        dataptr = iband * scan.nvals + detoffset;

        /*** Thermal Emissive Bands ***/
        if (BANDLIST[iband].sdsindex == TEB_1KM) {
            scale = ((float*) atts[RAD_SCALE].data)[ical];
            offset = ((float*) atts[RAD_OFFSET].data)[ical];

            for (ip = 0; ip < l1rec->npix; ip++) {
                value = 0.0; /* init to fill value */

                /* Apply calibration for valid values */
                SI_val = scan.allbands[dataptr + ip];
                if (SI_val < MIN_BAD_SI) {

                    /* Convert radiance from Watts/m^2/micrometer/steradian */
                    /* to mW/cm2/um/sr, as expected by l2gen */
                    value = (SI_val - offset) * scale / 10.0;

                    /* Terra adjustments */
                    if (l1rec->l1file->sensorID == MODIST) {
                        value *= (1.0 - radoff[iteb][idet / (ndets / 10)]);
                        if (l1rec->mside == 1)
                            value *= (1.0 + mfact[iteb]);
                    }                        /* Aqua adjustments */
                    else {
                        /* channel 22 detector balancing */
                        if (ical == 1)
                            value /= l1_input->ch22detcor[line % 10];
                        /* channel 23 detector balancing */
                        if (ical == 2)
                            value /= l1_input->ch23detcor[line % 10];
                    }
                }

                /* populate l1rec->Ltir */
                ipb = ip * l1rec->l1file->nbandsir + iteb; /* band varies fastest */
                l1rec->Ltir[ipb] = value;

            } // ip
            iteb++;
        }// TEB

            /*** Cirrus Band ***/
        else if (BANDLIST[iband].sdsindex == CIR_1KM) {
            scale = ((float*) atts[REFL_SCALE].data)[ical];
            offset = ((float*) atts[REFL_OFFSET].data)[ical];

            for (ip = 0; ip < l1rec->npix; ip++) {
                value = BAD_FLT; /* init to fill value */

                /* Apply calibration for valid values */
                SI_val = scan.allbands[dataptr + ip];
                if (SI_val < MIN_BAD_SI)

                    /* Normalize reflectance by solar zenith angle */
                    value = (SI_val - offset) * scale
                        / cos(l1rec->solz[ip] / OEL_RADEG);

                /* populate l1rec->rho_cirrus */
                l1rec->rho_cirrus[ip] = value;

            } // ip
        }// CIR

            /*** Reflective Solar Bands ***/
        else {
            l1rec->Fo[irsb] = Fobar[irsb] * l1rec->fsol;
            scale = ((float*) atts[REFL_SCALE].data)[ical];
            offset = ((float*) atts[REFL_OFFSET].data)[ical];

            for (ip = 0; ip < l1rec->npix; ip++) {
                value = BAD_FLT; /* init to fill value */

                /* skip night data for visible bands */
                if (SOLZNIGHT < l1rec->solz[ip])
                    continue;

                /* Apply calibration for valid values */
                SI_val = scan.allbands[dataptr + ip];
                if (SI_val < MIN_BAD_SI)

                    /* Convert from reflectance to radiance */
                    value = (SI_val - offset) * scale * l1rec->Fo[irsb] / OEL_PI;

                    /* Flag any sentinel values */
                else
                    switch ((int32_t) SI_val) {
                    case 65527: /* not in earth view: SECTOR_ROTATION_SI */
                        l1rec->navfail[ip] = 1;
                        break;
                    case 65529: /* saturation: TEB_OR_RSB_GT_MAX_SI  */
                    case 65533: /* saturation: SATURATED_DETECTOR_SI */
                        value = 1000.0;
                        if (irsb < 13)
                            l1rec->hilt[ip] = 1;
                        break;
                    default: /* all others */
                        break;
                    }

                /* populate l1rec->Lt */
                ipb = ip * l1rec->l1file->nbands + irsb; /* band varies fastest */
                l1rec->Lt[ipb] = value;

            } // ip
            irsb++;
        } // RSB

    } // iband

    /* Convert IR bands to brightness temperature */
    radiance2bt(l1rec, l1_input->resolution);

    return SUCCESS;
}

/***********************************************************************/

/**
 *
 * @return
 */
int closel1_modis() {
    int32_t i;

    /* Free memory allocated for global variables */
    for (i = 0; i < REF_NUM_SDS; i++)
        free_sds(&(ref[i].sds));
    for (i = 0; i < GEO_NUM_SDS; i++)
        free_sds(&(geo[i].sds));
    for (i = 0; i < L1B_NUM_SDS; i++)
        free_sds(&(l1b[i].sds));

    /* Free any memory allocated for input files */
    if (file_geo.sds) {
        for (i = 0; i < file_geo.n_sds; i++)
            free_sds(&(file_geo.sds[i]));
        free(file_geo.sds);
    }
    if (file_1km.sds) {
        for (i = 0; i < file_1km.n_sds; i++)
            free_sds(&(file_1km.sds[i]));
        free(file_1km.sds);
    }
    if (file_hkm.sds) {
        for (i = 0; i < file_hkm.n_sds; i++)
            free_sds(&(file_hkm.sds[i]));
        free(file_hkm.sds);
    }
    if (file_qkm.sds) {
        for (i = 0; i < file_qkm.n_sds; i++)
            free_sds(&(file_qkm.sds[i]));
        free(file_qkm.sds);
    }

    /* End HDF access */
    if (file_geo.id) SDend(file_geo.id);
    if (file_1km.id) SDend(file_1km.id);
    if (file_hkm.id) SDend(file_hkm.id);
    if (file_qkm.id) SDend(file_qkm.id);
    file_geo.id = file_1km.id = file_hkm.id = file_qkm.id = 0;

    /* TO DO: free scan structure */

    return SUCCESS;
}

/***********************************************************************/

int readl1_lonlat_modis(filehandle *l1file, int32_t line, l1str *l1rec) {
    static float *lonscan = NULL; /* storage buffer arrays */
    static float *latscan = NULL; /* (at output resolution) */

    int32_t ndets = l1b[RSB_250].ndets; /* output resolution dimensions*/
    int32_t npixl = l1b[RSB_250].nframes * l1b[RSB_250].nsf;
    int32_t nvals = ndets * npixl;
    int32_t resolution = l1_input->resolution;

    double *newdata = NULL; /* interpolation buffer */
    static int32_t prev_scan = -1;
    int32_t i, iscan, idet;
    size_t detoffset, nbytes;

    /* Allocate space to hold lon/lat data for one scan at output resolution */
    if (lonscan == NULL) {
        TRYMEM(__FILE__, __LINE__, (lonscan = malloc(nvals * sizeof (*lonscan))));
        TRYMEM(__FILE__, __LINE__, (latscan = malloc(nvals * sizeof (*latscan))));
    }

    /* Scan and detector indices for this line */
    iscan = line / ndets;
    idet = line % ndets;

    /* Load next scan as needed */
    if (iscan != prev_scan) {
        prev_scan = iscan;

        /* Read values */
        read_sds_1scan(&geo[GEO_LON].sds,
                iscan,
                geo[GEO_LON].scandim,
                geo[GEO_LON].ndets);
        read_sds_1scan(&geo[GEO_LAT].sds,
                iscan,
                geo[GEO_LAT].scandim,
                geo[GEO_LAT].ndets);

        /* If native resolution, copy values into storage buffers */
        if (geo[GEO_LON].resolution == resolution) {
            memcpy(lonscan, geo[GEO_LON].sds.data, nvals * sizeof (*lonscan));
            memcpy(latscan, geo[GEO_LAT].sds.data, nvals * sizeof (*latscan));
        }
            /* Otherwise, interpolate first  */
        else {
            TRYMEM(__FILE__, __LINE__,
                    (newdata = malloc(nvals * sizeof (*newdata))));
            modis_interp(geo[GEO_LON], 0, resolution, NULL, newdata);
            LOOPI lonscan[i] = (float) newdata[i];
            modis_interp(geo[GEO_LAT], 0, resolution, NULL, newdata);
            LOOPI latscan[i] = (float) newdata[i];
            free(newdata);
        }

    } // loaded new scan

    /* Copy specified line from storage buffers into l1rec structure */
    detoffset = idet * npixl + spixl;
    nbytes = l1rec->npix * sizeof (float);
    memcpy(l1rec->lon, &lonscan[detoffset], nbytes);
    memcpy(l1rec->lat, &latscan[detoffset], nbytes);

    return SUCCESS;
}
