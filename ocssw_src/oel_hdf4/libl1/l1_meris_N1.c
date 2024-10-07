/**
 *  @file l1_meris_N1.h
 *  @brief MERIS reader
 *  @author Paul Lyon
 *  @author Naval Research Laboratory, Stennis Space Center, MS
 */

/* ============================================================================ */
/* module l1_meris_N1.c - functions to read MERIS   Reduced RES L2 for MSL12    */
/* NOTE!! THIS IS A TEMPORARY READER FOR L2 FILES, INTO L1 SEADAS RECORD        */
/* Written By:  Paul E. Lyon NRL, Oct. 2006.                                    */
/*                                                                              */
/* ============================================================================ */

#include "l1_meris_N1.h"
#include "epr_api.h"
#include "epr_field.h"
#include <math.h>
#include "smile.h"
#include <libnav.h>

#include "l1.h"

#include <stdbool.h>

#define MERIS_NBANDS 15

#define MERIS_BANDINFO_FILENAME      "band_info_meris.txt"

#define MERIS_WAVELENGTH_FR_FILENAME "central_wavelen_fr.txt"
#define MERIS_WAVELENGTH_RR_FILENAME "central_wavelen_rr.txt"
#define MERIS_SUN_FLUX_FR_FILENAME   "sun_spectral_flux_fr.txt"
#define MERIS_SUN_FLUX_RR_FILENAME   "sun_spectral_flux_rr.txt"

// # of detectors in a full resolution scan line
#define MERIS_FR_DETECTORS 3700
// # of detectors in a full resolution scan line
#define MERIS_RR_DETECTORS 925

// max line length in a config table file
#define MERIS_LINE_MAX 1024

/* ------------------------------------------------------------------------ */
/* These L1 flag values were derrived from the MERIS specification          */
/*                                                                          */
/*  http://earth.esa.int/pub/ESA_DOC/ENVISAT/Vol11_Meris_5b.pdf             */
/*  page 64,  Table 11.4.1.7.4.2-3 MER_RR__1P - Quality Flag Coding         */
/*                                                                          */
/* dshea                                                                    */
/* ------------------------------------------------------------------------ */
#define MERIS_L1FLAG_COSMETIC   0x01
#define MERIS_L1FLAG_DUPLICATED 0x02
#define MERIS_L1FLAG_GLINT      0x04
#define MERIS_L1FLAG_SUSPECT    0x08
#define MERIS_L1FLAG_LAND       0x10
#define MERIS_L1FLAG_BRIGHT     0x20
#define MERIS_L1FLAG_COASTLINE  0x40
#define MERIS_L1FLAG_INVALID    0x80

/* ------------------------------------------------------------------------ */
/* These L2 flag values were derrived from the MERIS specification          */
/*                                                                          */
/*  http://earth.esa.int/pub/ESA_DOC/ENVISAT/Vol11_Meris_5b.pdf             */
/*  page 105,  Table 11.5.1.7.4.8-2 Description of the Flag Coding          */
/*                                                                          */
/* dshea                                                                    */
/* ------------------------------------------------------------------------ */
#define MERIS_L2FLAG_WHITE_SCATTER   0x000001
#define MERIS_L2FLAG_PRESSURE_CONF   0x000002
#define MERIS_L2FLAG_HIGH_GLINT      0x000004
#define MERIS_L2FLAG_DDV             0x000008
#define MERIS_L2FLAG_MEDIUM_GLINT    0x000010
#define MERIS_L2FLAG_ICE_HAZE        0x000020
#define MERIS_L2FLAG_CASE2_Y         0x000040
#define MERIS_L2FLAG_CASE2_ANOM      0x000080
#define MERIS_L2FLAG_CASE2_S         0x000100
#define MERIS_L2FLAG_ABSOA_DUST      0x000200
#define MERIS_L2FLAG_OOADB           0x000400
#define MERIS_L2FLAG_SUSPECT         0x000800
#define MERIS_L2FLAG_COSMETIC        0x001000
#define MERIS_L2FLAG_COASTLINE       0x002000
#define MERIS_L2FLAG_PCD_19          0x004000
#define MERIS_L2FLAG_PCD_18          0x008000
#define MERIS_L2FLAG_PCD_17          0x010000
#define MERIS_L2FLAG_PCD_16          0x020000
#define MERIS_L2FLAG_PCD_15          0x040000
#define MERIS_L2FLAG_PCD_14          0x080000
#define MERIS_L2FLAG_PCD_1_13        0x100000
#define MERIS_L2FLAG_WATER           0x200000
#define MERIS_L2FLAG_CLOUD           0x400000
#define MERIS_L2FLAG_LAND            0x800000

static EPR_SProductId *fileID = NULL;
static int spix = 0;
static int file_npix;
static double fileStartTime; // first scan time in unix seconds
static double time_interval; // scan interval in seconds

static int fullResolution = 0;
static int numDetectors;
static float *detectorWL;
static float *detectorE0;

/**
 * @brief opens a MERIS file for reading to load into L1 record
 * @param[in]   file   file handle to MERIS file
 *
 * @author Paul E. Lyon NRL, Oct. 2006.
 */

int
openl1_meris_N1(filehandle * file) {
    const char *fltime;
    const char *fname;
    char monthstr[10];
    static char months_list[12][4] ={"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL",
        "AUG", "SEP", "OCT", "NOV", "DEC"};
    unsigned int source_w, source_h;
    int minute, hour, month, i;
    float fsec;
    const EPR_SRecord *record;
    const EPR_SField *field;
    EPR_SBandId *band_id = NULL;
    char *sph_names[] ={"NUM_BANDS", "LINE_TIME_INTERVAL", "NUM_SLICES"};
    int nf = 3;

    // Initialize the API
    epr_init_api(e_log_debug, NULL, NULL);

    // Open the N1 input file 

    fileID = epr_open_product(file->name);
    if (fileID == NULL) {
        fprintf(stderr, "-E- %s line %d: epr_open_product(%s) failed.\n",
                __FILE__, __LINE__, file->name);
        return (1);
    }

    band_id = epr_get_band_id(fileID, "l1_flags");
    if (band_id == NULL) {
        fprintf(stderr, "-E- %s line %d - Can not find \"l1_flags\" in Meris N1 file\n",
                __FILE__, __LINE__);
        return (1);
    }

    // Get pixel and scan dimensions 

    source_h = fileID->scene_height;
    source_w = fileID->scene_width;
    file_npix = source_w;
    if (want_verbose) {
        printf("MERIS Level1B Npix  :%u Nscans:%u\n", source_w, source_h);
    } // want_verbose

    // get specific product header (SPH)

    record = epr_get_sph(fileID);
    if (record == NULL) {
        fprintf(stderr, "-E- %s line %d: epr_get_sph(fileID) failed.\n",
                __FILE__, __LINE__);
        return (1);
    }


    // get first scan time 

    field = epr_get_field(record, "FIRST_LINE_TIME");
    if (field == NULL) {
        fprintf(stderr,
                "-E- %s line %d: epr_get_field(record,FIRST_LINE_TIME) failed.\n",
                __FILE__, __LINE__);
        return (1);
    }

    int year, day;
    fltime = epr_get_field_elem_as_str(field);
    sscanf(fltime, "%02d-%3s-%04d %02d:%02d:%f", &day, monthstr, &year,
            &hour, &minute, &fsec);
    monthstr[4] = '\0';
    month = 1;
    for (i = 0; i < 12; i++)
        if (strncmp(monthstr, months_list[i], 3) == 0)
            month = i + 1;
    fileStartTime = ymds2unix(year, month, day, hour * 3600.0 + minute * 60.0 + fsec);

    field = epr_get_field(record, "LINE_TIME_INTERVAL");
    if (field == NULL) {
        fprintf(stderr,
                "-E- %s line %d: epr_get_field(record,FIRST_LINE_TIME) failed.\n",
                __FILE__, __LINE__);
        return (1);
    }
    time_interval = (double) epr_get_field_elem_as_uint(field, 0); // microseconds;
    time_interval /= 1e6; // change to seconds

    if (want_verbose)
        printf("MERIS time interval %lf\n", time_interval);

    // dump a few more items for infomational purposes

    for (i = 0; i < nf; i++) {
        field = epr_get_field(record, sph_names[i]);
        if (field == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: epr_get_field(record,%s) failed.\n",
                    __FILE__, __LINE__, sph_names[i]);
            return (1);
        }
        if (want_verbose)
            printf("\t%s = %d\n", sph_names[i],
                (int) epr_get_field_elem_as_uint(field, 0));
    }

    // get the START_PIXEL and END_PIXEL that were injected into the MPH
    // by the extractor.  START_PIXEL and END_PIXEL are 1 based.  The epr
    // api reverses the pixels on a line, so spix = source_w - epix.
    record = epr_get_mph(fileID);
    if (record) {
        field = epr_get_field(record, "PRODUCT");
        if (field) {
            fname = epr_get_field_elem_as_str(field);
            if (strncmp(fname, "MER_RR", 6) == 0) {
                fullResolution = 0;
                numDetectors = MERIS_RR_DETECTORS;
                strcpy(file->spatialResolution, "1.2 km");
            } else {
                fullResolution = 1;
                numDetectors = MERIS_FR_DETECTORS;
                strcpy(file->spatialResolution, "300 m");
            }
        } else {
            printf("-E- PRODUCT field not found in the MPH header.");
            exit(1);
        }

        field = epr_get_field(record, "END_PIXEL");
        if (field) {
            int tmp_epix = epr_get_field_elem_as_uint(field, 0);
            field = epr_get_field(record, "START_PIXEL");
            if (field) {
                int tmp_spix = epr_get_field_elem_as_uint(field, 0);
                if ((tmp_spix < source_w) && (tmp_spix < tmp_epix)) {
                    spix = source_w - tmp_epix;
                    source_w = tmp_epix - tmp_spix + 1; // spix and epix are inclusive
                    if (want_verbose)
                        printf("OBPG Extract - spix=%d, npix=%d\n",
                            spix + 1, (int) source_w);
                } // spix, epix reasonable
            } // spix good
        } // epix good
    } else { // found the MPH
        printf("-E- MPH header not found.");
        exit(1);
    }

    // define number of input products

    file->nbands = MERIS_NBANDS;
    file->npix = (int) source_w;
    file->nscan = (int) source_h;

    return (LIFE_IS_GOOD);
}

/**
 * @brief reads 1 scan line from MERIS file, loads l1rec
 * @param[in]   file   file handle to MERIS file
 * @param[in]   scan   scan number to read
 * @apram[out]  l1rec  output l1rec
 *
 * @author Paul E. Lyon NRL, Oct. 2006.
 * W. Robinson, SAIC, 22 May 2012  account for msec going to next day
 */

int readl1_meris_N1(filehandle *file, int32_t scan, l1str *l1rec) {
    static int firstCall = 1;
    int err_code;

    static char *names[] = {
        "radiance_1", "radiance_2", "radiance_3", "radiance_4", // 412.7, 442.6, 489.9, 509.8
        "radiance_5", "radiance_6", "radiance_7", "radiance_8", // 559.7, 619.6, 664.5, 680.8
        "radiance_9", "radiance_10", "radiance_11", "radiance_12", // 708.3, 753.4, 761.5, 778.4
        "radiance_13", "radiance_14", "radiance_15" // 864.9, 884.9, 900.0
    };
    int32_t npix;
    int32_t nbands;
    int32_t ip, ib, ipb;
    EPR_SBandId *band_id = NULL;
    EPR_SRaster *temp_raster;
    epr_uint flag;
    static char *invalid_flag;
    static char *land_flag;


    if (firstCall) {
        firstCall = 0;
        invalid_flag = malloc(sizeof (char) * file->npix);
        land_flag = malloc(sizeof (char) * file->npix);

        if (l1_input->rad_opt != 0) {
            // setup the smile correction
            char *tmp_str;
            char line[MERIS_LINE_MAX];
            char path[FILENAME_MAX] = "";
            char filename[FILENAME_MAX] = "";
            FILE *fp = NULL;
            int count;
            int dummy;
            int detector;
            float *floatp;

            if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
                printf("OCDATAROOT environment variable is not defined.\n");
                exit(1);
            }

            strcpy(path, tmp_str);
            strcat(path, "/");
            strcat(path, sensorId2SensorDir(l1rec->l1file->sensorID));
            strcat(path, "/cal/");

            /*-------------------------------------------------------------*/
            detectorWL = (float*) malloc(sizeof (float)*MERIS_NBANDS * numDetectors);
            if (!detectorWL) {
                printf("-E- %s line %d : error allocating memory for detectorWL.\n",
                        __FILE__, __LINE__);
                exit(1);
            }

            /*-------------------------------------------------------------*/
            // read the central wavelength table
            strcpy(filename, path);
            if (fullResolution)
                strcat(filename, MERIS_WAVELENGTH_FR_FILENAME);
            else
                strcat(filename, MERIS_WAVELENGTH_RR_FILENAME);
            if ((fp = fopen(filename, "r")) == NULL) {
                fprintf(stderr,
                        "-E- %s line %d: unable to open %s for reading\n", __FILE__, __LINE__, filename);
                exit(1);
            }

            // discard the first line of column labels
            fgets(line, MERIS_LINE_MAX, fp);

            floatp = detectorWL;
            for (detector = 0; detector < numDetectors; detector++) {
                if (!fgets(line, MERIS_LINE_MAX, fp)) {
                    fprintf(stderr,
                            "-E- %s line %d: unable to read detector %d from file %s\n", __FILE__, __LINE__, detector, filename);
                    exit(1);
                }

                /* I really should read all of the bands defined by MERIS_NBANDS, but I am just reading
             15 bands for now.  dshea */
                count = sscanf(line, "%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", &dummy, floatp, floatp + 1,
                        floatp + 2, floatp + 3, floatp + 4, floatp + 5, floatp + 6, floatp + 7, floatp + 8, floatp + 9,
                        floatp + 10, floatp + 11, floatp + 12, floatp + 13, floatp + 14);
                if (count != 16) {
                    fprintf(stderr,
                            "-E- %s line %d: unable to read whole detector %d line from file %s, count = %d\n",
                            __FILE__, __LINE__, detector, filename, count);
                    exit(1);
                }

                floatp += MERIS_NBANDS;
            } // for detector
            fclose(fp);


            /*-------------------------------------------------------------*/
            // read solar flux table
            detectorE0 = (float*) malloc(sizeof (float)*MERIS_NBANDS * numDetectors);
            if (!detectorE0) {
                printf("-E- %s line %d : error allocating memory for detectorE0.\n",
                        __FILE__, __LINE__);
                exit(1);
            }

            strcpy(filename, path);
            if (fullResolution)
                strcat(filename, MERIS_SUN_FLUX_FR_FILENAME);
            else
                strcat(filename, MERIS_SUN_FLUX_RR_FILENAME);
            if ((fp = fopen(filename, "r")) == NULL) {
                fprintf(stderr,
                        "-E- %s line %d: unable to open %s for reading\n", __FILE__, __LINE__, filename);
                exit(1);
            }

            // discard the first line of column labels
            fgets(line, MERIS_LINE_MAX, fp);

            floatp = detectorE0;
            for (detector = 0; detector < numDetectors; detector++) {
                if (!fgets(line, MERIS_LINE_MAX, fp)) {
                    fprintf(stderr,
                            "-E- %s line %d: unable to read detector %d from file %s\n", __FILE__, __LINE__, detector, filename);
                    exit(1);
                }

                /* I really should read all of the bands defined by MERIS_NBANDS, but I am just reading
             15 bands for now.  dshea */
                count = sscanf(line, "%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", &dummy, floatp, floatp + 1,
                        floatp + 2, floatp + 3, floatp + 4, floatp + 5, floatp + 6, floatp + 7, floatp + 8, floatp + 9,
                        floatp + 10, floatp + 11, floatp + 12, floatp + 13, floatp + 14);
                if (count != 16) {
                    fprintf(stderr,
                            "-E- %s line %d: unable to read whole detector %d line from file %s, count = %d\n",
                            __FILE__, __LINE__, detector, filename, count);
                    exit(1);
                }

                floatp += MERIS_NBANDS;
            } // for detector
            fclose(fp);

            strcpy(filename, path);
            strcat(filename, MERIS_BANDINFO_FILENAME);
            smile_init(MERIS_NBANDS, numDetectors, filename, detectorWL, detectorE0);

        } // if input->rad_opt

    } // first call

    nbands = file->nbands;
    npix = file->npix;

    l1rec->scantime = fileStartTime + time_interval*scan;
    int16_t syear, sday;
    double secs;
    unix2yds(l1rec->scantime, &syear, &sday, &secs);

    int32_t yr = syear;
    int32_t dy = sday;
    int32_t msec = (int32_t) (secs * 1000.0);
    double esdist = esdist_(&yr, &dy, &msec);

    l1rec->fsol = pow(1.0 / esdist, 2);
    // read L1 flags and set invalid_flag if suspect or invalid set

    band_id = epr_get_band_id(fileID, "l1_flags");
    if (band_id == NULL) {
        printf("-E- %s line %d: Error finding band_id:l1_flags\n", __FILE__,
                __LINE__);
        return (1);
    }
    temp_raster = epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
    if (temp_raster == NULL) {
        printf("-E- %s line %d: Error allocating raster space.\n", __FILE__,
                __LINE__);
        return (1);
    }
    err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
    for (ip = 0; ip < npix; ip++) {
        flag = epr_get_pixel_as_uint(temp_raster, spix + ip, 0);
        if ((flag & MERIS_L1FLAG_SUSPECT) || (flag & MERIS_L1FLAG_INVALID)) {
            invalid_flag[ip] = 1;
        } else {
            invalid_flag[ip] = 0;
        }
        if (flag & MERIS_L1FLAG_LAND){
            land_flag[ip] = 1;
        } else{
            land_flag[ip] = 0;
        }
    }

    epr_free_raster(temp_raster);


    // read latitude

    band_id = epr_get_band_id(fileID, "latitude");
    if (band_id == NULL) {
        printf("-E- %s line %d: Error finding band_id:latitude\n", __FILE__,
                __LINE__);
        return (1);
    }
    temp_raster = epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
    if (temp_raster == NULL) {
        printf("-E- %s line %d: Error allocating raster space.\n", __FILE__,
                __LINE__);
        return (1);
    }
    err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
    for (ip = 0; ip < npix; ip++)
        l1rec->lat[ip] = epr_get_pixel_as_float(temp_raster, spix + ip, 0);

    // read longitude

    band_id = epr_get_band_id(fileID, "longitude");
    if (band_id == NULL) {
        printf("-E- %s line %d: Error finding band_id:longitude\n", __FILE__,
                __LINE__);
        return (1);
    }
    err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
    for (ip = 0; ip < npix; ip++)
        l1rec->lon[ip] = epr_get_pixel_as_float(temp_raster, spix + ip, 0);
    epr_free_raster(temp_raster);

    // read sun zenith

    band_id = epr_get_band_id(fileID, "sun_zenith");
    if (band_id == NULL) {
        printf("-E- %s line %d: Error finding band_id:sun_zenith", __FILE__,
                __LINE__);
        return (1);
    }
    temp_raster = epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
    if (temp_raster == NULL) {
        printf("-E- %s line %d: Error allocating raster space.\n", __FILE__,
                __LINE__);
        return (1);
    }
    err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
    for (ip = 0; ip < npix; ip++)
        l1rec->solz[ip] = epr_get_pixel_as_float(temp_raster, spix + ip, 0);
    epr_free_raster(temp_raster);

    // read sun azimuth

    band_id = epr_get_band_id(fileID, "sun_azimuth");
    if (band_id == NULL) {
        printf("-E- %s line %d: Error finding band_id:sun_azimuth", __FILE__,
                __LINE__);
        return (1);
    }
    temp_raster = epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
    if (temp_raster == NULL) {
        printf("-E- %s line %d: Error allocating raster space.\n", __FILE__,
                __LINE__);
        return (1);
    }
    err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
    for (ip = 0; ip < npix; ip++)
        l1rec->sola[ip] = epr_get_pixel_as_float(temp_raster, spix + ip, 0);
    epr_free_raster(temp_raster);

    // read view zenith

    band_id = epr_get_band_id(fileID, "view_zenith");
    if (band_id == NULL) {
        printf("-E- %s line %d: Error finding band_id:view_zenith", __FILE__,
                __LINE__);
        return (1);
    }
    temp_raster = epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
    if (temp_raster == NULL) {
        printf("-E- %s line %d: Error allocating raster space.\n", __FILE__,
                __LINE__);
        return (1);
    }
    err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
    for (ip = 0; ip < npix; ip++)
        l1rec->senz[ip] = epr_get_pixel_as_float(temp_raster, spix + ip, 0);
    epr_free_raster(temp_raster);

    // read view azimuth

    band_id = epr_get_band_id(fileID, "view_azimuth");
    if (band_id == NULL) {
        printf("-E- %s line %d: Error finding band_id:view_azimuth", __FILE__,
                __LINE__);
        return (1);
    }
    temp_raster = epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
    if (temp_raster == NULL) {
        printf("-E- %s line %d: Error allocating raster space.\n", __FILE__,
                __LINE__);
        return (1);
    }
    err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
    for (ip = 0; ip < npix; ip++)
        l1rec->sena[ip] = epr_get_pixel_as_float(temp_raster, spix + ip, 0);
    epr_free_raster(temp_raster);

    // set pixnum and check for navigation failure

    for (ip = 0; ip < npix; ip++) {
        l1rec->pixnum[ip] = spix + ip;
        if (l1rec->lon[ip] < -181.0 || l1rec->lon[ip] > 181.0 ||
                l1rec->lat[ip] < -91.0 || l1rec->lat[ip] > 91.0)
            l1rec->navfail[ip] = 1;
    }

    // read detector_index
    band_id = epr_get_band_id(fileID, "detector_index");
    if (band_id == NULL) {
        printf("-E- %s line %d: Error finding band_id:detector_index\n", __FILE__,
                __LINE__);
        return (1);
    }
    temp_raster = epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
    if (temp_raster == NULL) {
        printf("-E- %s line %d: Error allocating raster space.\n", __FILE__,
                __LINE__);
        return (1);
    }
    err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
    for (ip = 0; ip < npix; ip++) {
        l1rec->pixdet[ip] = epr_get_pixel_as_uint(temp_raster, spix + ip, 0);
    }
    epr_free_raster(temp_raster);


    // read in data

    for (ib = 0; ib < nbands; ib++) {

        band_id = epr_get_band_id(fileID, names[ib]);
        if (band_id == NULL) {
            printf("-E- %s line %d: Error finding band_id:%s\n", __FILE__,
                    __LINE__, names[ib]);
            return (1);

        } else {
            temp_raster = epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
            if (temp_raster == NULL) {
                printf("-E- %s line %d: Error allocating raster space.\n",
                        __FILE__, __LINE__);
                return (1);
            }
            err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
            if (!err_code) {

                // copy to Lt record.  Note that this might actually be Rrs
                // if the input file is a MERIS L2 file

                for (ip = 0; ip < npix; ip++) {

                    ipb = ip * nbands + ib;
                    l1rec->Lt[ipb] = epr_get_pixel_as_float(temp_raster, spix + ip, 0);

                    // set the value of all of the bands to BAD_FLT if
                    // navfail has been flagged.
                    if (invalid_flag[ip])
                        l1rec->Lt[ipb] = BAD_FLT;
                    else
                        l1rec->Lt[ipb] /= 10.0; // units conversion

                    // mark negative input data as HILT
                    if (l1rec->Lt[ipb] < 0.0)
                        l1rec->hilt[ip] = 1;

                }
            }
            epr_free_raster(temp_raster);

        } // if band_id found

    } // for ib
    // radcor needs all Lts populated before running so...one more time around the block...
    if (l1_input->rad_opt != 0) {
        for (ip = spix; ip < npix; ip++) {
            // es not corrected f0, so setting 0 as 4 element - seemed silly to make a variable for it
            radcor(l1rec, ip, land_flag[ip], 0);
            for (ib = 0; ib < nbands; ib++) {
                ipb = ip * nbands + ib;
                l1rec->Lt[ipb] += l1rec->radcor[ipb];
            }
        }
    }
    l1rec->npix = file->npix;

    return (LIFE_IS_GOOD);
}

/**
 * @brief reads 1 scan line from MERIS file, loads l1rec
 * @param[in]   file   file handle to MERIS file
 * @param[in]   scan   scan number to read
 * @apram[out]  l1rec  output l1rec
 *
 * @author Paul E. Lyon NRL, Oct. 2006.
 */

int readl1_lonlat_meris_N1(filehandle *file, int32_t scan, l1str *l1rec) {
    int32_t npix;
    int32_t ip;
    static EPR_SBandId *band_id;
    static EPR_SRaster *temp_raster;

    npix = file->npix;

    // Time
    l1rec->scantime = fileStartTime + time_interval*scan;

    // read latitude

    band_id = epr_get_band_id(fileID, "latitude");
    if (band_id == NULL) {
        printf("-E- %s line %d: Error finding band_id:latitude\n", __FILE__,
                __LINE__);
        return (1);
    }
    temp_raster = epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
    if (temp_raster == NULL) {
        printf("-E- %s line %d: Error allocating raster space.\n", __FILE__,
                __LINE__);
        return (1);
    }
    epr_read_band_raster(band_id, 0, scan, temp_raster);
    for (ip = 0; ip < npix; ip++)
        l1rec->lat[ip] = epr_get_pixel_as_float(temp_raster, spix + ip, 0);

    // read longitude

    band_id = epr_get_band_id(fileID, "longitude");
    if (band_id == NULL) {
        printf("-E- %s line %d: Error finding band_id:longitude\n", __FILE__,
                __LINE__);
        return (1);
    }
    epr_read_band_raster(band_id, 0, scan, temp_raster);
    for (ip = 0; ip < npix; ip++)
        l1rec->lon[ip] = epr_get_pixel_as_float(temp_raster, spix + ip, 0);
    epr_free_raster(temp_raster);

    // read sun zenith

    band_id = epr_get_band_id(fileID, "sun_zenith");
    if (band_id == NULL) {
        printf("-E- %s line %d: Error finding band_id:sun_zenith", __FILE__,
                __LINE__);
        return (1);
    }
    temp_raster = epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
    if (temp_raster == NULL) {
        printf("-E- %s line %d: Error allocating raster space.\n", __FILE__,
                __LINE__);
        return (1);
    }
    epr_read_band_raster(band_id, 0, scan, temp_raster);
    for (ip = 0; ip < npix; ip++)
        l1rec->solz[ip] = epr_get_pixel_as_float(temp_raster, spix + ip, 0);
    epr_free_raster(temp_raster);

    return (LIFE_IS_GOOD);
}

/**
 * @brief closes MERIS file, loads l1rec
 * @param[in]   file   file handle to MERIS file
 *
 * @author Paul E. Lyon NRL, Oct. 2006.
 */

int
closel1_meris_N1(filehandle *file) {
    if (epr_close_product(fileID)) {
        fprintf(stderr,
                "-E- %s line %d: epr_close_product failed for file, %s.\n",
                __FILE__, __LINE__, file->name);
        return (1);
    }

    if (detectorWL)
        free(detectorWL);
    if (detectorE0)
        free(detectorE0);

    return (LIFE_IS_GOOD);
}

