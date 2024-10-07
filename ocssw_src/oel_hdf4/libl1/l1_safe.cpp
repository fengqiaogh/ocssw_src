/* =============================================================================================== */
/* module l1_safe.cpp - functions to read MERIS and OLCI SAFE formatted (coastal color) for MSL12  */
/* Written By:  Richard Healy (SAIC) July 29, 2015.                                                */
/* Modified By:  Don Shea (SAIC) April 2, 2022.                                                    */
/*                                                                                                 */
/* =============================================================================================== */

#include <netcdf.h>

#include "l1_safe.h"
#include "smile.h"
#include "math.h"
#include <string>
#include <vector>

#include <gsl/gsl_spline2d.h>

using namespace std;

#define NEW_CACHE_SIZE 10000000
#define NEW_CACHE_NELEMS 23
#define NEW_CACHE_PREEMPTION 1.0

#define BANDINFO_FILENAME_MERIS "band_info_meris.txt"
#define BANDINFO_FILENAME_OLCI "band_info_olci.txt"

#define ANGLE_FILL_VALUE (-999)

static int16_t npix;
static int64_t *scan_start_tai;
static double lastvalidtime;
static int lastvalidscan = 0;
static vector<string> radFilename;
static vector<int32_t> radFileID;
static vector<string> radVarname;
static string geoCoordinatesFilename, tieGeometriesFilename, instrumentFilename, timeCoordinatesFilename, qualityFlagsFilename;
static int32_t geoCoordinatesFileID, tieGeometriesFileID, instrumentFileID, timeCoordinatesFileID, qualityFileID;


int openl1_safe(filehandle * l1file) {
    size_t source_w, source_h, source_b;
    int32_t nscan;
    int xid, yid, bid, retval, sds_id;
    int i;
    size_t start[3], count[3];
    unsigned short orbit_number;
    size_t attrLength;
    
    string indirStr;
    string tmpStr = l1file->name;
    size_t pos = tmpStr.rfind('/');
    if(pos != string::npos) {
        indirStr = tmpStr.substr(0,pos+1);
    }

    instrumentFilename = indirStr + "instrument_data.nc";
    timeCoordinatesFilename = indirStr + "time_coordinates.nc";
    geoCoordinatesFilename = indirStr + "geo_coordinates.nc";
    tieGeometriesFilename = indirStr + "tie_geometries.nc";
    qualityFlagsFilename = indirStr + "qualityFlags.nc";

    // set vector sizes
    radFilename.resize(l1file->nbands);
    radFileID.resize(l1file->nbands);
    radVarname.resize(l1file->nbands);

    for (i = 0; i < l1file->nbands; i++) {
        string numStr = to_string(i+1);
        if(numStr.size() == 1)
            numStr = (string) "0" + numStr;
        if(l1file->sensorID == MERIS)
            radFilename[i] = indirStr + "M" + numStr + "_radiance.nc";
        else
            radFilename[i] = indirStr + "Oa" + numStr + "_radiance.nc";

        if(want_verbose)
            printf("SAFE rad=%s\n", radFilename[i].c_str());
    }

    // Open the netcdf4 input file
    if (nc_set_chunk_cache(NEW_CACHE_SIZE, NEW_CACHE_NELEMS,
            NEW_CACHE_PREEMPTION) != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_set_chunk_cache (%s) failed.\n",
                __FILE__, __LINE__, l1file->name);
        exit(EXIT_FAILURE);
    }
    retval = nc_open(instrumentFilename.c_str(), NC_NOWRITE, &instrumentFileID);
    if (retval != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
                __FILE__, __LINE__, instrumentFilename.c_str());
        exit(EXIT_FAILURE);
    }

    // Get pixel and scan dimensions
    DPTB(nc_inq_dimid(instrumentFileID, "rows", &yid));
    DPTB(nc_inq_dimid(instrumentFileID, "columns", &xid));
    DPTB(nc_inq_dimid(instrumentFileID, "bands", &bid));
    DPTB(nc_inq_dimlen(instrumentFileID, xid, &source_w));
    DPTB(nc_inq_dimlen(instrumentFileID, yid, &source_h));
    DPTB(nc_inq_dimlen(instrumentFileID, bid, &source_b));
    if (source_b != (size_t)l1file->nbands) {
        fprintf(stderr, "-E- %s line %d: num bands (%d) does not match expected SAFE num bands (%d).\n",
                __FILE__, __LINE__, (int) source_b, l1file->nbands);
        exit(EXIT_FAILURE);
    }

    DPTB(nc_get_att_ushort(instrumentFileID, NC_GLOBAL, "absolute_orbit_number", &orbit_number));
    DPTB(nc_inq_attlen(instrumentFileID, NC_GLOBAL, "product_name", &attrLength));
    char* productName = (char*) allocateMemory(attrLength + 1, "product_name");
    DPTB(nc_get_att_text(instrumentFileID, NC_GLOBAL, "product_name", productName));
    if (strstr(productName, "ME_1_FRG")) {
        strcpy(l1file->spatialResolution, "300 m");
    } else if (strstr(productName, "ME_1_RRG")) {
        strcpy(l1file->spatialResolution, "1.2 km");
    } else if (strstr(productName, "OL_1_EFR")) {
        strcpy(l1file->spatialResolution, "300 m");
    } else if (strstr(productName, "OL_1_ERR")) {
        strcpy(l1file->spatialResolution, "1.2 km");
    } else {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) illegal product = \"%s\".\n",
                __FILE__, __LINE__, l1file->name, productName);
        exit(EXIT_FAILURE);
    }
    free(productName);

    if (want_verbose) {
        printf("%s L1B SAFE Npix  :%d Nscans:%d\n", sensorId2SensorName(l1file->sensorID), (int) source_w,
                (int) source_h);
    } // want_verbose

    npix = (int32_t) source_w;
    nscan = (int32_t) source_h;

    l1file->orbit_number = orbit_number;
    l1file->nbands = source_b;
    l1file->npix = npix;
    l1file->nscan = nscan;

    scan_start_tai = (int64_t *) calloc(nscan, sizeof (int64_t));

    start[0] = 0;
    count[0] = nscan;

    retval = nc_open(timeCoordinatesFilename.c_str(), NC_NOWRITE, &timeCoordinatesFileID);
    if (retval != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
                __FILE__, __LINE__, timeCoordinatesFilename.c_str());
        exit(EXIT_FAILURE);
    }

    DPTB(nc_inq_varid(timeCoordinatesFileID, "time_stamp", &sds_id)); //microseconds since 1/1/2000
    DPTB(nc_get_vara_long(timeCoordinatesFileID, sds_id, start, count, (long*) scan_start_tai));

    retval = nc_open(geoCoordinatesFilename.c_str(), NC_NOWRITE, &geoCoordinatesFileID);
    if (retval != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
                __FILE__, __LINE__, geoCoordinatesFilename.c_str());
        exit(EXIT_FAILURE);
    }

    retval = nc_open(tieGeometriesFilename.c_str(), NC_NOWRITE, &tieGeometriesFileID);
    if (retval != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
                __FILE__, __LINE__, tieGeometriesFilename.c_str());
        exit(EXIT_FAILURE);
    }

    retval = nc_open(qualityFlagsFilename.c_str(), NC_NOWRITE, &qualityFileID);
    if (retval != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
                __FILE__, __LINE__, qualityFlagsFilename.c_str());
        return (1);
    }

    for (i = 0; i < l1file->nbands; i++) {
        // Open each of the netcdf4 Lt files for each band
        retval = nc_open(radFilename[i].c_str(), NC_NOWRITE, &radFileID[i]);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_open failed for file, %s.\n",
                    __FILE__, __LINE__, radFilename[i].c_str());
            exit(EXIT_FAILURE);
        }
        string numStr = to_string(i+1);
        if(numStr.size() == 1)
            numStr = (string) "0" + numStr;
        if(l1file->sensorID == MERIS)
            radVarname[i] = (string) "M" + numStr + "_radiance";
        else
            radVarname[i] = (string) "Oa" + numStr + "_radiance";
    }

    return (LIFE_IS_GOOD);
}

int readl1_safe(filehandle *file, int32_t scan, l1str *l1rec) {
    static int firstCall = 1;
    static double time_interval;
    static double tai93at2000;
    static int32_t *lon, *lat;
    static int32_t lonFillValue, latFillValue;
    static int longitudeVarID, latitudeVarID;
    static int solaVarID, solzVarID, senaVarID, senzVarID;
    static double scale_lon, scale_lat;
    static double scale_solz, scale_sola, scale_senz, scale_sena;
    static int detectorIndexVarID;
    static int16_t *detector_index;
    static uint32_t *qualityFlags;
    static vector<int> radVarID;
    static vector<float> radScale;
    static vector<float> radOffset;
    static vector<uint16_t> radFillValue;

    static float *lambda0, *solar_flux;
    static size_t detectors;
    static uint16_t *rad_data;
    static size_t num_tie_cols, num_tie_rows;
    static uint16_t tie_row_subsample, tie_col_subsample;

    static double *solz_xa, *solz_ya, *solz_za;
    static double *sola_x_xa, *sola_x_ya, *sola_x_za;
    static double *sola_y_xa, *sola_y_ya, *sola_y_za;
    static double *senz_xa, *senz_ya, *senz_za;
    static double *sena_x_xa, *sena_x_ya, *sena_x_za;
    static double *sena_y_xa, *sena_y_ya, *sena_y_za;
    static gsl_spline2d *solz_spline;
    static gsl_spline2d *sola_x_spline;
    static gsl_spline2d *sola_y_spline;
    static gsl_spline2d *senz_spline;
    static gsl_spline2d *sena_x_spline;
    static gsl_spline2d *sena_y_spline;
    static gsl_interp_accel *solz_xacc, *solz_yacc;
    static gsl_interp_accel *sola_x_xacc, *sola_x_yacc;
    static gsl_interp_accel *sola_y_xacc, *sola_y_yacc;
    static gsl_interp_accel *senz_xacc, *senz_yacc;
    static gsl_interp_accel *sena_x_xacc, *sena_x_yacc;
    static gsl_interp_accel *sena_y_xacc, *sena_y_yacc;

    float *detectorTmpFloat;

    int retval;

    int32_t ip, ib, ipb;
    int32_t nbands = l1rec->l1file->nbands;
    size_t start[3], count[3];

    int xid, yid;

    if (firstCall) {
        firstCall = 0;
        if (want_verbose)
            printf("file->nbands = %d, l1rec->nbands = %d\n",
                (int) file->nbands, (int) l1rec->l1file->nbands);

        // set the size of the vectors
        radVarID.resize(nbands);
        radScale.resize(nbands);
        radOffset.resize(nbands);
        radFillValue.resize(nbands);

        // This is the time offset between what the time utilities need (1/1/1970) and what
        // the SAFE netcdf provides (1/1/2000)
        tai93at2000 = yds2tai93(2000, 1, 0.0);

        DPTB(nc_inq_dimid(tieGeometriesFileID, "tie_columns", &xid));
        DPTB(nc_inq_dimlen(tieGeometriesFileID, xid, &num_tie_cols));
        DPTB(nc_inq_dimid(tieGeometriesFileID, "tie_rows", &yid));
        DPTB(nc_inq_dimlen(tieGeometriesFileID, yid, &num_tie_rows));
        DPTB(nc_get_att_ushort(tieGeometriesFileID, NC_GLOBAL, "ac_subsampling_factor", &tie_col_subsample));
        DPTB(nc_get_att_ushort(tieGeometriesFileID, NC_GLOBAL, "al_subsampling_factor", &tie_row_subsample));

        if (((num_tie_rows-1) * tie_row_subsample + 1) != (size_t)file->nscan) {
            printf("-E- %s line %d: Sanity check failed - tie_rows (%d) x tie_row_pts (%d) < nscan (%d) in file %s\n",
                    __FILE__, __LINE__, (int) num_tie_rows, tie_row_subsample, file->nscan, tieGeometriesFilename.c_str());
            exit(EXIT_FAILURE);
        }
        if (((num_tie_cols - 1) * tie_col_subsample + 1) != (size_t)file->npix) {
            printf("-E- %s line %d: Sanity check failed - tie_cols (%d) x tie_col_pts (%d) != npix (%d) in file %s\n",
                    __FILE__, __LINE__, (int) num_tie_cols, tie_col_subsample, file->nscan, tieGeometriesFilename.c_str());
            exit(EXIT_FAILURE);
        }

        lon = (int32_t *) calloc(npix, sizeof (int32_t));
        lat = (int32_t *) calloc(npix, sizeof (int32_t));

        detector_index = (int16_t*) calloc(npix, sizeof (int16_t));
        qualityFlags = (uint32_t*) calloc(npix, sizeof (uint32_t));

        DPTB(nc_inq_dimid(instrumentFileID, "detectors", &xid));
        DPTB(nc_inq_dimlen(instrumentFileID, xid, &detectors));

        detectorTmpFloat = (float*) allocateMemory(nbands * detectors * sizeof (float), "detectorTmpFloat");
        lambda0 = (float*) allocateMemory(nbands * detectors * sizeof (float), "lambda0");
        solar_flux = (float*) allocateMemory(nbands * detectors * sizeof (float), "solar_flux");

        rad_data = (unsigned short *) calloc(npix, sizeof (unsigned short)); //BYSCAN

        start[0] = 0;
        start[1] = 0;
        start[2] = 0;
        count[0] = nbands;
        count[1] = detectors;
        count[2] = 0;

        if (l1_input->rad_opt != 0) {
            retval = nc_inq_varid(instrumentFileID, "lambda0", &xid);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: nc_get_varid failed for file, %s  field, %s.\n",
                        __FILE__, __LINE__, instrumentFilename.c_str(), "lambda0");
                exit(EXIT_FAILURE);
            }
            retval = nc_get_vara_float(instrumentFileID, xid, start, count, detectorTmpFloat);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                        __FILE__, __LINE__, instrumentFilename.c_str(), "lambda0");
                exit(EXIT_FAILURE);
            }
            // swap axis of the array
            for (int32_t i = 0; i < nbands; i++)
                for (size_t j = 0; j < detectors; j++)
                    lambda0[j * nbands + i] = detectorTmpFloat[i * detectors + j];

            retval = nc_inq_varid(instrumentFileID, "solar_flux", &xid);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                        __FILE__, __LINE__, instrumentFilename.c_str(), "solar_flux");
                exit(EXIT_FAILURE);
            }
            retval = nc_get_vara_float(instrumentFileID, xid, start, count, detectorTmpFloat);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                        __FILE__, __LINE__, instrumentFilename.c_str(), "solar_flux");
                exit(EXIT_FAILURE);
            }
            // swap axis of the array
            for (int32_t i = 0; i < nbands; i++)
                for (size_t j = 0; j < detectors; j++)
                    solar_flux[j * nbands + i] = detectorTmpFloat[i * detectors + j];
            free(detectorTmpFloat);

            char *tmp_str;
            char filename[FILENAME_MAX];
            if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
                printf("OCDATAROOT environment variable is not defined.\n");
                exit(EXIT_FAILURE);
            }
            strcpy(filename, tmp_str);
            strcat(filename, "/");
            strcat(filename, sensorId2SensorDir(l1rec->l1file->sensorID));
            if(l1rec->l1file->sensorID == MERIS) {
                strcat(filename, "/cal/");
                strcat(filename, BANDINFO_FILENAME_MERIS);
            } else {
                strcat(filename, "/");
                strcat(filename, subsensorId2SubsensorDir(l1rec->l1file->subsensorID));
                strcat(filename, "/cal/");
                strcat(filename, BANDINFO_FILENAME_OLCI);
            }

            smile_init(nbands, detectors, filename, lambda0, solar_flux);
        } // if rad_opt

        DPTB(nc_inq_varid(geoCoordinatesFileID, "longitude", &longitudeVarID));
        DPTB(nc_get_att_double(geoCoordinatesFileID, longitudeVarID, "scale_factor", &scale_lon));
        lonFillValue = -2147483648;
        retval = nc_get_att_int(geoCoordinatesFileID, longitudeVarID, "_FillValue", &lonFillValue);
        DPTB(nc_inq_varid(geoCoordinatesFileID, "latitude", &latitudeVarID));
        DPTB(nc_get_att_double(geoCoordinatesFileID, latitudeVarID, "scale_factor", &scale_lat));
        latFillValue = -2147483648;
        retval = nc_get_att_int(geoCoordinatesFileID, latitudeVarID, "_FillValue", &latFillValue);

        DPTB(nc_inq_varid(tieGeometriesFileID, "SAA", &solaVarID));
        DPTB(nc_get_att_double(tieGeometriesFileID, solaVarID, "scale_factor", &scale_sola));
        DPTB(nc_inq_varid(tieGeometriesFileID, "SZA", &solzVarID));
        DPTB(nc_get_att_double(tieGeometriesFileID, solzVarID, "scale_factor", &scale_solz));

        DPTB(nc_inq_varid(tieGeometriesFileID, "OAA", &senaVarID));
        DPTB(nc_get_att_double(tieGeometriesFileID, senaVarID, "scale_factor", &scale_sena));
        DPTB(nc_inq_varid(tieGeometriesFileID, "OZA", &senzVarID));
        DPTB(nc_get_att_double(tieGeometriesFileID, senzVarID, "scale_factor", &scale_senz));

        DPTB(nc_inq_varid(instrumentFileID, "detector_index", &detectorIndexVarID));

        for (ib = 0; ib < nbands; ib++) {
            retval = nc_inq_varid(radFileID[ib], radVarname[ib].c_str(), &radVarID[ib]);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: nc_inq_varid failed for file=%s, field=%s.\n",
                        __FILE__, __LINE__, radFilename[ib].c_str(), radVarname[ib].c_str());
                exit(EXIT_FAILURE);
            }
            retval = nc_get_att_float(radFileID[ib], radVarID[ib], "scale_factor", &radScale[ib]);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: nc_get_att_float failed for file=%s, field=%s, attribute=scale_factor.\n",
                        __FILE__, __LINE__, radFilename[ib].c_str(), radVarname[ib].c_str());
                exit(EXIT_FAILURE);
            }
            retval = nc_get_att_float(radFileID[ib], radVarID[ib], "add_offset", &radOffset[ib]);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: nc_get_att_float failed for file=%s, field=%s, attribute=add_offset.\n",
                        __FILE__, __LINE__, radFilename[ib].c_str(), radVarname[ib].c_str());
                exit(EXIT_FAILURE);
            }
            retval = nc_get_att_ushort(radFileID[ib], radVarID[ib], "_FillValue", &radFillValue[ib]);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: nc_get_att_float failed for file=%s, field=%s, attribute=_FillValue.\n",
                        __FILE__, __LINE__, radFilename[ib].c_str(), radVarname[ib].c_str());
                exit(EXIT_FAILURE);
            }

        }

        int num_tie_points = num_tie_cols * num_tie_rows;
        int32_t *tmp_int = (int32_t*) allocateMemory(num_tie_points * sizeof(int32_t), "tmp_int array for sol and sen");
        uint32_t *tmp_uint = (uint32_t*) tmp_int;

        solz_xa = (double*) allocateMemory(num_tie_cols * sizeof(double), "solz_xa");
        solz_ya = (double*) allocateMemory(num_tie_rows * sizeof(double), "solz_ya");
        solz_za = (double*) allocateMemory(num_tie_points * sizeof(double), "solz_za");
        sola_x_xa = (double*) allocateMemory(num_tie_cols * sizeof(double), "sola_pos_xa");
        sola_x_ya = (double*) allocateMemory(num_tie_rows * sizeof(double), "sola_pos_ya");
        sola_x_za = (double*) allocateMemory(num_tie_points * sizeof(double), "sola_pos_za");
        sola_y_xa = (double*) allocateMemory(num_tie_cols * sizeof(double), "sola_neg_xa");
        sola_y_ya = (double*) allocateMemory(num_tie_rows * sizeof(double), "sola_neg_ya");
        sola_y_za = (double*) allocateMemory(num_tie_points * sizeof(double), "sola_neg_za");
        senz_xa = (double*) allocateMemory(num_tie_cols * sizeof(double), "senz_xa");
        senz_ya = (double*) allocateMemory(num_tie_rows * sizeof(double), "senz_ya");
        senz_za = (double*) allocateMemory(num_tie_points * sizeof(double), "senz_za");
        sena_x_xa = (double*) allocateMemory(num_tie_cols * sizeof(double), "sena_pos_xa");
        sena_x_ya = (double*) allocateMemory(num_tie_rows * sizeof(double), "sena_pos_ya");
        sena_x_za = (double*) allocateMemory(num_tie_points * sizeof(double), "sena_pos_za");
        sena_y_xa = (double*) allocateMemory(num_tie_cols * sizeof(double), "sena_neg_xa");
        sena_y_ya = (double*) allocateMemory(num_tie_rows * sizeof(double), "sena_neg_ya");
        sena_y_za = (double*) allocateMemory(num_tie_points * sizeof(double), "sena_neg_za");
        int index = 0;
        for(size_t i=0; i<num_tie_cols; i++) {
            solz_xa[i] = index;
            sola_x_xa[i] = index;
            sola_y_xa[i] = index;
            senz_xa[i] = index;
            sena_x_xa[i] = index;
            sena_y_xa[i] = index;
            index += tie_col_subsample;
        }
        index = 0;
        for(size_t i=0; i<num_tie_rows; i++) {
            solz_ya[i] = index;
            sola_x_ya[i] = index;
            sola_y_ya[i] = index;
            senz_ya[i] = index;
            sena_x_ya[i] = index;
            sena_y_ya[i] = index;
            index += tie_row_subsample;
        }
        DPTB(nc_get_var_uint(tieGeometriesFileID, solzVarID, tmp_uint));
        for(int i=0; i<num_tie_points; i++) {
            solz_za[i] = tmp_uint[i] * scale_solz;
        }
        DPTB(nc_get_var_int(tieGeometriesFileID, solaVarID, tmp_int));
        for(int i=0; i<num_tie_points; i++) {
            double angle = tmp_int[i] * scale_sola / RADEG;
            sola_x_za[i] = cos(angle);
            sola_y_za[i] = sin(angle);
        }
        DPTB(nc_get_var_uint(tieGeometriesFileID, senzVarID, tmp_uint));
        for(int i=0; i<num_tie_points; i++) {
            senz_za[i] = tmp_uint[i] * scale_senz;
        }
        DPTB(nc_get_var_int(tieGeometriesFileID, senaVarID, tmp_int));
        for(int i=0; i<num_tie_points; i++) {
            double angle = tmp_int[i] * scale_sena / RADEG;
            sena_x_za[i] = cos(angle);
            sena_y_za[i] = sin(angle);
        }
        const gsl_interp2d_type *splineType = gsl_interp2d_bilinear;
        //const gsl_interp2d_type *splineType = gsl_interp2d_bicubic;
        solz_spline = gsl_spline2d_alloc(splineType, num_tie_cols, num_tie_rows);
        sola_x_spline = gsl_spline2d_alloc(splineType, num_tie_cols, num_tie_rows);
        sola_y_spline = gsl_spline2d_alloc(splineType, num_tie_cols, num_tie_rows);
        senz_spline = gsl_spline2d_alloc(splineType, num_tie_cols, num_tie_rows);
        sena_x_spline = gsl_spline2d_alloc(splineType, num_tie_cols, num_tie_rows);
        sena_y_spline = gsl_spline2d_alloc(splineType, num_tie_cols, num_tie_rows);
        solz_xacc = gsl_interp_accel_alloc();
        solz_yacc = gsl_interp_accel_alloc();
        sola_x_xacc = gsl_interp_accel_alloc();
        sola_x_yacc = gsl_interp_accel_alloc();
        sola_y_xacc = gsl_interp_accel_alloc();
        sola_y_yacc = gsl_interp_accel_alloc();
        senz_xacc = gsl_interp_accel_alloc();
        senz_yacc = gsl_interp_accel_alloc();
        sena_x_xacc = gsl_interp_accel_alloc();
        sena_x_yacc = gsl_interp_accel_alloc();
        sena_y_xacc = gsl_interp_accel_alloc();
        sena_y_yacc = gsl_interp_accel_alloc();
        gsl_spline2d_init(solz_spline, solz_xa, solz_ya, solz_za, num_tie_cols, num_tie_rows);
        gsl_spline2d_init(sola_x_spline, sola_x_xa, sola_x_ya, sola_x_za, num_tie_cols, num_tie_rows);
        gsl_spline2d_init(sola_y_spline, sola_y_xa, sola_y_ya, sola_y_za, num_tie_cols, num_tie_rows);
        gsl_spline2d_init(senz_spline, senz_xa, senz_ya, senz_za, num_tie_cols, num_tie_rows);
        gsl_spline2d_init(sena_x_spline, sena_x_xa, sena_x_ya, sena_x_za, num_tie_cols, num_tie_rows);
        gsl_spline2d_init(sena_y_spline, sena_y_xa, sena_y_ya, sena_y_za, num_tie_cols, num_tie_rows);

        free(tmp_int);
    } // first call


    //      set time for this scan - if scan_start_time value not properly set, estimate from scene start time.
    if (scan_start_tai[scan] > 0) {
        l1rec->scantime = tai93_to_unix(scan_start_tai[scan] / 1000000.0 + tai93at2000);
        lastvalidtime = l1rec->scantime;
        if (scan > 0)
            time_interval = (scan_start_tai[scan] - lastvalidtime) / (scan - lastvalidscan);
        lastvalidscan = scan;
    } else {
        l1rec->scantime = lastvalidtime + (time_interval * (scan - lastvalidscan));
        lastvalidtime = l1rec->scantime;
    }

    // calculate fsol which is needed by radcor
    int16_t syear, sday;
    double secs;
    unix2yds(l1rec->scantime, &syear, &sday, &secs);
    int32_t yr = syear;
    int32_t dy = sday;
    int32_t msec = (int32_t) (secs * 1000.0);
    double esdist = esdist_(&yr, &dy, &msec);
    l1rec->fsol = pow(1.0 / esdist, 2);

    start[0] = scan;
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = npix;
    count[2] = 0;

    //until geolocation is read, set fill values -
    for (ip = 0; ip < npix; ip++) {
        l1rec->navfail[ip] = 0;
        l1rec->lon[ip] = ANGLE_FILL_VALUE;
        l1rec->lat[ip] = ANGLE_FILL_VALUE;

        l1rec->solz[ip] = ANGLE_FILL_VALUE;
        l1rec->sola[ip] = ANGLE_FILL_VALUE;
        l1rec->senz[ip] = ANGLE_FILL_VALUE;
        l1rec->sena[ip] = ANGLE_FILL_VALUE;
    }

    retval = nc_get_vara_int(geoCoordinatesFileID, longitudeVarID, start, count, lon);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, geoCoordinatesFilename.c_str(), "longitude");
        exit(EXIT_FAILURE);
    }
    retval = nc_get_vara_int(geoCoordinatesFileID, latitudeVarID, start, count, lat);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, geoCoordinatesFilename.c_str(), "lat");
        exit(EXIT_FAILURE);
    }
    retval = nc_get_vara_short(instrumentFileID, detectorIndexVarID, start, count, detector_index);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_short failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, instrumentFilename.c_str(), "detector_index");
        exit(EXIT_FAILURE);
    }
    for (ip = 0; ip < npix; ip++) {
        l1rec->pixnum[ip] = ip;
        l1rec->flags[ip] = 0;
        l1rec->pixdet[ip] = detector_index[ip];
        if(lon[ip] == lonFillValue || lat[ip] == latFillValue) {
            l1rec->lon[ip] = l1rec->lat[ip] = ANGLE_FILL_VALUE;
            l1rec->navfail[ip] = 1;
        } else {
            l1rec->lon[ip] = lon[ip] * scale_lon;
            l1rec->lat[ip] = lat[ip] * scale_lat;
        }
    }

    retval = nc_inq_varid(qualityFileID, "quality_flags", &xid);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_varid failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, qualityFlagsFilename.c_str(), "quality_flags");
        exit(FATAL_ERROR);
    }
    retval = nc_get_vara_uint(qualityFileID, xid, start, count, qualityFlags);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_uint failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, qualityFlagsFilename.c_str(), "quality_flags");
        exit(FATAL_ERROR);
    }

    // read in radiance data
    for (ib = 0; ib < nbands; ib++) {
        retval = nc_get_vara_ushort(radFileID[ib], radVarID[ib], start, count, rad_data); //BYSCAN
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, radFilename[ib].c_str(), radVarname[ib].c_str());
            exit(FATAL_ERROR);
        }
        // copy to Lt record.
        for (ip = 0; ip < npix; ip++) {
            ipb = ip * nbands + ib;
            if(rad_data[ip] == radFillValue[ib]) {
                l1rec->Lt[ipb] = BAD_FLT;
                l1rec->navfail[ip] = 1;
            } else {
                l1rec->Lt[ipb] = (rad_data[ip] * radScale[ib] + radOffset[ib]) / 10.; //BYSCAN

                // mark negative input data as HILT
                if (l1rec->Lt[ipb] < 0.0)
                    l1rec->hilt[ip] = 1;
            }
        }
    } // for ib

    for (ip = 0; ip < npix; ip++) {
        double x, y;
        l1rec->solz[ip] = gsl_spline2d_eval(solz_spline, ip, scan, solz_xacc, solz_yacc);
        x = gsl_spline2d_eval(sola_x_spline, ip, scan, sola_x_xacc, sola_x_yacc);
        y = gsl_spline2d_eval(sola_y_spline, ip, scan, sola_y_xacc, sola_y_yacc);
        l1rec->sola[ip] = atan2(y, x) * RADEG;
        l1rec->senz[ip] = gsl_spline2d_eval(senz_spline, ip, scan, senz_xacc, senz_yacc);
        x = gsl_spline2d_eval(sena_x_spline, ip, scan, sena_x_xacc, sena_x_yacc);
        y = gsl_spline2d_eval(sena_y_spline, ip, scan, sena_y_xacc, sena_y_yacc);
        l1rec->sena[ip] = atan2(y, x) * RADEG;

        // radcor needs all Lts populated before running so...one more time around the block...
        if (l1_input->rad_opt != 0) {
            uint32_t isLand = qualityFlags[ip] & 2147483648;
            // es corrected f0, so setting 1 as 4 element - seemed silly to make a variable for it
            radcor(l1rec, ip, isLand, 1);
            for (ib = 0; ib < nbands; ib++) {
                if(l1rec->navfail[ip] != 1) {
                    ipb = ip * nbands + ib;
                    l1rec->Lt[ipb] += l1rec->radcor[ipb];
                }
            }
        }
    }

    l1rec->npix = file->npix;
    l1rec->l1file->terrain_corrected = 1;

    return (LIFE_IS_GOOD);
}

int closel1_safe(filehandle *file) {
    int i;

    free(scan_start_tai);

    DPTB(nc_close(instrumentFileID));
    DPTB(nc_close(timeCoordinatesFileID));
    DPTB(nc_close(geoCoordinatesFileID));
    DPTB(nc_close(tieGeometriesFileID));
    DPTB(nc_close(qualityFileID));
    for (i = 0; i < file->nbands; i++) {
        DPTB(nc_close(radFileID[i]));
    }

    return (LIFE_IS_GOOD);
}

