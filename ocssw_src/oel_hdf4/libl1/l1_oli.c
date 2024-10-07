#include <stdlib.h>
#include <stdio.h>

#include "l1.h"
#include "l1_oli.h"

#include <tiffio.h>
#include <geotiff.h>
#include <xtiffio.h>
#include <geo_normalize.h>
#include <filetype.h>
#include <libnav.h>
#include <libgen.h>
#include "xcal.h"

static const int itemSize = 500;
static const int maxReflBands = 8;
static float *Fobar;
// following added by Sudipta for FPM based correction
static double      *xcal_factor = NULL;               /* xcal factors for each FPM and band */
static int tile_exist;
static uint32_t tileLength = 0;
static uint32_t tileWidth = 0;
static uint32_t imageWidth = 0;

typedef struct oli_struct {
    int32_t year, doy, msec;
    double sunAzimuth, sunElevation;
    double *lat, *lon;
    double *scale, *offset;
    double *refl_scale, *refl_offset;
    int line_num_cached;
    TIFF** tif; // file handle for each band
    GTIF* gtif; // geotiff handle for the first file
    GTIFDefn* defn; // geotiff definition structure for first file
    uint16_t* buf_tile_row; // buffer used to read a row of tile * 9 bands from TIFF file
    uint16_t* buf; // buffer used to read one scan line from TIFF file
} oli_t;

oli_t* createPrivateData(int numBands) {
    int i;

    oli_t* data = (oli_t*) calloc(1, sizeof (oli_t));
    if (data == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to allocate private data for OLI\n",
                __FILE__, __LINE__);
        exit(1);
    }

    data->scale = (double *) malloc(numBands * sizeof (double));
    data->offset = (double *) malloc(numBands * sizeof (double));
    if (data->scale == NULL || data->offset == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to allocate scale/offset data for OLI\n",
                __FILE__, __LINE__);
        exit(1);
    }

    data->refl_scale = (double *) malloc(numBands * sizeof (double));
    data->refl_offset = (double *) malloc(numBands * sizeof (double));
    if (data->refl_scale == NULL || data->refl_offset == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to allocate reflectance scale/offset data for OLI\n",
                __FILE__, __LINE__);
        exit(1);
    }

    for (i = 0; i < numBands; i++) {
        data->scale[i] = BAD_FLT;
        data->offset[i] = BAD_FLT;
        data->refl_scale[i] = BAD_FLT;
        data->refl_offset[i] = BAD_FLT;
    }

    data->tif = (TIFF**) calloc(numBands, sizeof (TIFF*));
    if (data->tif == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to allocate TIFF pointers for OLI\n",
                __FILE__, __LINE__);
        exit(1);
    }

    data->defn = (GTIFDefn*) malloc(sizeof (GTIFDefn));
    if (data->defn == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to allocate GEOTIFF definition structure for OLI\n",
                __FILE__, __LINE__);
        exit(1);
    }
    
    data->line_num_cached = -1; //

    return data;
}

void freePrivateData(oli_t* data) {
    free(data->scale);
    free(data->offset);
    free(data->refl_scale);
    free(data->refl_offset);
    free(data->tif);
    free(data->defn);
    free(data);
}

int readNextLine(FILE* fp, char* tag, int* i, char* val) {
    char* result;
    char line[itemSize];
    int count;

    result = fgets(line, itemSize, fp);
    if (result == NULL) {
        return 0;
    }
    trimBlanks(line);

    count = sscanf(line, "%s = %s", tag, val);
    if (count != 2) {
        // not found so return blank line
        tag[0] = 0;
        *i = 0;
        val[0] = 0;
        return 1;
    }

    // grab index if it exists
    result = strrchr(tag, '_');
    if (result) {
        result++;
        *i = atoi(result);
    } else {
        *i = 0;
    }

    // get rid of the quotes from val
    if (val[0] == '"')
        val[0] = ' ';
    count = strlen(val) - 1;
    if (val[count] == '"')
        val[count] = 0;
    trimBlanks(val);
    return 1;
}

int openl1_oli(filehandle *file) {
    int i;
    FILE *fp;
    char tag[itemSize];
    char val[itemSize];
    char dateStr[32];
    char timeStr[32];
    char fileName[FILENAME_MAX];

    oli_t* data = file->private_data = createPrivateData(maxReflBands);

    if (want_verbose)
        printf("OLI Level-1B %s\n", file->name);

    /* Open file */
    if ((fp = fopen(file->name, "r")) == NULL) {
        fprintf(stderr, "-E- %s line %d: unable open %s\n",
                __FILE__, __LINE__, file->name);
        exit(1);
    }

    int collectionNumber = 0;
    int dateFound = 0;
    int timeFound = 0;
    int numLinesFound = 0;
    int numSamplesFound = 0;
    int sunAzimuthFound = 0;
    int sunElevationFound = 0;
    int angleFileFound = 0;

    // loop metadata
    while (readNextLine(fp, tag, &i, val)) {

        // skip blank lines
        if (tag[0] == 0)
            continue;

        // get collection number
        if (!strcmp(tag, "COLLECTION_NUMBER")) {
            collectionNumber = atoi(val);
            if (want_verbose)
                printf("OLI Level-1B Collection %d\n", collectionNumber);

            // get date
        } else if (!strcmp(tag, "DATE_ACQUIRED")) {
            dateFound = 1;
            strcpy(dateStr, val);

            // get time
        } else if (!strcmp(tag, "SCENE_CENTER_TIME")) {
            timeFound = 1;
            strcpy(timeStr, val);

            // get num lines
        } else if (!strcmp(tag, "REFLECTIVE_LINES")) {
            numLinesFound = 1;
            file->nscan = atoi(val);

            // get num samples
        } else if (!strcmp(tag, "REFLECTIVE_SAMPLES")) {
            numSamplesFound = 1;
            file->npix = atoi(val);

            // get band file names
        } else if (!strncmp(tag, "FILE_NAME_BAND_", 15)) {
            i--;
            if ((i >= 0 && i < file->nbands) || i == 8) {
                // dirname might destroy input, so we pass it a copy
                char dir[FILENAME_MAX];
                strcpy(dir, file->name);
                strcpy(fileName, dirname(dir));
                strcat(fileName, "/");
                strcat(fileName, val);
                if (want_verbose)
                    printf("OLI Level-1B Band[%d]:%s\n", i, fileName);
                if (i != 8) {
                    data->tif[i] = XTIFFOpen(fileName, "r"); 
                    if (!data->tif[i]) {
                        fprintf(stderr, "-E- %s line %d: unable open TIFF file %s\n",
                                __FILE__, __LINE__, fileName);
                        exit(1);
                    }
                }
                else {
                    data->tif[i - 1] = XTIFFOpen(fileName, "r"); // open Band 9  
                    if (!data->tif[i -1 ]) {
                        fprintf(stderr, "-E- %s line %d: unable open TIFF file %s\n",
                                __FILE__, __LINE__, fileName);
                        exit(1);
                    }
                }
            }

            // get geofile name
        } else if (!strcmp(tag, "ANGLE_COEFFICIENT_FILE_NAME") || !strcmp(tag, "FILE_NAME_ANGLE_COEFFICIENT")) {
            angleFileFound = 1;
            if ((file->geofile == NULL) || (file->geofile[0] == 0)) {
                // dirname might destroy input, so we pass it a copy
                char dir[FILENAME_MAX];
                strcpy(dir, file->name);
                strcpy(fileName, dirname(dir));
                strcat(fileName, "/");
                strcat(fileName, val);
                file->geofile = strdup(fileName);
                if (want_verbose)
                    printf("OLI Level-1B Angle File:%s\n", fileName);
            }

            // get sun azimuth
        } else if (!strcmp(tag, "SUN_AZIMUTH")) {
            sunAzimuthFound = 1;
            data->sunAzimuth = atof(val);

            // get sun elevation
        } else if (!strcmp(tag, "SUN_ELEVATION")) {
            sunElevationFound = 1;
            data->sunElevation = atof(val);

            // get refl_scale
        } else if (!strncmp(tag, "REFLECTANCE_MULT_BAND_", 22)) {
            i--;
            if (i >= 0 && i < maxReflBands) {
                data->refl_scale[i] = atof(val);
            }

            // get refl_offset
        } else if (!strncmp(tag, "REFLECTANCE_ADD_BAND_", 21)) {
            i--;
            if (i >= 0 && i < maxReflBands) {
                data->refl_offset[i] = atof(val);
            }

        } else if (!strncmp(tag, "RADIANCE_MULT_BAND_", 19)) {
            i--;
            if (i >= 0 && i < maxReflBands) {
                data->scale[i] = atof(val);
            }

            // get offset
        } else if (!strncmp(tag, "RADIANCE_ADD_BAND_", 18)) {
            i--;
            if (i >= 0 && i < maxReflBands) {
                data->offset[i] = atof(val);
            }
        }

    } // while

    fclose(fp);

    if (!dateFound) {
        fprintf(stderr, "-E- %s line %d: Did not find DATE_ACQUIRED in %s\n",
                __FILE__, __LINE__, file->name);
        exit(EXIT_FAILURE);
    }
    if (!timeFound) {
        fprintf(stderr, "-E- %s line %d: Did not find SCENE_CENTER_TIME in %s\n",
                __FILE__, __LINE__, file->name);
        exit(EXIT_FAILURE);
    }
    if (!numLinesFound) {
        fprintf(stderr, "-E- %s line %d: Did not find REFLECTIVE_LINES in %s\n",
                __FILE__, __LINE__, file->name);
        exit(EXIT_FAILURE);
    }
    if (!numSamplesFound) {
        fprintf(stderr, "-E- %s line %d: Did not find REFLECTIVE_SAMPLES in %s\n",
                __FILE__, __LINE__, file->name);
        exit(EXIT_FAILURE);
    }
    if (!sunAzimuthFound) {
        fprintf(stderr, "-E- %s line %d: Did not find SUN_AZIMUTH in %s\n",
                __FILE__, __LINE__, file->name);
        exit(EXIT_FAILURE);
    }
    if (!sunElevationFound) {
        fprintf(stderr, "-E- %s line %d: Did not find SUN_ELEVATION in %s\n",
                __FILE__, __LINE__, file->name);
        exit(EXIT_FAILURE);
    }
    if (collectionNumber > 0) {
        if (!angleFileFound) {
            fprintf(stderr, "-E- %s line %d: Did not find ANGLE_COEFFICIENT_FILE_NAME in %s\n",
                    __FILE__, __LINE__, file->name);
            exit(EXIT_FAILURE);
        }
    }
    for (i = 0; i < maxReflBands; i++) {
        if (data->tif[i] == NULL) {
            fprintf(stderr, "-E- %s line %d: Did not find FILE_NAME_BAND_%d in %s\n",
                    __FILE__, __LINE__, i + 1, file->name);
            exit(EXIT_FAILURE);
        }
        if (data->scale[i] == BAD_FLT) {
            fprintf(stderr, "-E- %s line %d: Did not find RADIANCE_MULT_BAND_%d in %s\n",
                    __FILE__, __LINE__, i + 1, file->name);
            exit(EXIT_FAILURE);
        }
        if (data->offset[i] == BAD_FLT) {
            fprintf(stderr, "-E- %s line %d: Did not find RADIANCE_ADD_BAND_%d in %s\n",
                    __FILE__, __LINE__, i + 1, file->name);
            exit(EXIT_FAILURE);
        }
        if (data->refl_scale[i] == BAD_FLT) {
            fprintf(stderr, "-E- %s line %d: Did not find REFLECTANCE_MULT_BAND_%d in %s\n",
                    __FILE__, __LINE__, i + 1, file->name);
            exit(EXIT_FAILURE);
        }
        if (data->refl_offset[i] == BAD_FLT) {
            fprintf(stderr, "-E- %s line %d: Did not find REFLECTANCE_ADD_BAND_%d in %s\n",
                    __FILE__, __LINE__, i + 1, file->name);
            exit(EXIT_FAILURE);
        }
    }

    // allocate lat and lon storage
    data->lat = (double *) malloc(file->npix * sizeof (double));
    data->lon = (double *) malloc(file->npix * sizeof (double));
    if (data->lat == NULL || data->lon == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to allocate lat/lon data for OLI\n",
                __FILE__, __LINE__);
        exit(1);
    }

    // only need the GEO TIFF info from one file
    data->gtif = GTIFNew(data->tif[0]);
    if (!data->gtif) {
        fprintf(stderr, "-E- %s line %d: unable open GEOTIFF file %s\n",
                __FILE__, __LINE__, fileName);
        exit(1);
    }

    if (!GTIFGetDefn(data->gtif, data->defn)) {
        fprintf(stderr, "-E- %s line %d: unable populate GEOTIFF defn structure for %s\n",
                __FILE__, __LINE__, fileName);
        exit(1);
    }

    // allocate buffer to hold one scan line
    int size = TIFFScanlineSize(data->tif[0]);
    if (size != file->npix * 2) {
        fprintf(stderr, "-E- %s line %d: unexpected pixel data size in %s\n",
                __FILE__, __LINE__, fileName);
        exit(1);
    }
    data->buf = (uint16_t*) malloc(size);


    tile_exist = TIFFGetField(data->tif[0], TIFFTAG_TILELENGTH, &tileLength);
        
    if (tile_exist != 0) {    
        TIFFGetField(data->tif[0], TIFFTAG_TILEWIDTH, &tileWidth);
        TIFFGetField(data->tif[0], TIFFTAG_IMAGEWIDTH, &imageWidth);   
        // allocate buffer to hold a row of tile * 8 bands
        int size_tile_row =  maxReflBands * tileLength * size;
        data->buf_tile_row = (uint16_t*) malloc(size_tile_row);
    }

    // get date "2013-07-19"
    int year, month, day;
    sscanf(dateStr, "%d-%d-%d", &year, &month, &day);

    // get time "10:41:59.9930740Z"
    int hour, minute;
    double sec;
    sscanf(timeStr, "%d:%d:%lf", &hour, &minute, &sec);

    int isec = (int) sec;
    ymdhms2ydmsec(year, month, day, hour, minute, isec,
            &data->year, &data->doy, &data->msec);
    sec -= isec;
    data->msec += sec * 1000;

    if (want_verbose) {
        printf("OLI Start Time: %4d-%02d-%02d %03d %02d:%02d:%f\n",
                year, month, day, data->doy, hour, minute, sec + isec);

        printf("OLI file has %d bands, %d samples, %d lines\n",
                file->nbands, file->npix, file->nscan);
    }
    strcpy(file->spatialResolution, "30 m");
    file->terrain_corrected = 1;
    /*
     *  get the Fobar here to set up Fo
     */
    rdsensorinfo(file->sensorID, l1_input->evalmask, "Fobar", (void **) &Fobar);

    return 0;
}

int readl1_oli(filehandle *file, int recnum, l1str *l1rec, int lonlat) {
    int ip, ib, ipb;
    oli_t* data = (oli_t*) file->private_data;
    
    /* Sudipta new addition for OLI to extract SCA and Det numbers */
    static short  *sca_num;        /* SCA number                     */
    static short  *det_num;        /* Detector number                */
    int         isb, iw;
    static int firstCall = 1;
     /* End Sudipta new addition for OLI */

    //  set information about data
    l1rec->npix = file->npix;
    l1rec->scantime = yds2unix((int16_t) data->year, (int16_t) data->doy,
            (double) (data->msec / 1000.0));

    double esdist = esdist_(&data->year, &data->doy, &data->msec);
    double fsol = pow(1.0 / esdist, 2);

    //  get lat-lon
    for (ip = 0; ip < file->npix; ip++) {
        data->lon[ip] = ip;
        data->lat[ip] = recnum;
        GTIFImageToPCS(data->gtif, data->lon + ip, data->lat + ip);
    }

    if (!GTIFProj4ToLatLong(data->defn, file->npix, data->lon, data->lat)) {
        fprintf(stderr, "-E- %s line %d: unable reproject points for scan %d\n",
                __FILE__, __LINE__, recnum);
        exit(1);
    }

    for (ip = 0; ip < file->npix; ip++) {

        l1rec->pixnum[ip] = ip;

        if (isnan(data->lat[ip]))
            data->lat[ip] = -999.0;
        if (isnan(data->lon[ip]))
            data->lon[ip] = -999.0;
        l1rec->lat[ip] = data->lat[ip];
        l1rec->lon[ip] = data->lon[ip];

        if (l1rec->lon[ip] < -181.0 || l1rec->lon[ip] > 181.0 ||
                l1rec->lat[ip] < -91.0 || l1rec->lat[ip] > 91.0)
            l1rec->navfail[ip] = 1;

    }
    
    // read path angles if user supplied, else use best guess

    // new addition for OLI
    /*
     *  if required for that record, set up the geom_per_band storage
     */
    if ((l1_input->geom_per_band == 1) && (l1rec->geom_per_band == NULL)) {
        init_geom_per_band(l1rec);
        //      gm_p_b = l1rec->geom_per_band; // store this address so that it can be later destroyed in close()
    }

    if (firstCall) {
        firstCall = 0;
        // Sudipta - new addition to read the FPM xcal file and ingest the values
        for (iw = 0; iw < l1_input->xcal_nwave; iw++)
            if ((l1_input->xcal_opt[iw] & XCALOLI) != 0) {
                xcal_factor = get_fpm_xcal(l1_input->xcal_file);
                break;
            }
        /* printf("geom_per_band = %d\n",l1_input->geom_per_band);
        printf("xcal_opt = %d\n", l1_input->xcal_opt[0]);
        printf("xcal_factor = %lf\n", xcal_factor[0]);
        exit(0)*/;

        // allocate arrays for the sca and detector number arrays for all bands and pixels
        /* Sudipta new addition for OLI to extract SCA and Det numbers */
        if ((sca_num = (short *) calloc(file->npix * file->nbands, sizeof (short))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate SCA number array.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        if ((det_num = (short *) calloc(file->npix * file->nbands, sizeof (short))) == NULL) {
            free(sca_num);
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate detector number array.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        /* End Sudipta new addition for OLI */
                    
    }

    if (file->geofile == NULL || file->geofile[0] == 0) {
        if (get_oli_nom_angles(file->name, file->npix, file->nscan, recnum,
                l1rec->solz, l1rec->sola, l1rec->senz, l1rec->sena, sca_num,
                det_num, 0) == -1) {
            fprintf(stderr, "-E- %s line %d: Unable to compute geometries\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        // Sudipta added sca_num and det_num to get_oli_nom_angles
        // also added provision for band based geometry
        // band based geometry requested
        if (l1_input->geom_per_band) {
            if (get_oli_nom_angles(file->name, file->npix, file->nscan, recnum,
                    l1rec->geom_per_band->solz, l1rec->geom_per_band->sola,
                    l1rec->geom_per_band->senz, l1rec->geom_per_band->sena,
                    sca_num, det_num, l1_input->geom_per_band) == -1) {
                fprintf(stderr, "-E- %s line %d: Unable to compute geometries\n",
                        __FILE__, __LINE__);
                exit(1);
            }
        }
    } else {
        if (chk_oli_geo(file->geofile).type != FT_OLIL1B) {
            printf("Exit get_oli_angles with error\n");
            fprintf(stderr, "-E- %s line %d: geofile specified must be an Enhanced Meta-Data file (EMD format)\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        // nominal angles
        if (get_oli_angles(file->geofile, file->npix, file->nscan, recnum,
                l1rec->solz, l1rec->sola, l1rec->senz, l1rec->sena, 0) == -1) {
            printf("Exit get_oli_angles with error\n");
            fprintf(stderr, "-E- %s line %d: Unable to compute geometries\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        // band based geometry requested
        if (l1_input->geom_per_band) {
            if (get_oli_angles(file->geofile, file->npix, file->nscan, recnum,
                    l1rec->geom_per_band->solz, l1rec->geom_per_band->sola,
                    l1rec->geom_per_band->senz, l1rec->geom_per_band->sena,
                    l1_input->geom_per_band) == -1) {
                printf("Exit get_oli_angles with error\n");
                fprintf(stderr, "-E- %s line %d: Unable to compute geometries\n",
                        __FILE__, __LINE__);
                exit(1);
            }
        }
    }

    if (lonlat)
        return (LIFE_IS_GOOD);
    
    // int y = 0;
    int y_in_tile = 0;
    
    if (tile_exist != 0) {
        // y = recnum / tileLength;
        y_in_tile = recnum % tileLength;
        if (data->line_num_cached == -1 || recnum < data->line_num_cached || recnum >= (data->line_num_cached + tileLength)) {
            int y = recnum / tileLength;
            data->line_num_cached = y * tileLength;
            for (ib = 0; ib < maxReflBands; ib++) { // read Bands 1-7 && 9
                for (int x = 0; x < imageWidth; x += tileWidth) {
                    if (TIFFReadTile(data->tif[ib], (void*) (data->buf_tile_row + ib * imageWidth * tileLength + x * tileLength), x, y * tileLength, 0, 0) == -1) {
                        fprintf(stderr, "-E- %s line %d: Failed to read Lt, band %d, recnum %d\n",
                            __FILE__, __LINE__, ib, recnum);
                        exit(1);
                    }
                }
            }
        }
    } 
    
    for (ib = 0; ib < file->nbands; ib++) {

        if (tile_exist == 0) {
            if (TIFFReadScanline(data->tif[ib], (void*) data->buf, recnum, 0) == -1) {
                fprintf(stderr, "-E- %s line %d: Failed to read Lt, band %d, recnum %d\n",
                        __FILE__, __LINE__, ib, recnum);
                exit(1);
            }
        } else {
            for (int x = 0; x < imageWidth; x += tileWidth) {
                if ((imageWidth - x) >= tileWidth) {
                    _TIFFmemcpy((void*) (data->buf + x), (void*) (data->buf_tile_row + ib * imageWidth * tileLength + x * tileLength + y_in_tile * tileWidth), tileWidth * 2);
                }                    
                else
                    _TIFFmemcpy((void*) (data->buf + x), (void*) (data->buf_tile_row + ib * imageWidth * tileLength + x * tileLength + y_in_tile * tileWidth), (imageWidth - x) * 2);
            }                 
        }

        l1rec->Fo[ib] = Fobar[ib] * fsol;

        for (ip = 0; ip < file->npix; ip++) {
            ipb = ip * file->nbands + ib;
            if (data->buf[ip] == 0) {
                l1rec->Lt[ipb] = BAD_FLT; // I assume this is outside the projected tile
                l1rec->navwarn[ip] = 1; // so set navigation failure
            } else {
                l1rec->Lt[ipb] = (data->buf[ip] * data->refl_scale[ib] + data->refl_offset[ib]) * l1rec->Fo[ib] / PI;
                // New addition by Sudipta to scale by band specific FPM gain factor
                for (iw = 0; iw < l1_input->xcal_nwave; iw++) {
                    if ((bindex_get(l1_input->xcal_wave[iw]) == ib) &&
                            ((l1_input->xcal_opt[iw] & XCALOLI) != 0)) {
                        isb = (sca_num[ipb] - 1) * file->nbands + ib;
                        l1rec->Lt[ipb] *= xcal_factor[isb];
                        break;
                    }
                }
                // end Sudipta addition
            }
        }
    }

    // read Cirrus Band 9

    ib = 7;

    if (tile_exist == 0) {
        if (TIFFReadScanline(data->tif[ib], (void*) data->buf, recnum, 0) == -1) {
            fprintf(stderr, "-E- %s line %d: Failed to read Lt, band %d, recnum %d\n",
                        __FILE__, __LINE__, ib, recnum);
            exit(1);
        }
    } else {
            for (int x = 0; x < imageWidth; x += tileWidth) {
                if ((imageWidth - x) >= tileWidth) {
                    _TIFFmemcpy((void*) (data->buf + x), (void*) (data->buf_tile_row + ib * imageWidth * tileLength + x * tileLength + y_in_tile * tileWidth), tileWidth * 2);
                }                    
                else
                    _TIFFmemcpy((void*) (data->buf + x), (void*) (data->buf_tile_row + ib * imageWidth * tileLength + x * tileLength + y_in_tile * tileWidth), (imageWidth - x) * 2);
            }                 
    }

    for (ip = 0; ip < file->npix; ip++) {
        if (data->buf[ip] == 0)
            l1rec->rho_cirrus[ip] = BAD_FLT;
        else
            l1rec->rho_cirrus[ip] = (data->buf[ip] * data->refl_scale[ib] + data->refl_offset[ib])
            / cos(l1rec->solz[ip] / RADEG);

    }
    return 0;
}

int closel1_oli(filehandle *file) {
    int ib;
    oli_t* data = (oli_t*) file->private_data;

    // undo what open allocated
    free(data->lat);
    free(data->lon);
    data->lat = data->lon = NULL;

    free(data->buf);
    data->buf = NULL;

    if (tile_exist != 0) {
        free(data->buf_tile_row);
        data->buf_tile_row = NULL;
    }

    GTIFFree(data->gtif);
    data->gtif = NULL;

    for (ib = 0; ib < file->nbands; ib++) {
        XTIFFClose(data->tif[ib]);
    }
    freePrivateData(data);
    file->private_data = NULL;

    if (xcal_factor) free(xcal_factor);

    return 0;
}

