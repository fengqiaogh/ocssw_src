#include <stdlib.h>
#include <stdio.h>

#include <netcdf.h>
#include "l1.h"
#include "l1_l5tm.h"

#include <tiffio.h>
#include <geotiff.h>
#include <xtiffio.h>
#include <geo_normalize.h>
#include <libgen.h>

/* For Landsat 5 TM */

int get_l57tm_nom_angles( char *meta_filename, int32_t npix, int32_t nscan, int32_t iscan,
        float *solz, float *sola, float *senz, float *sena);
int get_l5tm_angles( char *emeta_filename, int32_t npix, int32_t nscan, int32_t iscan,
        float *solz, float *sola, float *senz, float *sena);
int32_t chk_l5tm_geo(char *fname);

static const int itemSize = 500;
static const int maxBands = 7;

typedef struct L5TM_struct {
    int32_t year, doy, msec;
    double sunAzimuth, sunElevation;
    double *lat, *lon;
    double *scale, *offset;
    /* double *refl_scale, *refl_offset; */ // not available for Landsat 5 TM
    TIFF** tif; // file handle for each band
    GTIF* gtif; // geotiff handle for the first file
    GTIFDefn* defn; // geotiff definition structure for first file
    uint8_t* buf; // buffer used to read one scan line from TIFF file
} l5tm_t;


l5tm_t* createPrivateData_l5tm(int numBands) {
    l5tm_t* data = (l5tm_t*)calloc(1, sizeof(l5tm_t));
    if(data == NULL) {
        fprintf(stderr,"-E- %s line %d: unable to allocate private data for L5TM\n",
            __FILE__,__LINE__);
        exit(1);
    }

    data->scale = (double *) malloc(numBands*sizeof(double) );
    data->offset = (double *) malloc(numBands*sizeof(double) );
    if(data->scale==NULL || data->offset==NULL) {
        fprintf(stderr,"-E- %s line %d: unable to allocate scale/offset data for L5TM\n",
            __FILE__,__LINE__);
        exit(1);
    }

    /* data->refl_scale = (double *) malloc(numBands*sizeof(double) );
    data->refl_offset = (double *) malloc(numBands*sizeof(double) );
    if(data->refl_scale==NULL || data->refl_offset==NULL) {
        fprintf(stderr,"-E- %s line %d: unable to allocate reflectance scale/offset data for L5TM\n",
            __FILE__,__LINE__);
        exit(1);
    } */

    data->tif = (TIFF**) calloc(numBands, sizeof(TIFF*) );
    if(data->tif==NULL) {
        fprintf(stderr,"-E- %s line %d: unable to allocate TIFF pointers for L5TM\n",
            __FILE__,__LINE__);
        exit(1);
    }

    data->defn = (GTIFDefn*) malloc(sizeof(GTIFDefn) );
    if(data->defn==NULL) {
        fprintf(stderr,"-E- %s line %d: unable to allocate GEOTIFF definition structure for L5TM\n",
            __FILE__,__LINE__);
        exit(1);
    }

    return data;
}

void freePrivateData_l5tm(l5tm_t* data) {
    free(data->scale);
    free(data->offset);
    /* free(data->refl_scale);
    free(data->refl_offset); */
    free(data->tif);
    free(data->defn);
    free(data);
}

void readNextLine_l5tm(FILE* fp, char* tag, int* i, char* val) {
    char* result;
    char line[itemSize];
    int count;

    result = fgets(line, itemSize, fp);
    if(result == NULL) {
        fprintf(stderr,"-E- %s line %d: unable to read all of the required metadata from L5TM file\n",
            __FILE__,__LINE__);
        exit(1);
    }
    trimBlanks(line);

    count = sscanf(line, "%s = %s", tag, val);
    if(count != 2) {
        // not found so return blank line
        tag[0] = 0;
        *i = 0;
        val[0] = 0;
        return;
    }

    // grab index if it exists
    result = strrchr(tag, '_');
    if(result) {
        result++;
        *i = atoi(result);
    } else {
        *i = 0;
    }

    // get rid of the quotes from val
    if(val[0] == '"')
        val[0] = ' ';
    count = strlen(val) - 1;
    if(val[count] == '"')
        val[count] = 0;
    trimBlanks(val);
}


int read_l5tm_angles(char *file, int32_t npix, int32_t nscan, int32_t iscan, 
                    float *solz, float *sola, float *senz, float *sena) {

    static int  firstCall = 1;
    static int  fileID;

    int    retval;
    int    xid, yid, varID=-1;
    size_t geo_npix, geo_nscan;
    size_t start[3], count[3];
 
    if (firstCall) {

        firstCall = 0;
  
        printf("Reading path angles from %s.\n",file);

        // Open the netcdf file
        retval = nc_open(file, NC_NOWRITE, &fileID);
        if (retval != NC_NOERR) {
            fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
                    __FILE__, __LINE__, file);
            return (-1);
        }

        // Get pixel and scan dimensions
        retval = nc_inq_dimid(fileID, "Pixels", &xid);
        retval = nc_inq_dimid(fileID, "Lines" , &yid);
        retval = nc_inq_dimlen(fileID, xid, &geo_npix);
        retval = nc_inq_dimlen(fileID, yid, &geo_nscan);

        if (retval != NC_NOERR) {
            fprintf(stderr, "-E- %s line %d: Error reading dimensions from %s.\n",
                    __FILE__, __LINE__, file);
            return (-1);
        }
        if (geo_npix != npix || geo_nscan != nscan) {
            fprintf(stderr, "-E- %s line %d: geofile dimensions (%zu,%zu) do not match image dimensions (%d,%d).\n",
                    __FILE__, __LINE__, geo_npix,geo_nscan,npix,nscan);
            return (-1);
        }
 
    }

    start[0] = iscan;
    start[1] = 0;
    count[0] = 1;
    count[1] = npix;

    retval = nc_inq_varid(fileID, "solz", &varID);
    if (retval == NC_NOERR)
        retval = nc_get_vara_float(fileID, varID, start, count, solz);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: Unable to read solz from %s.\n",
                __FILE__, __LINE__, file);
        return (-1);
    }

    retval = nc_inq_varid(fileID, "sola", &varID);
    if (retval == NC_NOERR)
        retval = nc_get_vara_float(fileID, varID, start, count, sola);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: Unable to read sola from %s.\n",
                __FILE__, __LINE__, file);
        return (-1);
    }

    retval = nc_inq_varid(fileID, "senz", &varID);
    if (retval == NC_NOERR)
        retval = nc_get_vara_float(fileID, varID, start, count, senz);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: Unable to read senz from %s.\n",
                __FILE__, __LINE__, file);
        return (-1);
    }

    retval = nc_inq_varid(fileID, "sena", &varID);
    if (retval == NC_NOERR)
        retval = nc_get_vara_float(fileID, varID, start, count, sena);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: Unable to read sena from %s.\n",
                __FILE__, __LINE__, file);
        return (-1);
    }        
 
    return(0);
}


int openl1_l5tm(filehandle *file) {
    int   i;
    FILE *fp;
    char tag[itemSize];
    char val[itemSize];
    char dateStr[32];
    char timeStr[32];
    char fileName[itemSize];

    l5tm_t* data = file->private_data = createPrivateData_l5tm(maxBands);

    if(want_verbose)
        printf("L5TM Level-1B %s\n", file->name );

    /* Open file */
    if ((fp = fopen(file->name, "r")) == NULL) {
        fprintf(stderr,"-E- %s line %d: unable open %s\n",
            __FILE__,__LINE__,file->name);
        exit(1);
    }

    int dateNeeded = 1;
    int timeNeeded = 1;
    int filesNeeded = 1;
    int numLinesNeeded = 1;
    int numSamplesNeeded = 1;
    int sunAzimuthNeeded = 1;
    int sunElevationNeeded = 1;
    int scaleNeeded = 1;
    int offsetNeeded = 1;
    /* int reflScaleNeeded = 1;
    int reflOffsetNeeded = 1; */

    // loop metadata
    while(dateNeeded ||
            timeNeeded  ||
            numLinesNeeded ||
            numSamplesNeeded ||
            filesNeeded ||
            sunAzimuthNeeded ||
            sunElevationNeeded ||
            /* reflScaleNeeded ||
            reflOffsetNeeded || */
            scaleNeeded ||
            offsetNeeded) {

        readNextLine_l5tm(fp, tag, &i, val);

        // skip blank lines
        if(tag[0] == 0)
            continue;

        // get date
        if(!strcmp(tag, "DATE_ACQUIRED")) {
            dateNeeded = 0;
            strcpy(dateStr, val);

        // get time
        } else if(!strcmp(tag, "SCENE_CENTER_TIME")) {
            timeNeeded = 0;
            strcpy(timeStr, val);

        // get num lines
        } else if(!strcmp(tag, "REFLECTIVE_LINES")) {
            numLinesNeeded = 0;
            file->nscan = atoi(val);

        // get num samples
        } else if(!strcmp(tag, "REFLECTIVE_SAMPLES")) {
            numSamplesNeeded = 0;
            file->npix = atoi(val);


        // get band file names
        } else if(!strncmp(tag, "FILE_NAME_BAND_", 15)) {
            if(i == maxBands)
                filesNeeded = 0;
            i--;
            if(i>=0 && i<maxBands) {
                // dirname might destroy input, so we pass it a copy
                char dir[FILENAME_MAX];
                strcpy(dir,file->name);
                strcpy(fileName, dirname(dir));
                strcat(fileName, "/");
                strcat(fileName, val);
                if(want_verbose)
                    printf("L5TM Level-1B Band[%d]:%s\n", i, fileName );
                data->tif[i] = XTIFFOpen(fileName, "r");
                if (!data->tif[i]) {
                    fprintf(stderr,"-E- %s line %d: unable open TIFF file %s\n",
                        __FILE__,__LINE__,fileName);
                    exit(1);
                }
            }

        // get sun azimuth
        } else if(!strcmp(tag, "SUN_AZIMUTH")) {
            sunAzimuthNeeded = 0;
            data->sunAzimuth = atof(val);

        // get sun elevation
        } else if(!strcmp(tag, "SUN_ELEVATION")) {
            sunElevationNeeded = 0;
            data->sunElevation = atof(val);

        // get refl_scale
        } /* else if(!strncmp(tag, "REFLECTANCE_MULT_BAND_", 22)) {
            if(i == maxReflBands)
                reflScaleNeeded = 0;
            i--;
            if(i>=0 && i<maxReflBands) {
                data->refl_scale[i] = atof(val);
            }

        // get refl_offset
        }  else if(!strncmp(tag, "REFLECTANCE_ADD_BAND_", 21)) {
            if(i == maxReflBands)
                reflOffsetNeeded = 0;
            i--;
            if(i>=0 && i<maxReflBands) {
                data->refl_offset[i] = atof(val);
            }

        } */
        // radiance scale 
        else if(!strncmp(tag, "RADIANCE_MULT_BAND_", 19)) {
            if(i == maxBands)
                scaleNeeded = 0;
            i--;
            if(i>=0 && i<maxBands) {
                data->scale[i] = atof(val);
            }

        // get offset
        } else if(!strncmp(tag, "RADIANCE_ADD_BAND_", 18)) {
            if(i == maxBands)
                offsetNeeded = 0;
            i--;
            if(i>=0 && i<maxBands) {
                data->offset[i] = atof(val);
            }
        }

    } // while

    fclose(fp);

    // allocate lat and lon storage
    data->lat = (double *) malloc(file->npix*sizeof(double) );
    data->lon = (double *) malloc(file->npix*sizeof(double) );
    if(data->lat==NULL || data->lon==NULL) {
        fprintf(stderr,"-E- %s line %d: unable to allocate lat/lon data for L5TM\n",
            __FILE__,__LINE__);
        exit(1);
    }

    // only need the GEO TIFF info from one file
    data->gtif = GTIFNew(data->tif[0]);
    if (!data->gtif) {
        fprintf(stderr,"-E- %s line %d: unable open GEOTIFF file %s\n",
            __FILE__,__LINE__,fileName);
        exit(1);
    }

    if(!GTIFGetDefn(data->gtif, data->defn)) {
        fprintf(stderr,"-E- %s line %d: unable populate GEOTIFF defn structure for %s\n",
            __FILE__,__LINE__,fileName);
        exit(1);
    }

    // allocate buffer to hold one scan line
    int size = TIFFScanlineSize(data->tif[0]);
    if(size != file->npix) /* For Landsat 5 TM the data type is Byte */
    {
        fprintf(stderr,"-E- %s line %d: unexpected pixel data size in %s\n",
            __FILE__,__LINE__,fileName);
        exit(1);
    }
    data->buf = (uint8_t*) malloc(size);

    // get date "2013-07-19"
    int year, month, day;
    sscanf(dateStr, "%d-%d-%d", &year, &month, &day);

    // get time "10:41:59.9930740Z"
    int hour, minute;
    double sec;
    sscanf(timeStr, "%d:%d:%lf", &hour, &minute, &sec);

    int isec = (int)sec;
    ymdhms2ydmsec(year, month, day, hour, minute, isec,
                       &data->year, &data->doy, &data->msec);
    sec -= isec;
    data->msec += sec * 1000;

    if(want_verbose) {
        printf("L5TM Start Time: %4d-%02d-%02d %03d %02d:%02d:%f\n",
                year, month, day, data->doy, hour, minute, sec);

        printf("L5TM file has %d bands, %d samples, %d lines\n",
                file->nbands, file->npix, file->nscan );
    }
    strcpy(file->spatialResolution,"30 m");

    return 0;
}

int readl1_l5tm( filehandle *file, int recnum, l1str *l1rec, int lonlat) {
    int ip, ib, ipb;
    l5tm_t* data = (l5tm_t*)file->private_data;

    //  set information about data
    l1rec->npix = file->npix;
    l1rec->l1file->sensorID = file->sensorID;
    l1rec->scantime = yds2unix((int16_t) data->year, (int16_t) data->doy, 
            (double) (data->msec / 1000.0 ));

    //  get lat-lon
    for (ip=0; ip<file->npix; ip++) {
        data->lon[ip] = ip;
        data->lat[ip] = recnum;
        GTIFImageToPCS(data->gtif, data->lon+ip, data->lat+ip);
    }

    if (!GTIFProj4ToLatLong(data->defn, file->npix, data->lon, data->lat) ) {
        fprintf(stderr,"-E- %s line %d: unable reproject points for scan %d\n",
            __FILE__,__LINE__,recnum);
        exit(1);
    }

    for (ip=0; ip<file->npix; ip++) {

        l1rec->pixnum[ip] = ip;

        if ( isnan(data->lat[ip]) )
            data->lat[ip] = -999.0;
        if ( isnan(data->lon[ip]) )
            data->lon[ip] = -999.0;
        l1rec->lat[ip] = data->lat[ip];
        l1rec->lon[ip] = data->lon[ip];
        
        /* printf("lat = %f, lon = %f\n", l1rec->lat[ip], l1rec->lon[ip]); */
        
        if (l1rec->lon[ip] < -181.0 || l1rec->lon[ip] > 181.0 ||
            l1rec->lat[ip] <  -91.0 || l1rec->lat[ip] >  91.0 )
        {
          l1rec->navfail[ip] = 1;
          printf("ERROR: lat = %f, lon = %f\n", l1rec->lat[ip], l1rec->lon[ip]);
        }

    }

    // read path angles if user supplied, else use best guess

    /* For Landsat 5 this IS the only option 
     * if (file->geofile == NULL || file->geofile[0] == 0)
    { */
    if (get_l57tm_nom_angles(file->name,file->npix,file->nscan,recnum,l1rec->solz,l1rec->sola,l1rec->senz,l1rec->sena) == -1) {
        fprintf(stderr, "-E- %s line %d: missing or incompatible geofile\n",
                __FILE__,__LINE__);
        exit(1);
    /* } */ 

    }
    /* not sure if for Landsat 5 enhanced metadata file (EMD) is even available so commented out
     * 
     * else {
        if (chk_l5tm_geo(file->geofile) != FT_L5TML1B) {
            printf("Exit get_oli_angles with error\n");
            fprintf(stderr, "-E- %s line %d: non-existent or incompatible geofile.\n-E- File specified must be an Enhanced Meta-Data file (EMD format)\n",
                    __FILE__,__LINE__);
            exit(1);
        }
        if ( get_oli_angles(file->geofile,file->npix,file->nscan,recnum,l1rec->solz,l1rec->sola,l1rec->senz,l1rec->sena) == -1) {
            printf("Exit get_oli_angles with error\n");
           fprintf(stderr, "-E- %s line %d: missing or incompatible geofile\n",
                    __FILE__,__LINE__);
            exit(1);
        }
    } */

    // exit for lonlat mode, read in only lon, lat, time and solz
    if (lonlat)
        return(LIFE_IS_GOOD);

    for(ib = 0; ib < file->nbands; ib++)
    {
        
        /* For Landsat 5 skip band 6 */
        if (ib == 5)
        {
            if(TIFFReadScanline(data->tif[ib+1], (void*)data->buf, recnum, 0) == -1) {
            fprintf(stderr, "-E- %s line %d: Failed to read Lt, band %d, recnum %d\n",
                    __FILE__,__LINE__, ib, recnum );
            exit(1);
            }
        }
        else
        {
            if(TIFFReadScanline(data->tif[ib], (void*)data->buf, recnum, 0) == -1) {
                fprintf(stderr, "-E- %s line %d: Failed to read Lt, band %d, recnum %d\n",
                        __FILE__,__LINE__, ib, recnum );
                exit(1);
            }
        }

        for (ip=0; ip<file->npix; ip++) {
            ipb = ip*file->nbands+ib;
            if(data->buf[ip] == 0) {
                l1rec->Lt[ipb] = BAD_FLT;   // I assume this is outside the projected tile
                l1rec->navfail[ip] = 1;     // so set navigation failure
            } else
            {
                /* For Landsat 5 skip band 6 */
                if (ib == 5)
                {
                    l1rec->Lt[ipb] = (data->buf[ip] * data->scale[ib+1] + data->offset[ib+1]) / 10.0;
                    /* printf("band = %d, scale = %f, offset = %f\n", ib+1, data->scale[ib+1], data->offset[ib+1]);
                    printf("buf = %d, final = %f\n", data->buf[ip], l1rec->Lt[ipb]); */
                }
                else
                {
                    l1rec->Lt[ipb] = (data->buf[ip] * data->scale[ib] + data->offset[ib]) / 10.0;
                    /* printf("band = %d, scale = %f, offset = %f\n", ib, data->scale[ib], data->offset[ib]);
                    printf("buf = %d, final = %f\n", data->buf[ip], l1rec->Lt[ipb]); */
                }
                
            }
        }
    }

    // read Cirrus Band 9

    /* ib = 8;
    if(TIFFReadScanline(data->tif[ib], (void*)data->buf, recnum, 0) == -1) {
        fprintf(stderr, "-E- %s line %d: Failed to read cirrus band, recnum %d\n",
                __FILE__,__LINE__, recnum );
        exit(1);
    }

    for (ip=0;ip<file->npix; ip++) {
        if(data->buf[ip] == 0)
            l1rec->rho_cirrus[ip] = BAD_FLT;
        else
            l1rec->rho_cirrus[ip] = (data->buf[ip] * data->refl_scale[ib] + data->refl_offset[ib])
                                  / cos(l1rec->solz[ip]/RADEG);
       
    } */

    return 0;
}

int closel1_l5tm(filehandle *file) {
    int ib;
    l5tm_t* data = (l5tm_t*) file->private_data;

    // undo what open allocated
    free(data->lat);
    free(data->lon);
    data->lat = data->lon = NULL;

    free(data->buf);
    data->buf = NULL;

    GTIFFree(data->gtif);
    data->gtif = NULL;

    for(ib=0; ib<file->nbands; ib++) {
        XTIFFClose(data->tif[ib]);
    }
    freePrivateData_l5tm(data);
    file->private_data = NULL;

    return 0;
}
