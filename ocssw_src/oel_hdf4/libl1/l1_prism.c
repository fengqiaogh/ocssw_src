/*
 * l1_prism.c
 *
 *  Portable Remote Imaging SpectroMeter (PRISM)
 *
 *  Created on: June 11, 2015
 *      Author: Rick Healy SAIC
 *              NASA-GSFC OBPG
 */
#include "l1.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "jplaeriallib.h"
#include "prism.h"
#include <libnav.h>

#define SKIP -9999

static const int maxBands = 285;
static double *lat, *lon;

prism_t* createPrivateData_pr(int numBands, int32_t nscan, int32_t npix) {

    prism_t* data = (prism_t*) calloc(1, sizeof (prism_t));
    if (data == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to allocate private data for prism\n",
                __FILE__, __LINE__);
        exit(1);
    }

    data->wave = (double *) malloc(numBands * sizeof (double));
    data->fwhm = (double *) malloc(numBands * sizeof (double));
    data->gain = (double *) malloc(numBands * sizeof (double));
    if (data->wave == NULL || data->fwhm == NULL || data->gain == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to allocate scale/offset data for prism\n",
                __FILE__, __LINE__);
        exit(1);
    }

    return data;
}

int openl1_prism(filehandle *file) {

    FILE *ptr;
    char *val0;
    char *inbasename, *infile;
    int i, pos;
    int numBands;
    char* result;
    char line[itemSize];
    int year, month, day, hour, minute, second;
    //float knts2mps = 0.51444444444; // knots to meters per second
    float ft2m = 0.3048; // feet to meters

    char projStr[1024];
    static char *dupline;
    int cnt, linelength;

    printf("file=%s\n", file->name);
    inbasename = getinbasename(file->name);
    printf("Basename=%s\n", inbasename);
    pos = strlen(inbasename);
    if (pos <= 0) {
        fprintf(stderr, "-E- %s line %d: Not a avalid prism file %s\n",
                __FILE__, __LINE__, file->name);
        return 1;
    }

    sscanf(inbasename, "prm%4d%2d%2dt%2d%2d%2d", &year, &month, &day, &hour, &minute, &second);

    prism_t* data = file->private_data = createPrivateData_pr(maxBands, file->nscan, file->npix);

    data->stime = ymds2unix(year, month, day, hour * 3600.0 + minute * 60.0 + second);

    printf("Date of prism flight: %s\n", unix2isodate(data->stime, 'G'));

    if ((ptr = fopen(file->name, "r")) == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to open %s\n",
                __FILE__, __LINE__, file->name);
        return 1;
    }

    int numLinesNeeded = 1;
    int numSamplesNeeded = 1;
    int bandsNeeded = 1;
    int utmZoneNeeded = 1;
    int eastingNeeded = 1;
    int northingNeeded = 1;
    int pixelSizeNeeded = 1;
    int interleaveNeeded = 1;
    int rotationAngNeeded = 1;
    int EndTimeNeeded = 1;
    int altNeeded = 1;

    // loop metadata
    while (numLinesNeeded ||
            numSamplesNeeded ||
            bandsNeeded ||
            pixelSizeNeeded ||
            utmZoneNeeded ||
            eastingNeeded ||
            northingNeeded ||
            rotationAngNeeded ||
            EndTimeNeeded ||
            altNeeded ||
            interleaveNeeded) {

        result = fgets(line, itemSize, ptr);
        if (result == NULL) {
            fprintf(stderr, "-E- %s line %d: unable to read all of the required metadata from prism file\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        trimBlanks(line);

        if ((val0 = checkTagLine(line, "lines"))) {
            numLinesNeeded = 0;
            file->nscan = atoi(val0);
            data->nscan = file->nscan;
            printf("lines=%d\n", data->nscan);
        }
        if ((val0 = checkTagLine(line, "samples"))) {
            numSamplesNeeded = 0;
            file->npix = atoi(val0);
            data->npix = file->npix;
            printf("samples=%d\n", data->npix);
        }
        if ((val0 = checkTagLine(line, "Alt"))) {
            altNeeded = 0;
            data->alt = atof(val0) * ft2m;
            printf("Altitude=%lf\n", data->alt);
        }
        if ((val0 = checkTagLine(line, "EndTime"))) {
            EndTimeNeeded = 0;
            sscanf(val0, "%2d%2d", &hour, &minute);
            printf("End hour=%d minute=%d\n", hour, minute);
            data->etime = ymds2unix(year, month, day, hour * 3600.0 + minute * 60.0);
            printf("End Time=%s\n", unix2isodate(data->etime, 'G'));
        }
        if ((val0 = checkTagLine(line, "bands"))) {
            bandsNeeded = 0;
            numBands = atoi(val0);
            data->numBands = numBands;
            if (numBands > maxBands) {
                fprintf(stderr, "-E- %s line %d: number of bands (%d) from prism file > maxBands (%d)\n",
                        __FILE__, __LINE__, numBands, maxBands);
                exit(1);
            }
            printf("numBands=%d\n", data->numBands);
        }
        if ((val0 = checkTagLine(line, "interleave"))) {
            interleaveNeeded = 0;
            if (strstr(val0, "bip")) {
                data->interleave = BIP;
            } else if (strstr(val0, "bil")) {
                data->interleave = BIL;
            } else {
                fprintf(stderr, "Interleave = %s is not supported\n", val0);
                exit(1);
            }
            printf("Interleave=%d\n", data->interleave);
        }
        if ((val0 = checkTagLine(line, "map info"))) {
            cnt = 0;
            linelength = strlen(line);
            if (dupline) free(dupline);
            if ((dupline = (char *) malloc(linelength * sizeof (char))) == NULL) {
                fprintf(stderr,
                        "-E- %s line %d: Memory allocation failure.\n",
                        __FILE__, __LINE__);
                return (0);
            }
            strcpy(dupline, line);
            result = strtok(dupline, ",");
            while (result) {
                switch (cnt) {
                case 3:
                    data->easting = atof(result);
                    eastingNeeded = 0;
                    break;
                case 4:
                    data->northing = atof(result);
                    northingNeeded = 0;
                    break;
                case 5:
                    data->pixelSize = atof(result);
                    pixelSizeNeeded = 0;
                    break;
                case 7:
                    data->utmZone = atoi(result);
                    utmZoneNeeded = 0;
                    break;
                default:
                    break;
                }
                printf(">>%d) %s\n", cnt, result);
                cnt++;
                result = strtok(NULL, ",");
            }
            if ((val0 = checknspTagLine(line, "rotation"))) {
                rotationAngNeeded = 0;
                data->rotation = atof(val0);
            } else {
                printf("Rotation angle expected in line: %s\n", val0);
                exit(-1);
            }
            printf("Rotation/easting/northing/pixsize/utmzone=%f/%lf/%lf/%lf/%d\n", data->rotation, data->easting, data->northing, data->pixelSize, data->utmZone);
        }
        //        if((val0=checkTagLine(line,"map info"))) {
        //            sscanf(val0,"%*[^,], %*f, %*f, %lf, %lf, %lf, %lf, %d, %*[^,], %*[^,], %*[^,], %s}",&data->easting,&data->northing,&data->pixelSize,&data->pixelSize,&data->utmZone,val1);
        //            if((val0=checknspTagLine(line,"rotation"))) {
        //                rotationAngNeeded = 0;
        //                rotation = atof(val0);
        //                if (rotation > 45)
        //                    data->eastbyscan = -1;
        //                 else if (rotation < -45)
        //                     data->eastbyscan = 1;
        //                 else
        //                     data->eastbyscan = 0;
        //            } else {
        //                printf("Rotation angle expected in line: %s\n",val0);
        //                exit(-1);
        //            }
        //            printf("Rotation/easting/northing/pixsize/utmzone=%f/%lf/%lf/%lf/%d\n",rotation,data->easting,data->northing,data->pixelSize,data->utmZone);
        //            pixelSizeNeeded = 0;
        //            northingNeeded = 0;
        //            eastingNeeded = 0;
        //            utmZoneNeeded = 0;
        //
        //        }


    }

    fclose(ptr);


    // Get the sensor and solar data
    data->sena = (double **) malloc(data->nscan * sizeof (double*));
    data->senz = (double **) malloc(data->nscan * sizeof (double*));
    data->solz = (double **) malloc(data->nscan * sizeof (double*));
    data->sola = (double **) malloc(data->nscan * sizeof (double*));
    data->utc = (double **) malloc(data->nscan * sizeof (double*));
    if (data->sena == NULL || data->senz == NULL || data->sola == NULL || data->solz == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to allocate sensor and solar angle data for prism\n",
                __FILE__, __LINE__);
        exit(1);
    }

    for (i = 0; i < data->nscan; i++) {
        data->sena[i] = (double *) malloc(data->npix * sizeof (double));
        data->senz[i] = (double *) malloc(data->npix * sizeof (double));
        data->sola[i] = (double *) malloc(data->npix * sizeof (double));
        data->solz[i] = (double *) malloc(data->npix * sizeof (double));
        data->utc[i] = (double *) malloc(data->npix * sizeof (double));
        if (data->sena[i] == NULL || data->senz[i] == NULL || data->sola[i] == NULL || data->solz[i] == NULL || data->utc[i] == NULL) {
            fprintf(stderr, "-E- %s line %d: unable to allocate sensor and solar angle data for prism\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    }

    for (i = 0; i < numBands; i++) {
        *(data->gain + i) = 1.0;
    }

    // free(infile);
    infile = malloc((pos + strlen("_rdn_ort")) * sizeof (char));
    strcpy(infile, inbasename);
    strcat(infile, "_rdn_ort");
    printf("Opening prism image file %s\n", infile);

    if ((data->av_fp = fopen(infile, "rb")) == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to open %s\n",
                __FILE__, __LINE__, infile);
        return 1;

    }

    PJ *pj;
    sprintf(projStr, "+proj=utm +zone=%d +ellps=WGS84 +datum=WGS84 +units=m +no_defs ",
            data->utmZone);

    // init the proj4 projections
    pj = proj_create_crs_to_crs(PJ_DEFAULT_CTX,
                               projStr,
                               "+proj=longlat +ellps=WGS84 +datum=WGS84",
                               NULL);
    if(pj == NULL) {
        printf("Error - prism first PROJ projection failed to init\n");
        exit(1);
    }
    data->pj = proj_normalize_for_visualization(PJ_DEFAULT_CTX, pj);
    if(data->pj == NULL) {
        printf("Error - prism visualization PROJ projection failed to init\n");
        exit(1);
    }
    proj_destroy(pj);

    lat = (double *) malloc(file->npix * sizeof (double));
    lon = (double *) malloc(file->npix * sizeof (double));

    return (0);

}

int readl1_prism(filehandle *file, int recnum, l1str *l1rec, int lonlat)
/*
 *  fill standard record with L1B line of data
 */ {
    static int firstCall = 1;
    double pos[3];
    float epos[3], sunpos[3];
    int16_t year, doy;
    double secondOfDay;
    float longitude, latitude, sunDist;
    int i;
    int npix = file->npix, ip;
    prism_t* data = (prism_t*) file->private_data;

    if (firstCall) {
        if (want_verbose)
            printf("file->nbands = %d\n", (int) file->nbands);
        firstCall = 0;

        for (ip = 0; ip < npix; ip++) {
            l1rec->pixnum[ip] = ip;
            l1rec->flags[ip] = 0;
        }
    }

    //  set information about data
    npix = file->npix;
    l1rec->npix = file->npix;
    l1rec->alt = data->alt;

    //    *(l1rec->year) = data->year;
    //    *(l1rec->day)  = data->doy;
    //    *(l1rec->msec) = data->msec + (data->emsec - data->msec)*recnum/(data->nscan-1);
    l1rec->scantime = data->stime + (data->etime - data->stime) * recnum / (data->nscan - 1);

    //  get lat-lon
    for (ip = 0; ip < npix; ip++) {
        //rotate pixel to project onto map
        lon[ip] = data->easting + ip * cos(deg2rad(data->rotation)) * data->pixelSize + recnum * sin(deg2rad(data->rotation)) * data->pixelSize; // starts in upper left corner
        lat[ip] = data->northing + ip * sin(deg2rad(data->rotation)) * data->pixelSize - recnum * cos(deg2rad(data->rotation)) * data->pixelSize;
    }

    prism_proj4_convert(data, npix, lon, lat);

    if (lat[npix / 2] > SKIP)
        latitude = lat[npix / 2];
    else {
        fprintf(stderr, "-E- %s line %d: Don't have sensor latitude for geometry calculation\n",
                __FILE__, __LINE__);
        exit(1);
    }

    if (lon[npix / 2] > SKIP)
        longitude = lon[npix / 2];
    else {
        fprintf(stderr, "-E- %s line %d: Don't have sensor longitude for geometry calculation\n",
                __FILE__, __LINE__);
        exit(1);
    }

    getPosVec(latitude, longitude, data->alt, pos); // get position vector of sensor
    unix2yds(l1rec->scantime, &year, &doy, &secondOfDay);

    int32_t iyear, idoy;
    iyear = (int32_t) year;
    idoy = (int32_t) doy;
    l_sun_(&iyear, &idoy, &secondOfDay, sunpos, &sunDist); // get position vector for the sun

    for (i = 0; i < 3; i++) {
        sunpos[i] *= 1.496e8; //convert to km for call to get_zenaz
        epos[i] = pos[i];
    }

    for (ip = 0; ip < npix; ip++) {

        l1rec->pixnum[ip] = ip;

        if (isnan(lat[ip])) lat[ip] = SKIP;
        if (isnan(lon[ip])) lon[ip] = SKIP;
        l1rec->lat[ip] = lat[ip];
        l1rec->lon[ip] = lon[ip];
        //printf("ip=%d scan=%d lat=%f lon=%f\n",ip,recnum,lat[ip],lon[ip]);

        if (l1rec->lon[ip] < -181.0 || l1rec->lon[ip] > 181.0 ||
                l1rec->lat[ip] < -91.0 || l1rec->lat[ip] > 91.0)
            l1rec->navfail[ip] = 1;

        get_zenaz(epos, lon[ip], lat[ip], &l1rec->senz[ip], &l1rec->sena[ip]);
        get_zenaz(sunpos, lon[ip], lat[ip], &l1rec->solz[ip], &l1rec->sola[ip]);

        //printf("RJH: %d %d senz=%f sena=%f solz=%f sola=%f\n",recnum, ip, l1rec->senz[ip],l1rec->sena[ip], l1rec->solz[ip],l1rec->sola[ip]);


    }

    readBinScanLine_float(l1rec->Lt, recnum, l1rec->npix, data->gain, l1rec->l1file->nbands, data->numBands, data->interleave, 0, data->av_fp);

    return (LIFE_IS_GOOD);
}

void prism_proj4_convert(prism_t *data, int numPoints, double *x, double *y) {
    int i;
    PJ_COORD c, c_out;
    
    // set default z and t
    c.xyzt.z = 0.0;
    c.xyzt.t = HUGE_VAL;
        
    for (i = 0; i < numPoints; i++) {
        c.xy.x = x[i];
        c.xy.y = y[i];
        c_out = proj_trans(data->pj, PJ_FWD, c);
        x[i] = c_out.xy.x;
        y[i] = c_out.xy.y;
    }
}

int closel1_prism(filehandle *file) {
    prism_t* data = (prism_t*) file->private_data;

    // undo what open allocated

    freePrivateData_pr(data);
    free(file->private_data);

    return 0;
}

void freePrivateData_pr(prism_t* data) {
    int k;
    free(data->gain);
    for (k = 0; k < data->nscan; k++) {
        free(data->sena[k]);
        free(data->senz[k]);
        free(data->sola[k]);
        free(data->solz[k]);
        free(data->utc[k]);
    }
    free(data->sena);
    free(data->senz);
    free(data->sola);
    free(data->solz);
    free(data->utc);
    gsl_interp_accel_free(data->spl_acc);

}

