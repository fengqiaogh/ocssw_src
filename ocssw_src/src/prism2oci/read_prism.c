/*
 * l1_aviris.c
 *
 *  Created on: May 18, 2015
 *      Author: Rick Healy SAIC
 *              NASA-GSFC OBPG
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <proj_api.h>
#include "timeutils.h"
#include "prism.h"
#include "jplaeriallib.h"
#include <math.h>
#include <genutils.h>
#define SKIP -9999
#define BIP 0
#define BIL 1
#define BSQ 2
#define MAXLINESZ 4290

static const int maxBands = 285;
static double *lat, *lon;

int close_prism(prism4ocia_t *data);
prism4ocia_t* open_prism(char *filename, prism4ocia_t **data);
int read_prism(prism4ocia_t *data, int32_t recnum);
char *checkTagLine(char *linein, char* tag);
int checkTagLine_i(char *linein, char* tag);
int readBinScanLine4Ocip_float(float *Lt, int32_t recnum, int32_t npix, double *gain, int nbands, int numBands, int interleave, int swap, FILE *ptr);
void trimBlanks(char* str);
double angular_distance(double lat1, double lon1, double lat2, double lon2);
double deg2rad(double deg);
void l_sun_(int *iyr, int *iday, double *sec, float *sunr, float *rs);

prism4ocia_t* createPrivateData_pr(int numBands, int32_t nscan, int32_t npix) {

    int i;

    prism4ocia_t* data = (prism4ocia_t*) calloc(1, sizeof (prism4ocia_t));
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

void freePrivateData_pr(prism4ocia_t* data) {
    int k;
    free(data->gain);
    free(data->sena);
    free(data->senz);
    free(data->sola);
    free(data->solz);
    free(data->Lt);
    free(data->lon);
    free(data->lat);
    free(data->utc);
    gsl_interp_accel_free(data->spl_acc);

}

prism4ocia_t *open_prism(char *filename, prism4ocia_t **data) {

    int16_t buffer[10];
    FILE *ptr;
    char tag[itemSize];
    char *val0;
    char val[itemSize];
    char *inbasename, *infile;
    static char *infile2;
    int i, j, k, status, pos;
    double *indata;
    float *indataf;
    float tmp;
    int numBands, num;
    char* result;
    char *line, val1[itemSize];
    int count;
    int year, month, day, hour, minute, second;
    float knts2mps = 0.51444444444; // knots to meters per second
    float ft2m = 0.3048; // feet to meters
    prism4ocia_t *temp;
    double sec;

    char projStr[1024];
    static char *dupline;
    int cnt, linelength, itmp;
    char *iptr;

    temp = *data;
    inbasename = getinbasename(filename);
    pos = strlen(inbasename);
    if (pos <= 0) {
        fprintf(stderr, "-E- %s line %d: Not a avalid prism file %s\n",
                __FILE__, __LINE__, filename);
        exit(-1);
    }

    cdata_(); //  initialize global FORTRAN common block data for l_sun call

    //prism4ocia_t* data = file->private_data = createPrivateData_pr(maxBands,file->nscan,file->npix);
    *data = (prism4ocia_t*) malloc(sizeof (prism4ocia_t));
    if (*data == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to allocate private data for Prism\n",
                __FILE__, __LINE__);
        exit(1);
    }
    line = (char *) malloc(MAXLINESZ * sizeof (char));
    (*data)->wave = (double *) malloc(maxBands * sizeof (double));
    (*data)->fwhm = (double *) malloc(maxBands * sizeof (double));
    (*data)->gain = (double *) malloc(maxBands * sizeof (double));
    if ((*data)->wave == NULL || (*data)->fwhm == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to allocate scale/offset data for aviris\n",
                __FILE__, __LINE__);
        exit(1);
    }

    temp = *data;

    sscanf(inbasename, "prm%4d%2d%2dt%2d%2d%2d", &year, &month, &day, &hour, &minute, &second);

    temp->month = month;
    temp->day = day;
    temp->hour = hour;
    temp->min = minute;
    temp->sec = second;

    ymdhms2ydmsec(year, month, day, hour, minute, second,
            &temp->year, &temp->doy, &temp->msec);


    temp->stime = ymds2unix(year, month, day, hour * 3600.0 + minute * 60.0 + second);

    printf("Date of prism flight: %s\n", unix2isodate(temp->stime, 'G'));

    if ((ptr = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to open %s\n",
                __FILE__, __LINE__, filename);
        exit(-1);
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
    int waveLengthNeeded = 1;
    int fwhmNeeded = 1;

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
            waveLengthNeeded ||
            fwhmNeeded ||
            interleaveNeeded) {

        result = fgets(line, MAXLINESZ, ptr);
        if (result == NULL) {
            fprintf(stderr, "-E- %s line %d: unable to read all of the required metadata from prism file\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        trimBlanks(line);

        if ((itmp = checkTagLine_i(line, "lines")) && itmp > 0) {
            numLinesNeeded = 0;
            temp->nscan = itmp; //strtol(val0,&iptr,10);
            printf("lines=%d\n", temp->nscan);
        }
        if ((itmp = checkTagLine_i(line, "samples")) && itmp > 0) {
            numSamplesNeeded = 0;
            temp->npix = itmp; //strtol(val0,&iptr,10);
            printf("samples=%d\n", temp->npix);
        }
        if ((val0 = checkTagLine(line, "Alt"))) {
            altNeeded = 0;
            temp->alt = atof(val0) * ft2m;
            printf("Altitude=%lf\n", temp->alt);
        }
        if ((val0 = checkTagLine(line, "EndTime"))) {
            EndTimeNeeded = 0;
            sscanf(val0, "%2d%2d", &hour, &minute);
            printf("End hour=%d minute=%d\n", hour, minute);
            temp->etime = ymds2unix(year, month, day, hour * 3600.0 + minute * 60.0);
            printf("End Time=%s\n", unix2isodate(temp->etime, 'G'));
        }
        if ((itmp = checkTagLine_i(line, "bands")) && itmp > 0) {
            bandsNeeded = 0;
            numBands = itmp;
            temp->numBands = numBands;
            if (numBands > maxBands) {
                fprintf(stderr, "-E- %s line %d: number of bands (%d) from prism file > maxBands (%d)\n",
                        __FILE__, __LINE__, numBands, maxBands);
                exit(1);
            }
            printf("numBands=%d\n", temp->numBands);
        }
        if ((val0 = checkTagLine(line, "interleave"))) {
            interleaveNeeded = 0;
            if (strstr(val0, "bip")) {
                temp->interleave = BIP;
            } else if (strstr(val0, "bil")) {
                temp->interleave = BIL;
            } else {
                fprintf(stderr, "Interleave = %s is not supported\n", val0);
                exit(1);
            }
            printf("Interleave=%d\n", temp->interleave);
        }
        if ((val0 = checkTagLine(line, "wavelength ="))) {
            waveLengthNeeded = 0;
            i = 0;
            result = strtok(line, "wavelength = {");
            trimBlanks(result);
            tmp = atof(result);

            while (i < maxBands && result) {
                temp->wave[i] = 1000.0 * tmp;
                printf("i=%d wave=%f\n", i, temp->wave[i]);
                result = strtok(NULL, ",");
                if (result) {
                    trimBlanks(result);
                    tmp = atof(result);
                    i++;
                }
            }


        }
        if ((val0 = checkTagLine(line, "fwhm"))) {
            fwhmNeeded = 0;
            i = 0;
            result = strtok(line, "fwhm = {");
            trimBlanks(result);
            tmp = atof(result);

            while (i < maxBands && result) {
                temp->fwhm[i] = 1000.0 * tmp;
                printf("i=%d fwhm=%f\n", i, temp->fwhm[i]);
                result = strtok(NULL, ",");
                if (result) {
                    trimBlanks(result);
                    tmp = atof(result);
                    i++;
                }
            }
        }
        if ((val0 = checkTagLine(line, "map info"))) {
            cnt = 0;
            linelength = strlen(line);
            if (dupline) free(dupline);
            if ((dupline = (char *) malloc(linelength * sizeof (char))) == NULL) {
                fprintf(stderr,
                        "-E- %s line %d: Memory allocation failure.\n",
                        __FILE__, __LINE__);
                exit(-1);
            }
            strcpy(dupline, line);
            result = strtok(dupline, ",");
            while (result) {
                switch (cnt) {
                case 3:
                    temp->easting = atof(result);
                    eastingNeeded = 0;
                    break;
                case 4:
                    temp->northing = atof(result);
                    northingNeeded = 0;
                    break;
                case 5:
                    temp->pixelSize = atof(result);
                    pixelSizeNeeded = 0;
                    break;
                case 7:
                    temp->utmZone = strtol(result, &iptr, 10);
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
                temp->rotation = atof(val0);
            } else {
                printf("Rotation angle expected in line: %s\n", val0);
                exit(-1);
            }
            printf("Rotation/easting/northing/pixsize/utmzone=%f/%lf/%lf/%lf/%d\n", temp->rotation, temp->easting, temp->northing, temp->pixelSize, temp->utmZone);
        }
        //        if((val0=checkTagLine(line,"map info"))) {
        //            sscanf(val0,"%*[^,], %*f, %*f, %lf, %lf, %lf, %lf, %d, %*[^,], %*[^,], %*[^,], %s}",&temp->easting,&temp->northing,&temp->pixelSize,&temp->pixelSize,&temp->utmZone,val1);
        //            if((val0=checknspTagLine(line,"rotation"))) {
        //                rotationAngNeeded = 0;
        //                rotation = atof(val0);
        //                if (rotation > 45)
        //                    temp->eastbyscan = -1;
        //                 else if (rotation < -45)
        //                     temp->eastbyscan = 1;
        //                 else
        //                     temp->eastbyscan = 0;
        //            } else {
        //                printf("Rotation angle expected in line: %s\n",val0);
        //                exit(-1);
        //            }
        //            printf("Rotation/easting/northing/pixsize/utmzone=%f/%lf/%lf/%lf/%d\n",rotation,temp->easting,temp->northing,temp->pixelSize,temp->utmZone);
        //            pixelSizeNeeded = 0;
        //            northingNeeded = 0;
        //            eastingNeeded = 0;
        //            utmZoneNeeded = 0;
        //
        //        }


    }

    fclose(ptr);


    // Get the sensor and solar temp

    temp->sena = (float *) malloc(temp->npix * sizeof (float));
    temp->senz = (float *) malloc(temp->npix * sizeof (float));
    temp->solz = (float *) malloc(temp->npix * sizeof (float));
    temp->sola = (float *) malloc(temp->npix * sizeof (float));
    temp->utc = (float *) malloc(temp->npix * sizeof (float));
    temp->lon = (double *) malloc(temp->npix * sizeof (double));
    temp->lat = (double *) malloc(temp->npix * sizeof (double));
    temp->Lt = (float *) malloc(temp->numBands * temp->npix * sizeof (float));
    if (temp->sena == NULL || temp->senz == NULL || temp->sola == NULL || temp->solz == NULL || temp->Lt == NULL || temp->lon == NULL || temp->lat == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to allocate sensor and solar angle data for AVIRIS\n",
                __FILE__, __LINE__);
        exit(1);
    }

    for (i = 0; i < numBands; i++) {
        temp->gain[i] = 1.0;
    }

    // free(infile);
    infile = malloc((pos + strlen("_rdn_ort")) * sizeof (char));
    strcpy(infile, inbasename);
    strcat(infile, "_rdn_ort");
    printf("Opening prism image file %s\n", infile);

    if ((temp->av_fp = fopen(infile, "rb")) == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to open %s\n",
                __FILE__, __LINE__, infile);
        exit(-1);

    }

    sprintf(projStr, "+proj=utm +zone=%d +ellps=WGS84 +datum=WGS84 +units=m +no_defs ",
            temp->utmZone);
    //    printf("RJH: projStr=%s\n",projStr);
    if (!(temp->pj_ortho = pj_init_plus(projStr))) {
        printf("Error - prism first projection failed to init\n");
        exit(1);
    }
    if (!(temp->pj_latlong = pj_latlong_from_proj(temp->pj_ortho))) {
        fprintf(stderr, "-E- %s line %d: prism latlon projection failed to init\n",
                __FILE__, __LINE__);
        exit(1);
    }

    lat = (double *) malloc(temp->npix * sizeof (double));
    lon = (double *) malloc(temp->npix * sizeof (double));

    return (temp);

}

int read_prism(prism4ocia_t *data, int recnum)
/*
 *  fill standard record with L1B line of data
 */ {
    int status, min_msec = 86401 * 1000, max_msec = -1;

    static double last_good_hour = 18;
    static int firstCall = 1;
    float rel_sec;
    double sec, pos[3];
    float epos[3], sunpos[3];
    int16_t year, doy, month, day, hour, minute;
    double dist;
    double secondOfDay;
    float longitude, latitude, sunDist;
    int i;
    int npix = data->npix, ip, ib, ipb;


    //  set information about data
    npix = data->npix;

    //    *(l1rec->year) = data->year;
    //    *(l1rec->day)  = data->doy;
    //    *(l1rec->msec) = data->msec + (data->emsec - data->msec)*recnum/(data->nscan-1);
    data->scantime = data->stime + (data->etime - data->stime) * recnum / (data->nscan - 1);
    unix2ymds(data->scantime, &year, &month, &day, &secondOfDay);
    data->year = year;
    data->month = month;
    data->day = day;
    data->hour = secondOfDay / 3600;
    data->min = (secondOfDay - data->hour * 3600) / 60;
    data->sec = secondOfDay - data->hour * 3600 - data->min * 60;
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
    unix2yds(data->scantime, &year, &doy, &secondOfDay);

    int32_t iyear, idoy;
    iyear = (int32_t) year;
    idoy = (int32_t) doy;
    l_sun_(&iyear, &idoy, &secondOfDay, sunpos, &sunDist); // get position vector for the sun

    for (i = 0; i < 3; i++) {
        sunpos[i] *= 1.496e8; //convert to km for call to get_zenaz
        epos[i] = pos[i];
    }

    for (ip = 0; ip < npix; ip++) {


        if (isnan(lat[ip])) lat[ip] = SKIP;
        if (isnan(lon[ip])) lon[ip] = SKIP;
        data->lat[ip] = lat[ip];
        data->lon[ip] = lon[ip];
        //printf("ip=%d scan=%d lat=%f lon=%f\n",ip,recnum,lat[ip],lon[ip]);


        get_zenaz(epos, lon[ip], lat[ip], &data->senz[ip], &data->sena[ip]);
        get_zenaz(sunpos, lon[ip], lat[ip], &data->solz[ip], &data->sola[ip]);

        //printf("RJH: %d %d senz=%f sena=%f solz=%f sola=%f\n",recnum, ip, l1rec->senz[ip],l1rec->sena[ip], l1rec->solz[ip],l1rec->sola[ip]);


    }

    readBinScanLine_float(data->Lt, recnum, data->npix, data->gain, data->numBands, data->numBands, data->interleave, 0, data->av_fp);

    return (0);
}

void prism4ocia_proj4_convert(prism4ocia_t *data, int32_t numPoints, double *x, double *y) {
    int i;
    if (pj_transform(data->pj_ortho, data->pj_latlong, numPoints, 1, x, y, NULL)) {
        fprintf(stderr, "-E- %s line %d: AVIRIS proj4 transformation blew up\n",
                __FILE__, __LINE__);
        exit(1);
    }
    for (i = 0; i < numPoints; i++) {
        x[i] *= RAD_TO_DEG;
        y[i] *= RAD_TO_DEG;
    }
}

int close_prism(prism4ocia_t* data) {
    int ib;

    // undo what open allocated

    freePrivateData_pr(data);
    free(data);

    return 0;
}

float getValidOrcaAngle(float *ang, int32_t npix, int32_t skip, float *fillangle) {
    int32_t i;
    float angle = *fillangle;

    for (i = 0; i < npix && ang[i] <= skip; i++)
        angle = ang[i];

    *fillangle = angle;
    return (angle);

}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::  Compute the angular distance                                  :*/

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

double angular_distance(double lat1, double lon1, double lat2, double lon2) {
    double theta, dist;
    theta = lon1 - lon2;
    dist = sin(deg2rad(lat1)) * sin(deg2rad(lat2)) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * cos(deg2rad(theta));
    dist = dist > 1 ? 1 : dist;
    dist = dist<-1 ? -1 : dist;
    return (acos(dist));
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::  This function converts decimal degrees to radians             :*/

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
double deg2rad(double deg) {
    return (deg * OEL_PI / 180);
}

void prism_proj4_convert(prism4ocia_t * data, int numPoints, double *x, double *y) {
    int i;
    if (pj_transform(data->pj_ortho, data->pj_latlong, numPoints, 1, x, y, NULL)) {
        fprintf(stderr, "-E- %s line %d: prism proj4 transformation blew up\n",
                __FILE__, __LINE__);
        exit(1);
    }
    for (i = 0; i < numPoints; i++) {
        x[i] *= RAD_TO_DEG;
        y[i] *= RAD_TO_DEG;
    }
}

