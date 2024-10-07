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
#include <proj.h>
#include "timeutils.h"
#include "aviris.h"
#include "jplaeriallib.h"
#include <math.h>
#include <genutils.h>
#include <libnav.h>


#define SKIP -9999


static const int maxBands = 224;

void freePrivateData_ocia(aviris4ocia_t* data) {
    free(data->wave);
    free(data->fwhm);
    free(data->gain);
    free(data->sena);
    free(data->senz);
    free(data->sola);
    free(data->solz);
    free(data->utc);
    gsl_interp_accel_free(data->spl_acc);

}

float getValidAngleB(float *ang, int32_t npix, int32_t skip, float *fillangle) {
    int32_t i;
    float angle = *fillangle;

    for (i = 0; i < npix && ang[i] <= skip; i++)
        angle = ang[i];

    *fillangle = angle;
    return (angle);

}

void aviris4ocia_proj4_convert(aviris4ocia_t *data, int32_t numPoints, double *x, double *y) {
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

void get_nav_data(char* navfile, int32_t nscans, int32_t npix, aviris4ocia_t* data) {

    FILE* ptr;
    double sec;
    char line[itemSize];
    float sunpos[3];
    float latitude, longitude;
    double secondOfDay;
    float sunDist;
    int i, k, ip;
    int16_t year, hour, minute, doy;
    char* result;
    int32_t iyear, idoy;
    double *range, *lon, *lat;
    double c0, c1, d0, d1, cov00, cov01, cov11, chisq, dth = 0.00087, th;
    float xlon, xlat, vel[3], pos[3], pview[3], att[3] = {0, 0, 0}, coef[10],
            *smat[3];
    float sena, senz, sola, solz;

    if (nscans < 2) {
        fprintf(stderr,
                "-E- %s line %d: Navigation needs at least 2 scan lines\n",
                __FILE__, __LINE__);
        exit(-1);
    }
    if ((ptr = fopen(navfile, "r")) == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to open %s\n", __FILE__,
                __LINE__, navfile);
        exit(-1);
    }


    range = (double*) malloc(nscans * sizeof (double));
    lon = (double*) malloc(nscans * sizeof (double));
    lat = (double*) malloc(nscans * sizeof (double));
    for (i = 0; i < 3; i++)
        smat[i] = (float*) malloc(3 * sizeof (float));
    i = 0;
    while ((result = fgets(line, itemSize, ptr)) && i < nscans) {
        if (result == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: unable to read all of the navigation data from AVIRIS nav file\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        trimBlanks(line);
        sscanf(line, "%f %f %lf %lf ", &data->utc[i], &data->alt[i], &lat[i],
                &lon[i]);
        range[i] = i;
        i++;
    }
    // smooth the digitized lon/lat's
    gsl_fit_linear(range, 1, lon, 1, nscans, &c0, &c1, &cov00, &cov01, &cov11,
            &chisq);
    gsl_fit_linear(range, 1, lat, 1, nscans, &d0, &d1, &cov00, &cov01, &cov11,
            &chisq);
    //Get the velocity vector
    nav_get_vel(d0, d1, c0, c1, vel);
    for (i = 0; i < nscans; i++) {
        longitude = c0 + c1 * i;
        latitude = d0 + d1 * i;
        hour = (int) (data->utc[i]);
        minute = (data->utc[i] - hour) * 60;
        sec = ((data->utc[i] - hour) * 60 - minute) * 60;
        data->scantime = ymds2unix(data->year, data->month, data->day,
                (hour * 3600. + minute * 60. + sec));
        data->hour = hour;
        data->min = minute;
        data->sec = sec;
        //printf("Date:utc=%f sec=%f minute=%d (utc-hour)x60=%d \n",data->utc[i],data->sec,data->min,data->hour );
        //getPosVecR(latitude,longitude, data->alt[i], pos); // get position vector of sensor
        unix2yds(data->scantime, &year, &doy, &secondOfDay);
        nav_get_pos(latitude, longitude, data->alt[i], pos);
        nav_gd_orient(pos, vel, att, smat, coef);
        iyear = (int32_t) year;
        idoy = (int32_t) doy;
        l_sun_(&iyear, &idoy, &secondOfDay, sunpos, &sunDist); // get position vector for the sun
        for (k = 0; k < 3; k++) {
            sunpos[k] *= 1.496e8; //convert to km for call to get_zenaz
            //epos[k]    = pos[k];
        }
        //get_zenaz  (epos, longitude, latitude, &data->senz[i*npix], &data->sena[i*npix]);
        //get_zenaz(sunpos, longitude, latitude, &data->solz[i*npix], &data->sola[i*npix]);
        for (ip = 0; ip < npix; ip++) {
            th = (ip - (npix - 1) / 2) * dth;
            pview[0] = 0;
            pview[1] = -sin(th);
            pview[2] = cos(th);
            nav_get_geonav(sunpos, pos, pview, coef, smat, &xlon, &xlat, &solz,
                    &sola, &senz, &sena);
            data->lon[i * npix + ip] = xlon;
            data->lat[i * npix + ip] = xlat;
            data->senz[i * npix + ip] = senz;
            data->sena[i * npix + ip] = sena;
            data->solz[i * npix + ip] = solz;
            data->sola[i * npix + ip] = sola;
            //			if (ip == 0 || ip == npix/2 || ip == npix-1)
            //				printf("RJH: %d %d %f senz=%f sena=%f solz=%f sola=%f\n", i, ip,
            //						th, data->senz[i * npix + ip],
            //						data->sena[i * npix + ip], data->solz[i * npix + ip],
            //						data->sola[i * npix + ip]);
        }
        //if (i % 100 == 0) printf("Calculated Geometry for scan = %d of %d\n",i,nscans);
    }
    data->have_nav = 1;
    return;
}

aviris4ocia_t *open_aviris(char *filename, char *imgfile, char *navfile, char *gainfile, aviris4ocia_t **data) {

    FILE *ptr;
    char tag[itemSize];
    char *val0;
    char val[itemSize];
    char *inbasename, *infile;
    int i, j, k, status, pos;
    double *indata, *elev;
    float *indataf;
    float rotation;
    int numBands, num;
    char* result;
    char line[itemSize];
    int16_t year, month, day, hour, minute;
    int32_t npix, nscans;
    aviris4ocia_t *temp;
    int isec = 0;
    double sec;
    float valf;

    char projStr[1024];

    cdata_(); //  initialize global FORTRAN common block data for l_sun call

    //    *data = createPrivateData_av(maxBands);
    *data = (aviris4ocia_t*) malloc(sizeof (aviris4ocia_t));
    if (*data == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to allocate private data for aviris\n",
                __FILE__, __LINE__);
        exit(1);
    }

    (*data)->wave = (float *) malloc(maxBands * sizeof (float));
    (*data)->fwhm = (float *) malloc(maxBands * sizeof (float));
    (*data)->gain = (double *) malloc(maxBands * sizeof (double));
    if ((*data)->wave == NULL || (*data)->fwhm == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to allocate scale/offset data for aviris\n",
                __FILE__, __LINE__);
        exit(1);
    }

    temp = *data;
    inbasename = getinbasename_av(filename);
    pos = strlen(inbasename);
    if (pos <= 0) {
        fprintf(stderr, "-E- %s line %d: Not a avalid AVIRIS file %s\n",
                __FILE__, __LINE__, filename);
        exit(-1);
    }

    sscanf(inbasename, "f%2hd%2hd%2hd", &year, &month, &day);

    if (year >= 92) year += 1900;
    else year += 2000;

    sec = 0;
    hour = 0;
    minute = 0;
    isec = (int) sec;

    temp->month = month;
    temp->day = day;
    temp->hour = hour;
    temp->min = minute;
    temp->sec = isec;
    temp->have_nav = 0;
    temp->have_gain = 0;

    ymdhms2ydmsec(year, month, day, hour, minute, isec,
            &temp->year, &temp->doy, &temp->msec);

    sec -= isec;
    temp->msec += sec * 1000;

    printf("Date of AVIRIS flight: Y-%d M-%d D-%d\n", year, month, day);

    if ((ptr = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to open %s\n",
                __FILE__, __LINE__, filename);
        exit(-1);
    }


    int numLinesNeeded = 1;
    int numSamplesNeeded = 1;
    int bandsNeeded = 1;
    int waveLengthNeeded = 1;
    int fwhmNeeded = 1;
    int utmZoneNeeded = 1;
    int eastingNeeded = 1;
    int northingNeeded = 1;
    int pixelSizeNeeded = 1;
    int interleaveNeeded = 1;
    int rotationAngNeeded = 1;

    // loop metadata
    result = fgets(line, itemSize, ptr); // Skip ENVI line and set result
    while ((numLinesNeeded ||
            numSamplesNeeded ||
            bandsNeeded ||
            fwhmNeeded ||
            pixelSizeNeeded ||
            utmZoneNeeded ||
            eastingNeeded ||
            northingNeeded ||
            waveLengthNeeded ||
            rotationAngNeeded ||
            interleaveNeeded) && result) {

        result = fgets(line, itemSize, ptr);
        if (!result) continue;
        //         if(result == NULL) {
        //            fprintf(stderr,"-E- %s line %d: unable to read all of the required metadata from AVIRIS file\n",
        //                    __FILE__,__LINE__);
        //            exit(1);
        //        }
        trimBlanks(line);

        if ((val0 = checkTagLine(line, "lines"))) {
            numLinesNeeded = 0;
            //nscans = (int)atoi(val0);
            sscanf(val0, "%d", &nscans);
            temp->nscans = nscans;
        }
        if ((val0 = checkTagLine(line, "samples"))) {
            numSamplesNeeded = 0;
            sscanf(val0, "%d", &npix);
            //npix = (int)atoi(val0);
            temp->npix = npix;
        }
        if ((val0 = checkTagLine(line, "bands"))) {
            bandsNeeded = 0;
            //numBands = (int)atoi(val0);
            sscanf(val0, "%d", &numBands);
            temp->numBands = numBands;
            if (numBands > maxBands) {
                fprintf(stderr, "-E- %s line %d: number of bands (%d) from AVIRIS file > maxBands (%d)\n",
                        __FILE__, __LINE__, numBands, maxBands);
                exit(1);
            }
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
        }

        if ((val0 = checkTagLine(line, "rotation angle"))) {
            rotationAngNeeded = 0;
            //rotation = atof(val0);
            rotation = checkTagLine_f(line, "rotation angle");
            temp->rotation = rotation;
            if (rotation > 45)
                temp->eastbyscan = -1;
            else //(rotation < -45)
                temp->eastbyscan = 1;
            //             else
            //                 temp->eastbyscan = 0;
        }

        if ((val0 = checkTagLine(line, "pixel size"))) {
            pixelSizeNeeded = 0;
            valf = checkTagLine_f(line, "pixel size");
            temp->pixelSize = valf;
        }
        if ((val0 = checkTagLine(line, "Northing"))) {
            northingNeeded = 0;
            temp->northing = checkTagLine_f(line, "Northing");

        }
        if ((val0 = checkTagLine(line, "Easting"))) {
            eastingNeeded = 0;
            temp->easting = checkTagLine_f(line, "Easting");
            //            temp->easting = atof(val0);
        }
        if ((val0 = checkTagLine(line, "UTM zone"))) {
            utmZoneNeeded = 0;
            //temp->utmZone = (int)atof(val0);
            sscanf(val0, "%d", &temp->utmZone);
        }

        if ((val0 = checkTagLine(line, "wavelength"))) {
            waveLengthNeeded = 0;
            i = 0;
            readWavInfo_jpl(ptr, tag, val);
            while (i < maxBands && strcmp(tag, "}")) {
                temp->wave[i] = atof(val);
                readWavInfo_jpl(ptr, tag, val);
                i++;
            }

            if (!strcmp(tag, "}") && i <= maxBands) {
                temp->wave[i] = atof(val);
                i++;
            } else { // if (i> maxBands) {

                fprintf(stderr, "-E- %s line %d: number of bands (%d) from AVIRIS file > maxBands (%d)\n",
                        __FILE__, __LINE__, numBands, maxBands);
                exit(1);
            }
            numBands = i - 1;
        }

        if ((val0 = checkTagLine(line, "data gain values"))) {
            i = 0;
            readWavInfo_jpl(ptr, tag, val);
            while (i < maxBands && strcmp(tag, "}")) {
                temp->gain[i] = atof(val);
                readWavInfo_jpl(ptr, tag, val);
                i++;
            }

            if (!strcmp(tag, "}") && i <= maxBands) {
                temp->gain[i] = atof(val);
                i++;
            } else { // if (i> AV_MAXBANDS) {

                fprintf(stderr, "-E- %s line %d: number of bands (%d) from AVIRIS file > AV_MAXBANDS (%d)\n",
                        __FILE__, __LINE__, numBands, maxBands);
                exit(1);
            }
            numBands = i - 1;
            temp->have_gain = 1;
        }

        if ((val0 = checkTagLine(line, "fwhm"))) {
            fwhmNeeded = 0;
            i = 0;
            readWavInfo_jpl(ptr, tag, val);
            while (i < maxBands && strcmp(tag, "}")) {
                temp->fwhm[i] = atof(val);
                readWavInfo_jpl(ptr, tag, val);
                i++;
            }
            if (!strcmp(tag, "}") && i < maxBands) {
                temp->fwhm[i] = atof(val);
                i++;
            } else {
                fprintf(stderr, "-E- %s line %d: number of bands (%d) from AVIRIS file > maxBands (%d)\n",
                        __FILE__, __LINE__, numBands, maxBands);
                exit(1);
            }
        }
    }

    fclose(ptr);
    if ((numLinesNeeded ||
            numSamplesNeeded ||
            bandsNeeded ||
            fwhmNeeded ||
            waveLengthNeeded)) {

        fprintf(stderr, "-E- %s line %d: unable to read all of the required metadata from AVIRIS file\n",
                __FILE__, __LINE__);
        exit(1);
    }

    temp->sena = (float *) malloc(nscans * npix * sizeof (float));
    temp->senz = (float *) malloc(nscans * npix * sizeof (float));
    temp->solz = (float *) malloc(nscans * npix * sizeof (float));
    temp->sola = (float *) malloc(nscans * npix * sizeof (float));
    temp->utc = (float *) malloc(nscans * npix * sizeof (float));
    temp->Lt = (float *) malloc(temp->numBands * npix * sizeof (float));
    temp->lon = (double *) malloc(nscans * npix * sizeof (double));
    temp->lat = (double *) malloc(nscans * npix * sizeof (double));
    temp->alt = (float *) malloc(nscans * sizeof (float));

    if (temp->sena == NULL || temp->senz == NULL || temp->sola == NULL || temp->solz == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to allocate sensor and solar angle data for AVIRIS\n",
                __FILE__, __LINE__);
        exit(1);
    }

    if (navfile[0] != '\0') get_nav_data(navfile, nscans, npix, temp);

    if (gainfile[0] != '\0') {
        if ((ptr = fopen(gainfile, "r")) == NULL) {
            fprintf(stderr, "-E- %s line %d: unable to open %s\n",
                    __FILE__, __LINE__, gainfile);
            exit(-1);

        }
        i = 0;
        while ((result = fgets(line, itemSize, ptr)) && i < numBands) {

            if (result == NULL) {
                fprintf(stderr, "-E- %s line %d: unable to read all of the gain data from AVIRIS gain file\n",
                        __FILE__, __LINE__);
                exit(1);
            }
            trimBlanks(line);

            sscanf(line, "%lf", &temp->gain[i]);
            i++;
        }
        temp->have_gain = 1;
    }

    if (navfile[0] == '\0') {
        // Get info about the WGS84 data
        // Get the lat/lon/elev
        numLinesNeeded = 1;
        numSamplesNeeded = 1;

        infile = malloc((pos + strlen("_obs.hdr")) * sizeof (char));
        strcpy(infile, inbasename);
        strcat(infile, "_obs.hdr");

        if ((ptr = fopen(infile, "r")) == NULL) {
            fprintf(stderr, "-E- %s line %d: unable to open %s\n",
                    __FILE__, __LINE__, infile);
            exit(-1);

        }
        // loop metadata
        while (numLinesNeeded ||
                numSamplesNeeded) {

            readNextLine_jpl(ptr, tag, val);
            // skip blank lines
            if (tag[0] == 0)
                continue;

            // get date
            if (!strcmp(tag, "lines")) {
                numLinesNeeded = 0;
                temp->wgs_nscan = (int) atof(val);
            }
            if (!strcmp(tag, "samples")) {
                numSamplesNeeded = 0;
                temp->wgs_npix = (int) atof(val);
            }
        }

        fclose(ptr);


        //Get the lat/lon/elev

        infile = malloc((pos + strlen("_lonlat_eph")) * sizeof (char));
        strcpy(infile, inbasename);
        strcat(infile, "_lonlat_eph");

        printf("Reading lon/lat/elev information from file %s\n", infile);

        // allocate lat and lon storage
        num = 6;
        //    data->lat  = (double *) malloc(data->wgs_nscan*sizeof(double) );
        //    data->lon  = (double *) malloc(data->wgs_nscan*sizeof(double) );
        elev = (double *) malloc(temp->wgs_nscan * sizeof (double));
        indata = (double *) malloc(temp->wgs_nscan * num * sizeof (double));

        //    if(temp->lat==NULL || temp->lon==NULL || temp->elev==NULL || indata==NULL) {
        if (elev == NULL || indata == NULL) {
            fprintf(stderr, "-E- %s line %d: unable to allocate lat/lon data for aviris\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        if ((ptr = fopen(infile, "rb")) == NULL) {
            fprintf(stderr, "-E- %s line %d: unable to open %s\n",
                    __FILE__, __LINE__, infile);
            exit(-1);

        }
        status = fread(indata, sizeof (double), temp->wgs_nscan*num, ptr);
        if (status != temp->wgs_nscan * num) {
            printf("Wrong data read size: want %d got %d in file %s\n", temp->wgs_nscan*num, status, infile);
            exit(1);
        }

        i = 0;

        gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, temp->wgs_nscan);

        double *dist = (double *) calloc(temp->wgs_nscan, sizeof (double));
        double *xlon = (double *) calloc(temp->wgs_nscan, sizeof (double));
        double *xlat = (double *) calloc(temp->wgs_nscan, sizeof (double));

        temp->spl_acc = gsl_interp_accel_alloc();
        temp->lon0 = indata[0];
        temp->lat0 = indata[1];
        temp->distmin = 999;
        temp->distmax = -999;
        while (i < temp->wgs_nscan * num) {
            j = i / num;
            xlon[j] = indata[i];
            xlat[j] = indata[i + 1];
            elev[j] = indata[i + 2];
            //dist[j] = pow((xlon[j]-temp->lon0),2.0) + pow((xlat[j]-temp->lat0),2.0) ;
            dist[j] = (double) angular_distance((double) xlat[j], (double) xlon[j], (double) temp->lat0, (double) temp->lon0);
            //printf("RJH: elev[%d]=%lf dist=%lf lat=%f lon=%f lat0=%lf lon0=%lf\n",j,elev[j],dist[j],xlat[j],xlon[j], temp->lat0, temp->lon0);
            i += num;
        }
        // Sort distances and corresponding elevation values
        gsl_sort2(dist, 1, elev, 1, temp->wgs_nscan);
        temp->distmin = dist[0];
        temp->distmax = dist[temp->wgs_nscan - 2];
        //    for (j=0;j<temp->wgs_nscan;j++)
        //        printf("sort:RJH: elev[%d]=%lf dist=%lf lat=%f lon=%f lat0=%lf lon0=%lf\n",j,elev[j],dist[j],xlat[j],xlon[j], temp->lat0, temp->lon0);

        // Initiate spline
        gsl_spline_init(spline, dist, elev, temp->wgs_nscan);

        temp->spline = spline;

        free(indata);
        fclose(ptr);

        // Get the sensor and solar data

        num = 10;
        indataf = (float *) malloc(npix * sizeof (float)*num);

        infile = malloc((pos + strlen("_obs_ort")) * sizeof (char));
        strcpy(infile, inbasename);
        strcat(infile, "_obs_ort");

        printf("Reading sensor and solar angles information from file %s\n", infile);

        if ((ptr = fopen(infile, "rb")) == NULL) {
            fprintf(stderr, "-E- %s line %d: unable to open %s\n",
                    __FILE__, __LINE__, infile);
            exit(-1);

        }
        for (j = 0; j < nscans; j++) {
            //      for (j=nscans-1;j>=nscans;j--) {
            status = fread(indataf, sizeof (float), num*npix, ptr);
            if (status != num * npix) {
                fprintf(stderr, "-E- %s line %d: AVIRIS Wrong sensor and solar data read size: want %d got %d\n",
                        __FILE__, __LINE__, npix*num, status);
                exit(1);
            }

            for (k = 0; k < npix; k++) {
                temp->sena[j * npix + k] = indataf[1 * npix + k];
                temp->senz[j * npix + k] = indataf[2 * npix + k];
                temp->sola[j * npix + k] = indataf[3 * npix + k];
                temp->solz[j * npix + k] = indataf[4 * npix + k];
                temp->utc [j * npix + k] = indataf[9 * npix + k];
                //            if (temp->utc [j*npix+k] > 0) printf("i=%f %f\n",indataf[8*npix+k],temp->utc [j*npix+k]);
            }
        }

        free(indataf);
        fclose(ptr);


        PJ *pj;
        // init the proj4 projections
        sprintf(projStr, "+proj=utm +zone=%d +ellps=WGS84 +datum=WGS84 +units=m +no_defs ",
                temp->utmZone);
        pj = proj_create_crs_to_crs(PJ_DEFAULT_CTX,
                                   projStr,
                                   "+proj=longlat +ellps=WGS84 +datum=WGS84",
                                   NULL);
        if(pj == NULL) {
            printf("Error - AVIRIS first PROJ projection failed to init\n");
            exit(1);
        }
        temp->pj = proj_normalize_for_visualization(PJ_DEFAULT_CTX, pj);
        if(temp->pj == NULL) {
            printf("Error - AVIRIS visualization PROJ projection failed to init\n");
            exit(1);
        }
        proj_destroy(pj);
        
        // Get the gain data
        free(infile);
    }

    if (!temp->have_gain) {
        infile = malloc((pos + strlen(".gain")) * sizeof (char));
        strcpy(infile, inbasename);
        strcat(infile, ".gain");

        printf("Reading gain information from file %s\n", infile);

        if ((ptr = fopen(infile, "r")) == NULL) {
            fprintf(stderr, "-E- %s line %d: unable to open %s\n  Will assume gains.\n",
                    __FILE__, __LINE__, infile);
            for (k = 0; k < numBands; k++) {
                if (k < 110) temp->gain[i] = 300;
                else if (k < 160) temp->gain[i] = 600;
                else temp->gain[i] = 1200;

            }
            temp->have_gain = 1;

        } else {

            temp->scale_factor = (float *) malloc(numBands * sizeof (float));
            if (temp->gain == NULL) {
                fprintf(stderr, "-E- %s line %d: unable to allocate gain data for AVIRIS\n",
                        __FILE__, __LINE__);
                exit(1);
            }
            for (i = 0; i < numBands && fscanf(ptr, "%lf %d", &temp->gain[i], &k); i++) {
                temp->scale_factor[i] = 1 / (temp->gain[i]);
                printf("gain: %lf %d\n", *(temp->gain + i), k);
            }

            temp->have_gain = 1;
            fclose(ptr);
        }
        free(infile);
    }
    infile = malloc(FILENAME_MAX * sizeof (char));
    if (imgfile[0] == '\0') {
        strcpy(infile, inbasename);
        strcat(infile, "_sc01_ort_img");
    } else {
        strcpy(infile, imgfile);
    }
    printf("Opening AVIRIS image file %s\n", infile);

    if ((temp->av_fp = fopen(infile, "rb")) == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to open %s\n",
                __FILE__, __LINE__, infile);
        exit(-1);

    }

    return temp;

}

int read_aviris(aviris4ocia_t* data, int32_t recnum)
/*
 *  fill standard record with L1B line of data
 */ {
    int status;

    static double last_good_hour = 18;
    static int firstCall = 1;
    static float validSenz, validSena, validSolz, validSola;
    double sec;
    int hour, minute;
    double dist;
    int isec;
    int32_t npix = data->npix, ip;
    static int swap;
    int32_t skip = SKIP;


    if (firstCall) {
        firstCall = 0;


        if (endianess() == 1)
            swap = 1;
        else
            swap = 0;

    }


    if (!data->have_nav) {
        int k = 0;
        //    while ((hour = (int)(data->utc[recnum][k])) <0 && k<npix) k++;
        while ((hour = (int) (data->utc[recnum * npix + k])) < 0 && k < npix) k++;

        if (hour < 0) hour = last_good_hour;

        last_good_hour = hour;

        //    minute = (data->utc[recnum][k] - hour)*60;
        //    sec = ((data->utc[recnum][k] - hour)*60 - minute)*60;
        minute = (data->utc[recnum * npix + k] - hour)*60;
        sec = ((data->utc[recnum * npix + k] - hour)*60 - minute)*60;
        //    printf("Date:utc=%f sec=%f minute=%d (utc-hour)*60=%f ",data->utc[recnum][k],sec,minute,(data->utc[recnum][k] - hour)*60 );
        isec = (int) sec;

        data->hour = hour;
        data->min = minute;
        data->sec = sec;

        ymdhms2ydmsec(data->year, data->month, data->day, hour, minute, isec,
                &data->year, &data->doy, &data->msec);

        //data->msec = sec * 1000;

        //data->scantime = yds2unix(data->year, data->doy, (double) (data->msec / 1.e3));
        data->scantime = ymds2unix(data->year, data->month, data->day, (hour * 3600 + minute * 60 + sec));
        //    printf("Date=%4d %d %d\n",*(l1rec->year),*(l1rec->day),*(l1rec->msec));

        for (ip = 0; ip < npix; ip++) {

            // Rotate the image
            data->lon[recnum * npix + ip] = data->easting + ip * cos(deg2rad(data->rotation)) * data->pixelSize - recnum * sin(deg2rad(data->rotation)) * data->pixelSize; // starts in upper left corner
            data->lat[recnum * npix + ip] = data->northing - ip * sin(deg2rad(data->rotation)) * data->pixelSize - recnum * cos(deg2rad(data->rotation)) * data->pixelSize;

            if (isnan(data->lat[recnum * npix + ip])) data->lat[recnum * npix + ip] = -999.0;
            if (isnan(data->lon[recnum * npix + ip])) data->lon[recnum * npix + ip] = -999.0;

            if (data->senz[recnum * npix + ip] <= skip) data->senz[recnum * npix + ip] = getValidAngleB(&data->senz[recnum * npix], npix, skip, &validSenz);
            if (data->sena[recnum * npix + ip] <= skip) data->sena[recnum * npix + ip] = getValidAngleB(&data->sena[recnum * npix], npix, skip, &validSena);
            if (data->solz[recnum * npix + ip] <= skip) data->solz[recnum * npix + ip] = getValidAngleB(&data->solz[recnum * npix], npix, skip, &validSolz);
            if (data->sola[recnum * npix + ip] <= skip) data->sola[recnum * npix + ip] = getValidAngleB(&data->sola[recnum * npix], npix, skip, &validSola);
        }

        aviris4ocia_proj4_convert(data, npix, &data->lon[recnum * npix], &data->lat[recnum * npix]);


        // find interpolated elevation from wgs-84 lat/lon
        ip = npix / 2;
        //        dist = (data->lon[recnum*npix+ip]-data->lon0)*(data->lon[recnum*npix+ip]-data->lon0) +
        //                (data->lat[recnum*npix+ip]-data->lat0)*(data->lat[recnum*npix+ip]-data->lat0);
        dist = (float) angular_distance((double) data->lat[recnum * npix + ip], (double) data->lon[recnum * npix + ip], (double) data->lat0, (double) data->lon0);
        if (dist > data->distmax) {
            printf("lat/lon > range of wgs coordinates - using altitude of nearest neighbor\n");
            dist = data->distmax;
        } else if (dist < data->distmin) {
            printf("lat/lon < range of wgs coordinates - using altitude of nearest neighbor\n");
            dist = data->distmin;
        }
        data->alt[recnum] = (float) gsl_spline_eval(data->spline, dist, data->spl_acc);

    } else {
        hour = (int) (data->utc[recnum]);
        minute = (data->utc[recnum] - hour) * 60;
        sec = ((data->utc[recnum] - hour) * 60 - minute) * 60;
        data->scantime = ymds2unix(data->year, data->month, data->day,
                (hour * 3600. + minute * 60. + sec));
        data->hour = hour;
        data->min = minute;
        data->sec = sec;

    }

    status = readBinScanLine4ocia_int2(data->Lt, recnum, data->npix, data->gain, data->numBands, data->numBands, data->interleave, swap, data->av_fp);

    if (status == data->npix * data->numBands)
        return (0);
    else
        return (1);

}

int close_aviris(aviris4ocia_t *data) {

    // undo what open allocated

    freePrivateData_ocia(data);
    fclose(data->av_fp);
    return 0;
}

int checkAvProcessFile(char *filename, char *hdrfile, char *imgfile, char *navfile, char *gainfile, int itemsize) {
    FILE *fh;
    char *result;
    int status = 0;
    char line[itemsize];
    char desc[itemsize];
    char fobj[itemsize];

    if ((fh = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to open %s\n",
                __FILE__, __LINE__, filename);
        exit(-1);
    }
    imgfile[0] = '\0';
    navfile[0] = '\0';
    gainfile[0] = '\0';
    if (fgets(line, itemsize, fh)) {
        sscanf(line, "%s", desc);
        if (strstr(desc, "AVIRIS_PREPROC_HDR")) {
            while ((result = fgets(line, itemsize, fh))) {

                if (result == NULL) {
                    fprintf(stderr, "-E- %s line %d: unable to read all of the gain data from AVIRIS gain file\n",
                            __FILE__, __LINE__);
                    exit(1);
                }
                trimBlanks(line);

                sscanf(line, "%s = %s", desc, fobj);
                if (strstr(desc, "hdrfile"))
                    strncpy(hdrfile, fobj, itemsize);
                if (strstr(desc, "navfile"))
                    strncpy(navfile, fobj, itemsize);
                if (strstr(desc, "gainfile"))
                    strncpy(gainfile, fobj, itemsize);
                if (strstr(desc, "imgfile"))
                    strncpy(imgfile, fobj, itemsize);
            }
            status = 1;
        }
    }
    fclose(fh);
    return status;
}



