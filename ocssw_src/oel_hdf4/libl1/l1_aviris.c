/*
 * l1_aviris.c
 *
 *  Created on: May 18, 2015
 *      Author: Rick Healy SAIC
 *              NASA-GSFC OBPG
 */
#include <stdlib.h>
#include <stdio.h>

#include "l1.h"
#include "jplaeriallib.h"
#include "l1_aviris.h"
#include <math.h>

#include <libnav.h>
#include <proj.h>

#define SKIP -9999
static const int maxBands = 224;

/*
 * RJH - 11/14/2016 - Added code for processing preprocessed hdr data for old Aviris files
 * RJH - 10/24/2016 - If gain file is missing, assume values are the same as every scene seems to be
 *					  Removed unused commented code
 * Todo: 7/6/2016 - rjh
 *       1) remove repeating and overlapping bands (wait for Bryan's rayleigh files)
 *       2) Add Atrem files for gas correction
 */
aviris_t* createPrivateData_aviris(int numBands) {

    aviris_t* data = (aviris_t*) calloc(1, sizeof (aviris_t));
    if (data == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to allocate private data for aviris\n",
                __FILE__, __LINE__);
        exit(1);
    }

    data->wave = (double *) malloc(numBands * sizeof (double));
    data->fwhm = (double *) malloc(numBands * sizeof (double));
    if (data->wave == NULL || data->fwhm == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to allocate scale/offset data for aviris\n",
                __FILE__, __LINE__);
        exit(1);
    }

    data->have_nav = 0;
    data->have_gain = 0;

    return data;
}

void freePrivateData_aviris(aviris_t* data) {
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

void get_aviris_nav_data(char* navfile, int32_t nscans, int32_t npix, aviris_t* data) {

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
    static double *range, *lon, *lat;
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

    if (!data->lon)
        data->lon = (double*) malloc(npix * nscans * sizeof (double));
    if (!data->lat)
        data->lat = (double*) malloc(npix * nscans * sizeof (double));
    if (!range) range = (double*) malloc(nscans * sizeof (double));
    if (!lon) lon = (double*) malloc(nscans * sizeof (double));
    if (!lat) lat = (double*) malloc(nscans * sizeof (double));
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
        sscanf(line, "%lf %f %lf %lf ", &data->utc[i], &data->alt[i], &lat[i],
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
        //getPosVecR(latitude,longitude, data->alt[i], pos); // get position vector of sensor
        unix2yds(ymds2unix(data->year, data->month, data->day,
                (hour * 3600. + minute * 60. + sec)), &year, &doy, &secondOfDay);
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

int openl1_aviris(filehandle *file) {

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
    int year, month, day, hour, minute;

    int isec = 0;
    double sec;

    char projStr[1024];

    aviris_t* data = file->private_data;

    cdata_(); //  initialize global FORTRAN common block data for l_sun call

    data->isnetcdf = 0;
    inbasename = getinbasename_av(data->hdrfile);
    pos = strlen(inbasename);
    if (pos <= 0) {
        fprintf(stderr, "-E- %s line %d: Not a avalid AVIRIS file %s\n",
                __FILE__, __LINE__, data->hdrfile);
        exit(-1);
    }

    data->wave = (double *) malloc(maxBands * sizeof (double));
    data->fwhm = (double *) malloc(maxBands * sizeof (double));
    data->gain = (double *) malloc(maxBands * sizeof (double));

    sscanf(inbasename, "f%2d%2d%2d", &year, &month, &day);

    if (year >= 92) year = year + 1900;
    else year = year + 2000;

    sec = 0;
    hour = 0;
    minute = 0;
    isec = (int) sec;

    data->month = month;
    data->day = day;

    ymdhms2ydmsec(year, month, day, hour, minute, isec,
            &data->year, &data->doy, &data->msec);

    sec -= isec;
    data->msec += sec * 1000;

    printf("Date of AVIRIS flight: Y-%d M-%d D-%d\n", year, month, day);

    if ((ptr = fopen(data->hdrfile, "r")) == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to open %s\n",
                __FILE__, __LINE__, data->hdrfile);
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
            file->nscan = atoi(val0);
            data->nscan = file->nscan;
        }
        if ((val0 = checkTagLine(line, "samples"))) {
            numSamplesNeeded = 0;
            file->npix = atoi(val0);
            data->npix = file->npix;
        }
        if ((val0 = checkTagLine(line, "bands"))) {
            bandsNeeded = 0;
            numBands = atoi(val0);
            data->numBands = numBands;
            if (numBands > AV_MAXBANDS) {
                fprintf(stderr, "-E- %s line %d: number of bands (%d) from AVIRIS file > maxBands (%d)\n",
                        __FILE__, __LINE__, numBands, AV_MAXBANDS);
                exit(1);
            }
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
        }

        if ((val0 = checkTagLine(line, "rotation angle"))) {
            rotationAngNeeded = 0;
            rotation = atof(val0);
            data->rotation = rotation;
            if (rotation > 45)
                data->eastbyscan = -1;
            else if (rotation < -45)
                data->eastbyscan = 1;
            else
                data->eastbyscan = 0;
        }

        if ((val0 = checkTagLine(line, "pixel size"))) {
            pixelSizeNeeded = 0;
            data->pixelSize = atof(val0);
        }
        if ((val0 = checkTagLine(line, "Northing"))) {
            northingNeeded = 0;
            data->northing = atof(val0);
        }
        if ((val0 = checkTagLine(line, "Easting"))) {
            eastingNeeded = 0;
            data->easting = atof(val0);
        }
        if ((val0 = checkTagLine(line, "UTM zone"))) {
            utmZoneNeeded = 0;
            data->utmZone = atoi(val0);
        }

        if ((val0 = checkTagLine(line, "wavelength"))) {
            waveLengthNeeded = 0;
            i = 0;
            readWavInfo_jpl(ptr, tag, val);
            while (i < AV_MAXBANDS && strcmp(tag, "}")) {
                data->wave[i] = atof(val);
                readWavInfo_jpl(ptr, tag, val);
                i++;
            }

            if (!strcmp(tag, "}") && i <= AV_MAXBANDS) {
                data->wave[i] = atof(val);
                i++;
            } else { // if (i> AV_MAXBANDS) {

                fprintf(stderr, "-E- %s line %d: number of bands (%d) from AVIRIS file > AV_MAXBANDS (%d)\n",
                        __FILE__, __LINE__, file->nbands, AV_MAXBANDS);
                exit(1);
            }
            numBands = i - 1;
        }
        if ((val0 = checkTagLine(line, "data gain values"))) {
            data->gain = (double *) malloc(AV_MAXBANDS * sizeof (double));
            i = 0;
            readWavInfo_jpl(ptr, tag, val);
            while (i < AV_MAXBANDS && strcmp(tag, "}")) {
                data->gain[i] = atof(val);
                readWavInfo_jpl(ptr, tag, val);
                i++;
            }

            if (!strcmp(tag, "}") && i <= AV_MAXBANDS) {
                data->gain[i] = atof(val);
                i++;
            } else { // if (i> AV_MAXBANDS) {

                fprintf(stderr, "-E- %s line %d: number of bands (%d) from AVIRIS file > AV_MAXBANDS (%d)\n",
                        __FILE__, __LINE__, file->nbands, AV_MAXBANDS);
                exit(1);
            }
            numBands = i - 1;
            data->have_gain = 1;
        }


        if ((val0 = checkTagLine(line, "fwhm"))) {
            fwhmNeeded = 0;
            i = 0;
            readWavInfo_jpl(ptr, tag, val);
            while (i < AV_MAXBANDS && strcmp(tag, "}")) {
                data->fwhm[i] = atof(val);
                readWavInfo_jpl(ptr, tag, val);
                i++;
            }
            if (!strcmp(tag, "}") && i < AV_MAXBANDS) {
                data->fwhm[i] = atof(val);
                i++;
            } else {
                fprintf(stderr, "-E- %s line %d: number of bands (%d) from AVIRIS file > AV_MAXBANDS (%d)\n",
                        __FILE__, __LINE__, file->nbands, AV_MAXBANDS);
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
    int nscans = file->nscan;
    int npix = file->npix;

    data->sena = (double *) malloc(nscans * npix * sizeof (double));
    data->senz = (double *) malloc(nscans * npix * sizeof (double));
    data->solz = (double *) malloc(nscans * npix * sizeof (double));
    data->sola = (double *) malloc(nscans * npix * sizeof (double));
    data->utc = (double *) malloc(nscans * npix * sizeof (double));
    data->alt = (float *) malloc(nscans * sizeof (float));

    if (data->sena == NULL || data->senz == NULL || data->sola == NULL || data->solz == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to allocate sensor and solar angle data for AVIRIS\n",
                __FILE__, __LINE__);
        exit(1);
    }

    if (data->navfile[0] != '\0') get_aviris_nav_data(data->navfile, nscans, npix, data);

    if (data->gainfile[0] != '\0') {
        if ((ptr = fopen(data->gainfile, "r")) == NULL) {
            fprintf(stderr, "-E- %s line %d: unable to open %s\n",
                    __FILE__, __LINE__, data->gainfile);
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

            sscanf(line, "%lf", &data->gain[i]);
            i++;
        }
        data->have_gain = 1;
    }

    if (data->navfile[0] == '\0') {
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
                data->wgs_nscan = atoi(val);
            }
            if (!strcmp(tag, "samples")) {
                numSamplesNeeded = 0;
                data->wgs_npix = atoi(val);
            }
        }

        fclose(ptr);

        //Get the lat/lon/elev

        infile = malloc((pos + strlen("_lonlat_eph")) * sizeof (char));
        strcpy(infile, inbasename);
        strcat(infile, "_lonlat_eph");

        printf("Reading lon/lat/elev information from file %s\n", infile);

        num = 6;
        elev = (double *) malloc(data->wgs_nscan * sizeof (double));
        indata = (double *) malloc(data->wgs_nscan * num * sizeof (double));

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
        status = fread(indata, sizeof (double), data->wgs_nscan*num, ptr);
        if (status != data->wgs_nscan * num) {
            printf("Wrong data read size: want %d got %d in file %s\n", data->wgs_nscan*num, status, infile);
            exit(1);
        }

        i = 0;

        gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, data->wgs_nscan);

        double *dist = (double *) calloc(data->wgs_nscan, sizeof (double));
        double *xlon = (double *) calloc(data->wgs_nscan, sizeof (double));
        double *xlat = (double *) calloc(data->wgs_nscan, sizeof (double));

        data->spl_acc = gsl_interp_accel_alloc();
        data->lon0 = indata[0];
        data->lat0 = indata[1];
        data->distmin = 999;
        data->distmax = -999;
        while (i < data->wgs_nscan * num) {
            j = i / num;
            xlon[j] = indata[i];
            xlat[j] = indata[i + 1];
            elev[j] = indata[i + 2];
            dist[j] = (double) angular_distance((double) xlat[j], (double) xlon[j], (double) data->lat0, (double) data->lon0);
            i += num;
        }
        i = num;
        // Sort distances and corresponding elevation values
        gsl_sort2(dist, 1, elev, 1, data->wgs_nscan);
        data->distmin = dist[0];
        data->distmax = dist[data->wgs_nscan - 2];

        // Initiate spline
        gsl_spline_init(spline, dist, elev, data->wgs_nscan);

        data->spline = spline;

        free(indata);
        fclose(ptr);


        // Get the sensor and solar data

        num = 10;
        indataf = (float *) malloc(file->npix * sizeof (float)*num);

        if (data->sena == NULL || data->senz == NULL || data->sola == NULL || data->solz == NULL || indataf == NULL) {
            fprintf(stderr, "-E- %s line %d: unable to allocate sensor and solar angle data for AVIRIS\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        infile = malloc((pos + strlen("_obs_ort")) * sizeof (char));
        strcpy(infile, inbasename);
        strcat(infile, "_obs_ort");

        printf("Reading sensor and solar angles information from file %s\n", infile);

        if ((ptr = fopen(infile, "rb")) == NULL) {
            fprintf(stderr, "-E- %s line %d: unable to open %s\n",
                    __FILE__, __LINE__, infile);
            exit(-1);

        }
        for (j = 0; j < file->nscan; j++) {
            status = fread(indataf, sizeof (float), num * file->npix, ptr);
            if (status != num * file->npix) {
                fprintf(stderr, "-E- %s line %d: AVIRIS Wrong sensor and solar data read size: want %d got %d\n",
                        __FILE__, __LINE__, file->npix*num, status);
                exit(1);
            }

            for (k = 0; k < file->npix; k++) {
                data->sena[j * file->npix + k] = indataf[1 * file->npix + k];
                data->senz[j * file->npix + k] = indataf[2 * file->npix + k];
                data->sola[j * file->npix + k] = indataf[3 * file->npix + k];
                data->solz[j * file->npix + k] = indataf[4 * file->npix + k];
                data->utc [j * file->npix + k] = indataf[9 * file->npix + k];
            }
        }

        free(indataf);
        fclose(ptr);

        PJ *pj;
        sprintf(projStr, "+proj=utm +zone=%d +ellps=WGS84 +datum=WGS84 +units=m +no_defs ",
                data->utmZone);

        // init the proj4 projections
        pj = proj_create_crs_to_crs(PJ_DEFAULT_CTX,
                                   projStr,
                                   "+proj=longlat +ellps=WGS84 +datum=WGS84",
                                   NULL);
        if(pj == NULL) {
            printf("Error - AVIRIS first PROJ projection failed to init\n");
            exit(1);
        }
        data->pj = proj_normalize_for_visualization(PJ_DEFAULT_CTX, pj);
        if(data->pj == NULL) {
            printf("Error - AVIRIS visualization PROJ projection failed to init\n");
            exit(1);
        }
        proj_destroy(pj);

        free(infile);
    } // no input navfile


    // Get the gain data
    if (!data->have_gain) {
        infile = malloc((pos + strlen(".gain")) * sizeof (char));
        strcpy(infile, inbasename);
        strcat(infile, ".gain");

        printf("Attempting to read gain information from file %s\n", infile);

        if ((ptr = fopen(infile, "r")) == NULL) {
            fprintf(stderr, "-E- %s line %d: unable to open %s\n  Will assume gains.\n",
                    __FILE__, __LINE__, infile);
            for (k = 0; k < numBands; k++) {
                if (k < 110) data->gain[i] = 300;
                else if (k < 160) data->gain[i] = 600;
                else data->gain[i] = 1200;

            }
        } else {

            if (data->gain == NULL) {
                fprintf(stderr, "-E- %s line %d: unable to allocate gain data for AVIRIS\n",
                        __FILE__, __LINE__);
                exit(1);
            }
            for (i = 0; i < numBands && fscanf(ptr, "%lf %d", (data->gain + i), &k); i++) {
                //        printf("gain: %lf %d\n",*(data->gain+i),k);
            }
            fclose(ptr);
        }

        data->have_gain = 1;
        free(infile);
    }
    infile = malloc(FILENAME_MAX * sizeof (char));
    if (data->imgfile[0] == '\0') {
        strcpy(infile, inbasename);
        strcat(infile, "_sc01_ort_img");
    } else {
        strcpy(infile, data->imgfile);
    }
    printf("Opening AVIRIS image file %s\n", infile);

    if ((data->av_fp = fopen(infile, "rb")) == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to open %s\n",
                __FILE__, __LINE__, infile);
        exit(-1);

    }

    return (0);

}

int readl1_aviris(filehandle *file, int32_t recnum, l1str *l1rec)
/*
 *  fill standard record with L1B line of data
 */ {

    static double last_good_hour = 18;
    static int firstCall = 1;
    double sec;
    int hour, minute;
    double dist;
    int npix = file->npix, ip, ib, ipb;
    static int swap, *ibndx, *ipbndx;
    static float *Lt;


    aviris_t* data = (aviris_t*) file->private_data;
    int ipb_av, ibl1;

    if (firstCall) {
        if (want_verbose)
            printf("l1file->nbands = %d\n", (int) file->nbands);
        firstCall = 0;

        for (ip = 0; ip < npix; ip++) {
            l1rec->pixnum[ip] = ip;
            l1rec->flags[ip] = 0;
        }

        if (endianess() == 1)
            swap = 1;
        else
            swap = 0;

        Lt = (float*) malloc(npix * data->numBands * sizeof (float));
        if (!Lt) {
            fprintf(stderr, "-E- %s line %d: Out of Memory allocating buffer array in readl1_aviris.\n",
                    __FILE__, __LINE__);
            exit(1);

        }
        ipbndx = (int*) malloc(npix * file->nbands * sizeof (int));
        if (!ipbndx) {
            fprintf(stderr, "-E- %s line %d: Out of Memory allocating buffer array in readl1_aviris.\n",
                    __FILE__, __LINE__);
            exit(1);

        }
        ibndx = (int*) malloc(file->nbands * sizeof (int));
        if (!ibndx) {
            fprintf(stderr, "-E- %s line %d: Out of Memory allocating buffer array in readl1_aviris.\n",
                    __FILE__, __LINE__);
            exit(1);

        }

        ibl1 = 0;
        for (ib = 0; ib < 31; ib++) {
            ibndx[ibl1] = ib;
            ibl1++;
        }
        for (ib = 33; ib < 95; ib++) {
            ibndx[ibl1] = ib;
            ibl1++;
        }
        for (ib = 97; ib < 159; ib++) {
            ibndx[ibl1] = ib;
            ibl1++;
        }
        for (ib = 161; ib < 206; ib++) {
            ibndx[ibl1] = ib;
            ibl1++;
        }
        file->fwhm = (float*) malloc(data->numBands * sizeof (float));
        if (!file->fwhm) {
            fprintf(stderr, "-E- %s line %d: Out of Memory allocating fwhm array in readl1_aviris.\n",
                    __FILE__, __LINE__);
            exit(1);

        }
        for (ib = 0; ib < ibl1; ib++) file->fwhm[ib] = data->fwhm[ibndx[ib]] / 1000.; // convert from nm to microns

        // create an index array to map the Lt's to the wavelength indexes actually used
        for (ip = 0; ip < npix; ip++) {
            ibl1 = 0;
            for (ib = 0; ib < 31; ib++) {
                ipb = ip * file->nbands + ibl1;
                switch (data->interleave) {
                case BIP:
                    ipb_av = ip * data->numBands + ib;
                    break;
                case BIL:
                    ipb_av = ib * npix + ip;
                    break;
                }
                ipbndx[ipb] = ipb_av;
                ibl1++;
            }
            for (ib = 33; ib < 95; ib++) {
                ipb = ip * file->nbands + ibl1;
                switch (data->interleave) {
                case BIP:
                    ipb_av = ip * data->numBands + ib;
                    break;
                case BIL:
                    ipb_av = ib * npix + ip;
                    break;
                }
                ipbndx[ipb] = ipb_av;
                ibl1++;
            }
            for (ib = 97; ib < 159; ib++) {
                ipb = ip * file->nbands + ibl1;
                switch (data->interleave) {
                case BIP:
                    ipb_av = ip * data->numBands + ib;
                    break;
                case BIL:
                    ipb_av = ib * npix + ip;
                    break;
                }
                ipbndx[ipb] = ipb_av;
                ibl1++;
            }
            for (ib = 161; ib < 206; ib++) {
                ipb = ip * file->nbands + ibl1;
                switch (data->interleave) {
                case BIP:
                    ipb_av = ip * data->numBands + ib;
                    break;
                case BIL:
                    ipb_av = ib * npix + ip;
                    break;
                }
                ipbndx[ipb] = ipb_av;
                ibl1++;
            }
        }

    }

    //  set information about data
    l1rec->npix = file->npix;
    l1rec->l1file->sensorID = file->sensorID;

    if (!data->have_nav) {

        int k = 0;
        while ((hour = (int) (data->utc[recnum * npix + k])) < 0 && k < npix) k++;

        if (hour < 0) hour = last_good_hour;

        last_good_hour = hour;

        minute = (data->utc[recnum * npix + k] - hour)*60;
        sec = ((data->utc[recnum * npix + k] - hour)*60 - minute)*60;

        l1rec->scantime = ymds2unix(data->year, data->month, data->day, (hour * 3600 + minute * 60 + sec));
        
        PJ_COORD c, c_out;

        // set default z and t
        c.xyzt.z = 0.0;
        c.xyzt.t = HUGE_VAL;

        for (ip = 0; ip < npix; ip++) {

            l1rec->pixnum[ip] = ip;

            c.xy.x = data->easting + ip * cos(deg2rad(data->rotation)) * data->pixelSize - recnum * sin(deg2rad(data->rotation)) * data->pixelSize; // starts in upper left corner
            c.xy.y = data->northing - ip * sin(deg2rad(data->rotation)) * data->pixelSize - recnum * cos(deg2rad(data->rotation)) * data->pixelSize;
            c_out = proj_trans(data->pj, PJ_FWD, c);
            
            if(isfinite(c_out.xy.x) && isfinite(c_out.xy.y)) {
                l1rec->lon[ip] = c_out.xy.x;
                l1rec->lat[ip] = c_out.xy.y;
            } else {
                l1rec->lon[ip] = -999.0;
                l1rec->lat[ip] = -999.0;
            }
            
            if (l1rec->lon[ip] < -181.0 || l1rec->lon[ip] > 181.0 ||
                    l1rec->lat[ip] < -91.0 || l1rec->lat[ip] > 91.0)
                l1rec->navfail[ip] = 1;

            if (data->senz[recnum * npix + ip] > SKIP)
                l1rec->senz[ip] = data->senz[recnum * npix + ip];
            else
                l1rec->senz[ip] = getValidAngle(&data->senz[recnum * npix], npix, SKIP);

            if (data->sena[recnum * npix + ip] > SKIP)
                l1rec->sena[ip] = data->sena[recnum * npix + ip];
            else
                l1rec->sena[ip] = getValidAngle(&data->sena[recnum * npix], npix, SKIP);

            if (data->solz[recnum * npix + ip] > SKIP)
                l1rec->solz[ip] = data->solz[recnum * npix + ip];
            else
                l1rec->solz[ip] = getValidAngle(&data->solz[recnum * npix], npix, SKIP);

            if (data->sola[recnum * npix + ip] > SKIP)
                l1rec->sola[ip] = data->sola[recnum * npix + ip];
            else
                l1rec->sola[ip] = getValidAngle(&data->sola[recnum * npix], npix, SKIP);


        }
        // find interpolated elevation from wgs-84 lat/lon
        ip = npix / 2;
        dist = (double) angular_distance((double) l1rec->lat[ip], (double) l1rec->lon[ip], (double) data->lat0, (double) data->lon0);
        if (dist > data->distmax) {
            printf("lat/lon > range of wgs coordinates - using altitude of nearest neighbor\n");
            dist = data->distmax;
        } else if (dist < data->distmin) {
            printf("lat/lon < range of wgs coordinates - using altitude of nearest neighbor\n");
            dist = data->distmin;
        }
        l1rec->alt = (float) gsl_spline_eval(data->spline, dist, data->spl_acc) / 1000.; // and convert to km

    } else { // have_nav
        l1rec->alt = data->alt[recnum];
        for (ip = 0; ip < npix; ip++) {
            l1rec->senz[ip] = data->senz[recnum * npix + ip];
            l1rec->sena[ip] = data->sena[recnum * npix + ip];
            l1rec->solz[ip] = data->solz[recnum * npix + ip];
            l1rec->sola[ip] = data->sola[recnum * npix + ip];
            l1rec->lon[ip] = data->lon[recnum * npix + ip];
            l1rec->lat[ip] = data->lat[recnum * npix + ip];
            if (l1rec->lon[ip] < -181.0 || l1rec->lon[ip] > 181.0 ||
                    l1rec->lat[ip] < -91.0 || l1rec->lat[ip] > 91.0)
                l1rec->navfail[ip] = 1;
        }
        hour = (int) (data->utc[recnum]);
        minute = (data->utc[recnum] - hour) * 60;
        sec = ((data->utc[recnum] - hour) * 60 - minute) * 60;
        l1rec->scantime = ymds2unix(data->year, data->month, data->day,
                (hour * 3600. + minute * 60. + sec));
    }
    readBinScanLine_int2(Lt, recnum, l1rec->npix, data->gain, data->numBands, data->interleave, swap, data->av_fp);

    for (ipb = 0; ipb < npix * file->nbands; ipb++)
        l1rec->Lt[ipb] = Lt[ipbndx[ipb]];

    return (LIFE_IS_GOOD);
}

int closel1_aviris(filehandle *file) {
    aviris_t* data = (aviris_t*) file->private_data;

    // undo what open allocated

    freePrivateData_aviris(data);
    free(file->private_data);

    return 0;
}


