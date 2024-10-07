
/* ============================================================================ */
/* module l1_meris_CC.c - functions to read MERIS CC (coastal color) for MSL12  */
/* Written By:  Joel M. Gales GSFC Futuretech, Sep. 2011.                       */
/*                                                                              */
/*  W. Robinson, SAIC, 22 May 2012 make repair to start, end and scan times     */
/* ============================================================================ */

#include "l1_meris_CC.h"

#include "l1.h"

#include <stdbool.h>
#include <stdlib.h>
#include <netcdf.h>
#include <math.h>

#define MERIS_NBANDS 15
#define MAXOCLIN 6700      /* max # lines */

static float *lon, *lat, *senz, *sena, *solz, *sola;
static int ncol;
static int32_t nrow;
static int16_t nline, npix;

static int32_t spix = 0;
static int year, day, msec;
static double time_interval;



/* ------------------------------------------------------------ */
/* navigation_meris() - generates navigation info interpolated  */
/*              for each pixel, using info from a MERIS_CC      */
/*              NETCDF file                                     */
/*                                                              */
/* ------------------------------------------------------------ */
int navigation_meris(int32_t fileID) {
    float *xctl;
    float *inlon, *inlat;
    float *insolz, *insola;
    float *insenz, *insena;
    float *yctl;
    float *in1, *in2;

    int i, j, k, l;
    int nlon, nlat;
    float out1[MAXOCLIN], out2[MAXOCLIN];
    float *spl_aux;
    float *x_ctl_ll, *y_ctl_ll, *z_ctl_ll;
    float *x_ctl_sol, *y_ctl_sol, *z_ctl_sol;
    float *x_ctl_sen, *y_ctl_sen, *z_ctl_sen;
    float * in_ptr[3][3], *out_ptr[3][2], *tout_ptr[3][2];

    float off_x, off_y, sub_x, sub_y;
    float *tlon, *tlat, *tsenz, *tsena, *tsolz, *tsola;

    int geo_id;

    inlon = (float *) calloc(ncol*nrow, sizeof (float));
    inlat = (float *) calloc(ncol*nrow, sizeof (float));
    insolz = (float *) calloc(ncol*nrow, sizeof (float));
    insola = (float *) calloc(ncol*nrow, sizeof (float));
    insenz = (float *) calloc(ncol*nrow, sizeof (float));
    insena = (float *) calloc(ncol*nrow, sizeof (float));

    xctl = (float *) calloc(ncol, sizeof (float));
    yctl = (float *) calloc(nrow, sizeof (float));
    in1 = (float *) calloc(nrow, sizeof (float));
    in2 = (float *) calloc(nrow, sizeof (float));

    /* read the data sets */
    nc_inq_varid(fileID, "longitude", &geo_id);
    nc_get_var_float(fileID, geo_id, inlon);
    nc_get_att_float(fileID, geo_id, "offset_x", &off_x);
    nc_get_att_float(fileID, geo_id, "offset_y", &off_y);
    nc_get_att_float(fileID, geo_id, "subsampling_x", &sub_x);
    nc_get_att_float(fileID, geo_id, "subsampling_y", &sub_y);

    nc_inq_varid(fileID, "latitude", &geo_id);
    nc_get_var_float(fileID, geo_id, inlat);
    nc_inq_varid(fileID, "sun_zenith", &geo_id);
    nc_get_var_float(fileID, geo_id, insolz);
    nc_inq_varid(fileID, "sun_azimuth", &geo_id);
    nc_get_var_float(fileID, geo_id, insola);
    nc_inq_varid(fileID, "view_zenith", &geo_id);
    nc_get_var_float(fileID, geo_id, insenz);
    nc_inq_varid(fileID, "view_azimuth", &geo_id);
    nc_get_var_float(fileID, geo_id, insena);


    /* define control point grid */
    if (want_verbose)
        printf("nrow,ncol = %d %d \n", nrow, ncol);
    for (i = 0; i < ncol; i++) {
        xctl[i] = ((float) i) * sub_x + off_x;
    }

    /* define output grid */
    nlon = npix; /* # of output pixs */
    nlat = nline;

    for (i = 0; i < nrow; i++) {
        yctl[i] = ((float) i) * sub_y + off_y;
    }


    /* Compute unit vectors from lon/lat of control points */
    x_ctl_ll = (float *) calloc(ncol*nrow, sizeof (float));
    y_ctl_ll = (float *) calloc(ncol*nrow, sizeof (float));
    z_ctl_ll = (float *) calloc(ncol*nrow, sizeof (float));

    x_ctl_sol = (float *) calloc(ncol*nrow, sizeof (float));
    y_ctl_sol = (float *) calloc(ncol*nrow, sizeof (float));
    z_ctl_sol = (float *) calloc(ncol*nrow, sizeof (float));

    x_ctl_sen = (float *) calloc(ncol*nrow, sizeof (float));
    y_ctl_sen = (float *) calloc(ncol*nrow, sizeof (float));
    z_ctl_sen = (float *) calloc(ncol*nrow, sizeof (float));


    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) {
            inlon[i * ncol + j] = inlon[i * ncol + j] / RADEG;
            inlat[i * ncol + j] = inlat[i * ncol + j] / RADEG;

            x_ctl_ll[i * ncol + j] = cos(inlat[i * ncol + j]) * cos(inlon[i * ncol + j]);
            y_ctl_ll[i * ncol + j] = cos(inlat[i * ncol + j]) * sin(inlon[i * ncol + j]);
            z_ctl_ll[i * ncol + j] = sin(inlat[i * ncol + j]);


            insola[i * ncol + j] = insola[i * ncol + j] / RADEG;
            insolz[i * ncol + j] = insolz[i * ncol + j] / RADEG;

            x_ctl_sol[i * ncol + j] = cos(insolz[i * ncol + j]) * cos(insola[i * ncol + j]);
            y_ctl_sol[i * ncol + j] = cos(insolz[i * ncol + j]) * sin(insola[i * ncol + j]);
            z_ctl_sol[i * ncol + j] = sin(insolz[i * ncol + j]);


            insena[i * ncol + j] = insena[i * ncol + j] / RADEG;
            insenz[i * ncol + j] = insenz[i * ncol + j] / RADEG;

            x_ctl_sen[i * ncol + j] = cos(insenz[i * ncol + j]) * cos(insena[i * ncol + j]);
            y_ctl_sen[i * ncol + j] = cos(insenz[i * ncol + j]) * sin(insena[i * ncol + j]);
            z_ctl_sen[i * ncol + j] = sin(insenz[i * ncol + j]);
        }
    }

    in_ptr[0][0] = x_ctl_ll;
    in_ptr[0][1] = y_ctl_ll;
    in_ptr[0][2] = z_ctl_ll;
    in_ptr[1][0] = x_ctl_sol;
    in_ptr[1][1] = y_ctl_sol;
    in_ptr[1][2] = z_ctl_sol;
    in_ptr[2][0] = x_ctl_sen;
    in_ptr[2][1] = y_ctl_sen;
    in_ptr[2][2] = z_ctl_sen;

    tlon = (float *) calloc(npix*nrow, sizeof (float));
    tlat = (float *) calloc(npix*nrow, sizeof (float));
    tsenz = (float *) calloc(npix*nrow, sizeof (float));
    tsena = (float *) calloc(npix*nrow, sizeof (float));
    tsolz = (float *) calloc(npix*nrow, sizeof (float));
    tsola = (float *) calloc(npix*nrow, sizeof (float));

    tout_ptr[0][0] = tlon;
    tout_ptr[0][1] = tlat;
    tout_ptr[1][0] = tsola;
    tout_ptr[1][1] = tsolz;
    tout_ptr[2][0] = tsena;
    tout_ptr[2][1] = tsenz;

    out_ptr[0][0] = lon;
    out_ptr[0][1] = lat;
    out_ptr[1][0] = sola;
    out_ptr[1][1] = solz;
    out_ptr[2][0] = sena;
    out_ptr[2][1] = senz;


    /* we now have all the info at each control point, so we
    can interpolate to all pixels, all lines */
    /* interpolate angles across each control point line */
    if (want_verbose)
        printf("Interpolating rows for longitude/azimuth\n");
    spl_aux = (float *) calloc(ncol, sizeof (float));

    for (i = 0; i < nrow; i++) {
        for (l = 0; l < 3; l++) {
            spline(xctl, in_ptr[l][0] + i*ncol, ncol, 1e30, 1e30, spl_aux);
            for (j = 0; j < nlon; j++)
                splint(xctl, in_ptr[l][0] + i * ncol, spl_aux, ncol,
                    (float) j, tout_ptr[l][0] + i * npix + j);

            spline(xctl, in_ptr[l][1] + i*ncol, ncol, 1e30, 1e30, spl_aux);
            for (j = 0; j < nlon; j++)
                splint(xctl, in_ptr[l][1] + i * ncol, spl_aux, ncol,
                    (float) j, tout_ptr[l][1] + i * npix + j);
        }
    }
    free(spl_aux);


    /* fill missing lines by interpolating columns */
    if (want_verbose)
        printf("Interpolating columns for longitude/azimuth\n");
    spl_aux = (float *) calloc(nrow, sizeof (float));

    for (i = 0; i < nlon; i++) {
        for (l = 0; l < 3; l++) {
            for (k = 0; k < nrow; k++) {
                in1[k] = *(tout_ptr[l][0] + k * npix + i);
                in2[k] = *(tout_ptr[l][1] + k * npix + i);
            }
            spline(yctl, in1, nrow, 1e30, 1e30, spl_aux);
            for (j = 0; j < nlat; j++)
                splint(yctl, in1, spl_aux, nrow, (float) j, (float *) &out1[j]);

            spline(yctl, in2, nrow, 1e30, 1e30, spl_aux);
            for (j = 0; j < nlat; j++)
                splint(yctl, in2, spl_aux, nrow, (float) j, (float *) &out2[j]);

            for (j = 0; j < nlat; j++) {
                *(out_ptr[l][0] + j * npix + i) = atan2(out2[j], out1[j]) * RADEG;
                if (l >= 1 && *(out_ptr[l][0] + j * npix + i) < 0) {
                    *(out_ptr[l][0] + j * npix + i) += 360;
                }
            }
        }
    }
    free(spl_aux);

    if (want_verbose)
        printf("Interpolating rows for latitude/zenith\n");
    spl_aux = (float *) calloc(ncol, sizeof (float));

    for (i = 0; i < nrow; i++) {
        for (l = 0; l < 3; l++) {
            spline(xctl, in_ptr[l][2] + i*ncol, ncol, 1e30, 1e30, spl_aux);
            for (j = 0; j < nlon; j++)
                splint(xctl, in_ptr[l][2] + i * ncol, spl_aux, ncol,
                    (float) j, tout_ptr[l][1] + i * npix + j);
        }
    }
    free(spl_aux);


    /* fill missing lines by interpolating columns */
    if (want_verbose)
        printf("Interpolating columns for latitude/zenith\n");
    spl_aux = (float *) calloc(nrow, sizeof (float));

    for (i = 0; i < nlon; i++) {
        for (l = 0; l < 3; l++) {
            for (k = 0; k < nrow; k++) {
                in1[k] = *(tout_ptr[l][1] + k * npix + i);
            }
            spline(yctl, in1, nrow, 1e30, 1e30, spl_aux);
            for (j = 0; j < nlat; j++)
                splint(yctl, in1, spl_aux, nrow, (float) j, (float *) &out1[j]);

            for (j = 0; j < nlat; j++) {
                *(out_ptr[l][1] + j * npix + i) = asin(out1[j]) * RADEG;
            }
        }
    }
    free(spl_aux);


    free(x_ctl_ll);
    free(y_ctl_ll);
    free(z_ctl_ll);
    free(x_ctl_sol);
    free(y_ctl_sol);
    free(z_ctl_sol);
    free(x_ctl_sen);
    free(y_ctl_sen);
    free(z_ctl_sen);

    free(inlon);
    free(inlat);
    free(insolz);
    free(insola);
    free(insenz);
    free(insena);
    free(xctl);
    free(yctl);
    free(in1);
    free(in2);

    free(tlon);
    free(tlat);
    free(tsolz);
    free(tsola);
    free(tsenz);
    free(tsena);

    return (0);
}

int
openl1_meris_CC(filehandle * file) {
    char *fltime;
    char *prodtype;
    char monthstr[10];
    static char months_list[12][4] ={"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL",
        "AUG", "SEP", "OCT", "NOV", "DEC"};
    size_t source_w, source_h;
    int32_t nscan;
    int fileID, xid, yid, retval;
    int minute, hour, month, sec, i;
    int year2, day2, msec2;
    float fsec, milisec;
    size_t att_len; // Change to size_t JMG 08/12/13
    double dblsec, starttime, stoptime;
    int16_t y16, m16, d16;
    size_t nctlpix = 0;
    size_t nctlscan = 0;

    // Open the CC input file 
    retval = nc_open(file->name, NC_NOWRITE, &fileID);
    if (retval != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
                __FILE__, __LINE__, file->name);
        return (1);
    }

    // Get pixel and scan dimensions
    retval = nc_inq_dimid(fileID, "x", &xid);
    retval = nc_inq_dimid(fileID, "y", &yid);
    retval = nc_inq_dimlen(fileID, xid, &source_w);
    retval = nc_inq_dimlen(fileID, yid, &source_h);

    if (want_verbose) {
        printf("MERIS CC Npix  :%d Nscans:%d\n", (int) source_w, (int) source_h);
    } // want_verbose
    npix = (int32_t) source_w;
    nscan = (int32_t) source_h;
    nline = nscan;

    // get year, day, msecs 
    retval = nc_inq_attlen(fileID, NC_GLOBAL, "start_date", &att_len);

    /* allocate required space before retrieving values */
    fltime = (char *) malloc(att_len + 1); /* + 1 for trailing null */

    /* get attribute values */
    retval = nc_get_att_text(fileID, NC_GLOBAL, "start_date", fltime);

    sscanf(fltime, "%02d-%3s-%04d %02d:%02d:%f", &day, monthstr, &year,
            &hour, &minute, &fsec);
    monthstr[4] = '\0';
    month = 1;
    for (i = 0; i < 12; i++)
        if (strncmp(monthstr, months_list[i], 3) == 0)
            month = i + 1;
    sec = (int) trunc((double) fsec);
    milisec = (float) (fsec - sec) * 1e3;
    ymdhms2ydmsec(year, month, day, hour, minute, sec,
            (int32_t *) & year, (int32_t *) & day, (int32_t *) & msec);
    msec += (int) milisec;

    dblsec = hour * 3600 + minute * 60 + fsec;
    y16 = year;
    m16 = 0; /* at this time, day is day-of-year */
    d16 = day;

    starttime = ymds2unix(y16, m16 + 1, d16, dblsec);


    retval = nc_inq_attlen(fileID, NC_GLOBAL, "stop_date", &att_len);

    /* allocate required space before retrieving values */
    fltime = (char *) malloc(att_len + 1); /* + 1 for trailing null */
    prodtype = (char *) malloc(att_len + 1); /* + 1 for trailing null */

    /* get attribute values */
    retval = nc_get_att_text(fileID, NC_GLOBAL, "stop_date", fltime);
    retval = nc_get_att_text(fileID, NC_GLOBAL, "product_type", prodtype);

    sscanf(fltime, "%02d-%3s-%04d %02d:%02d:%f", &day2, monthstr, &year2,
            &hour, &minute, &fsec);
    monthstr[4] = '\0';
    month = 1;
    for (i = 0; i < 12; i++)
        if (strncmp(monthstr, months_list[i], 3) == 0)
            month = i + 1;
    sec = (int) trunc((double) fsec);
    milisec = (float) (fsec - sec) * 1e3;
    ymdhms2ydmsec(year2, month, day2, hour, minute, sec,
            (int32_t *) & year2, (int32_t *) & day2, (int32_t *) & msec2);
    msec2 += (int) milisec;

    dblsec = hour * 3600 + minute * 60 + fsec;
    y16 = year2;
    m16 = 0; /* at this time, day is day-of-year */
    d16 = day2;

    stoptime = ymds2unix(y16, m16 + 1, d16, dblsec);
    time_interval = (stoptime - starttime) / (nscan - 1); /* in sec */

    file->sd_id = fileID;
    file->nbands = 15;
    file->npix = npix;
    file->nscan = nscan;
    if (strncmp(prodtype, "MER_RR", 6) == 0)
        strcpy(file->spatialResolution, "1.2 km");
    else
        strcpy(file->spatialResolution, "300 m");

    retval = nc_inq_dimid(file->sd_id, "tp_x", &xid);
    retval = nc_inq_dimlen(file->sd_id, xid, &nctlpix);
    ncol = nctlpix;

    retval = nc_inq_dimid(file->sd_id, "tp_y", &yid);
    retval = nc_inq_dimlen(file->sd_id, yid, &nctlscan);
    nrow = nctlscan;

    lon = (float *) calloc(npix*nline, sizeof (float));
    lat = (float *) calloc(npix*nline, sizeof (float));
    senz = (float *) calloc(npix*nline, sizeof (float));
    sena = (float *) calloc(npix*nline, sizeof (float));
    solz = (float *) calloc(npix*nline, sizeof (float));
    sola = (float *) calloc(npix*nline, sizeof (float));

    retval = navigation_meris(fileID);

    return (LIFE_IS_GOOD);
}

int readl1_meris_CC(filehandle *file, int32_t scan, l1str *l1rec)
/*
 * W. Robinson, SAIC, 22 May 2012  account for msec going to next day
 */ {
    static int firstCall = 1;

    // only look at SeaWiFS bands until we have rayleigh/aerosol tables
    // NOTE: the following names are part of the API and NOT the file.  See libmeris.
    static char *l1_names[] = {
        "radiance_1", "radiance_2", "radiance_3", "radiance_4", // 412.7, 442.6, 489.9, 509.8
        "radiance_5", "radiance_6", "radiance_7", "radiance_8", // 559.7, 619.6, 664.5, 680.8
        "radiance_9", "radiance_10", "radiance_11", "radiance_12", // 708.3, 753.4, 761.5, 778.4
        "radiance_13", "radiance_14", "radiance_15" // 864.9, 884.9, 900.0
    };
    int32_t nbands;
    int32_t ip, ib, ipb;
    size_t start[3], count[3];
    int band_id, retval;
    int16_t *rad_data;
    static float gain[15];
    double recsec, sec70;
    int16_t sh_year, sh_day;

    nbands = file->nbands;

    if (firstCall) {
        if (want_verbose)
            printf("file->nbands = %d\n", (int) file->nbands);
        firstCall = 0;

        for (ip = 0; ip < npix; ip++) l1rec->pixnum[ip] = ip;

        for (ib = 0; ib < nbands; ib++) {
            retval = nc_inq_varid(file->sd_id, l1_names[ib], &band_id);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: nc_inq_varid failed for file, %s  band, %s.\n",
                        __FILE__, __LINE__, file->name, l1_names[ib]);
                return (1);
            }

            retval = nc_get_att_float(file->sd_id, band_id, "scale_factor", &gain[ib]);
            if (retval != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: nc_get_att_float for file, %s  attribute, %s.\n",
                        __FILE__, __LINE__, file->name, "scale_factor");
                return (1);
            }
        }
    }

    /*  set time for this scan and account for overflow of msec to next day  */
    recsec = (double) (msec + time_interval * scan * 1e3) / 1.e3;
    sec70 = ymds2unix((int16_t) year, 1, (int16_t) day, recsec);
    unix2yds(sec70, &sh_year, &sh_day, &recsec);
    l1rec->scantime = sec70;

    for (ip = spix; ip < npix; ip++) {
        l1rec->lon[ip] = lon[scan * npix + ip];
        l1rec->lat[ip] = lat[scan * npix + ip];
        l1rec->solz[ip] = solz[scan * npix + ip];
        l1rec->sola[ip] = sola[scan * npix + ip];
        l1rec->senz[ip] = senz[scan * npix + ip];
        l1rec->sena[ip] = sena[scan * npix + ip];
    }


    // read in radiance data
    rad_data = (int16_t *) calloc(npix, sizeof (int16_t));

    start[0] = scan;
    start[1] = 0;
    count[0] = 1;
    count[1] = npix;

    for (ib = 0; ib < nbands; ib++) {

        retval = nc_inq_varid(file->sd_id, l1_names[ib], &band_id);
        retval = nc_get_vara_short(file->sd_id, band_id, start, count, rad_data);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_short failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, file->name, l1_names[ib]);
            return (1);
        }
        // copy to Lt record.
        for (ip = spix; ip < npix; ip++) {
            ipb = ip * nbands + ib;
            l1rec->Lt[ipb] = rad_data[ip] * gain[ib] / 10.0;

            // mark negative input data as HILT
            if (l1rec->Lt[ipb] < 0.0)
                l1rec->navfail[ip] = 1;
        }
    } // for ib
    free(rad_data);

    // disable smile correction
    l1_input->rad_opt = 0;

    l1rec->npix = file->npix;

    return (LIFE_IS_GOOD);
}

int
closel1_meris_CC(filehandle *file) {
    int retval;

    retval = nc_close(file->sd_id);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_close failed for file, %s.\n",
                __FILE__, __LINE__, file->name);
        return (1);
    }

    return (LIFE_IS_GOOD);
}

