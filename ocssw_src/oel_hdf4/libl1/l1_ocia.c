/* ============================================================================ */
/* module l1_ocia.c - functions to read OCIA (coastal color) for MSL12  */
/* Written By:  Rick Healy GSFC SAIC, July, 2016.                       */
/*                                                                              */
/* ============================================================================ */

#include "l1.h"

#include <netcdf.h>
#include <stdlib.h>

static int16_t nline, npix;

static int32_t spix = 0;
static double starttime;
static double *scan_start_tai, *scan_stop_tai;
static double time_interval;
static float *alt;

int openl1_ocia(filehandle * file) {

    size_t source_w, source_h, source_b;
    int32_t nscan;
    int fileID, xid, yid, band_id, retval, grp_id, sds_id;
    double stoptime;
    size_t start[3], count[3];

    // Open the netcdf4 input file
    retval = nc_open(file->name, NC_NOWRITE, &fileID);
    if (retval != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
                __FILE__, __LINE__, file->name);
        return (1);
    }

    // Get pixel and scan dimensions
    retval = nc_inq_dimid(fileID, "number_of_scans", &yid);
    retval = nc_inq_dimid(fileID, "number_of_pixels", &xid);
    retval = nc_inq_dimid(fileID, "number_of_bands", &band_id);
    retval = nc_inq_dimlen(fileID, xid, &source_w);
    retval = nc_inq_dimlen(fileID, yid, &source_h);
    retval = nc_inq_dimlen(fileID, band_id, &source_b);

    if (want_verbose) {
        printf("OCIA L1B Npix  :%d Nscans:%d\n", (int) source_w,
                (int) source_h);
    } // want_verbose
    npix = (int32_t) source_w;
    nscan = (int32_t) source_h;
    nline = nscan;

    file->sd_id = fileID;
    file->npix = npix;
    file->nscan = nscan;
    file->nbands = source_b;
    strcpy(file->spatialResolution, "750 m");

    scan_start_tai = (double *) calloc(nscan, sizeof (double));
    scan_stop_tai = (double *) calloc(nscan, sizeof (double));
    alt = (float *) calloc(nscan, sizeof (float));

    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    count[0] = nscan;
    count[1] = 0;
    count[2] = 0;

    nc_inq_ncid(file->sd_id, "scan_line_attributes", &grp_id);

    retval = nc_inq_varid(grp_id, "altitude", &sds_id);
    retval = nc_get_vara_float(grp_id, sds_id, start, count, alt);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "altitude");
        return (1);
    }

    retval = nc_inq_varid(grp_id, "scan_start_time", &sds_id);
    retval = nc_get_vara_double(grp_id, sds_id, start, count, scan_start_tai);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "scan_start_time");
        return (1);
    }

    int iscan = 0, i1, i2;
    starttime = 9999;
    while (iscan < nscan) {
        if (scan_start_tai[iscan] > 0) {
            if (starttime > scan_start_tai[iscan])starttime = scan_start_tai[iscan];
        }
        iscan++;
    }
    i1 = iscan;

    //starttime = scan_start_tai[0];

    //    nc_inq_varid(grp_id,"scan_end_time",&sds_id);
    //    nc_get_vara_double(grp_id, sds_id, start, count, scan_stop_tai);
    iscan = nscan - 1;
    stoptime = -9999;
    while (iscan > 0) {
        if (scan_start_tai[iscan] > 0) {
            if (stoptime < scan_start_tai[iscan])stoptime = scan_start_tai[iscan];
        }
        iscan--;
    }
    i2 = iscan;
    //stoptime = scan_start_tai[nscan];

    time_interval = (stoptime - starttime) / ((i2 - i1) + 1); /* in sec */

    /*
        retval = navigation_meris(fileID);
     */
    start[0] = 0;
    count[0] = source_b;
    return (LIFE_IS_GOOD);
}

int readl1_ocia(filehandle *file, int32_t scan, l1str *l1rec) {
    static int firstCall = 1;
    int retval;
    float *lon, *lat, *senz, *sena, *solz, *sola;

    int32_t ip, ib, ipb, Ibp;
    int32_t nbands = l1rec->l1file->nbands;
    size_t start[3], count[3];
    int band_id, grp_id;
    float *rad_data;

    if (firstCall) {
        if (want_verbose)
            printf("file->nbands = %d, l1rec->nbands = %d\n",
                (int) file->nbands, (int) l1rec->l1file->nbands);
        firstCall = 0;

        for (ip = 0; ip < npix; ip++) {
            l1rec->pixnum[ip] = ip;
            l1rec->flags[ip] = 0;
        }

    }


    //      set time for this scan - if scan_start_time value not properly set, estimate from scene start time.
    // printf("RJH: scan_start_tai=%f\n",scan_start_tai[scan]);
    //    if (scan_start_tai[scan] > 0){
    //        lastvalidtime = scan_start_tai[scan]+tai2unixoffset;
    //        lastvalidscan = scan;
    //        unix2yds(lastvalidtime, &scan_year, &scan_day, &msec);
    //    } else {
    //        unix2yds(lastvalidtime+(time_interval * (scan-lastvalidscan)),&scan_year,&scan_day, &msec);
    //    }
    //    msec -= leapseconds_since_1993(lastvalidtime-tai2unixoffset)*1000;
    //    l1rec->scantime = yds2unix(scan_year, scan_day, (double) (msec / 1.e3));
    l1rec->scantime = scan_start_tai[scan];

    retval = nc_inq_ncid(file->sd_id, "geolocation_data", &grp_id);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_inq_ncid failed for file, %s  group, %s.\n",
                __FILE__, __LINE__, file->name, "geolocation_data");
        return (1);
    }

    start[0] = scan;
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = npix;
    count[2] = 0;

    lon = (float *) calloc(npix, sizeof (float));
    lat = (float *) calloc(npix, sizeof (float));
    senz = (float *) calloc(npix, sizeof (float));
    sena = (float *) calloc(npix, sizeof (float));
    solz = (float *) calloc(npix, sizeof (float));
    sola = (float *) calloc(npix, sizeof (float));

    retval = nc_inq_varid(grp_id, "longitude", &band_id);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "longitude");
        return (1);
    }
    retval = nc_get_vara_float(grp_id, band_id, start, count, lon);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "longitude");
        return (1);
    }
    retval = nc_inq_varid(grp_id, "latitude", &band_id);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "latitude");
        return (1);
    }
    retval = nc_get_vara_float(grp_id, band_id, start, count, lat);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "latitude");
        return (1);
    }
    //until geolocation is read, set fill values -
    for (ip = spix; ip < npix; ip++) {
        l1rec->lon[ip] = lon[ip];
        l1rec->lat[ip] = lat[ip];
        l1rec->solz[ip] = -999; //solz[scan * npix + ip];
        l1rec->sola[ip] = -999; //sola[scan * npix + ip];
        l1rec->senz[ip] = -999; //senz[scan * npix + ip];
        l1rec->sena[ip] = -999; //sena[scan * npix + ip];
    }

    retval = nc_inq_varid(grp_id, "sensor_azimuth", &band_id);
    retval = nc_get_vara_float(grp_id, band_id, start, count, sena);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "sensor_azimuth");
        return (1);
    }
    retval = nc_inq_varid(grp_id, "sensor_zenith", &band_id);
    retval = nc_get_vara_float(grp_id, band_id, start, count, senz);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "sensor_zenith");
        return (1);
    }
    retval = nc_inq_varid(grp_id, "solar_azimuth", &band_id);
    retval = nc_get_vara_float(grp_id, band_id, start, count, sola);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "solar_azimuth");
        return (1);
    }
    retval = nc_inq_varid(grp_id, "solar_zenith", &band_id);
    retval = nc_get_vara_float(grp_id, band_id, start, count, solz);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "solar_zenith");
        return (1);
    }

    for (ip = spix; ip < npix; ip++) {
        l1rec->lon[ip] = lon[ip];
        l1rec->lat[ip] = lat[ip];
        l1rec->solz[ip] = solz[ip];
        l1rec->sola[ip] = sola[ip];
        l1rec->senz[ip] = senz[ip];
        l1rec->sena[ip] = sena[ip];
    }

    // read in radiance data
    rad_data = (float *) calloc(npix*nbands, sizeof (float));

    start[0] = 0;
    start[1] = scan;
    start[2] = 0;
    count[0] = nbands;
    count[1] = 1;
    count[2] = npix;

    nc_inq_ncid(file->sd_id, "observation_data", &grp_id);

    retval = nc_inq_varid(grp_id, "Lt", &band_id);
    retval = nc_get_vara_float(grp_id, band_id, start, count, rad_data);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "Lt_visnir");
        return (1);
    }

    for (ib = 0; ib < nbands; ib++) {

        // copy to Lt record.
        for (ip = spix; ip < npix; ip++) {
            ipb = ip * nbands + ib;
            Ibp = ib * npix + ip;
            l1rec->Lt[ipb] = rad_data[Ibp];

            // mark negative input data as HILT
            // navfail commented out for Lt < 0 - too limiting for hyperspectral - RJH
            if (l1rec->Lt[ipb] < 0.0) {
                l1rec->Lt[ipb] = 0.0001;
                //                l1rec->navfail[ip] = 1;
            }
        }
    } // for ib

    l1rec->alt = alt[scan] / 1000.; //convert to km

    free(rad_data);
    free(lat);
    free(lon);
    free(solz);
    free(sola);
    free(senz);
    free(sena);

    l1rec->l1file->sensorID = file->sensorID;
    l1rec->npix = file->npix;

    return (LIFE_IS_GOOD);
}

int closel1_ocia(filehandle *file) {
    int retval;

    retval = nc_close(file->sd_id);
    if (retval != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_close failed for file, %s.\n",
                __FILE__, __LINE__, file->name);
        return (1);
    }
    free(scan_start_tai);

    return (LIFE_IS_GOOD);
}

