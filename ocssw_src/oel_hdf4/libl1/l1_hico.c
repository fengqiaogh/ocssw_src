/*
 * l1_hico_hdf.c
 *
 */

#include "l1_hico.h"

#include <string.h>
#include "h5io.h"
#include "l1.h"

#define NUM_GEO_DATA 6

static char *geo_name[] = {"latitude", "longitude", "sensor_azimuth",
    "sensor_zenith", "solar_azimuth", "solar_zenith"};

// HICO file private information stored in filehandle

typedef struct hico_private_struct {
    int numBands; // number of bands in the hico file
    int* wave_ix; // array that maps the sensor bands into the hico file
    int32_t syear;
    int32_t sday;
    int32_t smsec;
    float lt_slope;
    float lt_intercept;
    uint16_t *rad_data;
    h5io_str* fileID;
    h5io_str* ds_id;
    h5io_str* geo_dat_id[NUM_GEO_DATA];
    int orientation; // is HICO orientation flipped
    float *fwhm;
} hico_private_t;

static hico_private_t* allocatePrivateData() {
    int i;
    hico_private_t* pData = (hico_private_t*) calloc(1, sizeof (hico_private_t));
    pData->lt_slope = 1.0;
    pData->fileID = (h5io_str*) calloc(1, sizeof (h5io_str));
    pData->ds_id = (h5io_str*) calloc(1, sizeof (h5io_str));
    for (i = 0; i < NUM_GEO_DATA; i++)
        pData->geo_dat_id[i] = (h5io_str*) calloc(1, sizeof (h5io_str));
    pData->orientation = 0;
    return pData;
}

static void freePrivateData(hico_private_t* pData) {
    int i;

    if (pData->rad_data)
        free(pData->rad_data);
    if (pData->wave_ix)
        free(pData->wave_ix);
    free(pData->fileID);
    free(pData->ds_id);
    for (i = 0; i < NUM_GEO_DATA; i++)
        free(pData->geo_dat_id[i]);
    free(pData);
}

int openl1_hico(filehandle * file) {
    h5io_str gid;
    int year, month, day, hour, minute;
    int ids;
    int ndim;
    int dims[3];
    int sto_len;
    int sec;
    char sdate[16];
    char edate[16];
    char stime[16];
    char etime[16];
    char g_path[100];
    H5T_class_t f_class;
    hid_t f_native_typ;
    float* wvls;
    hico_private_t* pData;
    int32_t *Lambda;
    int i;
    char orientationStr[10];

    pData = allocatePrivateData();
    file->private_data = pData;

    // get rid of the silly HDF5 error printing
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    
    if (h5io_openr(file->name, 0, pData->fileID)) {
        fprintf(stderr, "-E- %s Line %d: Failure to open %s\n",
                __FILE__, __LINE__, file->name);
        freePrivateData(pData);
        return 1;
    }

    // set Lt dataset pointer
    if (h5io_set_grp(pData->fileID, "products", &gid) != 0) {
        fprintf(stderr, "-E- %s Line %d: Unable to open group products\n",
                __FILE__, __LINE__);
        h5io_close(pData->fileID);
        freePrivateData(pData);
        return 1;
    }

    if (h5io_set_ds(&gid, "Lt", pData->ds_id) != 0) {
        fprintf(stderr, "-E- %s Line %d: Unable to open data set: Lt\n",
                __FILE__, __LINE__);
        h5io_close(&gid);
        h5io_close(pData->fileID);
        freePrivateData(pData);
        return 1;
    }
    h5io_close(&gid);

    // get number of wavelengths in file
    if (h5io_info(pData->ds_id, NULL, &f_class, &f_native_typ, &ndim, dims, &sto_len)
            != 0) {
        fprintf(stderr, "-E- %s Line %d:  Unable to get data set dimensions\n",
                __FILE__, __LINE__);
        h5io_close(pData->ds_id);
        h5io_close(pData->fileID);
        freePrivateData(pData);
        return 1;
    }
    pData->numBands = dims[2];
    // allocate space for wavelength array
    pData->fwhm = (float*) malloc(sizeof (float)*pData->numBands);
    if (pData->fwhm == NULL) {
        fprintf(stderr,
                "-E- %s Line %d: Unable to allocate data for wavelength bandwidth\n",
                __FILE__, __LINE__);
        h5io_close(pData->ds_id);
        h5io_close(pData->fileID);
        freePrivateData(pData);
        return 1;
    }
    wvls = (float*) malloc(sizeof (float) * pData->numBands);
    if (wvls == NULL) {
        fprintf(stderr,
                "-E- %s Line %d: Unable to allocate data for wavelengths\n",
                __FILE__, __LINE__);
        h5io_close(pData->ds_id);
        h5io_close(pData->fileID);
        freePrivateData(pData);
        return 1;
    }

    // Get the wavelength bandwidth list
    if (h5io_rd_attr(pData->ds_id, "fwhm", (void *) pData->fwhm) != 0) {
        fprintf(stderr,
                "-E- %s Line %d: Unable to read the wavelength bandwidth (fwhm) attribute\n",
                __FILE__, __LINE__);
        free(pData->fwhm);
        h5io_close(pData->ds_id);
        h5io_close(pData->fileID);
        freePrivateData(pData);
        return 1;
    }
    // Get the wavelength list
    if (h5io_rd_attr(pData->ds_id, "wavelengths", (void *) wvls) != 0) {
        fprintf(stderr,
                "-E- %s Line %d: Unable to read the wavelength attribute\n",
                __FILE__, __LINE__);
        free(wvls);
        h5io_close(pData->ds_id);
        h5io_close(pData->fileID);
        freePrivateData(pData);
        return 1;
    }

    // allocate space for wavelength index array
    pData->wave_ix = (int*) malloc(sizeof (int) * file->nbands);
    if (pData->wave_ix == NULL) {
        fprintf(stderr,
                "-E- %s Line %d: Unable to allocate data for wavelength index array\n",
                __FILE__, __LINE__);
        free(wvls);
        h5io_close(pData->ds_id);
        h5io_close(pData->fileID);
        freePrivateData(pData);
        return 1;
    }

    // fill up the wavelength index array
    rdsensorinfo(file->sensorID, l1_input->evalmask, "Lambda", (void **) &Lambda);
    for (i = 0; i < file->nbands; i++) {
        pData->wave_ix[i] = windex(Lambda[i], wvls, pData->numBands);
    }
    free(wvls);

    // Get the start and end date,time
    if (h5io_set_grp(pData->fileID,
            "metadata/FGDC/Identification_Information/Time_Period_of_Content",
            &gid) != 0) {
        fprintf(stderr,
                "-E- %s Line %d: Unable to open group Time_Period_of_Content\n",
                __FILE__, __LINE__);
        h5io_close(pData->ds_id);
        h5io_close(pData->fileID);
        freePrivateData(pData);
        return 1;
    }

    bzero(&sdate, sizeof (sdate));
    if (h5io_rd_attr(&gid, "Beginning_Date", (void *) &sdate) != 0) {
        fprintf(stderr,
                "-E- %s Line %d: Unable to read the Beginning_Date attribute:\n",
                __FILE__, __LINE__);
        h5io_close(&gid);
        h5io_close(pData->ds_id);
        h5io_close(pData->fileID);
        freePrivateData(pData);
        return 1;
    }

    bzero(&stime, sizeof (stime));
    if (h5io_rd_attr(&gid, "Beginning_Time", (void *) &stime) != 0) {
        fprintf(stderr,
                "-E- %s Line %d: Unable to read the Beginning_Date attribute:\n",
                __FILE__, __LINE__);
        h5io_close(&gid);
        h5io_close(pData->ds_id);
        h5io_close(pData->fileID);
        freePrivateData(pData);
        return 1;
    }

    bzero(&edate, sizeof (edate));
    if (h5io_rd_attr(&gid, "Ending_Date", (void *) &edate) != 0) {
        fprintf(stderr,
                "-E- %s Line %d: Unable to read the Beginning_Date attribute:\n",
                __FILE__, __LINE__);
        h5io_close(&gid);
        h5io_close(pData->ds_id);
        h5io_close(pData->fileID);
        freePrivateData(pData);
        return 1;
    }

    bzero(&etime, sizeof (etime));
    if (h5io_rd_attr(&gid, "Ending_Time", (void *) &etime) != 0) {
        fprintf(stderr,
                "-E- %s Line %d: Unable to read the Beginning_Date attribute:\n",
                __FILE__, __LINE__);
        h5io_close(&gid);
        h5io_close(pData->ds_id);
        h5io_close(pData->fileID);
        freePrivateData(pData);
        return 1;
    }
    h5io_close(&gid);

    // Parse the date/time strings
    sscanf(sdate, "%04d%02d%02d", &year, &month, &day);
    sscanf(stime, "%02d%02d%02d", &hour, &minute, &sec);
    ymdhms2ydmsec(year, month, day, hour, minute, sec, &pData->syear, &pData->sday, &pData->smsec);

    if (h5io_rd_attr(pData->ds_id, "scale_factor", (void *) &pData->lt_slope) != 0) {
        if (h5io_rd_attr(pData->ds_id, "slope", (void *) &pData->lt_slope) != 0) {
            fprintf(stderr,
                    "-E- %s Line %d: Unable to read the Lt slope attribute:\n",
                    __FILE__, __LINE__);
            h5io_close(pData->ds_id);
            h5io_close(pData->fileID);
            freePrivateData(pData);
            return 1;
        }
    }

    if (h5io_rd_attr(pData->ds_id, "add_offset", (void *) &pData->lt_intercept) != 0) {
        if (h5io_rd_attr(pData->ds_id, "intercept", (void *) &pData->lt_intercept) != 0) {
            fprintf(stderr,
                    "-E- %s Line %d: Unable to read the Lt intercept attribute:\n",
                    __FILE__, __LINE__);
            h5io_close(pData->ds_id);
            h5io_close(pData->fileID);
            freePrivateData(pData);
            return 1;
        }
    }

    //  Set to the geolocation datasets for the lat, lon, and view angles
    //  TODO: Get real with the navigation code ;)

    for (ids = 0; ids < NUM_GEO_DATA; ids++) {
        sprintf(g_path, "navigation/%s", geo_name[ids]);
        if (h5io_set_ds(pData->fileID, g_path, pData->geo_dat_id[ids]) != 0) {

            // look for latitude/longitude and latitudes/longitudes
            if(ids == 0 || ids == 1) {
                sprintf(g_path, "navigation/%ss", geo_name[ids]);
                if (h5io_set_ds(pData->fileID, g_path, pData->geo_dat_id[ids]) != 0) {
                    fprintf(stderr,
                            "-E- %s Line %d:  Unable to set ds # %d in geolocation file\n",
                            __FILE__, __LINE__, ids);
                    h5io_close(pData->ds_id);
                    h5io_close(pData->fileID);
                    freePrivateData(pData);
                    return 1;
                }
            }
        }
    }

    // Get the HICO orientation
    if (h5io_set_grp(pData->fileID, "metadata/HICO/Calibration", &gid) != 0) {
        fprintf(stderr,
                "-E- %s Line %d: Unable to open group metadata/HICO/Calibration\n",
                __FILE__, __LINE__);
        h5io_close(pData->ds_id);
        h5io_close(pData->fileID);
        freePrivateData(pData);
        return 1;
    }

    // #### note that there was a space in the attribute name
    // #### how annoying.
    bzero(&orientationStr, sizeof (orientationStr));
    if (h5io_rd_attr(&gid, "hico_orientation_from_quaternion", (void *) &orientationStr) != 0) {
        if (h5io_rd_attr(&gid, "hico_orientation_from_quaternion ", (void *) &orientationStr) != 0) {
            fprintf(stderr,
                    "-E- %s Line %d: Unable to read the hico_orientation_from_quaternion attribute:\n",
                    __FILE__, __LINE__);
            h5io_close(&gid);
            h5io_close(pData->ds_id);
            h5io_close(pData->fileID);
            freePrivateData(pData);
            return 1;
        }
    }
    h5io_close(&gid);

    if (strstr(orientationStr, "-XVV")) {
        pData->orientation = 1;
    } else {
        pData->orientation = 0;
    }

    file->npix = dims[1];
    file->nscan = dims[0];
    file->sensorID = HICO;
    strcpy(file->spatialResolution, "100 m");

    pData->rad_data = (uint16_t *) calloc(file->npix, sizeof (uint16_t));

    return (LIFE_IS_GOOD);
}

int readl1_hico(filehandle *file, int32_t scan, l1str *l1rec, int lonlat) {
    int igeo;
    int start[3], count[3];
    int32_t ib, ip, ipb;
    float *iptr;
    int32_t dmsec = 13;
    hico_private_t* pData = (hico_private_t*) file->private_data;

    l1rec->scantime = yds2unix(pData->syear, pData->sday, (double) ((pData->smsec + scan * dmsec) / 1.e3));

    // get the location and view information
    for (igeo = 0; igeo < NUM_GEO_DATA; igeo++) {

        // when lonlat is on, only need to get lon, lat and solar zenith
        if (lonlat && (igeo > 1 && igeo != 5))
            continue;

        switch (igeo) {
        case 0:
            iptr = l1rec->lat;
            break;
        case 1:
            iptr = l1rec->lon;
            break;
        case 2:
            iptr = l1rec->sena;
            break;
        case 3:
            iptr = l1rec->senz;
            break;
        case 4:
            iptr = l1rec->sola;
            break;
        case 5:
            iptr = l1rec->solz;
            break;
        }

        if (pData->orientation) {
            start[0] = file->nscan - scan - 1;
        } else {
            start[0] = scan;
        }
        start[1] = 0;
        count[0] = 1;
        count[1] = file->npix;
        if (h5io_rd_ds_slice(pData->geo_dat_id[igeo], start, count, (void *) iptr)
                != 0) {
            fprintf(stderr,
                    "-E- %s Line %d:  Failed to geonav scan %d of parm %d\n",
                    __FILE__, __LINE__, scan, igeo);
            return 1;
        }
    }
    
    // if latlon mode, skip the rest
    if (lonlat)
        return (LIFE_IS_GOOD);

    // read in radiance data
    if (pData->orientation) {
        start[0] = file->nscan - scan - 1;
    } else {
        start[0] = scan;
    }
    start[1] = 0;
    count[0] = 1;
    count[1] = file->npix;
    count[2] = 1;

    for (ib = 0; ib < l1rec->l1file->nbands; ib++) {
        start[2] = pData->wave_ix[ib];
        if (h5io_rd_ds_slice(pData->ds_id, start, count, (void *) pData->rad_data) != 0) {
            fprintf(stderr,
                    "-E- %s Line %d:  Failed to read scan %d of band %d\n",
                    __FILE__, __LINE__, scan, ib);
            return 1;
        }

        for (ip = 0; ip < file->npix; ip++) {
            ipb = ip * l1rec->l1file->nbands + ib;

            if (pData->rad_data[ip] > 0) {
                l1rec->Lt[ipb] = ((float) (pData->rad_data[ip]) * pData->lt_slope
                        + pData->lt_intercept) / 10.0;
            } else {
                l1rec->Lt[ipb] = BAD_FLT;
                //                l1rec->navfail[ip] = 1;
            }
        }
    }

    return (LIFE_IS_GOOD);
}

int closel1_hico(filehandle *file) {
    int i;
    hico_private_t* pData = (hico_private_t*) file->private_data;

    for (i = 0; i < NUM_GEO_DATA; i++) {
        if (h5io_close(pData->geo_dat_id[i]) != 0) {
            fprintf(stderr,
                    "-E- %s Line %d:  Failed to close navigation data %s\n",
                    __FILE__, __LINE__, geo_name[i]);
            return 1;
        }
    }

    if (h5io_close(pData->ds_id) != 0) {
        fprintf(stderr, "-E- %s Line %d:  Failed to close Lt data set\n",
                __FILE__, __LINE__);
        return 1;
    }

    if (h5io_close(pData->fileID) != 0) {
        fprintf(stderr, "-E- %s Line %d:  Failed to close HDF5 file %s\n",
                __FILE__, __LINE__, file->name);
        return 1;
    }

    freePrivateData(pData);
    file->private_data = NULL;

    return (LIFE_IS_GOOD);
}
