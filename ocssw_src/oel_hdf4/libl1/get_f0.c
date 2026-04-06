#include <stdio.h>
#include <stdlib.h>
#include "l1.h"
#include "netcdf.h"
#include "genutils.h"
#include "math.h"
static int luts_set = 0;
static float *wavelengths;
static float *irradiance;
static size_t wave_len;
double d_wave_step;
double starting_wave;
double ending_wave;

int set_solar_irradiance(const char *f0_filename) {
    if (luts_set)
        return 0;
    FILE *file;
    file = fopen(f0_filename, "r");
    if (file == NULL) {
        printf("-E- %s:%d - unable to open %s for reading\n", __FILE__, __LINE__, f0_filename);
        return 1;
    }
    fclose(file);
    int status, ncid, varid;
    status = nc_open(f0_filename, NC_NOWRITE, &ncid);
    if (status != NC_NOERR) {
        return 1;
    }
    status = nc_inq_varid(ncid, "wavelength", &varid);
    if (status != NC_NOERR) {
        return 1;
    }
    int ndims;
    int dimids[NC_MAX_DIMS];

    status = nc_inq_var(ncid, varid, NULL, NULL, &ndims, dimids, NULL);
    if (status != NC_NOERR) {
        return 1;
    }
    if (ndims != 1) {
        return 1;
    }
    status = nc_inq_dimlen(ncid, dimids[0], &wave_len);
    if (status != NC_NOERR) {
        return 1;
    }
    wavelengths = (float *)malloc(wave_len * sizeof(float));
    if (wavelengths == NULL) {
        return 1;
    }
    status = nc_get_var_float(ncid, varid, wavelengths);
    if (status != NC_NOERR) {
        return 1;
    }
    status = nc_inq_varid(ncid, "irradiance", &varid);
    if (status != NC_NOERR) {
        return 1;
    }

    status = nc_inq_var(ncid, varid, NULL, NULL, &ndims, dimids, NULL);
    if (status != NC_NOERR) {
        return 1;
    }
    if (ndims != 1) {
        return 1;
    }
    size_t temp_len;
    status = nc_inq_dimlen(ncid, dimids[0], &temp_len);
    if (status != NC_NOERR) {
        return 1;
    }
    if (temp_len != wave_len) {
        return 1;
    }
    irradiance = (float *)malloc(wave_len * sizeof(float));
    if (irradiance == NULL) {
        return 1;
    }
    status = nc_get_var_float(ncid, varid, irradiance);
    if (status != NC_NOERR) {
        return 1;
    }
    luts_set = 1;
    nc_close(ncid);
    // only works for a uniform wavelength grid
    d_wave_step = wavelengths[1] - wavelengths[0];
    starting_wave = wavelengths[0];
    ending_wave = wavelengths[wave_len - 1];
    for (int iw = 0; iw < wave_len; iw++) {
        irradiance[iw] *= 100.0; // W/m2/nm -> mW/cm^2/um
    }
    return 0;
}

int get_f0(float wave, float width, float *f0) {
    if (luts_set == 0)
        return 1;
    if (wave < starting_wave || wave > ending_wave) {
        printf( "-E- %s:%d: wavelength of %f outside irradiance table range.\n", __FILE__, __LINE__,
                wave);
        exit(EXIT_FAILURE);
    }
    size_t i = (wave - starting_wave) / d_wave_step;
    size_t imin = round(MAX(i - width / 2 / d_wave_step, 0));
    size_t imax = MIN(i + width / 2 / d_wave_step, wave_len - 1); // in case of width = 1 we want `imin == imax`
    float sum = 0.0;
    for (i = imin; i <= imax; i++)
        sum += irradiance[i];
    sum /= (imax - imin + 1);
    *f0 = sum;
    return 0;
}
