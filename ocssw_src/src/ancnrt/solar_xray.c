/***************************************************************** 
 * File:        solar_xray.c
 *
 * Modification history:        
 *****************************************************************/
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include "hdf5.h"

double ymds2unix(short year, short month, short day, double secs);
void unix2yds(double usec, short *year, short *day, double *secs);

int main(int argc, char *argv[]) {

    int32_t i, j, k;

    int16_t yr, mon, dy, doy;
    double secs;
    float flux[288 * 5];
    char junk[120], outfile[80], buffer[80];
    FILE *fp;

    hsize_t dimsf[1]; /* dataset dimensions */

    /* 
     * Data  and output buffer initialization. 
     */
    hid_t file, dataset; /* handles */
    hid_t dataspace;

    strcpy(junk, basename(argv[1]));

    memcpy(buffer, &junk[0], 4);
    buffer[4] = 0;
    yr = (int16_t) atoi(buffer);

    memcpy(buffer, &junk[4], 2);
    buffer[2] = 0;
    mon = (int16_t) atoi(buffer);

    memcpy(buffer, &junk[6], 2);
    buffer[2] = 0;
    dy = (int16_t) atoi(buffer);

    unix2yds(ymds2unix(yr, mon, dy, 0.0), &yr, &doy, &secs);

    sprintf(outfile, "%s%d%03d%s", "N", yr, doy, "00_XRAY_GOES_24h.h5");

    printf("Creating: %s\n", outfile);

    file = H5Fcreate(outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    dimsf[0] = 288 * 5;
    dataspace = H5Screate_simple(1, dimsf, NULL);


    fp = fopen(argv[1], "r");
    for (i = 0; i < 19; i++) fgets(junk, 120, fp);

    for (i = 0; i < 288 * 5; i++) {
        fgets(junk, 120, fp);
        memcpy(buffer, &junk[45], 15);
        buffer[15] = 0;
        flux[i] = (float) atof(buffer);
    }


    /*
     * Create a new dataset within the file using defined dataspace and
     * default dataset creation properties.
     */
    dataset = H5Dcreate1(file, "1min_xray_irradiance", H5T_NATIVE_FLOAT,
            dataspace, H5P_DEFAULT);

    /*
     * Write the data to the dataset using default transfer properties.
     */
    H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, flux);


    /*
     * Close/release resources.
     */
    H5Sclose(dataspace);
    H5Dclose(dataset);
    H5Fclose(file);


    return 0;
}
