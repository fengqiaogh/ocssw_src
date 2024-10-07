#include <stdint.h>
#include <string.h>
#include "hdf5.h"

// wget http://nomads.ncep.noaa.gov/pub/data/nccf/com/wave/prod/wave.yyyymmdd/nww3.t00z.grib.grib2
// wgrib2 -match HTSGW:surface:anl nww3.t00z.grib.grib2 -bin swh.dat
// swh swh.dat Nyyyydoyhr_SWH_NCEP_6h.h5

int main(int argc, char *argv[]) {
    hsize_t dimsf[2]; /* dataset dimensions */

    /* 
     * Data  and output buffer initialization. 
     */
    hid_t file, dataset; /* handles */
    hid_t dataspace;

    float swh[157 * 288], buf[288];

    int i, j;
    uint32_t hdr;

    FILE *fp;

    // Read in binary swh
    fp = fopen(argv[1], "rb");
    fread(&hdr, sizeof (uint32_t), 1, fp);
    fread(swh, 157 * 288 * sizeof (float), 1, fp);
    fclose(fp);

    for (i = 0; i < 157 * 288; i++) if (swh[i] > 1e20) swh[i] = -999.;

    // Flip about equator (if not already flipped)  JMG 03/18/14
    if (swh[39] == -999) {
        for (i = 0; i < 157 / 2; i++) {
            memcpy(buf, &swh[288 * i], 288 * sizeof (float));
            memcpy(&swh[288 * i], &swh[288 * (156 - i)], 288 * sizeof (float));
            memcpy(&swh[288 * (156 - i)], buf, 288 * sizeof (float));
        }
    }

    file = H5Fcreate(argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    dimsf[0] = 157;
    dimsf[1] = 288;
    dataspace = H5Screate_simple(2, dimsf, NULL);

    /*
     * Create a new dataset within the file using defined dataspace and
     * default dataset creation properties.
     */
    dataset = H5Dcreate1(file, "significant wave height", H5T_NATIVE_FLOAT,
            dataspace, H5P_DEFAULT);

    /*
     * Write the data to the dataset using default transfer properties.
     */
    H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, swh);

    /*
     * Close/release resources.
     */
    H5Sclose(dataspace);
    H5Dclose(dataset);
    H5Fclose(file);

    return 0;
}
