/***************************************************************** 
 * File:        tec.c
 *
 * Modification history:        
 *****************************************************************/
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include "hdf5.h"

int main(int argc, char *argv[]) {

    int32_t i, j, k, Y, M, D, H, offset;

    int16_t tec_in[73 * 71];
    char junk[120], outfile[80];
    FILE *fp;

    hsize_t dimsf[2]; /* dataset dimensions */

    /* 
     * Data  and output buffer initialization. 
     */
    hid_t file, dataset; /* handles */
    hid_t dataspace;

    strcpy(outfile, "N20");
    strcpy(junk, basename(argv[1]));

    memcpy(&outfile[3], &junk[9], 2);
    memcpy(&outfile[5], &junk[4], 3);
    memset(&outfile[8], 0, 1);
    if (strncmp(junk, "igrg", 4) == 0)
        strcat(outfile, "00_TEC_IGR_24h.h5");
    if (strncmp(junk, "igsg", 4) == 0)
        strcat(outfile, "00_TEC_IGS_24h.h5");
    if (strncmp(junk, "CODG", 4) == 0)
        strcat(outfile, "00_TEC_CODE_24h.h5");

    printf("Creating: %s\n", outfile);

    file = H5Fcreate(outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Switched 12/01/10 by JMG
    dimsf[0] = 71;
    dimsf[1] = 73;
    dataspace = H5Screate_simple(2, dimsf, NULL);


    fp = fopen(argv[1], "r");
    while (1) {
        if (fgets(junk, 120, fp) == NULL) {
            printf("\"END OF HEADER\" not found.\n");
            exit(1);
        }
        if (strstr(junk, "END OF HEADER") != NULL) break;
    }


    for (i = 0; i < 13; i++) {
        fgets(junk, 120, fp); // START OF TEC MAP
        if (strstr(junk, "START OF TEC MAP") == NULL) {
            printf("\"START OF TEC MAP\" expected but not found.");
            exit(1);
        }

        fgets(junk, 120, fp); // EPOCH OF CURRENT MAP
        if (strstr(junk, "EPOCH OF CURRENT MAP") == NULL) {
            printf("\"EPOCH OF CURRENT MAP\" expected but not found.");
            exit(1);
        }
        sscanf(junk, "%d%d%d%d", &Y, &M, &D, &H);

        if (H == 0 && i == 12) H = 24;

        for (k = 0; k < 71; k++) {
            offset = 73 * k;

            fgets(junk, 120, fp); // LAT/LON1/LON2/DLON/H

            fgets(junk, 120, fp);
            // Note: "%hd" for short int
            for (j = 0; j < 16; j++) sscanf(&junk[5 * j], "%hd", &tec_in[offset + j]);

            fgets(junk, 120, fp);
            for (j = 0; j < 16; j++) sscanf(&junk[5 * j], "%hd", &tec_in[offset + 16 + j]);

            fgets(junk, 120, fp);
            for (j = 0; j < 16; j++) sscanf(&junk[5 * j], "%hd", &tec_in[offset + 32 + j]);

            fgets(junk, 120, fp);
            for (j = 0; j < 16; j++) sscanf(&junk[5 * j], "%hd", &tec_in[offset + 48 + j]);

            fgets(junk, 120, fp);
            for (j = 0; j < 9; j++) sscanf(&junk[5 * j], "%hd", &tec_in[offset + 64 + j]);
        }

        fgets(junk, 120, fp); // END OF TEC MAP
        if (strstr(junk, "END OF TEC MAP") == NULL) {
            printf("\"END OF TEC MAP\" expected but not found.");
            exit(1);
        }

        /*
         * Create a new dataset within the file using defined dataspace and
         * default dataset creation properties.
         */

        sprintf(junk, "%s%02d%s", "tec_", H, "h");
        dataset = H5Dcreate1(file, junk, H5T_NATIVE_USHORT,
                dataspace, H5P_DEFAULT);

        /*
         * Write the data to the dataset using default transfer properties.
         */
        H5Dwrite(dataset, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL,
                H5P_DEFAULT, tec_in);

    }

    /*
     * Close/release resources.
     */
    H5Sclose(dataspace);
    H5Dclose(dataset);
    H5Fclose(file);


    return 0;
}
